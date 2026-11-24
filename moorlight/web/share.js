// Moorlight share page: a read-only viewer for a snapshot written by `moorlight share`.
// Models, map boxes and tour stops are embedded in the page; maps are contoured here (surface nets),
// so no server is needed. Drawing follows the main viewer (app.js): Coot colours, sticks, chicken wire.
import * as THREE from "three";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { edgeTable, triTable } from "three/addons/objects/MarchingCubes.js";

const $ = (s) => document.querySelector(s);
const DATA = JSON.parse($("#ml-data").textContent);
const REDUCED = matchMedia("(prefers-reduced-motion: reduce)").matches;
const ELEMENT_COLORS = { C: 0xe6c229, N: 0x3060ff, O: 0xff2020, S: 0xe0e000, P: 0xff8000, H: 0xdddddd,
  SE: 0xe0a000, FE: 0xe07000, MG: 0x20e020, ZN: 0x8080c0, CA: 0x40c040, NA: 0x8040ff, CL: 0x20f020, K: 0x8040c0,
  MN: 0x9c7ac7, PD: 0x7fc8d8, PT: 0xd0d0e0, CU: 0xc88033, NI: 0x50d050, CO: 0xf090a0, BR: 0xa62929, I: 0x940094,
  F: 0x90e050 };
const CARBON_BY_MOL = [0xe6c229, 0x60c0ff, 0x80e080, 0xff80c0, 0xc0a0ff, 0xffffff];
const CONTOUR = { pos: 0x00d000, neg: 0xe00000, map: 0x3d7fff };
const WATERS = ["HOH", "WAT", "DOD"];

const S = {
  mols: [],            // {data, group, xyz, hideH}
  maps: [],            // {info, values: Float32Array, visible, level, obj, drawn, vol}
  center: new THREE.Vector3(),
  zoom: 30, slab: 9,
  stop: -1,
  labels: new Map(),   // "mol:idx" -> div
  distances: [],
  radius: DATA.contour_radius || 10,
  style: "mesh",       // mesh | surface | slice
  planes: [],          // saved slice planes
};

// ------------------------------------------------------------------ three
const vp = $("#viewport");
const renderer = new THREE.WebGLRenderer({ antialias: true });
const DPR = Math.min(window.devicePixelRatio || 1, 2);
renderer.setPixelRatio(DPR);
vp.prepend(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x000000);
scene.fog = new THREE.Fog(0x000000, 20, 40);
const camera = new THREE.PerspectiveCamera(30, 1, 1, 1000);
const rot = new THREE.Group();     // rotation about the centre
const world = new THREE.Group();   // translation (-centre)
rot.add(world); scene.add(rot);
const overlay = new THREE.Group(); world.add(overlay);
scene.add(new THREE.AmbientLight(0xffffff, 1.1));
const keyLight = new THREE.DirectionalLight(0xffffff, 2.0);
keyLight.position.set(0.4, 0.6, 1);
scene.add(keyLight);
const lineMaterials = new Set();

function lineMat(width, opts = {}) {
  const m = new LineMaterial({ linewidth: width * DPR, worldUnits: false, fog: true, ...opts });
  m.resolution.set(vp.clientWidth * DPR, vp.clientHeight * DPR);
  lineMaterials.add(m);
  return m;
}

function resize() {
  const w = Math.max(1, vp.clientWidth), h = Math.max(1, vp.clientHeight);
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  for (const m of lineMaterials) m.resolution.set(w * DPR, h * DPR);
  redraw();
}
new ResizeObserver(resize).observe(vp);

function updateCamera() {
  camera.position.set(0, 0, S.zoom);
  camera.near = Math.max(0.5, S.zoom - S.slab);
  camera.far = S.zoom + S.slab;
  camera.updateProjectionMatrix();
  scene.fog.near = S.zoom - S.slab * 0.3;
  scene.fog.far = S.zoom + S.slab * 1.05;
  world.position.copy(S.center).multiplyScalar(-1);
}

const crossObj = (() => {
  const g = new THREE.BufferGeometry();
  const s = 0.5;
  g.setAttribute("position", new THREE.Float32BufferAttribute([-s,0,0, s,0,0, 0,-s,0, 0,s,0, 0,0,-s, 0,0,s], 3));
  const o = new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: 0xbbbbbb, fog: false }));
  rot.add(o);
  return o;
})();

let needsRender = true;
function render() {
  requestAnimationFrame(render);
  animate();
  if (!needsRender) return;
  needsRender = false;
  updateCamera();
  sliceFrame();
  renderer.render(scene, camera);
  updateLabels();
}
const redraw = () => (needsRender = true);

function status(text) { $("#status").textContent = text; }

// --------------------------------------------------------------- molecules
function carbonColor(k) {
  const own = S.mols[k]?.data.color;
  return own ? new THREE.Color(own).getHex() : CARBON_BY_MOL[k % CARBON_BY_MOL.length];
}
function elColor(el, k) {
  const mono = S.mols[k]?.data.mono;
  if (mono) return carbonColor(k);
  el = (el || "").toUpperCase();
  if (el === "C") return carbonColor(k);
  return ELEMENT_COLORS[el] ?? 0xff40ff;
}

const STICK_R = 0.13, BALL_R = 0.24, H_STICK_R = 0.07, H_BALL_R = 0.12, WATER_R = 0.28, ION_R = 0.6;
const sphereGeo = new THREE.SphereGeometry(1, 12, 8);
const cylGeo = new THREE.CylinderGeometry(1, 1, 1, 10, 1, true);
const solidMat = new THREE.MeshLambertMaterial({ color: 0xffffff });
const _m = new THREE.Matrix4(), _q = new THREE.Quaternion(), _v = new THREE.Vector3(), _s = new THREE.Vector3();
const _up = new THREE.Vector3(0, 1, 0), _c = new THREE.Color();

function buildMolecule(k) {
  const mol = S.mols[k], d = mol.data;
  for (const ch of [...mol.group.children]) { mol.group.remove(ch); ch.geometry.dispose(); }
  mol.group.visible = d.visible !== false;
  const isH = (i) => d.element[i] === "H";
  const keep = [];
  for (let b = 0; b < d.bonds.length / 2; b++) {
    if (mol.hideH && (isH(d.bonds[2*b]) || isH(d.bonds[2*b+1]))) continue;
    keep.push(b);
  }
  const bonded = new Uint8Array(d.atoms.length);
  for (const i of d.bonds) bonded[i] = 1;
  const lone = [];
  for (let i = 0; i < d.atoms.length; i++) if (!bonded[i] && !(mol.hideH && isH(i))) lone.push(i);
  if (d.style === "lines") {
    const pos = [], col = [];
    for (const b of keep) {
      const i = d.bonds[2*b], j = d.bonds[2*b+1], x = mol.xyz;
      const mx = (x[3*i] + x[3*j]) / 2, my = (x[3*i+1] + x[3*j+1]) / 2, mz = (x[3*i+2] + x[3*j+2]) / 2;
      pos.push(x[3*i], x[3*i+1], x[3*i+2], mx, my, mz, mx, my, mz, x[3*j], x[3*j+1], x[3*j+2]);
      _c.setHex(elColor(d.element[i], k)); col.push(_c.r, _c.g, _c.b, _c.r, _c.g, _c.b);
      _c.setHex(elColor(d.element[j], k)); col.push(_c.r, _c.g, _c.b, _c.r, _c.g, _c.b);
    }
    if (pos.length) {
      const g = new LineSegmentsGeometry(); g.setPositions(pos); g.setColors(col);
      mol.group.add(new LineSegments2(g, lineMat(2.5, { vertexColors: true })));
    }
    lone.length && mol.group.add(...stickMeshes(mol, k, [], lone));
  } else {
    mol.group.add(...stickMeshes(mol, k, keep, lone));
  }
  redraw();
}

function stickMeshes(mol, k, keep, lone) {
  const d = mol.data, x = mol.xyz;
  const P = (i) => new THREE.Vector3(x[3*i], x[3*i+1], x[3*i+2]);
  const isH = (i) => d.element[i] === "H";
  const isWater = (i) => WATERS.includes(d.atoms[i][3]);
  const nb = new Map();
  for (const b of keep) {
    const i = d.bonds[2*b], j = d.bonds[2*b+1];
    if (isH(i) || isH(j)) continue;
    (nb.get(i) || nb.set(i, []).get(i)).push(j);
    (nb.get(j) || nb.set(j, []).get(j)).push(i);
  }
  const cyl = [];
  for (const b of keep) {
    const i = d.bonds[2*b], j = d.bonds[2*b+1];
    const pi = P(i), pj = P(j), mid = pi.clone().add(pj).multiplyScalar(0.5);
    const r = isH(i) || isH(j) ? H_STICK_R : STICK_R;
    const ci = elColor(d.element[i], k), cj = elColor(d.element[j], k);
    cyl.push([pi, mid, r, ci], [mid, pj, r, cj]);
    const order = d.orders ? d.orders[b] : 1;
    if (order === 2 || order === 3 || order === 4) {
      // second line offset towards a neighbour (inside of rings), shortened at both ends, like Coot
      const axis = pj.clone().sub(pi);
      const c = (nb.get(i) || []).find((n) => n !== j) ?? (nb.get(j) || []).find((n) => n !== i);
      let off;
      if (c !== undefined) {
        const ref = (nb.get(i) || []).includes(c) ? pi : pj;
        off = P(c).sub(ref);
        off.sub(axis.clone().multiplyScalar(off.dot(axis) / axis.lengthSq()));
      }
      if (!off || off.lengthSq() < 1e-6) off = new THREE.Vector3(1, 0, 0).cross(axis);
      if (off.lengthSq() < 1e-6) off = new THREE.Vector3(0, 1, 0).cross(axis);
      off.normalize().multiplyScalar(0.27);
      const a = pi.clone().lerp(pj, 0.15).add(off), bb = pj.clone().lerp(pi, 0.15).add(off);
      const m2 = a.clone().add(bb).multiplyScalar(0.5);
      cyl.push([a, m2, r * 0.6, ci], [m2, bb, r * 0.6, cj]);
      if (order === 3) {
        const a3 = pi.clone().lerp(pj, 0.15).sub(off), b3 = pj.clone().lerp(pi, 0.15).sub(off);
        const m3 = a3.clone().add(b3).multiplyScalar(0.5);
        cyl.push([a3, m3, r * 0.6, ci], [m3, b3, r * 0.6, cj]);
      }
    }
  }
  const out = [];
  if (cyl.length) {
    const im = new THREE.InstancedMesh(cylGeo, solidMat, cyl.length);
    cyl.forEach(([p, q, r, col], n) => {
      _v.subVectors(q, p);
      const len = _v.length();
      _q.setFromUnitVectors(_up, _v.normalize());
      _m.compose(p.clone().add(q).multiplyScalar(0.5), _q, _s.set(r, len, r));
      im.setMatrixAt(n, _m);
      im.setColorAt(n, _c.setHex(col));
    });
    im.frustumCulled = false;
    out.push(im);
  }
  const loneSet = new Set(lone);
  const balls = [];
  if (d.style === "lines") balls.push(...lone);
  else for (let i = 0; i < d.atoms.length; i++) if (!(mol.hideH && isH(i))) balls.push(i);
  if (balls.length) {
    const im = new THREE.InstancedMesh(sphereGeo, solidMat, balls.length);
    balls.forEach((i, n) => {
      let r = isH(i) ? H_BALL_R : BALL_R;
      if (loneSet.has(i)) r = isWater(i) ? WATER_R : (isH(i) ? H_BALL_R : ION_R);
      _m.compose(P(i), _q.identity(), _s.set(r, r, r));
      im.setMatrixAt(n, _m);
      im.setColorAt(n, _c.setHex(elColor(d.element[i], k)));
    });
    im.frustumCulled = false;
    out.push(im);
  }
  return out;
}

function atomLabel(k, i) {
  const a = S.mols[k].data.atoms[i];
  return `${a[3]}${a[1]}${a[2]} ${a[4]}${a[5] ? ":" + a[5] : ""}${a[9] ? " [" + a[9] + "]" : ""}`;
}
function atomInfo(k, i) {
  const a = S.mols[k].data.atoms[i], x = S.mols[k].xyz;
  const p = [x[3*i], x[3*i+1], x[3*i+2]];
  const vals = S.maps.filter((m) => m.visible).map((m) => {
    const v = valueAt(m, p);
    return Number.isFinite(v) ? `${m.info.name} ${fmtSigned(v)}${unit(m)}` : null;
  }).filter(Boolean);
  const mol = S.mols.length > 1 ? `  [${S.mols[k].data.name}]` : "";
  return `${a[0]}/${atomLabel(k, i)}  occ ${a[6].toFixed(2)}  B ${a[7].toFixed(1)}${mol}` + (vals.length ? "\n" + vals.join(" · ") : "");
}
const fmtSigned = (v) => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(1);

// ------------------------------------------------------------ picking, labels
function pick(clientX, clientY, maxPx) {
  updateCamera(); scene.updateMatrixWorld();
  const rect = renderer.domElement.getBoundingClientRect();
  const v = new THREE.Vector3();
  let best = null, bd = maxPx * maxPx;
  S.mols.forEach((mol, k) => {
    if (!mol.group.visible) return;
    const el = mol.data.element;
    for (let i = 0; i < mol.data.atoms.length; i++) {
      if (mol.hideH && el[i] === "H") continue;
      v.set(mol.xyz[3*i], mol.xyz[3*i+1], mol.xyz[3*i+2]).applyMatrix4(world.matrixWorld).project(camera);
      if (v.z < -1 || v.z > 1) continue;
      const sx = (v.x + 1) / 2 * rect.width, sy = (1 - v.y) / 2 * rect.height;
      const dd = (sx - (clientX - rect.left)) ** 2 + (sy - (clientY - rect.top)) ** 2;
      if (dd < bd || (best && Math.abs(dd - bd) < 4 && v.z < best.z)) { bd = dd; best = { mol: k, idx: i, z: v.z }; }
    }
  });
  return best;
}
function toggleLabel(k, i) {
  const key = `${k}:${i}`;
  if (S.labels.has(key)) { S.labels.get(key).remove(); S.labels.delete(key); }
  else {
    const div = document.createElement("div");
    div.className = "label"; div.textContent = atomLabel(k, i);
    $("#labels").appendChild(div); S.labels.set(key, div);
  }
  redraw();
}
function clearLabels() {
  for (const div of S.labels.values()) div.remove();
  S.labels.clear();
  for (const d of S.distances) { d.div.remove(); overlay.remove(d.line); d.line.geometry.dispose(); }
  S.distances = [];
  redraw();
}
function addDistance(a, b) {
  const A = S.mols[a[0]].xyz, B = S.mols[b[0]].xyz;
  const pa = new THREE.Vector3(A[3*a[1]], A[3*a[1]+1], A[3*a[1]+2]);
  const pb = new THREE.Vector3(B[3*b[1]], B[3*b[1]+1], B[3*b[1]+2]);
  const g = new THREE.BufferGeometry().setFromPoints([pa, pb]);
  const line = new THREE.LineSegments(g, new THREE.LineDashedMaterial({ color: 0x99ffff, dashSize: 0.2, gapSize: 0.15 }));
  line.computeLineDistances();
  overlay.add(line);
  const div = document.createElement("div");
  div.className = "label dist"; div.textContent = `${pa.distanceTo(pb).toFixed(2)} Å`;
  $("#labels").appendChild(div);
  S.distances.push({ a, b, pa, pb, line, div });
}
function updateLabels() {
  const rect = renderer.domElement.getBoundingClientRect();
  const v = new THREE.Vector3();
  scene.updateMatrixWorld();
  for (const [key, div] of S.labels) {
    const [k, i] = key.split(":").map(Number);
    const x = S.mols[k].xyz;
    v.set(x[3*i], x[3*i+1], x[3*i+2]).applyMatrix4(world.matrixWorld).project(camera);
    div.style.display = v.z > -1 && v.z < 1 && S.mols[k].group.visible ? "" : "none";
    div.style.left = `${(v.x + 1) / 2 * rect.width}px`;
    div.style.top = `${(1 - v.y) / 2 * rect.height}px`;
  }
  for (const d of S.distances) {
    const pa = d.pa.clone().applyMatrix4(world.matrixWorld).project(camera);
    const pb = d.pb.clone().applyMatrix4(world.matrixWorld).project(camera);
    d.div.style.left = `${((pa.x + pb.x) / 2 + 1) / 2 * rect.width}px`;
    d.div.style.top = `${(1 - (pa.y + pb.y) / 2) / 2 * rect.height + 11}px`;
  }
}

// --------------------------------------------------------------------- maps
async function decodeMap(info) {
  const bin = atob(info.data);
  const bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
  const q = new Int16Array(await new Response(stream).arrayBuffer());
  const f = new Float32Array(q.length);
  for (let i = 0; i < q.length; i++) f[i] = q[i] * info.scale;
  return f;
}

const unit = (m) => (m.info.units === "sigma" ? "σ" : (m.info.unit_label ? " " + m.info.unit_label : ""));
// map value in the units the level is given in (σ of the whole map, or the map's own units)
const toLevelUnits = (m, v) => (m.info.units === "sigma" ? (v - m.info.mean) / m.info.rms : v);

// Stored grids: position = origin + G (i, j, k). "native" maps are the map's own grid over the unit cell
// (periodic, wrap); "box" maps are resampled on an orthogonal box.
const mat3 = (M, v) => [M[0]*v[0] + M[1]*v[1] + M[2]*v[2], M[3]*v[0] + M[4]*v[1] + M[5]*v[2], M[6]*v[0] + M[7]*v[1] + M[8]*v[2]];
const wrapIndex = (g, n) => ((g % n) + n) % n;

function valueAt(m, p) {
  const { n, origin: o, Ginv, wrap } = m.info;
  const f = mat3(Ginv, [p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
  if (!wrap && f.some((v, a) => v < 0 || v > n[a] - 1)) return NaN;
  const i = f.map((v, a) => (wrap ? Math.floor(v) : Math.min(n[a] - 2, Math.floor(v))));
  const t = f.map((v, a) => v - i[a]);
  const g = (a, d) => (wrap ? wrapIndex(i[a] + d, n[a]) : i[a] + d);
  const D = m.values, syz = n[1] * n[2];
  let s = 0;
  for (let c = 0; c < 8; c++) {
    const dx = c & 1, dy = (c >> 1) & 1, dz = (c >> 2) & 1;
    const w = (dx ? t[0] : 1 - t[0]) * (dy ? t[1] : 1 - t[1]) * (dz ? t[2] : 1 - t[2]);
    s += w * D[g(0, dx) * syz + g(1, dy) * n[2] + g(2, dz)];
  }
  return toLevelUnits(m, s);
}

// Marching cubes (three.js tables): vertices where the contour crosses grid edges, shared between neighbouring
// cells. Mesh style draws the triangle edges (Coot's chicken wire), surface style the shaded triangles.
const MC_CORNER = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
const MC_EDGE = [[0, 0, 0, 0], [1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 1, 1],
  [0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 2], [1, 0, 0, 2], [1, 1, 0, 2], [0, 1, 0, 2]];   // start point, axis
const MCX = MC_CORNER.map((q) => q[0]), MCY = MC_CORNER.map((q) => q[1]), MCZ = MC_CORNER.map((q) => q[2]);

function marchingCubes(m, c, r, level, sign, wire) {
  const { n, origin: o, G, Ginv, wrap } = m.info;
  const D = m.values;
  const off = m.info.units === "sigma" ? m.info.mean : 0, scale = m.info.units === "sigma" ? m.info.rms : 1;
  const L = level * scale;
  // grid cells covering the cube c +/- r (a parallelepiped on a non-orthogonal grid); triangles leaving the cube
  // are dropped, so the contour fills a box around the centre as in Coot
  const ci = mat3(Ginv, [c[0] - o[0], c[1] - o[1], c[2] - o[2]]);
  const lo = [], hi = [];
  for (let a = 0; a < 3; a++) {
    const ext = r * (Math.abs(Ginv[3*a]) + Math.abs(Ginv[3*a+1]) + Math.abs(Ginv[3*a+2]));
    let l = Math.floor(ci[a] - ext), h = Math.ceil(ci[a] + ext);
    if (!wrap) { l = Math.max(0, l); h = Math.min(n[a] - 1, h); }
    lo.push(l); hi.push(h);
  }
  const sx = hi[0] - lo[0], sy = hi[1] - lo[1], sz = hi[2] - lo[2];
  const empty = { verts: new Float32Array(0), index: new Uint32Array(0), segments: new Float32Array(0) };
  if (sx <= 0 || sy <= 0 || sz <= 0) return empty;
  // offset in the value array of each grid plane along each axis (wrapped on the periodic cell grid)
  const stride = [n[1] * n[2], n[2], 1];
  const [X, Y, Z] = [0, 1, 2].map((a) => {
    const t = new Int32Array(hi[a] - lo[a] + 1);
    for (let q = 0; q < t.length; q++) t[q] = (wrap ? wrapIndex(lo[a] + q, n[a]) : lo[a] + q) * stride[a];
    return t;
  });
  const py = sy + 1, pz = sz + 1;
  const edgeVert = new Int32Array((sx + 1) * py * pz * 3).fill(-1);
  const verts = [], inBox = [], tris = [];
  const val = (a, b, d) => sign * (D[X[a] + Y[b] + Z[d]] - off);
  function edgeVertex(a, b, d, axis) {
    const key = ((a * py + b) * pz + d) * 3 + axis;
    let v = edgeVert[key];
    if (v >= 0) return v;
    const v0 = val(a, b, d);
    const v1 = axis === 0 ? val(a + 1, b, d) : axis === 1 ? val(a, b + 1, d) : val(a, b, d + 1);
    const t = (L - v0) / (v1 - v0);
    const fi = lo[0] + a + (axis === 0 ? t : 0), fj = lo[1] + b + (axis === 1 ? t : 0), fk = lo[2] + d + (axis === 2 ? t : 0);
    const x = o[0] + G[0]*fi + G[1]*fj + G[2]*fk, y = o[1] + G[3]*fi + G[4]*fj + G[5]*fk, z = o[2] + G[6]*fi + G[7]*fj + G[8]*fk;
    v = verts.length / 3;
    verts.push(x, y, z);
    inBox.push(Math.abs(x - c[0]) <= r && Math.abs(y - c[1]) <= r && Math.abs(z - c[2]) <= r);
    edgeVert[key] = v;
    return v;
  }
  const ev = new Int32Array(12);
  for (let a = 0; a < sx; a++) for (let b = 0; b < sy; b++) for (let d = 0; d < sz; d++) {
    let cube = 0;
    for (let q = 0; q < 8; q++) if (val(a + MCX[q], b + MCY[q], d + MCZ[q]) < L) cube |= 1 << q;
    const bits = edgeTable[cube];
    if (!bits) continue;
    for (let e = 0; e < 12; e++) {
      if (bits & (1 << e)) { const E = MC_EDGE[e]; ev[e] = edgeVertex(a + E[0], b + E[1], d + E[2], E[3]); }
    }
    for (let t = cube * 16; triTable[t] !== -1; t += 3) {
      const i = ev[triTable[t]], j = ev[triTable[t + 1]], k = ev[triTable[t + 2]];
      if (inBox[i] && inBox[j] && inBox[k]) tris.push(i, j, k);
    }
  }
  const out = { verts: new Float32Array(verts), index: new Uint32Array(tris), segments: new Float32Array(0) };
  if (wire) {   // each triangle edge once
    const seen = new Set(), seg = [];
    const add = (i, j) => {
      const key = i < j ? i * 67108864 + j : j * 67108864 + i;
      if (seen.has(key)) return;
      seen.add(key);
      seg.push(verts[3*i], verts[3*i+1], verts[3*i+2], verts[3*j], verts[3*j+1], verts[3*j+2]);
    };
    for (let t = 0; t < tris.length; t += 3) { add(tris[t], tris[t+1]); add(tris[t+1], tris[t+2]); add(tris[t+2], tris[t]); }
    out.segments = new Float32Array(seg);
  }
  return out;
}

const mapColor = (m, sign) => (m.info.kind === "diff" ? (sign > 0 ? CONTOUR.pos : CONTOUR.neg) : (m.info.color || CONTOUR.map));
function surfaceMat(color) {
  return new THREE.MeshPhongMaterial({ color, transparent: true, opacity: 0.5, side: THREE.DoubleSide, depthWrite: false,
    shininess: 40, specular: 0x333333 });
}
function clearMapObj(m) {
  for (const ch of [...m.obj.children]) { m.obj.remove(ch); ch.geometry.dispose(); }
}
function contourMap(m) {
  const c = S.center.toArray();
  clearMapObj(m);
  const surface = S.style === "surface";
  m.style = S.style === "slice" ? "mesh" : S.style;
  for (const sign of m.info.kind === "diff" ? [1, -1] : [1]) {
    const k = sign > 0 ? 0 : 1;
    const mc = marchingCubes(m, c, S.radius, m.level, sign, !surface);
    if (surface) {
      if (!mc.index.length) continue;
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.BufferAttribute(mc.verts, 3));
      g.setIndex(new THREE.BufferAttribute(mc.index, 1));
      g.computeVertexNormals();
      m.obj.add(new THREE.Mesh(g, m.smats[k] || (m.smats[k] = surfaceMat(mapColor(m, sign)))));
    } else {
      if (!mc.segments.length) continue;
      const g = new LineSegmentsGeometry();
      g.setPositions(mc.segments);
      m.obj.add(new LineSegments2(g, m.mats[k] || (m.mats[k] = lineMat(1.1, { color: mapColor(m, sign) }))));
    }
  }
  m.drawn = { center: S.center.clone(), level: m.level, style: m.style, radius: S.radius };
  redraw();
}
let contourTimer = null;
function scheduleContour(now = false) {
  clearTimeout(contourTimer);
  contourTimer = setTimeout(() => {
    const sliced = S.style === "slice" ? activeMap() : null;
    for (const m of S.maps) {
      m.obj.visible = m.visible && m !== sliced;
      if (!m.visible || m === sliced) continue;
      const d = m.drawn;
      const want = S.style === "slice" ? "mesh" : S.style;
      if (d && d.center.distanceTo(S.center) < 0.5 && d.level === m.level && d.style === want && d.radius === S.radius) continue;
      contourMap(m);
    }
    updateSlice();
    redraw();
  }, now ? 0 : 120);
}

// -------------------------------------------------------------------- slice
// Map section as in Olex2: a plane through the focus that always faces the viewer (so it turns with the
// molecule), at an adjustable depth. It shows the active map (the first visible one): green +, red - (blue for
// 2mFo-DFc-like maps), stronger towards 3 x level, with contour lines at 1, 2, 4, 8 ... x level. The plane is dark
// glass: atoms behind it are dimmed, atoms in front stay bright.
const SLICE = { mesh: null, halo: null, map: null, failed: false, depth: 0, tiltX: 0, tiltY: 0, thick: 0, lock: null,
  off: false };   // off: live plane hidden (saved planes stay)   // lock: {c, Rt}
const SLICE_VS = `
uniform mat4 uPlane;
varying vec3 vModel;
varying vec2 vXY;
void main() {
  vModel = (uPlane * vec4(position, 1.0)).xyz;
  vXY = position.xy;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}`;
const SLICE_FS = `
precision highp float;
precision highp sampler3D;
uniform sampler3D uMap;
uniform mat3 uGinv;
uniform vec3 uOrigin, uN, uPos, uNeg, uNormal;
uniform float uLevel, uThick;
uniform int uSigned;
varying vec3 vModel;
varying vec2 vXY;
float sampleAt(vec3 p) {
  vec3 idx = uGinv * (p - uOrigin);
  return texture(uMap, (idx.zyx + 0.5) / uN.zyx).r;
}
void main() {
  // thin plane: the value on it; thick slice: the strongest value (sign kept) within +/- uThick / 2 of it
  float v = sampleAt(vModel);
  if (uThick > 0.01) {
    for (int s = 0; s < 13; s++) {
      float w = sampleAt(vModel + uNormal * uThick * (float(s) / 12.0 - 0.5));
      if (uSigned == 0 ? w > v : abs(w) > abs(v)) v = w;
    }
  }
  float a = (uSigned == 0 && v < 0.0) ? 0.0 : abs(v);
  vec3 hue = v >= 0.0 ? uPos : uNeg;
  float k = clamp(a / (3.0 * uLevel), 0.0, 1.0);
  vec3 col = hue * k;
  float alpha = 0.45 + 0.5 * k;
  float g = log2(max(a, 1e-6) / uLevel);
  float line = (1.0 - smoothstep(0.0, 1.5 * fwidth(g), abs(fract(g + 0.5) - 0.5))) * smoothstep(0.85, 0.95, a / uLevel);
  col = mix(col, min(hue * 1.2 + 0.3, vec3(1.0)), line);
  alpha = max(alpha, line);
  gl_FragColor = vec4(col * alpha, alpha);
}`;

function mapTexture(m) {
  if (m.tex) return m.tex;
  const [n0, n1, n2] = m.info.n;
  const half = new Uint16Array(m.values.length);
  for (let q = 0; q < half.length; q++) half[q] = THREE.DataUtils.toHalfFloat(toLevelUnits(m, m.values[q]));
  const tex = new THREE.Data3DTexture(half, n2, n1, n0);   // x = third grid axis (fastest), z = first
  tex.format = THREE.RedFormat;
  tex.type = THREE.HalfFloatType;
  tex.minFilter = tex.magFilter = THREE.LinearFilter;
  tex.wrapS = tex.wrapT = tex.wrapR = m.info.wrap ? THREE.RepeatWrapping : THREE.ClampToEdgeWrapping;
  tex.unpackAlignment = 1;
  tex.needsUpdate = true;
  return (m.tex = tex);
}
// a slice plane: unit square with the map shader and an outline; its matrix (model coordinates) places it
function makePlaneMesh(outlineColor) {
  const mat = new THREE.ShaderMaterial({
    vertexShader: SLICE_VS, fragmentShader: SLICE_FS, side: THREE.DoubleSide, transparent: true,
    premultipliedAlpha: true, depthWrite: false,
    uniforms: {
      uMap: { value: null }, uPlane: { value: new THREE.Matrix4() }, uGinv: { value: new THREE.Matrix3() },
      uOrigin: { value: new THREE.Vector3() }, uN: { value: new THREE.Vector3() }, uPos: { value: new THREE.Color() },
      uNeg: { value: new THREE.Color() }, uLevel: { value: 1 }, uSigned: { value: 1 },
      uNormal: { value: new THREE.Vector3(0, 0, 1) }, uThick: { value: 0 },
    },
  });
  const mesh = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), mat);
  mesh.matrixAutoUpdate = false;
  mesh.frustumCulled = false;
  mesh.renderOrder = 5;
  const frame = new LineSegmentsGeometry();
  frame.setPositions([-1, -1, 0, 1, -1, 0, 1, -1, 0, 1, 1, 0, 1, 1, 0, -1, 1, 0, -1, 1, 0, -1, -1, 0]);
  const outline = new LineSegments2(frame, lineMat(2, { color: outlineColor }));
  outline.frustumCulled = false;
  mesh.add(outline);
  world.add(mesh);
  return mesh;
}
function setPlaneMap(u, m) {
  u.uMap.value = mapTexture(m);
  u.uGinv.value.set(...m.info.Ginv);
  u.uOrigin.value.set(...m.info.origin);
  u.uN.value.set(...m.info.n);
  u.uPos.value.set(m.info.kind === "diff" ? CONTOUR.pos : (m.info.color || CONTOUR.map));
  u.uNeg.value.set(CONTOUR.neg);
  u.uSigned.value = m.info.kind === "diff" ? 1 : 0;
  u.uLevel.value = m.level;
}
function updateSlice() {
  const m = S.style === "slice" && !SLICE.off ? activeMap() : null;
  if (!m) {
    if (SLICE.mesh) SLICE.mesh.visible = false;
    if (SLICE.halo) SLICE.halo.visible = false;
    return;
  }
  if (!SLICE.mesh) {
    SLICE.mesh = makePlaneMesh(0xe6c229);
    // atoms the plane cuts through
    const total = S.mols.reduce((t, x) => t + x.data.atoms.length, 0);
    SLICE.halo = new THREE.InstancedMesh(sphereGeo, new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true,
      opacity: 0.45, depthWrite: false }), Math.max(1, total));
    SLICE.halo.count = 0;
    SLICE.halo.frustumCulled = false;
    world.add(SLICE.halo);
  }
  if (SLICE.map !== m) { SLICE.map = m; setPlaneMap(SLICE.mesh.material.uniforms, m); }
  SLICE.mesh.material.uniforms.uLevel.value = m.level;
  SLICE.mesh.visible = SLICE.halo.visible = true;
}
// where the live plane is: facing the viewer through the focus (or where it was fixed), tilted, at the chosen depth
// along its normal. Fills _Rt (axes as columns), _M (placement of the unit square), _n, _ex, _ey, _p0.
const _M = new THREE.Matrix4(), _T = new THREE.Matrix4(), _Rt = new THREE.Matrix4();
const _n = new THREE.Vector3(), _ex = new THREE.Vector3(), _ey = new THREE.Vector3(), _p0 = new THREE.Vector3();
function livePlacement() {
  const c = SLICE.lock ? SLICE.lock.c : S.center;
  if (SLICE.lock) _Rt.copy(SLICE.lock.Rt);
  else _Rt.makeRotationFromQuaternion(rot.quaternion).transpose();
  _Rt.multiply(_T.makeRotationX(THREE.MathUtils.degToRad(SLICE.tiltX))).multiply(_T.makeRotationY(THREE.MathUtils.degToRad(SLICE.tiltY)));
  const r = S.radius;
  _M.makeTranslation(c.x, c.y, c.z).multiply(_Rt).multiply(_T.makeTranslation(0, 0, SLICE.depth))
    .multiply(_T.makeScale(r, r, 1));
  _ex.setFromMatrixColumn(_Rt, 0); _ey.setFromMatrixColumn(_Rt, 1); _n.setFromMatrixColumn(_Rt, 2);
  _p0.copy(c).addScaledVector(_n, SLICE.depth);
  return r;
}
// per frame: place the live plane; atoms inside it get a halo; saved planes follow their map's level
function sliceFrame() {
  for (const pl of S.planes) {
    pl.mesh.visible = pl.visible;
    pl.mesh.material.uniforms.uLevel.value = pl.level ?? pl.map.level;
  }
  if (!SLICE.mesh || !SLICE.mesh.visible) return;
  const r = livePlacement();
  SLICE.mesh.matrix.copy(_M);
  SLICE.mesh.matrixWorldNeedsUpdate = true;
  const u = SLICE.mesh.material.uniforms;
  u.uPlane.value.copy(_M);
  u.uNormal.value.copy(_n).normalize();
  u.uThick.value = SLICE.thick;
  const half = Math.max(0.5, SLICE.thick / 2);
  let k = 0;
  for (const mol of S.mols) {
    if (!mol.group.visible) continue;
    const x = mol.xyz, el = mol.data.element;
    for (let i = 0; i < el.length; i++) {
      if (el[i] === "H") continue;
      const dx = x[3*i] - _p0.x, dy = x[3*i+1] - _p0.y, dz = x[3*i+2] - _p0.z;
      if (Math.abs(dx * _n.x + dy * _n.y + dz * _n.z) > half) continue;
      if (Math.abs(dx * _ex.x + dy * _ex.y + dz * _ex.z) > r || Math.abs(dx * _ey.x + dy * _ey.y + dz * _ey.z) > r) continue;
      _m.compose(_v.set(x[3*i], x[3*i+1], x[3*i+2]), _q.identity(), _s.set(0.42, 0.42, 0.42));
      SLICE.halo.setMatrixAt(k++, _m);
    }
  }
  SLICE.halo.count = k;
  SLICE.halo.instanceMatrix.needsUpdate = true;
}

// ------------------------------------------------------------ saved planes
// A saved plane stays in the molecule with its own outline colour and is shown in every map style; the live
// plane is free for the next one. Each viewer's planes are kept in this browser between visits.
const PLANE_COLORS = [0x4fd0ff, 0xff6bd6, 0xffa64d, 0xa0ff70, 0xffffff, 0xb59cff];
const PLANE_KEY = "moorlight-share-planes";
let planeCount = 0;
function addSavedPlane(pl) {
  pl.mesh = makePlaneMesh(pl.color);
  const M = new THREE.Matrix4().makeTranslation(pl.p0.x, pl.p0.y, pl.p0.z).multiply(pl.Rt)
    .multiply(new THREE.Matrix4().makeScale(pl.r, pl.r, 1));
  pl.mesh.matrix.copy(M);
  pl.mesh.matrixWorldNeedsUpdate = true;
  const u = pl.mesh.material.uniforms;
  u.uPlane.value.copy(M);
  setPlaneMap(u, pl.map);
  u.uNormal.value.setFromMatrixColumn(pl.Rt, 2).normalize();
  u.uThick.value = pl.thick;
  S.planes.push(pl);
}
function savePlane() {
  if (S.style !== "slice" || !SLICE.mesh || !SLICE.mesh.visible) {
    status(SLICE.off ? "Turn the live plane on to place a new one." : "Switch to Slice and place the plane first.");
    return;
  }
  const r = livePlacement();
  planeCount += 1;
  addSavedPlane({ name: `Plane ${planeCount}`, color: PLANE_COLORS[(planeCount - 1) % PLANE_COLORS.length],
    p0: _p0.clone(), Rt: _Rt.clone(), r, thick: SLICE.thick, map: SLICE.map, level: null, visible: true });
  renderPlaneList(); storePlanes(); redraw();
  status(`Plane ${planeCount} saved. Move the live plane to place the next one.`);
}
// planes the page was made with (spec "planes"); they come back on every visit
function loadPresetPlanes() {
  for (const p of DATA.planes || []) {
    const map = S.maps.find((m) => m.info.name === p.map) || S.maps.find((m) => m.info.visible) || S.maps[0];
    if (!map) return;
    const Rt = new THREE.Matrix4().makeBasis(new THREE.Vector3(...p.ex), new THREE.Vector3(...p.ey), new THREE.Vector3(...p.n));
    addSavedPlane({ name: p.name, color: p.color ?? PLANE_COLORS[S.planes.length % PLANE_COLORS.length],
      p0: new THREE.Vector3(...p.p0), Rt, r: p.r, thick: p.thick, map, level: p.level ?? null,
      visible: p.visible !== false, preset: true });
  }
}
// turn the view onto a saved plane: centred on it, facing it
function goToPlane(i) {
  const pl = S.planes[i];
  pl.visible = true;
  renderPlaneList();
  animateTo(pl.p0.clone(), new THREE.Quaternion().setFromRotationMatrix(pl.Rt.clone().transpose()));
  status(`${pl.name}${pl.level != null ? `, contoured at ${pl.level}` : ""}`);
}
// a parallel copy of a plane, moved along its normal by the copy offset; it belongs to the viewer
function copyPlane(i) {
  const pl = S.planes[i];
  const off = +$("#copy-offset").value || 0;
  const n = new THREE.Vector3().setFromMatrixColumn(pl.Rt, 2).normalize();
  planeCount += 1;
  addSavedPlane({ name: `${pl.name} ${off >= 0 ? "+" : "\u2212"}${Math.abs(off).toFixed(1)} Å`,
    color: PLANE_COLORS[(planeCount - 1) % PLANE_COLORS.length], p0: pl.p0.clone().addScaledVector(n, off),
    Rt: pl.Rt.clone(), r: pl.r, thick: pl.thick, map: pl.map, level: pl.level, visible: true });
  renderPlaneList(); storePlanes(); redraw();
  status(`Copied ${pl.name}, ${off.toFixed(1)} Å along its normal.`);
}
function removePlane(i) {
  const [pl] = S.planes.splice(i, 1);
  world.remove(pl.mesh);
  pl.mesh.geometry.dispose();
  renderPlaneList(); storePlanes(); redraw();
}
function editPlane(i) {        // the saved plane becomes the live one again
  const pl = S.planes[i];
  removePlane(i);
  if (S.style !== "slice") setStyle("slice");
  if (SLICE.off) setLivePlane(true);
  SLICE.lock = { c: pl.p0.clone(), Rt: pl.Rt.clone() };
  for (const b of [$("#slice-lock"), $("#mini-lock")]) { b.setAttribute("aria-pressed", "true"); b.textContent = "Release"; }
  setTilt("tiltX", 0); setTilt("tiltY", 0); setSliceDepth(0); setThickness(pl.thick);
  status(`${pl.name} is the live plane again.`);
}
function renderPlaneList() {
  const ul = $("#plane-list");
  ul.innerHTML = "";
  S.planes.forEach((pl, i) => {
    const li = document.createElement("li");
    li.className = "plane";
    const chip = "#" + pl.color.toString(16).padStart(6, "0");
    li.innerHTML = `<input type="checkbox" id="plane-vis-${i}"><label for="plane-vis-${i}"><span class="chip" ` +
      `style="border-color:${chip}" aria-hidden="true"></span></label>` +
      `<button type="button" class="go"></button>` +
      `<button type="button">Copy</button><button type="button">Edit</button>` +
      `<button type="button" aria-label="Delete">&times;</button>`;
    const go = li.querySelector(".go");
    go.textContent = pl.name + (pl.preset ? "" : " (yours)");
    go.title = `Look at ${pl.name}: ${pl.map.info.name}${pl.thick ? `, ${pl.thick.toFixed(1)} Å thick` : ""}`;
    go.addEventListener("click", () => goToPlane(i));
    const cb = li.querySelector("input");
    cb.checked = pl.visible;
    cb.addEventListener("change", () => { pl.visible = cb.checked; storePlanes(); redraw(); });
    const [, copy, edit, del] = li.querySelectorAll("button");
    copy.setAttribute("aria-label", `Copy ${pl.name} at the copy offset`);
    copy.addEventListener("click", () => copyPlane(i));
    edit.addEventListener("click", () => editPlane(i));
    del.setAttribute("aria-label", `Delete ${pl.name}`);
    del.addEventListener("click", () => removePlane(i));
    ul.appendChild(li);
  });
  $("#plane-list").hidden = !S.planes.length;
  $("#planes-empty").hidden = S.planes.length > 0;
  $("#copy-row").hidden = !S.planes.length;
}
function storePlanes() {
  try {
    localStorage.setItem(PLANE_KEY, JSON.stringify({ count: planeCount, planes: S.planes.filter((pl) => !pl.preset).map((pl) => ({ name: pl.name,
      color: pl.color, p0: pl.p0.toArray(), Rt: pl.Rt.toArray(), r: pl.r, thick: pl.thick, map: pl.map.info.name,
      level: pl.level, visible: pl.visible })) }));
  } catch (err) { /* storage blocked */ }
}
function loadPlanes() {
  let saved = null;
  try { saved = JSON.parse(localStorage.getItem(PLANE_KEY) || "null"); } catch (err) { saved = null; }
  if (!saved || !Array.isArray(saved.planes)) return;
  planeCount = saved.count || 0;
  for (const p of saved.planes) {
    const map = S.maps.find((m) => m.info.name === p.map);
    if (!map) continue;
    addSavedPlane({ name: p.name, color: p.color, p0: new THREE.Vector3(...p.p0), Rt: new THREE.Matrix4().fromArray(p.Rt),
      r: p.r, thick: p.thick, map, level: p.level ?? null, visible: p.visible !== false });
  }
}
function setTilt(axis, deg) {
  SLICE[axis] = Math.max(-90, Math.min(90, deg));
  const id = axis === "tiltX" ? "slice-tilt-x" : "slice-tilt-y";
  $(`#${id}`).value = String(SLICE[axis]);
  $(`#${id}-out`).textContent = `${Math.round(SLICE[axis])}°`;
  redraw();
}
function setLivePlane(on) {
  SLICE.off = !on;
  for (const b of [$("#slice-live"), $("#mini-live")]) {
    b.setAttribute("aria-pressed", String(on));
    b.textContent = on ? "Live plane on" : "Live plane off";
  }
  $("#app").classList.toggle("live-off", !on);
  updateSlice();
  redraw();
}
function setThickness(t) {
  SLICE.thick = Math.max(0, Math.min(4, t));
  $("#slice-thick").value = String(SLICE.thick);
  $("#slice-thick-out").textContent = SLICE.thick ? `${SLICE.thick.toFixed(1)} Å` : "thin";
  redraw();
}
function faceMe() {
  if (SLICE.lock) setSliceLock(false);
  setTilt("tiltX", 0); setTilt("tiltY", 0); setSliceDepth(0);
}
// best plane through the labelled atoms (3 or more): the plane is fixed there and the view turns to face it
function planeThroughLabels() {
  const pts = [];
  for (const key of S.labels.keys()) {
    const [k, i] = key.split(":").map(Number);
    const x = S.mols[k].xyz;
    pts.push(new THREE.Vector3(x[3*i], x[3*i+1], x[3*i+2]));
  }
  if (pts.length < 3) { status("Label three or more atoms first (tap them), then put the plane through them."); return; }
  const c = pts.reduce((a, p) => a.add(p), new THREE.Vector3()).multiplyScalar(1 / pts.length);
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of pts) {
    const d = [p.x - c.x, p.y - c.y, p.z - c.z];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += d[a] * d[b];
  }
  const { vectors } = eigenSym3(C);          // columns, ascending eigenvalues
  const n = new THREE.Vector3(...vectors.map((r) => r[0])).normalize();
  const ex = new THREE.Vector3(...vectors.map((r) => r[2])).normalize();
  const ey = new THREE.Vector3().crossVectors(n, ex).normalize();
  const Rt = new THREE.Matrix4().makeBasis(ex, ey, n);
  SLICE.lock = { c, Rt };
  for (const b of [$("#slice-lock"), $("#mini-lock")]) { b.setAttribute("aria-pressed", "true"); b.textContent = "Release"; }
  setTilt("tiltX", 0); setTilt("tiltY", 0); setSliceDepth(0);
  if (S.style !== "slice") setStyle("slice");
  animateTo(c, new THREE.Quaternion().setFromRotationMatrix(Rt.clone().transpose()));
  status(`Plane through ${pts.length} labelled atoms (rms out of plane ${Math.sqrt(Math.max(0, eigenSym3(C).values[0]) / pts.length).toFixed(2)} Å).`);
}
// eigen decomposition of a symmetric 3x3 matrix (Jacobi); returns values ascending and vectors as columns
function eigenSym3(A) {
  const a = A.map((r) => r.slice()), v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < 2; p++) for (let q = p + 1; q < 3; q++) off += a[p][q] * a[p][q];
    if (off < 1e-20) break;
    for (let p = 0; p < 2; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(a[p][q]) < 1e-30) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
      for (let k = 0; k < 3; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = cs * akp - sn * akq; a[k][q] = sn * akp + cs * akq;
      }
      for (let k = 0; k < 3; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = cs * apk - sn * aqk; a[q][k] = sn * apk + cs * aqk;
      }
      for (let k = 0; k < 3; k++) {
        const vkp = v[k][p], vkq = v[k][q];
        v[k][p] = cs * vkp - sn * vkq; v[k][q] = sn * vkp + cs * vkq;
      }
    }
  }
  const order = [0, 1, 2].sort((i, j) => a[i][i] - a[j][j]);
  return { values: order.map((i) => a[i][i]), vectors: v.map((r) => order.map((i) => r[i])) };
}
function setSliceLock(on) {
  SLICE.lock = on ? { c: S.center.clone(), Rt: new THREE.Matrix4().makeRotationFromQuaternion(rot.quaternion).transpose() } : null;
  for (const b of [$("#slice-lock"), $("#mini-lock")]) {
    b.setAttribute("aria-pressed", String(on));
    b.textContent = on ? "Release" : "Fix plane";
  }
  if (on) status("Plane fixed in the molecule: rotate to see it from the side; the depth slider moves it along its normal.");
  redraw();
}
renderer.debug.onShaderError = (gl, program, vs, fs) => {
  SLICE.failed = true;
  console.error("shader error", gl.getProgramInfoLog(program), gl.getShaderInfoLog(fs));
  status("The map slice is not available in this browser; back to mesh.");
  setTimeout(() => setStyle("mesh"), 0);
};
function setSliceDepth(z) {
  SLICE.depth = Math.max(-S.slab, Math.min(S.slab, z));
  for (const el of [$("#slice-depth"), $("#mini-depth-range")]) el.value = String(SLICE.depth);
  const txt = `${SLICE.depth >= 0 ? "+" : "−"}${Math.abs(SLICE.depth).toFixed(1)} Å`;
  $("#slice-depth-out").textContent = txt;
  $("#mini-depth-out").textContent = txt;
  redraw();
}

function setStyle(style) {
  if (style === "slice" && SLICE.failed) style = "mesh";
  S.style = style;
  $("#app").classList.toggle("slicing", style === "slice");
  document.querySelectorAll("#style-row button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.style === style)));
  try { localStorage.setItem("moorlight-share-style", style); } catch (err) { /* storage blocked */ }
  scheduleContour(true);
}

function renderMapList() {
  const ul = $("#map-list");
  ul.innerHTML = "";
  S.maps.forEach((m, k) => {
    const li = document.createElement("li");
    li.className = "map";
    const sw = m.info.kind === "diff"
      ? `<span class="sw" aria-hidden="true"><i style="background:var(--pos)"></i><i style="background:var(--neg)"></i></span>`
      : `<span class="sw" aria-hidden="true"><i style="background:var(--map)"></i></span>`;
    li.innerHTML = `<div class="map-head"><input type="checkbox" id="map-vis-${k}"><label for="map-vis-${k}">${sw}` +
      `<span class="name"></span></label></div>` +
      `<div class="lvl"><button type="button" aria-label="Lower contour level">&minus;</button>` +
      `<input type="range" id="map-lvl-${k}" min="0" max="1000" step="1">` +
      `<button type="button" aria-label="Raise contour level">+</button><output for="map-lvl-${k}"></output></div>`;
    li.querySelector(".name").textContent = m.info.name;
    const cb = li.querySelector("input[type=checkbox]");
    cb.checked = m.visible;
    cb.addEventListener("change", () => { m.visible = cb.checked; updateLevels(); scheduleContour(true); });
    const [down, up] = li.querySelectorAll("button");
    down.addEventListener("click", () => changeLevel(m, -levelStep(m)));
    up.addEventListener("click", () => changeLevel(m, levelStep(m)));
    const range = li.querySelector("input[type=range]");
    range.setAttribute("aria-label", `Contour level of ${m.info.name}`);
    bindSlider(range, () => m);
    m.range = range;
    m.out = li.querySelector("output");
    ul.appendChild(li);
  });
  bindSlider($("#mini-level-range"), activeMap);
  updateLevels();
}
const levelStep = (m) => (m.info.kind === "diff" ? 0.5 : 0.25);
const fmtLevel = (m, v = m.level) => `${m.kind === "diff" ? "±" : ""}${v < 10 ? v.toFixed(1) : v.toFixed(0)}${unit(m)}`;
const activeMap = () => S.maps.find((x) => x.visible) || null;

// Sliders run on a log scale from lo (1 for difference maps, 0.3 otherwise) to the map's largest value, so the
// usual 3-6 level range gets most of the travel and the strongest peaks are still in reach.
function levelRange(m) {
  const lo = m.info.kind === "diff" ? 1 : 0.3;
  return [lo, Math.max(lo * 4, m.maxAbs || lo * 20)];
}
function levelToPos(m, v) {
  const [lo, hi] = levelRange(m);
  return Math.round(1000 * Math.log(Math.min(hi, Math.max(lo, v)) / lo) / Math.log(hi / lo));
}
function posToLevel(m, p) {
  const [lo, hi] = levelRange(m);
  const v = lo * Math.pow(hi / lo, p / 1000);
  return v < 10 ? Math.round(v * 10) / 10 : Math.round(v);
}
function bindSlider(range, getMap) {
  range.addEventListener("input", () => {
    const m = getMap();
    if (!m) return;
    setLevel(m, posToLevel(m, +range.value), range);
  });
  range.addEventListener("change", () => scheduleContour(true));
}
function updateLevels(source = null) {
  for (const m of S.maps) {
    if (m.out) m.out.textContent = fmtLevel(m);
    if (m.range && m.range !== source) m.range.value = String(levelToPos(m, m.level));
  }
  S.maps.forEach((m, k) => { const cb = $(`#map-vis-${k}`); if (cb) cb.checked = m.visible; });
  const a = activeMap(), box = $("#mini-level");
  box.classList.toggle("none", !a);
  if (a) {
    $("#mini-level-name").textContent = a.info.name;
    $("#mini-level-out").textContent = fmtLevel(a);
    const r = $("#mini-level-range");
    r.setAttribute("aria-label", `Contour level of ${a.info.name}`);
    if (r !== source) r.value = String(levelToPos(a, a.level));
  }
}
function setLevel(m, v, source = null) {
  m.level = Math.max(0.1, +v.toFixed(2));
  if (SLICE.map === m && SLICE.mesh) { SLICE.mesh.material.uniforms.uLevel.value = m.level; redraw(); }
  if (!m.visible) m.visible = true;
  updateLevels(source);
  status(`${m.info.name}: contour at ${fmtLevel(m)}`);
  scheduleContour(!source);    // while a slider moves, recontour after a short pause
}
function changeLevel(m, step) { setLevel(m, m.level + step); }

// ------------------------------------------------------------ view and stops
let anim = null;
function animateTo(center, quat) {
  if (REDUCED) {
    S.center.copy(center);
    if (quat) rot.quaternion.copy(quat);
    redraw(); scheduleContour(true);
    return;
  }
  anim = { c0: S.center.clone(), c1: center.clone(), q0: rot.quaternion.clone(), q1: quat ? quat.clone() : null,
    t0: performance.now(), dur: 450 };
}
function animate() {
  if (!anim) return;
  const t = Math.min(1, (performance.now() - anim.t0) / anim.dur);
  const e = t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
  S.center.lerpVectors(anim.c0, anim.c1, e);
  if (anim.q1) rot.quaternion.slerpQuaternions(anim.q0, anim.q1, e);
  needsRender = true;
  if (t >= 1) { anim = null; scheduleContour(); }
}

function fitZoom(radius) {
  const half = Math.tan(THREE.MathUtils.degToRad(camera.fov / 2)) * Math.min(1, camera.aspect);
  return Math.min(400, Math.max(4, radius / half));
}

// the starting view of the page (also "Reset view")
function applyView() {
  const s = DATA.start;
  for (const m of S.maps) {
    const v = s.maps ? s.maps[m.info.name] : undefined;
    if (v === false || v === null) m.visible = false;
    else if (typeof v === "number") { m.visible = true; m.level = v; }
    else { m.visible = v === true ? true : m.info.visible; m.level = m.info.level; }
  }
  updateLevels();
  clearLabels();
  for (const [mk, i] of s.labels || []) if (S.mols[mk - 1]) toggleLabel(mk - 1, i);
  for (const [a, b] of s.distances || []) addDistance([a[0] - 1, a[1]], [b[0] - 1, b[1]]);
  S.radius = s.contour_radius || DATA.contour_radius || 10;
  for (const m of S.maps) m.drawn = null;
  S.zoom = fitZoom(s.radius || 10);
  S.slab = s.slab || s.radius || 10;
  SLICE.depth = 0;
  if (SLICE.lock) setSliceLock(false);
  setTilt("tiltX", 0); setTilt("tiltY", 0);
  syncDepthRange();
  const q = s.orient ? new THREE.Quaternion(...s.orient) : null;
  animateTo(new THREE.Vector3(...s.center), q);
  scheduleContour(true);
  status("");
}

// ---------------------------------------------------------------- structure
function renderModelList() {
  const ul = $("#mol-list");
  ul.innerHTML = "";
  S.mols.forEach((mol, k) => {
    const li = document.createElement("li");
    li.className = "mol";
    const chip = "#" + carbonColor(k).toString(16).padStart(6, "0");
    li.innerHTML = `<input type="checkbox" id="mol-vis-${k}"><label for="mol-vis-${k}">` +
      `<span class="chip" style="background:${chip}" aria-hidden="true"></span><span class="name"></span></label>`;
    li.querySelector(".name").textContent = mol.data.name;
    const cb = li.querySelector("input");
    cb.checked = mol.group.visible;
    cb.addEventListener("change", () => setMolVisible(k, cb.checked));
    ul.appendChild(li);
  });
  updateAtomsButton();
}
function setMolVisible(k, on) {
  S.mols[k].group.visible = on;
  const cb = $(`#mol-vis-${k}`);
  if (cb) cb.checked = on;
  updateAtomsButton();
  redraw();
}
function updateAtomsButton() {
  $("#mini-atoms").setAttribute("aria-pressed", String(S.mols.some((m) => m.group.visible)));
}
function toggleStructure() {
  const show = !S.mols.some((m) => m.group.visible);
  S.mols.forEach((m, k) => setMolVisible(k, show));
}

// ------------------------------------------------------------------ input
const el = renderer.domElement;
const pointers = new Map();
let gesture = null, lastTap = null;
el.addEventListener("contextmenu", (e) => e.preventDefault());
el.addEventListener("pointerdown", (e) => {
  e.preventDefault();
  try { el.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
  pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  if (pointers.size === 1) {
    gesture = { t0: performance.now(), x0: e.clientX, y0: e.clientY, moved: false, multi: false,
      button: e.button, pan: e.ctrlKey || e.metaKey || e.shiftKey, touch: e.pointerType !== "mouse" };
    anim = null;
  } else if (gesture) { gesture.multi = true; gesture.moved = true; }
});
el.addEventListener("pointermove", (e) => {
  const p = pointers.get(e.pointerId);
  if (!p || !gesture) return;
  const prev = { x: p.x, y: p.y };
  if (pointers.size >= 2) {
    const pts = [...pointers.entries()];
    const geom = (sub) => {
      const ps = pts.map(([id, q]) => (id === e.pointerId ? sub : q));
      const cx = ps.reduce((s, q) => s + q.x, 0) / ps.length, cy = ps.reduce((s, q) => s + q.y, 0) / ps.length;
      const d = Math.hypot(ps[0].x - ps[1].x, ps[0].y - ps[1].y);
      return { cx, cy, d };
    };
    const before = geom(prev), after = geom({ x: e.clientX, y: e.clientY });
    p.x = e.clientX; p.y = e.clientY;
    if (before.d > 0 && after.d > 0) S.zoom = Math.min(400, Math.max(4, S.zoom * before.d / after.d));
    pan(after.cx - before.cx, after.cy - before.cy);
    redraw();
    return;
  }
  p.x = e.clientX; p.y = e.clientY;
  if (!gesture.moved && Math.hypot(e.clientX - gesture.x0, e.clientY - gesture.y0) < (gesture.touch ? 8 : 3)) return;
  gesture.moved = true;
  const dx = e.clientX - prev.x, dy = e.clientY - prev.y;
  if (gesture.button === 2) {
    S.zoom = Math.min(400, Math.max(4, S.zoom * Math.exp(-dy * 0.01)));
    redraw();
  } else if (gesture.button === 1 || gesture.pan) {
    pan(dx, dy);
  } else {
    const axis = new THREE.Vector3(dy, dx, 0);
    const ang = axis.length() * 0.008;
    if (ang > 0) {
      rot.quaternion.premultiply(new THREE.Quaternion().setFromAxisAngle(axis.normalize(), ang));
      redraw();
    }
  }
});
function pan(dx, dy) {
  const h = 2 * S.zoom * Math.tan(THREE.MathUtils.degToRad(camera.fov / 2));
  const s = h / el.clientHeight;
  S.center.add(new THREE.Vector3(-dx * s, dy * s, 0).applyQuaternion(rot.quaternion.clone().invert()));
  redraw(); scheduleContour();
}
function endPointer(e) {
  if (!pointers.has(e.pointerId)) return;
  pointers.delete(e.pointerId);
  if (pointers.size || !gesture) return;
  const g = gesture;
  gesture = null;
  if (g.moved || g.multi || e.type === "pointercancel" || performance.now() - g.t0 > 600) return;
  const p = pick(e.clientX, e.clientY, g.touch ? 24 : 12);
  if (!p) return;
  const now = performance.now();
  if (lastTap && lastTap.mol === p.mol && lastTap.idx === p.idx && now - lastTap.t < 400) {
    toggleLabel(p.mol, p.idx);   // undo the first tap's label toggle
    const x = S.mols[p.mol].xyz;
    animateTo(new THREE.Vector3(x[3*p.idx], x[3*p.idx+1], x[3*p.idx+2]), null);
    lastTap = null;
    return;
  }
  lastTap = { ...p, t: now };
  toggleLabel(p.mol, p.idx);
  status(atomInfo(p.mol, p.idx));
}
el.addEventListener("pointerup", endPointer);
el.addEventListener("pointercancel", endPointer);
el.addEventListener("wheel", (e) => {
  e.preventDefault();
  if (e.ctrlKey) { S.zoom = Math.min(400, Math.max(4, S.zoom * Math.exp(e.deltaY * 0.01))); redraw(); return; }
  const m = S.maps.find((x) => x.visible);
  if (m) changeLevel(m, (e.shiftKey ? 0.5 : 0.1) * (e.deltaY > 0 ? -1 : 1));
}, { passive: false });
window.addEventListener("keydown", (e) => {
  if (e.target.closest && e.target.closest("input, select")) return;
  if (e.key === "f" || e.key === "F") setFull(!$("#app").classList.contains("full"));
  else if ((e.key === "l" || e.key === "L") && S.style === "slice") setSliceLock(!SLICE.lock);
  else if (e.key === "a" || e.key === "A") toggleStructure();
  else if ((e.key === "s" || e.key === "S") && S.style === "slice") savePlane();
});

document.querySelectorAll("#style-row button").forEach((b) => b.addEventListener("click", () => setStyle(b.dataset.style)));
try {
  const st = localStorage.getItem("moorlight-share-style");
  if (st === "surface" || st === "slice") S.style = st;
} catch (err) { /* storage blocked */ }
document.querySelectorAll("#style-row button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.style === S.style)));
$("#app").classList.toggle("slicing", S.style === "slice");
for (const el of [$("#slice-depth"), $("#mini-depth-range")]) el.addEventListener("input", () => setSliceDepth(+el.value));
for (const b of [$("#slice-lock"), $("#mini-lock")]) b.addEventListener("click", () => setSliceLock(!SLICE.lock));
$("#slice-tilt-x").addEventListener("input", (e) => setTilt("tiltX", +e.target.value));
$("#slice-tilt-y").addEventListener("input", (e) => setTilt("tiltY", +e.target.value));
$("#slice-face").addEventListener("click", faceMe);
for (const b of [$("#slice-live"), $("#mini-live")]) b.addEventListener("click", () => setLivePlane(SLICE.off));
$("#slice-thick").addEventListener("input", (e) => setThickness(+e.target.value));
$("#copy-offset").addEventListener("input", (e) => {
  const v = +e.target.value;
  $("#copy-offset-out").textContent = `${v >= 0 ? "+" : "\u2212"}${Math.abs(v).toFixed(1)} Å`;
});
for (const b of [$("#slice-save"), $("#mini-save")]) b.addEventListener("click", savePlane);
$("#slice-atoms").addEventListener("click", planeThroughLabels);
$("#mini-atoms").addEventListener("click", toggleStructure);
$("#reset").addEventListener("click", applyView);
$("#clear").addEventListener("click", () => { clearLabels(); status(""); });
$("#slab-down").addEventListener("click", () => { S.slab = Math.max(1, S.slab / 1.3); syncDepthRange(); redraw(); });
$("#slab-up").addEventListener("click", () => { S.slab = Math.min(200, S.slab * 1.3); syncDepthRange(); redraw(); });
function syncDepthRange() {
  for (const el of [$("#slice-depth"), $("#mini-depth-range")]) { el.min = String(-S.slab); el.max = String(S.slab); }
  setSliceDepth(SLICE.depth);
}
$("#hyd").addEventListener("click", (e) => {
  const show = S.mols.some((m) => m.hideH);
  S.mols.forEach((m, k) => { m.hideH = !show; buildMolecule(k); });
  e.currentTarget.setAttribute("aria-pressed", String(show));
});
const svgIcon = (d) => `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" ` +
  `stroke-linejoin="round" aria-hidden="true"><path d="${d}"/></svg>`;
const ICON_FULL = svgIcon("M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5");
const ICON_MENU = svgIcon("M4 6h16M4 12h16M4 18h16");
$("#fullscreen").innerHTML = ICON_FULL;
function setFull(on) {
  $("#app").classList.toggle("full", on);
  const b = $("#fullscreen");
  b.setAttribute("aria-pressed", String(on));
  b.setAttribute("aria-label", on ? "Show menus" : "Full screen view");
  b.innerHTML = on ? ICON_MENU : ICON_FULL;
  try { localStorage.setItem("moorlight-share-full", on ? "1" : "0"); } catch (err) { /* storage blocked */ }
}
$("#fullscreen").addEventListener("click", () => setFull(!$("#app").classList.contains("full")));
try { if (localStorage.getItem("moorlight-share-full") === "1") setFull(true); } catch (err) { /* storage blocked */ }

// ------------------------------------------------------------------- start
async function start() {
  DATA.models.forEach((d, k) => {
    const mol = { data: d, group: new THREE.Group(), xyz: Float32Array.from(d.xyz), hideH: !d.hydrogens };
    world.add(mol.group);
    S.mols.push(mol);
    buildMolecule(k);
  });
  $("#hyd").setAttribute("aria-pressed", String(S.mols.some((m) => !m.hideH)));
  resize();
  render();
  if (DATA.maps.length && typeof DecompressionStream === "undefined") {
    $("#loading").textContent = "This browser cannot unpack the maps. Use a current Safari, Chrome or Firefox.";
  } else {
    for (const info of DATA.maps) {
      const m = { info, values: await decodeMap(info), visible: info.visible, level: info.level, kind: info.kind,
        obj: new THREE.Group(), mats: [], smats: [], drawn: null, vol: null, tex: null };
      let mx = 0;
      for (const v of m.values) { const a = Math.abs(toLevelUnits(m, v)); if (a > mx) mx = a; }
      m.maxAbs = mx;
      world.add(m.obj);
      S.maps.push(m);
    }
    $("#loading").hidden = true;
  }
  renderMapList();
  renderModelList();
  loadPresetPlanes();
  loadPlanes();
  renderPlaneList();
  applyView();
}
start().catch((err) => {
  $("#loading").hidden = false;
  $("#loading").textContent = `Could not show the snapshot: ${err.message}`;
  console.error(err);
});
