// Moorlight share page: a read-only viewer for a snapshot written by `moorlight share`.
// Models, map boxes and tour stops are embedded in the page; maps are contoured here (surface nets),
// so no server is needed. Drawing follows the main viewer (app.js): Coot colours, sticks, chicken wire.
import * as THREE from "three";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";

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
  maps: [],            // {info, values: Float32Array, visible, level, obj, lastCenter, lastLevel}
  center: new THREE.Vector3(),
  zoom: 30, slab: 9,
  stop: -1,
  labels: new Map(),   // "mol:idx" -> div
  distances: [],
  radius: DATA.contour_radius || 10,
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
    div.style.display = v.z > -1 && v.z < 1 ? "" : "none";
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

// Surface nets: one vertex per grid cell the contour passes through (mean of its edge crossings), joined to the
// vertices of face-neighbours whose shared face the contour crosses. Drawn as lines this is Coot-like chicken wire.
const CORNERS = [...Array(8).keys()].map((c) => [c & 1, (c >> 1) & 1, (c >> 2) & 1]);
const CUBE_EDGES = [];
for (let a = 0; a < 8; a++) for (const bit of [1, 2, 4]) if (!(a & bit)) CUBE_EDGES.push([a, a | bit]);
const FACE = [0xaa, 0xcc, 0xf0];   // corners on the +x, +y, +z face of a cell

function surfaceNet(m, c, r, level, sign) {
  const { n, origin: o, G, Ginv, wrap } = m.info;
  const D = m.values;
  const off = m.info.units === "sigma" ? m.info.mean : 0, scale = m.info.units === "sigma" ? m.info.rms : 1;
  const L = level * scale;
  // grid cells covering the cube c +/- r (a parallelepiped on a non-orthogonal grid); vertices outside the cube are
  // dropped below, so the contour fills a box around the centre as in Coot
  const ci = mat3(Ginv, [c[0] - o[0], c[1] - o[1], c[2] - o[2]]);
  const lo = [], hi = [];
  for (let a = 0; a < 3; a++) {
    const ext = r * (Math.abs(Ginv[3*a]) + Math.abs(Ginv[3*a+1]) + Math.abs(Ginv[3*a+2]));
    let l = Math.floor(ci[a] - ext), h = Math.ceil(ci[a] + ext);
    if (!wrap) { l = Math.max(0, l); h = Math.min(n[a] - 1, h); }
    lo.push(l); hi.push(h);
  }
  const sx = hi[0] - lo[0], sy = hi[1] - lo[1], sz = hi[2] - lo[2];
  if (sx <= 0 || sy <= 0 || sz <= 0) return new Float32Array(0);
  // offset in the value array of each grid plane along each axis (wrapped on the periodic cell grid)
  const stride = [n[1] * n[2], n[2], 1];
  const [X, Y, Z] = [0, 1, 2].map((a) => {
    const t = new Int32Array(hi[a] - lo[a] + 1);
    for (let q = 0; q < t.length; q++) t[q] = (wrap ? wrapIndex(lo[a] + q, n[a]) : lo[a] + q) * stride[a];
    return t;
  });
  const CX = CORNERS.map((q) => q[0]), CY = CORNERS.map((q) => q[1]), CZ = CORNERS.map((q) => q[2]);
  const vid = new Int32Array(sx * sy * sz).fill(-1);
  const masks = new Uint8Array(sx * sy * sz);
  const verts = [];
  const val = new Float64Array(8);
  for (let a = 0; a < sx; a++) for (let b = 0; b < sy; b++) {
    let cell = (a * sy + b) * sz;
    for (let d = 0; d < sz; d++, cell++) {
      let mask = 0;
      for (let q = 0; q < 8; q++) {
        const v = sign * (D[X[a + CX[q]] + Y[b + CY[q]] + Z[d + CZ[q]]] - off);
        val[q] = v;
        if (v > L) mask |= 1 << q;
      }
      if (mask === 0 || mask === 255) continue;
      let ax = 0, ay = 0, az = 0, cnt = 0;
      for (const [p, q] of CUBE_EDGES) {
        if (((mask >> p) & 1) === ((mask >> q) & 1)) continue;
        const t = (L - val[p]) / (val[q] - val[p]);
        ax += CX[p] + t * (CX[q] - CX[p]); ay += CY[p] + t * (CY[q] - CY[p]); az += CZ[p] + t * (CZ[q] - CZ[p]);
        cnt++;
      }
      const fi = lo[0] + a + ax / cnt, fj = lo[1] + b + ay / cnt, fk = lo[2] + d + az / cnt;
      const px = o[0] + G[0]*fi + G[1]*fj + G[2]*fk;
      const py = o[1] + G[3]*fi + G[4]*fj + G[5]*fk;
      const pz = o[2] + G[6]*fi + G[7]*fj + G[8]*fk;
      if (Math.abs(px - c[0]) > r || Math.abs(py - c[1]) > r || Math.abs(pz - c[2]) > r) continue;
      vid[cell] = verts.length / 3;
      masks[cell] = mask;
      verts.push(px, py, pz);
    }
  }
  const seg = [];
  const strides = [sy * sz, sz, 1], lim = [sx, sy, sz];
  for (let a = 0; a < sx; a++) for (let b = 0; b < sy; b++) for (let d = 0; d < sz; d++) {
    const cell = (a * sy + b) * sz + d;
    const v = vid[cell];
    if (v < 0) continue;
    const idx = [a, b, d];
    for (let ax = 0; ax < 3; ax++) {
      if (idx[ax] + 1 >= lim[ax]) continue;
      const fm = masks[cell] & FACE[ax];
      if (fm === 0 || fm === FACE[ax]) continue;
      const w = vid[cell + strides[ax]];
      if (w < 0) continue;
      seg.push(verts[3*v], verts[3*v+1], verts[3*v+2], verts[3*w], verts[3*w+1], verts[3*w+2]);
    }
  }
  return new Float32Array(seg);
}

function clearMapObj(m) {
  for (const ch of [...m.obj.children]) { m.obj.remove(ch); ch.geometry.dispose(); }
}
function contourMap(m) {
  const c = S.center.toArray();
  clearMapObj(m);
  const signs = m.info.kind === "diff" ? [1, -1] : [1];
  for (const sign of signs) {
    const seg = surfaceNet(m, c, S.radius, m.level, sign);
    if (!seg.length) continue;
    const g = new LineSegmentsGeometry();
    g.setPositions(seg);
    const color = m.info.kind === "diff" ? (sign > 0 ? CONTOUR.pos : CONTOUR.neg) : (m.info.color || CONTOUR.map);
    m.obj.add(new LineSegments2(g, m.mats[sign > 0 ? 0 : 1] || (m.mats[sign > 0 ? 0 : 1] = lineMat(1.1, { color }))));
  }
  m.lastCenter = S.center.clone(); m.lastLevel = m.level;
  redraw();
}
let contourTimer = null;
function scheduleContour(now = false) {
  clearTimeout(contourTimer);
  contourTimer = setTimeout(() => {
    for (const m of S.maps) {
      m.obj.visible = m.visible;
      if (!m.visible) continue;
      if (m.lastCenter && m.lastCenter.distanceTo(S.center) < 0.5 && m.lastLevel === m.level) continue;
      contourMap(m);
    }
    redraw();
  }, now ? 0 : 120);
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
    li.innerHTML = `<input type="checkbox" id="map-vis-${k}"><label for="map-vis-${k}">${sw}<span class="name"></span></label>` +
      `<div class="lvl"><button type="button" aria-label="Lower contour level">&minus;</button><output></output>` +
      `<button type="button" aria-label="Raise contour level">+</button></div>`;
    li.querySelector(".name").textContent = m.info.name;
    const cb = li.querySelector("input");
    cb.checked = m.visible;
    cb.addEventListener("change", () => { m.visible = cb.checked; scheduleContour(true); });
    const [down, up] = li.querySelectorAll("button");
    down.addEventListener("click", () => changeLevel(m, -levelStep(m)));
    up.addEventListener("click", () => changeLevel(m, levelStep(m)));
    m.out = li.querySelector("output");
    ul.appendChild(li);
  });
  updateLevels();
}
const levelStep = (m) => (m.info.kind === "diff" ? 0.5 : 0.25);
function updateLevels() {
  for (const m of S.maps) if (m.out) m.out.textContent = `${m.kind === "diff" ? "±" : ""}${m.level.toFixed(1)}${unit(m)}`;
  S.maps.forEach((m, k) => { const cb = $(`#map-vis-${k}`); if (cb) cb.checked = m.visible; });
}
function changeLevel(m, step) {
  m.level = Math.max(0.1, +(m.level + step).toFixed(2));
  if (!m.visible) m.visible = true;
  updateLevels();
  status(`${m.info.name}: contour at ${m.info.kind === "diff" ? "±" : ""}${m.level.toFixed(2)}${unit(m)}`);
  scheduleContour(true);
}

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

function applyStop(k) {
  const stops = DATA.stops;
  S.stop = stops.length ? Math.max(0, Math.min(stops.length - 1, k)) : -1;
  const s = S.stop >= 0 ? stops[S.stop] : DATA.start;
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
  for (const m of S.maps) m.lastCenter = null;
  S.zoom = fitZoom(s.radius || 10);
  S.slab = s.slab || s.radius || 10;
  const q = s.orient ? new THREE.Quaternion(...s.orient) : null;
  animateTo(new THREE.Vector3(...s.center), q);
  scheduleContour(true);
  renderStop(s);
}

function renderStop(s) {
  const n = DATA.stops.length;
  $("#stop-pos").textContent = n ? `${S.stop + 1} / ${n}` : "";
  $("#stop-count").textContent = n ? `Stop ${S.stop + 1} of ${n}` : "View";
  $("#stop-name").textContent = s.name || "";
  const box = $("#stop-text");
  box.innerHTML = "";
  for (const para of String(s.text || "").split(/\n\s*\n/)) {
    if (!para.trim()) continue;
    const p = document.createElement("p");
    p.textContent = para.trim();
    box.appendChild(p);
  }
  $("#prev").disabled = !n || S.stop <= 0;
  $("#next").disabled = !n || S.stop >= n - 1;
  try { if (n) history.replaceState(null, "", `#stop-${S.stop + 1}`); } catch (e) { /* not allowed in some frames */ }
  status("");
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
  if (e.target.closest && e.target.closest("input")) return;
  if (e.key === "ArrowRight") applyStop(S.stop + 1);
  else if (e.key === "ArrowLeft") applyStop(S.stop - 1);
});

$("#prev").addEventListener("click", () => applyStop(S.stop - 1));
$("#next").addEventListener("click", () => applyStop(S.stop + 1));
$("#reset").addEventListener("click", () => applyStop(S.stop));
$("#clear").addEventListener("click", () => { clearLabels(); status(""); });
$("#slab-down").addEventListener("click", () => { S.slab = Math.max(1, S.slab / 1.3); redraw(); });
$("#slab-up").addEventListener("click", () => { S.slab = Math.min(200, S.slab * 1.3); redraw(); });
$("#hyd").addEventListener("click", (e) => {
  const show = S.mols.some((m) => m.hideH);
  S.mols.forEach((m, k) => { m.hideH = !show; buildMolecule(k); });
  e.currentTarget.setAttribute("aria-pressed", String(show));
});
$("#collapse").addEventListener("click", (e) => {
  const panel = $("#panel");
  const collapsed = panel.classList.toggle("collapsed");
  e.currentTarget.textContent = collapsed ? "Show panel" : "Hide panel";
  e.currentTarget.setAttribute("aria-expanded", String(!collapsed));
});

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
        obj: new THREE.Group(), mats: [], lastCenter: null, lastLevel: null };
      world.add(m.obj);
      S.maps.push(m);
    }
    $("#loading").hidden = true;
  }
  renderMapList();
  const h = /^#stop-(\d+)$/.exec(location.hash || "");
  applyStop(h ? parseInt(h[1], 10) - 1 : (DATA.stops.length ? 0 : -1));
}
start().catch((err) => {
  $("#loading").hidden = false;
  $("#loading").textContent = `Could not show the snapshot: ${err.message}`;
  console.error(err);
});
