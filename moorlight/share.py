"""Shareable snapshot: one self-contained HTML viewer page with models, maps and preset slice planes.

No server is needed to view it. The page contours the maps itself (three.js from a CDN), so it can be sent as a
file or published as a claude.ai artifact and opened on a phone.

    moorlight share spec.json -o view.html        # build the page
    moorlight share spec.json --check             # numbers for writing stop texts: peaks, values at atoms
    moorlight share spec.json --check --at "-/102 N9" --near "-/102" --radius 4

The spec (JSON; paths relative to the spec file):

    {"title": "Cage 3.5 ps",                      # page name, two to four words
     "description": "one line shown under the title",
     "models": [{"path": "dark.pdb", "name": "dark"}],
     "maps": [{"path": "3p5ps.ccp4", "name": "3.5 ps light - dark", "kind": "diff",
               "units": "absolute", "unit_label": "SD", "level": 4.0, "visible": true}],
     "view": {"center": "-/102", "radius": 8, "orient": "face", "labels": ["-/102 N9"]},
     "planes": [{"name": "Azo 102", "atoms": ["-/102 C8", "-/102 N9", "-/102 N10", "-/102 C11"],
                 "thickness": 0.5, "map": "3.5 ps light - dark"}]}

See skills/moorlight-share/SKILL.md for every field.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import html
import itertools
import json
import os
import re
import sys

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .maps import SymMapper, detect_columns, check_columns, guess_kind
from .model import Molecule
from .session import parse_atom_spec

WEB = os.path.join(os.path.dirname(__file__), "web")
MAX_BYTES = 15 * 1024 * 1024    # artifact pages must stay at or below 16 MB
NATIVE_LIMIT = 4_000_000        # maps up to this many grid points are kept on their own grid
NEIGHBOUR_COLOR = "#8a9099"
_TAGGED = re.compile(r"^(.*?)\s*\[([^\]]*)\]\s*$")


# --------------------------------------------------------- neighbour copies
def _shift_tag(t):
    parts = [("+" if v > 0 else "-") + (str(abs(v)) if abs(v) > 1 else "") + "abc"[a] for a, v in enumerate(t) if v]
    return "".join(parts)


class Neighbours:
    """Symmetry and lattice copies of a model's residues that come within `margin` Å of the model (the model itself
    excluded), as Coot draws symmetry atoms. Residues are copied whole, so the copies keep their bonds."""

    def __init__(self, mol: Molecule, margin=5.0):
        st = mol.st
        cell = st.cell
        sg = st.find_spacegroup() or gemmi.SpaceGroup("P 1")
        M = np.array(cell.orth.mat.tolist())
        Minv = np.array(cell.frac.mat.tolist())
        heavy = np.array([e != "H" for e in mol.elements])
        F = mol.xyz @ Minv.T
        groups = {}
        for i, inf in enumerate(mol.info):
            if heavy[i]:
                groups.setdefault(inf[:4], []).append(i)
        groups = [np.array(g) for g in groups.values()]
        tree = cKDTree(mol.xyz[heavy])
        lo, hi = mol.xyz[heavy].min(0) - margin, mol.xyz[heavy].max(0) + margin
        corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]) @ Minv.T
        fmin, fmax = corners.min(0), corners.max(0)
        self.orig, self.tags, xyz = [], [], []
        for op in sg.operations():
            R, T = np.array(op.rot) / op.DEN, np.array(op.tran) / op.DEN
            Fo = F @ R.T + T
            ident = op.triplet() == "x,y,z"
            tmin = np.floor(fmin - Fo.max(0)).astype(int)
            tmax = np.ceil(fmax - Fo.min(0)).astype(int)
            for t in itertools.product(*(range(a, b + 1) for a, b in zip(tmin, tmax))):
                if ident and not any(t):
                    continue
                X = (Fo + np.array(t)) @ M.T
                tag = ("" if ident else op.triplet() + " ") + _shift_tag(t)
                for g in groups:
                    if np.isfinite(tree.query(X[g], distance_upper_bound=margin)[0]).any():
                        self.orig.extend(g.tolist())
                        self.tags.extend([tag] * len(g))
                        xyz.append(X[g])
        self.xyz = np.concatenate(xyz) if xyz else np.zeros((0, 3))
        self.index = {(o, t): k for k, (o, t) in enumerate(zip(self.orig, self.tags))}
        # bonds of the model between atoms of the same copy
        pos = {}
        for k, (o, t) in enumerate(zip(self.orig, self.tags)):
            pos[(o, t)] = k
        bonds, orders = [], []
        for (i, j), order in zip(mol.bonds.tolist(), mol.orders.tolist()):
            for t in set(self.tags):
                a, b = pos.get((i, t)), pos.get((j, t))
                if a is not None and b is not None:
                    bonds.append((a, b))
                    orders.append(order)
        self.bonds, self.orders = bonds, orders
        self.mol = mol

    def to_json(self, number, name):
        m = self.mol
        atoms = []
        for o, t in zip(self.orig, self.tags):
            i = m.info[o]
            atoms.append(list(i[:6]) + [round(i[6], 2), round(i[7], 1), i[8], t])
        return {"id": number, "name": name, "xyz": np.round(self.xyz, 3).ravel().tolist(),
                "element": [m.elements[o] for o in self.orig], "atoms": atoms,
                "bonds": [v for b in self.bonds for v in b], "orders": self.orders, "visible": True,
                "color": NEIGHBOUR_COLOR, "mono": True, "hydrogens": False, "style": "lines", "neighbour": True}


# ------------------------------------------------------------------ atoms
class Atoms:
    """Atom-spec lookup over what the page draws.

    Molecule numbers follow the page: the spec's models (1, 2, ...), then their neighbour copies. A spec with a tag
    in brackets, 'S/14 O [+a]', names the copy of that atom one lattice vector a away."""

    def __init__(self, mols, neighbours):
        self.mols = mols
        self.nb = neighbours              # model number -> Neighbours (or absent)
        self.nb_number = {k: len(mols) + j + 1 for j, k in enumerate(sorted(neighbours))}

    def resolve(self, spec, first_conformer=True):
        """Atom spec -> (page molecule number, [atom indices], xyz array). A residue spec gives its heavy atoms."""
        tag = None
        if isinstance(spec, str):
            m_ = _TAGGED.match(spec)
            if m_:
                spec, tag = m_.group(1), m_.group(2).strip()
        d = parse_atom_spec(spec)
        k = d["mol"] or 1
        if not 1 <= k <= len(self.mols):
            raise ValueError(f"atom spec {spec!r}: no model {k} (there are {len(self.mols)})")
        m = self.mols[k - 1]
        res = [i for i, key in enumerate(m.keys) if key[0] == d["chain"] and key[1] == d["resseq"] and key[2] == d["icode"]]
        if not res:
            raise ValueError(f"atom spec {spec!r}: no residue {d['chain']}/{d['resseq']}{d['icode']} in {m.name}; "
                             f"residues: {self.residue_list(m)}")
        if d["atom"] and d["atom"] == m.info[res[0]][3] and not any(m.keys[i][3] == d["atom"] for i in res):
            d["atom"] = None                        # 'A/85 LYS': a residue name, not an atom
        idx = [i for i in res if (m.keys[i][3] == d["atom"] if d["atom"] else m.elements[i] != "H")
               and (d["altloc"] is None or m.keys[i][4] in (d["altloc"], ""))]
        if not idx:
            names = sorted({m.keys[i][3] for i in res})
            raise ValueError(f"atom spec {spec!r}: no such atom; atoms of that residue: {' '.join(names)}")
        if d["altloc"] is None and first_conformer:
            alts = sorted({m.keys[i][4] for i in idx} - {""})
            idx = [i for i in idx if m.keys[i][4] in ("", alts[0] if alts else "")]
        if not tag:
            return k, idx, m.xyz[idx]
        nb = self.nb.get(k)
        if nb is None:
            raise ValueError(f"atom spec {spec!r} [{tag}]: the page draws no neighbour copies of model {k}")
        out = [nb.index.get((i, tag)) for i in idx]
        if any(o is None for o in out):
            have = sorted({t for (o, t) in nb.index if o in idx})
            raise ValueError(f"atom spec {spec!r}: no copy [{tag}] on the page; copies drawn: "
                             + (", ".join(f"[{t}]" for t in have) or "none"))
        return self.nb_number[k], out, nb.xyz[out]

    def point(self, spec):
        """xyz list, an atom spec (residue: centroid) or a list of specs (centroid of all)."""
        if _is_xyz(spec):
            return np.asarray(spec, float)
        if isinstance(spec, (list, tuple)):
            return np.mean([self.point(s) for s in spec], axis=0)
        return self.resolve(spec)[2].mean(0)

    def xyz(self, specs):
        """All atom positions of a list of specs (residues give all their heavy atoms)."""
        out = [self.resolve(s)[2] for s in specs]
        return np.concatenate(out) if out else np.zeros((0, 3))

    def name(self, number, i):
        """Spec of an atom on the page, as the specs are written."""
        if number <= len(self.mols):
            a, tag, k = self.mols[number - 1].info[i], "", number
        else:
            k = next(kk for kk, nn in self.nb_number.items() if nn == number)
            nb = self.nb[k]
            a, tag = nb.mol.info[nb.orig[i]], nb.tags[i]
        spec = f"{a[0]}/{a[1]}{a[2]} {a[4]}" + (f":{a[5]}" if a[5] else "")
        if len(self.mols) > 1:
            spec = f"{k}#{spec}"
        return spec + (f" [{tag}]" if tag else "")

    def nearest(self):
        """f(xyz) -> (spec, distance) of the nearest heavy atom drawn on the page (models and neighbour copies)."""
        pts, ref = [], []
        for k, m in enumerate(self.mols, 1):
            for i, e in enumerate(m.elements):
                if e != "H":
                    pts.append(m.xyz[i])
                    ref.append((k, i))
        for k, nb in self.nb.items():
            for i in range(len(nb.orig)):
                pts.append(nb.xyz[i])
                ref.append((self.nb_number[k], i))
        tree = cKDTree(np.array(pts))

        def f(xyz):
            d, j = tree.query(xyz)
            return self.name(*ref[j]), float(d)
        return f

    def near_model(self):
        """f(xyz) -> the symmetry/lattice image of xyz next to the models (for listing each peak once)."""
        mappers = [SymMapper(m.st, m.xyz, np.array([e != "H" for e in m.elements])) for m in self.mols]

        def f(xyz):
            return np.asarray(min((sm.nearest(xyz) for sm in mappers), key=lambda r: r[1])[0])
        return f

    @staticmethod
    def residue_list(m, limit=40):
        seen = []
        for key, inf in zip(m.keys, m.info):
            r = f"{key[0]}/{key[1]}{key[2]} {inf[3]}"
            if r not in seen:
                seen.append(r)
        return ", ".join(seen[:limit]) + (" ..." if len(seen) > limit else "")


def _is_xyz(v):
    return isinstance(v, (list, tuple)) and len(v) == 3 and all(isinstance(x, (int, float)) for x in v)


# ------------------------------------------------------------------- maps
def read_map(spec, base):
    """Map spec -> (gemmi.FloatGrid over the unit cell, kind). CCP4/MRC maps, or MTZ with columns 'F,PHI'."""
    p = _path(spec["path"], base)
    low = p.lower()
    if low.endswith((".ccp4", ".map", ".mrc")):
        g = gemmi.read_ccp4_map(p, setup=True).grid
        return g, spec.get("kind", "diff")
    if low.endswith(".mtz"):
        mtz = gemmi.read_mtz_file(p)
        if spec.get("columns"):
            f, phi = [c.strip() for c in spec["columns"].split(",")[:2]]
            check_columns(mtz, f, phi)
            kind = guess_kind(f)
        else:
            found = detect_columns(mtz)
            if not found:
                raise ValueError(f"{p}: no map columns found; give \"columns\": \"F,PHI\"")
            f, phi, kind, _ = found[0]
        sample = float(spec.get("sample_rate", 3.0))
        return mtz.transform_f_phi_to_map(f, phi, sample_rate=sample), spec.get("kind", kind)
    raise ValueError(f"{p}: unknown map format (use .ccp4/.map/.mrc or .mtz)")


class PageMap:
    """A map as the page stores it: values on a grid with positions xyz = origin + G @ (i, j, k).

    native: the map's own grid over the unit cell, periodic, nothing resampled (values as in Coot).
    box: resampled on an orthogonal box (for cells too large to ship whole)."""

    def __init__(self, spec, base, box_lo, box_hi):
        grid, kind = read_map(spec, base)
        whole = np.array(grid, copy=False)
        self.name = spec.get("name") or os.path.basename(spec["path"])
        self.kind = kind if kind in ("diff", "2fofc") else "diff"
        self.mean, self.rms = float(whole.mean()), float(whole.std())
        self.units = spec.get("units", "sigma")
        if self.units not in ("sigma", "absolute"):
            raise ValueError(f"map {self.name!r}: units must be 'sigma' or 'absolute'")
        self.unit_label = spec.get("unit_label") or ("σ" if self.units == "sigma" else "")
        self.level = float(spec.get("level", 3.0 if self.kind == "diff" else 1.5))
        self.visible = bool(spec.get("visible", True))
        self.color = spec.get("color")
        self.source = os.path.basename(spec["path"])
        mode = spec.get("grid", "auto")
        if mode == "auto":
            mode = "native" if whole.size <= NATIVE_LIMIT else "box"
        self.mode = mode
        self.cell = np.array(grid.unit_cell.orth.mat.tolist())
        if mode == "native":
            self.arr = np.array(whole, dtype=np.float32)
            self.n = np.array(self.arr.shape)
            self.G = self.cell / self.n            # column a: one grid step along axis a
            self.origin = np.zeros(3)
            self.wrap = True
        elif mode == "box":
            uc = grid.unit_cell
            step = float(spec.get("step") or min(uc.a / grid.nu, uc.b / grid.nv, uc.c / grid.nw) / 2)
            lo, hi = np.asarray(box_lo, float), np.asarray(box_hi, float)
            self.n = np.maximum(2, np.ceil((hi - lo) / step).astype(int) + 1)
            self.arr = np.zeros(tuple(self.n), dtype=np.float32)
            tr = gemmi.Transform()
            tr.mat.fromlist([[step, 0, 0], [0, step, 0], [0, 0, step]])
            tr.vec.fromlist([float(v) for v in lo])
            grid.interpolate_values(self.arr, tr, order=3)   # tricubic: trilinear on a new grid flattens peaks
            self.G = np.eye(3) * step
            self.origin = lo
            self.wrap = False
        else:
            raise ValueError(f"map {self.name!r}: grid must be 'auto', 'native' or 'box'")
        self.Ginv = np.linalg.inv(self.G)
        self._peaks = {}

    def to_level_units(self, v):
        """Map value -> the page's level units (σ of the whole map, or the map's own units)."""
        return (v - self.mean) / self.rms if self.units == "sigma" else v

    def value_at(self, xyz):
        """Trilinear interpolation on the stored grid, in level units (the page does the same); NaN outside a box."""
        f = self.Ginv @ (np.asarray(xyz, float) - self.origin)
        if not self.wrap and (np.any(f < 0) or np.any(f > self.n - 1)):
            return float("nan")
        i0 = np.floor(f).astype(int) if self.wrap else np.minimum(np.floor(f).astype(int), self.n - 2)
        t = f - i0
        v = 0.0
        for c in range(8):
            d = np.array([c & 1, (c >> 1) & 1, (c >> 2) & 1])
            idx = i0 + d
            if self.wrap:
                idx = idx % self.n
            v += np.prod(np.where(d, t, 1 - t)) * self.arr[tuple(idx)]
        return float(self.to_level_units(v))

    def image_near(self, xyz, c):
        """The lattice image of xyz closest to c (the stored grid is periodic only in native mode)."""
        if not self.wrap:
            return np.asarray(xyz, float)
        f = np.linalg.solve(self.cell, np.asarray(xyz, float) - c)
        return c + self.cell @ (f - np.round(f))

    def peaks(self, level, negative=False):
        """Grid points that are the extreme of their 3x3x3 neighbourhood beyond `level` (level units).
        One finder for every list: [(xyz, value)], strongest first."""
        key = (round(level, 4), negative)
        if key not in self._peaks:
            from scipy.ndimage import maximum_filter
            a = self.to_level_units(self.arr)
            a = -a if negative else a
            mx = maximum_filter(a, size=3, mode="wrap" if self.wrap else "nearest")
            idx = np.argwhere((a == mx) & (a > level))
            val = a[tuple(idx.T)]
            order = np.argsort(-val)
            xyz = self.origin + idx[order] @ self.G.T
            self._peaks[key] = [(p, float(-v if negative else v)) for p, v in zip(xyz, val[order])]
        return self._peaks[key]

    def to_json(self):
        scale = max(float(np.abs(self.arr).max()), 1e-6) / 32000.0
        q = np.round(self.arr / scale).astype("<i2")
        data = base64.b64encode(gzip.compress(q.tobytes(), compresslevel=9, mtime=0)).decode()
        return {"name": self.name, "kind": self.kind, "units": self.units, "unit_label": self.unit_label,
                "level": self.level, "visible": self.visible, "color": self.color, "mean": self.mean, "rms": self.rms,
                "grid": self.mode, "n": [int(v) for v in self.n], "origin": [round(float(v), 4) for v in self.origin],
                "G": [round(float(v), 7) for v in self.G.ravel()], "Ginv": [round(float(v), 7) for v in self.Ginv.ravel()],
                "wrap": self.wrap, "scale": scale, "source": self.source, "data": data}


# ------------------------------------------------------------- orientation
def screen_axes(xyz, mode="face"):
    """Rows: screen x, y, z (towards the viewer) in model coordinates, from the atoms' principal axes.

    face: their plane faces the viewer, longest axis horizontal. edge: their plane is seen edge-on (its normal points
    up on screen), which shows out-of-plane motion."""
    c = xyz - xyz.mean(0)
    w, v = np.linalg.eigh(c.T @ c)          # ascending
    e1, e2, e3 = v[:, 2], v[:, 1], v[:, 0]  # longest, middle, shortest
    if mode == "edge":
        e2, e3 = e3, e2
    return np.array([e1, e2, np.cross(e1, e2)])


def _quat(R):
    """Quaternion [x, y, z, w] (three.js order) of a rotation matrix."""
    t = np.trace(R)
    if t > 0:
        s = 2 * np.sqrt(t + 1)
        q = [(R[2, 1] - R[1, 2]) / s, (R[0, 2] - R[2, 0]) / s, (R[1, 0] - R[0, 1]) / s, s / 4]
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2 * np.sqrt(1 + R[0, 0] - R[1, 1] - R[2, 2])
        q = [s / 4, (R[0, 1] + R[1, 0]) / s, (R[0, 2] + R[2, 0]) / s, (R[2, 1] - R[1, 2]) / s]
    elif R[1, 1] > R[2, 2]:
        s = 2 * np.sqrt(1 + R[1, 1] - R[0, 0] - R[2, 2])
        q = [(R[0, 1] + R[1, 0]) / s, s / 4, (R[1, 2] + R[2, 1]) / s, (R[0, 2] - R[2, 0]) / s]
    else:
        s = 2 * np.sqrt(1 + R[2, 2] - R[0, 0] - R[1, 1])
        q = [(R[0, 2] + R[2, 0]) / s, (R[1, 2] + R[2, 1]) / s, s / 4, (R[1, 0] - R[0, 1]) / s]
    q = np.array(q)
    return [round(float(v), 6) for v in q / np.linalg.norm(q)]


def _rot_from_quat(q):
    x, y, z, w = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


# ------------------------------------------------------------------ build
class Share:
    def __init__(self, spec, base):
        self.spec, self.base = spec, base
        self.mols = []
        for ms in spec.get("models") or []:
            m = Molecule(_path(ms["path"], base), name=ms.get("name"))
            m.color, m.hydrogens = ms.get("color"), ms.get("hydrogens", False)
            m.style = ms.get("style")
            m.visible = ms.get("visible", True)
            m.neighbours = ms.get("neighbours", True)
            self.mols.append(m)
        if not self.mols:
            raise ValueError("the spec needs at least one model")
        lo, hi = self.region()
        margin = float(spec.get("neighbour_margin", 5.0))
        neighbours = {}
        for k, m in enumerate(self.mols, 1):
            if m.neighbours and m.st.cell.is_crystal():
                nb = Neighbours(m, margin)
                if len(nb.orig):
                    neighbours[k] = nb
        self.atoms = Atoms(self.mols, neighbours)
        self.maps = [PageMap(s, base, lo, hi) for s in spec.get("maps") or []]
        names = [x.name for x in self.maps]
        if len(set(names)) != len(names):
            raise ValueError(f"map names must be unique: {names}")
        # the page opens on one view; older specs with "stops" open on their first stop
        first = spec.get("view") or (spec.get("stops") or [None])[0]
        self.view = self.stop(first, 0) if first else None
        self.stops = [self.view] if self.view else []
        self.planes = []
        for k, p in enumerate(spec.get("planes") or []):
            self.planes.append(self.plane(p, k))
            for o in p.get("copies") or []:      # parallel copies at these offsets (Å) along the normal
                name = f"{self.planes[-1]['name']} {o:+g} \u00c5"
                self.planes.append(self.plane({**p, "name": name, "offset": float(p.get("offset", 0)) + float(o),
                                               "color": None}, k))

    def region(self):
        """Box for maps that are resampled (grid 'box'): all models plus a margin, or a cube around a centre."""
        r = self.spec.get("region") or {}
        if "center" in r:
            c, rad = self.point(r["center"]), float(r.get("radius", 12.0))
            return c - rad, c + rad
        xyz = np.concatenate([m.xyz for m in self.mols])
        margin = float(r.get("margin", 5.0))
        return xyz.min(0) - margin, xyz.max(0) + margin

    def point(self, spec):
        if hasattr(self, "atoms"):
            return self.atoms.point(spec)
        return Atoms(self.mols, {}).point(spec)

    def stop(self, s, k):
        where = f"stop {k + 1} ({s.get('name', '')!r})"
        try:
            center = self.atoms.point(s["center"])
            frame = s.get("frame")
            if frame is None and not _is_xyz(s["center"]):
                frame = s["center"] if isinstance(s["center"], (list, tuple)) else [s["center"]]
            mode = s.get("orient", "face")
            axes, orient, plane = None, None, None
            if isinstance(mode, (list, tuple)):
                orient = [float(v) for v in mode]
                axes = _rot_from_quat(orient)
            elif mode != "keep" and frame:
                pts = self.atoms.xyz(frame)
                if len(pts) < 3 and s.get("frame") is None and len(pts) == 1:
                    # a lone atom (a metal): use it and the atoms within 2.6 A, its coordination
                    near_pts = [xyz for _, _, xyz in self._drawn() if 0.1 < np.linalg.norm(xyz - pts[0]) <= 2.6]
                    if len(near_pts) >= 2:
                        pts = np.vstack([pts] + near_pts)
                if len(pts) >= 3:
                    axes = screen_axes(pts, mode)
                    orient = _quat(axes)
                    plane = (pts.mean(0), axes[2] if mode == "face" else axes[1],
                             "towards the viewer" if mode == "face" else "up on screen")
            maps = {}
            for name, v in (s.get("maps") or {}).items():
                if name not in {x.name for x in self.maps}:
                    raise ValueError(f"no map named {name!r}; maps: {[x.name for x in self.maps]}")
                maps[name] = v
            labels = []
            for spec in s.get("labels") or []:
                number, idx, _ = self.atoms.resolve(spec)
                labels.extend([number, i] for i in idx)
            dists = []
            for a, b in s.get("distances") or []:
                na, ia, _ = self.atoms.resolve(a)
                nb, ib, _ = self.atoms.resolve(b)
                dists.append([[na, ia[0]], [nb, ib[0]]])
        except (ValueError, KeyError) as e:
            raise ValueError(f"{where}: {e}") from None
        radius = float(s.get("radius", 8.0))
        out = {"name": s.get("name") or f"stop {k + 1}", "text": s.get("text", ""),
               "center": [round(float(v), 3) for v in center], "orient": orient, "radius": radius,
               "slab": float(s.get("slab", radius)), "maps": maps, "labels": labels, "distances": dists,
               "_center_spec": s["center"], "_plane": plane}
        if "contour_radius" in s:
            out["contour_radius"] = float(s["contour_radius"])
        return out

    def _drawn(self):
        """(page molecule number, index, xyz) of every heavy atom drawn on the page."""
        for k, m in enumerate(self.mols, 1):
            for i, e in enumerate(m.elements):
                if e != "H":
                    yield k, i, m.xyz[i]
        for k, nb in self.atoms.nb.items():
            for i in range(len(nb.orig)):
                yield self.atoms.nb_number[k], i, nb.xyz[i]

    def plane(self, p, k):
        """Preset slice plane: the best plane through `atoms` (three or more), or `center` with a `normal` (a vector,
        or two atom specs from the first to the second). `axis` (two atom specs) sets its horizontal direction;
        `offset` moves it along the normal; `thickness`, `radius`, `map`, `level`, `color`, `visible`."""
        where = f"plane {k + 1} ({p.get('name', '')!r})"
        try:
            rms = None
            if p.get("atoms"):
                pts = self.atoms.xyz(p["atoms"])
                if len(pts) < 3:
                    raise ValueError("\"atoms\" needs three or more atoms")
                ex, ey, n = screen_axes(pts, "face")
                centre = self.atoms.point(p["center"]) if "center" in p else pts.mean(0)
                rms = float(np.sqrt(np.mean(((pts - pts.mean(0)) @ n) ** 2)))
            elif "center" in p and "normal" in p:
                centre = self.atoms.point(p["center"])
                nv = p["normal"]
                if _is_xyz(nv):
                    n = np.asarray(nv, float)
                elif isinstance(nv, (list, tuple)) and len(nv) == 2:
                    n = self.atoms.point(nv[1]) - self.atoms.point(nv[0])
                else:
                    raise ValueError("\"normal\" is a vector [x, y, z] or two atom specs")
                n = n / np.linalg.norm(n)
                ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
                ex = np.cross(ref, n)
                ex /= np.linalg.norm(ex)
                ey = np.cross(n, ex)
            else:
                raise ValueError("give \"atoms\" (three or more) or \"center\" and \"normal\"")
            if p.get("axis"):
                a, b = p["axis"]
                d = self.atoms.point(b) - self.atoms.point(a)
                d -= n * (d @ n)
                if np.linalg.norm(d) > 1e-6:
                    ex = d / np.linalg.norm(d)
                    ey = np.cross(n, ex)
            name = p.get("map")
            if name is not None and name not in {x.name for x in self.maps}:
                raise ValueError(f"no map named {name!r}; maps: {[x.name for x in self.maps]}")
        except (ValueError, KeyError) as e:
            raise ValueError(f"{where}: {e}") from None
        p0 = centre + n * float(p.get("offset", 0.0))
        color = p.get("color")
        if isinstance(color, str):
            color = int(color.lstrip("#"), 16)
        r3 = lambda v: [round(float(t), 4) for t in v]
        return {"name": p.get("name") or f"Preset {k + 1}", "p0": r3(p0), "ex": r3(ex), "ey": r3(ey), "n": r3(n),
                "r": float(p.get("radius", self.spec.get("contour_radius", 10.0))),
                "thick": float(p.get("thickness", 0.0)), "map": name, "level": p.get("level"), "color": color,
                "visible": bool(p.get("visible", True)), "_rms": rms}

    def payload(self):
        mols = []
        for k, m in enumerate(self.mols, 1):
            j = m.to_json()
            mols.append({"id": k, "name": m.name, "xyz": j["xyz"], "element": j["element"], "atoms": j["atoms"],
                         "bonds": j["bonds"], "orders": j["orders"], "visible": m.visible, "color": m.color,
                         "hydrogens": bool(m.hydrogens), "style": m.style or "sticks"})
        for k, nb in sorted(self.atoms.nb.items()):
            mols.append(nb.to_json(self.atoms.nb_number[k], f"{self.mols[k - 1].name}, neighbouring cells"))
        planes = [{k: v for k, v in pl.items() if not k.startswith("_")} for pl in self.planes]
        if self.view:
            start = {k: v for k, v in self.view.items() if not k.startswith("_")}
        else:
            c = np.concatenate([m.xyz for m in self.mols]).mean(0)
            start = {"name": "", "text": "", "center": [round(float(v), 3) for v in c], "orient": None,
                     "radius": 12.0, "slab": 12.0, "maps": {}, "labels": [], "distances": []}
        return {"title": self.spec.get("title", "Moorlight"), "description": self.spec.get("description", ""),
                "contour_radius": float(self.spec.get("contour_radius", 10.0)), "models": mols,
                "maps": [x.to_json() for x in self.maps], "planes": planes, "start": start}

    def html(self, standalone=False):
        data = json.dumps(self.payload(), separators=(",", ":")).replace("</", "<\\/")
        with open(os.path.join(WEB, "share.html")) as fh:
            page = fh.read()
        with open(os.path.join(WEB, "share.js")) as fh:
            js = fh.read()
        page = (page.replace("{{TITLE}}", html.escape(self.spec.get("title", "Moorlight")))
                .replace("{{DESCRIPTION}}", html.escape(self.spec.get("description", "")))
                .replace("{{SCRIPT}}", js)
                .replace("{{DATA}}", data))
        if standalone:
            page = ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">\n"
                    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
                    + page + "\n</html>\n")
        return page

    # ------------------------------------------------------------- check
    def _peak_lines(self, x, center, radius, level, near, plane=None, limit=8):
        """Peaks of one map within radius of center (their images next to it), as text lines."""
        u = x.unit_label or x.units
        lines = []
        for neg in (False, True):
            found = []
            for xyz, v in x.peaks(level, negative=neg):
                p = x.image_near(xyz, center)
                if np.linalg.norm(p - center) <= radius:
                    found.append((p, v))
            head = f"    {'+' if not neg else '-'} beyond {'-' if neg else '+'}{level:g} {u}: {len(found)}"
            lines.append(head + (", strongest:" if found else ""))
            for p, v in found[:limit]:
                nm, d = near(p)
                h = ""
                if plane is not None:
                    h = f"  h {float((p - plane[0]) @ plane[1]):+.2f} A"
                lines.append(f"      {v:+7.2f} {u} at ({p[0]:.2f}, {p[1]:.2f}, {p[2]:.2f})  {d:.2f} A from {nm}{h}")
        return lines

    def report(self, at=(), near_specs=(), near_radius=5.0, n_peaks=12, stop_peaks=8, maps_only=(), frame=None):
        """Text for writing stop texts. maps_only: names of the maps to list (default all); frame: atom specs whose
        plane gives h for the --near lists."""
        out = []
        p = out.append
        near = self.atoms.nearest()
        to_model = self.atoms.near_model()
        p("Peaks are grid points that are the extreme of their 3x3x3 neighbourhood; one finder serves every list. "
          "Values at atoms are trilinear on the stored grid, as the page shows them on a tap. Atom names with a tag in "
          "brackets ([+a], [-b+c], ...) are neighbour copies drawn on the page; use the same text in specs.")
        p("\nmodels: " + "; ".join(f"{k}# {m.name}: {len(m.keys)} atoms ({Atoms.residue_list(m, 12)})"
                                  for k, m in enumerate(self.mols, 1)))
        for k, nb in sorted(self.atoms.nb.items()):
            p(f"neighbour copies of model {k} (page molecule {self.atoms.nb_number[k]}): {len(nb.orig)} heavy atoms in residues "
              f"within {self.spec.get('neighbour_margin', 5.0)} A of the model; copies: "
              + ", ".join(f"[{t}]" for t in sorted(set(nb.tags))))
        for x in self.maps:
            if maps_only and x.name not in maps_only:
                continue
            u = x.unit_label or x.units
            vals = x.to_level_units(x.arr)
            grid = (f"own grid {x.n.tolist()} over the cell, nothing resampled" if x.mode == "native"
                    else f"resampled on a box {x.n.tolist()} x {x.G[0, 0]:.2f} A (tricubic)")
            p(f"\nmap {x.name!r} ({x.source}, {x.kind}, {grid}): level {x.level:g} {u}; whole map mean {x.mean:.3f} "
              f"rms {x.rms:.3f}; stored values {vals.min():.2f} to {vals.max():.2f} {u}")
            for neg in ([False, True] if x.kind == "diff" else [False]):
                pk = x.peaks(x.level, negative=neg)
                listed = []          # symmetry/lattice copies of a peak are listed once, at the image next to the model
                for xyz, v in pk:
                    pos = to_model(xyz)
                    if any(np.linalg.norm(pos - q) < 0.5 for q, *_ in listed):
                        continue
                    listed.append((pos, v))
                    if len(listed) == n_peaks:
                        break
                p(f"  {'negative' if neg else 'positive'} peaks beyond {'-' if neg else '+'}{x.level:g} {u}: {len(pk)} on the "
                  f"stored grid" + (f"; strongest {len(listed)}, each once, next to the model:" if listed else ""))
                for pos, v in listed:
                    nm, d = near(pos)
                    p(f"    {v:+7.2f} {u} at ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})  {d:.2f} A from {nm}")
        for s in self.stops:
            c = np.array(s["center"])
            plane = s["_plane"]
            p(f"\nstarting view: centre {s['_center_spec']!r} = ({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f}), radius {s['radius']:g}, "
              + (f"orientation from the frame atoms; h = height above their plane, positive {plane[2]}" if plane is not None
                 else "orientation " + ("given" if s["orient"] else "kept from the previous stop")))
            for x in self.maps:
                if not x.wrap and not np.isfinite(x.value_at(c)):
                    p(f"  WARNING: centre outside the box of map {x.name!r}")
            for number, i in [tuple(l) for l in s["labels"]]:
                xyz = self._atom_xyz(number, i)
                vals = ", ".join(f"{x.name} {x.value_at(xyz):+.2f}" for x in self.maps)
                h = f"; h {float((xyz - plane[0]) @ plane[1]):+.2f} A" if plane is not None else ""
                p(f"  label {self.atoms.name(number, i)}: {vals}{h}")
            for x in self.maps:
                if x.kind != "diff":
                    continue
                lv = s["maps"].get(x.name)
                shown = (lv is not False and lv is not None) if x.name in s["maps"] else x.visible
                lv = float(lv) if isinstance(lv, (int, float)) and not isinstance(lv, bool) else x.level
                if maps_only and x.name not in maps_only:
                    continue
                p(f"  {x.name} ({'shown' if shown else 'hidden'} at this stop, level {lv:g}), within {s['radius']:g} A:")
                out.extend(self._peak_lines(x, c, s["radius"], lv, near, plane, limit=stop_peaks))
        for pl in self.planes:
            p0, n = np.array(pl["p0"]), np.array(pl["n"])
            ex, ey = np.array(pl["ex"]), np.array(pl["ey"])
            half = max(0.5, pl["thick"] / 2)
            p(f"\nplane {pl['name']!r}: through ({p0[0]:.2f}, {p0[1]:.2f}, {p0[2]:.2f}), normal ({n[0]:.3f}, {n[1]:.3f}, "
              f"{n[2]:.3f}), {2 * pl['r']:g} A square, thickness {pl['thick']:g} A"
              + (f", atoms {pl['_rms']:.2f} A rms out of plane" if pl["_rms"] is not None else ""))
            near_atoms = []
            for number, i, xyz in self._drawn():
                d = xyz - p0
                if abs(d @ n) <= half and abs(d @ ex) <= pl["r"] and abs(d @ ey) <= pl["r"]:
                    near_atoms.append((self.atoms.name(number, i), float(d @ n)))
            p(f"  atoms within {half:g} A of it: {len(near_atoms)}" + (": " + ", ".join(
                f"{nm} {h:+.2f}" for nm, h in near_atoms[:20]) + (" ..." if len(near_atoms) > 20 else "") if near_atoms else ""))
            for x in self.maps:
                if x.kind != "diff" or (pl["map"] and x.name != pl["map"]) or (maps_only and x.name not in maps_only):
                    continue
                lv = float(pl["level"]) if pl["level"] is not None else x.level
                u = x.unit_label or x.units
                for neg in (False, True):
                    found = []
                    for xyz, v in x.peaks(lv, negative=neg):
                        q = x.image_near(xyz, p0)
                        d = q - p0
                        if abs(d @ n) <= half and abs(d @ ex) <= pl["r"] and abs(d @ ey) <= pl["r"]:
                            found.append((q, v, float(d @ n)))
                    p(f"  {x.name}: {len(found)} {'negative' if neg else 'positive'} peaks beyond "
                      f"{'-' if neg else '+'}{lv:g} {u} in the slice" + (", strongest: " + "; ".join(
                          f"{v:+.1f} {near(q)[1]:.2f} A from {near(q)[0]} (h {h:+.2f})" for q, v, h in found[:stop_peaks])
                          if found else ""))
        fplane = None
        if frame:
            pts = self.atoms.xyz(frame)
            if len(pts) >= 3:
                fplane = (pts.mean(0), screen_axes(pts, "face")[2], "on the side the plane normal points to")
        for spec in near_specs:
            c = self.atoms.point(_maybe_xyz(spec))
            p(f"\nnear {spec!r} = ({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f}), within {near_radius:g} A"
              + (f"; h = height above the plane of {frame}" if fplane is not None else "") + ":")
            p("  atoms: " + ", ".join(f"{nm} {d:.2f}" for nm, d in self.atoms_near(c, near_radius)))
            for x in self.maps:
                if x.kind == "diff" and not (maps_only and x.name not in maps_only):
                    p(f"  {x.name}:")
                    out.extend(self._peak_lines(x, c, near_radius, x.level, near, fplane, limit=max(12, stop_peaks)))
        for spec in at:
            number, idx, xyzs = self.atoms.resolve(spec, first_conformer=False)
            p(f"\nvalues at {spec!r}:")
            for i, xyz in zip(idx, xyzs):
                vals = ", ".join(f"{x.name} {x.value_at(xyz):+.2f}" for x in self.maps)
                p(f"  {self.atoms.name(number, i):<22} {vals}")
        return "\n".join(out)

    def atoms_near(self, c, r):
        """[(spec, distance)] of the heavy atoms drawn on the page within r of c, nearest first."""
        out = []
        for k, m in enumerate(self.mols, 1):
            for i, e in enumerate(m.elements):
                d = float(np.linalg.norm(m.xyz[i] - c))
                if e != "H" and d <= r:
                    out.append((self.atoms.name(k, i), d))
        for k, nb in self.atoms.nb.items():
            for i in range(len(nb.orig)):
                d = float(np.linalg.norm(nb.xyz[i] - c))
                if d <= r:
                    out.append((self.atoms.name(self.atoms.nb_number[k], i), d))
        return sorted(out, key=lambda t: t[1])

    def _atom_xyz(self, number, i):
        if number <= len(self.mols):
            return self.mols[number - 1].xyz[i]
        k = next(kk for kk, nn in self.atoms.nb_number.items() if nn == number)
        return self.atoms.nb[k].xyz[i]


def _maybe_xyz(spec):
    """'1.0, -2.5, 3' -> [1.0, -2.5, 3.0]; anything else unchanged."""
    try:
        v = [float(t) for t in spec.replace(",", " ").split()]
        return v if len(v) == 3 else spec
    except ValueError:
        return spec


def _path(p, base):
    p = os.path.expanduser(p)
    return p if os.path.isabs(p) else os.path.join(base, p)


# ------------------------------------------------------------------- init
def _fmt(v):
    return ("+" if v >= 0 else "\u2212") + f"{abs(v):.1f}"


def init_spec(models, maps, title=None, units="sigma", unit_label=None, level=None, min_peak=None, n_planes=5,
              base="."):
    """A ready-to-build spec from files: a starting view of the whole structure and one preset slice plane per
    group of strong features of the first map.

    Features are peaks of either sign beyond `min_peak` within 1.5 A of a model atom. Residues with more than five
    heavy atoms are groups of their own (the plane goes through their atoms around the strongest feature); smaller
    ones (metals, waters, counter-ions) are grouped by name (a metal's plane goes through its coordination)."""
    spec = {"title": title or "Structure and maps", "description": "", "models": [], "maps": [], "planes": []}
    for p in models:
        name = os.path.splitext(os.path.basename(p))[0]
        spec["models"].append({"path": p, "name": name})
    for k, m in enumerate(maps):
        path, _, name = m.partition("=")
        d = {"path": path, "name": name or os.path.splitext(os.path.basename(path))[0], "kind": "diff",
             "units": units, "visible": k == 0}
        if unit_label:
            d["unit_label"] = unit_label
        if level is not None:
            d["level"] = level
        spec["maps"].append(d)
    sh = Share(spec, base)
    if not sh.maps:
        return spec, sh
    x = sh.maps[0]
    lv = x.level
    thr = float(min_peak if min_peak is not None else 2.5 * lv)
    u = x.unit_label or ("\u03c3" if x.units == "sigma" else "")
    sep = "" if u == "\u03c3" else " "
    to_model = sh.atoms.near_model()
    heavy = [(k, i) for k, m in enumerate(sh.mols, 1) for i, e in enumerate(m.elements) if e != "H"]
    tree = cKDTree(np.array([sh.mols[k - 1].xyz[i] for k, i in heavy]))
    feats = []                                   # (|v|, v, pos, model, atom index, distance)
    for neg in (False, True):
        seen = []
        for xyz, v in x.peaks(thr, negative=neg):
            pos = to_model(xyz)
            if any(np.linalg.norm(pos - q) < 0.5 for q in seen):
                continue
            seen.append(pos)
            d, j = tree.query(pos)
            if d <= 1.5:
                feats.append((abs(v), v, pos, *heavy[j], float(d)))
    feats.sort(key=lambda f: -f[0])
    # group: big residues alone, small ones by residue name
    groups = {}
    for f in feats:
        m = sh.mols[f[3] - 1]
        inf = m.info[f[4]]
        rkey = inf[:4]
        nheavy = len({i2[4] for i2, i3 in zip(m.info, m.elements) if i2[:4] == rkey and i3 != "H"})   # altlocs once
        gkey = (f[3], rkey) if nheavy > 5 else (f[3], "name", inf[3])
        groups.setdefault(gkey, []).append(f)
    # the strongest groups; shown with residues of their own (guests, ligands) first, then groups of small residues
    ranked = sorted(groups.items(), key=lambda kv: -kv[1][0][0])[:n_planes]
    ranked.sort(key=lambda kv: (kv[0][1] == "name", -kv[1][0][0]))

    def atom_spec(k, i):
        a = sh.mols[k - 1].info[i]
        sp = f"{a[0]}/{a[1]}{a[2]} {a[4]}" + (f":{a[5]}" if a[5] else "")
        return sp if len(sh.mols) == 1 else f"{k}#{sp}"

    def atom_name(k, i, residue=False):
        a = sh.mols[k - 1].info[i]
        if a[4] == a[3]:                         # one-atom residue named like its atom (PD, ZN)
            return f"{a[3]} {a[1]}{a[2]}"
        return f"{a[4]} ({a[3]} {a[1]}{a[2]})" if residue else a[4]

    def feature_text(f, residue=False):
        where = "on" if f[5] < 0.45 else f"{f[5]:.1f} \u00c5 from"
        return f"{_fmt(f[1])}{sep}{u} {where} {atom_name(f[3], f[4], residue)}"

    # starting view: the whole structure (without waters), face-on
    allxyz = np.concatenate([m.xyz[[i for i, inf in enumerate(m.info) if inf[3] not in ("HOH", "WAT", "DOD")
                                    and m.elements[i] != "H"]] for m in sh.mols])
    c = allxyz.mean(0)
    rad = float(np.percentile(np.linalg.norm(allxyz - c, axis=1), 90)) + 2.0
    spec["view"] = {"center": [round(float(v), 2) for v in c], "radius": round(rad, 1),
                    "orient": _quat(screen_axes(allxyz, "face"))}
    # one preset plane per group of strong features: through the atoms around the strongest one
    drawn = list(sh._drawn())
    for gkey, fs in ranked:
        k = gkey[0]
        m = sh.mols[k - 1]
        top = fs[0]
        inf = m.info[top[4]]
        res_idx = [i for i, i2 in enumerate(m.info) if i2[:4] == inf[:4] and m.elements[i] != "H"
                   and m.info[i][5] in ("", m.info[top[4]][5])]
        hot = m.xyz[top[4]]
        if len(res_idx) >= 3:
            part = [i for i in res_idx if np.linalg.norm(m.xyz[i] - hot) <= 2.6] if len(res_idx) > 12 else res_idx
            atoms = [atom_spec(k, i) for i in part]
            name = f"{inf[3]} {inf[1]}{inf[2]}"
        else:                          # one-atom residue (a metal): the plane of its coordination
            near_atoms = [(nb, i) for nb, i, xyz in drawn if 0.1 < np.linalg.norm(xyz - hot) <= 2.6]
            if len(near_atoms) < 2:
                continue
            atoms = [atom_spec(k, top[4])] + [sh.atoms.name(nb, i) for nb, i in near_atoms]
            name = f"{inf[3]} {inf[1]}{inf[2]} coordination"
        if len(atoms) < 3:
            continue
        spec["planes"].append({"name": name, "atoms": atoms, "map": x.name,
                               "level": thr if top[0] > 3 * thr else lv})
    return spec, sh


def init_main(argv):
    ap = argparse.ArgumentParser(prog="moorlight share init", description=
                                 "Write a ready-to-build share spec from model and map files: an overview and one stop "
                                 "per group of strong features of the first map, with draft texts that state the "
                                 "numbers. Then edit the texts and run `moorlight share spec.json -o page.html`.")
    ap.add_argument("files", nargs="+", help="models (.pdb/.cif) and maps (.ccp4/.map/.mrc/.mtz); name a map with "
                                             "path=name. The first map sets the preset planes.")
    ap.add_argument("-o", "--out", required=True, help="spec to write (.json)")
    ap.add_argument("--title", help="page name, two to four words")
    ap.add_argument("--units", choices=["sigma", "absolute"], default="sigma")
    ap.add_argument("--unit-label", help="text after levels, e.g. SD for maps in noise SD (with --units absolute)")
    ap.add_argument("--level", type=float, help="default contour level (default 3 for difference maps)")
    ap.add_argument("--min", type=float, dest="min_peak", help="features from this level on (default 2.5 x level)")
    ap.add_argument("--planes", "--stops", type=int, default=5, dest="planes", help="preset planes (default 5)")
    ap.add_argument("--build", metavar="PAGE", help="also build the page")
    a = ap.parse_args(argv)
    out_dir = os.path.dirname(os.path.abspath(a.out))
    models, maps = [], []
    for f in a.files:
        p = f.partition("=")[0].lower()
        (models if p.endswith((".pdb", ".ent", ".cif", ".mmcif")) else maps).append(
            os.path.relpath(os.path.abspath(f.partition("=")[0]), out_dir) + (("=" + f.partition("=")[2]) if "=" in f else ""))
    try:
        spec, sh = init_spec(models, maps, a.title, a.units, a.unit_label, a.level, a.min_peak, a.planes, out_dir)
        with open(a.out, "w") as fh:
            json.dump(spec, fh, indent=1, ensure_ascii=False)
        sh = Share(spec, out_dir)
        check = os.path.splitext(a.out)[0] + ".check.txt"
        with open(check, "w") as fh:
            fh.write(sh.report())
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print(f"wrote {a.out} ({len(spec['planes'])} preset planes: " + ", ".join(pl["name"] for pl in spec["planes"])
          + f") and {check}")
    if a.build:
        page = sh.html()
        with open(a.build, "w") as fh:
            fh.write(page)
        print(f"wrote {a.build}: {len(page.encode()) / 1e6:.2f} MB")
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "init":
        return init_main(argv[1:])
    # specs may start with '-' (chain '-'): '--at -/102' -> '--at=-/102'
    fixed, k = [], 0
    while k < len(argv):
        if argv[k] in ("--at", "--near", "--frame") and k + 1 < len(argv):
            fixed.append(f"{argv[k]}={argv[k + 1]}")
            k += 2
        else:
            fixed.append(argv[k])
            k += 1
    ap = argparse.ArgumentParser(prog="moorlight share", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec", help="share spec (JSON)")
    ap.add_argument("-o", "--out", help="output page (.html); default: next to the spec")
    ap.add_argument("--check", action="store_true", help="print peaks and values for writing the stop texts; no page")
    ap.add_argument("--at", action="append", default=[], metavar="SPEC",
                    help="with --check: map values at these atoms or residues (repeatable)")
    ap.add_argument("--near", action="append", default=[], metavar="SPEC",
                    help="with --check: peaks of every difference map around an atom, residue or 'x,y,z' (repeatable)")
    ap.add_argument("--radius", type=float, default=5.0, help="with --near: radius in A (default 5)")
    ap.add_argument("--peaks", type=int, default=12, help="with --check: peaks listed per map and sign")
    ap.add_argument("--stop-peaks", type=int, default=8, help="with --check: peaks listed per stop, map and sign")
    ap.add_argument("--map", action="append", default=[], dest="maps", metavar="NAME",
                    help="with --check: list only this map (repeatable)")
    ap.add_argument("--frame", metavar="SPECS", help="with --near: atoms (comma-separated) whose plane gives h")
    ap.add_argument("--standalone", action="store_true",
                    help="full HTML document for opening as a local file (artifact pages must not have it)")
    a = ap.parse_args(fixed)
    with open(a.spec) as fh:
        spec = json.load(fh)
    try:
        sh = Share(spec, os.path.dirname(os.path.abspath(a.spec)))
        if a.check:
            frame = [t.strip() for t in a.frame.split(",")] if a.frame else None
            print(sh.report(a.at, a.near, a.radius, a.peaks, a.stop_peaks, a.maps, frame))
            return 0
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    out = a.out or os.path.splitext(a.spec)[0] + ".html"
    page = sh.html(standalone=a.standalone)
    with open(out, "w") as fh:
        fh.write(page)
    size = len(page.encode())
    print(f"wrote {out}: {size / 1e6:.2f} MB, {len(sh.mols)} model(s) + {len(sh.atoms.nb)} neighbour set(s), "
          f"{len(sh.maps)} map(s) ({', '.join(x.mode for x in sh.maps)}), {len(sh.planes)} preset plane(s)")
    if size > MAX_BYTES:
        print(f"warning: larger than {MAX_BYTES / 1e6:.0f} MB; an artifact page must stay at or below 16 MB. "
              "Use fewer maps, or grid 'box' with a smaller region or a coarser step.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
