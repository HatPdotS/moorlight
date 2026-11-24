"""Shareable snapshot: one self-contained HTML viewer page with models, maps and preset slice planes.

No server is needed to view it. The page contours the maps itself (three.js from a CDN), so it can be sent as a
file or published as a claude.ai artifact and opened on a phone.

    moorlight share init model.pdb maps.mtz -o spec.json --build view.html    # draft spec (and page) from files
    moorlight share spec.json -o view.html        # build the page
    moorlight share spec.json --check             # what the page shows: peaks, planes, values at atoms
    moorlight share spec.json --check --at "A/85 SG" --near "A/301" --radius 4

The spec (JSON; paths relative to the spec file):

    {"title": "Cage 3.5 ps",                      # page name, two to four words
     "description": "one line shown under the title",
     "models": [{"path": "dark.pdb", "name": "dark"}],
     "maps": [{"path": "3p5ps.ccp4", "name": "3.5 ps light - dark", "kind": "diff",
               "units": "absolute", "unit_label": "SD", "level": 4.0, "visible": true}],
     "view": {"center": "-/102", "radius": 8, "orient": "face", "labels": ["-/102 N9"]},
     "planes": [{"name": "Azo 102", "atoms": ["-/102 C8", "-/102 N9", "-/102 N10", "-/102 C11"],
                 "thickness": 0.5, "map": "3.5 ps light - dark"},
                {"name": "Pd site", "coordination": "A/1 PD"}]}

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

from .maps import DEFAULT_VISIBLE, SymMapper, detect_columns, check_columns, guess_kind
from .model import Molecule
from .session import parse_atom_spec

WEB = os.path.join(os.path.dirname(__file__), "web")
MAX_BYTES = 15 * 1024 * 1024    # artifact pages must stay at or below 16 MB
POINT_BUDGET = 8_000_000        # grid points stored for all maps together (~1.4 bytes each on the page)
NEIGHBOUR_COLOR = "#8a9099"
WATERS = ("HOH", "WAT", "DOD", "H2O")
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
                tag = (("" if ident else op.triplet() + " ") + _shift_tag(t)).strip()
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


def coordination_cutoff(element):
    """Default metal-ligand cutoff (A): 2.6, or 2.8 for alkali and alkaline-earth ions (Ca-O up to ~2.7)."""
    return 2.8 if element.upper() in ("NA", "K", "RB", "CS", "CA", "SR", "BA") else 2.6


def _plane_rms(pts):
    """rms distance of points from their best plane."""
    c = pts - pts.mean(0)
    return float(np.sqrt(max(np.linalg.eigvalsh(c.T @ c)[0], 0.0) / len(pts)))


def _is_xyz(v):
    return isinstance(v, (list, tuple)) and len(v) == 3 and all(isinstance(x, (int, float)) for x in v)


# ------------------------------------------------------------------- maps
def read_map(spec, base):
    """Map spec -> (gemmi.FloatGrid over the unit cell, kind or None). CCP4/MRC maps (a map that covers part of
    the cell has NaN elsewhere), or MTZ with columns 'F,PHI'."""
    p = _path(spec["path"], base)
    low = p.lower()
    if low.endswith((".ccp4", ".map", ".mrc")):
        g = gemmi.read_ccp4_map(p, setup=True).grid
        return g, spec.get("kind")
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


def guess_map_kind(path, values):
    """'2fofc' or 'diff' from the file name (2mFo-DFc, FWT, emd_... / mFo-DFc, DELFWT, diff, minus, omit, ...),
    else from the values. Density maps have almost nothing below -2.5 rms (0.000-0.004 % of the grid in 2mFo-DFc
    maps); difference maps have noise there (0.4-1.4 %, also light - dark maps with strong positive features)."""
    name = os.path.basename(path).lower()
    if re.search(r"2m?fo|(?<!del)fwt|^emd", name):
        return "2fofc"
    if re.search(r"fo-?fc|mfo-?dfc|delfwt|diff|fo-?fo|minus|omit|polder|delta", name):
        return "diff"
    v = values[:: max(1, values.size // 1_000_000)].astype(np.float64)
    return "2fofc" if np.mean(v < v.mean() - 2.5 * v.std()) < 0.001 else "diff"


class PageMap:
    """A map as the page stores it: values on a grid with positions xyz = origin + G @ (i, j, k).

    native: the map's own grid over the unit cell, periodic, nothing resampled (values as in Coot).
    block: the part of that grid that covers the models and a margin, nothing resampled.
    box: resampled on an orthogonal box (when even the block is too large for the page)."""

    def __init__(self, spec, base, box_lo, box_hi, budget=POINT_BUDGET):
        grid, kind = read_map(spec, base)
        whole = np.array(grid, copy=False)
        self.name = spec.get("name") or os.path.basename(spec["path"])
        finite = np.isfinite(whole)
        self.coverage = float(finite.mean())
        vals = whole[finite] if self.coverage < 1 else whole.ravel()
        if not vals.size:
            raise ValueError(f"map {self.name!r}: no values")
        self.mean, self.rms = float(vals.mean()), float(vals.std())
        if self.coverage < 1:                  # gemmi fills the cell outside the map's extent with NaN
            whole = np.where(finite, whole, np.float32(self.mean))
        kind = spec.get("kind") or kind or guess_map_kind(spec["path"], vals)
        self.kind = {"density": "2fofc", "em": "2fofc"}.get(kind, kind)
        if self.kind not in ("diff", "2fofc"):
            raise ValueError(f"map {self.name!r}: kind must be 'diff' or '2fofc' (or 'density')")
        self.units = spec.get("units", "sigma")
        if self.units not in ("sigma", "absolute"):
            raise ValueError(f"map {self.name!r}: units must be 'sigma' or 'absolute'")
        self.unit_label = spec.get("unit_label") or ("σ" if self.units == "sigma" else "")
        self.level = float(spec.get("level", 3.0 if self.kind == "diff" else 1.5))
        self.visible = bool(spec.get("visible", True))
        self.color = spec.get("color")
        self.source = os.path.basename(spec["path"])
        self.cell = np.array(grid.unit_cell.orth.mat.tolist())
        n_cell = np.array(whole.shape)
        g_cell = self.cell / n_cell                     # column a: one grid step along axis a
        lo, hi = np.asarray(box_lo, float), np.asarray(box_hi, float)
        corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
        g = np.linalg.solve(g_cell, corners.T).T        # the region in grid coordinates
        b_lo, b_hi = np.floor(g.min(0)).astype(int), np.ceil(g.max(0)).astype(int)
        block = int(np.prod(b_hi - b_lo + 1))
        mode = spec.get("grid", "auto")
        if mode == "auto":                              # smallest exact choice that fits the page
            if block < 0.5 * whole.size and block <= budget:
                mode = "block"
            elif whole.size <= budget:
                mode = "native"
            else:
                mode = "block" if block <= budget else "box"
        self.mode = mode
        self.filled = 0.0                               # fraction of stored points outside the map's extent
        if mode == "native":
            self.arr = np.array(whole, dtype=np.float32)
            self.n = n_cell
            self.G = g_cell
            self.origin = np.zeros(3)
            self.wrap = True
            self.filled = 1 - self.coverage
        elif mode == "block":
            idx = np.ix_(*[np.arange(a, b + 1) % n for a, b, n in zip(b_lo, b_hi, n_cell)])
            self.arr = np.array(whole[idx], dtype=np.float32)
            self.n = np.array(self.arr.shape)
            self.G = g_cell
            self.origin = g_cell @ b_lo
            self.wrap = False
            self.filled = float(1 - finite[idx].mean())
        elif mode == "box":
            uc = grid.unit_cell
            step = float(spec.get("step") or max(min(uc.a / grid.nu, uc.b / grid.nv, uc.c / grid.nw) / 2,
                                                 (np.prod(hi - lo) / budget) ** (1 / 3)))
            self.n = np.maximum(2, np.ceil((hi - lo) / step).astype(int) + 1)
            self.arr = np.zeros(tuple(self.n), dtype=np.float32)
            tr = gemmi.Transform()
            tr.mat.fromlist([[step, 0, 0], [0, step, 0], [0, 0, step]])
            tr.vec.fromlist([float(v) for v in lo])
            if self.coverage < 1:
                grid = gemmi.FloatGrid(whole, grid.unit_cell, grid.spacegroup)
            grid.interpolate_values(self.arr, tr, order=3)   # tricubic: trilinear on a new grid flattens peaks
            self.G = np.eye(3) * step
            self.origin = lo
            self.wrap = False
        else:
            raise ValueError(f"map {self.name!r}: grid must be 'auto', 'native', 'block' or 'box'")
        # the page stores int16 steps of 1/100 of the map's rms (wider for maps with very high peaks); the values
        # here are rounded the same way, so every number --check prints is the one the page shows
        self.scale = max(self.rms / 100, float(np.abs(self.arr).max()) / 32000, 1e-12)
        self._raw = self.arr                   # peaks are found before rounding, which would make plateaus
        self.arr = (np.round(self.arr / self.scale) * self.scale).astype(np.float32)
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
            a = self.to_level_units(self._raw)
            a = -a if negative else a
            mx = maximum_filter(a, size=3, mode="wrap" if self.wrap else "nearest")
            hit = (a == mx) & (a > level)
            if not self.wrap:                  # the edge of a block or box is not a peak
                hit[[0, -1]] = hit[:, [0, -1]] = hit[:, :, [0, -1]] = False
            idx = np.argwhere(hit)
            val = self.to_level_units(self.arr[tuple(idx.T)]) * (-1 if negative else 1)   # the page's values
            order = np.argsort(-val)
            xyz = self.origin + idx[order] @ self.G.T
            self._peaks[key] = [(p, float(-v if negative else v)) for p, v in zip(xyz, val[order])]
        return self._peaks[key]

    def to_json(self):
        scale = self.scale
        q = np.ascontiguousarray(np.round(self.arr / scale).astype("<i2")).ravel().view(np.uint8).reshape(-1, 2)
        # high bytes first, then low bytes: the high bytes are mostly alike and compress ~15 % better apart
        data = base64.b64encode(gzip.compress(q[:, 1].tobytes() + q[:, 0].tobytes(), compresslevel=9, mtime=0)).decode()
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
        map_specs = spec.get("maps") or []
        self.maps = [PageMap(s, base, lo, hi, POINT_BUDGET / max(1, len(map_specs))) for s in map_specs]
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
            base = self.planes[-1]["name"]
            for o in p.get("copies") or []:      # parallel copies at these offsets (Å) along the normal
                name = f"{base} {o:+g} \u00c5"
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
        """Preset slice plane: the best plane through `atoms` (three or more), through an atom and the atoms within
        `cutoff` of it (`coordination`: a metal's coordination plane), or `center` with a `normal` (a vector, or two
        atom specs from the first to the second). `axis` (two atom specs) sets its horizontal direction; `offset`
        moves it along the normal; `thickness`, `radius`, `map`, `level`, `color`, `visible`."""
        where = f"plane {k + 1} ({p.get('name', '')!r})"
        try:
            rms, pts, centre, ligands, shift = None, None, None, None, 0.0
            if p.get("coordination"):
                number, idx, xyz = self.atoms.resolve(p["coordination"])
                if len(idx) != 1:
                    raise ValueError("\"coordination\" is one atom, e.g. \"A/1 PD\"")
                cut = float(p.get("cutoff") or coordination_cutoff(self._element(number, idx[0])))
                ligands = [(nb, i, float(np.linalg.norm(x - xyz[0]))) for nb, i, x in self._drawn()
                           if 0.1 < np.linalg.norm(x - xyz[0]) <= cut]
                if len(ligands) < 2:
                    raise ValueError(f"fewer than two atoms within {cut:g} A of {p['coordination']}; raise \"cutoff\"")
                ligands.sort(key=lambda t: t[2])
                lx = np.array([self._atom_xyz(nb, i) for nb, i, _ in ligands])
                use = list(range(len(ligands)))
                if len(ligands) >= 5 and p.get("equatorial", True):
                    # octahedral and higher: the four ligands that lie flattest with the metal (equatorial plane)
                    use = list(min(itertools.combinations(use, 4), key=lambda q: _plane_rms(np.vstack([xyz, lx[list(q)]]))))
                pts = np.vstack([xyz, lx[use]])
                centre = xyz[0]
                ligands = [(self.atoms.name(nb, i), d, j in use) for j, (nb, i, d) in enumerate(ligands)]
            elif p.get("atoms"):
                pts = self.atoms.xyz(p["atoms"])
                if len(pts) < 3:
                    raise ValueError("\"atoms\" needs three or more atoms (one atom: use \"coordination\")")
            if pts is not None:
                ex, ey, n = screen_axes(pts, "face")
                if "center" in p:
                    centre = self.atoms.point(p["center"])
                elif centre is None:
                    centre = pts.mean(0)
                rms = _plane_rms(pts)
                shift = float((centre - pts.mean(0)) @ n)
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
                raise ValueError("give \"atoms\" (three or more), \"coordination\" (one atom), or \"center\" and "
                                 "\"normal\"")
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
                "visible": bool(p.get("visible", True)), "_rms": rms, "_ligands": ligands, "_shift": shift}

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
        for neg in ((False, True) if x.kind == "diff" else (False,)):
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
        """What the page shows, in numbers. maps_only: names of the maps to list (default all); frame: atom specs whose
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
            grid = {"native": f"own grid {x.n.tolist()} over the cell, nothing resampled",
                    "block": f"block {x.n.tolist()} of its own grid around the models, nothing resampled",
                    "box": f"resampled on a box {x.n.tolist()} x {x.G[0, 0]:.2f} A (tricubic)"}[x.mode]
            if x.units == "sigma":     # the level on the map's own scale too
                on_scale = (f"±{x.level * x.rms:.3g}" if x.kind == "diff" else f"{x.mean + x.level * x.rms:.3g}") + " in map units"
            else:
                on_scale = f"{(x.level - (0 if x.kind == 'diff' else x.mean)) / x.rms:.3g} x rms"
            p(f"\nmap {x.name!r} ({x.source}, {x.kind}, {grid}): level {x.level:g} {u} = {on_scale}; whole map mean "
              f"{x.mean:.3g} rms {x.rms:.3g}; stored values {vals.min():.2f} to {vals.max():.2f} {u}")
            if x.filled > 0.001:
                p(f"  WARNING: the map file covers {100 * x.coverage:.0f} % of the cell; {100 * x.filled:.0f} % of the "
                  "stored grid lies outside it and shows the map's mean")
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
              + (f", atoms {pl['_rms']:.2f} A rms about their best plane" if pl["_rms"] is not None else "")
              + (f"; moved {pl['_shift']:+.2f} A along the normal to pass through the centre" if abs(pl["_shift"]) > 0.01
                 else ""))
            near_atoms = []
            for number, i, xyz in self._drawn():
                d = xyz - p0
                if abs(d @ n) <= half and abs(d @ ex) <= pl["r"] and abs(d @ ey) <= pl["r"]:
                    near_atoms.append((self.atoms.name(number, i), float(d @ n)))
            if pl["_ligands"]:
                inp = [f"{nm} {d:.2f} A" for nm, d, used in pl["_ligands"] if used]
                outp = [f"{nm} {d:.2f} A" for nm, d, used in pl["_ligands"] if not used]
                p("  coordination" + (" (equatorial plane)" if outp else "") + ": " + ", ".join(inp)
                  + ("; out of the plane: " + ", ".join(outp) if outp else ""))
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
            p("  atoms (distance" + (", h" if fplane is not None else "") + "): " + ", ".join(
                f"{nm} {d:.2f}" + (f" h {float((xyz - fplane[0]) @ fplane[1]):+.2f}" if fplane is not None else "")
                for nm, d, xyz in self.atoms_near(c, near_radius, with_xyz=True)))
            for x in self.maps:
                if not (maps_only and x.name not in maps_only):
                    p(f"  {x.name}:")
                    out.extend(self._peak_lines(x, c, near_radius, x.level, near, fplane, limit=max(12, stop_peaks)))
        for spec in at:
            number, idx, xyzs = self.atoms.resolve(spec, first_conformer=False)
            p(f"\nvalues at {spec!r}:")
            for i, xyz in zip(idx, xyzs):
                vals = ", ".join(f"{x.name} {x.value_at(xyz):+.2f}" for x in self.maps)
                p(f"  {self.atoms.name(number, i):<22} {vals}")
        return "\n".join(out)

    def atoms_near(self, c, r, with_xyz=False):
        """[(spec, distance)] (or (spec, distance, xyz)) of the heavy atoms drawn on the page within r of c,
        nearest first."""
        out = []
        for number, i, xyz in self._drawn():
            d = float(np.linalg.norm(xyz - c))
            if d <= r:
                out.append((self.atoms.name(number, i), d, xyz))
        out.sort(key=lambda t: t[1])
        return out if with_xyz else [t[:2] for t in out]

    def _element(self, number, i):
        if number <= len(self.mols):
            return self.mols[number - 1].elements[i]
        nb = self.atoms.nb[next(kk for kk, nn in self.atoms.nb_number.items() if nn == number)]
        return nb.mol.elements[nb.orig[i]]

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
def _map_entries(path, name, base):
    """Spec entries for one map file. An MTZ gives every map a refinement program writes for viewing (2mFo-DFc and
    mFo-DFc), each with its columns; `name` then prefixes the map names."""
    stem = os.path.splitext(os.path.basename(path))[0]
    if not path.lower().endswith(".mtz"):
        return [{"path": path, "name": name or stem}]
    if ":" in name:                                  # file.mtz=name:F,PHI
        name, _, cols = name.partition(":")
        return [{"path": path, "columns": cols, "name": name or stem}]
    found = detect_columns(gemmi.read_mtz_file(_path(path, base)))
    if not found:
        raise ValueError(f"{path}: no map columns found; give them as {path}=name:F,PHI")
    shown = [f for f in found if f[0] in DEFAULT_VISIBLE] or found[:1]
    return [{"path": path, "columns": f"{f},{phi}", "kind": kind,
             "name": (f"{name} {label}" if len(shown) > 1 else name) if name else label}
            for f, phi, kind, label in shown]


def _residue_spec(sh, k, i, atom=False):
    a = sh.mols[k - 1].info[i]
    sp = f"{a[0]}/{a[1]}{a[2]}" + (f" {a[4]}" + (f":{a[5]}" if a[5] else "") if atom else "")
    return sp if len(sh.mols) == 1 else f"{k}#{sp}"


def _feature_planes(sh, x, thr, n):
    """Preset planes at the strongest features of difference map x: peaks of either sign beyond thr within 1.5 A of
    a model atom. A residue with more than five heavy atoms is a group of its own (the plane goes through its atoms
    around the strongest feature); smaller ones (metals, waters, ions) are grouped by name (a metal's plane is its
    coordination)."""
    to_model = sh.atoms.near_model()
    heavy = [(k, i) for k, m in enumerate(sh.mols, 1) for i, e in enumerate(m.elements) if e != "H"]
    tree = cKDTree(np.array([sh.mols[k - 1].xyz[i] for k, i in heavy]))
    size = {}                                    # (model, residue) -> heavy atom names (altlocs once)
    for k, i in heavy:
        size.setdefault((k, tuple(sh.mols[k - 1].info[i][:4])), set()).add(sh.mols[k - 1].info[i][4])
    feats = []                                   # (|v|, v, model, atom index)
    for neg in (False, True):
        seen = []
        for xyz, v in x.peaks(thr, negative=neg):
            pos = to_model(xyz)
            if any(np.linalg.norm(pos - q) < 0.5 for q in seen):
                continue
            seen.append(pos)
            d, j = tree.query(pos)
            if d <= 1.5:
                feats.append((abs(v), v, *heavy[j]))
    feats.sort(key=lambda f: -f[0])
    groups = {}
    for f in feats:
        inf = sh.mols[f[2] - 1].info[f[3]]
        rkey = tuple(inf[:4])
        big = len(size[(f[2], rkey)]) > 5
        if big or inf[3] not in WATERS:
            groups.setdefault((f[2], rkey) if big else (f[2], "name", inf[3]), []).append(f)
    ranked = sorted(groups.items(), key=lambda kv: -kv[1][0][0])[:n]
    ranked.sort(key=lambda kv: (kv[0][1] == "name", -kv[1][0][0]))   # residues of their own first
    planes = []
    for gkey, fs in ranked:
        top = fs[0]
        k, m = top[2], sh.mols[top[2] - 1]
        inf = m.info[top[3]]
        name = f"{inf[3]} {inf[1]}{inf[2]}"
        res_idx = [i for i, i2 in enumerate(m.info) if i2[:4] == inf[:4] and m.elements[i] != "H"
                   and i2[5] in ("", inf[5])]
        if len(res_idx) >= 3:
            hot = m.xyz[top[3]]
            part = [i for i in res_idx if np.linalg.norm(m.xyz[i] - hot) <= 2.6] if len(res_idx) > 12 else res_idx
            if len(part) >= 3:
                planes.append({"name": name, "atoms": [_residue_spec(sh, k, i, atom=True) for i in part],
                               "level": thr if top[0] > 3 * thr else x.level})
        elif sum(1 for *_, xyz in sh._drawn() if 0.1 < np.linalg.norm(xyz - m.xyz[top[3]]) <= 2.6) >= 2:
            planes.append({"name": f"{name} coordination", "coordination": _residue_spec(sh, k, top[3], atom=True),
                           "level": thr if top[0] > 3 * thr else x.level})
    return planes


def _ligand_planes(sh, n):
    """Preset planes through the ligands (residues that are neither polymer nor water, largest first), then the
    coordination of single-atom residues (metals, ions) that have two or more atoms within 2.6 A."""
    ligands, ions = [], []
    drawn = None
    for k, m in enumerate(sh.mols, 1):
        residues = {}
        for i, (inf, e) in enumerate(zip(m.info, m.elements)):
            if e != "H" and inf[3] not in WATERS and inf[5] in ("", "A"):
                residues.setdefault(tuple(inf[:4]), []).append(i)
        for rkey, idx in residues.items():
            tab = gemmi.find_tabulated_residue(rkey[3])
            if tab is not None and (tab.is_amino_acid() or tab.is_nucleic_acid()):
                continue
            name = f"{rkey[3]} {rkey[1]}{rkey[2]}"
            if len(idx) >= 3:
                ligands.append((len(idx), {"name": name, "atoms": [_residue_spec(sh, k, idx[0])]}))
            elif len(idx) == 1:
                drawn = drawn if drawn is not None else np.array([xyz for *_, xyz in sh._drawn()])
                if np.sum((np.linalg.norm(drawn - m.xyz[idx[0]], axis=1) > 0.1)
                          & (np.linalg.norm(drawn - m.xyz[idx[0]], axis=1) <= 2.6)) >= 2:
                    ions.append({"name": f"{name} coordination", "coordination": _residue_spec(sh, k, idx[0], atom=True)})
    return ([p for _, p in sorted(ligands, key=lambda t: -t[0])] + ions)[:n]


def _requested_plane(sh, text):
    """'NAME=SPEC, SPEC, ...' -> a plane: one atom gives its coordination plane, more give their best plane."""
    name, _, specs = text.partition("=") if "=" in text else ("", "", text)
    items = [t.strip() for t in specs.split(",") if t.strip()]
    if not items:
        raise ValueError(f"--plane {text!r}: give NAME=SPEC, SPEC, ... (atom or residue specs)")
    try:
        count = sum(len(sh.atoms.resolve(t)[1]) for t in items)
    except ValueError as e:
        raise ValueError(f"--plane {text!r}: {e}") from None
    d = {"name": name.strip() or ", ".join(items)}
    d.update({"coordination": items[0]} if count == 1 else {"atoms": items})
    return d


def init_spec(models, maps, title=None, units="sigma", unit_label=None, level=None, min_peak=None, n_planes=None,
              planes=(), base=".", description=""):
    """A ready-to-build spec from files: a starting view of the whole structure and preset slice planes.

    Planes: the requested ones ('NAME=SPEC, SPEC, ...'), then up to n_planes automatic ones (default 5, or 0 when
    planes are requested): at the strongest features of the first difference map, or, without a difference map or
    features, through the ligands and the coordination of metals."""
    spec = {"title": title or "Structure and maps", "description": description or "", "models": [], "maps": [],
            "planes": []}
    for m in models:
        p, _, name = m.partition("=")
        spec["models"].append({"path": p, "name": name or os.path.splitext(os.path.basename(p))[0]})
    entries = []
    for m in maps:
        path, _, name = m.partition("=")
        entries.extend(_map_entries(path, name, base))
    names = [e["name"] for e in entries]
    for e in entries:
        if names.count(e["name"]) > 1:
            e["name"] = f"{os.path.splitext(os.path.basename(e['path']))[0]} {e['name']}"
        e["units"] = units
        if unit_label:
            e["unit_label"] = unit_label
    spec["maps"] = entries
    sh = Share(spec, base)
    kinds = set()
    any_diff = any(x.kind == "diff" for x in sh.maps)
    for e, x in zip(entries, sh.maps):
        e["kind"] = x.kind
        e["level"] = level if level is not None and (x.kind == "diff" or not any_diff) else x.level
        e["visible"] = x.kind not in kinds                  # the first map of each kind
        kinds.add(x.kind)
        x.level = e["level"]
    diff = next((x for x in sh.maps if x.kind == "diff"), None)
    shown = diff or (sh.maps[0] if sh.maps else None)
    for text in planes:
        spec["planes"].append(_requested_plane(sh, text))
    n_auto = (0 if planes else 5) if n_planes is None else n_planes
    auto = []
    if diff is not None and n_auto:
        auto = _feature_planes(sh, diff, float(min_peak if min_peak is not None else 2.5 * diff.level), n_auto)
    if not auto and n_auto:
        auto = _ligand_planes(sh, n_auto)
    spec["planes"] += auto
    if shown is not None:
        for pl in spec["planes"]:
            pl["map"] = shown.name
    # starting view: the whole structure (without waters), face-on
    allxyz = np.concatenate([m.xyz[[i for i, inf in enumerate(m.info) if inf[3] not in WATERS
                                    and m.elements[i] != "H"]] for m in sh.mols])
    c = allxyz.mean(0)
    rad = float(np.percentile(np.linalg.norm(allxyz - c, axis=1), 90)) + 2.0
    spec["view"] = {"center": [round(float(v), 2) for v in c], "radius": round(rad, 1),
                    "orient": _quat(screen_axes(allxyz, "face"))}
    return spec


def init_main(argv):
    ap = argparse.ArgumentParser(prog="moorlight share init", description=
                                 "Write a ready-to-build share spec from model and map files: a starting view of the "
                                 "whole structure and preset slice planes, plus the --check report next to it. Then "
                                 "edit the spec if needed and run `moorlight share spec.json -o page.html`.")
    ap.add_argument("files", nargs="+", help="models (.pdb/.cif) and maps (.ccp4/.map/.mrc/.mtz). Name a file with "
                                             "path=name; an MTZ gives its 2mFo-DFc and mFo-DFc maps (path=name:F,PHI "
                                             "picks columns)")
    ap.add_argument("-o", "--out", required=True, help="spec to write (.json)")
    ap.add_argument("--title", help="page name, two to four words")
    ap.add_argument("--description", help="one line under the title")
    ap.add_argument("--plane", action="append", default=[], metavar="NAME=SPECS",
                    help="a preset plane: 'Azo=-/102 C8, -/102 N9, -/102 N10' (best plane through the atoms), "
                         "'Heme=A/201' (a residue's heavy atoms) or 'Zn site=A/301 ZN' (one atom: its coordination "
                         "plane). Repeatable")
    ap.add_argument("--planes", type=int, default=None, dest="planes",
                    help="automatic preset planes (default 5, or 0 when --plane is given)")
    ap.add_argument("--units", choices=["sigma", "absolute"], default="sigma")
    ap.add_argument("--unit-label", help="text after levels, e.g. SD for maps in noise SD (with --units absolute)")
    ap.add_argument("--level", type=float, help="contour level of the difference maps (default 3 sigma)")
    ap.add_argument("--min", type=float, dest="min_peak", help="features from this level on (default 2.5 x level)")
    ap.add_argument("--build", metavar="PAGE", help="also build the page")
    a = ap.parse_args(argv)
    out_dir = os.path.dirname(os.path.abspath(a.out))
    os.makedirs(out_dir, exist_ok=True)
    models, maps = [], []
    for f in a.files:
        path, eq, name = f.partition("=")
        if not os.path.exists(path):
            print(f"error: no file {path}", file=sys.stderr)
            return 1
        rel = os.path.relpath(os.path.abspath(path), out_dir) + (eq + name)
        (models if path.lower().endswith((".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".cif.gz")) else maps).append(rel)
    try:
        spec = init_spec(models, maps, a.title, a.units, a.unit_label, a.level, a.min_peak, a.planes, a.plane, out_dir,
                         a.description)
        with open(a.out, "w") as fh:
            json.dump(spec, fh, indent=1, ensure_ascii=False)
        sh = Share(spec, out_dir)
        check = os.path.splitext(a.out)[0] + ".check.txt"
        with open(check, "w") as fh:
            fh.write(sh.report())
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print(f"wrote {a.out}: {len(spec['maps'])} map(s) (" + ", ".join(f"{m['name']} [{m['kind']}]" for m in spec["maps"])
          + f"), {len(spec['planes'])} preset plane(s) (" + ", ".join(pl["name"] for pl in spec["planes"])
          + f"); and {check}")
    if a.build:
        page = sh.html()
        with open(a.build, "w") as fh:
            fh.write(page)
        print(f"wrote {a.build}: {len(page.encode()) / 1e6:.2f} MB")
    return 0


def list_model(path):
    """What a model holds, in spec form: chains with their residue ranges, ligands, metal and other ions with their
    coordination, waters and residues with alternate conformers."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    out = [f"{path}: " + (f"{st.spacegroup_hm}, cell " + " ".join(f"{v:g}" for v in st.cell.parameters)
                          if st.cell.is_crystal() else "no crystal cell") + f"; {len(st)} model(s), the first listed"]
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate()
    ligands, ions, waters, alts = {}, {}, 0, []
    for ch in st[0]:
        poly = [r for r in ch if (t := gemmi.find_tabulated_residue(r.name)) is not None
                and (t.is_amino_acid() or t.is_nucleic_acid())]
        if poly:
            out.append(f"  chain {ch.name}: {len(poly)} polymer residues, {ch.name}/{poly[0].seqid.num}"
                       f"{poly[0].seqid.icode.strip()} {poly[0].name} to {ch.name}/{poly[-1].seqid.num}"
                       f"{poly[-1].seqid.icode.strip()} {poly[-1].name}")
        for r in ch:
            spec = f"{ch.name}/{r.seqid.num}{r.seqid.icode.strip()}"
            if any(a.altloc != "\0" for a in r):
                alts.append(spec)
            if r in poly:
                continue
            if r.name in WATERS:
                waters += 1
                continue
            heavy = [a for a in r if a.element.name != "H"]
            if len(heavy) == 1:
                a = heavy[0]
                cut = coordination_cutoff(a.element.name)
                near = set()
                for m in ns.find_atoms(a.pos, "\0", radius=cut):   # marks sit in the cell; measure to the nearest image
                    d = st.cell.find_nearest_image(a.pos, m.pos, gemmi.Asu.Any).dist()
                    if m.to_cra(st[0]).atom.element.name != "H" and 0.1 < d <= cut:
                        near.add((m.chain_idx, m.residue_idx, m.atom_idx, round(d, 3)))
                ions.setdefault((r.name, a.name, cut), []).append(f"{spec} ({len(near)})")
            else:
                ligands.setdefault((r.name, len(heavy)), []).append(spec)
    for (name, n), specs in ligands.items():
        out.append(f"  ligand {name} ({n} heavy atoms): " + ", ".join(specs[:40]) + (" ..." if len(specs) > 40 else ""))
    for (name, atom, cut), specs in ions.items():
        out.append(f"  {name} (atom {atom}; atoms within {cut:g} A): " + ", ".join(specs[:40])
                   + (" ..." if len(specs) > 40 else ""))
    out.append(f"  waters: {waters}")
    if alts:
        out.append(f"  alternate conformers in {len(alts)} residues: " + ", ".join(alts[:30])
                   + (" ..." if len(alts) > 30 else "") + " (specs without :B take the first conformer)")
    return "\n".join(out)


def list_main(argv):
    ap = argparse.ArgumentParser(prog="moorlight share list", description=
                                 "List chains, residue ranges, ligands, metals and waters of models, as atom specs.")
    ap.add_argument("models", nargs="+")
    a = ap.parse_args(argv)
    try:
        for m in a.models:
            print(list_model(m))
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "init":
        return init_main(argv[1:])
    if argv and argv[0] == "list":
        return list_main(argv[1:])
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
    ap.add_argument("--check", action="store_true", help="print what the page shows (peaks, planes, values at atoms); no page")
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
