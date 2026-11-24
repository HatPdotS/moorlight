"""Shareable snapshot: one self-contained HTML page with models, map boxes and tour stops.

No server is needed to view it. The page contours the maps itself (three.js from a CDN), so it can be sent as a
file or published as a claude.ai artifact and opened on a phone.

    moorlight share spec.json -o view.html        # build the page
    moorlight share spec.json --check             # numbers for writing stop texts: peaks, values at atoms
    moorlight share spec.json --check --at "-/102 N1" --at "A/1"

The spec (JSON; paths relative to the spec file):

    {"title": "Cage 3.5 ps",                      # page name, two to four words
     "description": "one line shown under the title",
     "models": [{"path": "dark.pdb", "name": "dark"}],
     "maps": [{"path": "3p5ps.ccp4", "name": "3.5 ps light - dark", "kind": "diff",
               "units": "absolute", "unit_label": "SD", "level": 4.0, "visible": true}],
     "contour_radius": 10,
     "stops": [{"name": "Pd cage", "text": "...", "center": "A/1 PD", "radius": 8,
                "orient": "face", "frame": ["A/1 PD", "A/2 PD", "A/3 PD"],
                "maps": {"3.5 ps light - dark": 4.0}, "labels": ["A/1 PD"], "distances": [["A/1 PD", "A/2 PD"]]}]}

See skills/moorlight-share/SKILL.md for every field.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import html
import json
import os
import sys

import gemmi
import numpy as np

from .maps import SymMapper, detect_columns, check_columns, guess_kind
from .model import Molecule
from .session import parse_atom_spec

WEB = os.path.join(os.path.dirname(__file__), "web")
MAX_BYTES = 15 * 1024 * 1024    # artifact pages must stay at or below 16 MB


# ------------------------------------------------------------------ atoms
class Atoms:
    """Atom-spec lookup over the spec's models (molecule numbers start at 1, in spec order)."""

    def __init__(self, mols):
        self.mols = mols

    def resolve(self, spec, first_conformer=True):
        """Atom spec -> (molecule number, [atom indices]). A residue spec gives its heavy atoms."""
        d = parse_atom_spec(spec)
        k = d["mol"] or 1
        if not 1 <= k <= len(self.mols):
            raise ValueError(f"atom spec {spec!r}: no molecule {k} (there are {len(self.mols)})")
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
        return k, idx

    def point(self, spec):
        """xyz list, an atom spec (residue: centroid) or a list of specs (centroid of all)."""
        if isinstance(spec, (list, tuple)) and len(spec) == 3 and all(isinstance(v, (int, float)) for v in spec):
            return np.asarray(spec, float)
        if isinstance(spec, (list, tuple)):
            return np.mean([self.point(s) for s in spec], axis=0)
        k, idx = self.resolve(spec)
        return self.mols[k - 1].xyz[idx].mean(0)

    def xyz(self, specs):
        """All atom positions of a list of specs (residues give all their heavy atoms)."""
        out = []
        for s in specs:
            k, idx = self.resolve(s)
            out.append(self.mols[k - 1].xyz[idx])
        return np.concatenate(out) if out else np.zeros((0, 3))

    def name(self, k, i):
        a = self.mols[k - 1].info[i]
        spec = f"{a[0]}/{a[1]}{a[2]} {a[4]}" + (f":{a[5]}" if a[5] else "")
        return spec if len(self.mols) == 1 else f"{k}#{spec}"

    def nearest(self):
        """f(xyz) -> (spec, distance, position) of the nearest heavy atom over symmetry and lattice copies of each model;
        position is the image of xyz next to the model."""
        mappers = [(k, SymMapper(m.st, m.xyz, np.array([e != "H" for e in m.elements])))
                   for k, m in enumerate(self.mols, 1)]

        def f(xyz):
            best = min(((k, *sm.nearest(xyz)) for k, sm in mappers), key=lambda r: r[2])
            k, pos, d, i = best
            return self.name(k, i), float(d), np.asarray(pos)
        return f

    @staticmethod
    def residue_list(m, limit=40):
        seen = []
        for key, inf in zip(m.keys, m.info):
            r = f"{key[0]}/{key[1]}{key[2]} {inf[3]}"
            if r not in seen:
                seen.append(r)
        return ", ".join(seen[:limit]) + (" ..." if len(seen) > limit else "")


# ------------------------------------------------------------------- maps
def read_map(spec, base):
    """Map spec -> (gemmi.FloatGrid, kind). CCP4/MRC maps, or MTZ with columns 'F,PHI' (default: first detected)."""
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


def sample_box(grid, lo, hi, step):
    """Resample a (periodic) map on an orthogonal grid from lo to hi (Å). Returns (array[x, y, z], n)."""
    n = np.maximum(2, np.ceil((hi - lo) / step).astype(int) + 1)
    arr = np.zeros(tuple(n), dtype=np.float32)
    tr = gemmi.Transform()
    tr.mat.fromlist([[step, 0, 0], [0, step, 0], [0, 0, step]])
    tr.vec.fromlist([float(v) for v in lo])
    grid.interpolate_values(arr, tr)
    return arr, n


def grid_spacing(grid):
    uc = grid.unit_cell
    return min(uc.a / grid.nu, uc.b / grid.nv, uc.c / grid.nw)


class BoxMap:
    """A map resampled on the region box, with the level convention of the page."""

    def __init__(self, spec, base, lo, hi, step=None):
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
        self.step = float(step or spec.get("step") or grid_spacing(grid))
        self.lo = np.asarray(lo, float)
        self.arr, self.n = sample_box(grid, self.lo, np.asarray(hi, float), self.step)
        self.source = os.path.basename(spec["path"])

    def to_level_units(self, v):
        """Map value -> the page's level units (σ of the whole map, or the map's own units)."""
        return (v - self.mean) / self.rms if self.units == "sigma" else v

    def value_at(self, xyz):
        """Trilinear interpolation on the box (same as the page)."""
        f = (np.asarray(xyz, float) - self.lo) / self.step
        i0 = np.clip(np.floor(f).astype(int), 0, self.n - 2)
        t = np.clip(f - i0, 0, 1)
        v = 0.0
        for c in range(8):
            d = np.array([c & 1, (c >> 1) & 1, (c >> 2) & 1])
            w = np.prod(np.where(d, t, 1 - t))
            v += w * self.arr[tuple(i0 + d)]
        return float(v)

    def inside(self, xyz):
        f = (np.asarray(xyz, float) - self.lo) / self.step
        return bool(np.all(f >= 0) and np.all(f <= self.n - 1))

    def peaks(self, level, negative=False, center=None, radius=None):
        """Local extrema beyond `level` (page units) in the box, optionally within radius of center: [(xyz, value)]."""
        from scipy.ndimage import maximum_filter
        a = self.to_level_units(self.arr)
        a = -a if negative else a
        mx = maximum_filter(a, size=3, mode="nearest")
        idx = np.argwhere((a == mx) & (a > level))
        xyz = self.lo + idx * self.step
        val = a[tuple(idx.T)]
        if center is not None:
            keep = np.linalg.norm(xyz - center, axis=1) <= radius
            xyz, val = xyz[keep], val[keep]
        order = np.argsort(-val)
        return [(xyz[j], float(-val[j] if negative else val[j])) for j in order]

    def to_json(self):
        scale = max(float(np.abs(self.arr).max()), 1e-6) / 32000.0
        q = np.round(self.arr / scale).astype("<i2")
        data = base64.b64encode(gzip.compress(q.tobytes(), compresslevel=9, mtime=0)).decode()
        return {"name": self.name, "kind": self.kind, "units": self.units, "unit_label": self.unit_label,
                "level": self.level, "visible": self.visible, "color": self.color, "mean": self.mean, "rms": self.rms,
                "n": [int(v) for v in self.n], "origin": [round(float(v), 4) for v in self.lo], "step": self.step,
                "scale": scale, "source": self.source, "data": data}


# ------------------------------------------------------------- orientation
def orientation(xyz, mode="face"):
    """Quaternion [x, y, z, w] (three.js order) that lays the atoms' longest axis horizontally on screen.

    face: their plane faces the viewer. edge: their plane is seen edge-on (shortest axis vertical on screen),
    which shows out-of-plane motion."""
    c = xyz - xyz.mean(0)
    if len(c) < 2:
        return None
    w, v = np.linalg.eigh(c.T @ c)          # ascending
    e1, e2, e3 = v[:, 2], v[:, 1], v[:, 0]  # longest, middle, shortest
    if mode == "edge":
        e2, e3 = e3, e2
    R = np.array([e1, e2, np.cross(e1, e2)])  # rows: screen x, y, z in model coordinates (right-handed)
    return _quat(R)


def _quat(R):
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


# ------------------------------------------------------------------ build
class Share:
    def __init__(self, spec, base):
        self.spec, self.base = spec, base
        self.mols = []
        for i, ms in enumerate(spec.get("models") or []):
            m = Molecule(_path(ms["path"], base), name=ms.get("name"))
            m.color, m.hydrogens = ms.get("color"), ms.get("hydrogens", False)
            m.style = ms.get("style")
            m.visible = ms.get("visible", True)
            self.mols.append(m)
        if not self.mols:
            raise ValueError("the spec needs at least one model")
        self.atoms = Atoms(self.mols)
        lo, hi = self.region()
        self.maps = [BoxMap(s, base, lo, hi, s.get("step")) for s in spec.get("maps") or []]
        names = [x.name for x in self.maps]
        if len(set(names)) != len(names):
            raise ValueError(f"map names must be unique: {names}")
        self.stops = [self.stop(s, k) for k, s in enumerate(spec.get("stops") or [])]

    def region(self):
        r = self.spec.get("region") or {}
        if "center" in r:
            c, rad = self.atoms.point(r["center"]), float(r.get("radius", 12.0))
            return c - rad, c + rad
        xyz = np.concatenate([m.xyz for m in self.mols])
        margin = float(r.get("margin", 5.0))
        return xyz.min(0) - margin, xyz.max(0) + margin

    def stop(self, s, k):
        where = f"stop {k + 1} ({s.get('name', '')!r})"
        try:
            center = self.atoms.point(s["center"])
            frame = s.get("frame")
            if frame is None and not (isinstance(s["center"], (list, tuple)) and len(s["center"]) == 3
                                      and all(isinstance(v, (int, float)) for v in s["center"])):
                frame = s["center"] if isinstance(s["center"], (list, tuple)) else [s["center"]]
            mode = s.get("orient", "face")
            if isinstance(mode, (list, tuple)):
                orient = [float(v) for v in mode]
            elif mode == "keep" or not frame:
                orient = None
            else:
                pts = self.atoms.xyz(frame)
                orient = orientation(pts, mode) if len(pts) >= 3 else None
            maps = {}
            for name, v in (s.get("maps") or {}).items():
                if name not in {x.name for x in self.maps}:
                    raise ValueError(f"no map named {name!r}; maps: {[x.name for x in self.maps]}")
                maps[name] = v
            labels = [[k2, i] for spec in s.get("labels") or [] for k2, idx in [self.atoms.resolve(spec)] for i in idx]
            dists = []
            for a, b in s.get("distances") or []:
                ka, ia = self.atoms.resolve(a)
                kb, ib = self.atoms.resolve(b)
                dists.append([[ka, ia[0]], [kb, ib[0]]])
        except (ValueError, KeyError) as e:
            raise ValueError(f"{where}: {e}") from None
        radius = float(s.get("radius", 8.0))
        return {"name": s.get("name") or f"stop {k + 1}", "text": s.get("text", ""),
                "center": [round(float(v), 3) for v in center], "orient": orient, "radius": radius,
                "slab": float(s.get("slab", radius)), "maps": maps, "labels": labels, "distances": dists,
                "_center_spec": s["center"]}

    def payload(self):
        mols = []
        for k, m in enumerate(self.mols, 1):
            j = m.to_json()
            mols.append({"id": k, "name": m.name, "xyz": j["xyz"], "element": j["element"], "atoms": j["atoms"],
                         "bonds": j["bonds"], "orders": j["orders"], "visible": m.visible, "color": m.color,
                         "hydrogens": bool(m.hydrogens), "style": m.style or "sticks"})
        stops = [{k: v for k, v in s.items() if not k.startswith("_")} for s in self.stops]
        if stops:
            start = stops[0]
        else:
            c = np.concatenate([m.xyz for m in self.mols]).mean(0)
            start = {"name": "", "text": "", "center": [round(float(v), 3) for v in c], "orient": None,
                     "radius": 12.0, "slab": 12.0, "maps": {}, "labels": [], "distances": []}
        return {"title": self.spec.get("title", "Moorlight"), "description": self.spec.get("description", ""),
                "contour_radius": float(self.spec.get("contour_radius", 10.0)), "models": mols,
                "maps": [x.to_json() for x in self.maps], "stops": stops, "start": start}

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
    def report(self, at=(), n_peaks=12):
        out = []
        p = out.append
        near = self.atoms.nearest()
        p(f"models: " + "; ".join(f"{k}# {m.name}: {len(m.keys)} atoms ({Atoms.residue_list(m, 12)})"
                                  for k, m in enumerate(self.mols, 1)))
        for x in self.maps:
            u = x.unit_label or x.units
            box = x.to_level_units(x.arr)
            p(f"\nmap {x.name!r} ({x.source}, {x.kind}): level {x.level} {u}; whole map mean {x.mean:.3f} rms {x.rms:.3f}; "
              f"in the region min {box.min():.2f} max {box.max():.2f} {u}; box {x.n.tolist()} x {x.step:.2f} A")
            lv = x.level
            signs = [False, True] if x.kind == "diff" else [False]
            for neg in signs:
                pk = x.peaks(lv, negative=neg)
                listed = []          # symmetry/lattice copies of a peak are listed once, at the image next to the model
                for xyz, v in pk:
                    nm, d, pos = near(xyz)
                    if any(np.linalg.norm(pos - q) < 0.5 for q, *_ in listed):
                        continue
                    listed.append((pos, v, nm, d))
                    if len(listed) == n_peaks:
                        break
                p(f"  {'negative' if neg else 'positive'} peaks beyond {'-' if neg else '+'}{lv} {u}: {len(pk)} in the region box"
                  + (f"; strongest {len(listed)} unique:" if listed else ""))
                for pos, v, nm, d in listed:
                    p(f"    {v:+7.2f} {u} at ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})  {d:.2f} A from {nm}")
        for s in self.stops:
            c = np.array(s["center"])
            p(f"\nstop {s['name']!r}: centre {s['_center_spec']!r} = ({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f}), radius {s['radius']}, "
              f"orient {'set' if s['orient'] else 'kept from the previous stop'}")
            for x in self.maps:
                if not x.inside(c):
                    p(f"  WARNING: centre outside the box of map {x.name!r}")
            vis = [x for x in self.maps if s["maps"].get(x.name, x.visible) not in (False, None)]
            for k, i in [tuple(l) for l in s["labels"]]:
                xyz = self.mols[k - 1].xyz[i]
                vals = ", ".join(f"{x.name}: {x.to_level_units(x.value_at(xyz)):+.2f}" for x in vis)
                p(f"  label {self.atoms.name(k, i)}: {vals}")
            for x in vis:
                if x.kind != "diff":
                    continue
                u = x.unit_label or x.units
                lv = s["maps"].get(x.name)
                lv = float(lv) if isinstance(lv, (int, float)) and not isinstance(lv, bool) else x.level
                for neg in (False, True):
                    pk = x.peaks(lv, negative=neg, center=c, radius=s["radius"])
                    if pk:
                        p(f"  {x.name}: {len(pk)} {'negative' if neg else 'positive'} peaks beyond {lv} {u} within "
                          f"{s['radius']} A; strongest: " + "; ".join(
                              f"{v:+.1f} at {near(xyz)[1]:.2f} A from {near(xyz)[0]}" for xyz, v in pk[:6]))
        for spec in at:
            k, idx = self.atoms.resolve(spec, first_conformer=False)
            p(f"\nvalues at {spec!r}:")
            for i in idx:
                xyz = self.mols[k - 1].xyz[i]
                vals = ", ".join(f"{x.name}: {x.to_level_units(x.value_at(xyz)):+.2f}"
                                 + ("" if x.inside(xyz) else " (outside box)") for x in self.maps)
                p(f"  {self.atoms.name(k, i):<18} {vals}")
        return "\n".join(out)


def _path(p, base):
    p = os.path.expanduser(p)
    return p if os.path.isabs(p) else os.path.join(base, p)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="moorlight share", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec", help="share spec (JSON)")
    ap.add_argument("-o", "--out", help="output page (.html); default: next to the spec")
    ap.add_argument("--check", action="store_true", help="print peaks and values for writing the stop texts; no page")
    ap.add_argument("--at", action="append", default=[], metavar="SPEC",
                    help="with --check: map values at these atoms or residues (repeatable)")
    ap.add_argument("--peaks", type=int, default=12, help="with --check: peaks listed per map and sign")
    ap.add_argument("--standalone", action="store_true",
                    help="full HTML document for opening as a local file (artifact pages must not have it)")
    a = ap.parse_args(argv)
    with open(a.spec) as fh:
        spec = json.load(fh)
    try:
        sh = Share(spec, os.path.dirname(os.path.abspath(a.spec)))
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    if a.check:
        print(sh.report(a.at, a.peaks))
        return 0
    out = a.out or os.path.splitext(a.spec)[0] + ".html"
    page = sh.html(standalone=a.standalone)
    with open(out, "w") as fh:
        fh.write(page)
    size = len(page.encode())
    print(f"wrote {out}: {size / 1e6:.2f} MB, {len(sh.mols)} model(s), {len(sh.maps)} map(s), {len(sh.stops)} stop(s)")
    if size > MAX_BYTES:
        print(f"warning: larger than {MAX_BYTES / 1e6:.0f} MB; an artifact page must stay at or below 16 MB. "
              "Use fewer maps, a smaller region or a coarser step.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
