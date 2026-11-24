---
name: moorlight-share
description: Make a self-contained, phone-friendly Moorlight page (models + maps + a short tour) that needs no server, for sending as a file or publishing as a claude.ai artifact. Use when the user wants to look at a structure and its maps from a phone, a Remote Control session or a link, or to share a guided view of a result.
---

# Shareable Moorlight snapshots

`moorlight share` turns a small JSON spec into one HTML page holding:
- the models, plus the atoms of neighbouring cells and symmetry copies around them;
- the maps;
- a list of tour stops.

The browser contours the maps itself, with chicken wire as in Coot. On the page you can:
- rotate, pinch-zoom and pan;
- switch maps and change their contour levels;
- tap an atom to get its label and the map values there;
- step through the stops.

Nothing runs on the cluster after the page is written.

## Workflow

0. **Fast start**: `moorlight share init MODEL "MAP=name" ["MAP2=name" ...] -o spec.json [--title "Cage 3.5 ps"]
   [--units absolute --unit-label SD] [--level 4] [--min 10] [--stops 5] [--build page.html]`.
   - Takes about 10 s and writes `spec.json` plus `spec.check.txt` (the full `--check` output).
   - The spec holds an overview and one stop per group of strong features of the first map: peaks beyond `--min`
     (default 2.5 × level) within 1.5 Å of an atom.
   - Residues with more than five heavy atoms get a stop each (guests and ligands, shown first, planar ones
     edge-on). Smaller residues (metals, waters, counter-ions) are grouped by name and shown at the strongest one.
   - Each stop has its view, level and labels, and a draft text that states the numbers.
   - Usually you only reword the texts, add what the user's notes say about reliability, drop or reorder stops, and
     build. `--build` gives a page right away; the drafts are factual and fine for a first look.
1. **Write or edit the spec** (`spec.json`, format below). Paths are relative to the spec file.
2. **Get the numbers**: `moorlight share spec.json --check [--at SPEC ...] [--near SPEC --radius R ...]`.
   It prints:
   - for each map, its strongest peaks of each sign, each listed once next to the model, with the nearest atom;
   - for each stop:
     - the values of every map at its labels;
     - the peaks of every difference map within its radius, with coordinates, the nearest drawn atom and the
       height above the frame plane (`h`);
   - for `--at`, the values at atoms or residues;
   - for `--near`, the atoms and the peaks around an atom, residue or `"x, y, z"`; add
     `--frame "A/1 N1, A/1 N2, A/1 N3"` for heights above that plane. Use this to probe without adding stops.

   Options: `--map NAME` lists only that map; `--peaks N` sets the whole-map lists and `--stop-peaks N` the
   per-stop lists (default 8). A label in another residue also gets its `h`, which tells which side of the
   frame plane that residue lies on.

   Every statement in a stop text should come from these numbers.
3. **Build**: `moorlight share spec.json -o page.html`. The command prints the size. An artifact page must stay at
   or below 16 MB. A small-molecule cell map takes about 0.4 MB.
4. **Hand it over**:
   - Publish it as an artifact (Artifact tool: `file_path` = the page, `icon` "molecule", a one-sentence
     `description`). The page is written in the artifact form, without `<!doctype>`/`<html>`/`<head>`/`<body>`.
   - Or build with `--standalone` for a file someone opens locally.

   The page cannot be rendered on a cluster node without a browser. Say so if you could not look at it.

## What the numbers mean

- **Maps are not resampled.** A map whose cell grid has up to 4 million points is stored on its own grid over the
  whole unit cell, and the page wraps it periodically. Values and peak heights are then the same as in Coot.
  Larger maps are resampled on a box around the models (`"grid": "box"`; tricubic, step half the grid spacing),
  and `--check` says so.
- **Peaks** are grid points that are the extreme of their 3×3×3 neighbourhood. One finder serves every list, so
  the same peak has the same value everywhere. Within a stop, each peak is shown at its lattice image nearest the
  stop centre.
- **Values at atoms** are trilinear interpolation on the stored grid. A tap on the page shows the same number.
- **`h`** is the signed distance from the best plane through the stop's frame atoms. With `orient: face` it is
  positive towards the viewer; with `edge`, positive is up on the screen. Use it for statements like "0.7 Å out of
  the guest plane".
- **Neighbouring cells**: by default the page also draws, in grey lines, the copies of every residue that come
  within `neighbour_margin` Å (default 5) of the model, as Coot draws symmetry atoms. `--check` names them with
  a tag, `S/14 O [-a]` (one lattice vector a back) or `[x,-y,z+1/2 +b]` (operator, then shift). Specs accept the
  same text, so a stop can label or centre on a copy. Set `"neighbours": false` on a model to leave them out.

## Spec

```json
{
  "title": "Cage 3.5 ps",
  "description": "Light minus dark at 3.5 ps on the dark structure; contours in noise SD.",
  "models": [{"path": "dark.pdb", "name": "dark", "hydrogens": false}],
  "maps": [{"path": "3p5ps_corrected.ccp4", "name": "3.5 ps light - dark", "kind": "diff",
            "units": "absolute", "unit_label": "SD", "level": 4.0, "visible": true}],
  "contour_radius": 10,
  "stops": [
    {"name": "Azo bond, guest 102", "text": "Two to four sentences, with numbers.",
     "center": "-/102", "radius": 7, "orient": "edge", "frame": ["-/102 C8", "-/102 N9", "-/102 N10", "-/102 C11"],
     "maps": {"3.5 ps light - dark": 4.0}, "labels": ["-/102 N9", "-/102 N10", "S/14 O [-a]"],
     "distances": [["-/102 N9", "-/102 N10"]]}
  ]
}
```

- **title**: the page's name, two to four words ("Cage 3.5 ps"). Never a name followed by a dash or colon and an
  explanation; the explanation goes in `description`.
- **models**:
  - `path` (PDB/mmCIF) and `name`;
  - optional `color` (carbon colour, CSS);
  - `hydrogens` (default false);
  - `style` (`sticks` default, or `lines`);
  - `visible`;
  - `neighbours` (default true).
- **maps**:
  - `path`: a CCP4/MRC map, or an MTZ with `columns: "F,PHI"`.
  - `name`: unique; stops refer to it.
  - `kind`: `diff` (green + / red − contours) or `2fofc` (one blue contour).
  - `units` and `level`, `visible`.
  - `grid`: `auto` (default), `native` or `box`. For `box`, also `step` and the top-level `region`.
  - `units: "sigma"` (default): levels are multiples of the whole map's rms, as in Coot's default.
  - `units: "absolute"`: levels are in the map's own units. Use this for maps already scaled to noise SD. Their
    rms includes signal, so σ contours would mislead. `unit_label` is the text after the number ("SD").
- **contour_radius**: how far from the view centre maps are contoured, in Å (default 10, as Coot's map radius). A
  stop may set its own.
- **neighbour_margin**: Å around the model within which neighbour copies are drawn (default 5).
- **region** (only for `grid: box`): `{"margin": 5}` around all models, or `{"center": SPEC, "radius": 12}` for a
  cube.
- **stops** (optional; without stops the page opens on the whole model):
  - `center`: atom spec, residue spec (centroid), a list of specs (centroid) or `[x, y, z]`.
  - `radius`: Å of context that fits the screen (default 8). `slab`: depth half-thickness (default `radius`).
  - `orient`: `face` (default; the frame atoms' plane faces the viewer, longest axis horizontal), `edge` (that
    plane seen edge-on, its normal pointing up, which shows out-of-plane motion), `keep` (the previous stop's
    rotation) or a quaternion `[x, y, z, w]`. `frame` lists the atoms that define the plane; the default is the
    centre's atoms. A lone atom (a metal) uses itself and the atoms within 2.6 Å (its coordination). Otherwise
    fewer than three atoms keep the rotation.
  - `maps`: `{map name: level}` shows that map at that level, `false` hides it. Maps not named go back to their
    defaults from `maps`, so every stop looks the same however the viewer got there.
  - `labels`: atom specs to label. `distances`: pairs of atom specs.
  - `contour_radius`: overrides the global one for this stop.
  - `text`: plain text. A blank line starts a new paragraph.

**Atom specs**:
- `A/85 CA`, `A/85 LYS CA:B` (altloc B), `2#A/85 SG` (model 2 in spec order), `A/85` (residue: its heavy atoms,
  or their centroid as a centre).
- A neighbour copy adds its tag: `S/14 O [-a]`.
- Chain names are as in the file, for example `-`.
- Without an altloc, an atom spec takes the first conformer; `--at` lists every conformer.
- An unknown atom gives an error that lists the atoms of that residue. An unknown copy lists the copies drawn.
- On the command line, `--at -/102` and `--near -/102` work as written.

## The page

- **Theme:** one deliberate dark look (Coot black) in every theme.
- **Layout:**
  - On an upright phone the panel sits under the view.
  - On a phone turned sideways, and on wide screens, the panel sits beside the view.
  - The button at the top right of the view (or F) shows the view alone, with a small floating ◀ ▶ bar for the
    stops. The page remembers that choice.
- **Map styles:**
  - Mesh: marching cubes, Coot's chicken wire.
  - Surface: shaded and transparent.
  - Slice: an Olex2-style section through the focus that faces the viewer and turns with the molecule. It shows
    the first visible map in green/red, with contour lines at 1, 2, 4, 8 … × level, and has a depth slider.
- **Levels:** every map has a log-scale level slider.

## Writing the tour

- Start with an overview stop (the whole molecule, `radius` large enough), then one stop per finding.
- One finding per stop, two to four sentences. Name the map and the level, and give the numbers from `--check`
  ("−14.0 SD on N9; +12.2 SD 0.7 Å from it, h −1.0 Å").
- The page is read on a phone; keep it short.
- Show each map only where it is needed and say which one is on.
- Report what the map shows. Interpretation that the user has not established does not belong in the texts.
- Check every stop in the `--check` output: are the centre and labels right, and is the strongest feature where
  the text says?
