---
name: moorlight-share
description: Make a self-contained, phone-friendly Moorlight viewer page (models + maps + preset slice planes) that needs no server, for sending as a file or publishing as a claude.ai artifact. Use when the user wants to look at a structure and its maps from a phone, a Remote Control session or a link, or asks for preset planes through a structure.
---

# Shareable Moorlight viewer pages

`moorlight share` turns a small JSON spec into one HTML page holding:
- the models, plus the atoms of neighbouring cells and symmetry copies around them;
- the maps;
- a starting view;
- preset slice planes.

The browser draws everything itself:
- contours as marching-cubes mesh (Coot's chicken wire) or as shaded surfaces;
- Olex2-style map slices.

On the page you can:
- rotate, zoom and pan;
- switch maps and set their levels with sliders;
- tap atoms for labels and map values;
- hide the structure;
- move, tilt, thicken, fix, save and copy slice planes.

It is a viewer, not a presentation: there are no stops or texts.

## Workflow

1. **Fast start**: `moorlight share init MODEL "MAP=name" ["MAP2=name" ...] -o spec.json [--title "Cage 3.5 ps"]
   [--units absolute --unit-label SD] [--level 4] [--min 10] [--planes 5] [--build page.html]`.
   - Takes about 10 s and writes `spec.json` plus `spec.check.txt` (the full `--check` output).
   - The spec holds a starting view of the whole structure and one preset plane per group of strong features of
     the first map (peaks beyond `--min`, default 2.5 × level, within 1.5 Å of an atom):
     - a large residue gets a plane through its atoms around the strongest feature;
     - a metal gets the plane of its coordination;
     - waters and other small residues are grouped by name.
2. **Edit the spec** (format below). When the user asks for a plane ("a plane through the azo N and the ring
   carbons", "perpendicular to that, through Pd1", "a stack 1 Å apart"), add it to `planes`, check it, build, and
   republish to the same artifact.
3. **Check**: `moorlight share spec.json --check [--at SPEC ...] [--near SPEC --radius R [--frame "A, B, C"]]`.
   It prints:
   - each map's strongest peaks of each sign, each listed once next to the model, with the nearest atom;
   - for each plane:
     - its centre and normal;
     - how flat its atoms are (rms);
     - the atoms within the slice, with their heights `h`;
     - the peaks in the slice.
   - with `--at`, the values at atoms;
   - with `--near`, the atoms and peaks around a point.

   Options: `--map NAME` lists only that map; `--peaks N` and `--stop-peaks N` set the list lengths.
4. **Build**: `moorlight share spec.json -o page.html`. The command prints the size. An artifact page must stay at
   or below 16 MB; a small-molecule cell map takes about 0.4 MB.
5. **Hand it over**:
   - Publish it as an artifact (Artifact tool: `file_path` = the page, `icon` "molecule", a one-sentence
     `description`). The page is in the artifact form, without `<!doctype>`/`<html>`/`<head>`/`<body>`.
   - Republishing the same file path keeps the link.
   - Or build with `--standalone` for a file someone opens locally.

   The page cannot be rendered on a cluster node without a browser. Say so if you could not look at it.

## What the numbers mean

- **Maps are not resampled.** A map whose cell grid has up to 4 million points is stored on its own grid over the
  unit cell and wraps periodically. Values and peak heights match Coot. Larger maps are resampled on a box around
  the models (`"grid": "box"`, tricubic).
- **Peaks** are grid points that are the extreme of their 3×3×3 neighbourhood. One finder serves every list.
- **Values at atoms** are trilinear interpolation on the stored grid, the same number a tap shows.
- **`h`** is the signed distance from a plane (plane lists) or from the frame atoms' plane (`--near --frame`).
- **Neighbour copies** are drawn in grey for residues within `neighbour_margin` Å (default 5) of the model. They
  are named with a tag, `S/14 O [-a]`, and specs accept the same text. Set `"neighbours": false` on a model to
  leave them out.

## Spec

```json
{
  "title": "Cage 3.5 ps",
  "description": "Light minus dark at 3.5 ps on the dark structure; contours in noise SD.",
  "models": [{"path": "dark.pdb", "name": "dark"}],
  "maps": [{"path": "3p5ps_corrected.ccp4", "name": "3.5 ps light - dark", "kind": "diff",
            "units": "absolute", "unit_label": "SD", "level": 4.0, "visible": true}],
  "view": {"center": ["A/1", "A/2", "A/3"], "radius": 13, "labels": ["-/102 N9"]},
  "planes": [
    {"name": "Guest 102 azo", "atoms": ["-/102 C8", "-/102 N9", "-/102 N10", "-/102 C11"],
     "map": "3.5 ps light - dark", "level": 4, "thickness": 0.5, "radius": 6, "copies": [-1, 1]},
    {"name": "Across the azo bond", "center": "-/102 N9", "normal": ["-/102 N9", "-/102 N10"]}
  ]
}
```

- **title**: the page's name, two to four words ("Cage 3.5 ps"). Never a name followed by a dash or colon and an
  explanation; the explanation goes in `description`.
- **models**:
  - `path` (PDB/mmCIF) and `name`;
  - optional `color` (carbon colour, CSS);
  - `hydrogens` (default false);
  - `style` (`sticks` or `lines`);
  - `visible`;
  - `neighbours` (default true).
- **maps**:
  - `path`: a CCP4/MRC map, or an MTZ with `columns: "F,PHI"`.
  - `name`: unique; planes refer to it.
  - `kind`: `diff` (green + / red −) or `2fofc` (blue).
  - `units` and `level`, `visible`.
  - `grid`: `auto`, `native` or `box`; for `box`, also `step` and a top-level `region`.
  - `units: "sigma"` (default): levels are multiples of the whole map's rms.
  - `units: "absolute"`: levels are in the map's own units. Use this for maps already in noise SD, whose rms
    includes signal. `unit_label` is the text after the number ("SD").
- **contour_radius**: Å around the view centre that is contoured (default 10, Coot's map radius).
- **view** (optional; default the whole structure) — the starting view and what Reset view returns to:
  - `center`: atom spec, residue spec, a list of specs (centroid) or `[x, y, z]`.
  - `radius`: Å that fit the screen. `slab`: depth half-thickness.
  - `orient`: `face` (default; the `frame` atoms' plane faces the viewer, by default the centre's atoms), `edge`,
    or a quaternion `[x, y, z, w]`.
  - `maps`: `{name: level}` or `false`. `labels` and `distances`: atom specs.
- **planes** (preset slice planes; the viewer sees them in a list and can look at, hide, copy, edit or delete
  them):
  - Placement: `atoms` (three or more) gives their best plane through their centroid. `center` + `normal` gives a
    plane through that point; `normal` is a vector or two atom specs (from the first to the second). `center`
    may also be given with `atoms` to move the plane's centre.
  - `axis`: two atom specs that set the plane's horizontal direction.
  - `offset`: Å along the normal.
  - `copies`: extra parallel planes at these offsets from it (Å, e.g. `[-1, 1]`), named "<name> +1 Å".
  - `thickness`: Å. A thick slice shows the strongest value within it, with its sign kept.
  - `radius`: half the square's side (default `contour_radius`).
  - `map`: the map shown (default the first visible).
  - `level`: fixed contour level for this plane (default: follows the map's slider).
  - `color`: outline colour, CSS hex. `visible`: default true.

**Atom specs**:
- `A/85 CA`, `A/85 LYS CA:B` (altloc B), `2#A/85 SG` (model 2), `A/85` (residue: its heavy atoms).
- A neighbour copy adds its tag: `S/14 O [-a]`. Chain names are as in the file, e.g. `-`.
- An unknown atom gives an error that lists the atoms of that residue.

## The page

- **Maps:** visibility and a log-scale level slider for each map.
- **Maps drawn as:** Mesh, Surface or Slice. In slice mode:
  - the live plane (yellow outline) faces the viewer and turns with the molecule until it is fixed;
  - sliders for depth, thickness and two tilts;
  - buttons: Fix plane, Through labelled atoms (best plane through three or more tapped atoms), Face me, and
    Save plane (S).
- **Planes:** the preset planes and the viewer's own saved ones, each with its outline colour.
  - Tap a name to look at that plane, centred and face-on.
  - The checkbox hides it.
  - Copy adds a parallel plane at the copy offset (slider, −5 to +5 Å).
  - Edit makes it the live plane; × removes it.
  - The viewer's own planes and copies are kept in their browser. Presets come with the page.
- **Structure:** a checkbox for each model and for the neighbouring cells, plus hydrogens, clear labels and Reset
  view.
- **Depth bar** (right edge of the view): the structure from front (top) to back.
  - Its window is the visible slab around the focus. Drag the window to move the focus into or out of the screen.
  - Drag the window's ends to change the slab; tap the bar to jump.
  - The yellow tick is the live slice. Keys: `,` `.` move, `[` `]` change the slab.
- **Full screen** (button or F): level and depth sliders, Fix plane, Save, and an Atoms toggle (A).
- **Theme:** one deliberate dark look (Coot black).
- **Layout:** on an upright phone the panel sits under the view, otherwise beside it.
