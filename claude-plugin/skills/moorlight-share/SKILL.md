---
name: moorlight-share
description: Make a self-contained, phone-friendly viewer page for a structure and its maps (models + density/difference maps + preset slice planes) that needs no server, for publishing as a claude.ai artifact or sending as a file. Works for small molecules, proteins and EM maps; inputs PDB/mmCIF and CCP4/MRC maps or MTZ map coefficients. Use when the user wants to look at a structure or map from a phone, a Remote Control session or a link, or asks for planes or slices through a structure.
---

# Shareable structure viewer pages

`moorlight-share` (the same as `moorlight share`) turns model and map files into one HTML page. The page holds:
- the models, plus the atoms of neighbouring cells and symmetry copies around them (Coot's symmetry atoms);
- the maps;
- a starting view;
- preset slice planes.

The browser draws everything itself:
- contours, as marching-cubes mesh (Coot's chicken wire);
- Olex2-style map slices.

On the page the viewer can rotate, zoom and pan, switch maps and set their levels with sliders, and tap atoms for
labels and map values. They can hide the structure, and move, tilt, thicken, fix, save and copy slice planes. It is a
viewer, not a presentation: there are no stops or texts.

**Running it.** Use `moorlight-share ...`. If that is not on PATH, use `moorlight share ...` or
`python -m moorlight share ...` from an environment where the `moorlight` package is installed (it needs gemmi, numpy
and scipy). There is no server and no port.

## Workflow

0. **Look up the residues**: `moorlight-share list MODEL` prints, as atom specs:
   - the chains with their residue ranges;
   - the ligands, grouped by name;
   - metals and ions, with the number of atoms around each;
   - the water count;
   - the residues with alternate conformers.

   Use it instead of grepping the file.
1. **Draft**: `moorlight-share init MODEL MAPS... -o spec.json [--title "Trypsin benzamidine"] [--description "..."] [--plane NAME=SPECS ...] --build page.html`
   - Takes about 10 s. It writes `spec.json`, the page, and `spec.check.txt`, the full `--check` report.
   - **Models:** `.pdb`, `.cif`, `.ent`, also gzipped. `model.pdb=name` names one (default: the file name).
   - **Maps:** `.ccp4`, `.map` or `.mrc`, or an `.mtz` of map coefficients.
     - An MTZ gives every map a refinement program writes for viewing: 2mFo-DFc and mFo-DFc (FWT/DELFWT,
       2FOFCWT/FOFCWT, ...).
     - TorchRef difference refinement (`difference_data.mtz`): the light − dark map dF on the dark phases,
       q-weighted (`dF*W_Q,PHDELWT`; also W_InVa-weighted and unweighted), and the extrapolated light map
       (FWT/PHWT there is 2Fext − Fc, not 2mFo-DFc). Its `*_dark-sf.mtz` gives the dark 2mFo-DFc and mFo-DFc.
     - `file.mtz=name:F,PHI` picks two columns; `F*W` weights the amplitudes by a weight column
       (`file.mtz=name:dF*W_Q,PHDELWT`). `file.ccp4=name` names a map.
   - **Map kind** (`2fofc` = density in blue, `diff` = green +/red −) is taken from the MTZ columns or the file name.
     Otherwise it is guessed from the values: a density map has almost nothing below −2.5 rms. `init` prints the
     kinds; correct a wrong one in the spec.
   - **Levels:** density 1.5 σ, difference 3 σ. `--level` sets the difference maps' level.
     - Maps already in noise units: `--units absolute --unit-label SD`.
   - **Planes:**
     - `--plane "Benzamidine=A/1"` is the best plane through a residue's heavy atoms.
     - `--plane "Ring=A/401 C1, A/401 N1, A/401 C2, A/401 C3"` goes through these atoms.
     - `--plane "Zn site=A/301 ZN"` (one atom) is that atom's coordination plane.
     - The option is repeatable.
     - In the page every plane shows the map in use (the first visible one) at its slider level, so switching
       maps or levels changes all planes. A plane's `map` and `level` (`init` sets the first difference map) are
       only its start: they apply while that map is the one in use and its level has not been moved.
     - Without `--plane`, `init` adds up to `--planes` (default 5) automatic planes:
       - at the strongest features of the first difference map (peaks beyond `--min`, default 1.5 × its level,
         within 1.5 Å of an atom); a residue present in several models (dark and light) gives one plane;
       - with no difference map or no such features: through the ligands, the most buried first (a bound
         chromophore or inhibitor before surface lipids and additives), and the coordination of metal ions.
     - With `--planes 0` and no `--plane`, there are none; add them to the spec later.
2. **Edit the spec** when the user asks for something else: another plane, a different centre, colour or level.
   - When they ask for a plane ("through the azo N and the ring carbons", "perpendicular to that, through Pd1",
     "a stack 1 Å apart"), add it to `planes`, check, build, and republish to the same artifact.
3. **Check**: `moorlight-share spec.json --check [--at "SPEC, SPEC"] [--near SPEC --radius R [--frame "A, B, C"]]`.
   - It prints:
     - each map's grid and strongest peaks of each sign, each listed once next to the model with its nearest atom;
     - for each plane:
       - its centre and normal;
       - its atoms' rms about their best plane, and how far the plane was moved to pass through a given `center`;
       - a metal's coordination with distances;
       - the atoms in the slice with their heights `h`;
       - the peaks in the slice.

       Three atoms always give rms 0. Use four or more to judge flatness.
     - `--at`: values of every map at these atoms or residues (comma-separated, or repeat `--at`).
     - `--near`: the atoms around a point, one line per residue, and every map's peaks there. With
       `--frame "A, B, C"`, it also prints each atom's and peak's height `h` above those atoms' plane.
   - `--map NAME` limits it to one map.
   - `spec.check.txt` from `init` describes the spec as `init` wrote it; after edits, run `--check` again.
   - **Read the warnings.** A map that covers only part of the cell shows its mean outside its extent.
4. **Build**: `moorlight-share spec.json -o page.html`. It prints the size and each map's grid. An artifact page must
   stay at or below 16 MB. Typical sizes:

   | Structure | Page size |
   |---|---|
   | Small-molecule cell | 0.4 MB per map |
   | Protein, 2 maps | 3–5 MB |
   | Atomic-resolution protein, 2 maps | about 8 MB |

5. **Hand it over**:
   - Publish it as an artifact: Artifact tool with `file_path` = the page, `icon` "molecule", a one-sentence
     `description`. The page is in the artifact form, without `<!doctype>`/`<html>`/`<head>`/`<body>`.
   - Republishing the same file path keeps the link.
   - Or build with `--standalone` for a file to send (e-mail, messenger, a paper's supplement or a Zenodo
     deposit). The page fetches nothing: three.js (0.170.0, MIT) is written into it and it uses system fonts, so it
     opens offline in any current browser with WebGL 2.
     - Desktop and Android browsers run it from the file. On an iPhone, a file preview (Mail, WhatsApp, Files) may
       not run it; send iPhone users a link instead (the artifact, or the file on a web server).
     - Its "About this file" section lists the input files with SHA-256 checksums, the Moorlight version and
       commit, the build time and the spec (file names only).
   - A machine without a browser cannot render the page. Say so if you could not look at it.

## What the numbers mean

- **Maps are not resampled.** Each map is stored on its own grid, so values and peak heights match Coot. The page
  keeps them to 1/100 of the map's rms.
  - `native`: the whole unit cell, periodic. Used for small cells.
  - `block`: the part of the grid around the models (5 Å margin). Used when that is less than half the cell, as for
    large cells and EM boxes.
  - `box`: tricubic resampling. Used only if a block would not fit the page: about 8 million grid points for all
    maps together.
  - Set `"grid"` on a map to force one; `box` also takes a `step`.
- **Peaks** are grid points that are the extreme of their 3×3×3 neighbourhood. One finder serves every list.
- **Values at atoms** are trilinear interpolation on the stored grid: the same number a tap on the page shows.
- **`h`** is the signed distance from a plane (plane lists) or from the frame atoms' plane (`--near --frame`).
- **σ** (`units: "sigma"`) is the rms of the whole map. For a light − dark map whose rms includes signal, use the
  map's own noise units (`"units": "absolute"`) if it has them.
- **Neighbour copies** are drawn in grey for residues within 5 Å (`neighbour_margin`) of the model, for models with a
  crystal cell (not for EM models).
  - They are named with their symmetry operator and lattice shift in brackets: `A/219 O [-x,y+1/2,-z+1/2 -b]`, or
    `S/14 O [-a]` in P1. Specs accept the same text.
  - Set `"neighbours": false` on a model to leave them out.

## Spec

```json
{
  "title": "Trypsin benzamidine",
  "description": "PDB-REDO 3PTB: 2mFo-DFc at 1.5 σ, mFo-DFc at ±3 σ.",
  "models": [{"path": "3ptb_final.pdb", "name": "3ptb"}],
  "maps": [
    {"path": "3ptb_final.mtz", "columns": "FWT,PHWT", "name": "2mFo-DFc", "kind": "2fofc", "level": 1.5},
    {"path": "3ptb_final.mtz", "columns": "DELFWT,PHDELWT", "name": "mFo-DFc", "kind": "diff", "level": 3}
  ],
  "view": {"center": "A/1", "radius": 12, "labels": ["A/195 OG"]},
  "planes": [
    {"name": "Benzamidine", "atoms": ["A/1"], "map": "mFo-DFc", "copies": [-1, 1]},
    {"name": "Ca site", "coordination": "A/480 CA", "map": "2mFo-DFc"},
    {"name": "Across Ser195", "center": "A/195 OG", "normal": ["A/195 CB", "A/195 OG"], "thickness": 1}
  ]
}
```

- **title**: the page's name, two to four words. Never a name followed by a dash or colon and an explanation; the
  explanation goes in `description`.
- **models**:
  - `path` (PDB/mmCIF) and `name`;
  - optional `color` (carbon colour, CSS);
  - `hydrogens` (default false);
  - `style` (`sticks` or `lines`);
  - `visible`;
  - `neighbours` (default true).
- **maps**:
  - `path`: a CCP4/MRC map, or an MTZ with `columns: "F,PHI"` (or `"F*W,PHI"`: amplitudes times a weight column).
  - `name`: unique; planes refer to it.
  - `kind`: `diff`, or `2fofc` (also `density`).
  - `units`, `level`, `visible`, `color`.
  - `grid`: `auto`, `native`, `block` or `box`; for `box`, also `step` and a top-level `region`.
  - `units: "sigma"` (default): levels are multiples of the whole map's rms.
  - `units: "absolute"`: levels are in the map's own units (e.g. maps already in noise SD). `unit_label` is the text
    after the number ("SD").
  - `sample_rate`: MTZ grid spacing d_min / rate (default 3).
- **contour_radius**: Å around the view centre that is contoured (default 10, Coot's map radius).
- **view** (optional; default the whole structure) — the starting view and what Reset view returns to:
  - `center`: atom spec, residue spec, a list of specs (centroid) or `[x, y, z]`.
  - `radius`: Å that fit the screen. `slab`: depth half-thickness.
  - `orient`: `face` (default; the `frame` atoms' plane faces the viewer, by default the centre's atoms), `edge`,
    or a quaternion `[x, y, z, w]`.
  - `maps`: `{name: level}` or `false`. `labels` and `distances`: atom specs.
- **planes** (the viewer sees them in a list and can look at, hide, copy, edit or delete them):
  - **Placement**, one of:
    - `atoms` (three or more, or residues) gives their best plane through their centroid;
    - `coordination` (one atom) gives a plane through it and its ligands, centred on it:
      - ligands are the atoms within `cutoff` Å: default 2.6, or 2.8 for Na, K, Ca, Sr and Ba;
      - with five or more ligands (octahedral and up), the plane goes through the four that lie flattest with the
        metal, its equatorial plane, and `--check` names the ones out of the plane;
      - `"equatorial": false` uses all ligands;
    - `center` + `normal` gives a plane through that point; `normal` is a vector or two atom specs (from the first to
      the second).

    `center` may also be given with `atoms` to move the plane's centre.
  - `axis`: two atom specs that set the plane's horizontal direction.
  - `toward`: an atom spec (or point) on the side the normal should point to. A best plane's normal otherwise
    points to an arbitrary side; `offset` and `copies` count along the normal, and `--check` prints it.
  - `offset`: Å along the normal.
  - `copies`: extra parallel planes at these offsets from it (Å, e.g. `[-1, 1]`), named "<name> +1 Å".
  - `thickness`: Å. A thick slice shows the strongest value within it, with its sign kept.
  - `radius`: half the square's side (default `contour_radius`).
  - `map`: the map the plane was made for (default the first visible). The page shows the map in use instead
    when the viewer switches maps.
  - `level`: start level for this plane on its `map`. Choosing the plane in the list sets the map to it; any
    slider or scroll change on that map then moves all its planes.
  - `color`: outline colour, CSS hex. `visible`: default true.

**Atom specs**:
- `A/85 CA`, `A/85 LYS CA`, `2#A/85 SG` (model 2), `A/85` (residue: its heavy atoms), `A/52A` (insertion code).
- **Alternate conformers:** a spec without an altloc takes the first conformer (A); add `:B` for another
  (`A/85 CA:B`). `--at` lists every conformer.
- **Neighbour copies:** add the copy's tag in brackets, as `--check` prints it.
- Chain names are as in the file (`-` is a valid chain).
- An unknown atom gives an error that lists the atoms of that residue; an unknown residue lists the residues.

## The page

- **Maps:** visibility and a log-scale level slider for each map. Under each map: its rms and the level on the
  map's own scale (σ maps: in map units; absolute maps: as a multiple of the rms). `--check` prints the same.
- **Maps drawn as:** Mesh or Slice. Switching or changing maps never changes this mode. Slice mode shows slices
  only (no contours), of the map in use (the first visible); with every map off the planes are shown empty. In
  slice mode:
  - the live plane (yellow outline) faces the viewer and turns with the molecule until it is fixed;
  - sliders for depth, thickness and two tilts;
  - buttons: Fix plane, Through labelled atoms (best plane through three or more tapped atoms), Face me, and
    Save plane (S).
- **Planes:** the preset planes and the viewer's own saved ones, each with its outline colour. They are drawn only
  in slice mode, with the map in use at its level.
  - Tap a name to look at that plane, centred and face-on.
  - The checkbox hides it.
  - Copy adds a parallel plane at the copy offset (slider, −5 to +5 Å).
  - Edit makes it the live plane; × removes it.
  - The viewer's own planes and copies are kept in their browser. Presets come with the page.
- **Structure:** a checkbox for each model and for the neighbouring cells, plus hydrogens, Hide in front (residues
  in front of the centre residue that cover it are hidden with their density; on by default, key C), clear labels
  and Reset view.
- **Depth bar** (right edge of the view): the structure from front (top) to back.
  - Its window is the layer you see. Drag it to look deeper; drag its ends to make the layer thinner or thicker;
    tap the bar to jump.
  - The rotation centre (white dot) never moves with it. Double-tap an atom to make that atom the centre.
  - The live slice rides in the middle of the window (yellow tick). Keys: `,` `.` move, `[` `]` change the
    thickness.
- **Zoom** only magnifies: the camera is orthographic and never enters the structure.
- **Full screen** (button or F): level and depth sliders, Fix plane, Save, and an Atoms toggle (A).
- **Theme:** one deliberate dark look (Coot black).
- **Layout:** on an upright phone the panel sits under the view, otherwise beside it.
