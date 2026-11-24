---
name: moorlight-share
description: Make a self-contained, phone-friendly Moorlight page (models + maps + a short tour) that needs no server, for sending as a file or publishing as a claude.ai artifact. Use when the user wants to look at a structure and its maps from a phone, a Remote Control session or a link, or to share a guided view of a result.
---

# Shareable Moorlight snapshots

`moorlight share` turns a small JSON spec into one HTML page. The page holds the models, the maps (resampled on a
box around the models) and a list of tour stops. The browser contours the maps itself, with chicken wire as in
Coot. The viewer can rotate, pinch-zoom, switch maps, change contour levels, tap atoms for labels and map values,
and step through the stops. Nothing runs on the cluster after the page is written.

## Workflow

1. **Write the spec** (`spec.json`, format below). Paths are relative to the spec file.
2. **Get the numbers**: `moorlight share spec.json --check [--at "<atom or residue>" ...]`.
   This prints:
   - for each map, its strongest peaks of each sign, each listed once at its image next to the model, with the
     nearest atom;
   - for each stop, the values at its labelled atoms and the peaks within its radius;
   - with `--at`, the values at the atoms you ask for.

   Every statement in a stop text should come from these numbers.
3. **Build**: `moorlight share spec.json -o page.html`. The command prints the size. An artifact page must stay at
   or below 16 MB. One map over a 30 Å box at 0.3 Å is about 3 MB.
4. **Hand it over**:
   - Publish it as an artifact (Artifact tool: `file_path` = the page, `icon` "molecule", a one-sentence
     `description`). The page is written in the artifact form, without `<!doctype>`/`<html>`/`<head>`/`<body>`.
   - Or build with `--standalone` for a file someone opens locally.

   The page cannot be rendered on a cluster node without a browser. Say so if you could not look at it.

## Spec

```json
{
  "title": "Cage 3.5 ps",
  "description": "Light minus dark at 3.5 ps on the dark structure; contours in noise SD.",
  "models": [{"path": "dark.pdb", "name": "dark", "hydrogens": false}],
  "maps": [{"path": "3p5ps_corrected.ccp4", "name": "3.5 ps light - dark", "kind": "diff",
            "units": "absolute", "unit_label": "SD", "level": 4.0, "visible": true}],
  "region": {"margin": 5},
  "contour_radius": 10,
  "stops": [
    {"name": "Azo bond, guest 102", "text": "Two to four sentences, with numbers.",
     "center": "-/102", "radius": 7, "orient": "edge", "frame": ["-/102 N9", "-/102 N10", "-/102 C8"],
     "maps": {"3.5 ps light - dark": 4.0}, "labels": ["-/102 N9", "-/102 N10"],
     "distances": [["-/102 N9", "-/102 N10"]]}
  ]
}
```

- **title**: the page's name, two to four words ("Cage 3.5 ps"). Never a name followed by a dash or colon and an
  explanation; the explanation goes in `description`.
- **models**: `path` (PDB/mmCIF), `name`, optional `color` (carbon colour, CSS), `hydrogens` (default false),
  `style` (`sticks` default, or `lines`), `visible`.
- **maps**: `path` (CCP4/MRC map, or MTZ with `columns: "F,PHI"`), `name` (unique; stops refer to it), `kind`
  (`diff`: green + / red − contours, or `2fofc`: one blue contour), `units` and `level`, `visible`, optional `step`
  (Å, default the map's grid spacing).
  - `units: "sigma"` (default): levels are multiples of the whole map's rms, as in Coot's default.
  - `units: "absolute"`: levels are in the map's own units. Use this for maps that are already scaled to noise SD.
    Their rms includes signal, so σ contours would mislead. `unit_label` is the text after the number ("SD").
- **region**: the map box, by default all models plus `margin` Å (default 5). For a smaller page use
  `{"center": "<atom spec or [x, y, z]>", "radius": 12}`.
- **contour_radius**: how far from the view centre maps are contoured, in Å (default 10, as Coot's map radius).
- **stops** (optional; without stops the page opens on the whole model):
  - `center`: atom spec, residue spec (centroid), a list of specs (centroid) or `[x, y, z]`.
  - `radius`: Å of context that fits the screen (default 8). `slab`: depth half-thickness (default `radius`).
  - `orient`: `face` (default; the frame atoms' plane faces the viewer, longest axis horizontal), `edge` (that
    plane seen edge-on, which shows out-of-plane motion), `keep` (the previous stop's rotation) or a quaternion
    `[x, y, z, w]`. `frame` lists the atoms that define the plane; the default is the centre's atoms. A single
    atom keeps the rotation.
  - `maps`: `{map name: level}` shows that map at that level, `false` hides it. Maps not named go back to their
    defaults from `maps`, so every stop looks the same however the viewer got there.
  - `labels`: atom specs to label. `distances`: pairs of atom specs.
  - `text`: plain text. A blank line starts a new paragraph.

**Atom specs**: `A/85 CA`, `A/85 LYS CA:B` (altloc B), `2#A/85 SG` (molecule 2 in spec order), `A/85` (residue:
its heavy atoms, or their centroid as a centre). Chain names are as in the file, for example `-` for a blank-ish
chain. An unknown atom gives an error that lists the atoms of that residue.

## Writing the tour

- Start with an overview stop (the whole molecule, `radius` large enough), then one stop per finding.
- One finding per stop, two to four sentences. Name the map and the level, and give the numbers from `--check`
  ("−13 SD on N9, +11.5 SD 0.9 Å from it"). The page is read on a phone; keep it short.
- Show each map only where it is needed and say which one is on.
- Report what the map shows. Interpretation that the user has not established does not belong in the texts.
- Check every stop in the `--check` output: are the centre and labels right, and is the strongest feature where
  the text says?
