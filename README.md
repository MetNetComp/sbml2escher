# sbml2escher

`sbml2escher.py` draws metabolic networks with nodes on a grid.
Cofactors may also be placed off the grid. Inorganic species are ignored.
Exchange reactions are placed on the outer boundary. Key pathways, such as
glucose to pyruvate, are drawn as straight lines and highlighted in red.

When flux sampling is available, reactions with zero median flux are not drawn.
The layout aims to minimize the total length of paths connecting the key
pathways to all exchange reactions. These connecting paths are shown in orange.

The criteria for selecting key pathways are provisional. The drawing criteria
are also provisional, and we welcome all ideas for improvement.

## Installation

Run these commands from the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-layout.txt -r requirements-escher.txt
```

The required Python packages are NumPy, SciPy, Matplotlib, and jsonschema.
Optimization uses HiGHS through SciPy; no separate solver executable or MATLAB
installation is required. A browser is needed to view the interactive output.
The Escher JavaScript viewer is bundled under `vendor/escher/`.

Run all examples from the project root with the virtual environment active.
Always specify `--model` to select an input file.

## Directory structure

```text
./
  README.md
  sbml2escher.py                  Main command-line entry point
  *.py                           Supporting Python modules
  requirements-layout.txt        Numerical and plotting dependencies
  requirements-escher.txt        Escher schema-validation dependency
  requirements-layout-lock.txt   Recorded layout dependency versions
  data/bacteria/                 Chemical-role reference catalog
  major_pathway.json             Canonical backbone reference
  pathway_output/                Additional canonical backbone references
  vendor/escher/                 Bundled viewer, schema, license, and provenance
  models/                        Input SBML models
  results/                       Saved drawing results
    manifest.json                Model-to-result mapping and verification record
```

## Input models: `models/`

This folder contains seven SBML models. The program also accepts BiGG/COBRA
model JSON through `--model`. Escher map JSON is an output format and cannot
be used as an input model.

| File | Description |
| --- | --- |
| [e_coli_core.xml](models/e_coli_core.xml) | E. coli core metabolic model |
| [iML1515.xml](models/iML1515.xml) | E. coli K-12 MG1655 model from BiGG |
| [iMM904.xml](models/iMM904.xml) | Saccharomyces cerevisiae model from BiGG |
| [Human-GEM.xml](models/Human-GEM.xml) | Human1, Human-GEM v1.19.0 |
| [iNESC2DN_CMap_Annotated.xml](models/iNESC2DN_CMap_Annotated.xml) | NeuroMap CMap with original CellDesigner annotations |
| [iNESC2DN_CTMap_Annotated.xml](models/iNESC2DN_CTMap_Annotated.xml) | NeuroMap CTMap with original CellDesigner annotations |
| [iNESC2DN_ComboMap_Annotated.xml](models/iNESC2DN_ComboMap_Annotated.xml) | NeuroMap ComboMap with original CellDesigner annotations |

The SBML reader targets metabolic models with fixed stoichiometry. Dynamic
rules, events, initial assignments, and variable stoichiometry are unsupported.
Models keep their original attribution and usage conditions.

When standard chemical identifiers are missing, the reader checks structured
SBML notes (names, abbreviations, formulas, and VMH links) and CellDesigner
species names. Matches use exact reference identifiers or names, including
explicit legacy D/L aliases. Ambiguous matches stay unresolved and can be set
with `--identity-map`. The evidence is saved in `id_mapping.json`.

Explicitly named exchange reactions can include a CellDesigner source/sink
symbol. Verified symbols are hidden only in the drawing, while the source model
keeps every participant and bound. These symbols trigger structural mode when
sampling is requested. CellDesigner coordinates and colors are not used to
select pathways or position nodes.

## Example commands

### Draw a structural map without flux sampling

```bash
python sbml2escher.py \
  --model models/e_coli_core.xml \
  --no-sampling \
  --coverage 100 \
  --output-dir results/runs/e_coli_core_structural
```

This displays all eligible reactions. Reactions with no visible participants
after display filtering are excluded. Arrows indicate directions permitted by
the source bounds: reversible reactions have arrows in both directions and
zero-bound reactions have no arrows. No flux or activity ranking is calculated.

A larger-model example:

```bash
python sbml2escher.py \
  --model models/iML1515.xml \
  --no-sampling \
  --output-dir results/runs/iML1515_structural
```

To calculate structural maps for every included model in Bash:

```bash
for model in models/*.xml models/*.sbml; do
  [ -f "$model" ] || continue
  name=$(basename "$model")
  python sbml2escher.py \
    --model "$model" \
    --no-sampling \
    --coverage 100 \
    --output-dir "results/runs/${name%.*}_structural"
done
```

These examples use `results/runs/` so the supplied result folders remain
available for comparison.

### View a drawing as an image

Open a PNG file to view the complete drawing.
Use the SVG version to zoom in without losing image quality.

| Model | PNG | SVG |
| --- | --- | --- |
| E. coli core | [PNG](results/e_coli_core/layout.png) | [SVG](results/e_coli_core/layout.svg) |
| iML1515 | [PNG](results/iML1515/layout.png) | [SVG](results/iML1515/layout.svg) |
| iMM904 | [PNG](results/iMM904/layout.png) | [SVG](results/iMM904/layout.svg) |
| Human1 | [PNG](results/Human-GEM/layout.png) | [SVG](results/Human-GEM/layout.svg) |
| iNESC2DN CMap | [PNG](results/iNESC2DN_CMap/layout.png) | [SVG](results/iNESC2DN_CMap/layout.svg) |
| iNESC2DN CTMap | [PNG](results/iNESC2DN_CTMap/layout.png) | [SVG](results/iNESC2DN_CTMap/layout.svg) |
| iNESC2DN ComboMap | [PNG](results/iNESC2DN_ComboMap/layout.png) | [SVG](results/iNESC2DN_ComboMap/layout.svg) |

The iNESC2DN drawings use SBML notes and CellDesigner identity information to
identify key pathways and cofactors. Red shows the matched pathways; orange
shows verified connections to exchange reactions. Unmatched compounds and
unreachable exchanges are listed in the result reports.

### Use a drawing in Escher

Open the [E. coli core Escher preview](results/e_coli_core/escher_preview.html)
in a browser to explore the map.

To load this map into another Escher session, download
[escher_map.json](results/e_coli_core/escher_map.json)
and select it as the map file.

To display the sampled reaction values, also download
[escher_flux_medians.json](results/e_coli_core/escher_flux_medians.json)
and load it as reaction data.

These values are sampling estimates, not measured fluxes.
See the [report](results/e_coli_core/report.json) for diagnostics.

GitHub file pages do not run the HTML previews. Download the repository,
then open `results/e_coli_core/escher_preview.html` in a browser.
To use the interactive layout instead, open `results/e_coli_core/index.html`.

### Discover pathways without drawing

```bash
python sbml2escher.py \
  --model models/iMM904.xml \
  --discover-only \
  --output-dir results/runs/iMM904_discovery
```

This writes `pathway_extraction.json` and `id_mapping.json` without sampling or
layout calculation.

### Draw a map using flux sampling

Sampling is enabled when `--no-sampling` is omitted:

```bash
python sbml2escher.py \
  --model models/e_coli_core.xml \
  --samples 2000 \
  --chains 4 \
  --thinning 100 \
  --burn-in 10000 \
  --seed 0 \
  --output-dir results/runs/e_coli_core_sampled
```

The model's original bounds are used without growth-floor or objective
conditioning. Reaction selection prioritizes identified backbones and then
absolute median flux. Zero-median reactions are excluded. Sampling diagnostics
can fail to meet convergence thresholds; such output is explicitly provisional.
If sampling is unavailable or fails, the program automatically falls back to
structural mode and records the reason in `report.json`.

To reuse samples generated by the preceding command:

```bash
python sbml2escher.py \
  --model models/e_coli_core.xml \
  --reuse-samples results/runs/e_coli_core_sampled \
  --coverage 50 \
  --output-dir results/runs/e_coli_core_sampled_50
```

Reuse requires the saved raw samples and an exactly matching parsed-model hash.
It cannot be combined with `--no-sampling`. Use a separate output directory for
structural runs rather than one containing sampling files.

### Inspect available options

```bash
python sbml2escher.py --help
```

Useful options include `--coverage` (0–100), `--time-limit` (the grid MILP time
limit, not a total runtime limit), `--max-milp-nodes`, `--grid-cols`, `--grid-rows`,
`--glucose-id`, `--identity-map`, and `--classification`. Identity-map files map
source metabolite IDs to reference roles; classification files map IDs to
`inorganic`, `coenzyme`, or `organic`.

## Saved drawings: `results/`

This folder contains drawings for the seven included models. Some results use
flux sampling; others show structure only. The table lists the mode and the
number of displayed reactions for each model.

| Result folder | Mode | Displayed / source reactions |
| --- | --- | ---: |
| [e_coli_core](results/e_coli_core/index.html) | Sampled; mixing checks passed | 76 / 95 |
| [iML1515](results/iML1515/index.html) | Sampled; convergence checks not met | 1,635 / 2,712 |
| [iMM904](results/iMM904/index.html) | Sampled; convergence checks not met | 848 / 1,577 |
| [Human-GEM](results/Human-GEM/index.html) | Structural; no sampling | 12,737 / 12,971 |
| [iNESC2DN_CMap](results/iNESC2DN_CMap/index.html) | Structural; no sampling | 1,834 / 1,915 |
| [iNESC2DN_CTMap](results/iNESC2DN_CTMap/index.html) | Structural; no sampling | 1,976 / 2,056 |
| [iNESC2DN_ComboMap](results/iNESC2DN_ComboMap/index.html) | Structural; no sampling | 2,222 / 2,305 |

[results/manifest.json](results/manifest.json) lists the model files, calculation
settings, validation results, and file hashes. Raw flux samples for E. coli core
are available in [flux_samples.npz](results/e_coli_core/flux_samples.npz).
Other supplied sampled results do not include raw samples. To use
`--reuse-samples`, use the E. coli core samples or generate samples with the
sampling command above.

Typical files in each result folder are:

| File | Purpose |
| --- | --- |
| `index.html` | Interactive layout; open in a browser |
| `layout.png`, `layout.svg` | Complete static drawing |
| `orange_routes.png`, `orange_routes.svg` | Backbone and EX-connector drawing, when provided by that result version |
| `escher_preview.html` | Browser-based Escher viewer |
| `escher_map.json` | Escher-format map |
| `source_model.json` | Parsed source model, including original reaction participants and bounds |
| `display_model.json` | Display projection after filtering |
| `coordinates.csv`, `edges.csv` | Saved node coordinates and reaction incidences |
| `report.json` | Layout, sampling, validation, and limitation report |
| `pathway_extraction.json`, `id_mapping.json` | Pathway matches and supporting chemical-role mappings |
| `reaction_ranking.csv` | Reaction selection and ordering |
| `exchange_connectors.json` | EX connections and unreachable exchanges, when available |

The display projection omits classified inorganic participants and verified
source/sink symbols. It must not be used as an FBA model. The original parsed model is retained in
`source_model.json`.

## Display conventions and limitations

- Red identifies recognized primary pathways. In results with connector
  highlighting, orange identifies shared undirected connections to exchange
  reactions; other branches are gray. Sampled results may also use flux colors.
- Ordinary metabolite and reaction nodes occupy distinct integer-grid points.
  Local cofactor copies normally sit off-grid near their reactions. Cofactors
  required to connect pathways to exchanges can instead become shared grid nodes.
- Exchange reactions are placed on the outer boundary. Layout constraints can
  differ between saved results; see each `report.json`.
- Large layouts use heuristic placement. Global optimality is not assumed, and
  edge crossings and label overlap are not minimized. Dense maps are best
  inspected through their interactive or SVG views.
- Connectivity and permitted directions do not establish steady-state
  feasibility or physiological pathway activity.
- The iNESC2DN files retain their CellDesigner annotations. Chemical-role matching
  uses a limited reference catalog, so some compounds remain unresolved. The
  reports list matching evidence and exchange reactions with no verified path.
  Restoring annotations does not prove physiological activity.
- Human1 and iMM904 also have partial pathway recognition. Consult the individual
  extraction reports rather than interpreting unmatched pathways as absent.
