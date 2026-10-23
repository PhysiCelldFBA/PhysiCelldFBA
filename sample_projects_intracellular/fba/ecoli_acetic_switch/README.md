# ecoli_acetic_switch

A dFBA model of the *E. coli* overflow/acetate-switch behavior: agents run flux balance analysis
against the core *E. coli* model (`config/Ecoli_core.xml`) and switch between glucose and acetate
uptake as glucose is depleted (`reintroduce_nutrients_function`, `inject_density` handle periodic
nutrient refeeding for longer runs). This is the paper's reference scenario — the genome-scale
(`iML1515`) and literature-parameter variants were trimmed in favor of the single scenario the
manuscript figures are generated from.

## Build & run

From the repository root:
```bash
make ecoli-acetic-switch-sample
make
./ecoli-dfba
python scripts/plot_ecoli_figure_a.py --output-folder output
python scripts/plot_ecoli_figure_b.py --output-folder output
```

## Contents

- `config/Ecoli_core.xml` — the *E. coli* core metabolic model, loaded via `sbml_filename` in
  `PhysiCell_settings.xml`.
- `config/PhysiCell_settings.xml` — the only configuration; no alternate scenarios.
- `config/cells_reduced.csv` — initial cell positions.
- `custom_modules/custom.cpp`, `custom.h` — cell setup and the acetate-switch coupling logic.
- `scripts/plot_ecoli_figure_a.py` — reproduces the paper's Figure A (acetate-switch time series:
  glucose/acetate/O₂/CO₂ on the left axis, total biomass on the right, 0–10 h).
- `scripts/plot_ecoli_figure_b.py` — reproduces the paper's Figure B (spatial glucose/acetate
  panels at 0.5 h / 2 h / 4.5 h). Both are self-contained CLI scripts (stdlib +
  numpy/matplotlib/scipy only); `--output-folder` is required (no implicit default).

## Figure reproducibility

| Manuscript figure | Panels | Script | XML input | Expected output |
|---|---|---|---|---|
| Main Figure 3 (acetate-switch dynamics) | Single time course, 0–6 h: glucose, oxygen, acetate, CO₂, total biomass | `python scripts/plot_ecoli_figure_a.py --output-folder <sim_output_dir>` | `config/PhysiCell_settings.xml` | `figure_A.png`, `figure_A.svg` |
| Supplementary Figure S5 (glucose and acetate fields) | Snapshots at 0.5 h, 2 h and 4.5 h; glucose and acetate rows | `python scripts/plot_ecoli_figure_b.py --output-folder <sim_output_dir>` | `config/PhysiCell_settings.xml` | `figure_B.png`, `figure_B.svg` |

Both scripts require `--output-folder` (no default). Figure A is drawn on a fixed 0–6 h axis, so the
simulation must cover at least 6 h of simulated time; Figure B needs output up to 4.5 h. Outputs
default to `analysis/results/ecoli/<tag>/` and can be redirected with `--results-dir`.

## Citation

See `ALL_CITATIONS.txt`.
