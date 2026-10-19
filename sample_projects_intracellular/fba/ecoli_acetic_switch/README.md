# ecoli_acetic_switch

A dFBA model of the *E. coli* overflow/acetate-switch behavior: agents run flux balance analysis
against the core *E. coli* model (`config/Ecoli_core.xml`) and switch between glucose and acetate
uptake as glucose is depleted (`reintroduce_nutrients_function`, `inject_density` handle periodic
nutrient refeeding for longer runs). This is the paper's reference scenario — the genome-scale
(`iML1515`) and literature-parameter variants were trimmed in favor of the single scenario the
manuscript figures are generated from.

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

## Build & run

From this directory:
```bash
make
./ecoli-dfba
python scripts/plot_ecoli_figure_a.py --output-folder output
python scripts/plot_ecoli_figure_b.py --output-folder output
```

Or from the repository root: `make ecoli-acetic-switch-sample && make`.

## Citation

See `ALL_CITATIONS.txt`.
