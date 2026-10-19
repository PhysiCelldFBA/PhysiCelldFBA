# metabolic_driven_motility

Demonstrates dFBA-driven chemotaxis-like behavior: each *E. coli* agent (core model,
`config/Ecoli_core.xml`) can switch its FBA objective from biomass production to ATP production
(`custom_atp_optimization`), letting a cell trade growth for directed migration toward a
glucose-richer region of the microenvironment (`metabolic_bound_migration_rule`).

## Contents

- `config/Ecoli_core.xml` — core *E. coli* model used by all scenario configs.
- `config/PhysiCell_settings.xml` — default configuration.
- `config/PhysiCell_settings_motility.xml`, `PhysiCell_settings_motility_proliferation.xml`,
  `PhysiCell_settings_no_motility.xml`, `PhysiCell_settings_multiple_sources.xml` — alternate
  scenarios (motility only, motility + proliferation, proliferation without motility, multiple
  glucose sources); select by passing the file as `argv[1]`.
- `config/cells.csv` — initial cell positions, shared by all scenarios.
- `custom_modules/custom.cpp`, `custom.h` — the objective-switching logic
  (`custom_atp_optimization`), secretion-node setup (`setup_secretion_nodes`), and the migration
  rule.
- `scripts/composite_figure.py`, `figure_generation.py`, `dynamic_analysis.py`,
  `fluxes_analysis.py`, `dynamics_analysis.ipynb` — figure and flux-dynamics analysis.
- `scripts/plot_df_journal_style_ext.py` — shared plotting style helper (also used by
  `ecoli_acetic_switch`/`dfba_unit_test`).

## Known issue

Every config variant in this project (including the default `PhysiCell_settings.xml`) sets a
substrate initial-condition `<filename>` of `./config/glucose.csv` — in the default file this is
hardcoded to a developer workstation path
(`/home/mruscone/Desktop/github/PhysiCelldFBA/config/glucose.csv`) — but no `glucose.csv` is
present anywhere in this repository. As committed, none of the five configs here will load
correctly; a `glucose.csv` initial-condition file needs to be supplied (and the hardcoded path
corrected) before running. Flagged here for visibility; not fixed as part of this cleanup pass.

## Build & run

From this directory (once `config/glucose.csv` is supplied, see above):
```bash
make
./ecoli-dfba-motility config/PhysiCell_settings_motility.xml
```

Or from the repository root: `make metabolic_driven_motility && make`.

## Citation

See `ALL_CITATIONS.txt`.
