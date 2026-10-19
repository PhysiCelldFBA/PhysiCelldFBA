# cancer_core_metabolism

A 3-D breast-cancer spheroid model in which each cell runs dFBA against a trimmed core model of
MCF7 metabolism (Zielinski et al., 2017), coupled to oxygen and nutrient gradients in the
microenvironment. Used to study how spatial oxygen/nutrient availability shapes flux distributions
across the spheroid (see `scripts/oxygen_dependency_experiment.ipynb`,
`scripts/reduced_cost_analysis.ipynb`).

## Contents

- `config/MCF7_Zielinski_2017_core_trimmed.sbml` — trimmed MCF7 core metabolic model (SBML), loaded
  via `sbml_filename` in both `PhysiCell_settings.xml` and `PhysiCell_settings_pfba.xml`.
- `config/PhysiCell_settings.xml` — default configuration (FBA).
- `config/PhysiCell_settings_pfba.xml` — alternate configuration running parsimonious FBA (pFBA)
  instead of plain FBA; pass it explicitly as `argv[1]` to use it.
- `config/cells_column.csv` — initial cell positions loaded by both settings files. The
  `scripts/cells_column.csv` copy is a stale output of `scripts/generate_cell_column.py`; only the
  `config/` copy is read by the simulation.
- `config/dfba_model.yaml` — standalone dFBA solver configuration kept for reference, not currently
  loaded by the XML pipeline.
- `MCF7_Zielinski_2017_exchanges.csv` / `.ods` — exchange-reaction bounds for the MCF7 core model,
  referenced when preparing/trimming the SBML model (see `scripts/preprocess_dfba_model.ipynb`).
- `custom_modules/custom.cpp`, `custom.h` — cell type setup, spheroid layout
  (`create_cell_sphere_positions`), and the dFBA coupling step.
- `scripts/` — preprocessing (`preprocess_dfba_model.ipynb`), spheroid cell-position generation
  (`create_cell_positions_spheroid.ipynb`, `generate_cell_column.py`), and figure generation
  (`cancer_tissue_figures.ipynb`, `figure_generation.py`, `plot_cancer_tissue_spatial*.py`,
  `plot_cells_spatial.py`, `plot_exchange_profiles.py`, `reduced_cost_analysis.ipynb`).

## Build & run

From this directory (standalone):
```bash
make
./cancer-core-metabolism                      # default FBA config
./cancer-core-metabolism config/PhysiCell_settings_pfba.xml   # pFBA variant
```

Or from the repository root:
```bash
make cancer-core-metabolism
make
```

## Citation

See `ALL_CITATIONS.txt`. The metabolic model is derived from Zielinski DC et al.,
"Systems biology analysis of drivers underlying hallmarks of cancer cell metabolism", Sci. Rep. 7,
41241 (2017).
