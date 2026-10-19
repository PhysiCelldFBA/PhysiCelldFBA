# crossfeeding

A two-species syntrophic crossfeeding model: *Clostridium beijerinckii*
(`config/sbml_models/C_beijerinckii_model_iCB925_trimmed.xml`) and *Methanosarcina barkeri*
(`config/sbml_models/M_barkeri_model_iMG746_trimmed.xml`) each run dFBA with species-specific
post-update coupling (`post_update_intracellular_c_beijerinckii`,
`post_update_intracellular_m_barkeri`), exchanging metabolites (e.g. H2, acetate) through the
shared microenvironment. Includes scheduled glucose refeeding
(`apply_glucose_dirichlet_refeed_schedule`, `apply_initial_glucose_refeed_setup`) to sustain
longer runs.

This project is actively developed on the `crossfeeding_dev` branch; this directory reflects its
state on `development` as of this branch's sync.

## Contents

- `config/sbml_models/` — the two trimmed genome-scale models actually loaded via `sbml_filename`
  in `PhysiCell_settings*.xml`. The `yaml/` subfolder (`c_beijerinckii.yaml`, `m_barkeri.yaml`,
  `crossfeeding.yaml`, `models.yaml`) is not referenced by any config/script in this project as of
  this writing — kept for reference pending confirmation of whether it's still needed.
- `config/PhysiCell_settings.xml` / `PhysiCell_settings_crossfeeding.xml` — simulation
  configuration; currently identical copies of the same crossfeeding scenario.
- `config/200cells.csv` — initial cell positions for both species.
- `config/rules.csv` — PhysiCell cell-behavior rules.
- `custom_modules/custom.cpp`, `custom.h` — species-specific dFBA coupling, Dirichlet-node
  substrate control, and the cross-feeding exchange logic.
- `scripts/paper_figures/` — self-contained figure-reproduction pipeline (own `README.md`,
  `requirements.txt`, and `run_paper_figures.py`) for the manuscript's crossfeeding panels. See
  that README for prerequisites (including the `pctk` MultiCellDS reader) and usage.

## Build & run

From this directory:
```bash
make
./crossfeeding config/PhysiCell_settings_crossfeeding.xml
```

Or from the repository root: `make crossfeeding && make`.

Figures are regenerated from a simulation output folder via
`scripts/paper_figures/run_paper_figures.py` — see `scripts/paper_figures/README.md`.

## Citation

See `ALL_CITATIONS.txt`.
