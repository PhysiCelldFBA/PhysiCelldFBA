# crossfeeding

A two-species syntrophic crossfeeding model: *Clostridium beijerinckii*
(`config/sbml_models/C_beijerinckii_model_iCB925_trimmed.xml`) and *Methanosarcina barkeri*
(`config/sbml_models/M_barkeri_model_iMG746_trimmed.xml`) each run dFBA with species-specific
post-update coupling (`post_update_intracellular_c_beijerinckii`,
`post_update_intracellular_m_barkeri`), exchanging metabolites (e.g. H2, acetate) through the
shared microenvironment.

This project is actively developed on the `crossfeeding_dev` branch; this directory reflects its
state on `development` as of this branch's sync.

## Build & run

From the repository root:
```bash
make crossfeeding
make
./crossfeeding
```

Figures are regenerated from a simulation output folder via
`scripts/paper_figures/run_paper_figures.py` — see `scripts/paper_figures/README.md`.

## Contents

- `config/sbml_models/` — the two trimmed genome-scale models, loaded via `sbml_filename`
  in `PhysiCell_settings.xml`.
- `config/PhysiCell_settings.xml` — simulation configuration (the only config; a byte-identical
  `PhysiCell_settings_crossfeeding.xml` duplicate was removed).
- `config/200cells.csv` — initial cell positions for both species.
- `config/rules.csv` — PhysiCell cell-behavior rules.
- `custom_modules/custom.cpp`, `custom.h` — species-specific dFBA coupling, Dirichlet-node
  substrate control, and the cross-feeding exchange logic.
- `scripts/paper_figures/` — self-contained figure-reproduction pipeline (own `README.md`,
  `requirements.txt`, and `run_paper_figures.py`) for the manuscript's crossfeeding panels. See
  that README for prerequisites (including the `pctk` MultiCellDS reader) and usage.

## Citation

See `ALL_CITATIONS.txt`.
