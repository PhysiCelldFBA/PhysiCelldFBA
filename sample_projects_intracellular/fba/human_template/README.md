# cancer_metabolism

A cancer-cell dFBA model built on a trimmed, cell-type-specific reduction of the genome-scale
human metabolic reconstruction Recon2.2 (`config/Recon2.2_FBC_trimmed_final.xml`). Each agent
solves this model against the local microenvironment every intracellular update, coupling cell
metabolism to oxygen/nutrient gradients.

**Note:** this project is the base for `human_template` (see
`../human_template/README.md`), the template sample project for building new
human/mammalian-metabolism PhysiCelldFBA cases.

## Contents

- `config/Recon2.2_FBC_trimmed_final.xml` — trimmed Recon2.2 SBML model, loaded via
  `sbml_filename` in `PhysiCell_settings.xml`.
- `config/iSIM.xml` — an earlier/alternate SBML model kept in the repository for reference; it is
  not loaded by the current `PhysiCell_settings.xml` (which points at `Recon2.2_FBC_trimmed_final.xml`).
- `config/PhysiCell_settings.xml` — simulation configuration (microenvironment, cell definitions,
  dFBA model path).
- `config/cells.csv` — initial cell positions.
- `custom_modules/custom.cpp`, `custom.h` — cell type setup, tissue layout, and the
  `post_update_intracellular` dFBA coupling step, including periodic nutrient reintroduction
  (`reintroduce_nutrients_function`).
- `CMakeLists.txt` — CMake build definition (alternative to the Makefile below).

## Build & run

From this directory (standalone):
```bash
make
./cancer-metabolism
```

Or from the repository root:
```bash
make cancer-metabolism-sample
make
```

## Citation

See `ALL_CITATIONS.txt`. The metabolic model is derived from Recon 2.2 (Swainston N et al.,
"Recon 2.2: from reconstruction to model of human metabolism", Metabolomics 12, 109, 2016).
