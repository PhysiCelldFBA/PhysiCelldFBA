# human_template

**Template sample project for PhysiCelldFBA.** Start here when building a new human/mammalian
dFBA case: copy this directory, swap in your own SBML model and cell layout, and adjust
`custom_modules/custom.cpp` as needed. Out of the box it runs a working example: a cancer-cell
dFBA model built on a trimmed reduction of the genome-scale human metabolic reconstruction Recon2.2
(`config/Recon2.2_FBC_trimmed_final.xml`). Each agent solves this model against the local
microenvironment every intracellular update, coupling cell metabolism to oxygen/nutrient
gradients.

## Using this as a template

1. Copy this directory to `sample_projects_intracellular/fba/<your_project>/`.
2. Replace `config/Recon2.2_FBC_trimmed_final.xml` with your SBML model, and update the
   `sbml_filename` entry in `config/PhysiCell_settings.xml` to point at it.
3. Update `config/cells.csv` (initial cell positions) and the microenvironment substrates in
   `PhysiCell_settings.xml` for your scenario.
4. Adjust `custom_modules/custom.cpp`/`custom.h` — cell type setup, tissue layout, and the
   `post_update_intracellular` dFBA coupling step — for your model's exchange reactions and
   behavior rules.
5. Set `PROGRAM_NAME` in `Makefile`, and the target name in `CMakeLists.txt`
   (`set(project ...)`), to your project's name.
6. Add a build rule for your project in `sample_projects/Makefile-default` (mirrored in the
   top-level `Makefile`) following the `human-template-sample` rule as a model, and add an
   `add_subdirectory(...)` line to the top-level `CMakeLists.txt` if you use CMake.

The settings file also contains disabled optional inputs: a cell-rules ruleset (`rules.csv`) and
MATLAB-format initial and Dirichlet conditions (`initial.mat`, `dirichlet.mat`). They are switched
off in this template and the files are not included. To use them, supply the files and set
`enabled="true"` on the corresponding block in `PhysiCell_settings.xml`.

## Build & run

From the repository root:
```bash
make human-template-sample
make
./human-template
```

## Contents

- `config/Recon2.2_FBC_trimmed_final.xml` — trimmed Recon2.2 SBML model, loaded via
  `sbml_filename` in `PhysiCell_settings.xml`.
- `config/PhysiCell_settings.xml` — simulation configuration (microenvironment, cell definitions,
  dFBA model path).
- `config/cells.csv` — initial cell positions.
- `custom_modules/custom.cpp`, `custom.h` — cell type setup, tissue layout, and the
  `post_update_intracellular` dFBA coupling step, including periodic nutrient reintroduction
  (`reintroduce_nutrients_function`).
- `Makefile` — build rules; `PROGRAM_NAME` sets the executable name.
- `main.cpp` — program entry point; copied to the repository root by the build rule.
- `CMakeLists.txt` — CMake build definition (alternative to the Makefile). It requires the dFBA
  addon's external dependencies (COIN-OR Clp and libSBML), which are downloaded at configure time
  by `addons/dFBA/CMakeLists.txt`.
- `README.md` — this file.

## Citation

See the PhysiCell core citation in the repository root `README.md`/`CITATION.txt`. The metabolic
model shipped as the default demo is derived from Recon 2.2 (Swainston N et al., "Recon 2.2: from
reconstruction to model of human metabolism", Metabolomics 12, 109, 2016).
