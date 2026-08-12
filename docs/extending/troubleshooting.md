# Troubleshooting

## Dependency download does not run

The Makefile checks for:

```text
addons/dFBA/ext/coin-or/include/coin/CoinPackedMatrix.hpp
```

If the header exists, automatic setup is skipped. If a package directory is non-empty but incomplete, `setup_fba.py` can also skip it. Inspect the Coin-OR and libSBML sentinel headers described in [Installation](../getting-started/installation.md).

## Download fails

Check network access, write permission for `addons/dFBA/ext/`, and the platform key in `beta/fba_packages.json`. The automatic manifest currently covers Linux x86-64, Windows 64-bit, macOS Intel, and macOS Apple Silicon.

## Compiler cannot find SBML or Coin headers

Confirm the dependencies are unpacked beneath:

```text
addons/dFBA/ext/libsbml/
addons/dFBA/ext/coin-or/
```

Then check that the selected sample Makefile uses `$(EXT_DIR)/libsbml/include` and `$(EXT_DIR)/coin-or/include`. A Makefile from an older non-dFBA project will not contain those paths.

## Linker reports missing CLP, CoinUtils, or libSBML symbols

The current runtime links the Coin-OR static archives and libSBML. Verify the corresponding files exist under the dependency `lib/` directories and that your compiler/architecture matches the downloaded packages.

## SBML model does not load

Common causes:

- invalid XML or fatal libSBML diagnostics;
- unsupported/missing FBC bounds;
- a machine-specific or incorrect `sbml_filename`;
- a file copied to a different project location;
- incompatible reaction IDs after SBML export.

Validate the file independently, then run from the repository root so relative `./config/...` paths resolve consistently.

## “substrate ... was not found in the microenvironment”

The `exchange substrate="..."` name does not exactly match a `<microenvironment_setup><variable name="...">` entry. Names are case-sensitive. Create the BioFVM field or correct the mapping.

## “exchange reaction not found in model”

The `fba_flux` ID does not exist in the loaded SBML model. Do not use the reaction's display name. Inspect the exact SBML `reaction id` after any COBRApy/import-export transformations.

## Objective reaction assertion or initialization failure

Confirm `objective_reaction` exists and the model contains a valid objective. Also check `max_growth_rate` is non-negative and expressed as h⁻¹ for the current implementation.

## Infeasible solutions

Test the same model and bounds outside PhysiCell. Frequent causes include:

- maintenance or objective demands that cannot be met;
- incorrect exchange sign conventions;
- a required nutrient omitted from the microenvironment;
- `Vmax = 0` on a required uptake;
- low local concentration relative to `Km`;
- a death-trigger lower bound that makes the network infeasible;
- over-aggressive model trimming.

An infeasible solve sets metabolic growth and reaction fluxes to zero and can activate metabolic death if enabled.

## Negative concentrations or mass mismatch

Check in this order:

1. `intracellular_dt = dt_diffusion`.
2. Concentrations are in mM.
3. Fluxes are in mmol/gDW/h.
4. `cell_density` and `reference_volume` have the intended physical meaning.
5. Uptake reactions use negative flux.
6. Voxel size and initial substrate mass match the analytical calculation.
7. No custom module overwrites net export rates after the dFBA update.

Use the mass-conservation unit test with the new model before debugging a large population.

## Growth is about 60-fold too fast or slow

This usually indicates confusion between h⁻¹ and min⁻¹, or between fmol/pg/min and mmol/gDW/h. XML unit labels do not convert values. The current growth update treats the objective as h⁻¹ and explicitly divides by 60.

## Simulation is unexpectedly slow

The dominant work can be one LP solve per cell per diffusion timestep. Reduce the problem systematically:

- use fewer cells and a shorter end time for debugging;
- validate a reduced metabolic reconstruction;
- remove blocked reactions;
- use a coarser spatial model only when scientifically acceptable;
- profile before changing the synchronization rule.

Do not increase `intracellular_dt` independently to gain speed; that mode is not supported because it introduces stale uptake bounds.

## Studio writes an unexpected dFBA block

Check that the Studio version supports the same PhysiCelldFBA XML structure. Save to a new file, compare the generated `intracellular` block with [the documented structure](../configuration/physicell-xml.md), and verify the timestep and reaction IDs manually.
