# Add a new metabolic model

Use a staged workflow. A genome-scale model that optimizes in COBRApy is not automatically ready for a spatial, mass-conserving agent simulation.

## 1. Define the scientific contract

Write down:

- the organism or cell type;
- the objective and its interpretation;
- extracellular metabolites that must diffuse;
- expected uptake/secretion directions;
- medium and boundary conditions;
- how metabolic failure affects the phenotype;
- measurements available for calibration and validation.

This prevents the mapping from being driven only by reaction names.

## 2. Validate the SBML model independently

Load the model in an established constraint-based tool and test each relevant medium:

- rich/replete conditions;
- depletion of each limiting substrate;
- secretion and re-uptake regimes;
- maintenance requirements;
- alternate objectives if used.

Record the objective, key exchange fluxes, and expected yield. These values become reference checks for the coupled model.

## 3. Create a project

Copy a current supported dFBA sample with the closest structure:

- the unit test for a single-species growth model;
- the colony for spatial microbial growth;
- cancer core metabolism for a large metabolite panel;
- cross-feeding for multiple SBML models;
- metabolic motility for a custom objective/behavior.

Rename the executable, output directory, cell definitions, and configuration files.

## 4. Add dependencies to the build

Retain the dFBA compile definition, include paths, objects, libraries, and dependency sentinel from the sample Makefile. The first build should be able to invoke `beta/setup_fba.py` automatically.

## 5. Define BioFVM substrates

Add only extracellular metabolites that participate in the spatial coupling. For each one, select:

- a stable name;
- diffusion coefficient and decay rate;
- initial condition;
- boundary condition;
- physical units.

The same name must appear in the cell's `exchange/@substrate` mapping.

## 6. Add the dFBA cell definition

Create the `settings`, `transport_model`, and `growth_model` blocks described in [PhysiCell XML](../configuration/physicell-xml.md). Keep:

```text
intracellular_dt = dt_diffusion
```

Use relative SBML paths and reaction IDs copied exactly from the model.

The [helper scripts](../configuration/helper-scripts.md) can produce a first YAML/XML draft for large exchange panels. Review every generated parameter.

## 7. Validate one cell

Run one cell in a small, closed domain with a known substrate mass:

1. compare its FBA objective with the independent solve;
2. compare exchange directions and magnitudes;
3. verify no concentration becomes negative;
4. compare final biomass with a yield calculation;
5. test nutrient depletion and infeasibility.

Do not scale to thousands of agents until this test passes.

## 8. Add space and mechanics

Increase complexity gradually:

1. introduce boundary sources;
2. create a small population;
3. inspect spatial concentrations and single-cell fluxes;
4. enable growth/division and death;
5. increase domain, cell count, and model size.

At every stage, preserve a small regression scenario that can be rerun after parameter or code changes.

## 9. Document reproducibility

For a publishable use case, record:

- SBML source, version, license, and preprocessing;
- all changed bounds and objective reactions;
- extracellular mappings and kinetic sources;
- timestep, voxel, domain, boundary, and cell settings;
- random seed and initial positions;
- exact run and analysis commands;
- expected qualitative and quantitative outputs.
