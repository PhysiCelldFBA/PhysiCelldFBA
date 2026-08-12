# Framework overview

PhysiCelldFBA joins three modelling layers:

| Layer | Responsibility |
| --- | --- |
| PhysiCell | Off-lattice cells, mechanics, volume, division, death, and behaviors |
| BioFVM | Spatial diffusion, decay, boundary conditions, and extracellular substrate concentrations |
| PhysiCelldFBA | Per-cell SBML model, exchange constraints, FBA solution, growth, and metabolite feedback |

## One model instance per cell

A cell definition with `<intracellular type="dfba">` provides a template. During initialization and cell creation, PhysiCelldFBA loads the referenced SBML model and keeps an independent metabolic-model state for each agent. This is what allows genetically identical cells to obtain different flux distributions when they occupy different microenvironments.

Multiple cell definitions can point to:

- the same SBML model with different kinetic or phenotype parameters; or
- different SBML models, as in a multispecies community.

## Bidirectional coupling

The coupling is not a one-way lookup:

1. The cell samples local BioFVM concentrations.
2. Those concentrations constrain metabolic uptake.
3. FBA returns an objective value and exchange fluxes.
4. The objective changes cell growth or a custom behavior.
5. Exchange fluxes become BioFVM source/sink terms.
6. Diffusion and neighboring cells change the environment before the next solve.

Consequently, temporal adaptations and spatial metabolic niches can emerge without hard-coding a phenotype switch. The acetate-switch and colony examples use the same principle in temporal and spatial settings.

## Optimization problem

For a stoichiometric matrix `N`, reaction-flux vector `v`, objective coefficients `c`, and reaction bounds, each cell solves:

```text
maximize    Z = cᵀv
subject to  Nv = 0
            v_lower ≤ v ≤ v_upper
```

The objective usually represents biomass synthesis. PhysiCelldFBA also exposes a custom optimization hook: the metabolism-driven motility example temporarily optimizes ATP maintenance and converts surplus ATP flux into migration speed.

## What is spatially explicit

Cells move continuously and interact mechanically without being fixed to a lattice. BioFVM still resolves concentrations on diffusion voxels. At each update a cell reads the voxel containing its position, but its position and mechanics remain off-lattice. This division of responsibilities supports dense colonies and tissues while retaining efficient reaction-diffusion calculations.

## Where the implementation lives

The dFBA runtime is under `addons/dFBA/src/`. Sample projects combine it with their own `main.cpp`, `custom_modules/`, configuration, metabolic models, and analysis scripts under `sample_projects_intracellular/fba/`. See the [add-on code guide](../addon/index.md).
