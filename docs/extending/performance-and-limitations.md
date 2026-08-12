# Performance and limitations

## Computational cost

Each metabolic cell owns an independent constraint-based model and solves a linear program at every supported dFBA/diffusion update. Runtime therefore grows with:

- number of metabolic agents;
- number of reactions and metabolites;
- simulation duration divided by `dt_diffusion`;
- difficulty of the changing LPs;
- cost of BioFVM diffusion and cell mechanics.

Genome-scale models can also create substantial per-cell memory pressure because the model and solver state are copied.

## Practical scaling strategy

1. Validate one cell with the full model.
2. Remove blocked and irrelevant reactions while preserving reference fluxes.
3. Benchmark tens, hundreds, and then thousands of cells.
4. Profile LP solves, diffusion, mechanics, and output independently.
5. Reduce output frequency during scale tests.
6. Keep a small full-output validation run for scientific checks.

BioFVM and PhysiCell already use parallel implementations for major tasks, and cell-level optimizations are conceptually independent. Actual scaling depends on the solver calls, memory bandwidth, model size, and build/platform.

!!! warning

    Do not trade correctness for speed by setting a coarser `intracellular_dt` than `dt_diffusion`. The current supported coupling requires equality to keep exchange bounds synchronized with the available voxel mass.

## Biological limitations

### Objective assumption

FBA requires an objective. Biomass maximization is useful for many growth systems but does not represent every cellular strategy. Maintenance, stress, migration, persistence, and competing functions can require custom or multi-objective formulations.

### Alternative optima

Different flux distributions can produce the same optimal objective. If alternative solutions differ in exchange flux, they can change extracellular gradients and downstream population behavior even though growth is unchanged.

### No explicit regulation by default

Environmental constraints cause metabolic adaptation, but the generic dFBA layer does not model transcriptional, signalling, or enzyme-capacity regulation. Such mechanisms require additional constraints or coupling to other intracellular modules.

### Parameter uncertainty

Transport kinetics, density, maintenance, diffusion, and boundary concentrations are often incompletely measured. A spatially detailed output is not automatically quantitatively predictive. Calibrate against independent data and report sensitivity.

### Steady-state intracellular assumption

Each FBA solve assumes intracellular pseudo-steady state. Fast intracellular transients, metabolite pools, and enzyme dynamics are not represented unless the model is extended.

## Numerical limitations

- Local uptake uses the concentration and amount in the current BioFVM voxel.
- A smaller voxel changes the mass available to a cell per update.
- Growth is continuous but cell division is discrete.
- Solver tolerances and degeneracy can affect small fluxes.
- An infeasible model produces zero fluxes and can feed into death logic.
- Custom behaviors can break accounting if they change secretion or intracellular outputs inconsistently.

## Responsible interpretation

Report model provenance, constraints, units, timestep, voxel size, solver, random seed, and calibration data. Distinguish robust qualitative behavior from parameter-sensitive quantitative predictions, particularly for genome-scale tissue models.
