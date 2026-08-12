# The dFBA update cycle

Each dFBA update has three conceptual stages: constrain, optimize, and feed back.

## 1. Sense the microenvironment and constrain exchanges

For every configured `<exchange>`, the cell:

1. finds the BioFVM substrate by the exact `substrate` name;
2. reads its concentration in the cell's current voxel;
3. evaluates the configured Michaelis-Menten uptake limit from `Km` and `Vmax`;
4. calculates the maximum mass that can physically be removed from the voxel during the update;
5. applies the more restrictive value as the lower bound of the corresponding SBML `fba_flux`.

This combines biological transport kinetics with a numerical mass-availability constraint.

## 2. Solve the metabolic model

Coin-OR CLP solves the cell's constrained FBA problem. An optimal solution provides:

- the objective flux, normally biomass production;
- a flux for each configured exchange reaction; and
- the complete intracellular flux distribution stored in the dFBA solution.

An infeasible solution produces zero metabolic growth and can activate the configured metabolic-death behavior. An unknown solver status is treated as an error.

## 3. Update the cell and extracellular fields

PhysiCelldFBA converts the objective from a specific growth rate into a cell-volume multiplier. It scales each exchange flux by the current cell dry weight and converts hours to minutes before assigning a BioFVM net export rate. Negative exchange fluxes act as uptake; positive fluxes act as secretion.

The cell then schedules its next solve and BioFVM advances the environment, closing the feedback loop.

## Supported timestep rule

!!! danger "Keep the dFBA and diffusion timesteps equal"

    The supported configuration is `intracellular_dt = dt_diffusion`. An independently coarser dFBA timestep is not currently a supported simulation mode.

For example:

```xml
<overall>
    <dt_diffusion units="min">0.01</dt_diffusion>
</overall>
...
<intracellular type="dfba">
    <settings>
        <sbml_filename>./config/model.xml</sbml_filename>
        <intracellular_dt units="min">0.01</intracellular_dt>
    </settings>
    ...
</intracellular>
```

Why synchronization matters:

- Uptake bounds are computed from the substrate available in the current voxel.
- BioFVM concentrations can change every diffusion step through diffusion and other cells.
- Reusing a metabolic solution across several diffusion steps makes its bounds stale.
- A stale uptake rate can request substrate that is no longer present, producing artifacts or violating the intended local mass balance.

If `intracellular_dt` is omitted, the C++ class defaults it to the PhysiCell diffusion timestep. Keeping the explicit values equal is still recommended because it makes the configuration auditable.

## Scheduling detail

PhysiCell checks intracellular models on diffusion updates. The dFBA object tracks `next_dfba_run` and reports whether it needs an update. After a solve, the next time is set from the configured dFBA timestep. With the supported equal-timestep rule, every metabolic cell is re-constrained and solved against the contemporaneous BioFVM field.
