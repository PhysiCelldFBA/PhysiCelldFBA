# Metabolism-driven motility

<div class="key-path">
  <strong>Directory:</strong> <code>sample_projects_intracellular/fba/metabolic_driven_motility/</code><br>
  <strong>Make target:</strong> <code>metabolic_driven_motility</code><br>
  <strong>Executable:</strong> <code>ecoli-dfba-motility</code>
</div>

## Biological question

Can a metabolic model regulate a cell behavior other than growth? This example connects ATP-production capacity to migration speed and uses a hysteretic glucose rule to switch between motility and biomass objectives.

## Run a scenario

```bash
make metabolic_driven_motility
make
./ecoli-dfba-motility config/PhysiCell_settings_motility_proliferation.xml
```

Available configurations include:

| File | Intended comparison |
| --- | --- |
| `PhysiCell_settings_no_motility.xml` | Growth without metabolic motility |
| `PhysiCell_settings_motility.xml` | Motility-focused setup |
| `PhysiCell_settings_motility_proliferation.xml` | Migration followed by growth recovery |
| `PhysiCell_settings_multiple_sources.xml` | Multiple nutrient sources |

## Custom objective

The generic dFBA objective is biomass. The example's custom module can instead call:

```cpp
dfba->optimize_for_objective("R_ATPM", 1.0);
```

It subtracts a basal ATP requirement, maps the remaining ATP capacity to a speed fraction, and sets PhysiCell's migration speed. When the cell switches back to growth mode, the normal biomass objective is restored.

Project-specific user parameters include:

- `ecoli_vmax`: maximum migration speed;
- `basal_atp_flux`: ATP reserved for maintenance;
- `motility_cost_at_vmax`: ATP surplus required for maximum speed;
- `phi_atp_hill`: response-shape exponent.

These are custom-module parameters, not generic `dFBAIntracellular` fields.

## Reported regimes

The interaction among glucose availability, uptake affinity, and ATP capacity produces:

1. **Migration and growth recovery:** low glucose `Km` permits efficient uptake; the cell moves toward glucose and switches back to biomass optimization after crossing the upper threshold.
2. **Persistent motility:** a higher `Km` provides enough ATP for motion but does not reach the glucose-flux condition needed for growth mode.
3. **Metabolic quiescence:** low glucose supply combined with weak uptake produces little ATP, minimal motion, and no meaningful growth.

<div class="figure-placeholder">
  <strong>Figure placeholder — metabolism-driven motility</strong>
  Cell trajectories, migration speed, growth, ATP production, motility demand, and glucose uptake
  <code>docs/assets/figures/final/metabolism-driven-motility.png</code>
</div>

## Extending the pattern

The same mechanism can connect a metabolic objective or flux to another phenotype:

1. implement a project-specific pre/post intracellular hook or custom optimization;
2. solve the alternative objective without permanently losing the original objective;
3. normalize the metabolic signal to a well-defined behavior range;
4. store diagnostic values in `custom_data`;
5. verify energy and mass accounting when the behavior consumes a metabolic resource.

This example is qualitative rather than a detailed chemotaxis model. Calibrate ATP costs, thresholds, and migration laws before making quantitative biological claims.
