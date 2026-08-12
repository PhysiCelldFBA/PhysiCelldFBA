# dFBA parameter reference

This page documents the XML consumed by `dFBAIntracellular`. Numerical unit strings are descriptive; the C++ code does not perform conversions from arbitrary units.

## Settings

| Location | Parameter | Required | Expected units | Constraint/default | Effect |
| --- | --- | --- | --- | --- | --- |
| `settings` | `sbml_filename` | Yes | Path | Must exist and be readable | SBML-FBC model loaded for the cell |
| `settings` | `intracellular_dt` | Recommended explicitly | min | Defaults to `diffusion_dt`; supported rule is equality with `dt_diffusion` | Interval used for dFBA scheduling, mass limit, growth, and uptake constants |
| `settings` | `time_step` | Deprecated fallback | min | Read only when `intracellular_dt` is absent | Backward-compatible name; do not use in new files |

## Exchange mapping

Each `transport_model` must contain at least one `exchange`.

| Location | Parameter | Required | Expected units | Constraint/default | Effect |
| --- | --- | --- | --- | --- | --- |
| `exchange` attribute | `substrate` | Yes | BioFVM name | Must match exactly | Selects the extracellular density and voxel concentration |
| `exchange` | `fba_flux` | Yes | SBML reaction ID | Reaction must exist | Selects the exchange reaction whose lower bound and output flux are coupled |
| `exchange` | `Km` | Yes | mM | Must be greater than 0 | Half-saturation concentration in the kinetic uptake limit |
| `exchange` | `Vmax` | Yes | mmol/gDW/h | Must be non-negative | Maximum uptake ceiling before the local-mass limit |

The C++ runtime uses `Vmax` directly as an FBA flux bound in mmol/gDW/h. Some older XML files label the same numerical field as fmol/pgDW/min; those units differ by a factor of 60 and the parser does not convert them. Treat the runtime interpretation as authoritative and correct legacy labels/values deliberately.

## Growth model

All four fields are required.

| Location | Parameter | Required | Expected units | Constraint | Effect |
| --- | --- | --- | --- | --- | --- |
| `growth_model` | `cell_density` | Yes | g/mL = pg/µm³ | Positive physical value | Converts solid cell volume into dry mass |
| `growth_model` | `reference_volume` | Yes | µm³ | Positive and consistent with phenotype volume | Reference for cell-volume/division behavior |
| `growth_model` | `max_growth_rate` | Yes | h⁻¹ | Non-negative | Upper bound placed on `objective_reaction` |
| `growth_model` | `objective_reaction` | Yes | SBML reaction ID | Reaction must exist | Standard FBA objective interpreted as specific growth |

The implementation converts an objective in h⁻¹ to the minute-based PhysiCell update by dividing by 60.

## Metabolism-dependent death

The entire `death_model` block is optional. Omitting it or setting `enabled="false"` disables metabolic-dependent death.

| Location | Parameter | Required when enabled | Expected units | Default | Effect |
| --- | --- | --- | --- | --- | --- |
| `death_model` attribute | `enabled` | Yes | Boolean | False when omitted at the block level | Enables parsing and application of the metabolic death rule |
| `death_model` | `death_type` | No | `apoptosis` or `necrosis` | `apoptosis` | PhysiCell death process whose rate is changed |
| `death_model` | `death_trigger_flux` | Functionally yes | SBML reaction ID | Missing value disables the rule | Reaction constrained/monitored for metabolic viability |
| `death_model` | `death_flux_threshold` | No | Model flux units | `1e-6` | Lower-bound threshold applied to the trigger reaction |
| `death_model` | `death_rate_increase` | No | min⁻¹ | `0.01` | Increase applied to the selected PhysiCell death rate when metabolism is infeasible |

## Relationships that must hold

```text
intracellular_dt = dt_diffusion
exchange/@substrate ∈ BioFVM variable names
fba_flux ∈ SBML reaction IDs
objective_reaction ∈ SBML reaction IDs
Km > 0
Vmax ≥ 0
max_growth_rate ≥ 0
```

When death is enabled, `death_trigger_flux` must also exist in the model and `death_type` must be either apoptosis or necrosis.

## Parameters outside the generic dFBA block

Individual examples can add `custom_data` or `user_parameters` that are not part of the generic add-on. For example, metabolism-driven motility defines ATP allocation, maximum speed, and response-shape parameters in its custom module. Document such fields with the example that implements them rather than treating them as universal PhysiCelldFBA parameters.
