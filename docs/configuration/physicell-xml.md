# Configure dFBA in PhysiCell XML

A dFBA model is attached to a cell definition by placing `<intracellular type="dfba">` inside that definition's `<phenotype>`. Each cell type can use its own model and coupling parameters.

## Complete example

```xml
<cell_definition name="ecoli" ID="0">
    <phenotype>
        <!-- cycle, death, volume, mechanics, motility, secretion, ... -->

        <intracellular type="dfba">
            <settings>
                <sbml_filename>./config/Ecoli_core.xml</sbml_filename>
                <intracellular_dt units="min">0.01</intracellular_dt>
            </settings>

            <transport_model>
                <exchange substrate="glucose">
                    <fba_flux>R_EX_glc__D_e</fba_flux>
                    <Km units="mM">0.02</Km>
                    <Vmax units="mmol/gDW/h">8.0</Vmax>
                </exchange>
                <exchange substrate="oxygen">
                    <fba_flux>R_EX_o2_e</fba_flux>
                    <Km units="mM">0.002</Km>
                    <Vmax units="mmol/gDW/h">18.0</Vmax>
                </exchange>
                <exchange substrate="CO2">
                    <fba_flux>R_EX_co2_e</fba_flux>
                    <Km units="mM">1.0</Km>
                    <Vmax units="mmol/gDW/h">0.0</Vmax>
                </exchange>
            </transport_model>

            <growth_model>
                <cell_density units="g/ml">1.04</cell_density>
                <reference_volume units="micron^3">1.3</reference_volume>
                <max_growth_rate units="1/h">0.8</max_growth_rate>
                <objective_reaction>R_BIOMASS_Ecoli_core_w_GAM</objective_reaction>
            </growth_model>

            <death_model enabled="false">
                <death_type>necrosis</death_type>
                <death_trigger_flux>R_ATPM</death_trigger_flux>
                <death_flux_threshold>8.39</death_flux_threshold>
                <death_rate_increase units="1/min">1.67e-5</death_rate_increase>
            </death_model>
        </intracellular>
    </phenotype>
</cell_definition>
```

## Define matching extracellular substrates

Every `exchange/@substrate` must already exist under `<microenvironment_setup>`:

```xml
<microenvironment_setup>
    <variable name="glucose" units="mM" ID="0">
        <physical_parameter_set>
            <diffusion_coefficient units="micron^2/min">100000</diffusion_coefficient>
            <decay_rate units="1/min">0</decay_rate>
        </physical_parameter_set>
        <initial_condition units="mM">10</initial_condition>
        <Dirichlet_boundary_condition units="mM" enabled="false">0</Dirichlet_boundary_condition>
    </variable>
</microenvironment_setup>
```

Names are case-sensitive. If `glucose` is configured in the dFBA block but only `D-Glucose` exists in BioFVM, initialization stops with a missing-substrate error.

## Keep the timesteps synchronized

Set:

```xml
<overall>
    <dt_diffusion units="min">0.01</dt_diffusion>
</overall>
```

and:

```xml
<intracellular_dt units="min">0.01</intracellular_dt>
```

to the same value. See [The dFBA update cycle](../concepts/update-cycle.md) for why a separate coarse intracellular timestep is not supported.

## Settings

`sbml_filename` is resolved relative to the directory from which the simulation is launched. Running from the repository root makes `./config/model.xml` the usual choice. Avoid machine-specific absolute paths in committed configurations.

`intracellular_dt` is retained in the schema but should equal `dt_diffusion`. If omitted, the code defaults it to the current diffusion timestep.

## Transport model

Add one `exchange` for each metabolite exchanged with BioFVM. You do not need to expose every intracellular metabolite as a diffusion field.

- `substrate` maps to the BioFVM name.
- `fba_flux` maps to an SBML reaction ID.
- `Km` controls concentration sensitivity and must be greater than zero.
- `Vmax` controls maximum uptake and must be non-negative.

The code applies both this kinetic limit and the amount physically available in the cell's voxel.

## Growth model

- `cell_density` converts solid cell volume to dry mass.
- `reference_volume` is a geometric reference used in volume/division behavior.
- `max_growth_rate` caps the metabolic objective reaction.
- `objective_reaction` identifies that reaction in the SBML model.

The objective is interpreted in h⁻¹ by the volume update. Use [the parameter reference](parameters.md) for units and constraints.

## Optional metabolic death

With `enabled="true"`, PhysiCelldFBA can increase an apoptosis or necrosis rate when the constrained metabolic problem cannot satisfy the configured trigger. The trigger reaction must exist in the SBML model. If metabolic death is not needed, omit the block or set `enabled="false"`.

## Multiple metabolic cell types

Repeat the block within each cell definition. The cross-feeding example assigns different SBML models and exchange mappings to two species while both interact through the same BioFVM fields. Reaction IDs can differ between their models; the shared contract is the extracellular substrate name.

## Configuration checklist

- [ ] `sbml_filename` is inside `<settings>` and points to a real file.
- [ ] `intracellular_dt` equals `dt_diffusion`.
- [ ] Every exchange substrate exists in `<microenvironment_setup>`.
- [ ] Every `fba_flux` and `objective_reaction` exists in the SBML model.
- [ ] `Km > 0` and `Vmax ≥ 0`.
- [ ] Flux and growth units follow the expected conventions.
- [ ] `reference_volume` is consistent with the cell's phenotype volume.
- [ ] Any death trigger refers to a valid reaction and PhysiCell death model.
