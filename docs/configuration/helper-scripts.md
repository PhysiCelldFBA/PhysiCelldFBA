# Configuration helper scripts

Two optional Python utilities turn an SBML model into a starting PhysiCelldFBA configuration:

```text
SBML model → generate_dfba_yaml.py → YAML draft
YAML draft → dfba_configurator.py → new PhysiCell XML
```

These tools accelerate repetitive mapping work. They do not replace biological review of units, kinetics, medium composition, model reduction, or boundary conditions.

## Install Python dependencies

The scripts import COBRApy, pandas, PyYAML, lxml, and `physicell-settings`. Install them in an isolated Python environment according to your platform and project dependency policy.

## Generate a YAML draft

Current command-line form:

```bash
python3 beta/generate_dfba_yaml.py \
    path/to/model.xml \
    config/dfba_model.yaml \
    --cell_type my_cell \
    --cell_volume 2494 \
    --ex_prefix R_EX
```

The two positional arguments are the input SBML model and output YAML file. Optional arguments control the cell-type name, reference volume, and exchange-reaction prefix.

The script:

1. loads the model with COBRApy;
2. selects the first reaction with a non-zero objective coefficient;
3. solves the model and records the optimal objective as `max_growth_rate`;
4. runs flux-variability analysis;
5. removes blocked proton/water exchanges;
6. creates BioFVM substrate and dFBA exchange entries;
7. writes a YAML configuration.

Defaults such as 5 mM initial concentration, 50,000 µm²/min diffusion, `Km = 0.001 mM`, and names derived from SBML metabolite names are only starting values.

!!! warning "Review generated Vmax values"

    Generated `Vmax` values are derived from flux-variability bounds, not measured transporter kinetics. Confirm the units and replace them with defensible values before interpreting a simulation biologically.

## Generate a PhysiCell XML

Current command-line form:

```bash
python3 beta/dfba_configurator.py \
    --config config/dfba_model.yaml \
    --sbml-folder . \
    --output config/PhysiCell_settings_generated.xml \
    --cell-prefix dfba_
```

The configurator validates required YAML sections and model paths, creates BioFVM substrates and cell definitions through `physicell-settings`, injects a dFBA intracellular block into each generated cell definition, and writes a new XML file.

!!! note "The current configurator creates a new configuration"

    Despite its name, the present script does not accept an arbitrary existing PhysiCell XML and update it in place. It starts from a new `PhysiCellConfig`. If your project already contains mechanics, cycle, boundary, or custom settings, transfer or merge those settings deliberately after generation.

## YAML shape

```yaml
substrates:
  - name: glucose
    diffusion_coefficient: 100000
    decay_rate: 0
    initial_condition: 10

models:
  ecoli:
    settings:
      sbml_path: config/Ecoli_core.xml
      intracellular_dt: 0.01
    growth_model:
      cell_density: 1.04
      reference_volume: 1.3
      nuclear_volume: 0
      max_growth_rate: 0.8
      objective_reaction: R_BIOMASS_Ecoli_core_w_GAM
    death_model:
      enabled: false
    exchanges:
      - substrate: glucose
        fba_flux: R_EX_glc__D_e
        Km: 0.02
        Vmax: 8.0
```

Keep `intracellular_dt` equal to the diffusion timestep in the generated XML. The YAML generator does not inspect a target PhysiCell configuration to enforce that rule.

Also review generated unit attributes: the current runtime interprets `max_growth_rate` as h⁻¹ and `Vmax` as mmol/gDW/h, regardless of the strings written into the XML.

## Validate the result

Before a full simulation:

1. Open the generated XML and verify model paths.
2. Check every substrate/reaction mapping.
3. Correct units and parameter values.
4. Add project-specific domain, boundary, mechanics, cycle, save, and initial-condition settings.
5. Run one cell in a small domain.
6. Compare its objective and exchange fluxes with an independent FBA solve.
7. Confirm local mass conservation.
