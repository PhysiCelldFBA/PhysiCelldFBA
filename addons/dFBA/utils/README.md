# utils

Helper scripts for building dFBA configuration files for PhysiCelldFBA. They are not part of the
PhysiCell build; the sample projects under `sample_projects_intracellular/fba/` are the source of
truth for each run.

| Script | Purpose |
|---|---|
| `generate_dfba_yaml.py` | Reads an SBML model and writes a dFBA YAML definition (substrates, exchanges, growth and death parameters). |
| `dfba_configurator.py` | Reads a dFBA YAML definition and writes a PhysiCell XML configuration with the dFBA blocks. |

## Requirements

- `generate_dfba_yaml.py`: Python 3, `cobra`, `pandas`, `pyyaml`.
- `dfba_configurator.py`: Python 3, `lxml`, `pandas`, `pyyaml`, and `physicell-settings` (the script imports `physicell_config`).

## Generate a YAML from an SBML model

```bash
python addons/dFBA/utils/generate_dfba_yaml.py <model.xml> <output.yaml> --cell_type <name> --cell_volume <um^3> [--ex_prefix R_EX]
```

The generator's substrate list, diffusion coefficients and initial concentrations are generic
defaults. Reconcile them with the microenvironment of the project before using the file.

## Build a PhysiCell XML from a YAML

```bash
python addons/dFBA/utils/dfba_configurator.py --output <config.xml> --config <models.yaml> --sbml-folder <sbml_dir> [--cell-prefix <prefix>]
```

`--sbml-folder` resolves the relative `sbml_path` entries in the YAML. The YAML can contain one or
several models; each model key becomes a dFBA-enabled cell definition.

The configurator writes a new PhysiCell XML containing only the substrates, the cell definitions
(with their volumes) and the dFBA block. The rest of the settings (domain, motility, rules) are not
generated and must be set by hand.

### Known limitations

- `sbml_path` must be given under `settings:` in each model (not directly under the model key).
- The `physicell-settings` library requires a positive `nuclear_volume`. Prokaryote models use 0, so
  give a small positive placeholder in the YAML and set the nuclear volume to 0 in the XML afterwards.
- The XML it writes labels `Vmax` as `fmol/pg DW cell/min` and `max_growth_rate` as `1/min`. The
  dFBA code stores these unit strings but does not use them in any calculation, so the numbers in
  the YAML are used as given (the sample project XMLs label them `mmol/gDW cell/h` and `1/h`).
