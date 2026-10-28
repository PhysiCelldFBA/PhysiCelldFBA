# metabolic_driven_motility

Demonstrates dFBA-driven chemotaxis-like behavior: each *E. coli* agent (core model,
`config/Ecoli_core.xml`) can switch its FBA objective from biomass production to ATP production
(`custom_atp_optimization`), letting a cell trade growth for directed migration toward a
glucose-richer region of the microenvironment (`metabolic_bound_migration_rule`).

## Build & run

From the repository root:
```bash
make metabolic_driven_motility
make
./ecoli-dfba-motility config/PhysiCell_settings_motility_proliferation.xml   # growth-permissive -> output_motility_prolif/
./ecoli-dfba-motility config/PhysiCell_settings_motility.xml                 # motility only     -> output_motility/
```

The project ships two scenario configs, which differ in the glucose K_m (0.02 and 0.5 mM),
`glucose_threshold_low` and `motility_cost_at_vmax`. Always pass one of them as `argv[1]`: without an
argument the program loads the repository's own `config/PhysiCell_settings.xml`, not a scenario of
this project. SVG snapshots are enabled in both.

## Contents

- `config/Ecoli_core.xml` — core *E. coli* model used by all scenario configs.
- `config/PhysiCell_settings_motility_proliferation.xml` — growth-permissive scenario (motility +
  proliferation, glucose K_m 0.02 mM); writes to `output_motility_prolif/`.
- `config/PhysiCell_settings_motility.xml` — motility-only scenario (glucose K_m 0.5 mM); writes to
  `output_motility/`.
- `config/Ecoli_core_motility_proliferation.yaml`, `config/Ecoli_core_motility.yaml` — dFBA definitions
  of the two scenario configs above; inputs for `addons/dFBA/utils/dfba_configurator.py` (see
  `addons/dFBA/utils/README.md`).
- `config/cells.csv` — initial cell positions, shared by all scenarios.
- `custom_modules/custom.cpp`, `custom.h` — the objective-switching logic
  (`custom_atp_optimization`), secretion-node setup (`setup_secretion_nodes`), and the migration
  rule.
- `scripts/composite_figure.py` — **generates the paper's figure** comparing the two scenarios
  (growth-permissive and motility-only); standalone CLI, outputs
  `fig_representative_scenarios.svg/png`.
- `scripts/dynamic_analysis.py` — standalone CLI; population-dynamics and phase-portrait figures
  plus a cells-over-glucose animation. Reads `output_motility_prolif/` by default; pass
  `--output-dir` for the other scenario.
- `scripts/fluxes_analysis.py` — standalone CLI; exchange-flux plots.
- `scripts/dynamics_analysis.ipynb` — exploratory notebook; reads simulation output from
  `../output_motility_prolif/` relative to its own location (i.e. `<repo_root>/scripts/` after
  `make metabolic_driven_motility`).

## Figure reproducibility

Manuscript Figure 7 compares two scenarios that differ only in the glucose half-saturation constant
K_m: a growth-permissive case (low K_m) and a motility-only case (higher K_m). Each runs from its own
XML, which writes to its own `<folder>`, so the two outputs can be kept side by side.

```bash
./ecoli-dfba-motility config/PhysiCell_settings_motility_proliferation.xml   # -> output_motility_prolif/
./ecoli-dfba-motility config/PhysiCell_settings_motility.xml                 # -> output_motility/

python scripts/composite_figure.py \
    --scenario output_motility_prolif::"Growth-permissive (K_m 0.02 mM)" \
    --scenario output_motility::"Motility only (K_m 0.5 mM)" \
    --results-dir results/representative_scenarios
```

| Manuscript figure | Panels | Script | XML inputs | Expected output |
|---|---|---|---|---|
| Main Figure 7 (parameter-dependent motility and proliferation under glucose limitation) | A–B: cell trajectories; C–D: migration speed and growth rate; E–F: ATP and glucose flux (growth-permissive in the left column, motility-only in the right) | `scripts/composite_figure.py` | `config/PhysiCell_settings_motility_proliferation.xml` (A, C, E); `config/PhysiCell_settings_motility.xml` (B, D, F) | `results/representative_scenarios/fig_representative_scenarios.{svg,png}`, `representative_scenarios_metrics.csv` |

Both simulations use `phi_atp_hill = 1` (the Hill coefficient stated in the Figure 7 legend). The
`--runs-dir` auto-discovery in `composite_figure.py` only matches the author's `hill*_job*/output`
sweep naming, so the scenarios are passed explicitly as above.

## Citation

See `ALL_CITATIONS.txt`.
