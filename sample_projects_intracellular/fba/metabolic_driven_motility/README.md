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
./ecoli-dfba-motility config/PhysiCell_settings_motility.xml
```

## Contents

- `config/Ecoli_core.xml` — core *E. coli* model used by all scenario configs.
- `config/PhysiCell_settings.xml` — default configuration.
- `config/PhysiCell_settings_motility.xml`, `PhysiCell_settings_motility_proliferation.xml`,
  `PhysiCell_settings_no_motility.xml`, `PhysiCell_settings_multiple_sources.xml` — alternate
  scenarios (motility only, motility + proliferation, proliferation without motility, multiple
  glucose sources); select by passing the file as `argv[1]`.
- `config/cells.csv` — initial cell positions, shared by all scenarios.
- `custom_modules/custom.cpp`, `custom.h` — the objective-switching logic
  (`custom_atp_optimization`), secretion-node setup (`setup_secretion_nodes`), and the migration
  rule.
- `scripts/composite_figure.py` — **generates the paper's figure** comparing representative
  scenarios (no motility/no growth, motility/no growth, motility + growth recovery); standalone
  CLI, outputs `fig_representative_scenarios.svg/png`.
- `scripts/dynamic_analysis.py` — standalone CLI; population-dynamics and phase-portrait figures
  plus a cells-over-glucose animation.
- `scripts/fluxes_analysis.py` — standalone CLI; exchange-flux plots.
- `scripts/dynamics_analysis.ipynb` — exploratory notebook; reads simulation output from
  `../output/` relative to its own location (i.e. `<repo_root>/scripts/` after
  `make metabolic_driven_motility`).

## Figure reproducibility

Manuscript Figure 7 compares two scenarios that differ only in the glucose half-saturation constant
K_m: a growth-permissive case (low K_m) and a motility-only case (higher K_m). Each runs from its own
XML, which writes to its own `<folder>`, so the two outputs can be kept side by side.

```bash
./ecoli-dfba-motility config/PhysiCell_settings_motility_proliferation.xml   # -> output_motility_proliferation_lit/
./ecoli-dfba-motility config/PhysiCell_settings_motility.xml                 # -> output_motility_lit/

python scripts/composite_figure.py \
    --scenario output_motility_proliferation_lit::"Growth-permissive (K_m 0.02 mM)" \
    --scenario output_motility_lit::"Motility only (K_m 0.5 mM)" \
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
