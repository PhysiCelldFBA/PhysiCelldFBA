# dfba_unit_test

A minimal single/few-agent dFBA setup used as a regression/unit test for the dFBA intracellular
module, not a biological case study. It checks two things:

1. **Mass conversion** — a single agent on the *E. coli* core model (`config/Ecoli_core.xml`),
   run with the default `config/PhysiCell_settings.xml` (no args needed — this scenario IS the
   default config), analyzed in `scripts/dfba_analysis.ipynb`, which renders the mass-conversion
   figure used in the paper (`dfba_unit_test_mass_conversion.svg`).
2. **Voxel-sharing correctness** — several agents drawing on the same microenvironment voxel must
   not drive its substrate negative, and an agent that cannot take up a substrate must not dilute
   the uptake available to others. Driven by `scripts/run_voxel_sharing_test.sh`, which runs each
   of `config/voxel_sharing_test/Unit_test_voxel_sharing_{n1,n4,n8,mixed}.xml` and checks the
   output with `scripts/check_voxel_nonnegative.py` and `scripts/check_voxel_draw_equal.py`.

## Build & run

From the repository root:
```bash
make dfba_unit_test
make
./dfba-unit-test                                    # mass-conversion test (default config)
bash scripts/run_voxel_sharing_test.sh ./dfba-unit-test   # voxel-sharing suite
```

## Contents

- `config/PhysiCell_settings.xml` — default settings (used by plain `./dfba-unit-test` with no
  args); this is the mass-conversion test scenario: a 16 × 16 × 16 µm domain with 2 µm voxels, all
  time steps 0.001 min (the grid and time steps of the other *E. coli* samples), 240 min.
- `config/voxel_sharing_test/Unit_test_voxel_sharing_{n1,n4,n8,mixed}.xml` — voxel-sharing
  scenarios (n agents per voxel, or a mixed layout), passed explicitly to the binary by
  `run_voxel_sharing_test.sh`.
- `config/Ecoli_core.xml` — the *E. coli* core metabolic model used by every scenario above.
- `config/Ecoli_core.yaml` — dFBA definition (substrates, cell definition and dFBA block) of the
  default `config/PhysiCell_settings.xml`; input for `addons/dFBA/utils/dfba_configurator.py` (see
  `addons/dFBA/utils/README.md`). The voxel-sharing XMLs have no YAML of their own.
- `custom_modules/custom.cpp`, `custom.h` — single-agent dFBA coupling for the unit test.
- `scripts/dfba_analysis.ipynb` — plots mass-conversion results and saves the paper figure.
- `scripts/plot_df_journal_style_ext.py` — shared plotting style helper (also used by
  `ecoli_acetic_switch`).
- `scripts/check_voxel_nonnegative.py`, `scripts/check_voxel_draw_equal.py` — pass/fail checks for
  the voxel-sharing scenarios.

## Figure reproducibility

| Manuscript figure | Panels | Script | XML input | Expected output |
|---|---|---|---|---|
| Main Figure 2 (mass-conservation validation) | A: glucose and CO₂ over time; B: biomass against the analytical prediction, with cell count | `scripts/dfba_analysis.ipynb` | `config/PhysiCell_settings.xml` (default; no arguments) | `<output_folder>/dfba_unit_test_mass_conversion.svg` |

## Citation

See `ALL_CITATIONS.txt`.
