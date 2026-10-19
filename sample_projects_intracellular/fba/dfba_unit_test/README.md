# dfba_unit_test

A minimal single/few-agent dFBA setup used as a regression/unit test for the dFBA intracellular
module, not a biological case study. It checks two things:

1. **Mass conversion** — a single agent on the *E. coli* core model (`config/Ecoli_core.xml`),
   validated against `config/Ecoli_core_noGAM.xml` (no growth-associated maintenance) and
   `config/Unit_test_mass_conversion.xml`, analyzed in `scripts/dfba_analysis.ipynb`.
2. **Voxel-sharing correctness** — several agents drawing on the same microenvironment voxel must
   not drive its substrate negative, and an agent that cannot take up a substrate must not dilute
   the uptake available to others. Driven by `scripts/run_voxel_sharing_test.sh`, which runs each
   of `config/Unit_test_voxel_sharing_{n1,n4,n8,mixed}.xml` and checks the output with
   `scripts/check_voxel_nonnegative.py` and `scripts/check_voxel_draw_equal.py`.

## Contents

- `config/PhysiCell_settings.xml` — default settings (used by plain `./dfba-unit-test` with no args).
- `config/Unit_test_voxel_sharing_{n1,n4,n8,mixed}.xml` — voxel-sharing scenarios (n agents per
  voxel, or a mixed layout), passed explicitly to the binary by `run_voxel_sharing_test.sh`.
- `config/Unit_test_mass_conversion.xml` — mass-conversion test scenario.
- `config/Ecoli_core.xml`, `config/Ecoli_core_noGAM.xml` — the two metabolic model variants
  compared for mass-conversion correctness.
- `custom_modules/custom.cpp`, `custom.h` — single-agent dFBA coupling for the unit test.
- `scripts/dfba_analysis.ipynb` — plots mass-conversion results.
- `scripts/plot_df_journal_style_ext.py` — shared plotting style helper (also used by
  `ecoli_acetic_switch`).
- `scripts/check_voxel_nonnegative.py`, `scripts/check_voxel_draw_equal.py` — pass/fail checks for
  the voxel-sharing scenarios.

## Build & run

From this directory:
```bash
make
./dfba-unit-test                                    # default config
./dfba-unit-test config/Unit_test_mass_conversion.xml
bash scripts/run_voxel_sharing_test.sh ./dfba-unit-test   # voxel-sharing suite
```

Or from the repository root: `make dfba_unit_test && make`.

## Citation

See `ALL_CITATIONS.txt`.
