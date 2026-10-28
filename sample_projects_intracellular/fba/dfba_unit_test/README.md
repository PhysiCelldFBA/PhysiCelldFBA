# dfba_unit_test

A minimal single/few-agent dFBA setup used as a regression/unit test for the dFBA intracellular
module, not a biological case study. It checks four things:

1. **Mass conversion** — a single agent on the *E. coli* core model (`config/Ecoli_core.xml`),
   run with the default `config/PhysiCell_settings.xml` (no args needed — this scenario IS the
   default config), analyzed in `scripts/dfba_analysis.ipynb`, which renders the mass-conversion
   figure used in the paper (`dfba_unit_test_mass_conversion.svg`).
2. **Voxel-sharing correctness** — several agents drawing on the same microenvironment voxel must
   not drive its substrate negative, and an agent that cannot take up a substrate must not dilute
   the uptake available to others. Driven by `scripts/run_voxel_sharing_test.sh`, which runs each
   of `config/voxel_sharing_test/Unit_test_voxel_sharing_{n1,n4,n8,mixed}.xml` and checks the
   output with `scripts/check_voxel_nonnegative.py` and `scripts/check_voxel_draw_equal.py`.
3. **Update order and thread count** — twelve cells of four types in one voxel, created in
   different orders and run on 1 and 4 threads, must give the same result. Also run by
   `scripts/run_voxel_sharing_test.sh`: `scripts/gen_update_order_configs.py` derives the configs
   from `config/voxel_sharing_test/Unit_test_voxel_sharing_mixed.xml` and
   `scripts/check_update_order.py` compares the runs.
4. **Voxel size and time step** — one *E. coli* core cell in a closed 16 µm box, run at
   dx ∈ {0.8, 1.6, 3.2} µm and dt ∈ {0.005, 0.01, 0.02} min (the diffusion and the dFBA time step
   kept equal). While the uptake cap is not limiting, glucose consumption, cell volume and growth
   rate must agree across dx and dt; with a non-diffusing glucose pool the cap must act and empty the
   cell's voxel without overdrawing it. The same runs check the per-step conservation ledger, the
   carbon balance and the elemental balance of the model. Driven by `scripts/run_dxdt_test.sh`, which
   builds its configs from the default `config/PhysiCell_settings.xml`
   (`scripts/gen_dxdt_configs.py`) and checks them with `scripts/check_dxdt.py` and
   `scripts/elemental_check.py`.

## Build & run

From the repository root:
```bash
make dfba_unit_test
make
./dfba-unit-test                                          # mass-conversion test (default config)
bash scripts/run_voxel_sharing_test.sh ./dfba-unit-test   # voxel-sharing and update-order tests
bash scripts/run_dxdt_test.sh ./dfba-unit-test            # voxel-size and time-step test
```

The check scripts need Python with `numpy` and `scipy`; `scripts/check_dxdt.py` and
`scripts/elemental_check.py` also need `cobra`. Set `PYTHON` to choose the interpreter. Each test
writes its runs to `output_*` folders in the repository root and moves earlier results aside instead
of deleting them. `run_dxdt_test.sh` runs its single-thread cases in parallel (`JOBS`, default 4).

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
- `custom_modules/dfba_conservation_ledger.h` — per-step conservation ledger of the extracellular
  substrates, called from `main.cpp`; each run writes `conservation_ledger.csv` to its output folder.
- `scripts/dfba_analysis.ipynb` — plots mass-conversion results and saves the paper figure.
- `scripts/plot_df_journal_style_ext.py` — shared plotting style helper (also used by
  `ecoli_acetic_switch`).
- `scripts/check_voxel_nonnegative.py`, `scripts/check_voxel_draw_equal.py` — pass/fail checks for
  the voxel-sharing scenarios.
- `scripts/gen_update_order_configs.py`, `scripts/check_update_order.py` — configs and pass/fail
  check of the update-order and thread-count test (run from `run_voxel_sharing_test.sh`).
- `scripts/run_dxdt_test.sh`, `scripts/gen_dxdt_configs.py`, `scripts/check_dxdt.py`,
  `scripts/elemental_check.py` — the voxel-size and time-step test: runner, config generator,
  pass/fail checks, and the elemental-balance check of the SBML model.

## Figure reproducibility

| Manuscript figure | Panels | Script | XML input | Expected output |
|---|---|---|---|---|
| Main Figure 2 (mass-conservation validation) | A: glucose and CO₂ over time; B: biomass against the analytical prediction, with cell count | `scripts/dfba_analysis.ipynb` | `config/PhysiCell_settings.xml` (default; no arguments) | `<output_folder>/dfba_unit_test_mass_conversion.svg` |

## Citation

See `ALL_CITATIONS.txt`.
