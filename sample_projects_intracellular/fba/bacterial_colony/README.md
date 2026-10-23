# bacterial_colony

A dFBA-driven *E. coli* colony growing on a shared glucose substrate, built on the *E. coli* core
metabolic model (`config/Ecoli_core.xml`). Each agent runs its own flux-balance optimization every
`post_update_intracellular` step, so colony-scale growth and substrate depletion emerge from
per-cell metabolism rather than a pre-defined growth law. Cells are placed procedurally in a
close-packed disc (`setup_tissue` → `create_cell_disc_positions`), not loaded from a CSV.

## Build & run

From the PhysiCelldFBA repository root:
```bash
make bacterial-colony
make
./bacterial-colony
```

Output (SVG snapshots, MultiCellDS `.mat`/`.xml` frames) is written to `./output/`.

## Contents

- `config/Ecoli_core.xml` — core *E. coli* metabolic model (SBML), loaded via `sbml_filename` in
  `PhysiCell_settings.xml`.
- `config/PhysiCell_settings.xml` — main simulation configuration; sets the microenvironment, cell
  definitions, and the `sbml_filename` loaded by each agent.
- `custom_modules/custom.cpp`, `custom.h` — cell type setup (`create_cell_types`), initial colony
  layout (`setup_tissue`, `create_cell_disc_positions`), and per-step dFBA coupling
  (`post_update_intracellular`).
- `scripts/bacterial_colony_figures.ipynb` — **generates the paper's Figure 4** end-to-end from a
  simulation output folder (per the notebook's own header).
- `scripts/animate_svg.sh` — standalone utility; stitches per-timestep SVG snapshots into an MP4
  animation via inkscape/ffmpeg.

## Figure reproducibility

| Manuscript figure | Panels | Script | XML input | Expected output |
|---|---|---|---|---|
| Main Figure 4A (total biomass over time) | Linear y-axis | `scripts/bacterial_colony_figures.ipynb` | `config/PhysiCell_settings.xml` | `fig_folder/panel_A_biomass.{svg,png}` |
| Supplementary Figure S6 (exponential-to-sub-exponential transition) | Same data, log y-axis, dashed line at t* | same notebook | same | `fig_folder/panel_A_biomass_log.{svg,png}` |
| Main Figure 4B (extracellular O₂, glucose, acetate, CO₂ at 9 h) | Spatial fields | same notebook | same | `fig_folder/panel_B_substrates.{svg,png}` |
| Main Figure 4C (growth rate, 1/h) | Single-cell spatial map | same notebook | same | `fig_folder/panel_C_growth_rate.{svg,png}` |
| Main Figure 4D (glucose uptake) | Single-cell spatial map | same notebook | same | `fig_folder/panel_D_glucose_uptake.{svg,png}` |
| Main Figure 4E (acetate exchange) | Single-cell spatial map | same notebook | same | `fig_folder/panel_E_acetate_exchange.{svg,png}` |

The notebook also writes `panel_A_biomass_cellcount_linear` and `panel_A_growth_transition_diagnostic`.
Neither appears as a separate panel in the manuscript; the diagnostic contains the breakpoint fit
described in the Supplementary Figure S6 legend. `fig_folder/` is created relative to the notebook
kernel's working directory. All panels come from one simulation run and one notebook pass.

## Citation

See `ALL_CITATIONS.txt` for the PhysiCell core citation.
