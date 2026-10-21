# bacterial_colony

A dFBA-driven *E. coli* colony growing on a shared glucose substrate, built on the *E. coli* core
metabolic model (`config/Ecoli_core.xml`). Each agent runs its own flux-balance optimization every
`post_update_intracellular` step, so colony-scale growth and substrate depletion emerge from
per-cell metabolism rather than a pre-defined growth law. Cells are placed procedurally in a
close-packed disc (`setup_tissue` → `create_cell_disc_positions`), not loaded from a CSV.

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

## Build & run

From this directory (standalone):
```bash
make
./bacterial-colony
```

Or from the PhysiCelldFBA repository root, to populate the top-level build:
```bash
make bacterial-colony
make
./bacterial-colony
```

Output (SVG snapshots, MultiCellDS `.mat`/`.xml` frames) is written to `./output/`.

## Citation

See `ALL_CITATIONS.txt` for the PhysiCell core citation.
