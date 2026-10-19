# bacterial_colony

A dFBA-driven *E. coli* colony growing on a shared glucose substrate, built on the
`iML1515_trimmed` genome-scale metabolic model. Each agent runs its own flux-balance
optimization every `post_update_intracellular` step, so colony-scale growth and substrate
depletion emerge from per-cell metabolism rather than a prescribed growth law.

## Contents

- `config/Ecoli_core.xml` — core *E. coli* metabolic model (SBML) used by the dFBA intracellular model.
- `config/PhysiCell_settings.xml` — main simulation configuration; sets the microenvironment, cell
  definitions, and the `sbml_filename` loaded by each agent.
- `config/dfba_config.yaml` — standalone dFBA solver configuration kept for reference; not currently
  wired into this project's XML pipeline (the model path is set directly via `sbml_filename` in
  `PhysiCell_settings.xml`).
- `custom_modules/custom.cpp`, `custom.h` — cell type setup (`create_cell_types`), initial colony
  layout (`setup_tissue`, `create_cell_disc_positions`), and per-step dFBA coupling
  (`post_update_intracellular`).
- `scripts/bacterial_colony_figures.ipynb` — notebook that reproduces the colony figures from a
  simulation output folder.
- `scripts/plot_bacteria_colony.py` — standalone plotting helper used by the notebook.
- `scripts/animate_svg.sh` — stitches per-timestep SVG snapshots into an animation.

## Build & run

From this directory (standalone):
```bash
make
./bacterial_colony
```

Or from the PhysiCelldFBA repository root, to populate the top-level build:
```bash
make bacterial-colony
make
./bacterial-colony
```

Output (SVG snapshots, MultiCellDS `.mat`/`.xml` frames) is written to `./output/`.

## Citation

See `ALL_CITATIONS.txt` for the PhysiCell core citation. This project additionally depends on the
`iML1515` genome-scale *E. coli* model (Monk et al., 2017).
