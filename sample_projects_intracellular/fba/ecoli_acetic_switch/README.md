# ecoli_acetic_switch

A dFBA model of the *E. coli* overflow/acetate-switch behavior: agents run flux balance analysis
against the core *E. coli* model (`config/Ecoli_core.xml`) or the genome-scale `iML1515.xml`,
and switch between glucose and acetate uptake as glucose is depleted
(`reintroduce_nutrients_function`, `inject_density` handle periodic nutrient refeeding for
longer runs).

## Contents

- `config/Ecoli_core.xml` — core model, used by the default and `*_static_condition*` configs.
- `config/iML1515.xml`, `config/Ecoli_large_scale.xml` — genome-scale model and its large-scale
  scenario config, used by `Ecoli_large_scale.xml`'s `sbml_filename` for bigger simulations.
- `config/PhysiCell_settings.xml` — default configuration.
- `config/PhysiCell_settings_static_condition.xml` /
  `PhysiCell_settings_static_condition_literature_params.xml` — fixed-condition scenarios with
  literature-derived parameters; pass as `argv[1]` to select.
- `config/cells_reduced.csv`, `config/cells_large_scale.csv` — initial positions for the reduced
  (core model) and large-scale (iML1515) scenarios respectively.
- `custom_modules/custom.cpp`, `custom.h` — cell setup and the acetate-switch coupling logic.
- `scripts/ecoli_acetate_switch.ipynb`, `dfba_analysis.py`, `fluxes_analysis.py` — flux/behavior
  analysis.
- `scripts/dFBA_analysis_large_scale.ipynb` — analysis for the `iML1515` large-scale scenario.
- `scripts/plot_ecoli_figure_a.py`, `plot_ecoli_figure_b.py` — manuscript figure generation.
- `scripts/plot_df_journal_style_ext.py` — shared plotting style helper (also used by
  `dfba_unit_test`).

## Build & run

From this directory:
```bash
make
./ecoli-dfba                                              # default (core model)
./ecoli-dfba config/PhysiCell_settings_static_condition.xml
./ecoli-dfba config/Ecoli_large_scale.xml                 # genome-scale iML1515 scenario
```

Or from the repository root: `make ecoli-acetic-switch-sample && make`.

## Citation

See `ALL_CITATIONS.txt`. The genome-scale model is `iML1515` (Monk JM et al., "iML1515, a
knowledgebase that computes Escherichia coli traits", Nat. Biotechnol. 35, 904-908, 2017).
