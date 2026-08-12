# E. coli acetate switch

<div class="key-path">
  <strong>Directory:</strong> <code>sample_projects_intracellular/fba/ecoli_acetic_switch/</code><br>
  <strong>Make target:</strong> <code>ecoli-acetic-switch-sample</code><br>
  <strong>Executable:</strong> <code>ecoli-dfba</code>
</div>

## Biological question

Can repeated metabolic optimization reproduce the *E. coli* acetate switch without an explicit regulatory rule? The model must move from aerobic glucose growth to acetate overflow and finally to acetate consumption after glucose depletion.

## Run

```bash
make ecoli-acetic-switch-sample
make
./ecoli-dfba config/PhysiCell_settings_static_condition.xml
```

Additional XML files explore other and larger-scale conditions. Always inspect the selected model path, initial-cell CSV, output folder, and timestep values before launching a long run. The current default `PhysiCell_settings.xml` contains an obsolete machine-specific SBML path and dFBA XML layout, so use the explicit portable configuration above until that file is corrected.

## Important files

| File | Role |
| --- | --- |
| `config/Ecoli_core.xml` | Core *E. coli* metabolic reconstruction |
| `config/iML1515.xml` | Genome-scale alternative |
| `config/PhysiCell_settings_static_condition.xml` | Portable acetate-switch setup |
| `config/PhysiCell_settings.xml` | Older default; currently requires correction before use |
| `scripts/ecoli_acetate_switch.ipynb` | Analysis notebook |
| `scripts/dfba_analysis.py` and `fluxes_analysis.py` | Output/flux analysis |

The dFBA block couples glucose, oxygen, acetate, and CO₂. Substrate-specific `Km` and `Vmax` values give glucose and oxygen high uptake capacity, acetate an intermediate uptake capacity, and CO₂ a zero uptake capacity while permitting secretion.

## Reported result

The model produces three phases:

1. **Aerobic glucose growth:** rapid glucose consumption and negligible acetate production.
2. **Oxygen-limited overflow:** glucose use continues while acetate is secreted.
3. **Acetate reutilization:** after glucose exhaustion, acetate becomes the carbon source and slower growth continues.

The switch follows from the feasible flux space under the evolving glucose, oxygen, and acetate fields. No separate acetate-switch rule is imposed.

<div class="figure-placeholder">
  <strong>Figure placeholder — E. coli acetate switch</strong>
  Glucose, oxygen, acetate, CO₂, and biomass across the three metabolic phases
  <code>docs/assets/figures/final/ecoli-acetate-switch.png</code>
</div>

## Interpretation

This example demonstrates temporal metabolic adaptation in a changing but shared environment. It is a useful template for diauxic growth, overflow metabolism, and by-product reuse.

When changing kinetic parameters, distinguish between:

- a model that cannot use acetate because its SBML bounds/pathways forbid it;
- uptake that is prevented by `Vmax = 0`;
- uptake that is weak because local acetate is much smaller than `Km`;
- a feasible acetate pathway that is not selected by the current objective.
