# Cancer core metabolism near a vessel

<div class="key-path">
  <strong>Directory:</strong> <code>sample_projects_intracellular/fba/cancer_core_metabolism/</code><br>
  <strong>Make target:</strong> <code>cancer-core-metabolism</code><br>
  <strong>Executable:</strong> <code>cancer-core-metabolism</code>
</div>

## Biological question

How does distance from a nutrient source organize metabolic state in a tumour-like cell population? MCF7 cells carrying the same reconstruction are arranged next to a fixed vessel that supplies oxygen and nutrients.

## Run

```bash
make cancer-core-metabolism
make
./cancer-core-metabolism config/PhysiCell_settings_wide_domain.xml
```

The directory also includes default and tiny configurations. Use the tiny setup to check model loading and mappings before a manuscript-scale simulation.

## Important files

| File | Role |
| --- | --- |
| `config/MCF7_Zielinski_2017_core.sbml` | MCF7-specific metabolic reconstruction |
| `config/MCF7_Zielinski_2017_core_trimmed.sbml` | Reduced reconstruction |
| `config/dfba_model.yaml` | Substrate/model mapping source |
| `config/PhysiCell_settings_tiny.xml` | Small test configuration |
| `config/PhysiCell_settings_wide_domain.xml` | Vessel-distance scenario |
| `scripts/plot_cancer_tissue_spatial.py` | Spatial analysis |
| `scripts/plot_exchange_profiles.py` | Distance-dependent flux profiles |
| `scripts/cancer_tissue_figures.ipynb` | Figure workflow |

The model is constrained by MCF7 exometabolomic measurements and exchanges a broad metabolite panel, including oxygen, glucose, lactate, amino acids, ammonium, phosphate, ornithine, and CO₂.

## Reported result

Distance from the vessel produces three zones:

| Approximate vessel distance | Reported state |
| --- | --- |
| 0–80 µm | Proliferative cells with active biomass production |
| 80–160 µm | Hypoxic transition with progressively reduced growth |
| Beyond 160 µm | Necrotic core where metabolic fluxes fall below maintenance |

The zonation appears coherently across multiple uptake and secretion fluxes, not only oxygen and glucose. Even within the proliferative zone, individual cells can show different exchange activity.

<div class="figure-placeholder">
  <strong>Figure placeholder — cancer metabolic zonation</strong>
  Vessel geometry, cell state/growth, and multi-metabolite flux profiles versus distance
  <code>docs/assets/figures/final/cancer-metabolic-zonation.png</code>
</div>

## Interpretation

The zones are not assigned with distance thresholds. They emerge from independently diffusing nutrients constraining a shared intracellular network. This example is therefore a template for coupling context-specific mammalian reconstructions to tissue-scale transport.

Transport parameters for metabolites without measured kinetics are phenomenological. Do not interpret individual flux magnitudes as patient-specific predictions without additional calibration and validation.
