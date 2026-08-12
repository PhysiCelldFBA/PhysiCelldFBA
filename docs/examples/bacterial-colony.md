# Growing E. coli colony

<div class="key-path">
  <strong>Directory:</strong> <code>sample_projects_intracellular/fba/bacterial_colony/</code><br>
  <strong>Make target:</strong> <code>bacterial-colony</code><br>
  <strong>Executable:</strong> <code>bacterial-colony</code>
</div>

## Biological question

Can spatial structure alone create heterogeneous metabolism in an initially homogeneous, clonal colony? This model couples colony growth and mechanics to oxygen, glucose, acetate, and CO₂ fields.

## Run

```bash
make bacterial-colony
make
./bacterial-colony config/PhysiCell_settings.xml
```

!!! danger "Check timestep synchronization"

    Before treating a run as supported PhysiCelldFBA behavior, verify that `intracellular_dt` in the selected cell definition equals `dt_diffusion`. Some in-progress colony configurations in the repository use different values and should be reconciled before a definitive reproduction run.

## Important files

| File | Role |
| --- | --- |
| `config/Ecoli_core.xml` | Reduced metabolic model |
| `config/iML1515_trimmed.xml` | Trimmed genome-scale model |
| `config/PhysiCell_settings*.xml` | Model/configuration variants |
| `config/bacterial_colonies_n1.csv` | Initial cell placement |
| `scripts/plot_bacteria_colony.py` | Spatial plotting |
| `scripts/preprocess_dfba_model.ipynb` | Model preprocessing |
| `scripts/trimming_model.ipynb` | Reduction workflow |

## Reported result

During the first four hours, biomass grows approximately exponentially. As the colony expands, nutrient penetration becomes limiting and total growth decelerates.

At approximately nine hours:

- peripheral cells have access to oxygen and glucose and grow rapidly;
- intermediate cells show progressively reduced growth;
- a non-viable nutrient-limited core appears;
- peripheral cells secrete acetate through overflow metabolism;
- inner cells show reduced secretion or net acetate uptake.

This creates a spatial division of metabolic labor without defining multiple phenotypes in advance.

<div class="figure-placeholder">
  <strong>Figure placeholder — spatial metabolism in an E. coli colony</strong>
  Biomass trajectory, extracellular fields, growth rate, glucose uptake, and acetate exchange
  <code>docs/assets/figures/final/ecoli-colony-spatial-metabolism.png</code>
</div>

## Parameters that shape the pattern

- Domain and voxel size determine gradient resolution.
- Boundary glucose and oxygen determine resource supply.
- Initial colony geometry controls diffusion distance.
- Uptake kinetics determine the transition from replete to limited regimes.
- The ATP-maintenance/death rule controls emergence of the non-viable core.
- Metabolic-model trimming can change alternative pathways and acetate use.

Validate any trimmed network outside PhysiCell over both glucose-rich and acetate-consuming conditions before comparing spatial results.
