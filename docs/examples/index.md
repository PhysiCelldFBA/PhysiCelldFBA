# Examples and use cases

The supplied projects progress from a quantitative one-cell test to spatial single-species, multicellular, multispecies, and behavior-coupled models.

| Example | Main question | Metabolic model(s) | Coupling demonstrated |
| --- | --- | --- | --- |
| [Mass-conservation test](unit-test.md) | Does substrate consumption predict the correct dry biomass? | *E. coli* core | Units, growth, and local mass balance |
| [E. coli acetate switch](acetate-switch.md) | Can a metabolic lifestyle switch emerge from changing resources? | *E. coli* core / iML1515 variants | Temporal adaptation and by-product reuse |
| [E. coli colony](bacterial-colony.md) | Can one clonal population form spatial metabolic niches? | *E. coli* core / trimmed iML1515 | Gradients, heterogeneous growth, and acetate exchange |
| [Cancer core metabolism](cancer-core-metabolism.md) | How does vessel distance organize tumour-cell metabolism? | MCF7-specific reduced reconstruction | Many diffusing metabolites and metabolic death |
| [Microbial cross-feeding](crossfeeding.md) | Can distinct species form a syntrophic community? | iCB925 and iMG746 | Multiple models and interspecies exchange |
| [Metabolism-driven motility](metabolic-driven-motility.md) | Can metabolic capacity control a behavior other than growth? | *E. coli* core | Alternative objective, ATP allocation, and chemotaxis |

## Common project layout

```text
sample_projects_intracellular/fba/<example>/
├── Makefile
├── main.cpp
├── config/
│   ├── PhysiCell_settings*.xml
│   └── metabolic model(s)
├── custom_modules/
└── scripts/
```

Some projects contain several configurations, preprocessing notebooks, or figure-generation scripts. The default XML is a starting point; manuscript-scale runs can use a named alternative and may require substantially more resources.

## Common workflow

```bash
# Select a sample from the repository root.
make <sample-target>

# Compile the selected root project.
make

# Run with the default root configuration.
./<executable>

# Or pass a specific configuration.
./<executable> config/<settings-file>.xml
```

!!! warning "Sample targets replace the root project"

    Selecting a sample copies its files into the repository root. Preserve unrelated work or use a separate worktree first.

## Reproducibility checklist

Before comparing a run with a figure or reported value:

- record the Git commit and chosen XML file;
- keep `intracellular_dt = dt_diffusion`;
- record the SBML model and any preprocessing;
- retain the random seed and initial-cell CSV;
- record boundary conditions and simulation duration;
- preserve the copied XML in the output directory;
- use the analysis script associated with the same example version.

The result summaries on these pages describe the manuscript scenarios. They are scientific reference points, not bit-for-bit regression assertions across platforms or changed configurations.
