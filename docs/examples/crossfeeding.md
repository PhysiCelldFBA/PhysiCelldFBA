# Microbial cross-feeding

<div class="key-path">
  <strong>Directory:</strong> <code>sample_projects_intracellular/fba/crossfeeding/</code><br>
  <strong>Make target:</strong> <code>crossfeeding</code><br>
  <strong>Executable:</strong> <code>crossfeeding</code>
</div>

## Biological question

Can two species carrying different genome-scale metabolic networks form a spatially structured syntrophic community through emergent exchange?

The model combines:

- *Clostridium beijerinckii* (CB; iCB925), which ferments supplied glucose and releases hydrogen and acetate;
- *Methanosarcina barkeri* (MB; iMG746), which can use hydrogen and CO₂ through hydrogenotrophic methanogenesis or acetate through acetoclastic methanogenesis, producing methane.

## Run

```bash
make crossfeeding
make
./crossfeeding config/PhysiCell_settings_crossfeeding.xml
```

The root sample target also makes the cross-feeding XML the default `config/PhysiCell_settings.xml`.

## Configuration pattern

The XML contains two cell definitions with separate `<intracellular type="dfba">` blocks. Each references its own SBML model and reaction identifiers, while both map to shared BioFVM fields such as hydrogen, acetate, CO₂, and methane.

This is the central multi-model pattern:

```text
CB SBML reactions ─┐
                   ├─ shared BioFVM substrates
MB SBML reactions ─┘
```

The extracellular name is the common key; the two SBML exchange IDs do not have to match.

## Reported result

Across the 50-hour scenario:

- CB grows rapidly under direct glucose supply;
- CB releases hydrogen and acetate;
- MB grows more slowly on CB-derived substrates;
- MB consumes hydrogen and CO₂ and increasingly uses acetate;
- methane production rises;
- diffusion and consumption produce spatial microgradients and niche separation.

The relative use of hydrogenotrophic and acetoclastic pathways is not prescribed. It follows from each cell's locally feasible optimum.

<div class="figure-placeholder">
  <strong>Figure placeholder — syntrophic microbial cross-feeding</strong>
  Community schematic, population growth, exchange rates, substrate totals, and spatial gradients
  <code>docs/assets/figures/final/microbial-crossfeeding.png</code>
</div>

## What this example teaches

- How to assign different SBML models to different cell definitions.
- How a secreted product of one model becomes an uptake substrate for another.
- How to interpret species-specific fluxes and shared population-level mass.
- How spatial separation and diffusion can limit an otherwise feasible syntrophy.

When replacing either reconstruction, validate its medium, exchange directions, maintenance demands, and growth objective independently before running the coupled community.
