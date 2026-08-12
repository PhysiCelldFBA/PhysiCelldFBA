# Prepare an SBML metabolic model

PhysiCelldFBA reads constraint-based metabolic models through libSBML and its Flux Balance Constraints package. The safest starting point is an SBML model that already solves correctly in an independent FBA tool such as COBRApy.

## Required model information

Before connecting a model to PhysiCell, confirm that it contains:

- a stoichiometric network with valid reactants, products, and coefficients;
- finite or intentionally unbounded lower and upper flux bounds through SBML-FBC;
- an objective reaction, normally a biomass pseudo-reaction;
- an exchange reaction for every extracellular metabolite that PhysiCell will couple;
- unique reaction identifiers that are preserved when the model is exported;
- flux units consistent with mmol/gDW/h.

PhysiCelldFBA looks reactions up by identifier. Display names and biochemical annotations are useful for humans but do not replace stable IDs.

## Objective reaction

The `objective_reaction` in the PhysiCell XML must identify a reaction in the SBML model. At initialization, PhysiCelldFBA sets that reaction's upper bound to `max_growth_rate` and uses the model objective for the standard solve.

For biomass-driven growth, the optimal objective value is treated as a specific growth rate in h⁻¹. A different objective can be used through custom code; the [motility example](../examples/metabolic-driven-motility.md) temporarily optimizes `R_ATPM`.

## Exchange reactions

An extracellular field is connected to an SBML exchange through:

```xml
<exchange substrate="glucose">
    <fba_flux>R_EX_glc__D_e</fba_flux>
    <Km units="mM">0.02</Km>
    <Vmax units="mmol/gDW/h">8.0</Vmax>
</exchange>
```

Check all three names:

1. `glucose` must exactly match a BioFVM variable name.
2. `R_EX_glc__D_e` must exactly match an SBML reaction ID.
3. The exchanged species and reaction direction must represent the intended metabolite.

PhysiCelldFBA uses the conventional sign rule: uptake is negative and secretion is positive. The local concentration constrains the reaction's lower bound. A zero `Vmax` disables uptake through the PhysiCell coupling while still allowing positive secretion if the metabolic model permits it.

## Bounds and units

The implementation interprets:

- exchange and internal fluxes as mmol/gDW/h;
- the biomass objective as h⁻¹;
- extracellular concentrations as mM.

SBML metadata is not used to convert arbitrary flux units automatically. Models expressed per cell, per protein mass, per litre, or in µmol units must be converted before use.

## Large models

Every metabolic agent keeps an independent model and solves it repeatedly. A genome-scale reconstruction can therefore dominate memory and runtime. Recommended preparation:

1. Remove blocked reactions and unused exchange reactions.
2. Preserve reactions required under every simulated environmental condition.
3. Validate the reduced model against the original over relevant nutrient regimes.
4. Compare objective and exchange fluxes before using the reduced model in PhysiCell.
5. Begin with one cell and a short simulation.

The colony and cancer-core example directories contain preprocessing and trimming notebooks that illustrate project-specific workflows. Treat them as examples rather than universal reduction algorithms.

## Validation checklist

- [ ] The model loads without fatal libSBML errors.
- [ ] Independent FBA returns an optimal solution.
- [ ] Every configured exchange reaction exists.
- [ ] The objective reaction exists and has the intended units.
- [ ] Uptake is represented by negative flux.
- [ ] Bounds are physiologically plausible.
- [ ] The same environmental constraints produce the expected yield outside PhysiCell.
- [ ] A one-cell PhysiCelldFBA test conserves mass.

You can generate an initial YAML mapping with [`generate_dfba_yaml.py`](helper-scripts.md), but its inferred substrate names, concentrations, bounds, and objective must still be reviewed.
