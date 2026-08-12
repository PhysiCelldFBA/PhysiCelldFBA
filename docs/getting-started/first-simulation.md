# First simulation

The mass-conservation unit test is the shortest route from a fresh clone to a verified dFBA run. It follows one *E. coli* cell in a closed, glucose-limited domain and compares the final dry biomass with an independent FBA yield calculation.

## 1. Select and compile the project

From the repository root:

```bash
make dfba_unit_test
make
```

The first compilation downloads libSBML and Coin-OR CLP automatically if they are not already present. The resulting executable is `dfba-unit-test`.

!!! warning

    The sample target replaces the root project files. Use a clean checkout or worktree if the repository root currently contains another simulation project.

## 2. Run the validation configuration

```bash
./dfba-unit-test config/unit_test_mass_conversion_lit_params.xml
```

The configuration places one cell in a closed three-dimensional domain, sets glucose to 10 mM, and synchronizes:

```xml
<dt_diffusion units="min">0.01</dt_diffusion>
...
<intracellular_dt units="min">0.01</intracellular_dt>
```

The output directory is controlled by `<save><folder>` in the XML and is `output/` in this configuration.

## 3. Check the run

A successful run:

- reports that the SBML model was loaded;
- completes without an unknown solver status;
- writes initial, periodic, and final MultiCellDS snapshots;
- consumes the finite glucose pool without producing negative concentrations;
- stops biomass accumulation when the limiting substrate is exhausted.

The supplied analysis notebook is:

```text
sample_projects_intracellular/fba/dfba_unit_test/scripts/dfba_analysis.ipynb
```

After selecting the sample, it is also copied into the root `scripts/` directory.

## Expected result

The independent calculation uses the biomass yield of the same FBA model under the same constraints. For the manuscript configuration:

- without non-growth-associated ATP maintenance, the analytical biomass is approximately **4.08 pg dry weight** and the simulated biomass is approximately **4.09 pg**;
- with ATP maintenance, the corresponding values are approximately **3.87 pg** and **3.88 pg**.

Both comparisons have a relative error below 1%. Small differences can result from saved-output timing and numerical tolerances; a large mismatch points to a unit, timestep, substrate, or exchange-bound problem.

<figure markdown="span">
  ![Mass-conservation validation showing glucose depletion, carbon dioxide production, and simulated total biomass converging to the analytical prediction](../assets/figures/final/mass_conservation.png)
  <figcaption>Glucose consumption and CO₂ production accompany biomass accumulation; the simulated total biomass converges to the analytical prediction.</figcaption>
</figure>

## Next steps

- Read [Units and mass conservation](../concepts/units-and-mass-conservation.md) for the conversion behind this test.
- Read [PhysiCell XML](../configuration/physicell-xml.md) to understand the configuration block.
- Continue to the [E. coli acetate-switch example](../examples/acetate-switch.md) for an adaptive metabolic phenotype.
