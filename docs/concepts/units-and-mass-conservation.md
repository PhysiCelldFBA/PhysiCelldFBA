# Units and mass conservation

PhysiCelldFBA couples metabolic fluxes normalized by dry biomass to extracellular concentrations and geometrical cell volume. Correct conversions are therefore part of the model, not an implementation detail.

## Core quantities

| Quantity | Units used by the coupling |
| --- | --- |
| Extracellular concentration | mM = mmol/L |
| SBML metabolic/exchange flux | mmol/gDW/h |
| Cell and voxel volume | µm³ |
| Cellular dry mass | gDW |
| PhysiCell time | min |

Because 1 g/mL is numerically equal to 1 pg/µm³, a `cell_density` expressed in g/mL can be multiplied by a solid volume in µm³ to obtain dry mass in pg.

## Cell dry weight

PhysiCell stores total cell volume and a fluid fraction. PhysiCelldFBA estimates:

```text
V_solid = V_cell × (1 − fluid_fraction)
m_DW(pg) = V_solid(µm³) × cell_density(pg/µm³)
m_DW(g) = m_DW(pg) × 10⁻¹²
```

The assumption is a constant dry-mass density. It connects the geometric agent to fluxes reported per gram dry weight.

## Concentration-dependent uptake

For local concentration `[S]`, the kinetic uptake ceiling is:

```text
v_kin = Vmax × [S] / (Km + [S])
```

`Km` is a half-saturation concentration in mM. `Vmax` is the maximum uptake rate. A `Vmax` of zero prevents uptake but does not prevent the FBA solution from secreting through that exchange reaction.

## Local available-mass limit

A kinetic rate alone can remove more material than is present in a small or depleted voxel. The absolute amount available is:

```text
M_voxel(mmol) = [S](mmol/L) × V_voxel(µm³) × 10⁻¹⁵(L/µm³)
```

For cell dry weight `m_DW` and synchronized timestep `Δt` in minutes, the largest feasible uptake flux is:

```text
v_limit = 60 × M_voxel / (m_DW × Δt)
```

The SBML lower bound is the negative of the smaller positive limit:

```text
v_lower = −min(v_kin, v_limit)
```

This sign follows the common constraint-based convention that uptake is negative and secretion is positive.

## Exchange flux to BioFVM feedback

For an exchange flux `v_i`:

```text
E_i(mmol/min/cell) = v_i(mmol/gDW/h) × m_DW(g) / 60
```

The implementation converts this agent-level rate to BioFVM's volume-based representation before setting the net export rate. Positive `E_i` adds substrate; negative `E_i` removes it.

## Growth and volume

The biomass objective is interpreted as a specific growth rate `µ` in h⁻¹. Over a timestep `Δt` in minutes:

```text
V(t + Δt) = V(t) × exp[(µ / 60) × Δt]
```

This continuous biomass/volume update is distinct from cell division, which is a discrete agent event. The total mass can match an analytical yield even when the final number of agents does not represent fractional cells.

## Validation

The [mass-conservation test](../examples/unit-test.md) starts from a known glucose amount and compares final biomass with the yield predicted independently by the same metabolic model. Agreement below 1% validates the combined dry-weight scaling, exchange conversion, and local uptake limit for that controlled scenario.

!!! warning "Unit labels must agree with model semantics"

    Reaction bounds and objective fluxes are interpreted as mmol/gDW/h and h⁻¹ respectively. Changing an XML label does not perform a conversion. Confirm the actual units of every SBML model and numerical parameter before comparing models.
