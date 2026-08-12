# The `addons/dFBA` code

The add-on is a small C++ layer between PhysiCell's intracellular-model interface, libSBML, and Coin-OR CLP. Reading it from the cell-facing class inward is the easiest way to understand the implementation.

## Directory map

```text
addons/dFBA/
├── CMakeLists.txt
├── src/
│   ├── dfba_intracellular.{h,cpp}
│   ├── dfba_Model.{h,cpp}
│   ├── dfba_Reaction.{h,cpp}
│   ├── dfba_Metabolite.{h,cpp}
│   └── dfba_Solution.{h,cpp}
├── scripts/
│   └── reusable YAML examples
├── test/
│   └── low-level SBML and solver programs
└── ext/
    └── external headers and libraries
```

The runtime source currently includes Coin-OR CLP headers and links CLP/CoinUtils. libSBML parses SBML-FBC models. Other solver source trees that may be present under `ext/` are not selected merely by existing there.

## Core classes

### `dFBAIntracellular`

Files: `dfba_intracellular.h` and `dfba_intracellular.cpp`

This is the PhysiCell-facing class. It derives from `PhysiCell::Intracellular` and owns the coupling state for one cell:

- SBML file path and `dFBAModel`;
- dFBA timestep and next scheduled run;
- extracellular-substrate-to-reaction mappings;
- objective, density, reference volume, and growth cap;
- current growth and solution status;
- optional metabolic-death settings.

It parses the XML, constrains exchanges, solves the model, updates cell volume and BioFVM rates, and exposes helpers for retrieving fluxes or optimizing a temporary custom objective.

The copy constructor deep-copies the metabolic model, allowing daughter/new agents to maintain independent solver states instead of sharing one mutable linear program.

### `dFBAModel`

Files: `dfba_Model.h` and `dfba_Model.cpp`

`dFBAModel` owns the constraint-based network and CLP problem. It:

- loads SBML through libSBML;
- extracts metabolites, reactions, stoichiometry, bounds, and objectives;
- builds the sparse matrix used by CLP;
- changes reaction bounds as the microenvironment changes;
- solves the linear program;
- records fluxes and solver status;
- saves/restores the objective around custom optimizations.

Solver outcomes are normalized to `optimal`, `infeasible`, or `unknown`.

### `dFBAReaction`

Files: `dfba_Reaction.h` and `dfba_Reaction.cpp`

A reaction stores its ID, name, reversibility, bounds, objective coefficient, current flux, and metabolite stoichiometry. `dFBAModel` uses it both to construct the matrix and to expose flux/bound information after a solve.

### `dFBAMetabolite`

Files: `dfba_Metabolite.h` and `dfba_Metabolite.cpp`

This lightweight object stores a metabolite identifier and name. Model-level index maps connect it to reaction stoichiometry.

### `dFBASolution`

Files: `dfba_Solution.h` and `dfba_Solution.cpp`

The solution container holds the objective value, status string, and reaction-ID-to-flux map returned by optimization.

## Initialization lifecycle

1. PhysiCell encounters `<intracellular type="dfba">` in a cell definition.
2. `dFBAIntracellular` parses `settings`, `transport_model`, `growth_model`, and optional `death_model`.
3. Every BioFVM substrate name is resolved to a density index.
4. libSBML loads the model and reports fatal/error diagnostics.
5. Every configured exchange and objective reaction is validated.
6. `max_growth_rate` is applied to the objective-reaction upper bound.
7. Any death-trigger bound is configured.
8. The original objective state is saved for later custom-objective use.

Missing required blocks or identifiers stop initialization rather than silently producing an uncoupled cell.

## Per-update lifecycle

PhysiCell's cell container asks each intracellular model whether it needs an update. For dFBA cells that are due:

1. a pre-update hook can modify the cell;
2. `update_dfba_inputs` computes exchange bounds;
3. either the standard solve or `custom_optimization` runs;
4. `update_dfba_outputs` applies death, growth, and exchange feedback;
5. a post-update hook can consume the new metabolic state;
6. `next_dfba_run` advances.

The motility example uses the custom optimization hook to switch between biomass and ATP objectives without replacing the generic coupling machinery.

## Build integration

Sample Makefiles:

- define `ADDON_PHYSIDFBA`;
- add `addons/dFBA/src`, libSBML, and Coin-OR include paths;
- compile the five dFBA implementation objects;
- link the static CLP/CoinUtils and libSBML libraries;
- run `beta/setup_fba.py` when the Coin-OR sentinel header is absent.

The exact root Makefile changes when a sample is selected. If you add dFBA to a new project, copy the dFBA include, object, dependency, and link sections from a current supported sample rather than an older project.

## Extension points

- Add a new extracellular coupling by defining another `exchange`; no C++ change is needed.
- Add a new cell type/model by defining another dFBA cell definition.
- Use `pre_update_intracellular` or `post_update_intracellular` for project-specific logic around a solve.
- Use `custom_optimization` and `optimize_for_objective` when phenotype logic requires an objective other than biomass.

When adding behavior, preserve the equal dFBA/diffusion timestep rule and the input → solve → output order that protects mass consistency.
