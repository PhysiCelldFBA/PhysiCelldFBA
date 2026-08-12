# Glossary

**Agent-based model (ABM):** A model in which individual cells are represented as discrete agents with their own state and behavior.

**BioFVM:** PhysiCell's reaction-diffusion solver for extracellular substrates.

**Constraint-based metabolic model:** A stoichiometric network whose feasible reaction fluxes are defined by mass-balance equations and flux bounds.

**dFBA:** Dynamic flux balance analysis: repeated FBA solves as extracellular conditions and model constraints change over time.

**Exchange reaction:** An SBML reaction representing uptake from or secretion to the environment. PhysiCelldFBA maps selected exchanges to BioFVM substrates.

**FBA:** Flux balance analysis: optimization of a linear objective over steady-state stoichiometric constraints and reaction bounds.

**FBC:** The SBML Flux Balance Constraints package, used to encode reaction bounds and objectives.

**Genome-scale metabolic model (GEM/GSM):** A constraint-based reconstruction intended to cover most known metabolic reactions of an organism or cell type.

**gDW / pgDW:** Grams or picograms of cellular dry weight. Metabolic fluxes are normally normalized by gDW.

**Intracellular timestep:** The configured interval between intracellular model updates. In supported PhysiCelldFBA configurations it equals the diffusion timestep.

**Km:** The extracellular concentration at which the Michaelis-Menten uptake ceiling is half of `Vmax`.

**Local mass limit:** The maximum uptake compatible with the absolute amount of substrate in a cell's current voxel over one synchronized timestep.

**Objective reaction:** The reaction optimized by FBA, normally a biomass pseudo-reaction for growth.

**Off-lattice:** A cell representation in which agent positions are continuous rather than restricted to grid sites. BioFVM concentrations are still discretized on voxels.

**SBML:** Systems Biology Markup Language, the XML standard used here to encode metabolic models.

**Syntrophy:** A metabolic interaction in which products released by one organism support another, often making the community dependent on shared exchange.

**Vmax:** The maximum uptake ceiling in the concentration-dependent transport relation.

**Voxel:** A BioFVM volume element that stores local substrate concentrations.
