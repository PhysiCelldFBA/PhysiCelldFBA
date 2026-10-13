/*
 dFBAintracellular.cpp
 *
 *  Created on: 11 jun. 2019
 *      Author: mponce
 */

#include "dfba_Model.h"
#include "dfba_Reaction.h"  // Include the full definition

/* Default dFBAModel used to initialize the initial cell */
dFBAModel default_dFBAModel;


dFBAModel::dFBAModel()
{
    this->id = "none";
    this->is_initialized = false;
    this->handler = NULL;
    this->solution.status = "none";
}

dFBAModel::~dFBAModel() {
    for(dFBAReaction* rxn: this->reactions)
        delete rxn;

    for(dFBAMetabolite* met: this->metabolites)
        delete met;

    for(dFBAGene* gene: this->genes)
        delete gene;

    if (this->handler != nullptr)
        delete handler;
}

dFBAModel::dFBAModel(const dFBAModel& copy) {
    // Copy primitive members
    this->id = copy.id;
    this->is_initialized = copy.is_initialized;
    this->track_flux_map = copy.track_flux_map;

    // Deep copy metabolites
    for (auto original_metabolite : copy.metabolites) {
        this->metabolites.push_back(new dFBAMetabolite(*original_metabolite));
    }

    // Deep copy reactions
    for (auto original_reaction : copy.reactions) {
        this->reactions.push_back(new dFBAReaction(*original_reaction));
    }

    // Deep copy genes, so a knockout in one cell never reaches another
    for (auto original_gene : copy.genes) {
        this->genes.push_back(new dFBAGene(*original_gene));
    }

    // Copy maps
    this->metaboliteIndexer = copy.metaboliteIndexer;
    this->reactionsIndexer = copy.reactionsIndexer;
    this->geneIndexer = copy.geneIndexer;
    this->geneToReactions = copy.geneToReactions;

    // Copy solution (assuming it can be copied by value)
    this->solution = copy.solution;

    // Handle ClpSimplex problem copying (depends on how `ClpSimplex` needs to be cloned)
    if (copy.is_initialized) {
        this->initProblem();
    }

    // Copy message handler
    if (copy.handler != nullptr) {
        this->handler = new CoinMessageHandler(*copy.handler);
    } else {
        this->handler = nullptr;
    }
}



dFBAModel& dFBAModel::operator=(const dFBAModel& other) {
    if (this == &other) {
        return *this; // Handle self-assignment
    }

    // Clean up existing resources to prevent memory leaks
    for (auto met : metabolites) {
        delete met;
    }
    metabolites.clear();

    for (auto rxn : reactions) {
        delete rxn;
    }
    reactions.clear();

    for (auto gene : genes) {
        delete gene;
    }
    genes.clear();

    if (handler != nullptr) {
        delete handler;
    }

    // Copy primitive members
    this->id = other.id;
    this->is_initialized = other.is_initialized;
    this->track_flux_map = other.track_flux_map;

    // Deep copy metabolites
    for (auto original_metabolite : other.metabolites) {
        this->metabolites.push_back(new dFBAMetabolite(*original_metabolite));
    }

    // Deep copy reactions
    for (auto original_reaction : other.reactions) {
        this->reactions.push_back(new dFBAReaction(*original_reaction));
    }

    // Deep copy genes, so a knockout in one cell never reaches another
    for (auto original_gene : other.genes) {
        this->genes.push_back(new dFBAGene(*original_gene));
    }

    // Copy other maps
    this->metaboliteIndexer = other.metaboliteIndexer;
    this->reactionsIndexer = other.reactionsIndexer;
    this->geneIndexer = other.geneIndexer;
    this->geneToReactions = other.geneToReactions;

    // Copy solution (assuming it can be copied by value)
    this->solution = other.solution;

    // Handle ClpSimplex problem copying
    if (other.is_initialized) {
        this->initProblem(); 
        this->is_initialized = true;
    }

    // Copy message handler
    if (other.handler != nullptr) {
        this->handler = new CoinMessageHandler(*other.handler);
    } else {
        this->handler = nullptr;
    }

    return *this;
}



void dFBAModel::clear(){

    this->id = "none";
    this->is_initialized = false;
    this->track_flux_map = false;
    this->reactions.clear();
    this->metabolites.clear();
    this->genes.clear();
    this->reactionsIndexer.clear();
    this->metaboliteIndexer.clear();
    this->geneIndexer.clear();
    this->geneToReactions.clear();
    this->solution.clear();
    this->problem = ClpSimplex();
    this->handler = NULL;

    return;
}

const ClpSimplex* dFBAModel::getLpModel() const
{
    return &this->problem;
}

const int dFBAModel::getNumReactions()
{
    return this->reactions.size();
}

const int dFBAModel::getNumMetabolites()
{
    return this->metabolites.size();
}

bool dFBAModel::hasMetabolite(const std::string& id) const
{
    return this->metaboliteIndexer.find(id) != this->metaboliteIndexer.end();
}

dFBAMetabolite* dFBAModel::addMetabolite(const std::string& id)
{
    // Check if metabolite already exists by its ID
    auto it = this->metaboliteIndexer.find(id);
    if (it == this->metaboliteIndexer.end()) {
        // If metabolite does not exist, create it and add it to the model
        dFBAMetabolite* newMet = new dFBAMetabolite(id);
        this->metabolites.push_back(newMet);
        this->metaboliteIndexer[id] = this->metabolites.size() - 1;
        return newMet;
    } else {
        // If metabolite already exists, return the existing one
        return this->metabolites[it->second];
    }
}

dFBAMetabolite* dFBAModel::getMetabolite(const std::string& mId) const {
    auto it = this->metaboliteIndexer.find(mId);
    if (it != this->metaboliteIndexer.end()) {
        int idx = it->second;
        return this->metabolites[idx];
    }
    return nullptr;  // Return nullptr if metabolite doesn't exist
}



const std::vector<dFBAMetabolite*> dFBAModel::getListOfMetabolites() const
{
    return this->metabolites;
}

bool dFBAModel::hasReaction(std::string rId)
{
    std::map<std::string, int>::iterator it;
    it = this->reactionsIndexer.find(rId);

    return it != this->reactionsIndexer.end();
}

dFBAReaction* dFBAModel::getReaction(std::string rId)
{
    if (this->hasReaction(rId))
    {
        int idx = this->reactionsIndexer[rId];
        dFBAReaction* rxn = this->reactions[idx];
        return rxn;
    }
    return nullptr;
}


double dFBAModel::getReactionUpperBound(std::string rId) 
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
        return rxn->getUpperBound();
    else {
        std::cerr << "Error: Reaction " << rId << " not found in getReactionUpperBound" << std::endl;
        return 0.0;
    }
}

void dFBAModel::setReactionUpperBound(std::string rId, double upperBound)
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
    {
        // A knocked-out reaction must stay clamped at zero: remember the new
        // bound for when the responsible gene is restored, but leave the LP shut.
        if (rxn->isKnockedOut())
        {
            rxn->setKnockoutSavedUpperBound(upperBound);
            return;
        }
        rxn->setUpperBound(upperBound);
        if (this->is_initialized)
        {
            int colIdx = this->reactionsIndexer[rId];
            this->problem.setColumnUpper(colIdx, upperBound);
        }
    }
    else{
        std::cerr << "Reaction with ID " << rId << " not found in the model." << std::endl;
    }
}

double dFBAModel::getReactionLowerBound(std::string rId) 
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
        return rxn->getLowerBound();
    else {
        std::cerr << "Error: Reaction " << rId << " not found in getReactionLowerBound" << std::endl;
        return 0.0;
    }
}

void dFBAModel::setReactionLowerBound(std::string rId, double lowerBound)
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
    {
        // A knocked-out reaction must stay clamped at zero: remember the new
        // bound for when the responsible gene is restored, but leave the LP shut.
        if (rxn->isKnockedOut())
        {
            rxn->setKnockoutSavedLowerBound(lowerBound);
            return;
        }
        rxn->setLowerBound(lowerBound);
        if (this->is_initialized)
        {
            int colIdx = this->reactionsIndexer[rId];
            this->problem.setColumnLower(colIdx, lowerBound);
        }
    }
    else{
        std::cerr << "Reaction with ID " << rId << " not found in the model." << std::endl;
    }

}

void dFBAModel::addReaction(dFBAReaction* rxn)
{
    if (this->is_initialized) {
        std::cerr << "CRITICAL ERROR: Cannot add reactions after model initialization." << std::endl;
        exit(-1);
    }
    if (!this->hasReaction( rxn->getId() ))
    {
        this->reactions.push_back(rxn);
        this->reactionsIndexer[rxn->getId()] = this->reactions.size() - 1;
    }
}

const int dFBAModel::getReactionIndex(std::string rId)
{
    if (this->hasReaction(rId))
        return this->reactionsIndexer[rId];
    else
        return -1;
}

const std::vector<dFBAReaction*> dFBAModel::getListOfReactions() const
{
    return this->reactions;
}

std::vector<dFBAReaction*> dFBAModel::getListOfBoundaryReactions()
{
    std::vector<dFBAReaction*> listOfBoundarys;
    for(dFBAReaction* dFBAReaction: this->reactions)
    {
        if (dFBAReaction->getNumberOfMetabolites() == 1)
        {
            listOfBoundarys.push_back(dFBAReaction);        
        }
    }
    return listOfBoundarys;
}

std::vector<std::string> dFBAModel::getListOfBoundaryReactionIds()
{
    std::vector<std::string> listOfBoundaryIds;
    for(dFBAReaction* dFBAReaction: this->reactions)
    {
        if (dFBAReaction->getNumberOfMetabolites() == 1)
        {
            listOfBoundaryIds.push_back(dFBAReaction->getId());
        }
    }
    return listOfBoundaryIds;
}

void dFBAModel::readSBMLModel(const char* sbmlFileName)
{
    SBMLReader reader;
    SBMLDocument* document = reader.readSBML(sbmlFileName);

    // Check if the document was successfully read
    if (document == nullptr) {
    std::cerr << "Error: SBMLDocument is null for " << sbmlFileName << "\n";
    return;
    }

    // Optional: run consistency checks (can add more warnings)
    // document->checkConsistency();

    const unsigned nFatal = document->getNumErrors(LIBSBML_SEV_FATAL);
    const unsigned nError = document->getNumErrors(LIBSBML_SEV_ERROR);
    const unsigned nWarn  = document->getNumErrors(LIBSBML_SEV_WARNING);

    if (nFatal + nError > 0) {
        std::cerr << "Failed to read SBML: " << sbmlFileName << "\n";
        document->printErrors();  // show why it failed
        delete document;
        return;
    }

    // Non-fatal: show warnings but continue
    if (nWarn) {
        std::cerr << "Validation warning(s): " << nWarn << "\n";
        document->printErrors();
    }

    Model* model = document->getModel();

    // Check if the model was successfully retrieved
    if (model == nullptr) {
        std::cerr << "Error: Model is null in SBML file: " << sbmlFileName << std::endl;
        delete document;
        return;
    }

    ListOfSpecies* listOfSpecies = model->getListOfSpecies();
    ListOfReactions* listOfReactions = model->getListOfReactions();
    ListOfParameters* listOfParameters = model->getListOfParameters();

    this->id = model->getId();

    for (unsigned int i = 0; i < model->getNumSpecies(); i++)
    {
        Species* species = listOfSpecies->get(i);
        // Skipping boundary dFBAMetabolites
        if ( species->getBoundaryCondition() )
            continue;

        dFBAMetabolite* metabolite = new dFBAMetabolite(species->getId());
        metabolite->setName(species->getName());
        this->addMetabolite(metabolite->getId());
    }

    for(unsigned int i = 0; i < model->getNumReactions(); i++)
    {
        Reaction* sbml_reaction = listOfReactions->get(i);

        dFBAReaction* reaction = new dFBAReaction(sbml_reaction->getId());
        // std::cout << "Adding reaction with ID: " << reaction->getId() << std::endl;
        reaction->setName(sbml_reaction->getName());

        FbcReactionPlugin* rxnFbc = static_cast<FbcReactionPlugin*> (sbml_reaction->getPlugin("fbc"));
        if ( rxnFbc )
        {
            // Getting dFBAReaction's upper and lower bounds
            std::string lbId = rxnFbc->getLowerFluxBound();
            Parameter* lbParam = listOfParameters->get(lbId);
            if (lbParam) {
                double lb = lbParam->getValue();
                reaction->setLowerBound(lb);
            }
            else {
                std::cout << "WARNING: Lower bound parameter " << lbId << " not found for reaction " << sbml_reaction->getId() << std::endl;
            }

            std::string ubId = rxnFbc->getUpperFluxBound();
            Parameter* ubParam = listOfParameters->get(ubId);
            if (ubParam) {
                double ub = ubParam->getValue();
                reaction->setUpperBound(ub);
            }
            else {
                std::cout << "WARNING: Upper bound parameter " << ubId << " not found for reaction " << sbml_reaction->getId() << std::endl;
            }
        }
        else {
            std::cout << "WARNING: No FBC plugin found for reaction " << sbml_reaction->getId() << std::endl;
        }
        int numReactans = sbml_reaction->getNumReactants();
        for(int j = 0; j < numReactans; j++)
        {
            SpeciesReference* sbml_species = sbml_reaction->getReactant(j);
            double stoich_coef = -1. * sbml_species->getStoichiometry();

            if ( !this->hasMetabolite(sbml_species->getSpecies()) )
            {
                dFBAMetabolite* metabolite = new dFBAMetabolite(sbml_species->getSpecies());
                metabolite->setName(sbml_species->getName());
                this->addMetabolite(metabolite->getId());
            }
            
            const dFBAMetabolite* metabolite = this->getMetabolite(sbml_species->getSpecies());
            if (metabolite != nullptr)
                reaction->addMetabolite(metabolite->getId(), stoich_coef);
            else
                std::cout << "ERROR: dFBAMetabolite " << sbml_species->getSpecies() << " not found" << std::endl;
        }

        int numProducts = sbml_reaction->getNumProducts();
        for(int j = 0; j < numProducts; j++)
        {
            SpeciesReference* sbml_species = sbml_reaction->getProduct(j);
            double stoich_coef = sbml_species->getStoichiometry();

            if ( !this->hasMetabolite(sbml_species->getSpecies()) )
            {
                dFBAMetabolite* metabolite = new dFBAMetabolite(sbml_species->getSpecies());
                metabolite->setName(sbml_species->getName());
                this->addMetabolite(metabolite->getId());
            }
            const dFBAMetabolite* metabolite = this->getMetabolite(sbml_species->getSpecies());
            if (metabolite != nullptr)
                reaction->addMetabolite(metabolite->getId(), stoich_coef);
            else
                std::cout << "ERROR: dFBAMetabolite " << sbml_species->getSpecies() << " not found" << std::endl;
        }
        this->addReaction(reaction);
    }

    // Parse the gene products and the gene-protein-reaction rules. Genes come
    // from fbc:listOfGeneProducts; each reaction's rule comes from its
    // fbc:geneProductAssociation. Reactions without an association keep an empty
    // GPR and are therefore never disabled by a knockout.
    FbcModelPlugin* gene_plugin = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));
    if (gene_plugin)
    {
        for (unsigned int i = 0; i < gene_plugin->getNumGeneProducts(); i++)
        {
            GeneProduct* gene_product = gene_plugin->getGeneProduct(i);
            if (gene_product == nullptr)
                continue;

            int gene_index = this->addGene(gene_product->getId());
            if (gene_index >= 0)
            {
                this->genes[gene_index]->setName(gene_product->getName());
                this->genes[gene_index]->setLabel(gene_product->getLabel());
            }
        }

        for (unsigned int i = 0; i < model->getNumReactions(); i++)
        {
            Reaction* sbml_reaction = listOfReactions->get(i);
            dFBAReaction* reaction = this->getReaction(sbml_reaction->getId());
            if (reaction == nullptr)
                continue;

            FbcReactionPlugin* rxnFbc =
                static_cast<FbcReactionPlugin*>(sbml_reaction->getPlugin("fbc"));
            if (rxnFbc == nullptr || !rxnFbc->isSetGeneProductAssociation())
                continue;

            GeneProductAssociation* gpa = rxnFbc->getGeneProductAssociation();
            if (gpa == nullptr || gpa->getAssociation() == nullptr)
                continue;

            dFBAGPR gpr;
            int root = this->parseGPRAssociation(gpa->getAssociation(), gpr);
            if (root >= 0)
            {
                gpr.setRoot(root);
                reaction->setGPR(gpr);
            }
        }
    }

    this->buildGeneReactionIndex();

    // The following code is intended to extract the objective function from the sbml using
    // the FbcdFBAModelPlugin; then the coefficients are assigned to the corresponding dFBAReactions
    FbcModelPlugin* mplugin = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));
    if (mplugin) {
        ListOfObjectives* listOfObjectives = mplugin->getListOfObjectives();
        Objective* objective = (listOfObjectives != nullptr)
            ? mplugin->getObjective(listOfObjectives->getActiveObjective())
            : nullptr;
        if (objective == nullptr) {
            std::cerr << "WARNING: no active fbc objective in " << sbmlFileName
                      << "; all objective coefficients left at zero." << std::endl;
            delete document;
            return;
        }
        ListOfFluxObjectives* listOfFluxObjectives = objective->getListOfFluxObjectives();

        for (unsigned int i = 0; i < listOfFluxObjectives->getNumFluxObjectives(); i++)
        {
            FluxObjective* fluxObjective = listOfFluxObjectives->get(i);
            std::string rId = fluxObjective->getReaction();
            double objectiveCoefficient = fluxObjective->getCoefficient();
            dFBAReaction* dFBAReaction = this->getReaction(rId);
            if (dFBAReaction) {
                dFBAReaction->setObjectiveCoefficient(objectiveCoefficient);
            } else {
                std::cerr << "ERROR: dFBAReaction " << rId << " not found" << std::endl;
            }
        }
    }

    delete document;
}

void dFBAModel::initProblem()
{
    int n_rows = this->getNumMetabolites();
    int n_cols = this->getNumReactions();

    this->handler = new CoinMessageHandler(nullptr);
    // std::cout << "Initilizing LP problem n=" << n_rows << std::endl;
    this->handler->setLogLevel(0);
    this->problem.passInMessageHandler(this->handler);

    CoinPackedMatrix matrix;
    matrix.setDimensions(n_rows, 0);

    double* row_lb = new double[n_rows]; //the row lower bounds
    double* row_ub = new double[n_rows]; //the row upper bounds
    double* col_lb = new double[n_cols]; //the column lower bounds
    double* col_ub = new double[n_cols]; //the column upper bounds
    double* objective = new double[n_cols]; //the objective coefficients

    for(int i=0; i< n_rows; i++)
    {
        row_lb[i] = 0;
        row_ub[i] = 0;
    }
    for(dFBAReaction* rxn: this->reactions)
    {
        int col_idx = this->reactionsIndexer[rxn->getId()];
        //std::cout << "Adding reaction " << rxn->getId() << " with index: " << this->reactionsIndexer[rxn->getId()] << " to the LP problem" << std::endl;
        col_lb[col_idx] = rxn->getLowerBound();
        col_ub[col_idx] = rxn->getUpperBound();
        objective[col_idx] = rxn->getObjectiveCoefficient();

        const std::map<std::string, double>& metabolites = rxn->getMetabolites();

        CoinPackedVector col;
        for (auto it = metabolites.begin(); it != metabolites.end(); ++it) {
            const std::string& metId = it->first;  // Get metabolite ID
            double stoich_coeff = it->second;

            // Use the metabolite ID to get the index from metaboliteIndexer
            auto idx_it = this->metaboliteIndexer.find(metId);
            if (idx_it != this->metaboliteIndexer.end()) {
                int row_idx = idx_it->second;
                col.insert(row_idx, stoich_coeff);
            }
        }
        matrix.appendCol(col);
    }

    this->problem.loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
    this->problem.setOptimizationDirection(-1);

    this->problem.setPerturbation(50); // 50 is a standard default for perturbation. It means that the perturbation is always applied. 100 is the default value (automatic)

    delete[] col_lb;
    delete[] col_ub;
    delete[] row_lb;
    delete[] row_ub;
    delete[] objective;

    this->problem.setPrimalTolerance(1e-7);
    this->problem.setDualTolerance(1e-7);
    this->problem.scaling(3);
    this->problem.setMaximumIterations(10000);
    this->problem.setLogLevel(0);
    this->problem.initialSolve();

    this->is_initialized = true;
}

void dFBAModel::initModel(const char* sbmlFileName)
{
    this->readSBMLModel(sbmlFileName);
    std::cout << "SBML model correctly loaded: " << sbmlFileName << std::endl;
    this->initProblem();
}

void dFBAModel::writeProblem(const char *filename)
{
    this->problem.writeLp(filename);
}

dFBASolution dFBAModel::optimize()
{

    this->problem.dual();
/*     // If Dual fails, try Primal (more robust from scratch) --> if numerical issue occurs, try primal to make sure the cell is really dead
    if (this->problem.status() != 0) { 
        this->problem.primal();
    } */
    bool isOptimal = problem.isProvenOptimal();

    solution.fluxes.clear();

    std::string status;
    switch (this->problem.status()) {
        case 0:
            status = "optimal";
            break;
        case 1:
        case 4:
            status = "infeasible";
            break;
        default:
            status = "unknown";
            std::cerr << "Solver returned unknown status code: " << this->problem.status() << std::endl;
            break;
    }

    if (isOptimal)
    {
        const double *columnPrimal = this->problem.getColSolution();

        double fopt =  problem.getObjValue();

        // Column i of the LP is reactions[i]: addReaction() appends to both the
        // vector and the indexer together, so the index equals the position and
        // the string lookup this used to do was redundant.
        //
        // The string-keyed flux map is only populated on request. Building it
        // here meant clearing and re-inserting one node per reaction on every
        // solve -- for a cell solving every intracellular_dt, that was the
        // majority of the cost of optimize() and nothing in the addon reads it.
        const int n_rxn = (int)this->reactions.size();
        for(int i = 0; i < n_rxn; i++)
            this->reactions[i]->setFluxValue(columnPrimal[i]);

        if (this->track_flux_map)
        {
            for(int i = 0; i < n_rxn; i++)
                solution.fluxes[this->reactions[i]->getId()] = columnPrimal[i];
        }

        solution.objective_value = fopt;
        solution.status = status;

        // Debugging info
        //std::cout << "Optimal solution found: Objective = " << fopt << "\n";
    }
    else if (status == "infeasible")
    {
        solution.status = status;

        //std::cerr << "FBA optimization infeasible. Cell should die.\n";
        for (auto &reaction : this->reactions)
        {
            reaction->setFluxValue(0.0);
            if (this->track_flux_map)
                solution.fluxes[reaction->getId()] = 0.0;
        }
        this->solution.objective_value = 0.0;
    }
    else
    {
        solution.status = status;

        std::cerr << "ERROR: FBA optimization failed with unknown status!\n";
        std::cerr << "Reaction bounds at failure:\n";
        for (auto &reaction : this->reactions)
        {
            int idx = this->reactionsIndexer[reaction->getId()];
            double lb = reaction->getLowerBound();
            double ub = reaction->getUpperBound();
            std::cerr << reaction->getId() << ": [" << lb << ", " << ub << "]\n";
            reaction->setFluxValue(0.0);
        }
    }

    return solution;
}


bool dFBAModel::getSolutionStatus()
{
    if (this->is_initialized)
        return this->problem.isProvenOptimal();
    else
        return false;
}

double dFBAModel::getObjectiveValue()
{
    assert(this->is_initialized);
    if (this->problem.isProvenOptimal())
        return this->problem.getObjValue();
    else
        std::cout << "WARNING: Primal infeasible" << std::endl;
    return 0;
}

bool dFBAModel::isInitialized(){
    return this->is_initialized;
}

void dFBAModel::setObjectiveCoefficient(std::string rId, double coefficient)
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
    {
        rxn->setObjectiveCoefficient(coefficient);
        if (this->is_initialized)
        {
            int colIdx = this->reactionsIndexer[rId];
            this->problem.setObjectiveCoefficient(colIdx, coefficient);
        }
    }
    else
    {
        std::cerr << "Reaction with ID " << rId << " not found in setObjectiveCoefficient." << std::endl;
    }
}

void dFBAModel::clearObjective()
{
    for(dFBAReaction* rxn: this->reactions)
    {
        rxn->setObjectiveCoefficient(0.0);
        if (this->is_initialized)
        {
            int colIdx = this->reactionsIndexer[rxn->getId()];
            this->problem.setObjectiveCoefficient(colIdx, 0.0);
        }
    }
}

void dFBAModel::deepCopy(const dFBAModel& source) {
    // Clean up existing resources
    for (auto met : metabolites) delete met;
    metabolites.clear();
    for (auto rxn : reactions) delete rxn;
    reactions.clear();
    for (auto gene : genes) delete gene;
    genes.clear();
    if (handler != nullptr) { delete handler; handler = nullptr; }

    // Copy primitive members
    id = source.id;
    is_initialized = source.is_initialized;
    track_flux_map = source.track_flux_map;

    // Deep copy metabolites
    for (auto original_metabolite : source.metabolites)
        metabolites.push_back(new dFBAMetabolite(*original_metabolite));

    // Deep copy reactions
    for (auto original_reaction : source.reactions)
        reactions.push_back(new dFBAReaction(*original_reaction));

    // Deep copy genes, so a knockout in one cell never reaches another
    for (auto original_gene : source.genes)
        genes.push_back(new dFBAGene(*original_gene));

    // Copy indexers and solution
    metaboliteIndexer = source.metaboliteIndexer;
    reactionsIndexer = source.reactionsIndexer;
    geneIndexer = source.geneIndexer;
    geneToReactions = source.geneToReactions;
    solution = source.solution;

    // Copy saved objective state
    original_objective_reaction = source.original_objective_reaction;
    original_objective_coefficients = source.original_objective_coefficients;

    // Reinitialize LP problem if source was initialized
    if (source.is_initialized)
        initProblem();

    // Copy message handler
    if (source.handler != nullptr)
        handler = new CoinMessageHandler(*source.handler);
}

void dFBAModel::saveObjectiveState() {
    original_objective_coefficients.clear();
    for (const auto* reaction : reactions) {
        if (reaction->getObjectiveCoefficient() != 0.0) {
            original_objective_coefficients[reaction->getId()] = reaction->getObjectiveCoefficient();
        }
    }
}

void dFBAModel::restoreObjectiveState() {
    clearObjective();
    for (const auto& pair : original_objective_coefficients) {
        setObjectiveCoefficient(pair.first, pair.second);
    }
}

// =====================================================================
//  Genes, GPR rules and gene knockouts
// =====================================================================

bool dFBAModel::hasGene(const std::string& gId) const
{
    return this->geneIndexer.find(gId) != this->geneIndexer.end();
}

int dFBAModel::getGeneIndex(const std::string& gId) const
{
    std::map<std::string, int>::const_iterator it = this->geneIndexer.find(gId);
    if (it == this->geneIndexer.end())
        return -1;
    return it->second;
}

dFBAGene* dFBAModel::getGene(const std::string& gId) const
{
    int index = this->getGeneIndex(gId);
    if (index < 0)
        return nullptr;
    return this->genes[index];
}

int dFBAModel::addGene(const std::string& gId)
{
    int existing = this->getGeneIndex(gId);
    if (existing >= 0)
        return existing;

    dFBAGene* gene = new dFBAGene(gId);
    this->genes.push_back(gene);
    int index = (int)this->genes.size() - 1;
    this->geneIndexer[gId] = index;
    return index;
}

const int dFBAModel::getNumGenes()
{
    return (int)this->genes.size();
}

const std::vector<dFBAGene*> dFBAModel::getListOfGenes() const
{
    return this->genes;
}

std::string dFBAModel::getGeneReactionRule(const std::string& rId)
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn == nullptr)
        return "";
    return rxn->getGPR().toString(this->genes);
}

const std::vector<int>& dFBAModel::getReactionsOfGene(int gene_index) const
{
    static const std::vector<int> empty;
    if (gene_index < 0 || gene_index >= (int)this->geneToReactions.size())
        return empty;
    return this->geneToReactions[gene_index];
}

int dFBAModel::parseGPRAssociation(const FbcAssociation* association, dFBAGPR& gpr)
{
    if (association == nullptr)
        return -1;

    const int type_code = association->getTypeCode();

    if (type_code == SBML_FBC_GENEPRODUCTREF)
    {
        const GeneProductRef* ref = static_cast<const GeneProductRef*>(association);
        const std::string gene_id = ref->getGeneProduct();

        // A rule may reference a gene product that is missing from
        // listOfGeneProducts. Register it so the rule stays faithful to the file
        // rather than silently dropping a term.
        int gene_index = this->getGeneIndex(gene_id);
        if (gene_index < 0)
            gene_index = this->addGene(gene_id);

        return gpr.addGeneNode(gene_index);
    }

    if (type_code == SBML_FBC_AND || type_code == SBML_FBC_OR)
    {
        unsigned int num_children = 0;
        if (type_code == SBML_FBC_AND)
            num_children = static_cast<const FbcAnd*>(association)->getNumAssociations();
        else
            num_children = static_cast<const FbcOr*>(association)->getNumAssociations();

        std::vector<int> children;
        for (unsigned int i = 0; i < num_children; i++)
        {
            const FbcAssociation* child = nullptr;
            if (type_code == SBML_FBC_AND)
                child = static_cast<const FbcAnd*>(association)->getAssociation(i);
            else
                child = static_cast<const FbcOr*>(association)->getAssociation(i);

            int child_index = this->parseGPRAssociation(child, gpr);
            if (child_index >= 0)
                children.push_back(child_index);
        }

        if (children.empty())
            return -1;

        // An AND/OR with a single child is just that child.
        if (children.size() == 1)
            return children[0];

        const GPRNode::NodeType node_type =
            (type_code == SBML_FBC_AND) ? GPRNode::AND : GPRNode::OR;
        return gpr.addLogicalNode(node_type, children);
    }

    return -1;
}

void dFBAModel::buildGeneReactionIndex()
{
    this->geneToReactions.clear();
    this->geneToReactions.resize(this->genes.size());

    for (unsigned int r = 0; r < this->reactions.size(); r++)
    {
        dFBAReaction* rxn = this->reactions[r];
        if (!rxn->hasGPR())
            continue;

        std::vector<int> gene_indices = rxn->getGPR().getGeneIndices();
        for (unsigned int i = 0; i < gene_indices.size(); i++)
        {
            const int g = gene_indices[i];
            if (g >= 0 && g < (int)this->geneToReactions.size())
                this->geneToReactions[g].push_back((int)r);
        }
    }
}

bool dFBAModel::applyGeneStateToReaction(int reaction_index)
{
    if (reaction_index < 0 || reaction_index >= (int)this->reactions.size())
        return false;

    dFBAReaction* rxn = this->reactions[reaction_index];

    // Without a GPR the reaction is not under genetic control.
    if (!rxn->hasGPR())
        return false;

    const bool available = rxn->getGPR().eval(this->genes);

    // The bounds are applied directly here rather than through
    // setReactionLowerBound/setReactionUpperBound, because those deliberately
    // refuse to move the bounds of a reaction that is currently knocked out.
    const int colIdx = this->reactionsIndexer[rxn->getId()];

    if (!available && !rxn->isKnockedOut())
    {
        // Save the bounds the reaction currently holds, then clamp it shut.
        rxn->saveBoundsForKnockout();
        rxn->setKnockedOut(true);
        rxn->setLowerBound(0.0);
        rxn->setUpperBound(0.0);
        if (this->is_initialized)
        {
            this->problem.setColumnLower(colIdx, 0.0);
            this->problem.setColumnUpper(colIdx, 0.0);
        }
        return true;
    }

    if (available && rxn->isKnockedOut())
    {
        const double lb = rxn->getKnockoutSavedLowerBound();
        const double ub = rxn->getKnockoutSavedUpperBound();
        rxn->setKnockedOut(false);
        rxn->setLowerBound(lb);
        rxn->setUpperBound(ub);
        if (this->is_initialized)
        {
            this->problem.setColumnLower(colIdx, lb);
            this->problem.setColumnUpper(colIdx, ub);
        }
        return true;
    }

    return false;
}

int dFBAModel::setGeneActive(const std::string& gId, bool active)
{
    const int gene_index = this->getGeneIndex(gId);
    if (gene_index < 0)
    {
        std::cerr << "Gene with ID " << gId << " not found in the model." << std::endl;
        return 0;
    }

    dFBAGene* gene = this->genes[gene_index];
    if (gene->isActive() == active)
        return 0;

    gene->setActive(active);

    // Only the reactions that reference this gene can have changed.
    const std::vector<int>& affected = this->getReactionsOfGene(gene_index);
    int changed = 0;
    for (unsigned int i = 0; i < affected.size(); i++)
    {
        if (this->applyGeneStateToReaction(affected[i]))
            changed++;
    }

    return changed;
}

int dFBAModel::knockOutGene(const std::string& gId)
{
    return this->setGeneActive(gId, false);
}

int dFBAModel::restoreGene(const std::string& gId)
{
    return this->setGeneActive(gId, true);
}

bool dFBAModel::isGeneActive(const std::string& gId) const
{
    dFBAGene* gene = this->getGene(gId);
    if (gene == nullptr)
        return false;
    return gene->isActive();
}

std::vector<std::string> dFBAModel::getKnockedOutGenes() const
{
    std::vector<std::string> out;
    for (unsigned int i = 0; i < this->genes.size(); i++)
    {
        if (!this->genes[i]->isActive())
            out.push_back(this->genes[i]->getId());
    }
    return out;
}

std::vector<std::string> dFBAModel::getKnockedOutReactions() const
{
    std::vector<std::string> out;
    for (unsigned int i = 0; i < this->reactions.size(); i++)
    {
        if (this->reactions[i]->isKnockedOut())
            out.push_back(this->reactions[i]->getId());
    }
    return out;
}

bool dFBAModel::getTrackFluxMap() const
{
    return this->track_flux_map;
}

void dFBAModel::setTrackFluxMap(bool value)
{
    this->track_flux_map = value;
}

std::map<std::string, double> dFBAModel::getFluxDistribution() const
{
    // Built on demand from the reactions, which always hold the last solution.
    // Prefer this to dFBASolution::getFluxes(), which is only filled when flux
    // tracking is enabled.
    std::map<std::string, double> out;
    for (dFBAReaction* rxn : this->reactions)
        out[rxn->getId()] = rxn->getFluxValue();
    return out;
}
