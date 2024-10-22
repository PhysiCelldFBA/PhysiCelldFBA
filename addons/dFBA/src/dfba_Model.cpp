/*
 dFBAintracellular.cpp
 *
 *  Created on: 11 jun. 2019
 *      Author: mponce
 */

#include "dfba_Model.h"


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

    if (this->handler != NULL)
        delete handler;
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

bool dFBAModel::hasMetabolite(std::string mId)
{
    std::map<std::string, int>::iterator itr;
    itr = this->metaboliteIndexer.find(mId);
    return itr != this->metaboliteIndexer.end();
}

void dFBAModel::addMetabolite(dFBAMetabolite* met)
{
    if (!this->hasMetabolite( met->getId() ))
    {
        this->metabolites.push_back(met);
        this->metaboliteIndexer[met->getId()] = this->metabolites.size() - 1;
    }

}

const dFBAMetabolite* dFBAModel::getMetabolite(std::string mId)
{
    if (this->hasMetabolite(mId))
    {
        int idx = this->metaboliteIndexer[mId];
        return this->metabolites[idx];
    }
    return nullptr;
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
    return rxn->getUpperBound();
}

void dFBAModel::setReactionUpperBound(std::string rId, double upperBound)
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
    {
        rxn->setUpperBound(upperBound);
        int colIdx = this->reactionsIndexer[rId];
        this->problem.setColumnUpper(colIdx, upperBound);
    }
}

double dFBAModel::getReactionLowerBound(std::string rId) 
{
    dFBAReaction* rxn = this->getReaction(rId);
    return rxn->getLowerBound();
}

void dFBAModel::setReactionLowerBound(std::string rId, double lowerBound)
{
    dFBAReaction* rxn = this->getReaction(rId);
    if (rxn)
    {
        rxn->setLowerBound(lowerBound);
        int colIdx = this->reactionsIndexer[rId];
        this->problem.setColumnLower(colIdx, lowerBound);
    }

}

void dFBAModel::addReaction(dFBAReaction* rxn)
{
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
    Model* model = document->getModel();

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
        this->addMetabolite(metabolite);
    }

    for(unsigned int i = 0; i < model->getNumReactions(); i++)
    {
        Reaction* sbml_reaction = listOfReactions->get(i);

        dFBAReaction* reaction = new dFBAReaction(sbml_reaction->getId());
        reaction->setName(sbml_reaction->getName());

        FbcReactionPlugin* rxnFbc = static_cast<FbcReactionPlugin*> (sbml_reaction->getPlugin("fbc"));
        if ( rxnFbc )
        {
            // Getting dFBAReaction's upper and lower bounds
            const std::string lbId = rxnFbc->getLowerFluxBound();
            double lb = listOfParameters->get(lbId)->getValue();
            reaction->setLowerBound(lb);

            const std::string ubId = rxnFbc->getUpperFluxBound();
            double ub = listOfParameters->get(ubId)->getValue();
            reaction->setUpperBound(ub);
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
                this->addMetabolite(metabolite);
            }
            
            const dFBAMetabolite* metabolite = this->getMetabolite(sbml_species->getSpecies());
            if (metabolite != nullptr)
                reaction->addMetabolite(metabolite, stoich_coef);
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
                this->addMetabolite(metabolite);
            }
            const dFBAMetabolite* metabolite = this->getMetabolite(sbml_species->getSpecies());
            if (metabolite != nullptr)
                reaction->addMetabolite(metabolite, stoich_coef);
            else
                std::cout << "ERROR: dFBAMetabolite " << sbml_species->getSpecies() << " not found" << std::endl;
        }
        this->addReaction(reaction);
    }

    // The following code is intended to extract the objective function from the sbml using
    // the FbcdFBAModelPlugin; then the coefficients are assigned to the corresponding dFBAReactions
    FbcModelPlugin* mplugin = static_cast<FbcModelPlugin*>(model->getPlugin("fbc"));
    ListOfObjectives* listOfObjectives =  mplugin->getListOfObjectives();
    Objective* objective =  mplugin->getObjective(listOfObjectives->getActiveObjective());
    ListOfFluxObjectives* listOfFluxObjectives = objective->getListOfFluxObjectives();

    for(unsigned int i=0; i <listOfFluxObjectives->getNumFluxObjectives(); i++)
    {
        FluxObjective* fluxObjective = listOfFluxObjectives->get(i);
        std::string rId = fluxObjective->getReaction();
        double objectiveCoefficient = fluxObjective->getCoefficient();
        dFBAReaction* dFBAReaction = this->getReaction(rId);
        dFBAReaction->setObjectiveCoefficient(objectiveCoefficient);
    }

    delete document;
}

void dFBAModel::initProblem()
{
    int n_rows = this->getNumMetabolites();
    int n_cols = this->getNumReactions();

    handler = new CoinMessageHandler(nullptr);
    std::cout << "Initilizing LP problem n=" << n_rows << std::endl;
    handler->setLogLevel(0);
    problem.passInMessageHandler(handler);

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
        col_lb[col_idx] = rxn->getLowerBound();
        col_ub[col_idx] = rxn->getUpperBound();
        objective[col_idx] = rxn->getObjectiveCoefficient();

        const std::map<const dFBAMetabolite*, double> metabolites = rxn->getMetabolites();
        CoinPackedVector col;
        for(auto it=metabolites.begin(); it!=metabolites.end(); it++)
        {
            const dFBAMetabolite* met = it->first;
            double stoich_coeff = it->second;
            int row_idx = this->metaboliteIndexer[met->getId()];
            col.insert(row_idx, stoich_coeff);
        }
        matrix.appendCol(col);
    }

    this->problem.loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
    this->problem.setOptimizationDirection(-1);

    delete[] col_lb;
    delete[] col_ub;
    delete[] row_lb;
    delete[] row_ub;
    delete[] objective;

    this->is_initialized = true;
}

void dFBAModel::initModel(const char* sbmlFileName)
{
    this->readSBMLModel(sbmlFileName);
    std::cout << "SBML model correctly loeaded: " << sbmlFileName << std::endl;
    this->initProblem();
}

void dFBAModel::writeProblem(const char *filename)
{
    this->problem.writeLp(filename);
}

dFBASolution dFBAModel::optimize()
{
    std::cout << "Running FBA... " << std::endl;
    std::cout << "Status before " << this->problem.statusOfProblem() << std::endl;
    this->problem.initialSolve();
    this->problem.primal();
    std::cout << "Status after running " << this->problem.statusOfProblem() << std::endl;
    std::cout << "Before checking... " << std::endl;
    if ( problem.isProvenOptimal() )
    {
        std::cout << "Optimal solution found... ";
        const double *columnPrimal = this->problem.getColSolution();
        std::map<std::string,double> fluxes;
        std::map<std::string,double> reduced_costs;

        std::string status;

        double fopt =  problem.getObjValue();
        
        if (problem.status() == 0){
            status = "optimal";
        }
        else if (problem.status() == 1){
            status = "infeasible";
        }
        else{
            status = "unknown";
        }
        
        for(dFBAReaction* reaction: this->reactions)
        {
            int column_idx = this->reactionsIndexer[reaction->getId()];
            double flux = columnPrimal[column_idx];
            fluxes[reaction->getId()] = flux;
            reaction->setFluxValue(flux);
        }

        solution.objective_value = fopt;
        solution.status = status;
        solution.fluxes = fluxes;
    }
    else
    {
        std::cout << "huston... ";
        for(dFBAReaction* reaction: this->reactions)
        { reaction->setFluxValue(0.0); }
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


