#include <iostream>
#include <vector>

#include "../src/dfba_Model.h"

using namespace std;

int main (int argc, const char *argv[])
{

    if (argc != 2)
    {
        cout << endl << "Usage: readSBML filename" << endl << endl;
        return 1;
    }
    const char *sbml_fileame = argv[1];

    dFBAModel *model = new dFBAModel();
    std::cout << "Reading SBML model from: " << sbml_fileame << " ";
    model->readSBMLModel(sbml_fileame);
    std::cout << "Ok!" << std::endl;
    std::cout << "Model " << model->getId() << "has been correctly loaded" << std::endl;
    std::cout << "Initializing LP model: ";
    model->initProblem();
    std::cout << "Ok!" << std::endl;
    model->setReactionLowerBound("R_EX_glc__D_e", -100.0);
    model->setReactionLowerBound("R_EX_o2_e", 0.0);


    std::cout << "Testing FBA: ";
    dFBASolution solution = model->optimize();
    if ( solution.status == "optimal"){
        std::cout << "OPTIMAL SOLUTION FOUND" << std::endl;
        solution = model->getSolution();
        float fopt = solution.getObjectiveValue();
	    std::cout << "Objective value: " << fopt << std::endl;   
    }

    delete model;
    return 0;
}
