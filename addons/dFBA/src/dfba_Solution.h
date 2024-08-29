/*
 * dfba_Reaction.h
 *
 *  Created on: jun. 2022
 *      Author: mponce
 */

#ifndef SRC_SOLUTION_H_
#define SRC_SOLUTION_H_

#include <iostream>
#include <map>
#include <vector>


class dFBASolution
{
	private:
	    /** \brief Constraint-Based Model Class that contains a solution*/


    public:

        float objective_value;
        std::string status;
        std::map<std::string, float> fluxes;
        std::map<std::string, float> reduced_costs;

        dFBASolution();
        dFBASolution(float objective_value, std::string status, std::map<std::string,float> fluxes);
        ~dFBASolution();

        float getObjectiveValue(){ return this->objective_value; }
        std::string getStatus(){ return this->status; };
        std::map<std::string, float> getFluxes(){ return fluxes; };

};

#endif 