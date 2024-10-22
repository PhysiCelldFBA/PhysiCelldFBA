/*
 * dfba_Reaction.cpp
 *
 *  Created on: jun. 2022
 *      Author: mponce
 */

#include "dfba_Solution.h"

#include <iostream>
#include <vector>
#include <map>


dFBASolution::dFBASolution()
{
    this->objective_value = 0;
    this->status = "none";
}

dFBASolution::dFBASolution(double objective_value, std::string status, std::map<std::string,double> fluxes)
{
    this->objective_value = objective_value;
    this->status = status;
    this->fluxes = fluxes;
}

dFBASolution::~dFBASolution()
{

}
