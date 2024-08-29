/*
 * dfba_Reaction.cpp
 *
 *  Created on: 13 jun. 2019
 *      Author: mponce
 */

#include "dfba_Reaction.h"

#include <iostream>
#include <vector>
#include <map>

#include "FBA_metabolite.h"



dFBAReaction::dFBAReaction(std::string id)
{
    this->id = id;
    this->name = "";
    this->lowerBound = 0;
    this->upperBound = 1000;
    this->objectiveCoefficient = 0;
    this->fluxValue = 0;
}

dFBAReaction::~dFBAReaction()
{

}

const std::string & dFBAReaction::getId() const
{
    return this->id;
}

void dFBAReaction::setName(std::string name)
{
    this->name = name;
}

const std::string & dFBAReaction::getName() const
{
    return this->name;
}

void dFBAReaction::setLowerBound(double lowerBound)
{
    this->lowerBound = lowerBound;
}

double dFBAReaction::getLowerBound()
{
    return this->lowerBound;
}

void dFBAReaction::setUpperBound(double upperBound)
{
    this->upperBound = upperBound;
}

double dFBAReaction::getUpperBound()
{
    return this->upperBound;
}

void dFBAReaction::setObjectiveCoefficient(double ojectiveCoefficient)
{
    this->objectiveCoefficient = ojectiveCoefficient;
}

double dFBAReaction::getObjectiveCoefficient()
{
    return this->objectiveCoefficient;
}

void dFBAReaction::setFluxValue(double fluxValue)
{
    this->fluxValue = fluxValue;
}

double dFBAReaction::getFluxValue()
{
    return this->fluxValue;
}

int dFBAReaction::getNumberOfMetabolites()
{
    return this->metabolites.size();
}

const std::map < const dFBAMetabolite *, double >&dFBAReaction::getMetabolites() const
{
    return this->metabolites;
}

bool dFBAReaction::reversible()
{
    return lowerBound < 0;
}

bool dFBAReaction::hasMetabolite(std::string mId)
{
    std::map < std::string, const dFBAMetabolite *>::iterator it;
    it = this->idMetaboliteMap.find(mId);
    bool hasMetabolite = (it != this->idMetaboliteMap.end());
    return hasMetabolite;
}


void dFBAReaction::addMetabolite(const dFBAMetabolite * met, double stoich)
{
    if (this->hasMetabolite(met->getId()))
    {
        this->metabolites[met] += stoich;
    }
    else
    {
        this->idMetaboliteMap[met->getId()] = met;
        this->metabolites[met] = stoich;
    }
}

std::vector < std::string > dFBAReaction::getReactants()
{
    std::vector < std::string > reactants;
    for (auto itr = this->metabolites.begin();
            itr != this->metabolites.end(); ++itr)
    {

        const dFBAMetabolite *met = itr->first;
        double sotich = itr->second;
        if (sotich < 0)
        {
            std::string mId = met->getId();
            reactants.push_back(mId);
        }
    }

    return reactants;
}

std::vector < std::string > dFBAReaction::getProducts()
{
    std::vector < std::string > products;
    for (auto itr = this->metabolites.begin();
            itr != this->metabolites.end(); ++itr)
    {

        const dFBAMetabolite *met = itr->first;
        double sotich = itr->second;
        if (sotich > 0)
        {
            std::string mId = met->getId();
            products.push_back(mId);
        }
    }
    return products;
}

double dFBAReaction::getStoichCoefficient(std::string mId)
{
    double coefficient = 0;
    if (hasMetabolite(mId))
    {
        const dFBAMetabolite *met = this->idMetaboliteMap[mId];
        coefficient = this->metabolites[met];
    }
    return coefficient;
}

std::string dFBAReaction::getReactionString()
{
    std::vector < std::string > compounds;
    std::string reactionString = "";

    compounds = this->getReactants();
    for (unsigned int i = 0; i < compounds.size(); i++)
    {
        std::string & mId = compounds[i];
        const dFBAMetabolite *met = this->idMetaboliteMap[mId];
        // Change the sign (-) of the coeeficient for printing
        float coeff = -1. * this->getStoichCoefficient(mId);
        if (int (coeff) == coeff)
            reactionString += std::to_string(int (coeff));
        else
            reactionString += std::to_string(coeff);

        //reactionString += " " + met->getName();
        reactionString += " " + met->getId();

        if (i < compounds.size() - 1)
            reactionString += " + ";
    }
    if (reversible())
        reactionString += " <==> ";
    else
        reactionString += " --> ";

    compounds = getProducts();
    for (unsigned int i = 0; i < compounds.size(); i++)
    {
        const std::string & mId = compounds[i];
        const dFBAMetabolite *met = this->idMetaboliteMap[mId];
        float coeff = getStoichCoefficient(mId);
        if (int (coeff) == coeff)
            reactionString += std::to_string(int (coeff));
        else
            reactionString += std::to_string(coeff);

        reactionString += " " + met->getId();

        if (i < compounds.size() - 1)
            reactionString += " + ";
    }
    return reactionString;
}
