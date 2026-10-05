/*
 * dfba_Gene.cpp
 *
 *  Created on: 2026
 *      Author: mponce
 */

#include "dfba_Gene.h"

dFBAGene::dFBAGene(std::string id)
{
    this->id = id;
    this->name = "";
    this->label = "";
    this->active = true;
}

dFBAGene::~dFBAGene() {
}

dFBAGene::dFBAGene(const dFBAGene& copy)
{
    this->id = copy.id;
    this->name = copy.name;
    this->label = copy.label;
    this->active = copy.active;
}

const std::string& dFBAGene::getId() const
{
    return this->id;
}

void dFBAGene::setName(std::string value)
{
    this->name = value;
}

const std::string& dFBAGene::getName() const
{
    return this->name;
}

void dFBAGene::setLabel(std::string value)
{
    this->label = value;
}

const std::string& dFBAGene::getLabel() const
{
    return this->label;
}

bool dFBAGene::isActive() const
{
    return this->active;
}

void dFBAGene::setActive(bool value)
{
    this->active = value;
}
