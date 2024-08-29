/*
 * Metabolite.cpp
 *
 *  Created on: 13 jun. 2019
 *      Author: mponce
 */


#include "dfba_Metabolite.h"

#include <iostream>
#include <vector>
#include <map>



dFBAMetabolite::dFBAMetabolite(std::string id) {
   this->id = id;
}

dFBAMetabolite::~dFBAMetabolite() {

}

const std::string& dFBAMetabolite::getId() const {
  return this->id;
}

void dFBAMetabolite::setName(std::string value) {
  this->name = value;
}

const std::string& dFBAMetabolite::getName() const {
  return this->name;
}


