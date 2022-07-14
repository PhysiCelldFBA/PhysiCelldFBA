/*
 * dfba_Reaction.h
 *
 *  Created on: jun. 2022
 *      Author: mponce
 */

#ifndef SRC_REACTION_H_
#define SRC_REACTION_H_

#include <iostream>
#include <vector>
#include <map>

#include "dfba_Metabolite.h"

class dFBAReaction
{
	private:
		std::string id;
		std::string name;

		double lowerBound;
		double upperBound;
		double objectiveCoefficient;
		double fluxValue;

		std::map<const dFBAMetabolite*, double> metabolites;
		std::map<std::string, const dFBAMetabolite*> idMetaboliteMap;

	public:
		dFBAReaction(std::string id);
		~dFBAReaction();

		const std::string& getId() const;

		void setName(std::string name);
		const std::string& getName() const;

		int getNumberOfMetabolites();

		void setLowerBound(double lowerBound);
		double getLowerBound();

		void setUpperBound(double upperBound);
		double getUpperBound();

		void setObjectiveCoefficient(double ojectiveCoefficient);
		double getObjectiveCoefficient();

		void setFluxValue(double flux_value);
		double getFluxValue();

		bool reversible();
		bool hasMetabolite(std::string mId);
		void addMetabolite(const dFBAMetabolite* met, double stoich);

		std::vector<std::string> getReactants();
		std::vector<std::string> getProducts();
		double getStoichCoefficient(std::string mId);

		std::string getReactionString();

		const std::map<const dFBAMetabolite*, double>& getMetabolites() const;
};



#endif /* SRC_REACTION_H_ */
