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
#include <set>

#include "dfba_Metabolite.h"
#include "dfba_GPR.h"

// Forward declaration of dFBAModel
class dFBAModel;

class dFBAReaction
{
	private:
		std::string id;
		std::string name;

		double lowerBound;
		double upperBound;
		double objectiveCoefficient;
		double fluxValue;

		/** \brief gene-protein-reaction rule; empty when the reaction has none */
		dFBAGPR gpr;

		/** \brief true while the reaction is disabled because its GPR evaluates false */
		bool knockedOut;

		/** \brief bounds held while knocked out, restored when the GPR becomes true again */
		double knockoutSavedLowerBound;
		double knockoutSavedUpperBound;

		// Map of metabolite IDs to stoichiometric coefficients
    	std::map<std::string, double> metabolites;

	public:
		/** \brief Constructor */
		dFBAReaction(std::string id);

		/** \brief Destructor */
		~dFBAReaction();

		/** \brief Copy */
		dFBAReaction(const dFBAReaction& copy);

		const std::string& getId() const;

		void setName(std::string name);
		const std::string& getName() const;

		int getNumberOfMetabolites() const;

		void setLowerBound(double lowerBound);
		double getLowerBound();

		void setUpperBound(double upperBound);
		double getUpperBound();

		void setObjectiveCoefficient(double ojectiveCoefficient);
		double getObjectiveCoefficient() const;

		void setFluxValue(double flux_value);
		double getFluxValue();

		bool reversible() const;
		bool hasMetabolite(const std::string& mId) const;
		void addMetabolite(const std::string& mId, double stoich);

		std::vector<std::string> getReactants() const;
    	std::vector<std::string> getProducts() const;
		double getStoichCoefficient(const std::string& mId) const;

		std::string getReactionString(const dFBAModel& model) const;

		/** \brief Gene-protein-reaction rule */
		bool hasGPR() const;
		const dFBAGPR& getGPR() const;
		dFBAGPR& getGPR();
		void setGPR(const dFBAGPR& value);

		/**
		 * \brief Knockout bookkeeping.
		 *
		 * A knocked-out reaction has its LP bounds forced to zero. The bounds it
		 * held at the moment of the knockout are saved here so that restoring the
		 * responsible gene reinstates them exactly.
		 */
		bool isKnockedOut() const;
		void setKnockedOut(bool value);
		void saveBoundsForKnockout();
		double getKnockoutSavedLowerBound() const;
		double getKnockoutSavedUpperBound() const;
		void setKnockoutSavedLowerBound(double value);
		void setKnockoutSavedUpperBound(double value);

		const std::map<std::string, double>& getMetabolites() const;
};



#endif /* SRC_REACTION_H_ */
