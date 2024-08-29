/**
 * \brief Wrapper class
 *
 *
 * Created on 06/11/2019
 * M. Ponce-de-Leon, Barcelona Supercomputing Center
 */

#ifndef __FBA_model_h__
#define __FBA_model_h__

#include <iostream>
#include <map>
#include <vector>

#include <sbml/SBMLTypes.h>
#include <sbml/packages/fbc/common/FbcExtensionTypes.h>

#include <coin/CoinPackedMatrix.hpp>
#include <coin/CoinPackedVector.hpp>
#include <coin/ClpSimplex.hpp>

#include "dfba_Metabolite.h"
#include "dfba_Reaction.h"
#include "dfba_Solution.h"

LIBSBML_CPP_NAMESPACE_USE



class dFBAModel
{
private:
	/** \brief Constraint-Based Model Class to perform FBA*/

	std::string id;

		/** \brief vector of reaction objects*/
		std::vector<dFBAMetabolite*> metabolites;

	/** \brief map between metabolites' ids and metabolites' references **/
	std::map<std::string, int> metaboliteIndexer;

		/** \brief vector of reaction objects*/
		std::vector<dFBAReaction*> reactions;

	/** \brief map between reaction IDs and reaction references */
	std::map< std::string, int> reactionsIndexer;

		/** \brief solution  */
		dFBASolution solution;

	bool is_initialized = false;

public:

	/** \brief Constructor */
	FBA_model();

		/** \brief Constructor */
		dFBAModel();

		/** \brief Destructor */
		~dFBAModel();

	/** \brief a metabolite pointer using a string Id*/
	const FBA_metabolite* getMetabolite(std::string mId);

		/** \brief a metabolite pointer using a string Id*/
		const dFBAMetabolite* getMetabolite(std::string mId);

		/** \brief Add new metabolite to the model*/
		void addMetabolite(dFBAMetabolite* met);

		
		/** \brief Check if there is a reaction with a given ID*/
		bool hasReaction(std::string rId);
		
		/** \brief Get a reaction pointer using string ID*/
		dFBAReaction* getReaction(std::string rId);

		/** \brief Add new reaction to the model*/
		void addReaction(dFBAReaction* rxn);

	/** \brief Get the number of model reactions*/
	const int getNumReactions();

	/** \brief Get the number of model metabolites*/
	const int getNumMetabolites();

	/** \brief Get a metabolite pointer using s string Id*/
	const std::vector<FBA_metabolite*> getListOfMetabolites() const;

		/** \brief Get a metabolite pointer using s string Id*/
		const std::vector<dFBAMetabolite*> getListOfMetabolites() const;

		/** \brief Get the list of reaction pointers*/
		const std::vector<dFBAReaction*> getListOfReactions() const;

			/** \brief Get the list of reaction pointers*/
		std::vector<dFBAReaction*> getListOfBoundaryReactions();

	/** \brief Get the ClpSimplex model */
	const ClpSimplex* getLpModel() const;

	/** \brief Parse and read a metabolic model from a SBML file*/
	void readSBMLModel(const char* sbmlFileName);
	
	/** \brief Get the ClpSimplex model */
	void initLpModel();

	void initFBAmodel(const char* sbmlFileName);

	/** \brief Write LP problem in MPS format */
	void writeLp(const char *filename);

	/** \brief Run FBA using primal method */
	void runFBA();

		/** \brief Run FBA using primal method */
		dFBASolution optimize();

		/** \brief Get solution status */
		std::string getId() { return this->id; };

		/** \brief Get solution status */
		bool getSolutionStatus();

		/** \brief Get objective value */
		float getObjectiveValue();

		dFBASolution getSolution(){ return this->solution; }
};



#endif
