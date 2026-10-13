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
#include "dfba_Solution.h"
#include "dfba_Gene.h"
#include "dfba_GPR.h"

LIBSBML_CPP_NAMESPACE_USE

// Forward declaration of dFBAReaction
class dFBAReaction;

class dFBAModel
{
	private:
		/** \brief Model ID */
		std::string id;

		/** \brief vector of metabolite objects*/
		std::vector<dFBAMetabolite*> metabolites;

		/** \brief Saved objective reaction ID for restore */
		std::string original_objective_reaction;

		/** \brief Saved objective coefficients for restore */
		std::map<std::string, double> original_objective_coefficients;

		/** \brief map between metabolites' ids and metabolites' references **/
		std::map<std::string, int> metaboliteIndexer;

		/** \brief vector of gene objects, parsed from fbc:listOfGeneProducts */
		std::vector<dFBAGene*> genes;

		/** \brief map between genes' ids and their index in genes */
		std::map<std::string, int> geneIndexer;

		/**
		 * \brief for each gene index, the indices of the reactions whose GPR
		 * references it.
		 *
		 * Built once after parsing and immutable thereafter. A knockout only has
		 * to re-evaluate this subset instead of sweeping every reaction, which
		 * matters for genome-scale models: a typical iML1515 gene appears in a
		 * handful of the 2712 reactions.
		 */
		std::vector< std::vector<int> > geneToReactions;

		/** \brief vector of reaction objects*/
		std::vector<dFBAReaction*> reactions;

		/** \brief map between reaction IDs and reaction references */
		std::map< std::string, int> reactionsIndexer;

		/** \brief solution  */
		dFBASolution solution;

		/** \brief Coin CLP simplex model to encode the FBA problem**/
		ClpSimplex problem;
		
		CoinMessageHandler* handler;

		bool is_initialized = false;

		/**
		 * \brief whether optimize() also fills dFBASolution::fluxes.
		 *
		 * Off by default. Populating that string-keyed map costs one map
		 * insertion per reaction on every solve, which dominated optimize()
		 * for cells solving at every intracellular step. Nothing in the addon
		 * reads it; use getFluxDistribution() instead, which builds the same
		 * map on demand from the reactions.
		 */
		bool track_flux_map = false;

	public:

		/** \brief Constructor */
		dFBAModel();

		/** \brief Destructor */
		~dFBAModel();

		/** \brief Copy */
		dFBAModel(const dFBAModel& copy);

		/** \brief new assignment operator */
		dFBAModel& operator=(const dFBAModel& other);

		/** \brief Clear function */
		void clear();

		/** \brief Check if the model is initialized*/
		bool isInitialized();

		/** \brief Check if there is a metaboltie with a given ID*/
		bool hasMetabolite(const std::string& id) const;

		/** \brief a metabolite pointer using a string Id*/
		dFBAMetabolite* getMetabolite(const std::string& mId) const;

		/** \brief Add new metabolite to the model*/
		dFBAMetabolite* addMetabolite(const std::string& id);

		
		/** \brief Check if there is a reaction with a given ID*/
		bool hasReaction(std::string rId);
		
		/** \brief Get a reaction pointer using string ID*/
		dFBAReaction* getReaction(std::string rId);

		/** \brief Add new reaction to the model*/
		void addReaction(dFBAReaction* rxn);

		// ================ genes and GPR rules ================

		/** \brief Check if there is a gene with a given ID*/
		bool hasGene(const std::string& gId) const;

		/** \brief Get a gene pointer using a string ID, nullptr when absent*/
		dFBAGene* getGene(const std::string& gId) const;

		/** \brief Add a new gene to the model, returning its index*/
		int addGene(const std::string& gId);

		/** \brief Get the integer index of a gene, -1 when absent*/
		int getGeneIndex(const std::string& gId) const;

		/** \brief Get the number of model genes*/
		const int getNumGenes();

		/** \brief Get the list of gene pointers*/
		const std::vector<dFBAGene*> getListOfGenes() const;

		/** \brief Render a reaction's GPR rule as an infix string*/
		std::string getGeneReactionRule(const std::string& rId);

		/** \brief Indices of the reactions whose GPR references a given gene*/
		const std::vector<int>& getReactionsOfGene(int gene_index) const;

		/**
		 * \brief Knock out a gene and disable the reactions that this makes
		 * unavailable.
		 *
		 * Sets the gene inactive, then re-evaluates the GPR of every reaction
		 * that references it. A reaction whose rule now evaluates false has its
		 * current bounds saved and is clamped to [0, 0]. Returns the number of
		 * reactions disabled by this call. A no-op if the gene is unknown or
		 * already inactive.
		 */
		int knockOutGene(const std::string& gId);

		/**
		 * \brief Reactivate a gene and re-enable the reactions that this makes
		 * available again.
		 *
		 * Re-evaluates rather than blindly restoring: a reaction requiring
		 * "G_a and G_b" stays disabled when only one of the two is restored.
		 * Returns the number of reactions re-enabled.
		 */
		int restoreGene(const std::string& gId);

		/** \brief Set a gene active or inactive, applying the consequences*/
		int setGeneActive(const std::string& gId, bool active);

		/** \brief True when the gene exists and has not been knocked out*/
		bool isGeneActive(const std::string& gId) const;

		/** \brief Ids of every gene currently knocked out*/
		std::vector<std::string> getKnockedOutGenes() const;

		/** \brief Ids of every reaction currently disabled by a knockout*/
		std::vector<std::string> getKnockedOutReactions() const;

		/** \brief Get the integer index of a reaction*/
		const int getReactionIndex(std::string rId);
		
		/** \brief Get the upper bound of a reactions*/
		double getReactionUpperBound(std::string rId);
		
		/** \brief Set the upper bound of a reactions*/
		void setReactionUpperBound(std::string rId, double upperBound);
		
		/** \brief Get the upper bound of a reactions*/
		double getReactionLowerBound(std::string rId);
		
		/** \brief Set the lower bound of a reactions*/
		void setReactionLowerBound(std::string rId, double lowerBound);

		/** \brief Get the number of model reactions*/
		const int getNumReactions();

		/** \brief Get the number of model metabolites*/
		const int getNumMetabolites();

		/** \brief Get a metabolite pointer using s string Id*/
		const std::vector<dFBAMetabolite*> getListOfMetabolites() const;

		/** \brief Get the list of reaction pointers*/
		const std::vector<dFBAReaction*> getListOfReactions() const;

			/** \brief Get the list of reaction pointers*/
		std::vector<dFBAReaction*> getListOfBoundaryReactions();

		/** \brief Get the list of IDs of the boundary reactions*/
		std::vector<std::string> getListOfBoundaryReactionIds();

		/** \brief Get the ClpSimplex model */
		const ClpSimplex* getLpModel() const;

		/** \brief Parse and read a metabolic model from a SBML file*/
		void readSBMLModel(const char* sbmlFileName);

		/** \brief Build geneToReactions from the parsed GPR rules*/
		void buildGeneReactionIndex();
		
		/** \brief Initialize the ClpSimplex model */
		void initProblem();

		/** \brief Initialize the CBM model and the LP problem */
		void initModel(const char* sbmlFileName);

		/** \brief Write problem in LP format */
		void writeProblem(const char *filename);

		/** \brief Run FBA using primal method */
		dFBASolution optimize();

		/** \brief Get solution status */
		std::string getId() { return this->id; };

		/** \brief Get solution status */
		bool getSolutionStatus();

		/** \brief Get objective value */
		double getObjectiveValue();

		dFBASolution getSolution(){ return this->solution; }

		/** \brief Flux of every reaction from the last solve, built on demand */
		std::map<std::string, double> getFluxDistribution() const;

		/** \brief Whether optimize() also fills dFBASolution::fluxes (default false) */
		bool getTrackFluxMap() const;
		void setTrackFluxMap(bool value);

		void setObjectiveCoefficient(std::string rId, double coefficient);
		void clearObjective();

		// Deep copy from another model (cleans up existing data first)
		void deepCopy(const dFBAModel& source);

		// Save/restore objective function state
		void saveObjectiveState();
		void restoreObjectiveState();

	private:
		/**
		 * \brief Re-evaluate one reaction's GPR and apply the result to its bounds.
		 *
		 * Returns true when the reaction's enabled/disabled state changed.
		 */
		bool applyGeneStateToReaction(int reaction_index);

		/** \brief Recursively translate an fbc association into GPR nodes*/
		int parseGPRAssociation(const FbcAssociation* association, dFBAGPR& gpr);
};

#endif
