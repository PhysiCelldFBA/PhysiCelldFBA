/*
 * dfba_GPR.h
 *
 * Gene-protein-reaction (GPR) rule attached to a reaction.
 *
 * The rule is stored as an expression tree over AND / OR / gene-reference
 * nodes, mirroring the fbc:geneProductAssociation element of SBML and the GPR
 * class of cobrapy. Real genome-scale models nest these: iML1515 contains
 * rules of the form (g1 and g2) or (g3 and g4) at depth 4, so a general tree
 * is required rather than a flat two-level scheme.
 *
 * The tree is held in a flat std::vector of nodes that reference each other by
 * index, and genes are referenced by their index in the owning dFBAModel
 * rather than by pointer. This is deliberate: every PhysiCell cell holds a deep
 * copy of its dFBAModel, and a flat vector of indices is copied correctly by
 * the implicitly generated copy constructor, with no ownership to manage and no
 * way for two cells to end up sharing GPR state.
 */

#ifndef SRC_GPR_H_
#define SRC_GPR_H_

#include <string>
#include <vector>

#include "dfba_Gene.h"

struct GPRNode
{
	enum NodeType { GENE, AND, OR };

	NodeType type;
	/** \brief index of the gene in the model, valid when type == GENE */
	int gene_index;
	/** \brief indices of child nodes, used when type is AND or OR */
	std::vector<int> children;

	GPRNode() : type(GENE), gene_index(-1) {}
};

class dFBAGPR
{
	private:
		std::vector<GPRNode> nodes;
		/** \brief index of the root node, -1 when the rule is empty */
		int root;

	public:
		dFBAGPR();

		/** \brief True when the reaction carries no GPR rule */
		bool empty() const;

		void clear();

		/** \brief Append a gene-reference node, returning its index */
		int addGeneNode(int gene_index);

		/** \brief Append an AND/OR node over already-added children */
		int addLogicalNode(GPRNode::NodeType type, const std::vector<int>& children);

		void setRoot(int node_index);
		int getRoot() const;

		/**
		 * \brief Evaluate the rule against the current gene activity states.
		 *
		 * AND is true when every child is true, OR when any child is true. An
		 * empty rule evaluates to true: a reaction without a GPR is not under
		 * genetic control and can never be disabled by a knockout.
		 */
		bool eval(const std::vector<dFBAGene*>& genes) const;

		/** \brief Indices of every gene referenced by the rule, without duplicates */
		std::vector<int> getGeneIndices() const;

		/** \brief Number of nodes in the expression tree */
		int getNumNodes() const;

		/**
		 * \brief Render the rule as an infix string, e.g. "G_b0351 or (G_a and G_b)".
		 *
		 * Parentheses are emitted around a nested sub-expression whose operator
		 * differs from its parent, matching the convention used by cobrapy.
		 */
		std::string toString(const std::vector<dFBAGene*>& genes) const;

	private:
		bool evalNode(int node_index, const std::vector<dFBAGene*>& genes) const;
		void collectGenes(int node_index, std::vector<int>& out) const;
		std::string nodeToString(int node_index, const std::vector<dFBAGene*>& genes) const;
};

#endif /* SRC_GPR_H_ */
