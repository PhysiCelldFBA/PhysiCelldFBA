/*
 * dfba_Gene.h
 *
 * Gene product of a constraint-based model, parsed from the SBML fbc
 * listOfGeneProducts. A gene carries an "active" flag, true by default, which
 * is set to false by a knockout. Reaction availability is derived from these
 * flags through each reaction's gene-protein-reaction (GPR) rule.
 */

#ifndef SRC_GENE_H_
#define SRC_GENE_H_

#include <iostream>
#include <string>

class dFBAGene
{
	private:
		std::string id;
		std::string name;
		/** \brief fbc:label, usually the locus tag (e.g. b1854) */
		std::string label;
		/** \brief false once the gene has been knocked out */
		bool active;

	public:
		dFBAGene(std::string id);
		~dFBAGene();

		/** \brief Copy */
		dFBAGene(const dFBAGene& copy);

		const std::string& getId() const;

		void setName(std::string value);
		const std::string& getName() const;

		void setLabel(std::string value);
		const std::string& getLabel() const;

		/** \brief True unless the gene has been knocked out */
		bool isActive() const;
		void setActive(bool value);
};

#endif /* SRC_GENE_H_ */
