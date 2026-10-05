/*
 * dfba_GPR.cpp
 *
 *  Created on: 2026
 *      Author: mponce
 */

#include "dfba_GPR.h"

#include <algorithm>

dFBAGPR::dFBAGPR()
{
    this->root = -1;
}

bool dFBAGPR::empty() const
{
    return this->root < 0 || this->nodes.empty();
}

void dFBAGPR::clear()
{
    this->nodes.clear();
    this->root = -1;
}

int dFBAGPR::addGeneNode(int gene_index)
{
    GPRNode node;
    node.type = GPRNode::GENE;
    node.gene_index = gene_index;
    this->nodes.push_back(node);
    return (int)this->nodes.size() - 1;
}

int dFBAGPR::addLogicalNode(GPRNode::NodeType type, const std::vector<int>& children)
{
    GPRNode node;
    node.type = type;
    node.gene_index = -1;
    node.children = children;
    this->nodes.push_back(node);
    return (int)this->nodes.size() - 1;
}

void dFBAGPR::setRoot(int node_index)
{
    this->root = node_index;
}

int dFBAGPR::getRoot() const
{
    return this->root;
}

int dFBAGPR::getNumNodes() const
{
    return (int)this->nodes.size();
}

bool dFBAGPR::eval(const std::vector<dFBAGene*>& genes) const
{
    // A reaction without a GPR is not under genetic control.
    if (this->empty())
        return true;

    return this->evalNode(this->root, genes);
}

bool dFBAGPR::evalNode(int node_index, const std::vector<dFBAGene*>& genes) const
{
    if (node_index < 0 || node_index >= (int)this->nodes.size())
        return true;

    const GPRNode& node = this->nodes[node_index];

    if (node.type == GPRNode::GENE)
    {
        // A reference to a gene absent from the model cannot constrain the
        // reaction; treat it as present rather than silently disabling flux.
        if (node.gene_index < 0 || node.gene_index >= (int)genes.size())
            return true;

        return genes[node.gene_index]->isActive();
    }

    if (node.children.empty())
        return true;

    if (node.type == GPRNode::AND)
    {
        for (unsigned int i = 0; i < node.children.size(); i++)
        {
            if (!this->evalNode(node.children[i], genes))
                return false;
        }
        return true;
    }

    // GPRNode::OR
    for (unsigned int i = 0; i < node.children.size(); i++)
    {
        if (this->evalNode(node.children[i], genes))
            return true;
    }
    return false;
}

std::vector<int> dFBAGPR::getGeneIndices() const
{
    std::vector<int> out;
    if (!this->empty())
        this->collectGenes(this->root, out);

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

void dFBAGPR::collectGenes(int node_index, std::vector<int>& out) const
{
    if (node_index < 0 || node_index >= (int)this->nodes.size())
        return;

    const GPRNode& node = this->nodes[node_index];

    if (node.type == GPRNode::GENE)
    {
        if (node.gene_index >= 0)
            out.push_back(node.gene_index);
        return;
    }

    for (unsigned int i = 0; i < node.children.size(); i++)
        this->collectGenes(node.children[i], out);
}

std::string dFBAGPR::toString(const std::vector<dFBAGene*>& genes) const
{
    if (this->empty())
        return "";

    return this->nodeToString(this->root, genes);
}

std::string dFBAGPR::nodeToString(int node_index, const std::vector<dFBAGene*>& genes) const
{
    if (node_index < 0 || node_index >= (int)this->nodes.size())
        return "";

    const GPRNode& node = this->nodes[node_index];

    if (node.type == GPRNode::GENE)
    {
        if (node.gene_index < 0 || node.gene_index >= (int)genes.size())
            return "?";
        return genes[node.gene_index]->getId();
    }

    const std::string op = (node.type == GPRNode::AND) ? " and " : " or ";

    std::string out = "";
    for (unsigned int i = 0; i < node.children.size(); i++)
    {
        const int child_index = node.children[i];
        std::string child = this->nodeToString(child_index, genes);

        // Parenthesise a child whose operator differs from this node's, so the
        // rendered string reparses to the same tree.
        if (child_index >= 0 && child_index < (int)this->nodes.size())
        {
            const GPRNode& child_node = this->nodes[child_index];
            if (child_node.type != GPRNode::GENE && child_node.type != node.type)
                child = "(" + child + ")";
        }

        out += child;
        if (i + 1 < node.children.size())
            out += op;
    }

    return out;
}
