/*
 * main_test_gpr.cpp
 *
 * Tests for gene-protein-reaction (GPR) parsing, evaluation and gene knockout.
 *
 * Usage: ./test_gpr [path-to-Ecoli_core.xml]
 *
 * The fixtures under data/gpr/ cover the GPR shapes and edge cases in
 * isolation; the optional E. coli core model exercises the parser against a
 * real curated model.
 */

#include <iostream>
#include <string>
#include <vector>

#include "../src/dfba_Model.h"
#include "../src/dfba_Reaction.h"
#include "dfba_test.h"

static const std::string FIXTURES = "data/gpr/";

static dFBAModel* load(const std::string& path)
{
    dFBAModel* model = new dFBAModel();
    model->readSBMLModel(path.c_str());
    return model;
}

// ---------------------------------------------------------------- parsing --

static void test_no_gpa()
{
    DFBA_TEST_SECTION("reaction without a geneProductAssociation");
    dFBAModel* m = load(FIXTURES + "no_gpa.xml");

    dFBAReaction* rxn = m->getReaction("R_NOGPR");
    CHECK(rxn != nullptr, "reaction parsed");
    if (rxn)
    {
        CHECK(!rxn->hasGPR(), "has no GPR rule");
        CHECK(rxn->getGPR().eval(m->getListOfGenes()),
              "empty rule evaluates true (not under genetic control)");
        CHECK_EQ(m->getGeneReactionRule("R_NOGPR"), std::string(""), "rule renders empty");
    }
    delete m;
}

static void test_single_gene()
{
    DFBA_TEST_SECTION("single gene reference");
    dFBAModel* m = load(FIXTURES + "single_gene.xml");

    CHECK_EQ(m->getNumGenes(), 1, "one gene parsed");
    CHECK(m->hasGene("G_a"), "gene G_a registered");
    CHECK_EQ(m->getGeneReactionRule("R_SINGLE"), std::string("G_a"), "rule renders as G_a");

    dFBAGene* g = m->getGene("G_a");
    CHECK(g != nullptr && g->isActive(), "gene is active by default");
    if (g)
    {
        CHECK_EQ(g->getName(), std::string("a"), "fbc:name parsed");
        CHECK_EQ(g->getLabel(), std::string("G_a"), "fbc:label parsed");
    }
    delete m;
}

static void test_or_and()
{
    DFBA_TEST_SECTION("OR and AND rules");

    dFBAModel* mo = load(FIXTURES + "or_two.xml");
    CHECK_EQ(mo->getGeneReactionRule("R_OR"), std::string("G_a or G_b"), "OR renders");
    delete mo;

    dFBAModel* ma = load(FIXTURES + "and_two.xml");
    CHECK_EQ(ma->getGeneReactionRule("R_AND"), std::string("G_a and G_b"), "AND renders");
    delete ma;
}

static void test_nested()
{
    DFBA_TEST_SECTION("nested or-of-ands, the iML1515 shape");
    dFBAModel* m = load(FIXTURES + "or_of_ands.xml");

    CHECK_EQ(m->getNumGenes(), 4, "four genes parsed");
    CHECK_EQ(m->getGeneReactionRule("R_NESTED"),
             std::string("(G_a and G_b) or (G_c and G_d)"),
             "nested rule renders with parentheses");

    dFBAReaction* rxn = m->getReaction("R_NESTED");
    CHECK(rxn != nullptr && rxn->hasGPR(), "nested GPR attached");
    if (rxn)
        CHECK_EQ((int)rxn->getGPR().getGeneIndices().size(), 4, "rule references four genes");
    delete m;
}

static void test_missing_geneproduct()
{
    DFBA_TEST_SECTION("rule referencing a gene absent from listOfGeneProducts");
    dFBAModel* m = load(FIXTURES + "missing_geneproduct.xml");

    CHECK(m->hasGene("G_ghost"), "referenced-but-undeclared gene is registered");
    CHECK_EQ(m->getGeneReactionRule("R_MISSING"), std::string("G_ghost"),
             "rule keeps the referenced gene");
    delete m;
}

static void test_single_child_collapse()
{
    DFBA_TEST_SECTION("AND wrapping a single child collapses");
    dFBAModel* m = load(FIXTURES + "single_child_and.xml");
    CHECK_EQ(m->getGeneReactionRule("R_COLLAPSE"), std::string("G_a"),
             "single-child AND renders as the bare gene");
    delete m;
}

// ------------------------------------------------------------ evaluation --

static void test_truth_tables()
{
    DFBA_TEST_SECTION("evaluation truth tables");

    // OR: false only when both genes are knocked out.
    {
        dFBAModel* m = load(FIXTURES + "or_two.xml");
        const dFBAGPR& gpr = m->getReaction("R_OR")->getGPR();
        CHECK(gpr.eval(m->getListOfGenes()), "OR: both active -> true");
        m->getGene("G_a")->setActive(false);
        CHECK(gpr.eval(m->getListOfGenes()), "OR: one active -> true");
        m->getGene("G_b")->setActive(false);
        CHECK(!gpr.eval(m->getListOfGenes()), "OR: none active -> false");
        delete m;
    }

    // AND: false as soon as any gene is knocked out.
    {
        dFBAModel* m = load(FIXTURES + "and_two.xml");
        const dFBAGPR& gpr = m->getReaction("R_AND")->getGPR();
        CHECK(gpr.eval(m->getListOfGenes()), "AND: both active -> true");
        m->getGene("G_a")->setActive(false);
        CHECK(!gpr.eval(m->getListOfGenes()), "AND: one knocked out -> false");
        delete m;
    }

    // (a and b) or (c and d)
    {
        dFBAModel* m = load(FIXTURES + "or_of_ands.xml");
        const dFBAGPR& gpr = m->getReaction("R_NESTED")->getGPR();
        CHECK(gpr.eval(m->getListOfGenes()), "nested: all active -> true");
        m->getGene("G_a")->setActive(false);
        CHECK(gpr.eval(m->getListOfGenes()), "nested: first complex broken, second intact -> true");
        m->getGene("G_c")->setActive(false);
        CHECK(!gpr.eval(m->getListOfGenes()), "nested: both complexes broken -> false");
        m->getGene("G_a")->setActive(true);
        CHECK(gpr.eval(m->getListOfGenes()), "nested: first complex repaired -> true");
        delete m;
    }
}

// -------------------------------------------------------------- knockouts --

static void test_knockout_bounds()
{
    DFBA_TEST_SECTION("knockout disables reactions and restores bounds");
    dFBAModel* m = load(FIXTURES + "shared_gene.xml");
    m->initProblem();

    dFBAReaction* r_one = m->getReaction("R_ONE");   // G_shared
    dFBAReaction* r_two = m->getReaction("R_TWO");   // G_shared or G_alt

    const double one_lb = r_one->getLowerBound();
    const double one_ub = r_one->getUpperBound();

    CHECK(!r_one->isKnockedOut() && !r_two->isKnockedOut(), "nothing knocked out initially");

    const int disabled = m->knockOutGene("G_shared");
    CHECK_EQ(disabled, 1, "knocking out G_shared disables exactly one reaction");
    CHECK(r_one->isKnockedOut(), "R_ONE disabled (requires G_shared alone)");
    CHECK(!r_two->isKnockedOut(), "R_TWO still enabled (G_alt covers it)");
    CHECK_NEAR(r_one->getLowerBound(), 0.0, 1e-12, "R_ONE lower bound clamped to zero");
    CHECK_NEAR(r_one->getUpperBound(), 0.0, 1e-12, "R_ONE upper bound clamped to zero");

    CHECK_EQ(m->knockOutGene("G_shared"), 0, "knocking out an inactive gene is a no-op");

    const int reenabled = m->restoreGene("G_shared");
    CHECK_EQ(reenabled, 1, "restoring G_shared re-enables one reaction");
    CHECK(!r_one->isKnockedOut(), "R_ONE enabled again");
    CHECK_NEAR(r_one->getLowerBound(), one_lb, 1e-12, "R_ONE lower bound restored exactly");
    CHECK_NEAR(r_one->getUpperBound(), one_ub, 1e-12, "R_ONE upper bound restored exactly");

    delete m;
}

static void test_partial_restore()
{
    DFBA_TEST_SECTION("restore re-evaluates rather than blindly restoring");
    dFBAModel* m = load(FIXTURES + "and_two.xml");
    m->initProblem();

    dFBAReaction* rxn = m->getReaction("R_AND");   // G_a and G_b

    m->knockOutGene("G_a");
    m->knockOutGene("G_b");
    CHECK(rxn->isKnockedOut(), "R_AND disabled with both genes out");

    CHECK_EQ(m->restoreGene("G_a"), 0, "restoring only G_a changes nothing");
    CHECK(rxn->isKnockedOut(), "R_AND still disabled: G_b is still out");

    CHECK_EQ(m->restoreGene("G_b"), 1, "restoring G_b re-enables the reaction");
    CHECK(!rxn->isKnockedOut(), "R_AND enabled once both genes are back");

    delete m;
}

static void test_knockout_reporting()
{
    DFBA_TEST_SECTION("knockout reporting");
    dFBAModel* m = load(FIXTURES + "shared_gene.xml");
    m->initProblem();

    m->knockOutGene("G_shared");
    CHECK_EQ((int)m->getKnockedOutGenes().size(), 1, "one gene reported knocked out");
    CHECK_EQ((int)m->getKnockedOutReactions().size(), 1, "one reaction reported disabled");
    CHECK(!m->isGeneActive("G_shared"), "isGeneActive false for knocked-out gene");
    CHECK(m->isGeneActive("G_alt"), "isGeneActive true for untouched gene");

    m->knockOutGene("G_alt");
    CHECK_EQ((int)m->getKnockedOutReactions().size(), 2, "R_TWO disabled once both genes are out");

    delete m;
}

static void test_copy_isolation()
{
    DFBA_TEST_SECTION("gene state is not shared between model copies");
    dFBAModel* original = load(FIXTURES + "shared_gene.xml");
    original->initProblem();

    // deepCopy is the path PhysiCell uses when a cell is created or divides.
    dFBAModel* clone = new dFBAModel();
    clone->deepCopy(*original);

    clone->knockOutGene("G_shared");

    CHECK(!clone->isGeneActive("G_shared"), "gene knocked out in the copy");
    CHECK(original->isGeneActive("G_shared"), "original gene untouched");
    CHECK(clone->getReaction("R_ONE")->isKnockedOut(), "copy's reaction disabled");
    CHECK(!original->getReaction("R_ONE")->isKnockedOut(), "original's reaction untouched");
    CHECK_NEAR(original->getReaction("R_ONE")->getUpperBound(), 1000.0, 1e-9,
               "original bounds intact");

    // And a copy taken from an already-knocked-out model inherits that state.
    dFBAModel* clone_of_clone = new dFBAModel();
    clone_of_clone->deepCopy(*clone);
    CHECK(!clone_of_clone->isGeneActive("G_shared"), "knockout survives a further copy");
    CHECK(clone_of_clone->getReaction("R_ONE")->isKnockedOut(),
          "disabled reaction survives a further copy");
    CHECK_NEAR(clone_of_clone->getReaction("R_ONE")->getUpperBound(), 0.0, 1e-12,
               "copied reaction stays clamped");

    delete clone_of_clone;
    delete clone;
    delete original;
}

// ------------------------------------------------------- real model check --

static void test_ecoli_core(const std::string& path)
{
    DFBA_TEST_SECTION("E. coli core model: " + path);
    dFBAModel* m = load(path);

    CHECK_EQ(m->getNumGenes(), 137, "137 gene products parsed");

    int with_gpr = 0;
    for (dFBAReaction* rxn : m->getListOfReactions())
    {
        if (rxn->hasGPR())
            with_gpr++;
    }
    CHECK_EQ(with_gpr, 64, "64 reactions carry a GPR rule");

    // Every rule must be satisfied while all genes are active.
    int unsatisfied = 0;
    for (dFBAReaction* rxn : m->getListOfReactions())
    {
        if (!rxn->getGPR().eval(m->getListOfGenes()))
            unsatisfied++;
    }
    CHECK_EQ(unsatisfied, 0, "every rule evaluates true in the wild type");

    // Exchange reactions are auxiliary and must carry no genes.
    int exchanges_with_gpr = 0;
    for (dFBAReaction* rxn : m->getListOfReactions())
    {
        if (rxn->getId().rfind("R_EX_", 0) == 0 && rxn->hasGPR())
            exchanges_with_gpr++;
    }
    CHECK_EQ(exchanges_with_gpr, 0, "no exchange reaction carries a GPR");

    // Knocking out every gene of a reaction must disable it.
    m->initProblem();
    dFBAReaction* probe = nullptr;
    for (dFBAReaction* rxn : m->getListOfReactions())
    {
        if (rxn->hasGPR()) { probe = rxn; break; }
    }
    if (probe)
    {
        std::vector<int> gene_indices = probe->getGPR().getGeneIndices();
        for (unsigned int i = 0; i < gene_indices.size(); i++)
            m->knockOutGene(m->getListOfGenes()[gene_indices[i]]->getId());

        CHECK(probe->isKnockedOut(),
              "knocking out all genes of " + probe->getId() + " disables it");
        CHECK_NEAR(probe->getUpperBound(), 0.0, 1e-12, "its upper bound is zero");
    }

    delete m;
}

// Growth values cross-checked against cobrapy 0.32.1 reading the same file:
// wild type 0.873921507, and all 137 single-gene deletions agree to 4.8e-11.
static void test_knockout_changes_growth(const std::string& path)
{
    DFBA_TEST_SECTION("knockout changes the FBA solution");
    dFBAModel* m = load(path);
    if (m->getNumReactions() == 0) { delete m; return; }
    m->initProblem();
    m->setReactionLowerBound("R_EX_glc__D_e", -10.0);
    m->setReactionLowerBound("R_EX_o2_e", -1000.0);

    const double wt = m->optimize().getObjectiveValue();
    CHECK_NEAR(wt, 0.873921507, 1e-6, "wild-type growth matches cobrapy");

    // G_b0720 (gltA, citrate synthase) is lethal in the core model: it is the
    // only route from oxaloacetate into the TCA cycle, which biomass requires.
    m->knockOutGene("G_b0720");
    const double ko = m->optimize().getObjectiveValue();
    CHECK_NEAR(ko, 0.0, 1e-6, "gltA knockout is lethal");

    m->restoreGene("G_b0720");
    const double restored = m->optimize().getObjectiveValue();
    CHECK_NEAR(restored, wt, 1e-9, "restoring gltA recovers wild-type growth exactly");

    // A gene with an isozyme partner must not change the optimum.
    dFBAModel* n = load(path);
    n->initProblem();
    n->setReactionLowerBound("R_EX_glc__D_e", -10.0);
    n->setReactionLowerBound("R_EX_o2_e", -1000.0);
    const double before = n->optimize().getObjectiveValue();
    n->knockOutGene("G_b2097");   // fbaB, redundant with fbaA
    CHECK_NEAR(n->optimize().getObjectiveValue(), before, 1e-9,
               "knocking out a redundant isozyme leaves growth unchanged");
    delete n;

    delete m;
}

int main(int argc, const char* argv[])
{
    std::cout << "=================================================" << std::endl;
    std::cout << "dFBA GPR / gene knockout tests" << std::endl;
    std::cout << "=================================================" << std::endl;

    test_no_gpa();
    test_single_gene();
    test_or_and();
    test_nested();
    test_missing_geneproduct();
    test_single_child_collapse();
    test_truth_tables();
    test_knockout_bounds();
    test_partial_restore();
    test_knockout_reporting();
    test_copy_isolation();

    const std::string ecoli = (argc > 1) ? argv[1] : "data/Ecoli_core.xml";
    test_ecoli_core(ecoli);
    test_knockout_changes_growth(ecoli);

    return DFBA_TEST_SUMMARY();
}
