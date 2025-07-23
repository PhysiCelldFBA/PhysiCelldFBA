#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include "coin/ClpSimplex.hpp"
#include "coin/OsiClpSolverInterface.hpp"
#include "coin/CoinWarmStartBasis.hpp"


using namespace std;

int main() {
    // --- Cold start solve ---
    OsiClpSolverInterface solver;
    std::cout << "Reading LP file..." << std::endl;
    int status = solver.readLp("data/human_gem.lp");
    if (status != 0) {
        std::cerr << "Failed to read LP file." << std::endl;
        return 1;
    }

    std::cout << "Solving (cold start)..." << std::endl;
    auto start_cold = std::chrono::high_resolution_clock::now();
    solver.initialSolve();
    auto end_cold = std::chrono::high_resolution_clock::now();

    if (solver.isProvenOptimal()) {
        std::cout << "Objective value (cold): " << solver.getObjValue() << std::endl;
    } else {
        std::cout << "Cold start solve did not reach optimal solution." << std::endl;
    }

    double cold_time = std::chrono::duration<double>(end_cold - start_cold).count();
    std::cout << "Time (cold): " << cold_time << " seconds\n" << std::endl;

    // --- Get warm start basis ---
    const CoinWarmStart* warm = solver.getWarmStart();
    if (!warm) {
        std::cerr << "Warning: no warm start basis available." << std::endl;
        return 1;
    }
    // Clone it before passing to next solver
    CoinWarmStart* warm_clone = warm->clone();

    // --- Warm start solve ---
    OsiClpSolverInterface solver_warm;
    std::cout << "Reading LP file for warm start..." << std::endl;
    status = solver_warm.readLp("data/human_gem.lp");


    if (status != 0) {
        std::cerr << "Failed to read LP file for warm start." << std::endl;
        delete warm_clone;
        return 1;
    }

    // Set warm start basis (solver takes ownership, do NOT delete warm_clone after this)
    solver_warm.setWarmStart(warm_clone);

    std::cout << "Solving (warm start)..." << std::endl;
    auto start_warm = std::chrono::high_resolution_clock::now();
    solver_warm.initialSolve();
    auto end_warm = std::chrono::high_resolution_clock::now();

    if (solver_warm.isProvenOptimal()) {
        std::cout << "Objective value (warm): " << solver_warm.getObjValue() << std::endl;
    } else {
        std::cout << "Warm start solve did not reach optimal solution." << std::endl;
    }

    double warm_time = std::chrono::duration<double>(end_warm - start_warm).count();
    std::cout << "Time (warm): " << warm_time << " seconds" << std::endl;

    return 0;
}