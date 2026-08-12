# Installation

## Clone the repository

```bash
git clone https://github.com/PhysiCelldFBA/PhysiCelldFBA.git
cd PhysiCelldFBA
```

PhysiCelldFBA is distributed as a PhysiCell source tree. A sample target copies the selected project's `main.cpp`, `custom_modules/`, configuration files, analysis scripts, and Makefile into the repository root; a subsequent `make` compiles it.

## Build a dFBA sample

The mass-conservation example is the best installation check:

```bash
make dfba_unit_test
make
```

On the first build, the Makefile automatically runs:

```bash
python3 beta/setup_fba.py
```

The script downloads libSBML and Coin-OR CLP into `addons/dFBA/ext/`. Messages saying that a populated package directory is being skipped are normal on later builds.

!!! warning "Selecting a sample changes the root project"

    PhysiCell sample targets copy files into the root `config/` and `custom_modules/` directories and replace the root `main.cpp` and `Makefile`. Save unrelated work first. Use a clean checkout or a dedicated Git worktree when you want to preserve an existing root project.

## Confirm the executable

For the unit test, a successful build creates:

```text
./dfba-unit-test
```

The other documented targets and executables are:

| Example | Root Make target | Executable |
| --- | --- | --- |
| Mass-conservation test | `dfba_unit_test` | `dfba-unit-test` |
| E. coli acetate switch | `ecoli-acetic-switch-sample` | `ecoli-dfba` |
| E. coli colony | `bacterial-colony` | `bacterial-colony` |
| Cancer core metabolism | `cancer-core-metabolism` | `cancer-core-metabolism` |
| Microbial cross-feeding | `crossfeeding` | `crossfeeding` |
| Metabolism-driven motility | `metabolic_driven_motility` | `ecoli-dfba-motility` |

## Manual dependency recovery

If an interrupted download leaves an incomplete dependency directory, the setup script may consider it installed because the directory is non-empty. Inspect:

```bash
ls addons/dFBA/ext/coin-or/include/coin/CoinPackedMatrix.hpp
ls addons/dFBA/ext/libsbml/include/sbml/SBMLTypes.h
```

If either file is missing, remove only the incomplete package directory and rerun `python3 beta/setup_fba.py`. Do not remove the entire repository or a broad parent directory.

For compiler, linker, SBML, or solver errors, continue with [Troubleshooting](../extending/troubleshooting.md).
