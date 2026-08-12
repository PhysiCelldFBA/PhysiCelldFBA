# Requirements

PhysiCelldFBA uses the normal PhysiCell C++ build workflow plus two libraries for constraint-based modelling:

- **libSBML 5.20.4** reads SBML models and their Flux Balance Constraints (FBC) annotations.
- **Coin-OR CLP 1.17.10** solves the linear programs used for FBA.

You do not normally install these libraries by hand.

!!! success "Automatic first-build setup"

    Every supplied dFBA Makefile checks for `addons/dFBA/ext/coin-or/include/coin/CoinPackedMatrix.hpp`. If it is absent, the build calls `python3 beta/setup_fba.py`, which downloads both libSBML and Coin-OR CLP into `addons/dFBA/ext/`. Later builds reuse the existing directories.

## System tools

Before building, make sure the following commands are available:

```bash
python3 --version
g++ --version
make --version
```

The compiler must support C++11 and OpenMP. The sample Makefiles use `g++` by default. As in PhysiCell, you can override it with the `PHYSICELL_CPP` environment variable when an OpenMP-capable compiler has a different name.

## Prebuilt dependency platforms

The package manifest currently provides libSBML and Coin-OR archives for:

| Operating system | Architecture |
| --- | --- |
| Linux | x86-64 |
| Windows | 64-bit |
| macOS | Intel x86-64 |
| macOS | Apple Silicon arm64 |

The setup script detects the host platform and selects the corresponding entry from `beta/fba_packages.json`. Other platforms require a manual build of compatible libSBML and Coin-OR libraries with include and library paths matching the project Makefile.

## Network and disk access

The first dFBA build needs:

- Internet access to download the two dependency archives.
- Write access to `addons/dFBA/ext/`.
- Enough space for the unpacked headers and static libraries.

After the dependencies are present, compilation and simulation can run offline.

## Python for optional utilities

The C++ simulation does not use the Python configuration utilities at runtime. The optional YAML workflow under `beta/` needs packages including:

- COBRApy
- pandas
- PyYAML
- lxml
- `physicell-settings`

Install these only if you intend to generate YAML/configuration files with [the helper scripts](../configuration/helper-scripts.md). Individual analysis notebooks and plotting scripts can have additional scientific-Python requirements.

## Model requirement

A dFBA cell needs an SBML model containing usable FBC reaction bounds and an objective reaction. Exchange reaction identifiers referenced in the PhysiCell XML must exist in that model. See [Prepare an SBML model](../configuration/sbml.md).
