# Third-party software and data

PhysiCelldFBA itself is distributed under the BSD 3-Clause licence (`LICENSE`). It extends PhysiCell and BioFVM
(BSD 3-Clause; `licenses/PhysiCell.txt`, `licenses/BioFVM.txt`) and bundles pugixml and MaBoSS code under their own
licences (`licenses/pugixml.txt`, `licenses/MaBoSS.txt`).

## Libraries installed by `beta/setup_fba.py` (not distributed in this repository)

| Component | Licence | Use |
| --- | --- | --- |
| COIN-OR Clp 1.17.x | Eclipse Public License 2.0 (EPL-2.0) | linear programming solver |
| libSBML 5.20.x | GNU Lesser General Public License 2.1 or later (LGPL-2.1+) | reading SBML models |

The libraries are downloaded and built by the setup script; their licence files come with their sources.

## Metabolic models in the sample projects

The SBML files in `sample_projects_intracellular/fba/*/config/` are redistributed from public sources and keep the
terms of those sources; they are **not** covered by the BSD 3-Clause licence of the code.

| Model | Source | Terms |
| --- | --- | --- |
| E. coli core, iML1515 | BiGG Models (https://bigg.ucsd.edu) | BiGG Models licence (https://bigg.ucsd.edu/license): free for educational, research and non-profit use; commercial use requires permission from the University of California, San Diego. Cite the BiGG publication and the model paper. |
| iCB925 (C. beijerinckii) | BioModels MODEL1507180034 | CC0 public-domain dedication (notice in the SBML file) |
| iMG746 (M. barkeri), MCF7 core model | Supporting information of the original publications (Gonnerman et al. 2013; Zielinski et al. 2017) | terms of the publications; cite the original papers |

See the `dfba_sbml_models` repository for the exact files, their sources and the preparation scripts.
