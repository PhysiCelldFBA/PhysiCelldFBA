"""Elemental balance of the SBML models: for every non-exchange, non-biomass reaction, cobrapy check_mass_balance()
(element and charge imbalance from the species formulas). usage: elemental_check.py MODEL.xml ..."""
import sys, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
import cobra
for f in sys.argv[1:]:
    m = cobra.io.read_sbml_model(f)
    nof = [x.id for x in m.metabolites if not x.formula]
    ex = {r.id for r in m.exchanges} | {r.id for r in m.demands} | {r.id for r in m.sinks}
    obj = {r.id for r in m.reactions if r.objective_coefficient}
    bio = {r.id for r in m.reactions if 'biomass' in r.id.lower() or 'biomass' in (r.name or '').lower()}
    internal = [r for r in m.reactions if r.id not in ex and r.id not in obj and r.id not in bio]
    bad = {}; unk = 0
    for r in internal:
        if any(not x.formula for x in r.metabolites): unk += 1; continue
        mb = r.check_mass_balance()
        if mb: bad[r.id] = mb
    print(f"{f.split('/')[-1]}: reactions {len(m.reactions)}, exchanges {len(ex)}, biomass/objective {len(obj|bio)}, internal {len(internal)}")
    print(f"   species without formula {len(nof)}/{len(m.metabolites)}; internal reactions not checkable {unk}; checked {len(internal)-unk}; elementally unbalanced {len(bad)}")
    for k, v in list(bad.items())[:6]: print("     ", k, v)
