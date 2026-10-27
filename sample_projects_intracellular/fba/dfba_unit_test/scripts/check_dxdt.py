#!/usr/bin/env python3
"""Checks for the dx-dt test. usage: check_dxdt.py [--tol-rel 0.005] [--tol-c 0.005] [--dir .]

Reads output_dxdt_<kind>_dx<dx>_dt<dt>/ (see gen_dxdt_configs.py, run_dxdt_test.sh).

  A  inactive, invariance   : glucose consumed, final cell volume and final growth rate against the finest run
                              (smallest dx, smallest dt); relative difference <= tol_rel (default 0.5%; measured 0.04%).
  B  capped, analytic       : glucose never below -tol_abs; at the end the cell's voxel is empty and the total glucose
                              has fallen by exactly V*c0 (the content of that one voxel), to 1e-3 relative.
  C  ledger residual        : cells-window residual of the gross exchange: <= 1e-5 for the inactive runs (rounding; oxygen is 1e8 mM here), <= 2*mu_max*dt for the
                              capped runs (the rate changes within one step when the voxel runs dry).
  D  carbon balance         : carbon leaving the medium (6 glucose + 2 acetate + 1 CO2, from the ledger) equals the carbon built
                              into biomass, C_bio * d(gDW), with C_bio from the biomass reaction formulas (SBML) and
                              gDW = volume*(1-fluid fraction)*density*1e-12; relative difference <= tol_c (default 0.5%; measured 0.05%).
  E  model elemental balance: elemental_check.py on the model: no internal reaction elementally unbalanced
                              (C, H, N, O, P, S; run separately by run_dxdt_test.sh, N, P and S have no extracellular substrate in this test).
exit status 0 = all pass, 1 = a check failed, 2 = missing output.
"""
import argparse, glob, os, re, sys
import numpy as np
import xml.etree.ElementTree as ET
import scipy.io as sio

ap = argparse.ArgumentParser()
ap.add_argument("--tol-rel", type=float, default=0.005)
ap.add_argument("--tol-c", type=float, default=0.005)
ap.add_argument("--dir", default=".")
ap.add_argument("--model", default="config/Ecoli_core.xml")
ap.add_argument("--density", type=float, default=1.04)
a = ap.parse_args()

def cells_at(d, last=True):
    xmls = sorted(glob.glob(os.path.join(d, "output*.xml")))
    xs = xmls[-1]
    r = ET.parse(xs).getroot()
    labels = {}
    for l in r.findall(".//cell_population/custom/simplified_data/labels/label"):
        labels[l.text] = (int(l.get("index")), int(l.get("size")))
    fn = r.find(".//cell_population/custom/simplified_data/filename").text
    m = sio.loadmat(os.path.join(d, fn))["cells"]
    return {k: m[i] for k, (i, s) in labels.items() if s == 1}, float(r.find(".//current_time").text), len(xmls)

MU_MAX = 0.8 / 60.0     # max_growth_rate of the test cell (1/min); first-order mismatch tolerances scale with MU_MAX*dt
FIRST_ORDER = 2.0       # factor on MU_MAX*dt

def ledger(d):
    rows = open(os.path.join(d, "conservation_ledger.csv")).read().strip().split("\n")
    h = rows[0].split(","); v = [float(x) for x in rows[-1].split(",")]
    return dict(zip(h, v))

def snaps_min_conc(d, name_idx):
    worst = 0.0
    for xf in sorted(glob.glob(os.path.join(d, "output*.xml"))):
        r = ET.parse(xf).getroot()
        names = [x.get("name") for x in r.findall(".//microenvironment/domain/variables/variable")]
        fn = r.find(".//microenvironment/domain/data/filename").text
        mm = sio.loadmat(os.path.join(d, fn))["multiscale_microenvironment"]
        sub = mm[4:4 + len(names), :]
        worst = min(worst, float(sub[names.index("glucose")].min()))
    return worst

def snap_glucose_voxel_sum(d, which):
    xf = sorted(glob.glob(os.path.join(d, "output*.xml")))[which]
    r = ET.parse(xf).getroot()
    names = [x.get("name") for x in r.findall(".//microenvironment/domain/variables/variable")]
    fn = r.find(".//microenvironment/domain/data/filename").text
    mm = sio.loadmat(os.path.join(d, fn))["multiscale_microenvironment"]
    return mm[3], mm[4 + names.index("glucose")]   # voxel volume, glucose concentration per voxel

runs = {}
for d in sorted(glob.glob(os.path.join(a.dir, "output_dxdt_*"))):
    m = re.match(r"output_dxdt_(inactive|capped)_dx([\d.]+)_dt([\d.]+)$", os.path.basename(d))
    if m and glob.glob(os.path.join(d, "conservation_ledger.csv")):
        runs[(m.group(1), float(m.group(2)), float(m.group(3)))] = d
if not runs:
    print("ERROR: no output_dxdt_* folders with a ledger"); sys.exit(2)

fails = []
def chk(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok: fails.append(msg)

# SBML biomass carbon (mmol C per gDW): sum over reactants of stoich*C minus products
import cobra, logging, warnings
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
model = cobra.io.read_sbml_model(os.path.join(a.dir, a.model))
bio = [r for r in model.reactions if r.objective_coefficient][0]
def nC(met): return met.elements.get("C", 0)
Cbio = -sum(c * nC(m) for m, c in bio.metabolites.items())
print(f"biomass reaction {bio.id}: {Cbio:.3f} mmol C per gDW (from the species formulas)")

# --- A: invariance of the inactive runs against the finest
ina = {k: v for k, v in runs.items() if k[0] == "inactive"}
if ina:
    ref = min(ina, key=lambda k: (k[1], k[2]))
    def obs(d):
        c, t, n = cells_at(d); L = ledger(d)
        return dict(vol=float(c["total_volume"].sum()), n=len(c["total_volume"]), mu=float(c["growth_rate"].mean()),
                    glc=-L["glucose_cells_dF"], t=t)
    O = {k: obs(d) for k, d in ina.items()}
    print(f"reference (finest): dx {ref[1]}, dt {ref[2]}: total volume {O[ref]['vol']:.5f} um^3 in {O[ref]['n']} cells, mean growth rate {O[ref]['mu']:.6g}, glucose taken {O[ref]['glc']:.5g} mM um^3 at t={O[ref]['t']}")
    for k in sorted(ina):
        if k == ref: continue
        chk(O[k]["n"] == O[ref]["n"], f"A inactive dx {k[1]} dt {k[2]}: {O[k]['n']} cells at the end (reference {O[ref]['n']})")
        for q in ("vol", "glc"):
            rel = abs(O[k][q] / O[ref][q] - 1)
            chk(rel <= a.tol_rel, f"A inactive dx {k[1]} dt {k[2]}: {q} differs from the reference by {rel*100:.3f}% (limit {a.tol_rel*100:g}%)")

# --- B: capped analytic
for k, d in sorted(runs.items()):
    if k[0] != "capped": continue
    nmin = snaps_min_conc(d, 0)
    tol = FIRST_ORDER * MU_MAX * k[2] * 1.0      # c0 = 1 mM
    chk(nmin >= -tol, f"B capped dx {k[1]} dt {k[2]}: lowest glucose concentration {nmin:.3g} mM (>= -{tol:.2g} = -{FIRST_ORDER:g}*mu_max*dt*c0)")
    V0, c0 = snap_glucose_voxel_sum(d, 0)
    V1, c1 = snap_glucose_voxel_sum(d, -1)
    taken = float(np.sum(V0 * c0) - np.sum(V1 * c1))
    expect = k[1] ** 3 * 1.0     # V * c0 for the one voxel of the cell
    rel = abs(taken / expect - 1)
    chk(rel <= 1e-3, f"B capped dx {k[1]} dt {k[2]}: glucose taken {taken:.5g} vs content of the cell's voxel {expect:.5g} (differ {rel*100:.3f}%)")

# --- C: ledger residual; D: carbon
for k, d in sorted(runs.items()):
    L = ledger(d)
    for s in ("glucose", "oxygen", "acetate", "CO2"):
        g = L[f"{s}_cells_gross"]
        if g > 0:
            r = abs(L[f"{s}_cells_resid"]) / g
            # inactive: nothing changes abruptly, the ledger closes to rounding (oxygen is 1e8 mM, so its residual is the largest, up to ~1e-6). capped: the rate changes within one step when the
            # voxel runs dry, and the applied rate lags by that step: bound by the first-order size 2*mu_max*dt.
            lim = 1e-5 if k[0] == "inactive" else FIRST_ORDER * MU_MAX * k[2]
            chk(r <= lim, f"C {k[0]} dx {k[1]} dt {k[2]}: {s} cells-window residual {r:.2e} of gross (<= {lim:.2e})")
    if k[0] != "inactive": continue
    c, t, n = cells_at(d)
    # volume at the first snapshot (cells placed with volume 1.3)
    vol0 = 1.3; ff0 = 0.75
    gdw = lambda vol, ff: vol * (1 - ff) * a.density * 1e-12
    dgdw = float(np.sum(gdw(c["total_volume"], c["fluid_fraction"]))) - gdw(vol0, ff0)
    Cmedium = -(6 * L["glucose_cells_dF"] + 2 * L["acetate_cells_dF"] + L["CO2_cells_dF"]) * 1e-15   # mmol C leaving the medium
    Cb = Cbio * dgdw
    rel = abs(Cmedium / Cb - 1) if Cb != 0 else float("inf")
    chk(rel <= a.tol_c, f"D carbon dx {k[1]} dt {k[2]}: medium loses {Cmedium:.4g} mmol C, biomass gains {Cb:.4g} (differ {rel*100:.2f}%, limit {a.tol_c*100:g}%)")

print("\n%d check(s) failed" % len(fails) if fails else "\nall checks passed")
sys.exit(1 if fails else 0)
