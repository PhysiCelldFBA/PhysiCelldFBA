#!/usr/bin/env python3
"""Write the voxel-size / time-step configs for the dx-dt test (run_dxdt_test.sh).

Base: config/Unit_test_mass_conversion.xml (one e_coli_core cell in a closed 16 um box).
Only dx (= dy = dz), dt_diffusion and intracellular_dt (kept equal, the rule of the addon), the run length,
the thread count and the output folder change. The box size does not change, so the total glucose is the same at every dx.

  inactive_dx<dx>_dt<dt>  : glucose 10 mM, normal diffusion. The uptake cap should not act, so the outputs must agree across dx and dt.
  capped_dx<dx>_dt<dt>    : glucose 1 mM, glucose D = 0 (the cell's voxel is an isolated pool). The cap must act at the first step and the
                            voxel is emptied, never overdrawn.

usage: gen_dxdt_configs.py <config_dir> [max_time_min]
"""
import re, sys, os
cdir = sys.argv[1]
T = sys.argv[2] if len(sys.argv) > 2 else "120"
base = open(os.path.join(cdir, "Unit_test_mass_conversion.xml")).read()
DX = [0.8, 1.6, 3.2]
DT = [0.005, 0.01, 0.02]

def sub1(pat, rep, s, flags=re.S):
    out, n = re.subn(pat, rep, s, count=1, flags=flags)
    assert n == 1, pat
    return out

def make(kind, dx, dt):
    s = base
    for t in ("dx", "dy", "dz"):
        s = sub1(rf"<{t}>[^<]*</{t}>", f"<{t}>{dx}</{t}>", s)
    s = sub1(r"(<max_time units=\"min\">)[^<]*", rf"\g<1>{T}", s)
    s = sub1(r"(<dt_diffusion units=\"min\">)[^<]*", rf"\g<1>{dt}", s)
    s = sub1(r"<intracellular_dt>[^<]*</intracellular_dt>", f"<intracellular_dt>{dt}</intracellular_dt>", s)
    s = sub1(r"<omp_num_threads>[^<]*</omp_num_threads>", "<omp_num_threads>1</omp_num_threads>", s)
    s = sub1(r"(<full_data>\s*<interval units=\"min\">)[^<]*", r"\g<1>30", s)
    s = sub1(r"(<SVG>\s*<interval units=\"min\">)[^<]*", r"\g<1>30", s)
    s = sub1(r"(<SVG>\s*<interval units=\"min\">[^<]*</interval>\s*<enable>)true", r"\g<1>false", s)
    s = sub1(r"(<user_parameters>)", r'\g<1>\n        <dfba_conservation_ledger type="int" units="dimensionless" description="write conservation_ledger.csv">1</dfba_conservation_ledger>', s)
    name = f"dxdt_{kind}_dx{dx}_dt{dt}"
    s = sub1(r"(<save>\s*<folder>)[^<]*", rf"\g<1>output_{name}", s)
    if kind == "capped":
        s = sub1(r'(<variable name="glucose".*?<diffusion_coefficient units="micron\^2/min">)[^<]*', r"\g<1>0", s)
        s = sub1(r'(<variable name="glucose".*?<initial_condition units="mM">)[^<]*', r"\g<1>1", s)
    open(os.path.join(cdir, f"Unit_test_{name}.xml"), "w").write(s)
    return name

names = [make(k, dx, dt) for k in ("inactive", "capped") for dx in DX for dt in DT]
print("\n".join(names))
