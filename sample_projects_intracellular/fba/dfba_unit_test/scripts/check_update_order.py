#!/usr/bin/env python3
"""The result must not depend on the order of the cells or on the number of threads.

Compares every output_order_<order>_t<threads> folder with output_order_identity_t1: per cell type the sorted glucose/oxygen
exchange flux, growth rate and volume, and the whole glucose and oxygen fields, at every snapshot. Also requires that no concentration
goes below -tolerance.

Before the fixes in core/PhysiCell_cell_container.cpp (bit-packed std::vector<bool> written in a parallel loop: cells lose their
"due for update" flag and are never updated) and BioFVM/BioFVM_basic_agent.cpp (unsynchronised update of a voxel's density vector),
the 4-thread runs differ from the 1-thread run by up to 1 mmol/gDW/h and 0.1 mM and the shared voxel goes negative.

usage: check_update_order.py [tolerance] [negative_tolerance_mM]
  tolerance default 1e-9 (fluxes, mM, um^3); negative_tolerance default 1e-4 mM, the same as check_voxel_nonnegative.py (a growth-related
  overshoot of about 1e-5 mM is expected and is not an order effect)
exit status: 0 = all runs agree, 1 = not"""
import glob, os, sys
import numpy as np
import scipy.io as sio
import xml.etree.ElementTree as ET

def load(d):
    snaps = []
    for xf in sorted(glob.glob(f'{d}/output*.xml')):
        r = ET.parse(xf).getroot()
        lab = {l.text: int(l.get('index')) for l in r.findall('.//cell_population//labels/label')}
        c = sio.loadmat(xf.replace('.xml', '_cells.mat'))['cells']
        mm = sio.loadmat(os.path.join(d, r.find('.//microenvironment/domain/data/filename').text))['multiscale_microenvironment']
        names = [v.get('name') for v in r.findall('.//microenvironment/domain/variables/variable')]
        typ = c[lab['cell_type']].astype(int)
        cells = {k: {int(T): np.sort(c[lab[k], typ == T]) for T in np.unique(typ)}
                 for k in ('oxygen_flux', 'glucose_flux', 'growth_rate', 'total_volume')}
        snaps.append((cells, {n: mm[4 + i] for i, n in enumerate(names)}))
    return snaps

def main():
    tol = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-9
    ntol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
    ref_dir = 'output_order_identity_t1'
    if len(glob.glob(f'{ref_dir}/output*.xml')) < 3:
        print(f'ERROR: no output in {ref_dir}'); return 2
    ref = load(ref_dir)
    bad = 0
    for d in sorted(p for p in glob.glob('output_order_*') if os.path.isdir(p)):
        cur = load(d)
        if len(cur) != len(ref):
            print(f'{d}: {len(cur)} snapshots, reference has {len(ref)}'); bad += 1; continue
        dev = 0.0; low = 0.0
        for (c1, f1), (c2, f2) in zip(ref, cur):
            for k in c1:
                for T in c1[k]:
                    dev = max(dev, float(np.abs(c1[k][T] - c2[k][T]).max()))
            for n in ('glucose', 'oxygen'):
                dev = max(dev, float(np.abs(f1[n] - f2[n]).max()))
            low = min(low, float(f2['glucose'].min()), float(f2['oxygen'].min()))
        ok = dev <= tol and low >= -ntol
        bad += not ok
        print(f'{d:34s} max difference to identity_t1 {dev:9.3e}   lowest concentration {low:9.3e}   {"ok" if ok else "FAIL"}')
    print('PASS' if not bad else f'FAIL: {bad} run(s) differ from the 1-thread reference by more than {tol:g}')
    return 1 if bad else 0

if __name__ == '__main__':
    sys.exit(main())
