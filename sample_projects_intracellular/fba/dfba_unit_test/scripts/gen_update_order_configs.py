#!/usr/bin/env python3
"""Write the configs of the update-order test: 12 cells in one voxel, four types that differ in volume and glucose Vmax,
created in different orders (the row order of the cell CSV is the order of all_cells), each at 1 and 4 threads.

types: A 1.3 um3, Vmax glc 8 | B 2.6 um3, 8 | C 1.3 um3, 24 | N 1.3 um3, 0 (cannot take glucose); three cells of each type.
Base: config/voxel_sharing_test/Unit_test_voxel_sharing_mixed.xml (isolated voxels for glucose, 0.165 mM, dt 0.01 min, 0.1 min).
usage: gen_update_order_configs.py [config_dir] [n_random_orders] [seed]
Writes <config_dir>/Unit_test_order_<order>_t<threads>.xml and cells_order_<order>.csv."""
import os, random, re, sys
cdir = sys.argv[1] if len(sys.argv) > 1 else 'config'
nrand = int(sys.argv[2]) if len(sys.argv) > 2 else 3
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
src = open(f'{cdir}/voxel_sharing_test/Unit_test_voxel_sharing_mixed.xml').read()
m1 = re.search(r'( *)<cell_definition name="ecoli1" ID="0">.*?</cell_definition>\n', src, re.S)
m2 = re.search(r'( *)<cell_definition name="ecoli2" ID="1">.*?</cell_definition>\n', src, re.S)
assert m1 and m2 and m2.start() == m1.end(), 'unexpected layout of the mixed config'
types = {'A': (1.3, '8.0'), 'B': (2.6, '8.0'), 'C': (1.3, '24.0'), 'N': (1.3, '0.0')}
defs = ''
for i, (t, (vol, vmax)) in enumerate(types.items()):
    b = m1.group(0).replace('name="ecoli1" ID="0"', f'name="{t}" ID="{i}"')
    b = b.replace('<cell_adhesion_affinity name="ecoli1">', f'<cell_adhesion_affinity name="{t}">')
    b = re.sub(r'(<volume>\s*<total units="micron\^3">)[\d.]+', lambda m: m.group(1) + str(vol), b)
    b = re.sub(r'(<exchange substrate="glucose">.*?<Vmax units="mmol/gDW cell/h">)[\d.]+', lambda m: m.group(1) + vmax, b, flags=re.S)
    defs += b
cfg = src[:m1.start()] + defs + src[m2.end():]
assert '<cell_positions type="csv" enabled="false">' in cfg
cfg = cfg.replace('<cell_positions type="csv" enabled="false">', '<cell_positions type="csv" enabled="true">')
rows = [t for t in types for _ in range(3)]
rng = random.Random(seed)
orders = {'identity': rows[:], 'reversed': rows[::-1]}
for k in range(nrand):
    p = rows[:]; rng.shuffle(p); orders[f'rand{k}'] = p
for name, p in orders.items():
    with open(f'{cdir}/cells_order_{name}.csv', 'w', newline='\n') as f:      # LF: a CRLF file makes PhysiCell reject every cell type
        f.write('x,y,z,type\n' + ''.join(f'0,0,0,{t}\n' for t in p))
    c = cfg.replace('<filename>cells.csv</filename>', f'<filename>cells_order_{name}.csv</filename>')
    for th in (1, 4):
        c2 = c.replace('<omp_num_threads>1</omp_num_threads>', f'<omp_num_threads>{th}</omp_num_threads>')
        c2 = c2.replace('<folder>output_voxel_mixed</folder>', f'<folder>output_order_{name}_t{th}</folder>')
        open(f'{cdir}/Unit_test_order_{name}_t{th}.xml', 'w').write(c2)
print(' '.join(f'{k}={"".join(v)}' for k, v in orders.items()))
