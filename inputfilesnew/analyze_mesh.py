#!/usr/bin/env python3
"""Compare old vs new macro_in.e: block volumes, sideset areas, element quality."""
import numpy as np
from netCDF4 import Dataset
import sys

# local face numbering (side -> node indices) for exodus TET4
TET4_FACES = {1: [0, 1, 3], 2: [1, 2, 3], 3: [0, 3, 2], 4: [0, 2, 1]}

def tet_vol(p):
    return np.linalg.det(np.column_stack([p[1]-p[0], p[2]-p[0], p[3]-p[0]])) / 6.0

def analyze(path):
    d = Dataset(path)
    x = d.variables['coordx'][:]; y = d.variables['coordy'][:]; z = d.variables['coordz'][:]
    coords = np.column_stack([x, y, z])
    print(f"\n=== {path}")
    print(f"  nodes={len(x)}  bbox x[{x.min():.2f},{x.max():.2f}] y[{y.min():.2f},{y.max():.2f}] z[{z.min():.2f},{z.max():.2f}]")

    # block names
    names = [b''.join(c for c in row if c != b'--').decode(errors='ignore').strip('\x00')
             for row in d.variables['eb_names'][:]] if 'eb_names' in d.variables else []
    ss_names = [b''.join(c for c in row if c != b'--').decode(errors='ignore').strip('\x00')
                for row in d.variables['ss_names'][:]] if 'ss_names' in d.variables else []

    blocks = {}
    nblk = d.dimensions['num_el_blk'].size
    for i in range(1, nblk + 1):
        conn = d.variables[f'connect{i}'][:] - 1
        name = names[i-1] if names and names[i-1] else f'block{i}'
        elem_type = d.variables[f'connect{i}'].elem_type
        vols = np.array([tet_vol(coords[e]) for e in conn]) if conn.shape[1] == 4 else None
        blocks[i] = (name, conn)
        if vols is not None:
            zb = coords[conn].reshape(-1,3)[:,2]
            print(f"  block {i} '{name}' ({elem_type}): {len(conn)} elems, vol={vols.sum():.4g}, "
                  f"min|V|={np.abs(vols).min():.3g}, neg vols={np.sum(vols<=0)}, z[{zb.min():.2f},{zb.max():.2f}]")
        else:
            print(f"  block {i} '{name}' ({elem_type}): {len(conn)} elems")

    # global element id -> (block, local idx)
    counts = [len(blocks[i][1]) for i in sorted(blocks)]
    offsets = np.cumsum([0] + counts)

    if 'num_side_sets' in d.dimensions:
        nss = d.dimensions['num_side_sets'].size
        for i in range(1, nss + 1):
            els = d.variables[f'elem_ss{i}'][:] - 1
            sides = d.variables[f'side_ss{i}'][:]
            name = ss_names[i-1] if ss_names and ss_names[i-1] else f'ss{i}'
            area = 0.0
            zs = []
            for e, s in zip(els, sides):
                bi = np.searchsorted(offsets, e, side='right')
                conn = blocks[bi][1]
                nodes = conn[e - offsets[bi-1]]
                f = nodes[TET4_FACES[int(s)]] if len(nodes) == 4 else None
                if f is None:
                    continue
                p = coords[f]
                area += 0.5 * np.linalg.norm(np.cross(p[1]-p[0], p[2]-p[0]))
                zs.extend(coords[f][:,2])
            zs = np.array(zs) if zs else np.array([np.nan])
            print(f"  sideset {i} '{name}': {len(els)} faces, area={area:.4g}, z[{np.nanmin(zs):.2f},{np.nanmax(zs):.2f}]")
    d.close()

for p in sys.argv[1:]:
    analyze(p)
