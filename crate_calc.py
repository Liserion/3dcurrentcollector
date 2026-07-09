#!/usr/bin/env python3
"""
crate_calc.py — compute the correct BC currents for a desired C-rate.

Measures the cathode volume and boundary areas DIRECTLY from the exodus mesh
(never trust stale comments!), then applies:

    I_top(xC)    = x * (1-eps) * V_cathode / (3600 * A_top * (1 - t0))
    I_cat_cc(xC) = I_top * A_top / A_cat_cc          (balanced solid-side flux)

Definition: 1C = the current that fills the full stoichiometry range (cs: 0->1)
of the active material in exactly 3600 s (time unit of the nondim system = 1 s).
Note: runs starting at SOC=0.5 hold only half that capacity, so a 1C discharge
of the remaining window lasts ~1800 s.

Usage:
    python3 crate_calc.py macro_in.e --crate 1.0
    python3 crate_calc.py macro_in.e --crate 0.5 --eps 0.2 \
        --cathode-block cathode --top top --catcc cat_cc

Requires: numpy, netCDF4  (pip install numpy netCDF4)
"""
import argparse
import sys

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    sys.exit("needs netCDF4:  pip install netCDF4")

# Exodus II side -> face nodes for TET4 (1-based local node ids)
TET4_SIDES = {1: (0, 1, 3), 2: (1, 2, 3), 3: (0, 3, 2), 4: (0, 2, 1)}


def names(ds, var):
    if var not in ds.variables:
        return []
    out = []
    for row in ds.variables[var][:]:
        s = b"".join(row.compressed() if np.ma.isMaskedArray(row) else row)
        out.append(s.decode().strip("\x00").strip())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mesh", help="exodus mesh file (e.g. macro_in.e)")
    ap.add_argument("--crate", type=float, default=1.0, help="desired C-rate (default 1.0)")
    ap.add_argument("--eps", type=float, default=0.2, help="cathode porosity (default 0.2)")
    ap.add_argument("--cathode-block", default="cathode", help="cathode block name (default 'cathode')")
    ap.add_argument("--top", default="top", help="top sideset name (default 'top')")
    ap.add_argument("--catcc", default="cat_cc", help="current-collector sideset name (default 'cat_cc')")
    ap.add_argument("--ce0", type=float, default=0.0874891,
                    help="initial electrolyte concentration for t0 (default 0.0874891)")
    a = ap.parse_args()

    ds = Dataset(a.mesh)
    x = ds.variables["coordx"][:]
    y = ds.variables["coordy"][:]
    z = ds.variables["coordz"][:]
    P = np.column_stack([x, y, z])

    # ---- cathode volume ------------------------------------------------
    eb_names = names(ds, "eb_names")
    vol = 0.0
    blk_found = None
    for i, bn in enumerate(eb_names, start=1):
        if bn.lower() != a.cathode_block.lower():
            continue
        blk_found = bn
        conn = ds.variables[f"connect{i}"][:] - 1  # 0-based
        if conn.shape[1] != 4:
            sys.exit(f"block '{bn}' is not TET4 (has {conn.shape[1]} nodes/elem)")
        p0, p1, p2, p3 = (P[conn[:, k]] for k in range(4))
        vol += np.abs(np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), p3 - p0)).sum() / 6.0
    if blk_found is None:
        sys.exit(f"block '{a.cathode_block}' not found; blocks: {eb_names}")

    # ---- sideset areas ---------------------------------------------------
    ss_names = names(ds, "ss_names")
    # map global element id -> (block connectivity row)
    n_blocks = int(ds.dimensions["num_el_blk"].size)
    conns, offsets = [], [0]
    for i in range(1, n_blocks + 1):
        c = ds.variables[f"connect{i}"][:] - 1
        conns.append(c)
        offsets.append(offsets[-1] + c.shape[0])

    def elem_nodes(gid):  # gid 0-based global element id
        for b in range(n_blocks):
            if gid < offsets[b + 1]:
                return conns[b][gid - offsets[b]]
        raise IndexError(gid)

    def sideset_area(name):
        for j, sn in enumerate(ss_names, start=1):
            if sn.lower() != name.lower():
                continue
            elems = ds.variables[f"elem_ss{j}"][:] - 1
            sides = ds.variables[f"side_ss{j}"][:]
            area = 0.0
            for e, s in zip(elems, sides):
                nd = elem_nodes(int(e))
                f = [nd[k] for k in TET4_SIDES[int(s)]]
                v1, v2 = P[f[1]] - P[f[0]], P[f[2]] - P[f[0]]
                area += 0.5 * np.linalg.norm(np.cross(v1, v2))
            return area
        return None

    def constructed_area(name):
        """Mimic MOOSE's construct_side_list_from_node_list = true:
        every element side whose 3 nodes all belong to the nodeset is added
        (this is what the flux BCs and AreaPostprocessor actually integrate
        over, and it can exceed the raw geometric sideset — on curved
        surfaces, internal faces whose nodes all lie on the surface get
        included too, from both neighboring elements)."""
        ns_names_ = names(ds, "ns_names")
        ns_ids = ds.variables["ns_prop1"][:].tolist() if "ns_prop1" in ds.variables else []
        ss_ids = ds.variables["ss_prop1"][:].tolist() if "ss_prop1" in ds.variables else []
        # target nodeset: by name, or (gmsh exports leave nodesets unnamed)
        # by the ID of the equally-named sideset
        want_id = None
        for j, sn in enumerate(ss_names):
            if sn.lower() == name.lower() and j < len(ss_ids):
                want_id = ss_ids[j]
        for j, nn in enumerate(ns_names_, start=1):
            match = nn.lower() == name.lower() or (
                nn == "" and want_id is not None
                and j - 1 < len(ns_ids) and ns_ids[j - 1] == want_id)
            if not match:
                continue
            nodes = set((ds.variables[f"node_ns{j}"][:] - 1).tolist())
            area = 0.0
            for conn in conns:
                inset = np.isin(conn, list(nodes))
                for s, idx in TET4_SIDES.items():
                    mask = inset[:, idx].all(axis=1)
                    if not mask.any():
                        continue
                    f = conn[mask][:, idx]
                    v1 = P[f[:, 1]] - P[f[:, 0]]
                    v2 = P[f[:, 2]] - P[f[:, 0]]
                    area += 0.5 * np.linalg.norm(np.cross(v1, v2), axis=1).sum()
            return area
        return None

    def boundary_area(name):
        geo = sideset_area(name)
        eff = constructed_area(name)
        if geo is None and eff is None:
            sys.exit(f"boundary '{name}' found neither as sideset ({ss_names}) "
                     f"nor as nodeset ({names(ds, 'ns_names')})")
        use = eff if eff is not None else geo
        return use, geo, eff

    a_top, a_top_geo, a_top_eff = boundary_area(a.top)
    a_cc, a_cc_geo, a_cc_eff = boundary_area(a.catcc)

    # ---- currents --------------------------------------------------------
    t0 = 0.0107907 + 1.48837e-4 * a.ce0
    i_top_1c = (1 - a.eps) * vol / (3600.0 * a_top * (1 - t0))
    i_top = a.crate * i_top_1c
    i_cc = i_top * a_top / a_cc

    def fmt(v):
        return f"{v:.6e}" if v is not None else "     --     "

    print(f"mesh:            {a.mesh}")
    print(f"V_cathode      = {vol:.6e}")
    print(f"A_top          = {a_top:.6e}   (geometric sideset {fmt(a_top_geo)}, "
          f"MOOSE-constructed {fmt(a_top_eff)})")
    print(f"A_cat_cc       = {a_cc:.6e}   (geometric sideset {fmt(a_cc_geo)}, "
          f"MOOSE-constructed {fmt(a_cc_eff)})")
    print(f"                 using the MOOSE-constructed values — that is what the")
    print(f"                 flux BCs integrate over (construct_side_list_from_node_list)")
    print(f"                 ratio A_top/A_cat_cc = {a_top/a_cc:.5f}")
    print(f"porosity eps   = {a.eps},  t0 = {t0:.6f}")
    print(f"I_top(1C)      = {i_top_1c:.7f}")
    print("-" * 56)
    print(f"===>  {a.crate}C:   I_top = {i_top:.7f}    I_cat_cc = {i_cc:.7f}")
    print("-" * 56)
    print("input file values:")
    print(f"  flux_c    (top):     I = {i_top:.7f}")
    print(f"  flux_phi1 (cat_cc):  I = {i_cc:.7f}")
    print(f"  flux_phi2 (top):     I = {i_top:.7f}")
    print("or as MOOSE CLI overrides:")
    print(f"  BCs/flux_c/I={i_top:.7f} BCs/flux_phi1/I={i_cc:.7f} BCs/flux_phi2/I={i_top:.7f}")
    print("\nreminder: 1C fills the FULL cs range 0->1 in 3600 s; a run starting")
    print("at SOC=0.5 discharges its remaining window in ~1800 s at 1C.")


if __name__ == "__main__":
    main()
