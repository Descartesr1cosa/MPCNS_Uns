#!/usr/bin/env python3
"""Relative L2 differences of two MercuryUns global field files on the same mesh (fields_final.unsb).
usage: compare_states.py <geometry_dir (global.unsb)> ref.unsb test.unsb [test2.unsb ...]
prints, over the fluid cells, the relative L2 error of rho_H, rho_Na, p_H (from U) and over all faces of the induced flux B (name 'B')."""
import struct, sys, numpy as np
def unsb(path):
    d = open(path, 'rb').read(); assert d[:6] == b'UNSBIN'; n, = struct.unpack('<Q', d[8:16]); o = 16; R = {}
    dt = {0: 'i1', 1: '<i4', 2: '<i8', 3: '<f4', 4: '<f8'}
    for _ in range(n):
        nl, = struct.unpack('<I', d[o:o + 4]); o += 4; name = d[o:o + nl].decode(); o += nl
        t = d[o]; o += 4; cols, = struct.unpack('<I', d[o:o + 4]); o += 4; rows, = struct.unpack('<Q', d[o:o + 8]); o += 8
        a = np.frombuffer(d, dtype=dt[t], count=rows * cols, offset=o).copy(); o += rows * cols * a.itemsize; o += (8 - (rows * cols * a.itemsize) % 8) % 8
        R[name] = a.reshape(rows, cols) if cols > 1 else a
    return R
G = unsb(sys.argv[1] + '/global.unsb'); fluid = G['cell_material'] == 0
def fields(fn):
    F = unsb(fn)
    UH = F['f.U_H'].reshape(-1, 5)[fluid]; UN = F['f.U_Na'].reshape(-1, 5)[fluid]
    p = (5 / 3 - 1) * (UH[:, 4] - 0.5 * (UH[:, 1:4] ** 2).sum(1) / UH[:, 0])
    return UH[:, 0], UN[:, 0], p, np.asarray(F['f.B']).ravel()
ref = fields(sys.argv[2])
print("%-40s %10s %10s %10s %10s" % ("file", "rho_H", "rho_Na", "p_H", "B flux"))
for fn in sys.argv[3:]:
    t = fields(fn)
    print("%-40s" % fn.split('/')[-3] + " ".join("%10.2e" % (np.linalg.norm(a - b) / np.linalg.norm(b)) for a, b in zip(t, ref)))
