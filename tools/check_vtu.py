#!/usr/bin/env python3
"""Independent check of a .vtu written by uns::io::write_vtu (no VTK needed).

Parses the appended raw data, then recomputes every polyhedron's volume from its
(outward) face loops with the divergence theorem and compares with the CellData
array named by --volume.  Run:  check_vtu.py file.vtu [--volume volume]
"""
import re, struct, sys
import numpy as np

def read_vtu(path):
    raw = open(path, 'rb').read()
    head, _, _ = raw.partition(b'<AppendedData')
    app = raw[raw.index(b'_', raw.index(b'<AppendedData')) + 1:]
    text = head.decode()
    arrays = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="([^"]+)" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"/>', text):
        typ, name, nc, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        n = struct.unpack_from('<Q', app, off)[0]
        dt = {'Float64': '<f8', 'Int64': '<i8', 'UInt8': 'u1'}[typ]
        a = np.frombuffer(app, dtype=dt, count=n // np.dtype(dt).itemsize, offset=off + 8)
        arrays[name] = a.reshape(-1, nc) if nc > 1 else a
    return arrays

def main():
    path = sys.argv[1]
    volname = sys.argv[sys.argv.index('--volume') + 1] if '--volume' in sys.argv else 'volume'
    a = read_vtu(path)
    pts, faces, foffs, types = a['Points'], a['faces'], a['faceoffsets'], a['types']
    assert (types == 42).all()
    vols = []
    start = 0
    for end in foffs:
        k = start
        nf = faces[k]; k += 1
        loops = []
        for _ in range(nf):
            m = faces[k]; loops.append(faces[k + 1:k + 1 + m]); k += 1 + m
        assert k == end
        # Warped faces have no unique volume; the convention (docs/Conventions.md) triangulates
        # every face about its area-weighted fan centroid, V = (1/3) sum_tri a_tri . (c_f - x0).
        v = 0.0
        for L in loops:
            p = pts[L]; q = np.roll(p, -1, axis=0); p0 = p.mean(axis=0)
            w = 0.5 * np.linalg.norm(np.cross(p - p0, q - p0), axis=1)
            cf = ((p0 + p + q) / 3.0 * w[:, None]).sum(axis=0) / w.sum() if w.sum() > 0 else p0
            area = 0.5 * np.cross(p - cf, q - cf)
            v += np.dot(area.sum(axis=0), cf) / 3.0
        vols.append(v)
        start = end
    vols = np.array(vols)
    ref = a[volname]
    err = np.abs(vols - ref).max() / np.abs(ref).max()
    print(f'{path}: cells={len(vols)} min_vol={vols.min():.4g} max rel volume mismatch vs "{volname}" = {err:.2e}')
    sys.exit(0 if err < 1e-9 and (vols > 0).all() else 1)

main()
