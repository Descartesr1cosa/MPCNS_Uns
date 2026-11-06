#!/usr/bin/env python3
"""Convert MercuryUns / MarsUns output (per-rank VTU pieces, polyhedral cells) to a Tecplot BINARY file (.plt, format #!TDV112).

Why: the solvers write VTK_POLYHEDRON cells (face-based connectivity).  ParaView reads them, Tecplot does not ("missing connectivity").
The converter merges the rank pieces of one output step (duplicate nodes on rank interfaces are merged by exact coordinates) and writes
  zone "hex"  ZONETYPE=FEBRICK        all hexahedral cells (every unrefined cubed-sphere / box mesh: the whole Mars and Mercury meshes)
  zone "poly" ZONETYPE=FETETRAHEDRON  every other cell (hanging-node polyhedra of refined meshes), each cut into tetrahedra about its centroid
                                      and its face centres (extra nodes: x,y,z only); cell data are copied to the sub-tetrahedra.
If the mesh is purely hexahedral the file has the single FEBRICK zone.
All cell data of the VTU files are written cell-centred; vectors are split into name_x, name_y, name_z.  Node variables: x, y, z
(nondimensional, as in the solver).  Values are stored as 32-bit floats (--double for 64-bit).

usage:  python3 vtu2tec.py [config.yml]        (default: config.yml next to this script)
All settings are in the yml file (see config.yml here); relative paths in it are relative to the directory of the yml file:
  input_dir   directory with sol_<step>_rXXXXX.vtu (the solver's output_dir), optionally solution.pvd for the solution times
  output_dir  created if missing
  name        output file name WITHOUT extension (empty: sol_<step>).  One step: NAME.plt;  step: all -> NAME_<step>.plt for every step
  step        last | all | a step number (e.g. 15659 or 0015659)
  vars        list of cell variables to write (prefixes with * allowed, e.g. rho_*); empty/omitted = all
  double      true: 64-bit values (default false: 32-bit)
Tecplot: File > Load Data > Tecplot Data Loader (.plt).  Independent check: the file parses with the pure-python `tecio` package.
"""
import fnmatch, glob, os, re, struct, sys
import numpy as np

DT = {'Float64': '<f8', 'Float32': '<f4', 'Int64': '<i8', 'Int32': '<i4', 'UInt8': 'u1'}
SKIP = ('Points', 'connectivity', 'offsets', 'types', 'faces', 'faceoffsets')


def read_vtu(path):
    """Arrays of a vtu written by uns::io::write_vtu (appended raw data, UInt64 headers, no compression)."""
    raw = open(path, 'rb').read()
    i = raw.index(b'<AppendedData')
    head = raw[:i].decode()
    app = raw[raw.index(b'_', i) + 1:]
    out = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="([^"]+)" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"/>', head):
        typ, name, nc, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        n = struct.unpack_from('<Q', app, off)[0]
        a = np.frombuffer(app, dtype=DT[typ], count=n // np.dtype(DT[typ]).itemsize, offset=off + 8)
        out[name] = a.reshape(-1, nc) if nc > 1 else a
    return out


def cell_faces(piece):
    """list (per cell) of face loops (arrays of local node ids)"""
    faces, offs = piece['faces'], piece['faceoffsets']
    res, start = [], 0
    for end in offs:
        k = start
        nf = faces[k]; k += 1
        loops = []
        for _ in range(nf):
            m = faces[k]; loops.append(faces[k + 1:k + 1 + m]); k += 1 + m
        res.append(loops)
        start = end
    return res


def hex_order(loops, pts):
    """FEBRICK node order (8 ids) of a hexahedron given by its 6 quad loops, or None if the cell is not a hexahedron."""
    if len(loops) != 6 or any(len(L) != 4 for L in loops):
        return None
    adj = {}
    for L in loops:
        for i in range(4):
            a, b = int(L[i]), int(L[(i + 1) % 4])
            adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    if len(adj) != 8 or any(len(v) != 3 for v in adj.values()):
        return None
    bot = [int(v) for v in loops[0]]
    top = []
    for a in bot:
        rest = adj[a] - set(bot)
        if len(rest) != 1:
            return None
        top.append(rest.pop())
    if np.dot(np.cross(pts[bot[1]] - pts[bot[0]], pts[bot[3]] - pts[bot[0]]), pts[top[0]] - pts[bot[0]]) < 0:  # right-handed brick
        bot = [bot[0], bot[3], bot[2], bot[1]]
        top = [top[0], top[3], top[2], top[1]]
    return bot + top


def steps_in(indir):
    s = sorted({re.search(r'sol_(\d+)_r', f).group(1) for f in glob.glob(os.path.join(indir, 'sol_*_r*.vtu'))})
    if not s:
        sys.exit(f'no sol_<step>_rXXXXX.vtu files in {indir}')
    return s


def step_times(indir):
    t = {}
    pvd = os.path.join(indir, 'solution.pvd')
    if os.path.exists(pvd):
        for m in re.finditer(r'timestep="([^"]+)"[^>]*file="sol_(\d+)\.pvtu"', open(pvd).read()):
            t[m.group(2)] = float(m.group(1))
    return t


# ------------------------------------------------------------------------------------------------ Tecplot binary (TDV112)
def tstr(s):
    return b''.join(struct.pack('<i', ord(c)) for c in s) + struct.pack('<i', 0)


def write_plt(path, title, varnames, zones, double):
    """zones: list of dict(name, ztype 4|5, nodes (N,3), cellvars [arrays (E,)], conn (E,k) 0-based, time, strand)"""
    nvar = len(varnames)
    vt, dt = (2, '<f8') if double else (1, '<f4')
    out = bytearray(b'#!TDV112')
    out += struct.pack('<ii', 1, 0)                    # byte-order marker, file type FULL
    out += tstr(title) + struct.pack('<i', nvar)
    for n in varnames:
        out += tstr(n)
    for z in zones:
        out += struct.pack('<f', 299.0) + tstr(z['name'])
        out += struct.pack('<ii', -1, z['strand'])     # parent zone, strand id (-1 static)
        out += struct.pack('<d', z['time'])
        out += struct.pack('<ii', -1, z['ztype'])      # zone colour (default), zone type
        out += struct.pack('<i', 1) + b''.join(struct.pack('<i', 0 if v < 3 else 1) for v in range(nvar))  # var location: x,y,z nodal, rest cell-centred
        out += struct.pack('<ii', 0, 0)                # raw local 1-to-1 face neighbours, user-defined face neighbours
        out += struct.pack('<ii', len(z['nodes']), len(z['conn']))   # NumPts, NumElements
        out += struct.pack('<iii', 0, 0, 0)            # cell dimensions (unused for FE zones)
        out += struct.pack('<i', 0)                    # no auxiliary data
    out += struct.pack('<f', 357.0)                    # end of header
    with open(path, 'wb') as f:
        f.write(out)
        for z in zones:
            arrs = [z['nodes'][:, 0], z['nodes'][:, 1], z['nodes'][:, 2]] + list(z['cellvars'])
            f.write(struct.pack('<f', 299.0))
            f.write(struct.pack('<%di' % nvar, *([vt] * nvar)))
            f.write(struct.pack('<i', 0))              # no passive variables
            f.write(struct.pack('<i', 0))              # no variable sharing
            f.write(struct.pack('<i', -1))             # no connectivity sharing
            for a in arrs:
                f.write(struct.pack('<dd', float(np.min(a)) if len(a) else 0.0, float(np.max(a)) if len(a) else 0.0))
            for a in arrs:
                f.write(np.ascontiguousarray(a, dtype=dt).tobytes())
            f.write(np.ascontiguousarray(z['conn'], dtype='<i4').tobytes())


# ------------------------------------------------------------------------------------------------ polyhedron -> tetrahedra
def tets_of_cell(loops_g, xyz, extra, key2id, nbase):
    """tetrahedra (global node ids) of one polyhedron: centroid + face centre + edge; `extra` collects new node coordinates"""
    ids = np.unique(np.concatenate(loops_g))
    xc = xyz[ids].mean(axis=0)
    ic = nbase + len(extra); extra.append(xc)
    tets = []
    def add(a, b, c, d):
        pa, pb, pc, pd = (xc if i == ic else (extra[i - nbase] if i >= nbase else xyz[i]) for i in (a, b, c, d))
        v = np.dot(np.cross(pb - pa, pc - pa), pd - pa)
        tets.append([a, b, c, d] if v >= 0 else [a, c, b, d])
    for L in loops_g:
        L = [int(v) for v in L]
        if len(L) == 3:
            add(ic, L[0], L[1], L[2])
            continue
        key = tuple(sorted(L))
        icf = key2id.get(key)
        if icf is None:
            icf = nbase + len(extra); extra.append(xyz[L].mean(axis=0)); key2id[key] = icf
        for i in range(len(L)):
            add(ic, icf, L[i], L[(i + 1) % len(L)])
    return tets


def convert(indir, outdir, step, name, varsel, time, strand, double):
    files = sorted(glob.glob(os.path.join(indir, f'sol_{step}_r*.vtu')))
    pieces = [read_vtu(f) for f in files]
    allpts = np.concatenate([p['Points'] for p in pieces])
    xyz, inv = np.unique(allpts, axis=0, return_inverse=True)       # merge duplicate nodes of rank interfaces
    inv = inv.reshape(-1)
    offs = np.cumsum([0] + [len(p['Points']) for p in pieces])
    # ---- cell variables
    names, cols = [], []
    for key in pieces[0]:
        if key in SKIP or (varsel and not any(fnmatch.fnmatch(key, s) for s in varsel)):
            continue
        a = np.concatenate([p[key] for p in pieces])
        if a.ndim == 2:
            for c, suffix in zip(range(a.shape[1]), 'xyz'):
                names.append(f'{key}_{suffix}'); cols.append(a[:, c])
        else:
            names.append(key); cols.append(a)
    if varsel and not names:
        sys.exit('--vars matched no cell variable; available: ' + ', '.join(k for k in pieces[0] if k not in SKIP))
    # ---- classify cells
    hex_conn, hex_cells, poly_cells, poly_loops = [], [], [], []
    cell0 = 0
    for k, p in enumerate(pieces):
        gid = inv[offs[k]:offs[k + 1]]
        for ci, loops in enumerate(cell_faces(p)):
            h = hex_order(loops, p['Points'])
            if h is not None:
                hex_conn.append(gid[h]); hex_cells.append(cell0 + ci)
            else:
                poly_cells.append(cell0 + ci); poly_loops.append([gid[L] for L in loops])
        cell0 += len(p['faceoffsets'])
    zones = []
    def zone_nodes(conn):
        used, local = np.unique(conn, return_inverse=True)
        return used, local.reshape(conn.shape)
    if hex_cells:
        conn = np.array(hex_conn, dtype=np.int64)
        used, local = zone_nodes(conn)
        zones.append(dict(name='hex', ztype=5, nodes=xyz[used], conn=local, cellvars=[c[hex_cells] for c in cols],
                          time=time or 0.0, strand=1 if strand else -1))
    if poly_cells:
        extra, key2id, tets, owner = [], {}, [], []
        for ci, loops in zip(poly_cells, poly_loops):
            t = tets_of_cell(loops, xyz, extra, key2id, len(xyz))
            tets += t; owner += [ci] * len(t)
        allx = np.vstack([xyz, np.array(extra)])
        conn = np.array(tets, dtype=np.int64)
        used, local = zone_nodes(conn)
        zones.append(dict(name='poly (tetrahedra)', ztype=4, nodes=allx[used], conn=local, cellvars=[c[owner] for c in cols],
                          time=time or 0.0, strand=2 if strand else -1))
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name + '.plt')
    write_plt(path, name, ['x', 'y', 'z'] + names, zones, double)
    desc = ', '.join(f"{z['name']}: {len(z['conn'])} {'bricks' if z['ztype'] == 5 else 'tetrahedra'}/{len(z['nodes'])} nodes" for z in zones)
    print(f'{path}: {desc}; {len(names)} cell variables; {os.path.getsize(path) / 1e6:.1f} MB')


def load_config(path):
    """flat yml: PyYAML if available, otherwise a minimal parser (key: value, inline lists [a, b], '- item' lists)"""
    text = open(path).read()
    try:
        import yaml
        cfg = yaml.safe_load(text) or {}
    except ImportError:
        cfg, key = {}, None
        for raw in text.splitlines():
            line = raw.split('#', 1)[0].rstrip()
            if not line.strip():
                continue
            if line.lstrip().startswith('- ') and key is not None:
                if cfg.get(key) is None:
                    cfg[key] = []
                cfg[key].append(line.lstrip()[2:].strip().strip('"\''))
                continue
            k, _, v = line.partition(':')
            key, v = k.strip(), v.strip()
            if v == '':
                cfg[key] = None
            elif v.startswith('['):
                cfg[key] = [x.strip().strip('"\'') for x in v.strip('[]').split(',') if x.strip()]
            elif v.lower() in ('true', 'false'):
                cfg[key] = v.lower() == 'true'
            else:
                cfg[key] = v.strip('"\'')
    return cfg


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    cfgpath = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'config.yml'))
    if not os.path.exists(cfgpath):
        sys.exit(f'config file not found: {cfgpath}')
    cfg = load_config(cfgpath)
    base = os.path.dirname(cfgpath)
    path = lambda p: p if os.path.isabs(p) else os.path.normpath(os.path.join(base, p))
    for k in ('input_dir', 'output_dir'):
        if not cfg.get(k):
            sys.exit(f'{cfgpath}: "{k}" is required')
    indir, outdir = path(str(cfg['input_dir'])), path(str(cfg['output_dir']))
    step = str(cfg.get('step') or 'last')
    name = cfg.get('name') or None
    varsel = cfg.get('vars') or None
    if isinstance(varsel, str):
        varsel = [v.strip() for v in varsel.split(',') if v.strip()]
    double = bool(cfg.get('double', False))
    steps = steps_in(indir)
    if step == 'last':
        sel = [steps[-1]]
    elif step == 'all':
        sel = steps
    else:
        s = f'{int(step):07d}'
        if s not in steps:
            sys.exit(f'step {s} not found; available: {steps[0]} .. {steps[-1]} ({len(steps)} steps)')
        sel = [s]
    times = step_times(indir)
    for s in sel:
        out = (f'{name}_{s}' if len(sel) > 1 else name) if name else f'sol_{s}'
        convert(indir, outdir, s, out, varsel, times.get(s), len(sel) > 1, double)


if __name__ == '__main__':
    main()
