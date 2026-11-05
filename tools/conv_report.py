#!/usr/bin/env python3
"""Convergence table of MercuryUns pseudo-time (steady) runs.
usage: conv_report.py [--tol_H 1e-5 --tol_B 1e-4] run.log [run.log ...]
Parses the 'step ... res_H= res_Na= res_B=' lines and the final 'real' time (from `time mpirun ...`, if present).
For every log: first logged step with res_H < tol_H and res_B < tol_B, residuals at the end, wall time, wall time per step, and the
(interpolated) wall time to reach the tolerance."""
import re, sys, argparse
ap = argparse.ArgumentParser()
ap.add_argument('--tol_H', type=float, default=1e-5)
ap.add_argument('--tol_B', type=float, default=1e-4)
ap.add_argument('logs', nargs='+')
a = ap.parse_args()
def num(s, key):
    m = re.search(key + r'=([-+0-9.eE]+)', s)
    try:
        return float(m.group(1)) if m else float("nan")
    except ValueError:
        return float("nan")
print("%-24s %7s %9s %9s %9s %8s %8s | %8s %9s" % ("run", "steps", "res_H", "res_Na", "res_B", "wall[s]", "ms/step", "it->tol", "wall->tol"))
for fn in a.logs:
    txt = open(fn).read()
    rows = []
    for line in txt.splitlines():
        if line.startswith('step'):
            st = int(line.split()[1])
            rows.append((st, num(line, 'res_H'), num(line, 'res_Na'), num(line, 'res_B')))
    m = re.search(r'real\s+(\d+)m([\d.]+)s', txt)
    wall = int(m.group(1)) * 60 + float(m.group(2)) if m else float('nan')
    if not rows:
        print(fn, "no data"); continue
    last = rows[-1]
    hit = next((r for r in rows if r[0] > 0 and r[1] < a.tol_H and r[3] < a.tol_B), None)
    per = wall / last[0] * 1000 if last[0] else float('nan')
    print("%-24s %7d %9.2e %9.2e %9.2e %8.1f %8.1f | %8s %9s" % (fn.split('/')[-2] if '/' in fn else fn, last[0], last[1], last[2], last[3], wall, per,
          hit[0] if hit else '-', ("%.1f" % (hit[0] * per / 1000)) if hit else '-'))
