#!/usr/bin/env python3
"""Debug aid: print the extreme fluid cells (min p_H, max |J|, max |u|, hottest Na, max rho_Na, ...) of a MercuryUns dump.
usage: inspect_vtu_extremes.py <step as 7 digits, e.g. 0000060> <DATA directory>
(reads sol_<step>_r*.vtu with the reader of tools/check_vtu.py; run the solver with output_step=1 output_start=N to get the dumps)
"""
import sys, glob, numpy as np
import os
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'check_vtu.py')).read().split("def main")[0]
exec(src)
step=sys.argv[1]; d=sys.argv[2]
rows=[]
for f in sorted(glob.glob(f'{d}/sol_{step}_r*.vtu')):
    a=read_vtu(f); pts=a['Points']; faces=a['faces']; foffs=a['faceoffsets']
    start=0
    for ci,end in enumerate(foffs):
        k=start; nf=faces[k]; k+=1; ids=set()
        for _ in range(nf):
            m=faces[k]; ids.update(faces[k+1:k+1+m]); k+=1+m
        start=end
        c=pts[list(ids)].mean(axis=0)
        rows.append(dict(x=c,r=np.linalg.norm(c),mat=a['material'][ci],rH=a['rho_H'][ci],pH=a['p_H'][ci],rN=a['rho_Na'][ci],pN=a['p_Na'][ci],uH=np.linalg.norm(a['u_H'][ci]),uN=np.linalg.norm(a['u_Na'][ci]),B=np.linalg.norm(a['B_total'][ci]),Bi=np.linalg.norm(a['B_induced'][ci]),J=np.linalg.norm(a['J'][ci])))
fl=[r for r in rows if r['mat']==0]
def show(title,key,rev=False,n=4):
    print(title)
    for r in sorted(fl,key=key,reverse=rev)[:n]:
        print("  r=%.3f x=%s rhoH=%.3g pH=%.3g rhoNa=%.3g pNa=%.3g |uH|=%.3g |uNa|=%.3g |B|=%.3g |Bind|=%.3g |J|=%.3g"%(r['r'],np.round(r['x'],2),r['rH'],r['pH'],r['rN'],r['pN'],r['uH'],r['uN'],r['B'],r['Bi'],r['J']))
show("min pH",lambda r:r['pH'])
show("max |J|",lambda r:-r['J'])
show("max |uH|",lambda r:-r['uH'])
show("max |uNa|",lambda r:-r['uN'])
show("max Na temperature p/rho",lambda r:-(r['pN']/max(r['rN'],1e-30)))
show("max H temperature p/rho",lambda r:-(r['pH']/max(r['rH'],1e-30)))
show("min rhoNa",lambda r:r['rN'])
show("max rhoNa",lambda r:-r['rN'],n=5)
show("max pNa",lambda r:-r['pN'],n=3)
