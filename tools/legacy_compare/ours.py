import struct, numpy as np, sys
from scipy.spatial import ConvexHull
def unsb(path):
    d=open(path,'rb').read(); assert d[:6]==b'UNSBIN'; n,=struct.unpack('<Q',d[8:16]); o=16; R={}
    dt={0:'i1',1:'<i4',2:'<i8',3:'<f4',4:'<f8'}
    for _ in range(n):
        nl,=struct.unpack('<I',d[o:o+4]); o+=4; name=d[o:o+nl].decode(); o+=nl
        t=d[o]; o+=4; cols,=struct.unpack('<I',d[o:o+4]); o+=4; rows,=struct.unpack('<Q',d[o:o+8]); o+=8
        a=np.frombuffer(d,dtype=dt[t],count=rows*cols,offset=o).copy(); o+=rows*cols*a.itemsize; o+=(8-(rows*cols*a.itemsize)%8)%8
        R[name]=a.reshape(rows,cols) if cols>1 else a
    return R
G=unsb(sys.argv[1]+'/global.unsb'); F=unsb(sys.argv[2])
x=G['x']; e2n=G['e2n.idx'].reshape(-1,2); f2e_p,f2e_i,f2e_s=G['f2e.ptr'],G['f2e.idx'],G['f2e.sgn']; c2f_p,c2f_i=G['c2f.ptr'],G['c2f.idx']
mat=G['cell_material']; nc=len(mat); rho=F['f.U_H'].reshape(-1,5); rN=F['f.U_Na'].reshape(-1,5)
out=[]
for c in range(nc):
    if mat[c]!=0: continue
    nodes=set()
    for f in c2f_i[c2f_p[c]:c2f_p[c+1]]:
        for e in f2e_i[f2e_p[f]:f2e_p[f+1]]: nodes.update(e2n[e])
    P=x[sorted(nodes)]; vol=ConvexHull(P).volume; ctr=P.mean(0)
    U=rho[c]; r=U[0]; u=U[1:4]/r; p=(5/3-0+0)*0
    out.append([*ctr,vol,r,*u,0.6667*(U[4]-0.5*r*(u@u)),rN[c][0]])
M=np.array(out); np.save('ours_cells.npy',M)
r=np.linalg.norm(M[:,:3],axis=1); print(M.shape,'r range',r.min(),r.max(),'vol',M[:,3].sum())
print('mass H %.1f mass Na %.3f'%((M[:,3]*M[:,4]).sum(),(M[:,3]*M[:,9]).sum()),'min rhoH %.3g min pH %.3g'%(M[:,4].min(),M[:,8].min()))
