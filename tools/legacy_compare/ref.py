import struct, glob, numpy as np, sys
from rd import read_restart
def read_grid(p):
    d=open(p,'rb').read(); nb,=struct.unpack('<i',d[:4]); o=4; names=[]
    for _ in range(nb):
        n,=struct.unpack('<i',d[o:o+4]); o+=4; names.append(d[o:o+n].rstrip(b'\0').decode()); o+=n
    dims=[]
    for _ in range(nb): dims.append(struct.unpack('<3i',d[o:o+12])); o+=12
    blocks=[]
    for nm,(ni,nj,nk) in zip(names,dims):
        n=ni*nj*nk; a=np.frombuffer(d,dtype='<f8',count=3*n,offset=o).reshape(3,n); o+=24*n
        blocks.append((nm,(ni,nj,nk),a))
    return blocks
D=glob.glob('../DATA/Step*')[0]
cells=[]  # rows: x,y,z,vol,rhoH,uHx,uHy,uHz,pH,rhoNa
for r in range(56):
    step,t,F=read_restart(f'{D}/flow_field{r:04d}.bin')
    G=read_grid(f'geometry/grid  {r:2d}.grd'.replace('  ','  ',1)) if False else read_grid(glob.glob(f'geometry/grid *{r}.grd')[[g.split()[-1]=='%d.grd'%r for g in glob.glob(f'geometry/grid *{r}.grd')].index(True)])
    # blocks of this rank = those active in U_H (fluid) in file order
    UH=F['U_H'][3]; UN=F['U_Na'][3]
    gi=0
    act=[ib for ib,(lo,hi,a) in enumerate(F['B_xi'][3]) if a is not None]
    for gb,ib in zip(G,act):
        nm,(ni,nj,nk),a=gb
        if nm!='Fluid': continue
        lo,hi,uh=UH[ib]; _,_,un=UN[ib]
        X=a.reshape(3,ni,nj,nk) if False else a.reshape(3,nk,nj,ni).transpose(0,3,2,1)  # try Fortran order (i fastest)
        c=0.125*sum(X[:, i:ni-1+i if i==0 else ni, :, :][:, :, :, :] for i in [])  if False else None
        def cc(X):
            s=0
            for di in (0,1):
                for dj in (0,1):
                    for dk in (0,1): s=s+X[:,di:ni-1+di,dj:nj-1+dj,dk:nk-1+dk]
            return s/8
        C=cc(X)  # (3,ni-1,nj-1,nk-1)
        n=(ni-1,nj-1,nk-1)
        U=uh[4:4+n[0],4:4+n[1],4:4+n[2],:]; V=un[4:4+n[0],4:4+n[1],4:4+n[2],:]
        rho=U[...,0]; u=U[...,1:4]/rho[...,None]; p=0.6667*(U[...,4]-0.5*rho*(u**2).sum(-1))
        # volume ~ product of mean edge vectors
        e1=X[:,1:,:-1,:-1]-X[:,:-1,:-1,:-1]; e2=X[:,:-1,1:,:-1]-X[:,:-1,:-1,:-1]; e3=X[:,:-1,:-1,1:]-X[:,:-1,:-1,:-1]
        vol=np.abs(np.einsum('cijk,cijk->ijk',np.cross(e1,e2,axis=0),e3))
        m=np.stack([C[0],C[1],C[2],vol,rho,u[...,0],u[...,1],u[...,2],p,V[...,0]],-1).reshape(-1,10)
        cells.append(m)
M=np.concatenate(cells); np.save('ref_cells.npy',M); print(M.shape, 'time',t)
r=np.hypot(np.hypot(M[:,0],M[:,1]),M[:,2]); print('r range',r.min(),r.max(), 'sum vol',M[:,3].sum())
print('mass H %.1f mass Na %.3f'%((M[:,3]*M[:,4]).sum(),(M[:,3]*M[:,9]).sum()),'min rhoH %.3g min pH %.3g'%(M[:,4].min(),M[:,8].min()))
