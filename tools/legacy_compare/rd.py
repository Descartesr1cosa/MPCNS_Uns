import struct, sys, numpy as np
def read_restart(path):
    f=open(path,'rb'); assert f.read(8)==b'MPCNSRST'
    ver,step=struct.unpack('<ii',f.read(8)); t,=struct.unpack('<d',f.read(8))
    nb,nf=struct.unpack('<ii',f.read(8)); out={}
    for _ in range(nf):
        n,=struct.unpack('<i',f.read(4)); name=f.read(n).decode()
        loc,nc,ng=struct.unpack('<iii',f.read(12)); blocks=[]
        for ib in range(nb):
            lo=struct.unpack('<3i',f.read(12)); hi=struct.unpack('<3i',f.read(12)); act,=struct.unpack('<i',f.read(4))
            if act:
                shp=[h-l for l,h in zip(lo,hi)]+[nc]
                a=np.frombuffer(f.read(8*int(np.prod(shp))),dtype='<f8').reshape(shp)  # C order i,j,k,m
            else: a=None
            blocks.append((lo,hi,a))
        out[name]=(loc,nc,ng,blocks)
    return step,t,out
if __name__=='__main__':
    s,t,o=read_restart(sys.argv[1]); print(s,t)
    for k,(loc,nc,ng,b) in o.items(): print(k,loc,nc,ng,[(x[0],x[1]) for x in b][:2], [x[2] is not None for x in b][:2])
