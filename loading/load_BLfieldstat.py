import numpy as np

def load_BLfieldstat(filename, N1, N2, N3, N4=11, Nload=11):
    stat = {}
    with open(filename, 'rb') as binfile:
        stat['nsamp'] = np.fromfile(binfile, dtype=np.int32, count=1)
        stat['time_interv'] = np.fromfile(binfile, dtype=np.float32, count=1)
        stat['time_incurr'] = np.fromfile(binfile, dtype=np.float32, count=1)
        dumm = np.fromfile(binfile, dtype=np.float64, count=N1*N2*N3*Nload)
        shape = (N1,N2,N3,Nload)
        dumm = dumm.reshape(shape, order='F')
        names = ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'p', 'wst']
        for i in range(Nload):
            stat[names[i]] = dumm[:,:,:,i]
    return stat

