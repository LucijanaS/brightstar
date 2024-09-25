import numpy as np
from fourier import fourier

SI_c = 299792458
SI_h = 6.626070040e-34
SI_k = 1.38064852e-23

# Generates amplitude map from T map

def sbright(T,lam,ds):
    c = SI_c
    h = SI_h
    k = SI_k
    z = h*c/(lam*k*(T+2.725))  # to avoid overflow in exp(z)
    z = np.clip(z,0,20)
    return (ds/lam)**2 / (np.exp(z)-1)

def correldens(S,lam):
    N = S.shape[0]
    c = SI_c
    h = SI_h
    pcoh = np.sum(S)
    print('Phi = %9.2e' % pcoh)
    pflux = 2*c/lam * pcoh
    mag = -2.5*(np.log10(h*pflux)+22.44)
    print('AB = %5.2f' % mag)
    V = fourier(S)
    Vmax = V[N//2,N//2]
    V2 = abs(V/Vmax)**2
    f = pcoh * V2
    return f

