import numpy as np
import matplotlib.pyplot as pl
import corner

#pl.style.use('dark_background')

from globi import ndim, params_lo, params_hi, nwalkers

params = []
logprob = []

for c in range(nwalkers):
    params_c = np.loadtxt("chains/params_{}.txt".format(c))
    logprob_c = np.loadtxt("chains/logprob_{}.txt".format(c))
    if c in list(range(nwalkers)):
    #if c in [0,1,2,4,5,6,7]:
        params.append(params_c)
        logprob.append(logprob_c)
    
for p in range(ndim):
    pl.subplot(ndim+1,1,p+1)
    for c in range(len(params)):
        pl.plot(params[c][:,p])

pl.subplot(ndim+1,1,ndim+1)
for c in range(len(logprob)):
    pl.plot(logprob[c])

for c in range(len(params)):
    bi = params[c].shape[0]//2
    params[c] = params[c][bi:,:]

prange = []

for p in range(params[0].shape[1]):
    prange.append((params_lo[p],params_hi[p]))

samp = np.concatenate(params)
figure = corner.corner(samp,labels=['t_peri','theta_a','e'],
                       range=prange)


axes = np.array(figure.axes).reshape((ndim, ndim))


pl.show()
    
