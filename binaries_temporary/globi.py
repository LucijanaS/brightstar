import numpy as np

ndim = 2

log2nw = 3
nwalkers  = 2**log2nw
nsteps=10000

per = 29.13376
gtruth = [1121,3e-8]
#         tperi theta_a

params_lo = np.array([1000,1e-8])
params_hi = np.array([1000+per,5e-8])

np.random.seed(42)
