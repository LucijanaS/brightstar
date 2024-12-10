import numpy as np
import emcee
from scipy.stats import qmc

import matplotlib.pyplot as plt

from binaries_fitfunction import seeing_double

from globi import ndim, per, gtruth, params_lo, params_hi, log2nw, nwalkers, nsteps

from reel import mjd, u_array, v_array

D = len(mjd)

def model(param):
    tph = 2*np.pi * (mjd-param[0])/per
    theta_a = param[1]
    e = 0.764
    omega = 130*np.pi/180
    I = np.pi/3
    Omega = np.pi/4
    theta_A = theta_B = 5e-9
    mdata = 0*mjd
    for k in range(D):
        mdata[k] = seeing_double(theta_A, theta_B, theta_a, e, I,
                                omega, Omega,
                                u_array[k], v_array[k], tph[k])
    return mdata

noiselev = 0.02
sigma = 0*mjd + noiselev
noise = np.random.normal(0, noiselev, D)

data = model(gtruth) + noise


print(data.shape,'data points')
total_snr = np.sum((data/sigma)**2)**(1/2)
smax = np.max(data/sigma)
print('SNR: max = %i, total = %i' % (smax,total_snr))


#data[data>1] = 1
#dcolor = [str(item) for item in data]
plt.plot(u_array,v_array,'.')
plt.gca().set_aspect(1)
plt.show()

plt.plot(u_array,data,'.')
plt.show()

def log_prob(param):
    if (param < params_lo).any() or (param > params_hi).any():
        return -1e30
    model_data = model(param)
    G = np.sum(sigma**(-2) * data * model_data)
    W = np.sum(sigma**(-2) * model_data * model_data)
    return (1/2) * (G*G/W - np.log(W))

sobol = qmc.Sobol(d=ndim,scramble=False)
qran = sobol.random_base2(m=log2nw)
p0 = 0*qran
for k in range(ndim):
    p0[:,k] = params_lo[k] + qran[:,k] * (params_hi[k]-params_lo[k])

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
sampler.run_mcmc(p0,nsteps, progress=True)
parameter_values = sampler.get_chain(discard=0)
log_prob_values = sampler.get_log_prob(discard=0)
acceptance_fraction = sampler.acceptance_fraction

for i in range(nwalkers):
    np.savetxt("chains/params_{}.txt".format(i), parameter_values[:,i,:])
    np.savetxt("chains/logprob_{}.txt".format(i), log_prob_values[:,i])
    np.savetxt("chains/accep_{}.txt".format(i), sampler.acceptance_fraction)

print("Autocorrelation time:", sampler.get_autocorr_time())

