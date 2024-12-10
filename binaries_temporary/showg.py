import numpy as np
import matplotlib.pyplot as pl

from binaries_fitfunction import seeing_double

scl = 5e-9  # about 1 mas
B = 1/scl

N = 100
u = np.linspace(-B,B,N)
v = np.linspace(-B,B,N)
u,v = np.meshgrid(u,v)
V2 = 0*u

theta_A = theta_B = scl
theta_a = 5*scl
e = 0.7
I = np.pi/3
omega = Omega = np.pi/4


for tph in np.linspace(0,2*np.pi,100):
    pl.clf()
    for i in range(N):
        for j in range(N):
            U = u[i,j]
            V = v[i,j]
            V2[i,j] = seeing_double(theta_A, theta_B, theta_a, e, I,
                                    omega, Omega, U, V, tph)
    pl.gca().set_aspect(1)
    pl.contourf(u,v,V2)
    pl.pause(0.01)
