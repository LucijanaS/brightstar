import numpy as np
import matplotlib.pyplot as plt

def get_stuff(fname):
    lam = 425e-9
    all = np.loadtxt(fname,skiprows=1)
    mjd = all[:,0]
    u,v = all[:,2]/lam, all[:,3]/lam
    return mjd,u,v

mjd, u_array,v_array = get_stuff('binary_4_extract.txt')

if __name__ == '__main__':
    plt.plot(mjd,'.')
    plt.show()
    plt.plot(u,v,'.')
    plt.gca().set_aspect(1)
    plt.show()