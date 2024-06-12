import numpy as np
import matplotlib.pyplot as pl

from fourier import grids
from sources import blob, smooth, crescent, eclipsed_sphere
from visib import sbright, correldens
from graphics import draw

pl.style.use('dark_background')

N = 2048  # N is grid size (equal on ground and sky)
ds = 1e-10  # ds is grid spacing on the sky

lam = 500e-9 # wavelegnth in meters

sx, sy, x, y = grids(ds, N, lam)

Teff = 1e4
Tmap = Teff * eclipsed_sphere(sx, sy, 3e-9, 1/10)
Tmap = smooth(Tmap, 3)
draw(sx, sy, Tmap, 16, 'sky', title='$T_{\\rm eff}$')
pl.show()

S = sbright(Tmap, lam, ds)
f = correldens(S, lam)

draw(x, y, f, 16, 'ground', cmap='coolwarm', title='$\Phi\,|V|^2$')
pl.show()
