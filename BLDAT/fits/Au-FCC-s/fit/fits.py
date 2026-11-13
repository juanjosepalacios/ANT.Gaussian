"""Fit the single-orbital FCC-s Bethe-lattice set to the CRYSTAL23 6s DOS of fcc Au.

python3 fits.py   (needs numpy and scipy; run in this directory)
"""
import numpy as np
from scipy.optimize import minimize_scalar, brentq
from bethe1 import dos

d = np.array([list(map(float, l.split())) for l in open('../crystal/DOSS.DAT')
              if l.strip() and not l.startswith(('#', '@'))])
HA = 27.211386
E, g = d[:, 0]*HA, d[:, 1]/HA/2          # eV relative to E_F, states/eV per spin
eta, z = 0.1, 12
rhoF = np.interp(0, E, g)
w = (E >= -15) & (E <= 4)

def cost(t):
    return np.sqrt(np.sum((dos(E[w], 0, t, z, eta) - g[w])**2)/np.sum(g[w]**2))

fits = {'A DOS(E_F)': brentq(lambda t: dos(np.array([0.0]), 0, t, z, eta)[0] - rhoF, 0.5, 10),
        'B band bottom': 14.73/(2*np.sqrt(z-1)),
        'C L2 -15..4 eV': minimize_scalar(cost, bounds=(0.5, 10), method='bounded').x}
print('CRYSTAL23 per spin: DOS(E_F) = %.4f states/eV' % rhoF)
for lab, t in fits.items():
    print('%-15s |sss| = %.4f eV  band +-%.2f eV  DOS(E_F) = %.4f  misfit = %.3f'
          % (lab, t, 2*np.sqrt(z-1)*t, dos(np.array([0.0]), 0, t, z, eta)[0], cost(t)))
np.savetxt('dos_fit.dat', np.c_[E, g, dos(E, 0, fits['A DOS(E_F)'], z, eta)], fmt='%.6f',
           header='E-E_F (eV)   CRYSTAL23 DOS per spin   FCC-s Bethe DOS (eta 0.1 eV), states/eV')
