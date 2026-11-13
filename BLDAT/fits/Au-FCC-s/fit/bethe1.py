"""Single-orbital Bethe lattice with z neighbours, orthogonal: local DOS per spin."""
import numpy as np
def dos(E, es, t, z=12, eta=0.1):
    w = E + 1j*eta - es
    # directional self-energy: s = t^2/(w - (z-1)s), retarded root (Im s <= 0)
    a = (z-1); disc = np.sqrt(w*w - 4*a*t*t + 0j)
    s1 = (w - disc)/(2*a); s2 = (w + disc)/(2*a)
    s = np.where(s1.imag <= 0, s1, s2)
    G = 1.0/(w - z*s)
    return -G.imag/np.pi
