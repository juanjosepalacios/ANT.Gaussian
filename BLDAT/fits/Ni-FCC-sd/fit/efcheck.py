"""Fermi-level check: occupations below E_F (Bethe vs CRYSTAL23) and the rigid shift of the Bethe
DOS that best matches the CRYSTAL23 d and total DOS (0 = Fermi levels consistent)."""
import sys, json, numpy as np
tag = sys.argv[1]; sys.argv = ['fitsd.py', tag]; import fitsd as f
x = np.array(json.load(open(f'best_{tag}.json'))['x']); res, tot, EF, m = f.model(x)
E = f.EG - EF; C = f.C
for sp in range(2):
    ob = [np.trapz(np.where(E <= 0, y, 0), E) for y in (res[sp][0], res[sp][2])]
    Ec = C[sp][0]; oc = [np.trapz(np.where(Ec <= 0, C[sp][j], 0), Ec) for j in (1, 2)]
    print('spin %d  N_s %.3f / %.3f  N_d %.3f / %.3f' % (sp+1, ob[0], oc[0], ob[1], oc[1]))
sh = np.arange(-0.5, 0.5001, 0.01)
def mis(s):
    return sum(np.sum((np.interp(f.W, E+s, y)-np.interp(f.W, C[sp][0], C[sp][j]))**2)
               for sp in range(2) for y, j in ((res[sp][2], 2), (tot[sp], 3)))
v = [mis(s) for s in sh]; print('best rigid shift %.2f eV (misfit %.2f at 0, %.2f at best)' % (sh[np.argmin(v)], mis(0), min(v)))
