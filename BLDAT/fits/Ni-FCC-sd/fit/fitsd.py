"""Bounded fit of an s+d Bethe-lattice set (12 fcc directions) to the CRYSTAL23 LSDA DOS
of fcc Ni (a = 3.524 A) with the CRENBS sd basis + ECP. Usage: python3 fitsd.py <start>
The s+p+d model of ../../Fe-BCC-spd/fit/bethe9.py is used with the p shell decoupled
(no p hoppings, p on-site energy far above the window), which gives the s+d Green's function."""
import numpy as np, json, time, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../Fe-BCC-spd/fit'))
import bethe9 as b
from scipy.optimize import least_squares
H = 27.211386
def crystal(path):
    B = [[]]
    for l in open(path):
        if l.strip().startswith('&'): B.append([]); continue
        if not l.strip() or l.lstrip()[0] in '#@': continue
        B[-1].append([float(x) for x in l.split()])
    out = {}
    for sp, blk in enumerate([np.array(x) for x in B if x][:2]):
        g = 1 if sp == 0 else -1
        out[sp] = (blk[:, 0]*H,) + tuple(g*blk[:, j]/H for j in range(1, blk.shape[1]))
    return out
C = crystal(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../crystal/DOSS.win.DAT'))   # s, d, total
EG = np.arange(-13, 16.0001, 0.05); ETA = 0.1; NEL = 10.0
W = np.arange(-10, 4.0001, 0.05)
tgt = {sp: [np.interp(W, C[sp][0], C[sp][j]) for j in (1, 2)] for sp in range(2)}   # s, d
mC = float(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../crystal/moment')).read())
edd1 = 1.3812081086548713            # gauge: same as the Fe sets
names = ['es1','des','dedd','edt1','dedt','sss','sds','dds','ddp','ddd']
lo = np.array([-5, -1, 0, edd1-2, 0, -3, -3, -2.5, 0.05, -1.0])
hi = np.array([15, 4, 6, edd1+2, 6, -0.2, -0.1, -0.2, 2.0, -0.01])
def fcc_dirs():
    r = 1/np.sqrt(2); out = []
    for a, c in [(1, 1), (1, -1)]:
        for i, j in [(0, 1), (1, 2), (0, 2)]:
            v = np.zeros(3); v[i] = a*r; v[j] = c*r; out += [v, -v]
    return np.array(out)
dirs = fcc_dirs()
R = 0.5*27.2113834
# Fe FCC-sd (BL026.dat), with a Ni-like exchange splitting
fefccsd = dict(es1=2.8011, des=0.2, dedd=0.6, edt1=2.0911, dedt=0.6, sss=-1.1748, sds=-1.0804, dds=-0.2, ddp=0.677, ddd=-0.0101)
# Papacon-old (BL028.dat, orthogonal spd), spin up shifted so that edd = edd1; spin splittings and hoppings as there
pap = dict(es1=edd1+(1.00243-0.48986)*R, des=(0.99207-1.00243)*R,
           dedd=(0.54495-0.48986)*R, edt1=edd1+(0.49574-0.48986)*R, dedt=(0.54881-0.49574)*R,
           sss=-0.0895*R, sds=-0.0405*R, dds=-0.0322*R, ddp=0.0225*R, ddd=-0.00663*R)
starts = {
 'A': fefccsd,
 'B': pap,
 'C': dict(pap, sss=-1.1748, sds=-1.0804),
 'D': dict(fefccsd, dds=pap['dds'], ddp=pap['ddp'], ddd=pap['ddd'])}
tag = sys.argv[1]
x0 = np.clip(np.array([starts[tag][n] for n in names], float), lo+1e-6, hi-1e-6)
def toP(x):
    v = dict(zip(names, x))
    P = {'es': [v['es1'], v['es1']+v['des']], 'ep': [200.0, 200.0],
         'edd': [edd1, edd1+v['dedd']], 'edt': [v['edt1'], v['edt1']+v['dedt']]}
    for k in ['sss','sds','dds','ddp','ddd']: P[k] = [v[k]]*2
    for k in ['sps','pps','ppp','pds','pdp']: P[k] = [0.0]*2
    return P
def model(x):
    P = toP(x); res = {sp: b.pdos(P, sp, EG, dirs, eta=ETA)[:3] for sp in range(2)}
    tot = {sp: res[sp][0]+res[sp][2] for sp in range(2)}
    N = np.cumsum(tot[0]+tot[1])*0.05; EF = float(np.interp(NEL, N, EG))
    m = sum(np.trapz(np.where(EG <= EF, tot[sp]*(1-2*sp), 0), EG) for sp in range(2))
    return res, tot, EF, m
ncall = [0]; best = [1e9]
def resid(x):
    res, tot, EF, m = model(x); r = []
    for sp in range(2):
        s, p, d = res[sp]
        r += [3*(np.interp(W, EG-EF, s)-tgt[sp][0]), np.interp(W, EG-EF, d)-tgt[sp][1]]
    r = np.concatenate(r+[[5*(m-mC)]]); c = float(r@r); ncall[0] += 1
    if c < best[0]:
        best[0] = c
        json.dump({'x': x.tolist(), 'names': names, 'cost': c, 'EF': EF, 'moment': m, 'edd1': edd1, 'P': toP(x)},
                  open(f'best_{tag}.json', 'w'))
    if ncall[0] % 20 == 0: print(tag, 'call', ncall[0], 'cost %.3f best %.3f moment %.3f' % (c, best[0], m), flush=True)
    return r
if __name__ == '__main__':
    t = time.time(); r = resid(x0); print(tag, 'start cost %.3f, %.1f s per call' % (r@r, time.time()-t), flush=True)
    sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=1e-3, max_nfev=400)
    print(tag, 'done', sol.status, sol.message, 'cost %.3f' % (2*sol.cost), flush=True)
    print(tag, 'active bounds:', [n for n, a in zip(names, sol.active_mask) if a != 0], flush=True)
