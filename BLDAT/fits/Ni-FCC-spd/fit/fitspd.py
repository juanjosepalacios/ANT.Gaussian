"""Bounded fit of an s+p+d Bethe-lattice set (12 fcc directions) to the CRYSTAL23 LSDA DOS
of fcc Ni (a = 3.524 A) with the CRENBS sd basis + ECP plus the scaled STO-3G 4p shell.
Usage: python3 fitspd.py <start>"""
import numpy as np, json, time, sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '../../Fe-BCC-spd/fit'))
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
        out[sp] = (blk[:, 0]*H, g*blk[:, 1]/H, g*blk[:, 2]/H, g*blk[:, 3]/H, g*blk[:, 4]/H)
    return out
C = crystal(os.path.join(here, '../crystal/s0.6/DOSS.win.DAT'))   # s, p, d, total
EG = np.arange(-13, 16.0001, 0.05); ETA = 0.1; NEL = 10.0
W = np.arange(-10, 4.0001, 0.05)
tgt = {sp: [np.interp(W, C[sp][0], C[sp][j]) for j in (1, 3, 4)] for sp in range(2)}   # s, d, total
mC = float(open(os.path.join(here, '../crystal/moment')).read())
edd1 = 1.3812081086548713            # gauge: same as the Fe sets and Ni FCC-sd
names = ['es1','des','ep1','dep','dedd','edt1','dedt','sss','sps','pps','ppp','sds','pds','pdp','dds','ddp','ddd']
lo = np.array([-5,-1, -2,-1, 0, edd1-2, 0, -4, 0.05, 0.05, -2, -3, -3, -3, -2.5, 0.05, -1.0])
hi = np.array([15, 6, 30, 6, 7, edd1+2, 7, -0.1, 5, 6, 1.5, -0.05, 3, 3, -0.1, 2.5, -0.01])
def fcc_dirs():
    r = 1/np.sqrt(2); out = []
    for a, c in [(1, 1), (1, -1)]:
        for i, j in [(0, 1), (1, 2), (0, 2)]:
            v = np.zeros(3); v[i] = a*r; v[j] = c*r; out += [v, -v]
    return np.array(out)
dirs = fcc_dirs()
R = 0.5*27.2113834
nisd = json.load(open(os.path.join(here, '../../Ni-FCC-sd/fit/best_sd.json')))      # the fit chosen for FCC-sd
nisd = dict(zip(nisd['names'], nisd['x']))
fefccspd = dict(es1=8.293, des=-0.7858, ep1=23.9446, dep=0.9025, dedd=5.16, edt1=3.0574, dedt=6.045, sss=-1.1099, sps=2.0028,
                pps=5.6321, ppp=-0.0109, sds=-0.9546, pds=-2.74, pdp=-0.0274, dds=-0.1, ddp=0.05, ddd=-0.01)   # Fe FCC-spd
pfe = {k: fefccspd[k] for k in ['ep1', 'dep', 'sps', 'pps', 'ppp', 'pds', 'pdp']}
# Papacon-old (BL028.dat, orthogonal spd), spin up shifted so that edd = edd1; spin-averaged hoppings
pap = dict(es1=edd1+(1.00243-0.48986)*R, des=(0.99207-1.00243)*R, ep1=edd1+(1.42873-0.48986)*R, dep=(1.39786-1.42873)*R,
           dedd=(0.54495-0.48986)*R, edt1=edd1+(0.49574-0.48986)*R, dedt=(0.54881-0.49574)*R,
           sss=-0.0895*R, sps=0.1366*R, pps=0.2285*R, ppp=0.0300*R, sds=-0.0405*R, pds=-0.0518*R, pdp=0.0270*R,
           dds=-0.0322*R, ddp=0.0225*R, ddd=-0.00663*R)
ppap = {k: pap[k] for k in ['ep1', 'dep', 'sps', 'pps', 'ppp', 'pds', 'pdp']}
starts = {
 'A': dict(nisd, **pfe),
 'B': pap,
 'C': dict(fefccspd, des=0.0, dep=0.0, dedd=0.6, dedt=0.6),
 'D': dict(nisd, **ppap)}
tag = sys.argv[1]
x0 = np.clip(np.array([starts[tag][n] for n in names], float), lo+1e-6, hi-1e-6)
def toP(x):
    v = dict(zip(names, x))
    P = {'es': [v['es1'], v['es1']+v['des']], 'ep': [v['ep1'], v['ep1']+v['dep']],
         'edd': [edd1, edd1+v['dedd']], 'edt': [v['edt1'], v['edt1']+v['dedt']]}
    for k in ['sss','sps','pps','ppp','sds','pds','pdp','dds','ddp','ddd']: P[k] = [v[k]]*2
    return P
def model(x):
    P = toP(x); res = {sp: b.pdos(P, sp, EG, dirs, eta=ETA)[:3] for sp in range(2)}
    tot = {sp: sum(res[sp]) for sp in range(2)}
    N = np.cumsum(tot[0]+tot[1])*0.05; EF = float(np.interp(NEL, N, EG))
    m = sum(np.trapz(np.where(EG <= EF, tot[sp]*(1-2*sp), 0), EG) for sp in range(2))
    return res, tot, EF, m
ncall = [0]; best = [1e9]
def resid(x):
    res, tot, EF, m = model(x); r = []
    for sp in range(2):
        s, p, d = res[sp]
        r += [3*(np.interp(W, EG-EF, s)-tgt[sp][0]), np.interp(W, EG-EF, d)-tgt[sp][1],
              np.interp(W, EG-EF, tot[sp])-tgt[sp][2]]
    r = np.concatenate(r+[[5*(m-mC)]]); c = float(r@r); ncall[0] += 1
    if c < best[0]:
        best[0] = c
        json.dump({'x': x.tolist(), 'names': names, 'cost': c, 'EF': EF, 'moment': m, 'edd1': edd1, 'P': toP(x)},
                  open(f'best_{tag}.json', 'w'))
    if ncall[0] % 20 == 0: print(tag, 'call', ncall[0], 'cost %.3f best %.3f moment %.3f' % (c, best[0], m), flush=True)
    return r
if __name__ == '__main__':
    t = time.time(); r = resid(x0); print(tag, 'start cost %.3f, %.1f s per call' % (r@r, time.time()-t), flush=True)
    sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=1e-3, max_nfev=300)
    print(tag, 'done', sol.status, sol.message, 'cost %.3f' % (2*sol.cost), flush=True)
    print(tag, 'active bounds:', [n for n, a in zip(names, sol.active_mask) if a != 0], flush=True)
