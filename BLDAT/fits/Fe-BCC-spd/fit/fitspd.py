"""Bounded fit of an s+p+d Bethe-lattice set (8 bcc directions) to the CRYSTAL23 LSDA DOS
of bcc Fe with the CRENBS + 4p basis. Usage: python3 fitspd.py <start>"""
import numpy as np, bethe9 as b, json, time, sys
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
C = crystal('/home/palacios/ANT-runs/Fe-bcc-bulk-spd-crystal23-a2.87-20261106/s0.586/DOSS.win.DAT')
EG = np.arange(-13, 16.0001, 0.05); dirs = b.bcc_dirs(); ETA = 0.1
W = np.arange(-10, 4.0001, 0.05)
tgt = {sp: [np.interp(W, C[sp][0], C[sp][j]) for j in (1, 3, 4)] for sp in range(2)}   # s, d, total
mC = 2.81
edd1 = 1.3812081086548713            # gauge: same as BCC-sd
names = ['es1','des','ep1','dep','dedd','edt1','dedt','sss','sps','pps','ppp','sds','pds','pdp','dds','ddp','ddd']
lo = np.array([-5,-1, -2,-1, 0, edd1-2, 0, -4, 0.05, 0.05, -2, -3, -3, -3, -2.5, 0.05, -1.0])
hi = np.array([15, 6, 30, 6, 7, edd1+2, 7, -0.1, 5, 6, 1.5, -0.05, 3, 3, -0.1, 2.5, -0.01])
bccsd = dict(es1=4.329, des=4.0, dedd=4.915, edt1=2.738, dedt=5.548, sss=-1.862, sds=-1.280, dds=-0.200, ddp=0.459, ddd=-0.294)
starts = {
 'A': dict(bccsd, ep1=12.0, dep=1.0, sps=1.8, pps=1.8, ppp=0.26, pds=-0.48, pdp=0.34),
 'B': dict(bccsd, ep1=9.0, dep=0.5, sps=1.2, pps=1.5, ppp=-0.3, pds=-0.8, pdp=0.4),
 'C': dict(bccsd, ep1=15.0, dep=1.5, sps=2.5, pps=2.5, ppp=-0.5, pds=0.8, pdp=-0.4),
 'D': dict(es1=5.858, des=0.5, ep1=13.6, dep=0.5, dedd=0.9, edt1=edd1, dedt=0.9, sss=-1.979, sps=1.81, pps=1.81,
           ppp=0.26, sds=-1.091, pds=-0.476, pdp=0.342, dds=-1.095, ddp=1.36, ddd=-0.257)}
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
    N = np.cumsum(tot[0]+tot[1])*0.05; EF = float(np.interp(8.0, N, EG))
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
t = time.time(); r = resid(x0); print(tag, 'start cost %.3f, %.1f s per call' % (r@r, time.time()-t), flush=True)
sol = least_squares(resid, x0, bounds=(lo, hi), diff_step=1e-3, max_nfev=300)
print(tag, 'done', sol.status, sol.message, 'cost %.3f' % (2*sol.cost), flush=True)
print(tag, 'active bounds:', [n for n, a in zip(names, sol.active_mask) if a != 0], flush=True)
