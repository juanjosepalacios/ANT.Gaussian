import sys, json, numpy as np
sys.argv = ['x', 'D']
src = open('fitspd.py').read().split('ncall = [0]')[0]; exec(src)
d = json.load(open('best_D.json')); x = np.array(d['x'])
res, tot, EF, m = model(x)
print('EF %.3f eV  moment %.3f' % (EF, m))
out={}
lab = ['s', 'p', 'd', 'total']
for j in range(4):
    for sp in range(2):
        g = 1-2*sp
        y = res[sp][j] if j < 3 else tot[sp]
        Eo=np.arange(-10,6.001,0.1); out[f'{lab[j]}{sp}']=[np.round(g*np.interp(Eo,EG-EF,y),3).tolist(),np.round(g*np.interp(Eo,C[sp][0],C[sp][j+1]),3).tolist()]
        yb = np.interp(0, EG-EF, y); yc = np.interp(0, C[sp][0], C[sp][j+1])
        mask = (W >= -10) & (W <= 4)
        rel = np.linalg.norm(np.interp(W, EG-EF, y)-np.interp(W, C[sp][0], C[sp][j+1]))/np.linalg.norm(np.interp(W, C[sp][0], C[sp][j+1]))
        print('%-5s spin %d  DOS(EF) Bethe %.3f CRYSTAL %.3f   rel. L2 misfit %.2f' % (lab[j], sp, yb, yc, rel))
json.dump(out,open('cmp_spd_D.json','w'))
