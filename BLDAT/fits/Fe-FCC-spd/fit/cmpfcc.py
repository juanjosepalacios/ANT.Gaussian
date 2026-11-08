import sys, json, numpy as np
tag = sys.argv[1] if len(sys.argv) > 1 else 'A'
sys.argv = ['x', tag]
src = open('fitfcc.py').read().split('ncall = [0]')[0]; exec(src)
d = json.load(open(f'best_{tag}.json')); x = np.array(d['x'])
res, tot, EF, m = model(x)
print('%s: cost %.3f  moment %.3f (CRYSTAL %.3f)' % (tag, d['cost'], m, mC))
lab = ['s', 'p', 'd', 'total']
for j in range(4):
    row = []
    for sp in range(2):
        y = res[sp][j] if j < 3 else tot[sp]
        yb = np.interp(0, EG-EF, y); yc = np.interp(0, C[sp][0], C[sp][j+1])
        ref = np.interp(W, C[sp][0], C[sp][j+1])
        rel = np.linalg.norm(np.interp(W, EG-EF, y)-ref)/np.linalg.norm(ref)
        row.append((yb, yc, rel))
    print('%-6s up %.3f / %.3f  down %.3f / %.3f   misfit %.2f / %.2f' % (lab[j], row[0][0], row[0][1], row[1][0], row[1][1], row[0][2], row[1][2]))
