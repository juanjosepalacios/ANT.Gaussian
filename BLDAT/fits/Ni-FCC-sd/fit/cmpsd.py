import sys, json, numpy as np
tag = sys.argv[1] if len(sys.argv) > 1 else 'A'
sys.argv = ['fitsd.py', tag]; import fitsd as f
d = json.load(open(f'best_{tag}.json')); x = np.array(d['x'])
res, tot, EF, m = f.model(x)
print('%s: cost %.3f  moment %.3f (CRYSTAL %.3f)  E_F %.3f eV' % (tag, d['cost'], m, f.mC, EF))
C = f.C; W = f.W; EG = f.EG
for lab, j in [('s', 0), ('d', 2), ('total', None)]:
    row = []
    for sp in range(2):
        y = res[sp][j] if j is not None else tot[sp]; c = C[sp][1] if j == 0 else C[sp][2] if j == 2 else C[sp][3]
        yb = np.interp(0, EG-EF, y); yc = np.interp(0, C[sp][0], c)
        ref = np.interp(W, C[sp][0], c); rel = np.linalg.norm(np.interp(W, EG-EF, y)-ref)/np.linalg.norm(ref)
        row.append((yb, yc, rel))
    print('%-6s up %.3f / %.3f  down %.3f / %.3f   misfit %.2f / %.2f' % (lab, row[0][0], row[0][1], row[1][0], row[1][1], row[0][2], row[1][2]))
