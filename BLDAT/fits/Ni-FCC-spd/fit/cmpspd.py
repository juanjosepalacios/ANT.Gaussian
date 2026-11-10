import sys, json, numpy as np
tag = sys.argv[1] if len(sys.argv) > 1 else 'A'
sys.argv = ['fitspd.py', tag]; import fitspd as f
d = json.load(open(f'best_{tag}.json')); x = np.array(d['x'])
res, tot, EF, m = f.model(x); C = f.C; W = f.W; EG = f.EG
print('%s: cost %.3f  moment %.3f (CRYSTAL %.3f)  E_F %.3f eV' % (tag, d['cost'], m, f.mC, EF))
for j, lab in enumerate(['s', 'p', 'd', 'total']):
    row = []
    for sp in range(2):
        y = res[sp][j] if j < 3 else tot[sp]
        yb = np.interp(0, EG-EF, y); yc = np.interp(0, C[sp][0], C[sp][j+1])
        ref = np.interp(W, C[sp][0], C[sp][j+1]); rel = np.linalg.norm(np.interp(W, EG-EF, y)-ref)/np.linalg.norm(ref)
        row.append((yb, yc, rel))
    print('%-6s up %.3f / %.3f  down %.3f / %.3f   misfit %.2f / %.2f' % (lab, row[0][0], row[0][1], row[1][0], row[1][1], row[0][2], row[1][2]))
