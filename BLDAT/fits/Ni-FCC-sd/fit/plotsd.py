"""python3 plotsd.py C D  (numpy/scipy) writes dos_fit.npz; then a python with matplotlib: plotsd.py --plot"""
import sys, json, numpy as np
if sys.argv[1] != '--plot':
    tags = sys.argv[1:]; sys.argv = ['fitsd.py', tags[0]]; import fitsd as f
    out = {'tags': np.array(tags)}
    for sp in range(2):
        out[f'ref{sp}'] = np.array(f.C[sp][:4])
    for t in tags:
        x = np.array(json.load(open(f'best_{t}.json'))['x']); res, tot, EF, m = f.model(x)
        for sp in range(2): out[f'{t}{sp}'] = np.array([f.EG-EF, res[sp][0], res[sp][2], tot[sp]])
    np.savez('dos_fit.npz', **out); sys.exit()
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = np.load('dos_fit.npz'); fig, ax = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
for sp in range(2):
    for t in D['tags']:
        for i in range(3): ax[i, sp].plot(D[f'{t}{sp}'][0], D[f'{t}{sp}'][i+1], label=f'Bethe {t}')
    for i in range(3): ax[i, sp].plot(D[f'ref{sp}'][0], D[f'ref{sp}'][i+1], 'k', lw=0.8, label='CRYSTAL23')
    ax[0, sp].set_title('spin ' + ['up', 'down'][sp]); ax[2, sp].set_xlim(-10, 4); ax[2, sp].set_xlabel('E - E_F (eV)')
ax[0, 0].legend(); [a.set_ylabel(l + ' DOS (states/eV)') for a, l in zip(ax[:, 0], ['s', 'd', 'total'])]
plt.tight_layout(); plt.savefig('dos_fit.png', dpi=90)
