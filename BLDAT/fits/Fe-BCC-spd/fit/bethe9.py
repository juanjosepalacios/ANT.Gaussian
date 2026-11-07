"""Python version of ANT's orthogonal Bethe-lattice Green's function for s+p+d sets
(BetheLattice.F90 InitBetheLattice/SolveDysonBL). Orbital order = ANT AOT 0..8:
s, px, py, pz, z2, xz, yz, x2-y2, xy. Parameters in eV."""
import numpy as np
RYD = 0.5*27.2113834
KEYS = ['es','ep','edd','edt','sss','sps','pps','ppp','sds','pds','pdp','dds','ddp','ddd']

def _funcs(p):
    x,y,z = p.T; r2 = x*x+y*y+z*z
    return np.stack([np.ones_like(x), x, y, z, (3*z*z-r2)/(2*np.sqrt(3)), x*z, y*z, (x*x-y*y)/2, x*y], 1)

def rotmat(d):
    d = np.asarray(d, float); d = d/np.linalg.norm(d)
    z = np.array([0, 0, 1.0]); v = np.cross(z, d); c = z @ d
    if np.linalg.norm(v) < 1e-12:
        R = np.eye(3) if c > 0 else np.diag([1, -1, -1.0])
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx/(1+c)
    p = np.random.default_rng(0).normal(size=(300, 3))
    T, *_ = np.linalg.lstsq(_funcs(p), _funcs(p @ R), rcond=None)
    return T.T

def hop_z(P, sp):
    V = np.zeros((9, 9)); g = lambda k: P[k][sp]
    V[0,0] = g('sss'); V[1,1] = V[2,2] = g('ppp'); V[3,3] = g('pps')
    V[4,4] = g('dds'); V[5,5] = V[6,6] = g('ddp'); V[7,7] = V[8,8] = g('ddd')
    V[0,3] = g('sps'); V[3,0] = -g('sps')
    V[0,4] = V[4,0] = g('sds')
    V[1,5] = g('pdp'); V[5,1] = -g('pdp'); V[2,6] = g('pdp'); V[6,2] = -g('pdp')
    V[3,4] = g('pds'); V[4,3] = -g('pds')
    return V

def H0(P, sp):
    e = lambda k: P[k][sp]
    return np.diag([e('es'), e('ep'), e('ep'), e('ep'), e('edd'), e('edt'), e('edt'), e('edd'), e('edt')])

def bcc_dirs():
    out = []
    for v in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]:
        v = np.array(v, float)/np.sqrt(3); out += [v, -v]
    return np.array(out)

def greens(P, sp, E, dirs, eta=1e-3, tol=1e-7, maxit=5000):
    Vz = hop_z(P, sp); H = H0(P, sp)
    Vk = np.array([rotmat(d) @ Vz @ rotmat(d).T for d in dirs])
    nk = len(dirs); partner = np.arange(nk) ^ 1
    z = (E + 1j*eta)[:, None, None]*np.eye(9)
    Sig = np.broadcast_to(-1j*np.eye(9), (len(E), nk, 9, 9)).copy()
    VkH = np.conj(np.transpose(Vk, (0, 2, 1)))
    for it in range(maxit):
        ST = Sig.sum(1)
        A = z[:, None] - H - ST[:, None] + Sig[:, partner]
        new = 0.5*(Vk @ np.linalg.inv(A) @ VkH + Sig)
        err = np.abs(new - Sig).mean(axis=(1, 2, 3)).max(); Sig = new
        if err < tol: break
    return np.linalg.inv(z - H - Sig.sum(1)), it

def pdos(P, sp, E, dirs, eta=1e-3):
    G, it = greens(P, sp, E, dirs, eta)
    d = -np.imag(np.diagonal(G, axis1=1, axis2=2))/np.pi
    return d[:, 0], d[:, 1:4].sum(1), d[:, 4:].sum(1), it

def load_set(path, name):
    L = open(path).read().splitlines(); i = 0
    while not (L[i].strip() == 'BEGINSET' and L[i+1].strip() == name): i += 1
    P = {k: [0.0, 0.0] for k in KEYS}; unit = 1.0; sp = 0
    for l in L[i+2:]:
        t = l.split()
        if not t or t[0].startswith('!'): continue
        if t[0] == 'ENDSET': break
        if t[0] == 'RYD': unit = RYD
        if t[0] == 'ISPIN': sp = int(t[2]) - 1
        if len(t) == 3 and t[1] == '=' and t[0] in P: P[t[0]][sp] = float(t[2])*unit
    return P
