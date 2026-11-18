# C60 between the Al tips of Al-benzene-Al (sp-Papacon), a hexagon facing each apex.
# Basis: LANL2DZ (with its ECP on Al) everywhere except the outermost Al plane of each tip,
# which touches the Bethe lattice and gets the minimal CRENBS-type sp basis with the
# Christiansen-Ermler ECP10CE, the basis of the Papacon Al parameters.
# Usage: python3 build.py ../../Al-benzene-Al/Al-benzene-Al/sp-Papacon/Al-benzene-Al.gjf Al-C60-Al.gjf Al-C60-Al.xyz
import math, re, sys, itertools
dist = lambda p, q: sum((x-y)**2 for x, y in zip(p, q))**0.5
phi = (1+5**0.5)/2
base = [(0, 1, 3*phi), (1, 2+phi, 2*phi), (phi, 2, 2*phi+1)]
V = set()
for b in base:
    for perm in ((0,1,2), (1,2,0), (2,0,1)):            # even permutations
        p = [b[i] for i in perm]
        for s in itertools.product((1,-1), repeat=3):
            V.add(tuple(round(s[k]*p[k], 9) for k in range(3)))
V = [list(v) for v in V]
assert len(V) == 60, len(V)
edge = 1.44; sc = edge/2.0
V = [[c*sc for c in v] for v in V]
nb = [[j for j in range(60) if j != i and abs(dist(V[i], V[j])-edge) < 0.01] for i in range(60)]
assert all(len(n) == 3 for n in nb)
# find a hexagonal face: 6-cycle; take the six atoms nearest to the direction of a hexagon centre
# hexagon centres lie along C3 axes = directions of icosahedron face centres, e.g. (1,1,1)
n = [1/3**0.5]*3
top = sorted(range(60), key=lambda i: -sum(V[i][k]*n[k] for k in range(3)))[:6]
c = [sum(V[i][k] for i in top)/6 for k in range(3)]
r = sum(c[k]**2 for k in range(3))**0.5; n = [c[k]/r for k in range(3)]
# rotation taking n to z (Rodrigues)
z = (0, 0, 1); ax_ = (n[1]*z[2]-n[2]*z[1], n[2]*z[0]-n[0]*z[2], n[0]*z[1]-n[1]*z[0])
s = sum(a*a for a in ax_)**0.5; cth = n[2]; ax_ = [a/s for a in ax_]
def rot(v):
    kx, ky, kz = ax_; th = math.acos(cth); ct, st = math.cos(th), math.sin(th)
    dot = kx*v[0]+ky*v[1]+kz*v[2]
    cr = (ky*v[2]-kz*v[1], kz*v[0]-kx*v[2], kx*v[1]-ky*v[0])
    return [v[i]*ct + cr[i]*st + [kx, ky, kz][i]*dot*(1-ct) for i in range(3)]
V = [rot(v) for v in V]
zs = sorted(v[2] for v in V); print("hexagon planes at z = %.3f and %.3f" % (zs[0], zs[-1]), file=sys.stderr)
assert abs(zs[-1]-zs[-6]) < 1e-6 and abs(zs[0]-zs[5]) < 1e-6
src = sys.argv[1]
lines = open(src).read().split('\n')
idx = [i for i, l in enumerate(lines) if re.match(r'^\s*[A-Z][a-z]?\s+-?\d+\.\d+\s+-?\d', l)]
A = [[lines[i].split()[0]]+list(map(float, lines[i].split()[1:4])) for i in idx]
toptip, ben, bot = A[:20], A[20:32], A[32:52]
cz = sum(a[3] for a in ben if a[0] == 'C')/6
ax = (toptip[19][1]+bot[0][1])/2; ay = (toptip[19][2]+bot[0][2])/2
M = [['C', ax+v[0], ay+v[1], cz+v[2]] for v in V]
for tip, api, sgn in ((toptip, 19, 1), (bot, 0, -1)):
    ap = tip[api]; h0 = (max if sgn > 0 else min)(a[3] for a in M) - ap[3]
    for a in tip: a[3] += h0
    h = 0.0
    while min(dist(a[1:], (tip[api][1], tip[api][2], tip[api][3]+sgn*h)) for a in M) < 2.5: h += 1e-4
    for a in tip: a[3] += sgn*h
print("apex gap %.3f A, apex to hexagon plane %.3f A" % (toptip[19][3]-bot[0][3], toptip[19][3]-max(a[3] for a in M)), file=sys.stderr)
allA = toptip + M + bot
assert len(allA) == 100
contact = list(range(1, 11)) + list(range(91, 101))
inner_al = list(range(11, 21)) + list(range(81, 91))
crenbs_basis = """S    4  1.0
     9.351    -0.002138
     0.9572   -0.226407
     0.1586    0.752703
     0.05381   0.393166
P    4  1.0
     3.906    -0.008570
     1.576    -0.015353
     0.1977    0.485440
     0.05730   0.622006
****"""
crenbs_ecp = """ECP10CE 2 10
L=2 COMPONENT
  4
    2       1.221100        -.537981
    2       3.368100       -5.459756
    2       9.750000      -16.655343
    1      29.269300       -6.475215
L=0 COMPONENT
  5
    2       1.563100      -56.205213
    2       1.771200      149.689955
    2       2.062300      -91.454394
    1       3.358300        3.728949
    0       2.130000        3.037994
L=1 COMPONENT
  5
    2       1.823100       93.675606
    2       2.124900     -189.888968
    2       2.570500      110.248104
    1       1.757500        4.199596
    0       6.769300        5.003356"""
L = lambda xs: ' '.join(map(str, xs)) + ' 0'
out = ["%mem=8GB", "%NProcShared=12", "%Subst L101 ../../../../bin", "%Subst L502 ../../../../bin",
       "#p ublyp/gen pseudo=read guess=mix scf=(maxcycle=300,conver=5,dsymm)", "", "Al-C60-Al", "", "0,1"]
out += ['%-2s %12.5f %12.5f %12.5f' % tuple(a) for a in allA]
out += ["", L(contact[:10]), crenbs_basis, L(contact[10:]), crenbs_basis,
        L(inner_al[:10]), "LANL2DZ", "****", L(inner_al[10:]), "LANL2DZ", "****", "C 0", "LANL2DZ", "****", "",
        L(contact[:10]), crenbs_ecp, L(contact[10:]), crenbs_ecp,
        L(inner_al[:10]), "LANL2DZ", L(inner_al[10:]), "LANL2DZ", "", ""]
open(sys.argv[2], 'w').write('\n'.join(out))
open(sys.argv[3], 'w').write('100\n\n' + '\n'.join('%s %.4f %.4f %.4f' % tuple(a) for a in allA) + '\n')
