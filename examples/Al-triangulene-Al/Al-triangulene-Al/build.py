# Build [4]triangulene (C33H15, S = 3/2) and place it between the Al tips of
# Al-benzene-Al with its plane containing the transport axis (z).
import math, re, sys
math.dist = lambda p, q: sum((x-y)**2 for x, y in zip(p, q))**0.5
d, dCH, n = 1.42, 1.09, 4
a1 = (math.sqrt(3)*d, 0.0); a2 = (math.sqrt(3)*d/2, 1.5*d)
C = []
for i in range(n):
    for j in range(n-i):
        cx = i*a1[0]+j*a2[0]; cy = i*a1[1]+j*a2[1]
        for k in range(6):
            t = math.radians(30+60*k); p = (cx+d*math.cos(t), cy+d*math.sin(t))
            if all(math.dist(p,q) > 0.1 for q in C): C.append(p)
nb = [[j for j,q in enumerate(C) if j!=i and math.dist(p,q) < d+0.1] for i,p in enumerate(C)]
H = []
for i,p in enumerate(C):
    if len(nb[i]) == 2:
        mx = sum(C[j][0] for j in nb[i])/2; my = sum(C[j][1] for j in nb[i])/2
        vx, vy = p[0]-mx, p[1]-my; r = math.hypot(vx,vy)
        H.append((p[0]+dCH*vx/r, p[1]+dCH*vy/r))
print(len(C), "C,", len(H), "H", file=sys.stderr)
assert len(C) == 33 and len(H) == 15
# centre on the carbon centroid; the C3 median (a vertex to the opposite edge) is along the 2D y axis
gx = sum(p[0] for p in C)/33; gy = sum(p[1] for p in C)/33
mol = [("C",p[0]-gx,p[1]-gy) for p in C] + [("H",p[0]-gx,p[1]-gy) for p in H]
src = sys.argv[1]
lines = open(src).read().split('\n')
idx = [i for i,l in enumerate(lines) if re.match(r'^\s*[A-Z][a-z]?\s+-?\d+\.\d+\s+-?\d', l)]
A = [[lines[i].split()[0]]+list(map(float, lines[i].split()[1:4])) for i in idx]
top, ben, bot = A[:20], A[20:32], A[32:52]
cz = sum(a[3] for a in ben if a[0]=='C')/6
ax = (top[19][1]+bot[0][1])/2; ay = (top[19][2]+bot[0][2])/2
# molecule plane = y-z plane: 2D x -> y, 2D y (C3 median) -> z
M = [[e, ax, ay+u, cz+v] for e,u,v in mol]
def mind(apex, s):
    return min(math.dist(a[1:], apex) for a in M)
# move each tip rigidly along z until its apex is 2.50 A from the nearest molecule atom
for tip, ap, sgn in ((top, top[19], 1), (bot, bot[0], -1)):
    h = 0.0
    while mind((ap[1], ap[2], ap[3]+sgn*h), sgn) < 2.5: h += 1e-4
    for a in tip: a[3] += sgn*h
    print("tip", "top" if sgn > 0 else "bottom", "moved %+.3f A" % (sgn*h), file=sys.stderr)
print("apex gap %.3f A" % (top[19][3]-bot[0][3]), file=sys.stderr)
mc = min(math.dist(a[1:], t[1:]) for a in M if a[0]=='C' for t in (top[19], bot[0]))
print("nearest C-apex %.3f A" % mc, file=sys.stderr)
allA = top + M + bot
head = lines[:idx[0]]; tail = lines[idx[-1]+1:]
head = [("0,4" if l.strip() in ("0,1","0 1") else l) for l in head]
body = ['%-2s %12.5f %12.5f %12.5f' % tuple(a) for a in allA]
open(sys.argv[2], 'w').write('\n'.join(head+body+tail))
open(sys.argv[3], 'w').write('%d\n\n' % len(allA) + '\n'.join('%s %.4f %.4f %.4f' % tuple(a) for a in allA) + '\n')
