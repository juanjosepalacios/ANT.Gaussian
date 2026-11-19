# Replace the Al tips of an Al-molecule-Al input by Au tips of the same shape.
# Usage: python3 build.py ../Al-triangulene-Al/sp-Papacon-UHF/Al-triangulene-Al.gjf Au-triangulene-Au Au-triangulene-Au.gjf Au-triangulene-Au.xyz
# Each 20-atom Al pyramid (fcc, nearest neighbours 2.496 A) is scaled about its apex to Au
# (nearest neighbours 2.884 A, a = 4.078 A) and moved rigidly along z until the nearest atom of
# the molecule is 2.50 A from its apex, as in the Al inputs. The molecule is unchanged.
# Basis: the outermost plane of each tip (coupled to the Bethe lattice) gets the CRENBS minimal
# spd basis of Au with its 68-electron ECP, the basis of the 9-orbital Papacon Au set (as in
# examples/benzene/input_Electrodo_Au_3_benzene); the other Au atoms and C and H get LANL2DZ.
import re, sys, os
dist = lambda p, q: sum((x-y)**2 for x, y in zip(p, q))**0.5
src, name, outg, outx = sys.argv[1:5]
lines = open(src).read().split('\n')
head = lines[:lines.index(next(l for l in lines if re.match(r'^-?\d+,\d+$', l.strip())))+1]
route = next(l for l in head if l.startswith('#'))
chm = head[-1].strip()
A = [[l.split()[0]]+list(map(float, l.split()[1:4])) for l in lines
     if re.match(r'^\s*[A-Z][a-z]?\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*$', l)]
N = len(A); tip1, mol, tip2 = A[:20], A[20:N-20], A[N-20:]
assert all(a[0] == 'Al' for a in tip1+tip2) and not any(a[0] == 'Al' for a in mol)
s = 2.884/2.496
for tip, api, sgn in ((tip1, 19, 1), (tip2, 0, -1)):
    ap = list(tip[api][1:])
    for a in tip:
        a[0] = 'Au'
        for k in range(3): a[k+1] = ap[k] + s*(a[k+1]-ap[k])
    h0 = (max if sgn > 0 else min)(a[3] for a in mol) - tip[api][3]
    for a in tip: a[3] += h0
    h = 0.0
    while min(dist(a[1:], (tip[api][1], tip[api][2], tip[api][3]+sgn*h)) for a in mol) < 2.5: h += 1e-4
    for a in tip: a[3] += sgn*h
gap = tip1[19][3]-tip2[0][3]
print("%s: %d atoms, apex gap %.3f A" % (name, N, gap), file=sys.stderr)
allA = tip1 + mol + tip2
contact = list(range(1, 11)) + list(range(N-9, N+1))
inner = list(range(11, 21)) + list(range(N-19, N-9))
L = lambda xs: ' '.join(map(str, xs)) + ' 0'
bas = """S   3   1.00
      0.4409000             -0.7478000
      0.2626000              0.8176000
      0.0617000              0.8113000
P   3   1.00
      0.6642000             -0.0459000
      0.1073000              0.3778000
      0.0319000              0.7233000
D   4   1.00
      1.3517000              0.4787000
      0.5149000              0.4582000
      0.1994000              0.2203000
      0.0776000              0.0412000
****"""
ecp = """AU-ECP     4     68
G POTENTIAL
  6
2      0.96569997            -0.98690200
2      2.29040003            -9.68630600
2      6.43839979           -32.57246399
2     16.36459923          -118.53087616
2     55.98929977          -286.25463867
1    171.10139465           -51.27423096
S-G POTENTIAL
  8
2      1.22479999           -84.96929932
2      1.43359995           280.17764282
2      1.89289999          -354.64355469
2      2.69630003           334.33190918
2      3.98990011          -219.39260864
2      6.00549984           170.85874939
1     14.58259964            56.33872986
0     39.31069946             6.40929985
P-G POTENTIAL
  8
2      0.90640002            -6.91190100
2      1.03110003            72.12546539
2      1.32869995          -114.97533417
2      1.90369999           123.61668396
2      2.83640003           -71.13571167
2      4.58839989            79.01660156
1     11.19810009            51.28700638
0     32.29570007             5.65198278
D-G POTENTIAL
  8
2      0.29310000            -0.11365700
2      1.03009999            -2.40402508
2      2.01859999            62.58111572
2      3.43919992           255.44717407
2      5.15789986          -161.80647278
2      2.51489997          -175.77076721
1      8.06239986            63.62079239
0     41.02849960             6.75031996
F-G POTENTIAL
  8
2      0.43349999           -47.60779953
2      0.50970000           176.31802368
2      0.57959998          -168.67689514
2      0.77109998            54.49244690
2      2.59349990            37.71539307
2      8.39070034           -24.13967514
1      6.92939997            49.16802216
0     22.25889969             4.75958586"""
els = sorted(set(a[0] for a in mol))
out = ["%mem=8GB", "%NProcShared=12", "%Subst L101 ../../../../bin", "%Subst L502 ../../../../bin",
       route, "", name, "", chm]
out += ['%-2s %12.5f %12.5f %12.5f' % tuple(a) for a in allA]
out += ["", L(contact), bas, L(inner), "LANL2DZ", "****", ' '.join(els) + " 0", "LANL2DZ", "****", "",
        L(contact), ecp, L(inner), "LANL2DZ", "", ""]
open(outg, 'w').write('\n'.join(out))
open(outx, 'w').write('%d\n%s, apex gap %.3f A\n' % (N, name, gap) + '\n'.join('%s %.5f %.5f %.5f' % tuple(a) for a in allA) + '\n')
print("MOLBLOCKS 21,%d" % (N-20), file=sys.stderr)
