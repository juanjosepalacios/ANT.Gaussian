import sys, os
lam = float(sys.argv[1]); out = sys.argv[2]
s = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'base.d12')).read()
s = s.replace('CRENBS sd basis + ECP', 'CRENBS sd + STO-3G 4p (scale %.4f) + ECP' % lam)
s = s.replace('228 2\n', '228 3\n')
p = [(0.6309328384, -0.1215468600), (0.2328538976, 0.5715227604), (0.1028307363, 0.5498949471)]
blk = '0 2 3 0.0 1.0\n' + ''.join('%.10f %.10f\n' % (a*lam, c) for a, c in p)
s = s.replace('0 3 5 8.0 1.0', blk + '0 3 5 8.0 1.0')
open(out, 'w').write(s)
