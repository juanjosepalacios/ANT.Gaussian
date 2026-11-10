import sys, os
l1, l2 = float(sys.argv[1]), float(sys.argv[2]); out = sys.argv[3]
s = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'base.d12')).read()
s = s.replace('CRENBS sd basis + ECP', 'CRENBS sd + STO-3G 4p (scales %.4f %.4f) + ECP' % (l1, l2)).replace('228 2\n', '228 3\n')
p = [(0.6309328384*l1, -0.1215468600), (0.2328538976*l1, 0.5715227604), (0.1028307363*l2, 0.5498949471)]
s = s.replace('0 3 5 8.0 1.0', '0 2 3 0.0 1.0\n' + ''.join('%.10f %.10f\n' % x for x in p) + '0 3 5 8.0 1.0')
open(out, 'w').write(s)
