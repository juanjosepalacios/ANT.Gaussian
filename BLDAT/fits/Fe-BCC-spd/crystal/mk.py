import sys
lam=float(sys.argv[1]); out=sys.argv[2]
s=open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'base.d12')).read()
s=s.replace('sd basis + ECP from the ANT Fe wire','CRENBS sd + STO-3G 4p (scale %.4f) + ECP'%lam)
s=s.replace('226 2\n','226 3\n')
p=[(0.5921156814,-0.1215468600),(0.2185279254,0.5715227604),(0.0965042359,0.5498949471)]
blk='0 2 3 0.0 1.0\n'+''.join('%.10f %.10f\n'%(a*lam,c) for a,c in p)
s=s.replace('0 3 5 6.0 1.0',blk+'0 3 5 6.0 1.0')
open(out,'w').write(s)
