import sys,os
l1,l2=float(sys.argv[1]),float(sys.argv[2]); out=sys.argv[3]
s=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'base.d12')).read()
s=s.replace('sd basis + ECP from the ANT Fe wire','CRENBS sd + STO-3G 4p (scales %.4f %.4f) + ECP'%(l1,l2)).replace('226 2\n','226 3\n')
p=[(0.5921156814*l1,-0.1215468600),(0.2185279254*l1,0.5715227604),(0.0965042359*l2,0.5498949471)]
s=s.replace('0 3 5 6.0 1.0','0 2 3 0.0 1.0\n'+''.join('%.10f %.10f\n'%x for x in p)+'0 3 5 6.0 1.0')
open(out,'w').write(s)
