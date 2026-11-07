import numpy as np, json
H=27.211386
def rd(p):
    B=[[]]
    for l in open(p):
        if l.strip().startswith('&'): B.append([]); continue
        if not l.strip() or l.lstrip()[0] in '#@': continue
        B[-1].append([float(x) for x in l.split()])
    return [np.array(b) for b in B if b][:2]
out={}
for tag,p,cols in [('sd','../Fe-bcc-bulk-crystal23-a2.87-20261105/DOSS.DAT',{'s':1,'d':2}),('spd','s0.586/DOSS.DAT',{'s':1,'p':2,'d':3})]:
    a,b=rd(p); o={}
    for sp,blk in (('up',a),('dn',b)):
        E=blk[:,0]*H; o['E']=E
        for k,c in cols.items():
            y=abs(blk[:,c])/H; o[sp+'_'+k]=y
            m=E<=0; print(tag,sp,k,'occ=%.3f'%np.trapz(y[m],E[m]),'DOS(EF)=%.3f'%np.interp(0,E,y))
    out[tag]={k:[round(float(x),4) for x in v[::2]] for k,v in o.items()}
json.dump(out,open('dos_cmp.json','w'))
