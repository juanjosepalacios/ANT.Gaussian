from pathlib import Path
import os,re,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
n=(base/'src/numeric.F90').read_text();s=(base/'src/OneDLead.F90').read_text()
def routine(text,name):
 m=re.search(r'^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'[^\n]*',text,re.I|re.M|re.S)
 if not m:raise ValueError(name)
 return m[0]
f=(base/'tests/lead_onset_fixture.f90.in').read_text().replace('@SEARCH@',routine(n,'FindChargeOnset'))
f=f.replace('@WRAPPER@','\n'.join(routine(s,x) for x in ['LowerBoundCharge','FindLowerChargeBound']))
fc=os.environ.get('FC','pgf90')
if Path(fc).name in ('ifx','ifort'):
 flags=['-i8','-r8','-recursive','-check','bounds','-O0']
elif Path(fc).name=='gfortran':
 flags=['-fdefault-integer-8','-fcheck=all','-ffree-line-length-none','-O0']
else:
 flags=['-i8','-r8','-Mrecursive','-Mbounds','-O0']
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 subprocess.check_call(['./test'],cwd=d)
