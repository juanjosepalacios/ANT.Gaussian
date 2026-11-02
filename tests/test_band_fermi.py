#!/usr/bin/env python3
"""Extract production Bloch and filling routines; test analytic tight-binding bands."""
from pathlib import Path
import os,re,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
s=(p/'src/OneDLead.F90').read_text();num=(p/'src/numeric.F90').read_text()
def extract(s,name):
 m=re.search(r'^  (?:recursive )?(?:subroutine|real\(double\) function) '+name+r'\(.*?^  end (?:subroutine|function) '+name+r'[^\n]*',s,re.M|re.S|re.I)
 if not m:raise ValueError(name)
 return re.sub(r'#ifdef PGI[^\n]*\n.*?#endif[^\n]*\n','',m[0],flags=re.S)
f=(p/'tests/band_fermi_module.f90.in').read_text().split('program test_seed')[0]
f=f.replace('real*8 :: ChargeAcc=1d-2,FermiAcc=1d-2,Infty=10000d0','real*8 :: ChargeAcc=1d-7,FermiAcc=1d-7,Infty=10000d0,eta=0.02d0')
f=f.replace('@CHDIAG@',extract(num,'CHDiag')).replace('@LEADTYPE@',re.search(r'  type T1DLead\n.*?  end type T1DLead',s,re.S)[0])
f=f.replace('@ROUTINES@','\n'.join(extract(s,n) for n in ['LeadBandPoint','BandPopulation','BandSearchFailure','FindBandFermi']))
f+=(p/'tests/band_fermi_fixture.f90.in').read_text()
fc=os.environ.get('FC','pgf90')
flags=['-i8','-r8','-recursive','-check','bounds','-qopenmp','-O0'] if Path(fc).name in ('ifx','ifort') else ['-i8','-r8','-Mrecursive','-Mbounds','-mp','-O0']
libs=['-qmkl=sequential'] if Path(fc).name in ('ifx','ifort') else ['-llapack_ilp64','-lblas_ilp64']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90',*libs,'-o','test'],cwd=d)
 for threads in [1,8]:
  subprocess.check_call(['./test'],cwd=d,env=dict(os.environ,OMP_NUM_THREADS=str(threads)))

 bad=subprocess.run(['./test','invalid-eta'],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 if bad.returncode==0 or 'Invalid band-search tolerances/ETA' not in bad.stdout:raise AssertionError(bad.stdout)
 print('PASS: invalid broadening returns nonzero failure')
