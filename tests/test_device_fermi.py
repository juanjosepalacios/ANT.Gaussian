#!/usr/bin/env python3
"""Test production device Fermi dispatch and Brent with stateful charge callbacks."""
from pathlib import Path
import os,re,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
def extract(s,name):
 m=re.search(r'^  (?:subroutine|real ?\*8 function) '+name+r'\(.*?^  end (?:subroutine|function) '+name+r'[^\n]*',s,re.M|re.S|re.I)
 if not m:raise ValueError(name)
 return m[0]
f=(p/'tests/device_fermi_fixture.f90.in').read_text()
n=(p/'src/numeric.F90').read_text();s=(p/'src/device.F90').read_text()
f=f.replace('@NUMERIC@','\n'.join(extract(n,x) for x in ['BracketRoot','BrentRoot','MULLER','SECANT','Bisec']))
f=f.replace('@DEVICE@','\n'.join(extract(s,x) for x in ['SolveDeviceFermi','SolveDeviceFermiMuller','CompDensMat','CompDensMat_SOC','CompDensMat2','CompDensMat2_SOC']))
fc=os.environ.get('FC','pgf90')
flags=['-i8','-r8','-recursive','-check','bounds','-qopenmp','-O0'] if Path(fc).name in ('ifx','ifort') else ['-i8','-r8','-Mrecursive','-Mbounds','-mp','-O0']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 subprocess.check_call(['./test'],cwd=d)
 for mode in ['no-bracket','nonfinite','charge-failure']:
  bad=subprocess.run(['./test',mode],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
  if bad.returncode==0 or 'ERROR:' not in bad.stdout:raise AssertionError(mode+bad.stdout)
 print('PASS: bracket, nonfinite and unacceptable final charge fail explicitly')
