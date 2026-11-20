#!/usr/bin/env python3
"""Extract FullPolTraces from device.F90 and check FULLPOL against dense Pauli-matrix
traces on random SOC spinor models: sign convention, T and POL identities, spin rotation."""
from pathlib import Path
import os,re,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
def extract(s,name):
 m=re.search(r'^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'[^\n]*',s,re.M|re.S|re.I)
 if not m:raise ValueError(name)
 return m[0]
s=(p/'src/device.F90').read_text()
f=(p/'tests/fullpol_fixture.f90.in').read_text().replace('@DEVICE@',extract(s,'FullPolTraces'))
fc=os.environ.get('FC','pgf90')
flags=['-i8','-r8','-check','bounds','-O0'] if Path(fc).name in ('ifx','ifort') else ['-i8','-r8','-Mbounds','-O0']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 out=subprocess.run(['./test'],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 print(out.stdout,end='')
 if out.returncode!=0 or 'All FULLPOL checks passed' not in out.stdout:raise AssertionError('FULLPOL checks failed')
print('PASS: FULLPOL sign convention, dense Pauli traces, T and POL identities, spin rotation, |P| <= T')
