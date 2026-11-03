#!/usr/bin/env python3
"""Exercise the actual EMin/EMax spectral-weight widening policy with controlled integrals."""
from pathlib import Path
import re,os,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
s=(p/'src/device.F90').read_text()
r=re.search(r'^  subroutine CheckEnergyBounds\(.*?^  end subroutine CheckEnergyBounds',s,re.M|re.S)[0]
f=(p/'tests/spectral_windows_fixture.f90.in').read_text().replace('@CHECK@',r)
fc=os.environ.get('FC','pgf90')
flags=['-i8','-r8','-check','bounds','-O0'] if Path(fc).name in ('ifx','ifort') else ['-i8','-r8','-Mbounds','-O0']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 for mode in ['first','widen','excess','quadrature','nonfinite','tail','outside','cap']:
  subprocess.check_call(['./test',mode],cwd=d)
 bad=subprocess.run(['./test','failure'],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 if bad.returncode==0 or 'failed after widening EMin/EMax' not in bad.stdout:raise AssertionError(bad.stdout)
 print('PASS: [EMin,EMax] check, all widening reasons, widened bounds kept, explicit final failure')
 print('PASS: CHARGEACC-scaled tolerance (capped at 0.1), quadrature at a tenth of it, ETA tail added back only for states inside the window')
