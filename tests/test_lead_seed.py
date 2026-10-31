#!/usr/bin/env python3
"""Test the actual repeated-cell generalized spectrum and lower-bound expansion."""
from pathlib import Path
import os,re,shlex,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
s=(base/'src/OneDLead.F90').read_text();numeric=(base/'src/numeric.F90').read_text()
def routine(text,name):
 pattern=rf'^  (?:recursive )?(?:complex\(double\) function|subroutine) {name}\(.*?^  end (?:function|subroutine) {name}[^\n]*'
 match=re.search(pattern,text,re.I|re.M|re.S)
 if not match:raise ValueError(name)
 return re.sub(r'#ifdef PGI[^\n]*\n.*?#endif[^\n]*\n','',match[0],flags=re.S)
fixture=(base/'tests/lead_seed_fixture.f90.in').read_text()
fixture=fixture.replace('@CHDIAG@',routine(numeric,'CHDiag'))
fixture=fixture.replace('@LEADTYPE@',re.search(r'  type T1DLead\n.*?  end type T1DLead',s,re.S)[0])
fixture=fixture.replace('@ROUTINES@','\n'.join(routine(s,name) for name in ['L1D_H0','L1D_S0','SetInitialLeadLowerBound','FindEBounds']))
fc=os.environ.get('FC','pgf90')
if Path(fc).name in ('ifx','ifort'):
 flags=['-i8','-r8','-recursive','-check','bounds','-O0'];libs=['-qmkl=sequential']
else:
 flags=['-i8','-r8','-Mrecursive','-Mbounds','-O0'];libs=['-llapack_ilp64','-lblas_ilp64']
if 'ANT_TEST_LIBS' in os.environ:libs=shlex.split(os.environ['ANT_TEST_LIBS'])
with tempfile.TemporaryDirectory(prefix='ant-lead-seed-') as directory:
 p=Path(directory);(p/'test.f90').write_text(fixture)
 subprocess.check_call([fc,*flags,'test.f90',*libs,'-o','test'],cwd=directory)
 subprocess.check_call(['./test'],cwd=directory)
 invalid=subprocess.run(['./test','invalid-overlap'],cwd=directory,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 if invalid.returncode==0 or 'overlap is not positive definite' not in invalid.stdout:
  raise AssertionError(invalid.stdout)
 print('PASS: nonpositive overlap is rejected explicitly')
