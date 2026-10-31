#!/usr/bin/env python3
"""Compile and exercise the production Brent routine and bisection reference."""
from pathlib import Path
import os,re,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
s=(base/'src/numeric.F90').read_text()
a=s.index('  real*8 function BrentRoot(');b=s.index('  end function BrentRoot',a)+len('  end function BrentRoot')
c=s.index('  REAL *8 FUNCTION Bisec(');d=s.index('  END FUNCTION Bisec',c)+len('  END FUNCTION Bisec')
fixture='module roots\ncontains\n'+s[a:b]+'\n'+s[c:d]+'\nend module\n'
fixture+=(base/'tests/brent_fixture.f90.in').read_text()
fc=os.environ.get('FC','pgf90')
if Path(fc).name in ('ifx','ifort'):
 flags=['-i8','-r8','-recursive','-check','bounds','-O0']
elif Path(fc).name=='gfortran':
 flags=['-fdefault-integer-8','-fcheck=all','-O0']
else:flags=['-i8','-r8','-Mrecursive','-Mbounds','-O0']
with tempfile.TemporaryDirectory(prefix='ant-brent-') as directory:
 p=Path(directory);(p/'test.f90').write_text(fixture)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=directory)
 subprocess.check_call(['./test'],cwd=directory)
