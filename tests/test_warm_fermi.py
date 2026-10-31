#!/usr/bin/env python3
"""Test actual FindFermi with distinct electrode models and real bracket/Brent code."""
from pathlib import Path
import os,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
num=(base/'src/numeric.F90').read_text()
lead=(base/'src/OneDLead.F90').read_text()
def extract(text,start,end):
 a=text.index(start);b=text.index(end,a)+len(end);return text[a:b]
bracket=extract(num,'  subroutine BracketRoot(','  end subroutine BracketRoot')
# Observe the center passed to the production bracket routine, without changing logic.
bracket=bracket.replace('    first=.true.', '''    test_bracket_calls=test_bracket_calls+1
    test_centers(test_bracket_calls)=center
    first=.true.''',1)
brent=extract(num,'  real*8 function BrentRoot(','  end function BrentRoot')
counted=extract(lead,'  real(double) function CountedLeadCharge(','  end function CountedLeadCharge')
find=extract(lead,'  subroutine FindFermi(','  end subroutine FindFermi')
fixture=(base/'tests/warm_fermi_fixture.f90.in').read_text()
fixture=fixture.replace('@BRACKET@',bracket).replace('@BRENT@',brent).replace('@COUNTED@',counted).replace('@FIND@',find)
# Read production warm-start declarations, including their per-electrode dimensions.
a=lead.index('  real(double) :: WarmLeadFermi(')
b=lead.index('\n',lead.index('  real(double), parameter :: WarmFermiHalfWidth=',a))
fixture=fixture.replace('@HISTORY@',lead[a:b])
fc=os.environ.get('FC','pgf90')
if Path(fc).name in ('ifx','ifort'):
 flags=['-qopenmp','-i8','-r8','-recursive','-check','bounds','-O0']
elif Path(fc).name=='gfortran':
 flags=['-fopenmp','-fdefault-integer-8','-fcheck=all','-O0']
else:flags=['-mp','-i8','-r8','-Mrecursive','-Mbounds','-O0']
with tempfile.TemporaryDirectory(prefix='ant-warm-fermi-') as directory:
 p=Path(directory);(p/'test.f90').write_text(fixture)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=directory)
 subprocess.check_call(['./test'],cwd=directory)
