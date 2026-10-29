#!/usr/bin/env python3
# Extract production routines into a small analytic two-spin test fixture.
from pathlib import Path
import os, re, subprocess, tempfile
base = Path(__file__).resolve().parents[1]
source = (base / 'src/device.F90').read_text()
fixture = (base / 'tests/totcharge_fixture.f90.in').read_text()
routines = {
    'GPLUS0': ('  recursive subroutine gplus0(', '  end subroutine gplus0'),
    'DDOS': ('  recursive real*8 function DDOS(', '  end function DDOS'),
    'TOTCHARGE': ('  real*8 function TotCharge(', '  end function TotCharge'),
}
for marker, (start, end) in routines.items():
    i = source.index(start)
    j = source.index(end, i) + len(end)
    routine = re.sub(r'#ifdef PGI\n.*?#endif\n', '', source[i:j], flags=re.S)
    fixture = fixture.replace('@' + marker + '@', routine)
with tempfile.TemporaryDirectory(prefix='ant-totcharge-') as directory:
    work = Path(directory)
    (work / 'regression.f90').write_text(fixture)
    compiler = os.environ.get('FC', 'pgf90')
    if Path(compiler).name in ('ifx', 'ifort'):
        flags = ['-qopenmp', '-i8', '-r8', '-recursive', '-check', 'bounds', '-O0']
    else:
        flags = ['-mp', '-i8', '-r8', '-Mrecursive', '-Mbounds', '-O0']
    subprocess.check_call([compiler] + flags +
                          ['regression.f90', '-o', 'regression'], cwd=directory)
    subprocess.check_call(['./regression'], cwd=directory)
