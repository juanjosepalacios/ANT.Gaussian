#!/usr/bin/env python3
"""Exercise actual 1D charge, Dyson, lead accessors and LAPACK on model leads.

FC selects pgf90 (default), ifx/ifort, or gfortran. ANT_TEST_LIBS can
supply ABI-compatible LAPACK/BLAS link flags (production uses 8-byte integers).
"""
from pathlib import Path
import os
import re
import shlex
import subprocess
import tempfile

base = Path(__file__).resolve().parents[1]
source = (base / 'src/OneDLead.F90').read_text()
numeric = (base / 'src/numeric.F90').read_text()
fixture = (base / 'tests/onedlead_charge_fixture.f90.in').read_text()

def routine(text, name):
    pattern = rf'^  (?:recursive )?(?:real\(double\) function|complex\(double\) function|integer function|subroutine) {name}\(.*?^  end (?:function|subroutine) {name}[^\n]*'
    match = re.search(pattern, text, re.I | re.M | re.S)
    if not match:
        raise ValueError(name)
    return re.sub(r'#ifdef PGI[^\n]*\n.*?#endif[^\n]*\n', '', match[0], flags=re.S)

names = ['CompSelfEnergy1D', 'SolveDyson1D', 'CompGreensFunc',
         'CompGreensFuncWorkspace', 'ChargePoint1D', 'TotCharge', 'CreateXSC',
         'L1D_H0', 'L1D_V1', 'L1D_V2', 'L1D_S0', 'L1D_S1', 'L1D_S2']
code = '\n'.join(routine(source, name) for name in names)
# Test-only observation inside the real Green-function path, not an extra team.
code = code.replace('    ispin = spin', '''!$omp critical(test_team)
    max_team=max(max_team,omp_get_num_threads())
    if (omp_get_level()>1) max_nested=max(max_nested,omp_get_num_threads())
!$omp end critical(test_team)
    ispin = spin''')
fixture = fixture.replace('@ROUTINES@', code)
fixture = fixture.replace('@LEADTYPE@', re.search(r'  type T1DLead\n.*?  end type T1DLead', source, re.S)[0])
fixture = fixture.replace('@CINV@', routine(numeric, 'CInv'))
a = numeric.index('  SUBROUTINE gauleg(')
b = numeric.index('  END SUBROUTINE gauleg', a) + len('  END SUBROUTINE gauleg')
fixture = fixture.replace('@GAULEG@', numeric[a:b])
compiler = os.environ.get('FC', 'pgf90')
name = Path(compiler).name
if name in ('ifx', 'ifort'):
    flags = ['-qopenmp', '-i8', '-r8', '-recursive', '-check', 'bounds', '-O0']
    libs = ['-qmkl=sequential', '-DMKL_ILP64']
elif name == 'gfortran':
    flags = ['-fopenmp', '-fdefault-integer-8', '-frecursive', '-fcheck=all', '-O0']
    libs = ['-llapack', '-lblas']  # Override for an ILP64 installation.
else:
    flags = ['-mp', '-i8', '-r8', '-Mrecursive', '-Mbounds', '-O0']
    libs = ['-llapack_ilp64', '-lblas_ilp64']
libs = shlex.split(os.environ['ANT_TEST_LIBS']) if 'ANT_TEST_LIBS' in os.environ else libs
with tempfile.TemporaryDirectory(prefix='ant-onedlead-') as directory:
    work = Path(directory)
    (work / 'regression.f90').write_text(fixture)
    subprocess.check_call([compiler, *flags, 'regression.f90', *libs, '-o', 'regression'], cwd=work)
    subprocess.check_call(['./regression'], cwd=work,
                          env={**os.environ, 'OMP_STACKSIZE': '64M', 'MKL_NUM_THREADS': '1', 'OMP_DYNAMIC': 'FALSE'})
