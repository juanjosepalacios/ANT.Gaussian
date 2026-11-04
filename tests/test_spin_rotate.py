#!/usr/bin/env python3
"""Compile the production SpinRotate module against stub modules and check CompHROT
against spin projectors (1 +/- n.sigma)/2 for global and per-atom rotations."""
from pathlib import Path
import os,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
f=(p/'tests/spin_rotate_fixture.f90.in').read_text().replace('@SPINROTATE@',(p/'src/SpinRotate.F90').read_text())
fc=os.environ.get('FC','pgf90')
flags=['-i8','-r8','-check','bounds','-O0'] if Path(fc).name in ('ifx','ifort') else ['-i8','-r8','-Mbounds','-O0']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 subprocess.check_call(['./test'],cwd=d)
print('PASS: CompHROT rigid global rotation, (90,90) along +y, per-atom on-site and averaged hopping blocks, spin-diagonal overlap')
