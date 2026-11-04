#!/usr/bin/env python3
"""Compile the production Zeeman module against stub modules and check CompHZM
against muB*B.sigma for a global field and for per-atom ZMATOM fields."""
from pathlib import Path
import os,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
f=(p/'tests/zeeman_fixture.f90.in').read_text().replace('@ZEEMAN@',(p/'src/Zeeman.F90').read_text(encoding='latin-1'))
fc=os.environ.get('FC','pgf90')
flags=['-i8','-r8','-check','bounds','-O0'] if Path(fc).name in ('ifx','ifort') else ['-i8','-r8','-Mbounds','-O0']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f,encoding='latin-1')
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 subprocess.check_call(['./test'],cwd=d)
print('PASS: CompHZM global field, per-atom ZMATOM fields, unlisted atoms without field')
