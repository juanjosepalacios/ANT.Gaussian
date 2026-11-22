#!/usr/bin/env python3
"""Extract the WBL kernels (WBLLowdin, WBLGreen, WBLTransMat, WBLMoleculeOnly) from device.F90
and check them on a small random non-orthogonal device with energy-dependent leads: T_WBL(E_F)
equals the normal T(E_F), T_WBL(E) against a dense Loewdin reference and the Schur form, T_WBL
is the same for two AO representations of one orthogonal model (the AO z*S_MM form is not),
T_WBL stays below the channel count, and the molecule-only diagnostic.
Needs LAPACK/BLAS with 8-byte integers (production uses -i8)."""
from pathlib import Path
import os,re,subprocess,tempfile
p=Path(__file__).resolve().parents[1]
def extract(s,name):
 m=re.search(r'^  subroutine '+name+r'\(.*?^  end subroutine '+name+r'[^\n]*',s,re.M|re.S|re.I)
 if not m:raise ValueError(name)
 return re.sub(r'#ifdef PGI[^\n]*\n.*?#endif[^\n]*\n','',m[0],flags=re.S)
s=(p/'src/device.F90').read_text()
kernels='\n'.join(extract(s,k) for k in ('WBLLowdin','WBLGreen','WBLTransMat','WBLMoleculeOnly'))
f=(p/'tests/wbl_fixture.f90.in').read_text().replace('@DEVICE@',kernels)
fc=os.environ.get('FC','pgf90')
if Path(fc).name in ('ifx','ifort'):flags=['-i8','-r8','-check','bounds','-O0'];libs=['-qmkl=sequential']
else:flags=['-i8','-r8','-Mbounds','-O0'];libs=['-llapack_ilp64','-lblas_ilp64']
with tempfile.TemporaryDirectory() as d:
 Path(d,'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90',*libs,'-o','test'],cwd=d)
 out=subprocess.run(['./test'],cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 print(out.stdout,end='')
 if out.returncode!=0 or 'All WBL checks passed' not in out.stdout:raise AssertionError('WBL checks failed')
print('PASS: WBL E_F identity, Loewdin dense and Schur references, overlap invariance, channel bound, molecule-only diagnostic')
