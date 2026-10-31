from pathlib import Path
import os,re,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
s=(base/'src/OneDLead.F90').read_text();a=s.index('  recursive subroutine SolveDyson1D(');b=s.index('  end subroutine SolveDyson1D',a)+len('  end subroutine SolveDyson1D')
routine=re.sub(r'#ifdef PGI[^\n]*\n.*?#endif[^\n]*\n','',s[a:b],flags=re.S)
fixture=(base/'tests/dyson_mixing_fixture.f90.in').read_text().replace('@SOLVER@',routine)
fc=os.environ.get('FC','pgf90')
if Path(fc).name in ('ifx','ifort'):
 flags=['-i8','-r8','-recursive','-check','bounds','-O0'];libs=['-qmkl=sequential']
else:flags=['-i8','-r8','-Mrecursive','-Mbounds','-O0'];libs=['-llapack_ilp64','-lblas_ilp64']
with tempfile.TemporaryDirectory(prefix='ant-dyson-mix-') as directory:
 p=Path(directory);(p/'test.f90').write_text(fixture)
 subprocess.check_call([fc,*flags,'test.f90',*libs,'-o','test'],cwd=directory)
 subprocess.check_call(['./test'],cwd=directory)
