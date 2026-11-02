from pathlib import Path
import os,re,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
s=(base/'src/numeric.F90').read_text();device=(base/'src/device.F90').read_text()
m=re.search(r'^  subroutine AdaptiveGK15\(.*?^  end subroutine AdaptiveGK15',s,re.I|re.M|re.S)
f=(base/'tests/adaptive_density_fixture.f90.in').read_text().replace('@GK@',m[0])
def routine(name):
 m=re.search(r'^  (?:recursive )?subroutine '+name+r'\(.*?^  end subroutine '+name,device,re.I|re.M|re.S)
 return re.sub(r'#ifdef PGI[^\n]*\n.*?#endif[^\n]*\n','',m[0],flags=re.S)
f=f.replace('@ROUTINES@','\n'.join(routine(n) for n in ['gplus0','AdaptiveDensityContour','AdaptiveDensityBatch','AdaptiveDensityPoint']))
fc=os.environ.get('FC','pgf90');flags=['-qopenmp','-i8','-r8','-recursive','-check','bounds','-O0'] if Path(fc).name in ('ifx','ifort') else ['-mp','-i8','-r8','-Mrecursive','-Mbounds','-O0']
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'test.f90').write_text(f)
 subprocess.check_call([fc,*flags,'test.f90','-o','test'],cwd=d)
 subprocess.check_call(['./test'],cwd=d,env=dict(os.environ,OMP_STACKSIZE='128M',OMP_DYNAMIC='FALSE'))
