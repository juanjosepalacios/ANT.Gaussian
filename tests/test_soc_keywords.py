#!/usr/bin/env python3
"""Compile parameters.F90 against stub modules and check that SCFSOC and the SOC
modifiers (SOC_FAC_*, RCUT, SOC_VAL_*, SOCFACATOM, SOCVALATOM) stop the input
reading unless SOC is also given."""
from pathlib import Path
import os,subprocess,tempfile
base=Path(__file__).resolve().parents[1]
stubs='''module preproc
 integer,parameter::MaxAtm=10000
end module
module constants
 integer,parameter::keyw_len=20,name_len=100
end module
module g09Common
contains
 integer function GetNAtoms()
 GetNAtoms=10
 end function
end module
module util
contains
 subroutine upcase(s)
 character(*)::s
 integer::i,n
 do i=1,len(s)
 n=iachar(s(i:i))
 if(n>=97.and.n<=122)s(i:i)=achar(n-32)
 enddo
 end subroutine
end module
'''
driver='''program test
 use parameters
 integer::ios
 character(50)::name='case'
 open(10,file='case.ant')
 call read_parameters(6,10,ios,name)
 call write_parameters(6)
end program
'''
cases=[
 ('none','ALPHA = 0.02\n',True,'SOC        =  F'),
 ('soc-only','SOC\n',True,'SOC        =  T'),
 ('soc-with-modifiers','SOC\nSCFSOC\nSOC_FAC_P = 2.0\nRCUT = 1.5\nSOC_VAL_D = 0.05\nSOCFACATOM\n1\n3, 1.0, 2.0, 0.0, 1.0\nSOCVALATOM\n1\n4, 0.0, 0.1, 0.0\n',True,'SCFSOC     =  T'),
 ('soc-after-scfsoc','SCFSOC\nSOC\n',True,'SCFSOC     =  T'),
 ('scfsoc','SCFSOC\n',False,'ERROR: SCFSOC requires SOC'),
 ('soc_fac_p','SOC_FAC_P = 1.0\n',False,'ERROR: SOC_FAC_P requires SOC'),
 ('soc_fac_d','SOC_FAC_D = 2.0\n',False,'ERROR: SOC_FAC_D requires SOC'),
 ('soc_fac_f','SOC_FAC_F = 2.0\n',False,'ERROR: SOC_FAC_F requires SOC'),
 ('rcut','RCUT = 1.5\n',False,'ERROR: RCUT requires SOC'),
 ('soc_val_p','SOC_VAL_P = 0.1\n',False,'ERROR: SOC_VAL_P requires SOC'),
 ('soc_val_d','SOC_VAL_D = 0.1\n',False,'ERROR: SOC_VAL_D requires SOC'),
 ('soc_val_f','SOC_VAL_F = 0.1\n',False,'ERROR: SOC_VAL_F requires SOC'),
 ('socfacatom','SOCFACATOM\n1\n3, 1.0, 2.0, 0.0, 1.0\n',False,'ERROR: SOCFACATOM requires SOC'),
 ('socvalatom','SOCVALATOM\n1\n4, 0.0, 0.1, 0.0\n',False,'ERROR: SOCVALATOM requires SOC'),
 ('first-reported','RCUT = 1.5\nSCFSOC\n',False,'ERROR: RCUT requires SOC')]
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'stubs.f90').write_text(stubs);(p/'driver.f90').write_text(driver)
 fc=os.environ.get('FC','pgf90');intel=Path(fc).name in ('ifx','ifort')
 flags=['-fpp','-DG09ROOT','-i8','-r8','-free','-extend-source','-check','bounds'] if intel else ['-Mpreprocess','-DG09ROOT','-i8','-r8','-Mextend','-Mbounds']
 subprocess.run([fc,*flags,'stubs.f90',str(base/'src/parameters.F90'),'driver.f90','-o','test'],cwd=p,check=True)
 for name,inp,success,match in cases:
  (p/'case.ant').write_text(inp)
  result=subprocess.run(['./test'],cwd=p,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert (result.returncode==0)==success,(name,result.stdout)
  assert ''.join(match.split()) in ''.join(result.stdout.split()),(name,result.stdout)
  print('PASS:',name)
