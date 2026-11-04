#!/usr/bin/env python3
"""Compile parameters.F90 against stub modules and check that SPINROTATOM alone
switches on the spin rotation and that the removed SPINROT, THETA and PHI
keywords stop the input reading."""
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
 print *,'ROTATION',rot,NSpinRotAtom
end program
'''
cases=[
 ('none','ALPHA = 0.02\n',True,'ROTATION F 0'),
 ('spinrotatom','SPINROTATOM\n2\n5, 90.0, 0.0\n6, -30.0, 45.0\n',True,'ROTATION T 2'),
 ('spinrotatom-logged','SPINROTATOM\n1\n6, -30.0, 45.0\n',True,'6 -30.0000 45.0000'),
 ('spinrot','SPINROT\n',False,'Undefined keyword: SPINROT'),
 ('theta','THETA = 90.0\n',False,'Undefined keyword: THETA'),
 ('phi','PHI = 90.0\n',False,'Undefined keyword: PHI')]
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
