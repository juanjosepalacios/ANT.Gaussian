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
 print *, 'COUNTS',NatPC,NEmbed
end program
'''
cases=[
 ('new', 'TYPE1 = 1DLEAD\nTYPE2 = 1DLEAD\nNATPC1 = 1\nNATPC2 = 3\n',True,'COUNTS'),
 ('old','NEMBED1 = 1\nNEMBED2 = 3\nTYPE1 = 1DLEAD\nTYPE2 = 1DLEAD\n',True,'deprecated'),
 ('equal','NATPC1 = 3\nNEMBED1 = 3\nTYPE1 = 1DLEAD\n',True,'deprecated'),
 ('conflict','NEMBED1 = 1\nTYPE1 = 1DLEAD\nNATPC1 = 3\n',False,'conflicting'),
 ('missing','TYPE1 = 1DLEAD\n',False,'specify positive'),
 ('zero','TYPE1 = 1DLEAD\nNATPC1 = 0\n',False,'specify positive'),
 ('negative','TYPE1 = 1DLEAD\nNATPC1 = -2\n',False,'specify positive'),
 ('negative-sentinel','TYPE1 = 1DLEAD\nNATPC1 = -1\nNEMBED1 = 3\n',False,'specify positive'),
 ('wrong-type','NATPC1 = 3\n',False,'requires TYPE'),
 ('bethe','NEMBED1 = 3\nNEMBED2 = 2\n',True,'COUNTS')]
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'stubs.f90').write_text(stubs);(p/'driver.f90').write_text(driver)
 fc=os.environ.get('FC','pgf90');intel=Path(fc).name in ('ifx','ifort')
 flags=['-fpp','-DG09ROOT','-i8','-r8','-free','-extend-source','-check','bounds'] if intel else ['-Mpreprocess','-DG09ROOT','-i8','-r8','-Mextend','-Mbounds']
 subprocess.run([fc,*flags,'stubs.f90',str(base/'src/parameters.F90'),'driver.f90','-o','test'],cwd=p,check=True)
 outputs={}
 for name,inp,success,match in cases:
  (p/'case.ant').write_text(inp)
  result=subprocess.run(['./test'],cwd=p,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  assert (result.returncode==0)==success,(name,result.stdout)
  assert match in result.stdout,(name,result.stdout)
  outputs[name]=result.stdout
  print('PASS:',name)
 assert outputs['new'].split('COUNTS')[-1].split()==outputs['old'].split('COUNTS')[-1].split()
 print('PASS: legacy and NATPC inputs resolve to identical counts')
