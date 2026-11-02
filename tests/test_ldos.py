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
 ('off','',True,'LDOS       =  F'),
 ('valid','LDOS\nLDOS_BEG = 2\nLDOS_END = 4\n',True,'LDOS       =  T'),
 ('order','LDOS_END = 4\nLDOS_BEG = 2\nLDOS\n',True,'LDOS       =  T'),
 ('single','LDOS\nLDOS_BEG = 3\nLDOS_END = 3\n',True,'LDOS       =  T'),
 ('bounds-only','LDOS_BEG = 2\nLDOS_END = 4\n',True,'local DOS is disabled'),
 ('missing-both','LDOS\n',False,'requires both'),
 ('missing-end','LDOS\nLDOS_BEG = 2\n',False,'requires both'),
 ('missing-beg','LDOS\nLDOS_END = 4\n',False,'requires both'),
 ('zero','LDOS\nLDOS_BEG = 0\nLDOS_END = 4\n',False,'number of atoms'),
 ('reversed','LDOS\nLDOS_BEG = 4\nLDOS_END = 2\n',False,'number of atoms'),
 ('past-end','LDOS\nLDOS_BEG = 2\nLDOS_END = 11\n',False,'number of atoms')]
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
  assert ''.join(match.split()) in ''.join(result.stdout.split()),(name,result.stdout)
  outputs[name]=result.stdout
  print('PASS:',name)
