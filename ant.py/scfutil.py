import sys
import time
import math
import numpy as np
from gauopen import QCOpMat as qco
from gauopen import QCBinAr as qcb

AlphaDensity = "Alpha SCF Density Matrix"
BetaDensity = "Beta SCF Density Matrix"
AlphaFock = "Alpha Fock Matrix"
BetaFock = "Beta Fock Matrix"

# check for one option, returning a default value if it is not present.
def get1opt (options,key,defval):
  key = key.casefold()
  if options.get(key,"default") == "default":
    print("mjf default for",key)
    return (False,defval)
  else:
    print("mjf use opt for",key)
    sval = options.get(key,defval)
    if type(defval) == type(1):  val = int(sval)
    elif type(defval) == type(1.1):  val = float(sval)
    else:  val = sval
    return (True,val)

# form the Fock matrices (alpha, beta) from the merged density matrix
# using Gaussian and compute the electronic energy.
def do_gau_fock(Debug,ENR,H,PCur,baf):
  nao,junk = np.shape(H)
  naosq = nao*nao
  Pa = PCur[0:naosq].reshape(nao,nao)
  Pb = PCur[naosq:2*naosq].reshape(nao,nao)
  # have Gaussian compute integrals and the Fock matrices
  PaO = qco.OpMat(AlphaDensity,Pa,dimens=(nao,nao))
  PaO.compress()
  baf.addobj(PaO)
  PbO = qco.OpMat(BetaDensity,Pb,dimens=(nao,nao))
  PbO.compress()
  baf.addobj(PbO)
  # force uhf (separate alpha and beta densities) in Gaussian
  baf.iopcl = 1
  baf.update(dofock="density")
  Fa = baf.matlist[AlphaFock].expand()
  Fb = baf.matlist[BetaFock].expand()
  # compute the electronic energy
  Ea = 0.5*np.einsum('ls,ls',Pa,Fa+H)
  Eb = 0.5*np.einsum('ls,ls',Pb,Fb+H)
  E = Ea + Eb + ENR
  return(E,Fa,Fb)

# diagonalize one Fock matrix in the basis of the current orbitals to produce new orbitals.
def dodiag (Debug,label,cmo,F):
  Fmo = np.einsum('im,mn,jn->ij',cmo,F,cmo,optimize=True)
  if Debug >= 3:  printfl ("Fmo",label,Fmo)
  oe,evec = np.linalg.eigh(Fmo)
  d1,d2 = np.shape(evec)
  for i in range(d1):
    emax=abs(evec[1,i])
    mmax=0
    for m in range(1,d1):
      if abs(evec[m,i]) > emax:
        mmax = m
        emacs = abs(evec[m,i])
    if evec[mmax,i] < 0:  evec[:,i] = -evec[:,i]
  ncmo = np.einsum('pi,pm->im',evec,cmo)
  if Debug >= 3:
    printfl ("cmo",label,cmo)
    printfl ("ncmo",label,ncmo)
    printfl ("dcmo",label,ncmo-cmo)
  return (oe,ncmo)

# print which flushes output, so that output from running
# jobs comes out continuously.
def printfl (*args):
  print (*args)
  sys.stdout.flush()

# time stamp routine for tuning.
def prtime (label):
  printfl ("%-20s" % label,time.asctime(time.localtime(time.time())))

# Form the density matrices (full nao+nno square) from the orbitals together
# in a single vector.  The order is alpha,beta,protons with the latter
# two optional.
def dodens(nae,cmoa,nbe=-1,cmob=None,nump=-1,cmon=None):
  cmoao = cmoa[:nae,:]
  Pa = np.einsum('im,in->mn',cmoao,cmoao)
  naou,nao = np.shape(cmoa)
  naosq = nao*nao
  nptot = naosq
  if nbe >= 0:
    cmobo = cmob[:nbe,:]
    Pb = np.einsum('im,in->mn',cmobo,cmobo)
    nptot += naosq
  indn = nptot
  if nump > 0:
    cmono = cmon[:nump,:]
    Pn = np.einsum('im,in->mn',cmono,cmono)
    nnou,nno = np.shape(cmon)
    nptot += nno*nno
  else:
    nnou = 0
    nno = 0
  PCur = np.zeros(shape=(nptot,))
  PCur[0:naosq] = Pa.flatten()
  if nbe >=0:  PCur[naosq:2*naosq] = Pb.flatten()
  if nump > 0:  PCur[indn:] = Pn.flatten()
  return (PCur)

# Generate 1e matrcies and some orbirals in Gaussian for an initial guess.
# Guess specifies how to get the orbitals:
#  "none"     -- just generate 1e matrices, don't return any orbitals.
#  "init"     -- return Gaussian's usual initial guess.
#  "partial"  -- do a partially-converged SCF in Gaussian.
#  "full"     -- fully converge the SCF in Gaussian.
# 
def do_gau_setup (baf,Debug,Guess,logfile=None):
  miscr = ""
  if Debug >= 1:  miscr = miscr + " scf=conventional noraff"
  if Guess == "none" or Guess == "init":  miscr = miscr + " guess=only"
  elif Guess == "partial":  miscr = miscr + " scfcon=2"
  elif Guess != "full":
    printfl("Unrecgonized guess type",Guess)
    exit
  # run Gaussian to update the binary array file object, in this case
  # by doing an SCF.
  baf.update (dofock="SCF",toutput=logfile,miscroute=miscr)
  enr0 = baf.scalar("enucrep")
  eref = baf.scalar("etotal")
  printfl ("got nuc rep",enr0,"energy",eref)
  if Guess == "none":
    cmoa = None
    cmob = None
  else:
    cmoa = baf.matlist[qcb.NamMOA].array.reshape(baf.nbsuse,baf.nbasis)
    if qcb.NamMOB in baf.matlist:
      cmob = baf.matlist[qcb.NamMOB].array.reshape(baf.nbsuse,baf.nbasis)
    else:
      cmob = cmoa
  if Debug >= 1:
    cmoao = cmoa[:baf.nae,:]
    cmobo = cmob[:baf.nbe,:]
    Pa = np.einsum('im,in->mn',cmoao,cmoao)
    Pb = np.einsum('im,in->mn',cmobo,cmobo)
    Pt = Pa + Pb
    H = baf.matlist["CORE HAMILTONIAN ALPHA"].expand().reshape(baf.nbasis,baf.nbasis)
    eone = np.einsum('ls,ls',Pt,H)
    printfl("init eone",eone)
    if qcb.Nam2E in baf.matlist:
      R = baf.matlist[qcb.Nam2E].expand().reshape(baf.nbasis,baf.nbasis,baf.nbasis,baf.nbasis)
      pjp = 0.5 * np.einsum('mn,ls,mnls',Pt,Pt,R,optimize=True)
      pkpa = 0.5 * np.einsum('mn,ls,mlns',Pa,Pa,R,optimize=True)
      pkpb = 0.5 * np.einsum('mn,ls,mlns',Pb,Pb,R,optimize=True)
      pgp = pjp - pkpa -pkpb
      printfl ("init pjp",pjp,"pkp",pkpa,pkpb,"pgp",pgp)
  return (cmoa,cmob)

# compute changes in the (merged) density matrix and use these to
# check convergence and possibly extrapolate the density.
# This routines does simple 3 and 4 point extrapolation of the
# density, which is inferior to DIIS but is simpler.
def extr34 (nao,IExtp,icount,olde,E,sp12,sp22,Conv,PCur,DeltaP):
  tiny = 1.e-12
  small = 1.e-30
  iflag = 0
  Done = False
  icount = icount + 1
  loc1 = 2 - (icount % 2)
  loc2 = 3 - loc1
  A1 = np.copy(DeltaP[0,:])
  DeltaP[0,:] = np.copy(PCur)
  # nothing else to do in first cycle or immediately after extrapolation
  rmsdp = 0.0
  if icount > 1:
    # p(n) now in a3 and p(n-1) in a1; form p(n)-p(n-1) in a2.
    A2 = PCur - A1
    # find length dp1
    sp11 = np.dot(A2,A2)
    naosq = nao*nao
    if nao > 0:
      rmsdpa = math.sqrt(np.dot(A2[0:naosq],A2[0:naosq])/naosq)
      rmsdpb = math.sqrt(np.dot(A2[naosq:2*naosq],A2[naosq:2*naosq])/naosq)
    else:
      rmsdpa = 0.0
      rmsdpb = 0.0
    # test for convergence by.
    # printfl("rmsdpa,b",rmsdpa,rmsdpb)
    rmsdp = math.sqrt(sp11/A2.size)
    Done = (rmsdp < Conv) and (E <= olde)
    if not Done:
      dp1 = math.sqrt(sp11*0.5)
      DeltaP[loc1,:] = np.copy(A2)
      if icount > 3:
        A1 = np.copy(DeltaP[loc1,:])
        sp23 = sp12
        sp33 = sp22
        sp13 = np.dot(A1,A2)
        dp3 = math.sqrt(sp33*0.5)
        # read p(n-1)-p(n-2) into a1
        A1 = np.copy(DeltaP[loc2,:])
        sp12 = np.dot(A1,A2)
        sp22 = np.dot(A1,A1)
        dp2 = math.sqrt(sp22*0.5)
        # find cosine of angle between successive displacements
        denom = 2.0*dp1*dp2
        do34 = denom > small
        if do34: cosphi = sp12 / denom
        else:  cosphi = 0.0
        # find cosine of angle between dp3 and plane of dp1 and dp2
        denom = sp11*sp22 - sp12*sp12
        do34 = do34 and (dp3 > small) and (denom > small)
        # do not extrapolate unless 4 consecutive points are nearly coplanar
        if do34:
          x = (sp13*sp22-sp12*sp23) / denom
          y = (sp23*sp11-sp12*sp13) / denom
          arsqrt = (x*x*sp11+y*y*sp22+2.0*x*y*sp12)*0.5
          do34 = arsqrt > 0.0
          if do34:
            cospsi = math.sqrt(arsqrt)/dp3
            do34 = (abs(x) > small) and (cospsi > 0.99)
        if do34:
          # express vector dp1 as x*dp3(projected)+y*dp2
          y = -y/x
          x = 1.0/x
          # test if 2*2 matrix has real eigenvalues between -.95 and +.95
          xy = y*y + 4.0*x
          xy1 = abs(y) + math.sqrt(abs(xy))
          ok4pt = (xy >= 0.0) and (xy1 <= 1.9)
          if ok4pt:
            if IExtp == 2 or IExtp == 3:
              iflag = 4
            else:
              xxx = x / (1.0-x-y)
              yyy = (x+y) / (1.0-x-y)
              PCur = PCur + xxx*A1 + yyy*A2
              iflag = 2
              icount = 0
          if (abs(cosphi)>0.995) and ((iflag==4 and IExtp==3) or not ok4pt):
            if IExtp == 1 or IExtp == 2:
              iflag = 3
            else:
              x = dp1 / (dp2*cosphi-dp1)
              PCur = PCur + x*A2
              iflag = 1
              icount = 0
        if iflag == 1 or iflag == 2:  DeltaP[0,:] = np.copy(PCur)
  return(icount,iflag,sp12,sp22,rmsdp,np.copy(PCur),np.copy(DeltaP),Done)


# report the energy at an iteration or what extrapolation was done
# (in which case the energy is non-variational and not reported).
def report_iter (label,iter,iflag,E,ETn,olde,rmsdp):
  print (label,iter,end=" ")
  if iflag == 1:
    printfl ("3-point extrapolation")
  elif iflag == 2:
    printfl ("4-point extrapolation")
  else:
    if ETn is None:  print ("EE=",E,end="")
    else:  print ("EE=",E-ETn,"E=",E,end="")
    if iter == 1:  printfl()
    else: printfl (" delta-E=",E-olde,"rmsdp=",rmsdp)
    olde = E
  return (olde)
