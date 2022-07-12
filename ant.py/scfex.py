#!/usr/bin/env python3
#
# This script illustrates how to do an SCF calculation, using Gaussian to
# generate Fock matrices from densities each iteration.
#
# scf.py fname option,option,...
#
# fname -- name of Gaussian input file, parsed to provide the atom specifications
#          and AO basis.
#
# options -- debug=n        for print level, default 1.
#            oname=filename final orbitals and density will be stored in
#                           a new binary array file.
#            iextp=n        0 for 3 and 4 extrapolation during iterations,
#                           1 for only 4pt,  2 for no extrapolation, 3 for
#                           3pt only.  Default 0.
#            maxit=n        up to N iterations, default 1000.
#            guess=gtype    type of initial guess:
#                           "core" to diagonalize the core Hamiltonian here.
#                           "init" to have Gaussian generate its usual guess.
#                           "full" to do full SCF in Gaussian for guess.
#                           "partial" to do a few iterations in Gaussian.
#                           The default is core.
#            conv=value     convergence on rms-density change, default 1.e-8
#
# Michael Frisch, 2022
#
import sys
import os
import re
import math
import io
import numpy as np
from gauopen import QCOpMat as qco
from gauopen import QCBinAr as qcb
from gauopen import QCUtil as qcu
import scfutil as su
AlphaMOs = "Alpha MO Coefficients"
BetaMOs = "Beta MO Coefficients"

# threshold on overlap eigenvalues for linear independence
ThrEvl = 1.e-6

def main():

  # parse command line arguments.
  FName,Debug,Guess,IExtp,MaxIt,Conv,OName = parse_args(sys.argv)

  # parse the input file into a binary array file object
  baf = qcb.BinAr(inputfile=FName,lenint=4)

  breakpoint()
  
  # Generate the 1e matrices and initial guess orbitals
  if Guess is None or Guess == "core":  GuessG = "none"
  else:  GuessG = Guess
  cmoa,cmob = su.do_gau_setup(baf,Debug,GuessG)
  if Guess is None or Guess == "core":
    # Diagonalize the core Hamiltonian, using the orthonormal basis
    # provided by Gaussian.
    obas = baf.matlist[qcb.NamOBas].array.reshape(baf.nbsuse,baf.nbasis)
    H = baf.matlist[qcb.NamCore].expand().reshape(baf.nbasis,baf.nbasis)
    eiga,cmoa = su.dodiag (Debug,"core",obas,H)
    cmob = cmoa

  # Get the 1-particle matrices.
  nao = baf.nbasis
  nmo = baf.nbsuse
  nae = baf.nae
  nbe = baf.nbe
  su.printfl (nao,"ao basis functions")
  S = baf.matlist[qcb.NamOvlp].expand().reshape(nao,nao)
  H = baf.matlist[qcb.NamCore].expand().reshape(nao,nao)
  T = baf.matlist[qcb.NamKE].expand().reshape(nao,nao)
  V = H - T
  if Debug >= 2:
    qco.sqout("S",nao,nao,S.flatten(),0,0)
    qco.sqout("T",nao,nao,T.flatten(),0,0)
    qco.sqout("V",nao,nao,V.flatten(),0,0)
  Enr = baf.scalar("enucrep")
  print ("ENR",Enr)

  # initialize variables related to convergence and extapolation.
  sp12 = 0.0
  sp22 = 0.0
  icount = 0
  iflag = 0
  olde = 0.0
  E = 0.0

  # diagonalize the core Hamiltonian if that is to be used for the guess.
  PCur = su.dodens(nae,cmoa,nbe,cmob)
  ldens = PCur.size
  naosq = nao*nao

  # iterations
  su.prtime ("time before iter")
  DeltaP = np.zeros(shape=(3,ldens),dtype=np.float64)
  PO = np.copy(PCur)
  Done = False
  # update method requires MOs in the object even though these are
  # not used with dofock=density
  baf.addobj(qco.OpMat(AlphaMOs,cmoa))
  baf.addobj(qco.OpMat(BetaMOs,cmob))
  for iter0 in range(MaxIt):
    iter = iter0 + 1
    PN = np.copy(PCur)
#   extr34 both checks convergence and extrapolates the density matrices.
    icount,iflag,sp12,sp22,rmsdp,PCur,DeltaP,Done = su.extr34(nao,IExtp,icount,olde,E,sp12,sp22,Conv,PCur,DeltaP)
    if Done:  break
#   damping if we didn't extrapolate
    if iflag == 0:  PCur = 0.5*PO + 0.5*PN
    PO = np.copy(PCur)
    E,Fa,Fb = su.do_gau_fock(Debug,Enr,H,PCur,baf)
    olde = su.report_iter("iter",iter,iflag,E,None,olde,rmsdp)
    eiga,cmoa = su.dodiag (Debug,"alpha",cmoa,Fa)
    eigb,cmob = su.dodiag (Debug,"beta",cmob,Fb)
    PCur = su.dodens(nae,cmoa,nbe,cmob)
  su.prtime ("time after iter")
  if Debug > 0:
    su.printfl ("final eiga",eiga)
    qco.sqout ("Final MOA",nao,nao,cmoa.flatten(),0,0)

  # possibly generate a new binary array file with the results.
  if OName and OName is not "No":
    su.printfl("OName",OName)
    # add the SCF results
    baf.addobj(qco.OpMat(qcb.NamOrbEA,eiga,dimens=(nao,)))
    baf.addobj(qco.OpMat(qcb.NamOrbEB,eigb,dimens=(nao,)))
    baf.addobj(qco.OpMat(qcb.NamMOA,cmoa,dimens=(nao,nao)))
    baf.addobj(qco.OpMat(qcb.NamMOB,cmob,dimens=(nao,nao)))
    baf.scalar["ESCF"] = E
    baf.scalar["ETOTAL"] = E
    baf.scalar["SCF RMSDP"] = rmsdp
    bafout.write(OName)

# parse the input options
def parse_args (argv,):
  if len(argv) < 2:
    print ("usage:  scf.py input-file-name [option=value[,option=value,...]]")
    exit()
  MaxIt = 1000
  IExtp = 0
  Conv = 1.e-8
  FName = argv[1]
  if len(argv) < 3:
    Debug = 0
    Guess = None
    OName = "No"
  else:
    options = {}
    for i in range(2,len(argv)):
      options = qcu.string_to_dict(argv[i].casefold(),options)
    Debug = options.get("debug",0)
    if Debug is None:  Debug = 1
    else: Debug = int(Debug)
    OName = options.get("out","No")
    Guess = options.get("guess","default")
    IExtp = int(options.get("iextp",0))
    MaxIt = int(options.get("maxit",MaxIt))
    Conv = float(options.get("conv",Conv))
  su.printfl ("input",FName,"Guess",Guess,"IExtp",IExtp,"MaxIt",MaxIt,"Conv",Conv)
  return (FName,Debug,Guess,IExtp,MaxIt,Conv,OName)

main()
