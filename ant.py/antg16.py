#!/usr/bin/env python3
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
import gaudata
AlphaMOs = "Alpha MO Coefficients"
BetaMOs = "Beta MO Coefficients"

# threshold on overlap eigenvalues for linear independence
ThrEvl = 1.e-6

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


def main():

	# parse command line arguments.
	FName,Debug,Guess,IExtp,MaxIt,Conv,OName = parse_args(sys.argv)

	# parse the input file into a binary array file object
	baf = qcb.BinAr(inputfile=FName,lenint=4)

	print( baf.nbasis )
	
	gaudata.gaudata.gaudata_init( baf.natoms, baf.ian, baf.c, baf.ne, baf.nae, baf.nbe, baf.atmchg, baf.nbasis, baf.ibfatm, baf.ibftyp )

	
main()
