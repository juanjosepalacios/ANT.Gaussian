//**********************************************************
//*                                                        *
//*      ANT.G09-2.4.1 - system.c                          *
//*                                                        *
//**********************************************************
//*                                                        *
//*  Copyright (c) by                                      *
//*                                                        *
//*  Juan Jose Palacios (1)                                *
//*  David Jacob (2)                                       *
//*                                                        *
//* (1) Departamento de Fisica de la Materia Condensada    *
//*     Universidad Autonoma de Madrid                     *
//*     28049 Madrid (SPAIN)                               *
//* (2) Theory Department                                  *
//*     Max-Planck-Institute for Microstructure Physics    *
//*     Halle, 06120 (GERMANY)                             *
//*                                                        *
//**********************************************************

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//************************************************
//* C function to obtain an environment variable *
//************************************************
void 
cgetenv_( envval, len, envname  )
     int *len;
     char *envval;
     char *envname;
{
  strcpy(envval, getenv( envname ));
  *len=strlen( envval );
}
