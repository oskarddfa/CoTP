#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-2016, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#
#  SLEPc is free software: you can redistribute it and/or modify it under  the
#  terms of version 3 of the GNU Lesser General Public License as published by
#  the Free Software Foundation.
#
#  SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY
#  WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS
#  FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for
#  more details.
#
#  You  should have received a copy of the GNU Lesser General  Public  License
#  along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
DIRS     = sys eps svd pep nep mfn

LOCDIR   = src/

CC = g++
LINKERFLAG = -lm -lslepc -lpetsc -lX11 -lpthread -llapack -lblas
default: uebung2.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

uebung2.o:  uebung2_.cpp chkopts
	-${CLINKER} uebung2_.cpp ${SLEPC_EPS_LIB} ${PETSC_CCPPFLAGS}  ${SLEPC_INCLUDE} ${LINKERFLAG}
