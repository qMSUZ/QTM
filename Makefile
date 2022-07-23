#
#/***************************************************************************
# *   Copyright (C) 2017 -- 2018 by Marek Sawerwain                         *
# *                                         <M.Sawerwain@gmail.com>         *
# *                                                                         *
# *   Part of the Quantum Trajectory Method:                                *
# *   https://github.com/qMSUZ/QTM                                          *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free Software Foundation; either version 3 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU General Public License for more details.                          *
# *                                                                         *
# *   You should have received a copy of the GNU General Public License     *
# *   along with this program; if not, write to the                         *
# *   Free Software Foundation, Inc.,                                       *
# *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
# ***************************************************************************/
#

#
#
# Makefile for Quantum Trajectory Method with MPI and ZVODE
#
#


FGET 	= wget
FC      = gfortran
MGXX    = g++
CC      = mpic++
#CC      = g++
RUN     = mpirun
CFLAGS  = -I. -fpermissive -O2
FCFLAGS = -I. -O2

OUTPUTNAME_U = unitary-ex
OUTPUTNAME_TRIHAM = triham-ex
OUTPUTNAME_BDPC = bdpc-ex
OUTPUTNAME_JCM = jcm-ex
OUTPUTNAME_PC = pumped-cavity-ex

UNITARY_MAIN_SRC = unitary-mc-ex.cc rgen_lfsr113.cc
UNITARY_MAIN_SRC_DIRECT = unitary-mc-ex-direct.cc rgen_lfsr113.cc
TRIHAM_OBJ_SRC = triham-mc-ex.cc rgen_lfsr113.cc
TRIHAM_OBJ_SRC_DIRECT = triham-mc-ex-direct.cc rgen_lfsr113.cc
BDPC_OBJ_SRC = bdpc-ex.cc rgen_lfsr113.cc
BDPC_MAIN_SRC_DIRECT = bdpc-ex-direct.cc rgen_lfsr113.cc
JCM-MC-EX_SRC = jaynes-cummings-model-mc-ex.cc rgen_lfsr113.cc
JCM-MC-EX_SRC_DIRECT = jaynes-cummings-model-mc-ex-direct.cc rgen_lfsr113.cc
PC_MAIN_SRC = pumped-cavity-ex.cc rgen_lfsr113.cc

ZVODE_SRC = zvode.f zgesl.f zgefa.f zgbsl.f zgbfa.f

#UNITARY_MAIN_SRC=$(QTM_MAIN_SRC:.cc=.o)
QTM_OBJ_ZVODE_SRC=$(ZVODE_SRC:.f=.o)

LIBRARY=

all:

unitary-ex: $(UNITARY_MAIN_SRC) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(UNITARY_MAIN_SRC) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_U) -lblas -llapack -lgfortran $(LIBRARY)

unitary-ex-direct: $(UNITARY_MAIN_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(UNITARY_MAIN_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_U) -lblas -llapack -lgfortran $(LIBRARY)

		
triham-ex: $(TRIHAM_OBJ_SRC) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(TRIHAM_OBJ_SRC) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_TRIHAM) -lblas -llapack -lgfortran $(LIBRARY)

triham-ex-direct: $(TRIHAM_OBJ_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(TRIHAM_OBJ_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_TRIHAM) -lblas -llapack -lgfortran $(LIBRARY)
		
bdpc-ex: $(BDPC_OBJ_SRC) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(BDPC_OBJ_SRC) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_BDPC) -lblas -llapack -lgfortran $(LIBRARY)

bdpc-ex-direct: $(BDPC_MAIN_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(BDPC_MAIN_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_BDPC) -lblas -llapack  -lgfortran $(LIBRARY)

jcm-ex: $(JCM-MC-EX_SRC) $(QTM_OBJ_ZVODE_SRC)
	$(CC) $(JCM-MC-EX_SRC) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_JCM) -lblas -llapack -lgfortran $(LIBRARY)
	
jcm-ex-direct: $(JCM-MC-EX_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC)
	$(CC) $(JCM-MC-EX_SRC_DIRECT) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_JCM) -lblas -llapack -lgfortran $(LIBRARY)

#pumped-cavity-ex: $(PC_MAIN_SRC) $(QTM_OBJ_ZVODE_SRC)
#		$(CC) $(PC_MAIN_SRC) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME_PC) -lblas -llapack -lgfortran $(LIBRARY)

testlfsr113:
	$(MGXX) $(CFLAGS) -c rgen_lfsr113.cc
	$(MGXX) $(CFLAGS) -o testrgen_lfsr113 testrgen_lfsr113.cc rgen_lfsr113.o

misctest:
	$(MGXX) $(CFLAGS) -o misctest misctest.cc

.cc.o:
	$(CC) -c $(CFLAGS) $< -o $@

.f.o: 
	$(FC) -c $(FCFLAGS) $< -o $@

zvode_download:
	$(FGET) http://netlib.sandia.gov/linpack/zgbfa.f
	$(FGET) http://netlib.sandia.gov/linpack/zgbsl.f
	$(FGET) http://netlib.sandia.gov/linpack/zgefa.f
	$(FGET) http://netlib.sandia.gov/linpack/zgesl.f
	$(FGET) http://netlib.sandia.gov/ode/zvode.f
	
run:
	mpirun -n 2 qtm.exe

run4:
	mpirun -n 4 qtm.exe

partclean:
	rm -f unitary-mc-ex.o triham-mc-ex.o bdpc-ex.o jaynes-cummings-model-mc-ex.o $(OUTPUTNAME_U) $(OUTPUTNAME_TRIHAM) $(OUTPUTNAME_BDPC) $(OUTPUTNAME_JCM)
	
clean: 
	rm -f *.o $(OUTPUTNAME_U) $(OUTPUTNAME_TRIHAM) $(OUTPUTNAME_BDPC) $(OUTPUTNAME_JCM)

