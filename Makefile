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
#CC      = mpic++
CC      = g++
RUN     = mpirun
CFLAGS  = -I. -fpermissive
FCFLAGS = -I.

OUTPUTNAME = qtm.exe

#QTM_MAIN_SRC = qtm.cc Uniform.cc RNG.cc MLCG.cc ACG.cc
QTM_MAIN_SRC = unitary-mc-ex.cc rgen_lfsr113.cc
ZVODE_SRC = zvode.f zgesl.f zgefa.f zgbsl.f zgbfa.f

QTM_OBJ_MAIN_SRC=$(QTM_MAIN_SRC:.cc=.o)
QTM_OBJ_ZVODE_SRC=$(ZVODE_SRC:.f=.o)

all: compile_task

compile_task: $(QTM_OBJ_MAIN_SRC) $(QTM_OBJ_ZVODE_SRC)
		$(CC) $(QTM_OBJ_MAIN_SRC) $(QTM_OBJ_ZVODE_SRC) -o $(OUTPUTNAME) -lblas -lgfortran -lmsmpi  -lwsock32

testlfsr113:
	$(MGXX) $(CFLAGS) -c rgen_lfsr113.cc
	$(MGXX) $(CFLAGS) -o testrgen_lfsr113 testrgen_lfsr113.cc rgen_lfsr113.o

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
 
clean:
	rm -f *.o *.exe

