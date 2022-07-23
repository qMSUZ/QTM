#/***************************************************************************
# *   Copyright (C) 2019, 2021, 2022 by Marek Sawerwain                     *
# *                                         <M.Sawerwain@gmail.com>         *
# *                                         <M.Sawerwain@issi.uz.zgora.pl   *
# *                                                                         *
# *   Copyright (C) 2005 -- 2012 by Marek Sawerwain                         *
# *                                         <M.Sawerwain@gmail.com>         *
# *                                                                         *
# *   Part of the Quantum Computing Simulator:                              *
# *   https://github.com/qMSUZ/QCS                                          *
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

CC=gcc
CXX=g++

ifeq ($(DEBUG),1) 
CFLAGS=-g -I./include 
else
CFLAGS=-I./include 
endif

ifeq ($(PYTHON),1)
CFLAGS+=-DPYTHON_SCRIPT $(PYTHON_CFLAGS) -I/opt/intel/oneapi/intelpython/python3.9/include/python3.9 -fPIC
endif


SWIGCMD=swig
SWIGOPT=-DPYTHON_SCRIPT -python -I./include
PYTHON_LIB=`pkg-config python3 --libs`
PYTHON_CFLAGS=`pkg-config python3 --cflags`
QCS_PYTHON_OUT = _qcs.so

C_MAIN_SOURCES = src/qcs_misc.c \
	src/qcs_hash_table.c \
	src/qcs_complex.c \
	src/qcs_rand.c \
	src/qcs_matrix_and_vector.c \
	src/qcs_qubit.c \
	src/qcs_qubit_gates.c \
	src/qcs_gates.c \
	src/qcs_quantum_register.c \
	src/qcs_info.c
	
C_OBJ_MAIN_SOURCES=$(C_MAIN_SOURCES:.c=.o)
	

all: microex1


qcs_wrap.o: src/qcs_wrap.c
	$(CC) $(CFLAGS) -c $<

src/qcs_wrap.c: src/qcs.i
	$(SWIGCMD) $(SWIGOPT) $<

library: $(C_OBJ_MAIN_SOURCES) $(CC_OBJ_MAIN_SOURCES)
	ar rcs libqcs.a $(C_OBJ_MAIN_SOURCES) $(CC_OBJ_MAIN_SOURCES)
	@echo "+--------------------------------------------------------------------------------+"
	@echo "Static library has been created."

python_port: library qcs_wrap.o
	$(CC) $(LIB_OPT) -shared qcs_wrap.o libqcs.a -o $(QCS_PYTHON_OUT) $(PYTHON_LIB) -llapack -lblas -lgfortran -lm
	cp src/qcs.py qcs.py

ex-spectral-decomposition-test: examples_ansi_c/ex-spectral-decomposition-test.c library
	$(CC) -o ex-spectral-decomposition-test examples_ansi_c/ex-spectral-decomposition-test.c -I./include -L. $(CFLAGS) $(PYTHON_LIB)  -lqcs -lpython3 -llapack -lblas -lgfortran -lm

ex-rand-test: examples_ansi_c/ex-rand-test.c library
	$(CC) -o ex-rand-test examples_ansi_c/ex-rand-test.c -I./include -L. $(CFLAGS) $(PYTHON_LIB)  -lqcs -lpython3.9 -llapack -lblas -lgfortran -lm

ex-matrix-and-vector-test: examples_ansi_c/ex-matrix-and-vector-test.c library
	$(CC) -o ex-matrix-and-vector-test examples_ansi_c/ex-matrix-and-vector-test.c -I./include -L. $(CFLAGS) -lqcs -llapack -lblas -lgfortran -lm


microex1: examples_ansi_c/microexamples.c library
	$(CC) -D__microEX1__ -o microex1 examples_ansi_c/microexamples.c -I./include -L. $(CFLAGS) -lqcs -llapack -lblas -lgfortran -lm
#	$(CC) -o ex1 examples_ansi_c/ex1.c -I./include -L. $(CFLAGS) $(PYTHON_LIB)  -lqcs -lpython3.8 -llapack -lblas -lgfortran -lm

microex2: examples_ansi_c/microexamples.c library
	$(CC) -D__microEX2__ -o microex2 examples_ansi_c/microexamples.c -I./include -L. $(CFLAGS) -lqcs -llapack -lblas -lgfortran -lm
#	$(CC) -o ex1 examples_ansi_c/ex1.c -I./include -L. $(CFLAGS) $(PYTHON_LIB)  -lqcs -lpython3.8 -llapack -lblas -lgfortran -lm

microex3: examples_ansi_c/microexamples.c library
	$(CC) -D__microEX3__ -o microex3 examples_ansi_c/microexamples.c -I./include -L. $(CFLAGS) -lqcs -llapack -lblas -lgfortran -lm
#	$(CC) -o ex1 examples_ansi_c/ex1.c -I./include -L. $(CFLAGS) $(PYTHON_LIB)  -lqcs -lpython3.8 -llapack -lblas -lgfortran -lm

microex4: examples_ansi_c/microexamples.c library
	$(CC) -D__microEX4__ -o microex4 examples_ansi_c/microexamples.c -I./include -L. $(CFLAGS) -lqcs -llapack -lblas -lgfortran -lm
#	$(CC) -o ex1 examples_ansi_c/ex1.c -I./include -L. $(CFLAGS) $(PYTHON_LIB)  -lqcs -lpython3.8 -llapack -lblas -lgfortran -lm


clean:
	rm -f *.o src/qcs_wrap.c examples_ansi_c/*.o src/*.o libqcs.a qcs.py qcs_wrap.c qcs_warp.o _qcs.so *.pyc \
		microex1 microex2 microex3 microex4 \
		ex-rand-test ex-spectral-decomposition-test ex-matrix-and-vector-test

help:
	@echo ".:             The Quantum Computing Simulator -- make build system             :."
	@echo "+--------------------------------------------------------------------------------+"
	@echo "+ Basic examples:                                                                +"
	@echo "+      \"make library\" for building static library using default configuration.   +"
	@echo "+      \"make PYTHON=1 python_port\" for building the module for Python language   +"
	@echo "+                                                                                +"
	@echo "+ Default action:                                                                +"
	@echo "+      \"microex1\" for building the micro example 1                               +"
	@echo "+                                                                                +"
	@echo "+--------------------------------------------------------------------------------+"
