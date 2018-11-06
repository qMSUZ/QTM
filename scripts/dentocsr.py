#! /usr/bin/python

#/***************************************************************************\
# *   Copyright (C) 2018 by Marek Sawerwain                         	    *
# *                                     <M.Sawerwain@gmail.com>             *
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
#\***************************************************************************/
 
import numpy as np
import scipy.sparse as sp

def print_csr_matrix(csrmat, outputname, fileoutname):
    f = open(fileoutname, 'w')
    i=0;
    for d in csrmat.data:
        f.write(outputname + "values["+str(i)+"]=make_simpleComplex("+str(d.real)+","+str(d.imag)+");\n");
        i=i+1;
    f.write("\n");
    
    i=0;
    for d in csrmat.indices:
        f.write(outputname + "col_ind["+str(i)+"]="+str(d)+";\n");
        i=i+1;
    f.write("\n");
    
    i=0;
    for d in csrmat.indptr:
        f.write(outputname + "row_ptr["+str(i)+"]="+str(d)+";\n");
        i=i+1;
    f.write("\n");

    f.close()


def print_alpha_zero(vardata, outputname, fileoutname):
    f = open(fileoutname, 'w')
    i=0;
    for v in vardata:
        d=v[0][0]
        f.write(outputname + "["+str(i)+"]=make_simpleComplex("+str(d.real)+","+str(d.imag)+");\n");
        i=i+1;
    f.write("\n");

    f.close()


def trival_ex1():
    denmat = np.matrix([
    [0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 2, 0, 0],
    [0, 0, 0, 3, 0],
    [0, 0, 0, 0, 4]]);
    
    
    denmat2 = np.matrix([    
    [10, 0, 0, 0, -2,  0],
    [ 3, 9, 0, 0,  0,  3],
    [ 0, 7, 8, 7,  0,  0], 
    [ 3, 0, 8, 7,  5,  0],
    [ 0, 8, 0, 9,  9, 13],
    [ 0, 4, 0, 0,  2, -1]]);
        
        
    print('matrix = ', denmat);
    print('matrix = ', denmat2);
    
    
    csrmat = sp.csr_matrix(denmat);
    csrmat2 = sp.csr_matrix(denmat2);
    
    print('csrmat data = ', csrmat.data);
    print('    indices = ', csrmat.indices);
    print('     indptr = ', csrmat.indptr);
    
    print('csrmat data = ', csrmat2.data);
    print('    indices = ', csrmat2.indices);
    print('     indptr = ', csrmat2.indptr);
    
    print("-----------------------------------------");
    
    #print_csr_matrix(csrmat, "col1.");

trival_ex1()