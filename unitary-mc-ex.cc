/***************************************************************************
 *   Copyright (C) 2017 -- 2018 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *                                                                         *
 *   Part of the Quantum Trajectory Method:                                *
 *   https://github.com/qMSUZ/QTM                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <cstdio>

#include "complexnum.h"
#include "qtm.h" 

const size_t Ntrj = 150;
const size_t N = 120;
const size_t WAVEVECTOR_LEAD_DIM = 2;
const size_t WAVEVECTOR_LEAD_DIM_SQR = 4;

uMatrix< simpleComplex<double>, WAVEVECTOR_LEAD_DIM > c_ops[ 1 ];
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > collapse_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > expect_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > H;
simpleComplex<double> alpha[WAVEVECTOR_LEAD_DIM];

#include "qtm.cc"


int main(int argc, char *argv[])
{
	int r = 0;
	collapse_operator[0] = make_simpleComplex( 0.0, 0.0 );  collapse_operator[1] = make_simpleComplex( 0.05, 0.0 );
    collapse_operator[2] = make_simpleComplex( 0.05, 0.0 ); collapse_operator[3] = make_simpleComplex( 0.0, 0.0 );

    expect_operator[0] = make_simpleComplex( 1.0, 0.0); expect_operator[1] = make_simpleComplex( 0.0, 0.0);
    expect_operator[2] = make_simpleComplex( 0.0, 0.0); expect_operator[3] = make_simpleComplex(-1.0, 0.0);

    alpha[0] = make_simpleComplex( 1.0, 0.0);
    alpha[1] = make_simpleComplex( 0.0, 0.0);
	
	
	// effective Hamiltonian
	// Heff = (H - ((ih)/2.0) * sum(C^{+}_n C_n)) * -1.0i
	// i -- imaginary unity
    H[0] = make_simpleComplex( -0.00125, 0.0);    H[1] = make_simpleComplex( 0.0, -0.62831853);
    H[2] = make_simpleComplex( 0.0, -0.62831853); H[3] = make_simpleComplex( -0.00125, 0.0);
	
	c_ops[0].rows=2;
    c_ops[0].cols=2;
    c_ops[0].m = collapse_operator;
	
	r = mpi_main<N, Ntrj, WAVEVECTOR_LEAD_DIM, WAVEVECTOR_LEAD_DIM_SQR, 1>(argc, argv, 1, 0, 10, 1, 1);


	
	return r;
}