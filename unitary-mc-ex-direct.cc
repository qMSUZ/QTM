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

#define __USE_DENSE_EXPECT_OPERATORS
#define __USE_DENSE_COLLAPSE_OPERATORS
#define __USE_ADAMS

#include "complexnum.h"
#include "qtm.h" 

const size_t COLLAPSE_OPERATORS = 1;
const size_t Ntrj = 10;
const size_t N = 100;
const size_t WAVEVECTOR_LEAD_DIM = 2;
const size_t WAVEVECTOR_LEAD_DIM_SQR = 4;

uMatrix< simpleComplex<double>, WAVEVECTOR_LEAD_DIM > c_ops[ COLLAPSE_OPERATORS ];
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > collapse_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > expect_operator;
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > H;
simpleComplex<double> alpha[WAVEVECTOR_LEAD_DIM];

#include "qtm.cc"

extra_options opt;

int myfex_fnc_f1( long int *NEQ,
			double *T,
			dblcmplx *Y,
			dblcmplx *YDOT,
			dblcmplx *RPAR,
			long int *IPAR)
{
	simpleComplex<double> o0, o1, out0, out1;

	  o0.re=0.0;	  o0.im=0.0; o1.re=0.0;   o1.im=0.0;
	out0.re=0.0;	out0.im=0; out1.re=0.0; out1.im=0.0;

	o0 = Y[0] * H[0];
	o1 = Y[1] * H[1];

	out0 = o0 + o1;

	o0.re=0.0;   o0.im=0.0; o1.re=0.0;   o1.im=0.0;

	o0 = Y[0] * H[2];
	o1 = Y[1] * H[3];

	out1 = o0 + o1;

	YDOT[0] = out0;
	YDOT[1] = out1;

	return 0;
}

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
	// Heff = (H - ((ih)/2.0) * sum(C^{+}_n C_n))
	// i -- imaginary unity 
	H[0] = make_simpleComplex( -0.00125, 0.0);    H[1] = make_simpleComplex( 0.0, -0.62831853);
	H[2] = make_simpleComplex( 0.0, -0.62831853); H[3] = make_simpleComplex( -0.00125, 0.0);

	c_ops[0].rows=2;
	c_ops[0].cols=2;
	c_ops[0].m = collapse_operator;

	opt.type_output = OUTPUT_FILE_PYTHON_STYLE;
	//opt.type_output = OUTPUT_FILE;
	opt.state_of_trj_output = OUTPUT_STATE_OF_TRJ_FILE;
	opt.rnd_test_retry=10;
	opt.verbose_mode=1;
	opt.only_final_trj = 1;
	opt.ode_method = METADAMS;
	opt.tolerance = 1e-7;
	//opt.file_name = strdup("output-data.txt");
	opt.file_name = strdup("output-data-matplotfig.py");
	opt.fnc = &myfex_fnc_f1;


	r = mpi_main<N, Ntrj, WAVEVECTOR_LEAD_DIM, WAVEVECTOR_LEAD_DIM_SQR, COLLAPSE_OPERATORS>(argc, argv,
		0, 10, 
		1, 1, opt);

	return r;
}