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
#include <cmath>


#define __USE_SPARSE_CSR_EXPECT_OPERATORS 2
#define __USE_SPARSE_CSR_COLLAPSE_OPERATORS 2
//#define __USE_ADAMS METADAMS
#define __USE_BDF METBDF


#include "complexnum.h"
#include "complexnum.cc"

#include "qtm.h" 

#include "auxtools.cc"

const size_t COLLAPSE_OPERATORS = 3;
const size_t Ntrj = 100;
const size_t N = 100;

const size_t N0 = 8;
const size_t N1 = 8;
const size_t N2 = 8;

const size_t WAVEVECTOR_LEAD_DIM = N0*N1*N2;
const size_t WAVEVECTOR_LEAD_DIM_SQR = WAVEVECTOR_LEAD_DIM*WAVEVECTOR_LEAD_DIM;

uCSRMatrix< simpleComplex<double> > c_ops[3];

uCSRMatrix< simpleComplex<double> > co0;
uCSRMatrix< simpleComplex<double> > co1;
uCSRMatrix< simpleComplex<double> > co2;

uCSRMatrix< simpleComplex<double> > expect_operator;
uCSRMatrix< simpleComplex<double> > expect_operator1;
uCSRMatrix< simpleComplex<double> > expect_operator2;

uCSRMatrix< simpleComplex<double> > H;


//uCSRMatrix< simpleComplex<double> > co0(448, 513, 448);
//uCSRMatrix< simpleComplex<double> > co1(448, 513, 448);
//uCSRMatrix< simpleComplex<double> > co2(448, 513, 448);

//uCSRMatrix< simpleComplex<double> > expect_operator(448, 513, 448);
//uCSRMatrix< simpleComplex<double> > expect_operator1(448, 513, 448);
//uCSRMatrix< simpleComplex<double> > expect_operator2(448, 513, 448);

//uCSRMatrix< simpleComplex<double> > H(1197, 513, 1197);

simpleComplex<double> alpha[WAVEVECTOR_LEAD_DIM];

#include "qtm.cc"

extra_options opt;

int myfex_fnc_f1(	long int *NEQ,
			double *T,
			dblcmplx *Y,
			dblcmplx *YDOT,
			dblcmplx *RPAR,
			long int *IPAR)
{
    size_t i, j;

    for ( i=0; i < WAVEVECTOR_LEAD_DIM; i++)
    {
        YDOT[i].re = 0.0;
        YDOT[i].im = 0.0;
    }

    for ( i=0; i < WAVEVECTOR_LEAD_DIM; i++)
    {
        for ( j=H.row_ptr[i] ; j < H.row_ptr[i+1] ; j++)
        {
            YDOT[i]= YDOT[i] + (H.values[j] * Y[H.col_ind[j]]);
        }
    }

	return 0;
}


int prepare_matrices()
{

	int i;
	simpleComplex<double> m;
	simpleComplex<double> unity = make_simpleComplex(0.0,1.0);

	double K=1.0;
	double gamma0=0.1, gamma1=0.1, gamma2=0.4;
	double alphaval=sqrt(3);
	
	
	uMatrix< simpleComplex<double> > d1(N0);
	uMatrix< simpleComplex<double> > d2(N1);
	uMatrix< simpleComplex<double> > d3(N2);
	
	uMatrix< simpleComplex<double> > a0(N0*N1*N2), C0(N0*N1*N2), num0(N0*N1*N2);
	uMatrix< simpleComplex<double> > a1(N0*N1*N2), C1(N0*N1*N2), num1(N0*N1*N2);
	uMatrix< simpleComplex<double> > a2(N0*N1*N2), C2(N0*N1*N2), num2(N0*N1*N2);
	
	uMatrix< simpleComplex<double> > D(N0*N1*N2), Dtmp(N0*N1*N2);
	uMatrix< simpleComplex<double> > Hsys(N0*N1*N2), Heff(N0*N1*N2);
	
	simpleComplex<double> b0[N0];
	simpleComplex<double> b1[N1];
	simpleComplex<double> b2[N2];

	simpleComplex<double> bt0[N0*N1];
	simpleComplex<double> vacuum[N0*N1*N2];

	zero_matrix(D);
	zero_matrix(Dtmp);
	zero_matrix(Hsys);
	
	zero_matrix(a0);
	zero_matrix(a1);
	zero_matrix(a2);

	zero_matrix(C0);
	zero_matrix(C1);
	zero_matrix(C2);

	zero_matrix(num0);
	zero_matrix(num1);
	zero_matrix(num2);

	
	zero_matrix(d1); zero_matrix(d2); zero_matrix(d3);
	destroy_operator(d1);
	eye_of_matrix(d2);
	eye_of_matrix(d3);
	a0=tensor(d1,d2,d3);

	zero_matrix(d1); zero_matrix(d2); zero_matrix(d3);
	eye_of_matrix(d1);
	destroy_operator(d2);
	eye_of_matrix(d3);
	a1=tensor(d1,d2,d3);

	zero_matrix(d1); zero_matrix(d2); zero_matrix(d3);
	eye_of_matrix(d1);
	eye_of_matrix(d2);
	destroy_operator(d3);
	a2=tensor(d1,d2,d3);

	num0=dagger(a0)*a0;
	num1=dagger(a1)*a1;
	num2=dagger(a2)*a2;
	
	C0=sqrt(2.0*gamma0)*a0;
	C1=sqrt(2.0*gamma1)*a1;
	C2=sqrt(2.0*gamma2)*a2;
	
	co0 = uMatrix_to_uCSRMatrix(C0);	
	co1 = uMatrix_to_uCSRMatrix(C1);	
	co2 = uMatrix_to_uCSRMatrix(C2);	

	//dump_uCSRMatrix_to_file<double>(co0, "co0", "dump-triham-co0.txt");
	//dump_uCSRMatrix_to_file<double>(co1, "co1", "dump-triham-co1.txt");
	//dump_uCSRMatrix_to_file<double>(co2, "co2", "dump-triham-co2.txt");

	expect_operator  = uMatrix_to_uCSRMatrix( num0 );
	expect_operator1 = uMatrix_to_uCSRMatrix( num1 );
	expect_operator2 = uMatrix_to_uCSRMatrix( num2 );

	//dump_uCSRMatrix_to_file<double>(expect_operator,  "eo0", "dump-triham-eo0.txt");
	//dump_uCSRMatrix_to_file<double>(expect_operator1, "eo1", "dump-triham-eo1.txt");
	//dump_uCSRMatrix_to_file<double>(expect_operator2, "eo2", "dump-triham-eo2.txt");
	
	
	std_base_state<double, N0>(&b0[0], 0);
	std_base_state<double, N1>(&b1[0], 0);
	std_base_state<double, N2>(&b2[0], 0);
	
	tensor<double, N0, N1>(&b0[0], &b1[0], &bt0[0]);
	tensor<double, N0*N1, N2>(&bt0[0], &b2[0], &vacuum[0]);
	
	Dtmp = alphaval * dagger(a0) - alphaval*a0;
	exp_of_matrix(Dtmp, 10, D);
	
	/*
	vacuum=tensor(basis(N0,0),basis(N1,0),basis(N2,0))
	D=(alpha*a0.dag()-np.conj(alpha)*a0).expm()
	psi0=D*vacuum	
	*/
	
	mul_mat_vec(D, &vacuum[0], &alpha[0]);
	
	//dump_table_to_file<double, WAVEVECTOR_LEAD_DIM>(&alpha[0], "alpha", "dump-alpha-triham.txt");
	
	
	//H=unity*K*(a0*dagger(a1)*dagger(a2)-dagger(a0)*a1*a2);
	//H=1j*K*(a0*a1.dag()*a2.dag()-a0.dag()*a1*a2)
	//H=1.0im*K*(a0*dagger(a1)*dagger(a2) - dagger(a0)*a1*a2)

	Hsys = unity*K*(a0*dagger(a1)*dagger(a2)-dagger(a0)*a1*a2);
	
	//Heff = (H - ((1.0j)/2.0) * (C0.dag()*C0 + C1.dag()*C1 + C2.dag()*C2))	
	
	m.re=0;
	m.im=0.5;
	
	Heff = Hsys - m*(dagger(C0)*C0 + dagger(C1)*C1 + dagger(C2)*C2);
	H = uMatrix_to_uCSRMatrix(Heff);
	
	
//	#include "data-alpha-triham.txt"

//	#include "data-co0-triham.txt"	
//	#include "data-co1-triham.txt"	
//	#include "data-co2-triham.txt"

//	#include "data-h-triham.txt"
	
//	#include "data-e0-triham.txt"
//	#include "data-e1-triham.txt"
//	#include "data-e2-triham.txt"

	
	c_ops[0] = co0;
	c_ops[1] = co1;
	c_ops[2] = co2;

	
	m.re=0.0;
	m.im=-1.0;
	
	for(i=0;i<H._values_size;i++)
	{
		H.values[i] = H.values[i] * m;
	}
	
	//dump_uCSRMatrix_to_file<double>(H, "h", "dump-triham-h.txt");

	
	return 0;
}

int main(int argc, char *argv[])
{
	int r = 0;

	prepare_matrices();
	
	//opt.type_output = OUTPUT_FILE;
	opt.type_output = OUTPUT_FILE_PYTHON_STYLE;
	//opt.verbose_mode = 2;
	opt.verbose_mode = 0;
	opt.only_final_trj = 1;
	opt.ode_method = __USE_BDF;
	//opt.ode_method = __USE_ADAMS;
	opt.tolerance = 1e-7;
	//opt.file_name = strdup("output-data.txt");
	opt.file_name = strdup("output-data-matplotfig.py");
	opt.fnc = &myfex_fnc_f1;
	
	r = mpi_main<N, Ntrj, 
		WAVEVECTOR_LEAD_DIM, WAVEVECTOR_LEAD_DIM_SQR, COLLAPSE_OPERATORS>(argc, argv,
		0.0, 4.0,
		__USE_SPARSE_CSR_COLLAPSE_OPERATORS, __USE_SPARSE_CSR_EXPECT_OPERATORS, opt);
	
	return r;
}
