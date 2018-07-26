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
#include "qtm.h" 

const size_t COLLAPSE_OPERATORS = 3;
const size_t Ntrj = 300;
const size_t N = 100;
const size_t WAVEVECTOR_LEAD_DIM = 8*8*8;
const size_t WAVEVECTOR_LEAD_DIM_SQR = WAVEVECTOR_LEAD_DIM*WAVEVECTOR_LEAD_DIM;

uCSRMatrix<simpleComplex<double>, 448, 513, 448 > c_ops[ COLLAPSE_OPERATORS ];

uCSRMatrix<simpleComplex<double>, 448, 513, 448 > co0;
uCSRMatrix<simpleComplex<double>, 448, 513, 448 > co1;
uCSRMatrix<simpleComplex<double>, 448, 513, 448 > co2;


uCSRMatrix<simpleComplex<double>, 448, 513, 448 > expect_operator;
uCSRMatrix<simpleComplex<double>, 448, 513, 448 > expect_operator1;
uCSRMatrix<simpleComplex<double>, 448, 513, 448 > expect_operator2;

uCSRMatrix<simpleComplex<double>, 1197, 513, 1197> H;

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
	simpleComplex<double> moneimag;

	moneimag.re=0.0;
	moneimag.im=-1.0;

/*
	const int N0=8, N1=8, N2=8;
	double K=1.0;
	double gamma0=0.1, gamma1=0.1, gamma2=0.4;
	double alpha=sqrt(3);
	
	simpleComplex<double> unity = make_simpleComplex(0.0,1.0);
	
	uMatrix< simpleComplex<double>, N0 > d1;
	uMatrix< simpleComplex<double>, N1 > d2;
	uMatrix< simpleComplex<double>, N2 > d3;
	
	uMatrix< simpleComplex<double>, N0*N1*N2 > a0, C0, num0;
	uMatrix< simpleComplex<double>, N0*N1*N2 > a1, C1, num1;
	uMatrix< simpleComplex<double>, N0*N1*N2 > a2, C2, num2;
	
	uMatrix< simpleComplex<double>, N0*N1*N2 > H;
	
	destroy_operator(d1);
	eye_of_matrix(d2);
	eye_of_matrix(d3);
	a0=tensor(d1,d2,d3);

	eye_of_matrix(d1);
	destroy_operator(d2);
	eye_of_matrix(d3);
	a1=tensor(d1,d2,d3);

	eye_of_matrix(d1);
	destroy_operator(d2);
	eye_of_matrix(d3);
	a2=tensor(d1,d2,d3);

	num0=dagger(a0)*a0;
	num1=dagger(a1)*a1;
	num2=dagger(a2)*a2;

	
	C0=sqrt(2.0*gamma0)*a0;
	C1=sqrt(2.0*gamma1)*a1;
	C2=sqrt(2.0*gamma2)*a2;
*/
	
	/*
	vacuum=tensor(basis(N0,0),basis(N1,0),basis(N2,0))
	D=(alpha*a0.dag()-np.conj(alpha)*a0).expm()
	psi0=D*vacuum	
	*/
	
	//H=unity*K*(a0*dagger(a1)*dagger(a2)-dagger(a0)*a1*a2);

	/*
	Heff = (H - ((1.0j)/2.0) * (C0.dag()*C0 + C1.dag()*C1 + C2.dag()*C2))	
	
*/
	
	#include "data-alpha-triham.txt"

	#include "data-co0-triham.txt"	
	#include "data-co1-triham.txt"	
	#include "data-co2-triham.txt"

	#include "data-h-triham.txt"
	
	#include "data-e0-triham.txt"
	#include "data-e1-triham.txt"
	#include "data-e2-triham.txt"

	
	c_ops[0] = co0;
	c_ops[1] = co1;
	c_ops[2] = co2;

	for(i=0;i<H._values_size;i++)
	{
		H.values[i] = H.values[i] * moneimag;
	}
	
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