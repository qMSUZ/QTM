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

#define __USE_SPARSE_CSR_EXPECT_OPERATORS 2
//#define __USE_DENSE_COLLAPSE_OPERATORS
//#define __USE_ADAMS METADAMS
#define __USE_BDF METBDF

#include "complexnum.h"
#include "complexnum.cc"
#include "qtm.h" 

const size_t Ntrj = 1;
const size_t N = 600;
const size_t WAVEVECTOR_LEAD_DIM = 80;
const size_t WAVEVECTOR_LEAD_DIM_SQR = WAVEVECTOR_LEAD_DIM*WAVEVECTOR_LEAD_DIM;

uMatrix< simpleComplex<double> > c_ops[ 1 ]  = { {WAVEVECTOR_LEAD_DIM} };
uVector< simpleComplex<double> > collapse_operator;

uCSRMatrix< simpleComplex<double> > expect_operator;
uCSRMatrix< simpleComplex<double> > H;

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
	simpleComplex<double> m, mhalf, alphaval;

	uMatrix< simpleComplex<double> > a(WAVEVECTOR_LEAD_DIM), sigmaminus(WAVEVECTOR_LEAD_DIM), Hsys(WAVEVECTOR_LEAD_DIM), eops(WAVEVECTOR_LEAD_DIM);
	
	uMatrix< simpleComplex<double> > d(WAVEVECTOR_LEAD_DIM/2), big_id(WAVEVECTOR_LEAD_DIM/2);
	uMatrix< simpleComplex<double> > small_id(2), sigmam(2);
	
	simpleComplex<double> c[WAVEVECTOR_LEAD_DIM/2];
	simpleComplex<double> b[2];
	
	double g = 1.0;
	double delta = -0.1;

	m.re=0.0;
	m.im=-1.0;
	
	mhalf.re=0.0;
	mhalf.im=-0.5f;
	
	alphaval = make_simpleComplex(4.0, 0.0);
	
	zero_matrix(a);
	zero_matrix(sigmaminus);
	zero_matrix(Hsys);
	zero_matrix(eops);
	
	zero_matrix(d);
	zero_matrix(big_id);
	zero_matrix(small_id);
	zero_matrix(sigmam);
	
	destroy_operator( d );
	
	eye_of_matrix(small_id);
	eye_of_matrix(big_id);
	sigma_m_matrix(sigmam);
	
	a = tensor(d, small_id);
	sigmaminus = tensor(big_id, sigmam);
		
	Hsys = (delta * dagger(a) * a) + (g * (dagger(a) * sigmaminus  + a * dagger(sigmaminus )));
	eops = dagger(sigmaminus) * sigmaminus;
	
	coherent<double, WAVEVECTOR_LEAD_DIM/2>( &c[0], alphaval);
	std_base_state<double, 2>(&b[0], 1);
	
	tensor<double, WAVEVECTOR_LEAD_DIM/2, 2>( &c[0], &b[0], &alpha[0]);

	H = uMatrix_to_uCSRMatrix(Hsys);	
	expect_operator = uMatrix_to_uCSRMatrix(eops);
	
	for(i=0;i<H._values_size;i++)
	{
		H.values[i] = m * H.values[i];
	}

	return 0;
}


int main(int argc, char *argv[])
{
	int r = 0;

	
	prepare_matrices();
	
	opt.type_output = OUTPUT_FILE_PYTHON_STYLE;
	opt.only_final_trj = 1;
	//opt.ode_method = __USE_ADAMS;
	opt.ode_method = __USE_BDF;
	opt.tolerance = 1e-12;
	opt.file_name = strdup("output-data-matplotfig.py");
	opt.fnc = &myfex_fnc_f1;
	
	
	r = mpi_main<N, Ntrj, WAVEVECTOR_LEAD_DIM, WAVEVECTOR_LEAD_DIM_SQR, 0>(argc, argv,
	0.0, 35.0, 
	0, __USE_SPARSE_CSR_EXPECT_OPERATORS, opt);
	
	return r;
}