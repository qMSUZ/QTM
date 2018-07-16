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
#include "qtm.h" 

const size_t Ntrj = 1;
const size_t N = 600;
const size_t WAVEVECTOR_LEAD_DIM = 80;
const size_t WAVEVECTOR_LEAD_DIM_SQR = 80*80;

uMatrix< simpleComplex<double>, WAVEVECTOR_LEAD_DIM > c_ops[ 1 ];
uVector< simpleComplex<double>, WAVEVECTOR_LEAD_DIM_SQR > collapse_operator;

// template <typename T, size_t _S_ValueSize, size_t _S_RowPtr, size_t _S_ColInd>

uCSRMatrix<simpleComplex<double>, 40, 81, 40 > expect_operator;
uCSRMatrix<simpleComplex<double>, 158, 81, 158> H;

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
	
	#include "data-h-jcm.txt"
	#include "data-eops-jcm.txt"	
	#include "data-alpha-jcm.txt"
	
	return 0;
}

int main(int argc, char *argv[])
{
	int i, r = 0;
	simpleComplex<double> moneimag;
	
	moneimag.re=0.0;
	moneimag.im=-1.0;

	prepare_matrices();
	
	opt.type_output = OUTPUT_FILE_PYTHON_STYLE;
	opt.only_final_trj = 1;
	//opt.ode_method = __USE_ADAMS;
	opt.ode_method = __USE_BDF;
	opt.tolerance = 1e-12;
	opt.file_name = strdup("output-data-matplotfig.py");
	opt.fnc = &myfex_fnc_f1;
	

	for(i=0;i<H._values_size;i++)
	{
		H.values[i] = moneimag * H.values[i];
	}
	
	r = mpi_main<N, Ntrj, WAVEVECTOR_LEAD_DIM, WAVEVECTOR_LEAD_DIM_SQR, 0>(argc, argv, 1, 
	0.0, 35.0, 
	0, __USE_SPARSE_CSR_EXPECT_OPERATORS, opt);


	
	return r;
}