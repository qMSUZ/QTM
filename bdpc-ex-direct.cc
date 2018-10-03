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

#define __USE_DENSE_EXPECT_OPERATORS 1
#define __USE_DENSE_COLLAPSE_OPERATORS 1
//#define __USE_ADAMS METADAMS
#define __USE_BDF METBDF

#include "complexnum.h"
#include "qtm.h" 

const size_t COLLAPSE_OPERATORS = 2;
const size_t Ntrj = 200;
const size_t N = 100;
const size_t WAVEVECTOR_LEAD_DIM = 5;
const size_t WAVEVECTOR_LEAD_DIM_SQR = 25;

uMatrix< simpleComplex<double> > c_ops[ COLLAPSE_OPERATORS ] = { {WAVEVECTOR_LEAD_DIM}, {WAVEVECTOR_LEAD_DIM} };
uVector< simpleComplex<double> > co0(WAVEVECTOR_LEAD_DIM_SQR), co1(WAVEVECTOR_LEAD_DIM_SQR);
uVector< simpleComplex<double> > expect_operator(WAVEVECTOR_LEAD_DIM_SQR);
uVector< simpleComplex<double> > H(WAVEVECTOR_LEAD_DIM_SQR);
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
    // YDOT = H * Y;

    size_t i, k;

	for( i=0; i<WAVEVECTOR_LEAD_DIM ; i++)
    {
        YDOT[i].re = 0.0;
        YDOT[i].im = 0.0;
		for( k=0; k <  WAVEVECTOR_LEAD_DIM; k++)
        {
            YDOT[i] = YDOT[i] + H[ (i * WAVEVECTOR_LEAD_DIM) + k] * Y[k];
        } // for( k=0; k <  WAVEVECTOR_LEAD_DIM; k++)
    } // for( i=0; i<WAVEVECTOR_LEAD_DIM ; i++)

	return 0;
}



int main(int argc, char *argv[])
{
	int i, r = 0;
	
	simpleComplex<double> moneimag;
	
	moneimag.re=0.0;
	moneimag.im=-1.0;
	
	zero_vector(co0);
	zero_vector(co1);

	co0[ 1] = make_simpleComplex( 2.87059403, 0.0);
	co0[ 7] = make_simpleComplex( 4.05963301, 0.0);
	co0[13] = make_simpleComplex( 4.97201471, 0.0);
	co0[19] = make_simpleComplex( 5.74118806, 0.0);	

	co1[ 5] = make_simpleComplex( 0.69883624, 0.0);
	co1[11] = make_simpleComplex( 0.98830369, 0.0);
	co1[17] = make_simpleComplex( 1.21041988, 0.0);
	co1[23] = make_simpleComplex( 1.39767248, 0.0);	

	
	zero_vector( expect_operator );

    expect_operator[ 0] = make_simpleComplex( 0.0, 0.0);    
	expect_operator[ 6] = make_simpleComplex( 1.0, 0.0);
	expect_operator[12] = make_simpleComplex( 2.0, 0.0);
	expect_operator[18] = make_simpleComplex( 3.0, 0.0);
	expect_operator[24] = make_simpleComplex( 4.0, 0.0);
	
    alpha[0] = make_simpleComplex( 0.0, 0.0);
    alpha[1] = make_simpleComplex( 1.0, 0.0);
    alpha[2] = make_simpleComplex( 0.0, 0.0);
    alpha[3] = make_simpleComplex( 0.0, 0.0);
    alpha[4] = make_simpleComplex( 0.0, 0.0);
	
	
	// effective Hamiltonian
	// Heff = (H - ((ih)/2.0) * sum(C^{+}_n C_n))
	// i -- imaginary unity
	
	zero_vector( H );
    H[ 0] = make_simpleComplex( 0.0, -0.24418604651162792);    
	H[ 6] = make_simpleComplex( 1.0, -4.6085271317829459);
	H[12] = make_simpleComplex( 2.0, -8.9728682170542644);
	H[18] = make_simpleComplex( 3.0, -13.337209302325581 );
	H[24] = make_simpleComplex( 4.0, -16.480620155038761);
	
	for(i=0;i<WAVEVECTOR_LEAD_DIM_SQR;i++)
	{
		H[i] = moneimag * H[i];
	}
	
	
	//c_ops[0].rows=5;
    //c_ops[0].cols=5;
    //c_ops[0].m = co0;
	c_ops[0].recreate(5,5, co0);

	//c_ops[1].rows=5;
    //c_ops[1].cols=5;
    //c_ops[1].m = co1;
	c_ops[1].recreate(5,5, co1);
	
	//opt.type_output = OUTPUT_FILE;
	opt.type_output = OUTPUT_FILE_PYTHON_STYLE;
	opt.state_of_trj_output = 0;
	opt.rnd_test_retry=10;
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
		0.0, 0.8,
		__USE_DENSE_COLLAPSE_OPERATORS, __USE_DENSE_EXPECT_OPERATORS, opt);

	
	return r;
}