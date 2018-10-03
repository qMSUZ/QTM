/***************************************************************************
 *   Copyright (C) 2017 -- 2018 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *							   and Joanna Wiśniewska                       *
 *                                         <jwisniewska@wat.edu.pl>        *
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

#include <cmath>
#include "complexnum.h"

/* 
   Algorithm of matrix exponentiation is based on a work:
   "New Scaling and Squaring Algorithm for the Matrix Exponential", by
   Awad H. Al-Mohy and Nicholas J. Higham, August 2009.
   
   Implementation is based on the code:   
	https://github.com/cran/expm/blob/master/src/matexp_MH09.c
	by  Drew Schmidt (2013-2014) and Martin Maechle (2014)
   from 
	https://github.com/cran/expm
   
   
*/

extern "C" int zgemm_( 	char * 	TRANSA,
	char * 	TRANSB,
	long int *  	M,
	long int *  	N,
	long int *  	K,
	dblcmplx *  	ALPHA,
	dblcmplx *  	A,
	long int *  	LDA,
	dblcmplx *  	B,
	long int *  	LDB,
	dblcmplx *  	BETA,
	dblcmplx *  	C,
	long int *  	LDC );
		
template <typename T>
void matprod( uMatrix<simpleComplex<T> > &A, uMatrix<simpleComplex<T> > &B, uMatrix<simpleComplex<T> > &C)
{
	long int n = A.rows;
	
    dblcmplx one = make_simpleComplex(1.0, 0.0);
	dblcmplx zero = make_simpleComplex(0.0, 0.0);
	
    char trans = 'N';
	
    zgemm_(&trans, &trans, &n, &n, &n, &one, &A.m[0], &n, &B.m[0], &n, &zero, &C.m[0], &n);
}

extern "C" int zlacpy_ ( 	char *  	UPLO,
	long int *  	M,
	long int *  	N,
	dblcmplx *  	A,
	long int *  	LDA,
	dblcmplx *  	B,
	long int *  	LDB ); 	

	
template <typename T>	
void matcopy(uMatrix< simpleComplex<T> > &A, uMatrix< simpleComplex<T> > &B)
{
	long int n = A.rows;
	
	char uplo = 'A';

	zlacpy_(&uplo, &n, &n, &A.m[0], &n, &B.m[0], &n);
}



template <typename T>	
T matnorm_1(uMatrix<simpleComplex<T> > &m)
{
	T norm = 0.0;
	for (int j=0; j<m.rows; j++)
	{
		T tmp = 0;
		for (int i=0; i<m.cols; i++)
			tmp += simpleComplexAbs(m.m[i + j*m.rows]);
		if (tmp > norm)
			norm = tmp;
	}
	return norm;
}



template <typename T>	
long int matexp_scale_factor(uMatrix< simpleComplex<T> > &m)
{
	const long int NTHETA = 5;
    const T theta[] = {1.5e-2, 2.5e-1, 9.5e-1, 2.1e0, 5.4e0};
    const T x_1 = matnorm_1<T>(m);

    for (int i=0; i < NTHETA; i++) {
		if (x_1 <= theta[i])
			return 0;
    }

    int i = (int) ceil(log2(x_1/theta[4]));
    return 1 << i;
}

template <typename T>	
static void matpow_by_squaring( uMatrix< simpleComplex<T> > &A, long int b, uMatrix< simpleComplex<T> > &P)
{
	if (b == 1) {
		matcopy<T>(A, P);
		return;
	}
    
	eye_of_matrix(P);
	
    if (b == 0)
		return;

	uMatrix< simpleComplex<T> > TMP(A.rows, A.cols);
	
	zero_matrix( TMP );
	
    while (b) {
		if (b&1) { // P := P A
			matprod<T>(P, A, TMP);
			matcopy<T>(TMP, P);
		}

		b >>= 1;

		matprod<T>(A, A, TMP);
		matcopy<T>(TMP, A);
    }
}

const double matexp_pade_coefs[14] =
{
  1.0,
  0.5,
  0.12,
  1.833333333333333333333e-2,
  1.992753623188405797101e-3,
  1.630434782608695652174e-4,
  1.035196687370600414079e-5,
  5.175983436853002070393e-7,
  2.043151356652500817261e-8,
  6.306022705717595115002e-10,
  1.483770048404140027059e-11,
  2.529153491597965955215e-13,
  2.810170546219962172461e-15,
  1.544049750670308885967e-17
};

#define SGNEXP(x,pow) (x==0?(pow==0?1:0):(x>0?1:(pow%2==0?1:(-1))))

template <typename T>	
void matexp_pade_fillmats(const long int i,
			  uMatrix< simpleComplex<T> > &N, 
			  uMatrix< simpleComplex<T> > &D,
			  uMatrix< simpleComplex<T> > &B, 
			  uMatrix< simpleComplex<T> > &C)
{
  const double tmp = matexp_pade_coefs[i];
  const long int sgn = SGNEXP(-1, i);

    for (int j=0; j < N.rows*N.cols; j++)
	{
		simpleComplex<T> t_j = C.m[j];
		B.m[j] = t_j;
		
		t_j = t_j * tmp;
	
		N.m[j] = N.m[j] + t_j;
		D.m[j] = D.m[j] + sgn*t_j;
    }
}

extern "C" int zgesv_ 	( 	long int *  	N,
		long int *  	NRHS,
		dblcmplx *  	A,
		long int *  	LDA,
		long int *  	IPIV,
		dblcmplx *  	B,
		long int *  	LDB,
		long int *  	INFO 	);

template <typename T>
static void matexp_pade(const long int p, uMatrix< simpleComplex<T> > &A, uMatrix< simpleComplex<T> > &N)
{
	long int n = A.rows;
	long int ipiv[A.rows];
	
    long int i, info = 0;

	uMatrix< simpleComplex<T> > B(A.rows, A.cols);
	zero_matrix(B);

	uMatrix< simpleComplex<T> > C(A.rows, A.cols);
	zero_matrix(C);
	matcopy<T>(A, C);
	
	uMatrix< simpleComplex<T> > D(A.rows, A.cols);
	zero_matrix(D);

    for (i=0; i<N.size; i++)
	{
		N.m[i].re = 0.0;
		N.m[i].im = 0.0;
		D.m[i].re = 0.0;
		D.m[i].im = 0.0;
    }

    i = 0;
    while (i < (N.size) )
	{
		N.m[i].re = 1.0;
		N.m[i].im = 0.0;
		D.m[i].re = 1.0;
		D.m[i].im = 0.0;

		i += n+1;
    }

    for (i=1; i<=p; i++)
    {
		if (i > 1)
			matprod<T>(A, B, C);

		matexp_pade_fillmats<T>(i, N, D, B, C);
    }

	

    zgesv_(&n, &n, &D.m[0], &n, &ipiv[0], &N.m[0], &n, &info);
}


extern "C" int zscal_	( 	long int *	 	N,
		dblcmplx *  	ZA,
		dblcmplx *  	ZX,
		long int *  	INCX 	);

template <typename T>
void exp_of_matrix(uMatrix< simpleComplex<T> > &x, const long int p, uMatrix< simpleComplex<T> > &ret)
{	
	long int n = x.rows; 
	long int m = matexp_scale_factor<T>(x);
	
	if (m == 0) {
		matexp_pade<T>(p, x, ret);
		return;
	}

	long int nn = n*n, one = 1;
	simpleComplex<T> tmp ;
	tmp.re = 1.0 / ((double) m);
	tmp.im = 0.0;

	zscal_(&nn, &tmp, &x.m[0], &one);

	matexp_pade<T>(p, x, ret);

	matcopy<T>(ret, x);

	matpow_by_squaring<T>(x, m, ret);
}
