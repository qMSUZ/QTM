#ifndef __complexnum_h__
#define __complexnum_h__

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

#include <iostream>
#include <iomanip>

#include <vector>
#include <cmath>


using namespace std;

template <typename T>
struct simpleComplex {
	T re;
	T im;
};

typedef simpleComplex<double> dblcmplx;

template <typename T>
struct simpleComplex<T> make_simpleComplex (T r, T i);


/*
template <typename T>
struct simpleCmplxMatrix {
    unsigned int rows, cols;
    std::vector< simpleComplex<T> > m;
};
*/


template <typename T, size_t v_size>
struct uVector {

    unsigned int size;

    uVector() { size = v_size; };

    T m[ v_size ];

    T& operator[](const size_t idx) { return m[idx]; };
    const T& operator[](const size_t idx) const { return m[idx]; };

};

template <typename T, size_t SIZE>
struct uMatrix {

    unsigned int _size,  rows, cols;

    uMatrix() {
        _size = SIZE * SIZE;
         rows = SIZE;
         cols = SIZE;
    }

    uVector< T, SIZE * SIZE > m;


    // get size
    inline unsigned int size() const
    { return _size; }

    // function operator() for indexing:  for reading (r)-rows, (c)-cols
    inline const T & operator()(const unsigned int r, const unsigned int c) const
    { return m[r*cols + c]; }

    // function operator() for indexing:  for writing (r)-rows, (c)-cols
    inline T & operator()(const unsigned int r, const unsigned int c)
    { return m[r*cols + c]; }
};

template <typename T, size_t _S_ValueSize, size_t _S_RowPtr, size_t _S_ColInd>
struct uCSRMatrix {

    unsigned int _values_size,  _row_ptr, _col_ind;

    uCSRMatrix()
    {
        _values_size = _S_ValueSize;
        _row_ptr = _S_RowPtr;
        _col_ind = _S_ColInd;
    }

    uVector< T, _S_ValueSize > values;

    uVector< unsigned int, _S_RowPtr> row_ptr;
    uVector< unsigned int, _S_ColInd> col_ind;

};


template <typename T, size_t _S_ValueSize, size_t _S_RowPtr, size_t _S_ColInd, size_t SIZE>
struct uCSRMatrix<T, _S_ValueSize, _S_RowPtr, _S_ColInd> make_uCSRMatrix(  struct uMatrix< T, SIZE> m )
{
    struct uCSRMatrix<T, _S_ValueSize, _S_RowPtr, _S_ColInd> csr_m;

    int i, idx, rowptridx, row, lastrow, lastcol;

    idx=0;
    rowptridx=1;
    row=0;
    lastrow=0;
    lastcol=-1;


    csr_m.row_ptr[0]=0;
    for(i=0;i<m.rows * m.cols;i++)
    {

        if (m.m[i].re!=0.0)
        {
            csr_m.values[idx] = m.m[i];
            csr_m.col_ind[idx] = i % m.cols;

            row = i / m.cols; 
            if( row != lastrow )
            {
                    lastrow = row;
                    csr_m.row_ptr[rowptridx]=idx;
                    rowptridx++;
            }
            
            idx++;

        }
    }
    csr_m.row_ptr[rowptridx]=idx;

    csr_m._values_size = idx;
    csr_m._col_ind = idx;
    csr_m._row_ptr = SIZE+1;

    return csr_m;
}

template<typename T, size_t v_size, size_t _S_ValueSize, size_t _S_RowPtr, size_t _S_ColInd>
uVector< T, v_size> mulCSRMatByuVec(uCSRMatrix<T, _S_ValueSize, _S_RowPtr, _S_ColInd> m, uVector< T, v_size> x)
{
    size_t i, j;


    uVector< T, v_size>  y;

    for ( i=0; i < v_size; i++)
    {
        y[i] = make_simpleComplex(0.0,0.0);
    }

    for ( i=0; i < v_size; i++)
    {
        for ( j=m.row_ptr[i] ; j < m.row_ptr[i+1] ; j++)
        {
            y[i]= y[i] + (m.values[j] * x[m.col_ind[j]]);
        }
    }

    return y;
}

template <typename T, size_t SIZE>
void exp_of_matrix(uMatrix< simpleComplex<T>, SIZE> &x, const long int p, uMatrix< simpleComplex<T>, SIZE> &ret);

#if 0
template <typename T>
struct simpleCmplxMatrix {

    unsigned int rows, cols;

    std::vector< simpleComplex<T> > m;

    //constructors

    inline simpleCmplxMatrix() {};

    inline simpleCmplxMatrix(const unsigned int Rows, const unsigned int Cols)
      : rows(Rows), cols(Cols) { m.resize(rows*cols); }

    // get size
    inline unsigned int size() const
    { return m.size(); }

    // function operator() for indexing:  for reading (r)-rows, (c)-cols
    inline const simpleComplex<T> & operator()(const unsigned int r, const unsigned int c) const
    { return m[r*cols + c]; }

    // function operator() for indexing:  for writing (r)-rows, (c)-cols
    inline simpleComplex<T> & operator()(const unsigned int r, const unsigned int c)
    { return m[r*cols + c]; }
};
#endif

template <typename T>
struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b);

template <typename T>
struct simpleComplex<T> operator+(const T &a, const struct simpleComplex<T> &b);

template <typename T>
struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const T &b);



template <typename T>
struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b);

template <typename T>
struct simpleComplex<T> operator*(const T &a, const struct simpleComplex<T> &b);

template <typename T>
struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const T &b);


/*template <typename T>
const std::vector< struct simpleComplex<T> > operator*(const T &a, const std::vector< struct simpleComplex<T> > &v);
*/

template <typename T>
struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b);

template <typename T>
struct simpleComplex<T> operator/(const T &_a, const struct simpleComplex<T> &b);

template <typename T>
struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const T &_b);



double simpleComplexAbs (const struct simpleComplex<double> &a);
double simpleComplexMod (const struct simpleComplex<double> &a);
float simpleComplexMod (const struct simpleComplex<float> &a);



template <typename T>
struct simpleComplex<T> make_simpleComplex (T r, T i);

template <typename T>
struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a.re + b.re;
	t.im = a.im + b.im;

	return t;
}

template <typename T>
struct simpleComplex<T> operator-(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a.re - b.re;
	t.im = a.im - b.im;

	return t;
}

template <typename T>
struct simpleComplex<T> operator+(const T &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a + b.re;
	t.im = b.im;

	return t;
}

template <typename T>
struct simpleComplex<T> operator-(const T &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a - b.re;
	t.im = b.im;

	return t;
}

template <typename T>
struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = a.re + b;
	t.im = a.im ;

	return t;
}

template <typename T>
struct simpleComplex<T> operator-(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = a.re - b;
	t.im = a.im ;

	return t;
}

template <typename T>
struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = (a.re * b.re) - (a.im * b.im);
	t.im = (a.re * b.im) + (a.im * b.re);

	return t;
}

template <typename T>
struct simpleComplex<T> operator*(const T &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = (a * b.re);
	t.im = (a * b.im);

	return t;
}

template <typename T>
struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = (a.re * b);
	t.im = (a.im * b);

	return t;
}

template <typename T>
struct simpleComplex<T> operator*(const long int &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = ((T)a * b.re);
	t.im = ((T)a * b.im);

	return t;
}

template <typename T>
struct simpleComplex<T> operator*(const struct simpleComplex<T> &a, const long int &b)
{
	struct simpleComplex<T> t;

	t.re = (a.re * (T)b);
	t.im = (a.im * (T)b);

	return t;
}

template <typename T>
struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	T s =  (b.re * b.re) + (b.im * b.im);

	t.re = ((a.re * b.re) + (a.im * b.im)) / s;
	t.im = ((a.im * b.re) - (a.re * b.im)) / s;

	return t;
}

template <typename T>
struct simpleComplex<T> operator/(const T &_a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t, a;

	a.re = a;
	a.im = 0.0;

	T s =  (b.re * b.re) + (b.im * b.im);

	t.re = ((a.re * b.re) + (a.im * b.im)) / s;
	t.im = ((a.im * b.re) - (a.re * b.im)) / s;

	return t;
}

template <typename T>
struct simpleComplex<T> operator/(const struct simpleComplex<T> &a, const T &_b)
{
	struct simpleComplex<T> t, b;

    b.re = _b;
    b.im = 0.0;

	T s =  (b.re * b.re) + (b.im * b.im);

	t.re = ((a.re * b.re) + (a.im * b.im)) / s;
	t.im = ((a.im * b.re) - (a.re * b.im)) / s;

	return t;
}

double  simpleComplexAbs (const struct simpleComplex<double> &a)
{
	double f;

	f = sqrt( (a.re * a.re) + (a.im * a.im) );

	return f;
}


double  simpleComplexMod (const struct simpleComplex<double> &a)
{
	double f;

	f = sqrt( (a.re * a.re) + (a.im * a.im) );

	return f;
}

float simpleComplexMod (const struct simpleComplex<float> &a)
{
	float f;

	f = sqrtf( (a.re * a.re) + (a.im * a.im) );

	return f;
}

template <typename T>
simpleComplex<T> simpleComplexConj (const struct simpleComplex<T> &a)
{
	simpleComplex<T> tmp;
	
	tmp = a;
	
	tmp.im = -tmp.im;
	
	return tmp;
}


template <typename T>
struct simpleComplex<T> make_simpleComplex (T r, T i)
{
	simpleComplex<T> t;

	t.re = r ;
	t.im = i ;

	return t;
}


/* vector operators */

template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator+(const T &a, const uVector< struct simpleComplex<T>, SIZE > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T>,SIZE > vtmp;

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a + v[i];
    }

    return vtmp;
}


template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator+(const uVector< struct simpleComplex<T>, SIZE > &v, const T &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T>,SIZE > vtmp;

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] + a;
    }

    return vtmp;
}

template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator-(const T &a, const uVector< struct simpleComplex<T>, SIZE > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T>,SIZE > vtmp;

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a - v[i];
    }

    return vtmp;
}


template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator-(const uVector< struct simpleComplex<T>, SIZE > &v, const T &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T>,SIZE > vtmp;

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] - a;
    }

    return vtmp;
}


template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator+(const uVector< struct simpleComplex<T>, SIZE > &v1, const uVector< struct simpleComplex<T>, SIZE > &v2)
{
    unsigned int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v1.size;i++)
    {
        vtmp[i] = v1[i] + v2[i];
    }

    return vtmp;
}

template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator-(const uVector< struct simpleComplex<T>, SIZE > &v1, const uVector< struct simpleComplex<T>, SIZE > &v2)
{
    unsigned int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v1.size;i++)
    {
        vtmp[i] = v1[i] - v2[i];
    }

    return vtmp;
}


template <typename T, size_t SIZE>
uVector< struct simpleComplex<T>, SIZE > operator*(const T &a, const uVector< struct simpleComplex<T>, SIZE > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}

template <typename T, size_t SIZE>
uVector< struct simpleComplex<T>, SIZE > operator*(const simpleComplex<T> &a, const uVector< struct simpleComplex<T>, SIZE > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}


template <typename T, size_t SIZE>
uVector< struct simpleComplex<T>, SIZE > operator*(const uVector< struct simpleComplex<T>, SIZE > &v, const T &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}

template <typename T, size_t SIZE>
uVector< struct simpleComplex<T>, SIZE > operator*(const uVector< struct simpleComplex<T>, SIZE > &v, const struct simpleComplex<T> &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}


template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator/(const uVector< struct simpleComplex<T>, SIZE > &v, const T &a)
{
   unsigned  int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] / a;
    }

    return vtmp;
}

template <typename T, size_t SIZE>
const uVector< struct simpleComplex<T>, SIZE > operator/(const uVector< struct simpleComplex<T>, SIZE > &v, const struct simpleComplex<T> &a)
{
   unsigned  int i;

    uVector< struct simpleComplex<T>, SIZE > vtmp;
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] / a;
    }

    return vtmp;
}


template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator+(const struct uMatrix<T, SIZE> &m1, const struct uMatrix<T, SIZE> &m2)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m1.rows ; i++)
    {
        for( j=0; j < m1.cols ; j++)
        {
            vtmp(i,j) = m1(i,j) + m2(i,j);
        }
    }

    return vtmp;
}

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator-(const struct uMatrix<T, SIZE> &m1, const struct uMatrix<T, SIZE> &m2)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m1.rows ; i++)
    {
        for( j=0; j < m1.cols ; j++)
        {
            vtmp(i,j) = m1(i,j) - m2(i,j);
        }
    }
    return vtmp;
}


// scalar and matrix mul

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const T &s, const struct uMatrix<T, SIZE> &m)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const struct uMatrix<T, SIZE> &m, const T &s)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}


template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const simpleComplex<T> &s, const struct uMatrix<T, SIZE> &m)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const struct uMatrix<T, SIZE> &m, const simpleComplex<T> &s)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const double &s, const struct uMatrix<T, SIZE> &m)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const struct uMatrix<T, SIZE> &m, const double &s)
{
    unsigned int i, j;
    struct uMatrix<T, SIZE> vtmp;

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}



// matrix and matrix mul

template <typename T, size_t SIZE>
struct uMatrix<T, SIZE> operator*(const struct uMatrix<T, SIZE> &m1, const struct uMatrix<T, SIZE> &m2)
{
    unsigned int i, j, k;
    struct uMatrix<T, SIZE> vtmp;
    T sum;

    sum = make_simpleComplex(0.0, 0.0);

    for (i = 0; i < SIZE; i++)
    {
        for (j = 0; j < SIZE; j++)
        {
            for (k = 0; k < SIZE; k++)
            {
                sum = sum + m1(i,k)*m2(k,j);
            }
            vtmp(i,j) = sum;
            sum = make_simpleComplex(0.0, 0.0);
        }
    }

    return vtmp;
}

#if 0
// matrix and vector mult
template <typename T>
const std::vector< struct simpleComplex<T> > operator*(const struct simpleCmplxMatrix<T> &m, const std::vector< struct simpleComplex<T> > &v)
{
    unsigned int i, k;
    std::vector< struct simpleComplex<T> > vtmp( m.rows );

    for( i=0 ; i<m.cols ; i++)
    {
        vtmp[i] = make_simpleComplex(0.0, 0.0);
        for( k=0; k < m.rows ; k++)
        {
            vtmp[i] = vtmp[i] + m.m[ i * m.rows + k] * v[k];
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}
#endif

template <typename T, size_t SIZE>
const uVector< T, SIZE > operator*(const struct uMatrix<T, SIZE> &m, const uVector< T, SIZE > &v)
{
    unsigned int i, k;
    uVector< T, SIZE > vtmp ;

    for( i=0 ; i<m.cols ; i++)
    {
        vtmp[i] = make_simpleComplex(0.0, 0.0);
        for( k=0; k < m.rows ; k++)
        {
            vtmp[i] = vtmp[i] + m.m[ i * m.rows + k] * v[k];
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

#if 0
template <typename T>
struct simpleComplex<T> expect( int sx, int sy,  const std::vector< struct simpleComplex<T> > &m, const std::vector< struct simpleComplex<T> > &state)
{
        int i, j, k;

        struct simpleComplex<T> n1, n2;
        std::vector< struct simpleComplex<T> > tmp ( sx*sy );


        for(k=0; k < sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }


        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m[i*sy+k] * state[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;
}
#endif

template <typename T, size_t SIZE>
struct simpleComplex<T> expect( size_t sx, size_t sy,  const uVector< struct simpleComplex<T>, SIZE > &m, const uVector< struct simpleComplex<T>, SIZE > &state)
{
        size_t i, j, k;

        struct simpleComplex<T> n1, n2;

        uVector < struct simpleComplex<T>, SIZE*SIZE > tmp;


        for(k=0; k < sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }


        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m[i*sy+k] * state[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;
}

template <typename T, size_t SIZE>
struct simpleComplex<T> expect( size_t sx, size_t sy, uMatrix< simpleComplex<T>, SIZE > &m, const uVector< struct simpleComplex<T>, SIZE > &state)
{
        size_t i, j, k;

        struct simpleComplex<T> n1, n2;

        uVector < struct simpleComplex<T>, SIZE*SIZE > tmp;


        for(k=0; k < sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }


        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m.m[i*sy+k] * state[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;
}

#if 0
template <typename T>
struct simpleComplex<T> expect( int sx, int sy,  struct simpleComplex<T> m[], struct simpleComplex<T> state[] )
{
        int i, j, k;

        struct simpleComplex<T> n1, n2;
        struct simpleComplex<T> tmp [ sx*sy ];


        for(k=0; k < sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }


        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m[i*sy+k] * state[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;
}
#endif

#if 0
template <typename T>
struct simpleComplex<T> expect_cnv_denmat( int sx, int sy,  const std::vector< struct simpleComplex<T> > &m, const std::vector< struct simpleComplex<T> > &state)
{
        int i, j, k;

        struct simpleComplex<T> n1, n2;
        std::vector< struct simpleComplex<T> > tmp ( sx*sy );
        std::vector< struct simpleComplex<T> > matden ( sx*sy );


        for(k=0; k < sx*sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }

        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                n1 = state[i];
                n2 = state[j];
                n2.im = -n2.im;
                matden[ i*sy+j ] = n1 * n2;
            }
        }

        //print_Y( matden );

        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m[i*sy+k] * matden[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;

}
#endif

template <typename T, size_t SIZE1, size_t SIZE2 >
T expect_cnv_denmat( size_t sx, size_t sy,  const uVector< T, SIZE1 > &m, const uVector< T, SIZE2 > &state)
{
        size_t i, j, k;

        T n1, n2;

        //const unsigned int csize = sx*sy;

        uVector< T, SIZE1 > tmp ;
        uVector< T, SIZE1 > matden ;


        for(k=0; k < sx*sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }

        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                n1 = state[i];
                n2 = state[j];
                n2.im = -n2.im;
                matden[ i*sy+j ] = n1 * n2;
            }
        }

        //print_Y( matden );

        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m[i*sy+k] * matden[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;

}

template <typename T, size_t SIZE1, size_t SIZE2 >
T expect_cnv_denmat( size_t sx, size_t sy,  const uMatrix< T, SIZE2 > &m, const uVector< T, SIZE2 > &state)
{
        size_t i, j, k;

        T n1, n2;

        //const unsigned int csize = sx*sy;

        uVector< T, SIZE1 > tmp ;
        uVector< T, SIZE1 > matden ;


        for(k=0; k < sx*sy ; k++)
        {
            tmp[k] = make_simpleComplex(0.0, 0.0);
        }

        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                n1 = state[i];
                n2 = state[j];
                n2.im = -n2.im;
                matden[ i*sy+j ] = n1 * n2;
            }
        }

        //print_Y( matden );

        for(i=0 ; i<sx ; i++)
        {
            for(j=0; j < sy; j++)
            {
                for( k=0; k < sy ; k++)
                {
                    n1 = tmp[i*sy+j];

                    n2 = ( m.m[i*sy+k] * matden[k*sy+j] );

                    tmp[i*sy+j] = n1 + n2;
                }
            }
        }

        n1.re=0;
        n1.im=0;

        for(k=0; k < sy ; k++)
        {
           //cout << " i: " << k*sy + k;

           n1 = n1 + tmp[ k*sy + k ];
        }

        return n1;
}

template <typename T, size_t SIZE1, size_t SIZE2 >
T expect_cnv_denmat_ver2( size_t sx, size_t sy,  const uVector<  T, SIZE1 > &m, const uVector< T, SIZE2 > &state)
{
        size_t i, j, k;

        T n1, n2, md, md1, md2;

        n1.re=0;
        n1.im=0;

        for( k=0; k<SIZE2; k++)
        {
            for(i=0; i<SIZE2; i++)
            {
                md1 = state[k];
                md2 = state[i];
                md2.im = -md2.im;
                md = md1 * md2;

                n2 = m[i*sy+k];

                n1 = n1 + (md * n2);
            }
        }

        return n1;
}

template <typename T, size_t _S_ValueSize, size_t _S_RowPtr, size_t _S_ColInd, size_t v_size >
T expect_cnv_csrdenmat(const struct uCSRMatrix< T, _S_ValueSize, _S_RowPtr, _S_ColInd> &m, const uVector< T, v_size > &state)
{
    struct simpleComplex<double> r, tmp;
    uVector< simpleComplex<double>, v_size > vtmp;
    size_t i;

    for(i=0;i<v_size;i++)
        vtmp[i] = make_simpleComplex( 0.0,  0.0 );

    vtmp = mulCSRMatByuVec(m, state);

    for(i=0 ; i<v_size ; i++)
    {
        tmp.re =  state[i].re;
        tmp.im = -state[i].im;
        vtmp[i] = tmp * vtmp[i];
    }

    r.re=0;
    r.im=0;
    r = sum( vtmp );

    return r;
}

template <typename T, size_t SIZE1, size_t SIZE2 >
const uVector< struct simpleComplex<T>, SIZE2 >  mul_mat_vec ( size_t sx, size_t sy,  const uVector< struct simpleComplex<T>, SIZE1 > &m, const uVector< struct simpleComplex<T>, SIZE2 > &v)
{
    size_t i, k;
    uVector< struct simpleComplex<T>, SIZE2> vtmp ;

    for( i=0 ; i<sy ; i++)
    {
        vtmp[i] = make_simpleComplex(0.0, 0.0);
        for( k=0; k < sx ; k++)
        {
            vtmp[i] = vtmp[i] + m[ i * sy + k] * v[k];
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)


    return vtmp;

}

#if 0
template <typename T>
T norm( const std::vector< struct simpleComplex<T> > &m )
{
    int i;
    simpleComplex<T>  tmp, cnj;
    T v = 0.0;

    tmp.re = 0.0;
    tmp.im = 0.0;

    for (i = 0 ; i < m.size() ; i++)
    {
        cnj = m[i];
        cnj.im =-cnj.im;
        tmp = tmp + cnj * m[i];
    }

    v = sqrt(tmp.re * tmp.re + tmp.im * tmp.im);

    return v;
}
#endif

template <typename T, size_t SIZE>
T norm( const uVector< struct simpleComplex<T>, SIZE > &m )
{
    size_t i;
    simpleComplex<T>  tmp, cnj;
    T v = 0.0;

    tmp.re = 0.0;
    tmp.im = 0.0;

    for (i = 0 ; i < m.size ; i++)
    {
        cnj = m[i];
        cnj.im =-cnj.im;
        tmp = tmp + (cnj * m[i]);
    }

    v = sqrt(tmp.re * tmp.re + tmp.im * tmp.im);

    return v;
}


#if 0
template <typename T>
T normsqrt( const std::vector< struct simpleComplex<T> > &m )
{
    T v = 0.0;

    v = norm( m );

    return sqrt(v);
}
#endif

template <typename T, size_t SIZE>
T normsqrt( const uVector< struct simpleComplex<T>, SIZE > &m )
{
    T v = 0.0;

    v = norm( m );

    return sqrt(v);
}

#if 0
template <typename T>
struct simpleComplex<T> sum( const std::vector< struct simpleComplex<T> > &m )
{
    int i;
    simpleComplex<T>  tmp;

    tmp.re = 0.0;
    tmp.im = 0.0;

    for (i = 0 ; i < m.size() ; i++)
    {
        tmp = tmp + m[i];
    }

    return tmp;
}
#endif

template <typename T, size_t SIZE>
struct simpleComplex<T> sum( const uVector< struct simpleComplex<T>, SIZE > &m )
{
    size_t i;
    simpleComplex<T>  tmp;

    tmp.re = 0.0;
    tmp.im = 0.0;

    for (i = 0 ; i < m.size ; i++)
    {
        tmp = tmp + m[i];
    }

    return tmp;
}

#if 0
template <typename T>
void normalize( std::vector< struct simpleComplex<T> > &m )
{

    T v;
    v = normsqrt( m );

    m = m / v;
    return ;
}
#endif

template <typename T, size_t SIZE>
void normalize( uVector< struct simpleComplex<T>, SIZE > &m )
{

    T v;
    v = normsqrt( m );

    m = m / v;
    return ;
}

template <typename T, size_t SIZE>
void std_base_state( uVector< struct simpleComplex<T>, SIZE > &m, int state )
{
    size_t i;

    for (i = 0 ; i < m.size ; i++)
    {
		if( i != state ) {
			m[i].re = 0.0;
			m[i].im = 0.0;
		} else {
			m[i].re = 1.0;
			m[i].im = 0.0;
		}
    }
}

template <typename T, size_t SIZE1, size_t SIZE2>
uVector< struct simpleComplex<T>, SIZE1*SIZE2 > dagnotdag( const uVector< struct simpleComplex<T>, SIZE1*SIZE2 > &m )
{
    size_t i;

	uVector< struct simpleComplex<T>, SIZE1*SIZE2 > r;
	uMatrix< simpleComplex<T>, SIZE1> d,nd;
	
	d.cols = SIZE1;
	d.rows = SIZE2;
	d.m = m;

	nd.cols = SIZE1;
	nd.rows = SIZE2;
	nd.m = m;
	
	d = dagger(nd);
	
	d  = d * nd;
	
	r.size=SIZE1 * SIZE2;
	r = d.m;
	
	return r;
}


template <typename T, size_t SIZE>
void std_base_state( simpleComplex<T> *tbl, int state )
{
    size_t i;

    for (i = 0 ; i < SIZE ; i++)
    {
		if( i == state ) {
			tbl[i].re = 1.0;
			tbl[i].im = 0.0;
		} else {
			tbl[i].re = 0.0;
			tbl[i].im = 0.0;
		}
    }
}

template <typename T, size_t SIZE>
uVector< struct simpleComplex<T>, SIZE > coherent( struct simpleComplex<T> alpha )
{
	uVector< struct simpleComplex<T>, SIZE > x;
	uMatrix< struct simpleComplex<T>, SIZE > a, adag, tmpmat, D;

    size_t i;

    for(i=0; i<x.size; i++)
    {
        x[i] = make_simpleComplex(0.0, 0.0);
    }
	
	std_base_state(x, 0);
	
	destroy_operator( a );
	

	adag = dagger( a );
	
	tmpmat = (alpha * adag) - (simpleComplexConj(alpha) * a);
	
	exp_of_matrix(tmpmat, 10, D);
	
	return D * x;
}


//------------------------------------------------------------------------------------------------------------
// global output operator<< for simpleComplex<T>
//------------------------------------------------------------------------------------------------------------

template <typename T>
inline ostream& operator<< (ostream & output, const simpleComplex<T> & c)
{

#define width 6

  //fixed
  output << right << setprecision(width) << fixed << "(" << setw(width+3) << c.re << ", " << setw(width+3) << c.im << ")";

  //scientific
  //output << right << setprecision(width) << scientific << "(" << setw(width+7) <<  c.re << ", " << setw(width+7) << c.im << ")";

  //normal
  //output << right << setprecision(width) << "(" << setw(width+6) << c.re << ", " << setw(width+6) << c.im << ")";

#undef width

  return output;
}


//------------------------------------------------------------------------------------------------------------
//  global output operator<< for uMatrix<T>
//------------------------------------------------------------------------------------------------------------


template <typename T, size_t SIZE>
ostream& operator<< (ostream & output, const uMatrix<T, SIZE> &c)
{
  unsigned int i, j;
#define width 4
  simpleComplex<T> nc;

  for (i = 0; i < c.rows; ++i) {
    for (j = 0; j < c.cols-1; ++j)
      output << setw(width) << c(i,j) << ", ";

    if (i<c.rows-1)
      output << setw(width) << c(i,j) << endl;
    else
      output << setw(width) << c(i,j);
  }
#undef width  
  return output;
}

//------------------------------------------------------------------------------------------------------------
//  global output operator<< for uVector<T>
//------------------------------------------------------------------------------------------------------------


template <typename T, size_t SIZE>
ostream& operator<< (ostream & output, const uVector<T, SIZE> &c)
{
  unsigned int i, j;
#define width 4
  simpleComplex<T> nc;

  output << setw(width) << c[i] << " ";
  for (i = 1; i < c.size-1; ++i) {
      output << setw(width) << c[i] << ", ";
  }
      output << setw(width) << c[i] << endl;
#undef width  
  return output;
}


//------------------------------------------------------------------------------------------------------------
template <typename T, size_t SIZE>
inline void print_matrix(const uMatrix<T, SIZE> &A, string name)
{
  cout << name << "[" << A.rows << "x" << A.cols << "=" << A.size() << "]:\n" << A << endl;
}

template <typename T, size_t SIZE>
bool inverse_of_matrix(const uMatrix<T, SIZE> &M, uMatrix<T, SIZE> &A_inv)
{
  // check sizes
  unsigned int size_m = M.size();
  unsigned int size_i = A_inv.size();
  if (size_m != size_i || M.rows != A_inv.rows || M.cols != A_inv.cols) {
    const int width = log10(max(size_m, size_i)) + 1;
    const int width_rows = log10(max(M.rows, A_inv.rows)) + 1;
    const int width_cols = log10(max(M.cols, A_inv.cols)) + 1;
/*
    cerr << "Error: " << __FILE__ <<" (" << __LINE__ << "), "<< __FUNCTION__ << "(m1, m2):\n  "
        << "Different sizes of matrices m1 and m2.\n" << right
        << "    size(m1) = " << setw(width) << size_m << ": rows = " << setw(width_rows) <<     M.rows << ", cols = " << setw(width_cols) << M.cols << "\n"
        << "    size(m2) = " << setw(width) << size_i << ": rows = " << setw(width_rows) << A_inv.rows << ", cols = " << setw(width_cols) << A_inv.cols << endl;
*/		
    return false;
  }

  long int n = M.rows;       //size of the matrix M
  long int i, j, k;          //indices

  uMatrix<T, SIZE> A; //create copy of matrix M as A, because matrix A change its values

  for(i=0;i<M.size();i++)
    A.m[i]=M.m[i];

  uMatrix<T, SIZE> L, U;          // matrices
  uMatrix<T, SIZE> b, d, x;  // vectors as matrices with one column size (n,1)

  for(i=0;i<b.size();i++)
  {
    b.m[i] = make_simpleComplex(0.0, 0.0);
    d.m[i] = make_simpleComplex(0.0, 0.0);
    x.m[i] = make_simpleComplex(0.0, 0.0);
  }

  b.cols = 1;
  d.cols = 1;
  x.cols = 1;

  // simpleComplex numbers 0 and 1
  simpleComplex<T> zero = make_simpleComplex(0.0, 0.0);  // 0 + j*0
  simpleComplex<T> one  = make_simpleComplex(1.0, 0.0);  // 1 + j*0

  // forward elimination
  for (k=0; k < n-1; ++k) {
    for (i=k+1; i < n; ++i) {
      L(i,k) = A(i,k) / A(k,k);
      for (j=k+1; j < n; ++j)
        A(i,j) = A(i,j) - L(i,k)*A(k,j);
    }
  }

  // L matrix is a matrix of the elimination coefficient and the diagonal elements are 1.0
  for (i=0; i < n; ++i)
    L(i,i) = one;

  // U matrix is the upper triangular part of A
  for (j=0; j<n; ++j)
    for (i=0; i<=j; ++i)
      U(i,j) = A(i,j);

  // MAIN LOOP
  // compute columns of the inverse matrix A_inv
  for (k=0; k<n; ++k) {
    b(k,0) = one;
    d(0,0) = b(0,0);

    // Solve 'L*d = b' using the forward substitution
    for (i=1; i<n; ++i) {
      d(i,0) = b(i,0);

      for (j=0; j <= i-1; ++j)
        d(i,0) = d(i,0) - L(i,j)*d(j,0);
    }

    //Solve 'U*x = d' using the back substitution
    x(n-1,0) = d(n-1,0)/U(n-1,n-1);
    for (i = n-2; i>=0; --i) {
      x(i,0) = d(i,0);

      for (j=n-1; j>=i+1; --j)
        x(i,0) = x(i,0) - U(i,j)*x(j,0);

      x(i,0) = x(i,0)/U(i,i);
    }

    //fill the solutions 'x' into column 'k' of 'A_inv'
    for (i=0; i<n; ++i)
      A_inv(i,k) = x(i,0);

    b(k,0) = zero;
  }

  // solved successfully
  return true;
}

template <typename T, int SIZE>
void matprod(uMatrix<T, SIZE> &A, uMatrix<T, SIZE> &B, uMatrix<T, SIZE> &C);

template <typename T, int SIZE>	
void matcopy(uMatrix<T, SIZE> &A, uMatrix<T, SIZE> &B);

template <typename T, int SIZE>
static void matexp_pade(const long int p, uMatrix< simpleComplex<T>, SIZE> &A, uMatrix< simpleComplex<T>, SIZE> &N);


template <typename T, size_t SIZE>
void eye_of_matrix(uMatrix<T, SIZE> &m)
{
    size_t i, j;
    for( i=0 ; i < m.rows ; i++)
    {
        for( j=0 ; j < m.rows ; j++)
        {
            if(i == j)
                m(i,j) = make_simpleComplex(1.0, 0.0);
            else
                m(i,j) = make_simpleComplex(0.0, 0.0);
        }

    }
}

template <typename T, size_t SIZE>
void zero_matrix( uMatrix< T, SIZE > &m )
{
    size_t i;

    for(i=0; i<m.size(); i++)
    {
       m.m[i] = make_simpleComplex(0.0, 0.0);
    }
}

template <typename T, size_t SIZE>
void pauli_x_matrix( uMatrix< T, SIZE > &m )
{
    m.m[0] = make_simpleComplex(0.0, 0.0); m.m[1] = make_simpleComplex(1.0, 0.0);
    m.m[2] = make_simpleComplex(1.0, 0.0); m.m[3] = make_simpleComplex(0.0, 0.0);	
}

template <typename T, size_t SIZE>
void pauli_x_matrix( uVector< T, SIZE > &v)
{
    v[0] = make_simpleComplex(0.0, 0.0); v[1] = make_simpleComplex(1.0, 0.0);
    v[2] = make_simpleComplex(1.0, 0.0); v[3] = make_simpleComplex(0.0, 0.0);	
}

template <typename T, size_t SIZE>
void pauli_y_matrix( uVector< T, SIZE > &v)
{
    v[0] = make_simpleComplex(0.0, 0.0); v[1] = make_simpleComplex(0.0,-1.0);
    v[2] = make_simpleComplex(0.0, 1.0); v[3] = make_simpleComplex(0.0, 0.0);	
}

template <typename T, size_t SIZE>
void pauli_z_matrix( uVector< T, SIZE > &v)
{
    v[0] = make_simpleComplex(1.0, 0.0); v[1] = make_simpleComplex( 0.0, 0.0);
    v[2] = make_simpleComplex(0.0, 0.0); v[3] = make_simpleComplex(-1.0, 0.0);	
}

template <typename T, size_t SIZE>
void transpose_of_matrix( uMatrix< T, SIZE > &m )
{
    T a, b;
    size_t i,j;

    if(m.rows == m.cols)
    {
        for(i=0 ; i<m.rows ; i++)
        {
            for(j=i ; j<m.cols ; j++)
            {
                a = m(i,j);
                a.im = -a.im;
				
                b = m(j,i);
                b.im = -b.im;
				
                m(i, j) = b;
                m(j, i) = a;
            }
        }
    }
}

template <typename T, size_t SIZE>
uMatrix< T, SIZE > dagger( const uMatrix< T, SIZE > &_m )
{	
	uMatrix< T, SIZE > m;
	
	m = _m;
	
    T a, b;
    size_t i,j;

    for(i=0 ; i<m.rows ; i++)
    {
		for(j=i ; j<m.cols ; j++)
		{
			a = m(i,j);
			a.im = -a.im;
				
			b = m(j,i);
			b.im = -b.im;
				
            m(i, j) = b;
            m(j, i) = a;
            }
	}
	
	return m;
}


template <typename T, size_t SIZE>
void destroy_operator( uMatrix< T, SIZE > &m )
{
	size_t i, y;
	
	y=1;
	for(i=0; i < SIZE - 1; i++)
	{
		m(i, y) = make_simpleComplex( sqrt((double)y), 0.0 );
		y++;
	}
}

template <typename T, size_t SIZE>
void create_operator( uMatrix< T, SIZE > &m )
{
	size_t i, y;
	
	y=1;
	for(i=0; i < SIZE - 1; i++)
	{
		m(y, i) = make_simpleComplex( sqrt((double)y), 0.0 );
		y++;
	}
}


template <typename T, size_t SIZE>
void make_x_operator( uMatrix< T, SIZE > &m )
{
	size_t i, y;
	
	y=SIZE - 1;
	for(i=0; i < SIZE ; i++)
	{
		m(i, y) = make_simpleComplex( 1.0, 0.0 );
		y--;
	}
}

template <typename T, size_t SIZE>
void sigma_m_matrix( uMatrix< T, SIZE > &v)
{
    v(0,0) = make_simpleComplex(0.0, 0.0); v(0,1) = make_simpleComplex(0.0, 0.0);
    v(1,0) = make_simpleComplex(1.0, 0.0); v(1,1) = make_simpleComplex(0.0, 0.0);	
}


template <typename T, size_t SIZE>
T norminf(const uMatrix<T, SIZE> &m)
{
    size_t i,j;

    T t = 0.0;
    T r = 0.0;

    for(i=0; i<m.rows; i++)
    {
        t = 0.0;
        for(j=0; j<m.cols; j++)
        {
            t = t + sqrt( m(i,j).re * m(i,j).re + m(i,j).im * m(i,j).im );
        }
        if(t > r) r=t;
    }

    return r;

}

template <typename T, size_t SIZE>
T norminf( const uVector< T, SIZE > &v)
{
    size_t i;

    T t;
    T r = sqrt( v[0].re * v[0].re + v[0].im * v[0].im );;


    for(i=1; i<v.size; i++)
    {
        t = sqrt( v[i].re * v[i].re + v[i].im * v[i].im );

        if( t > r) r = t;
    }


    return r;

}

template <typename T, size_t SIZE>
void zero_vector( uVector< struct simpleComplex<T>, SIZE > &v)
{
    size_t i;

    for(i=0; i<v.size; i++)
    {
        v[i] = make_simpleComplex(0.0, 0.0);
    }
}

template <typename T, size_t SIZE1, size_t SIZE2>
uMatrix<simpleComplex<T>, SIZE1*SIZE2>  tensor(uMatrix<simpleComplex<T>, SIZE1> &m1, uMatrix<simpleComplex<T>, SIZE2> &m2)
{
	int x,y,i,j,ii,jj;
	uMatrix< simpleComplex<T>, SIZE1*SIZE2>  tmp;
    simpleComplex<T> num;

	for(i=0;i<tmp.rows;i++)
    {
		for(j=0;j<tmp.cols;j++)
		{
			tmp(i,j).re=0.0;
			tmp(i,j).im=0.0;
		}
	}

    ii=0;
    jj=0;
     for(x=0;x<m1.rows;x++)
     {
         for(y=0;y<m1.cols;y++)
         {
             ii=(x*m2.rows);
             jj=(y*m2.cols);
             for(i=0;i<m2.rows;i++)
             {
                 for(j=0;j<m2.cols;j++)
                 {
                     num.re=0.0;
                     num.im=0.0;

					 num = m1(x,y) * m2(i,j);
					 
					 tmp(ii,jj) = num;
					 
                     jj++;
                 }
                ii++;
                jj-=m2.cols;
             }
         }
     }
	return tmp;
}

template <typename T, size_t SIZE>
uMatrix<T, SIZE*SIZE*SIZE>  tensor(const uMatrix<simpleComplex<T>, SIZE> &m1, const uMatrix<simpleComplex<T>, SIZE> &m2, const uMatrix<simpleComplex<T>, SIZE> &m3)
{
	uMatrix<T, SIZE*SIZE*SIZE>  tmp;
	
	return tmp;
}


#endif // __complexnum_h__

