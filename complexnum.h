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

template <typename T>
struct uVector {

    unsigned int size;
    T *m;

    uVector()
	{
		size = 0;
		m = nullptr;
	};
	
	uVector(const unsigned int v_size)
	{
		size = v_size;
		m = new T[v_size];
	}

	uVector( const uVector &v )
	{
		unsigned int i;
		size = v.size;
		
		m = new T[v.size];
		
		for(i=0;i<v.size;i++)
		{
			m[i] = v.m[i];
		}
	}
	
	uVector& operator=( const uVector &v )
	{
		unsigned int i;
		
		if(m!=nullptr) delete [] m;
		
		size = v.size;
		
		m = new T[v.size];
		
		for(i=0;i<v.size;i++)
		{
			m[i] = v.m[i];
		}
		
		return *this;
	}
	
	void recreate( unsigned int l1, T *_m)
	{
		unsigned int i;
		
		if( m != nullptr ) delete [] m;

        size = l1;
		
		m = new T[ l1 ];
		
		for(i=0;i<size;i++)
		{
			m[i] = _m[i];
		}
	}	
	
	~uVector()
	{
		if( m!=nullptr) delete [] m;
	}
	

    T& operator[](const unsigned int idx) { return m[idx]; };
    const T& operator[](const unsigned int idx) const { return m[idx]; };

};

template <typename T>
struct uMatrix {

    unsigned int size;
	unsigned int rows;
	unsigned int cols;

    T *m;

    uMatrix() {
        size = 0;
        rows = 0;
        cols = 0;
		
		m = nullptr;
    }

    uMatrix(const unsigned int l) {
        size = l*l;
        rows  = l;
        cols  = l;
		
		m = new T[l*l];
    }

    uMatrix(const unsigned int l1, const unsigned int l2) {
        size = l1*l2;
        rows  = l1;
        cols  = l2;
		
		m = new T[l1*l2];
    }
	
	uMatrix( const uMatrix &_m )
	{
		unsigned int i;
		
        size = _m.rows * _m.cols;
        rows  = _m.rows;
        cols  = _m.cols;
		
		m = new T[ _m.rows * _m.cols ];
		
		for(i=0;i<_m.size;i++)
		{
			m[i] = _m.m[i];
		}		
	}

	uMatrix& operator=( const uMatrix &_m )
	{
		unsigned int i;
		
		if( m != nullptr ) delete [] m;

        size = _m.rows * _m.cols;
        rows  = _m.rows;
        cols  = _m.cols;
		
		m = new T[ _m.rows * _m.cols ];
		
		for(i=0;i<_m.size;i++)
		{
			m[i] = _m.m[i];
		}
		
		return *this;
	}

	void recreate( unsigned int l1, unsigned int l2, const uVector<T> &_m)
	{
		unsigned int i;
		
		if( m != nullptr ) delete [] m;

        size = l1 * l2;
        rows  = l1;
        cols  = l2;
		
		m = new T[ rows * cols ];
		
		for(i=0;i<size;i++)
		{
			m[i] = _m[i];
		}
	}
	
	~uMatrix()
	{
		if( m != nullptr ) delete [] m;
	}
	
    // function operator() for indexing:  for reading (r)-rows, (c)-cols
    inline const T & operator()(const unsigned int r, const unsigned int c) const
    { return m[r*cols + c]; }

    // function operator() for indexing:  for writing (r)-rows, (c)-cols
    inline T & operator()(const unsigned int r, const unsigned int c)
    { return m[r*cols + c]; }
};

template <typename T>
struct uCSRMatrix {

    unsigned int _values_size,  _row_ptr, _col_ind;
	
    T *values;
    unsigned int *row_ptr;
    unsigned int *col_ind;	

    uCSRMatrix()
    {
        _values_size = 0;
        _row_ptr = 0;
        _col_ind = 0;
		
		values = nullptr;
		row_ptr = nullptr;
		col_ind = nullptr;
    }

	
    uCSRMatrix(size_t _S_ValueSize, size_t _S_RowPtr, size_t _S_ColInd)
    {
        _values_size = _S_ValueSize;
        _row_ptr = _S_RowPtr;
        _col_ind = _S_ColInd;
		
		values = new T[_values_size];
		row_ptr = new unsigned int[_row_ptr];
		col_ind = new unsigned int[_col_ind];
    }
	
	uCSRMatrix(const uCSRMatrix &m)
	{
		unsigned int i;

        _values_size = m._values_size;
        _row_ptr = m._row_ptr;
        _col_ind = m._col_ind;
		
		values = new T[_values_size];
		row_ptr = new unsigned int[_row_ptr];
		col_ind = new unsigned int[_col_ind];
		
		for(i=0;i<_values_size;i++)
		{
			values[i] = m.values[i];
			col_ind[i] = m.col_ind[i];
		}
		
		for(i=0;i<_row_ptr;i++)
		{
			row_ptr[i] = m.row_ptr[i];
		}
	}
	
	uCSRMatrix& operator=( const uCSRMatrix &m )
	{
		unsigned int i;

		if(col_ind != nullptr) delete [] col_ind;
		if(row_ptr != nullptr) delete [] row_ptr;
		if(values != nullptr) delete [] values;
		
        _values_size = m._values_size;
        _row_ptr = m._row_ptr;
        _col_ind = m._col_ind;
		
		values = new T[_values_size];
		row_ptr = new unsigned int[_row_ptr];
		col_ind = new unsigned int[_col_ind];
		
		for(i=0;i<_values_size;i++)
		{
			values[i] = m.values[i];
			col_ind[i] = m.col_ind[i];
		}
		
		for(i=0;i<_row_ptr;i++)
		{
			row_ptr[i] = m.row_ptr[i];
		}
		
		return *this;
	}
	
	~uCSRMatrix()
	{
		if( col_ind != nullptr) delete [] col_ind;
		if( row_ptr != nullptr) delete [] row_ptr;
		if( values != nullptr) delete [] values;
	}

};

template<typename T>
void csr_matrix_scan(const uMatrix< simpleComplex<T> > &m, size_t &_values_size_ref, size_t &_row_ptr_ref, size_t &_col_ind_ref)
{
	int i;
	size_t _values_size = 0,  _row_ptr = 0, _col_ind = 0;
	
	simpleComplex<T> t;
	
	for(i=0;i<m.size;i++)
	{
		t = m.m[i];
		
		if( !((t.re ==0) && (t.im==0)) ) _values_size++;
	}
	
	_row_ptr = m.cols + 1;
	
	_values_size_ref = _values_size;
	_row_ptr_ref = _row_ptr;
	_col_ind_ref = _values_size;
}

template<typename T>
uCSRMatrix< simpleComplex<T> > uMatrix_to_uCSRMatrix(const uMatrix< simpleComplex<T> > &m)
{
	size_t _values_size = 0,  _row_ptr = 0, _col_ind = 0, v, r, i, j;

	simpleComplex<T> t;

	csr_matrix_scan(m, _values_size,  _row_ptr, _col_ind);

	uCSRMatrix< simpleComplex<T> > retcsrmat(_values_size,  _row_ptr, _col_ind);
	
	v=0; r=0;
	for(i=0;i<m.rows;i++)
	{
		for(j=0;j<m.cols;j++)
		{
			t = m(i,j);
			if( !((t.re ==0) && (t.im==0)) )
			{
				retcsrmat.values[v]=t;
				retcsrmat.col_ind[v]=j;
				v++;
				r++;
			}
		}
		retcsrmat.row_ptr[i+1]=r;
	}
	retcsrmat.row_ptr[0]=0;
	
	return retcsrmat;
}

template<typename T>
uVector< simpleComplex<T> > uMatrix_to_uVector(const uMatrix< simpleComplex<T> > &m)
{
	uVector< simpleComplex<T> > t;
	
	t.recreate( m.size, m.m);
	
	return t;
}

template<typename T>
uVector< simpleComplex<T> > mulCSRMatByuVec(const uCSRMatrix< simpleComplex<T> > &m, const uVector< simpleComplex<T> > &x)
{
    size_t i, j;


    uVector< simpleComplex<T> > y(x.size);

    for ( i=0; i < x.size; i++)
    {
        y[i] = make_simpleComplex(0.0,0.0);
    }

    for ( i=0; i < x.size; i++)
    {
        for ( j=m.row_ptr[i] ; j < m.row_ptr[i+1] ; j++)
        {
            y[i]= y[i] + (m.values[j] * x[m.col_ind[j]]);
        }
    }

    return y;
}

template <typename T>
void exp_of_matrix(uMatrix< simpleComplex<T> > &x, const long int p, uMatrix< simpleComplex<T> > &ret);

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

template <typename T>
const uVector< struct simpleComplex<T> > operator+(const T &a, const uVector< struct simpleComplex<T> > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp(v.size);

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a + v[i];
    }

    return vtmp;
}


template <typename T>
const uVector< struct simpleComplex<T> > operator+(const uVector< struct simpleComplex<T> > &v, const T &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] + a;
    }

    return vtmp;
}

template <typename T>
const uVector< struct simpleComplex<T> > operator-(const T &a, const uVector< struct simpleComplex<T> > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a - v[i];
    }

    return vtmp;
}


template <typename T>
const uVector< struct simpleComplex<T> > operator-(const uVector< struct simpleComplex<T> > &v, const T &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );

    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] - a;
    }

    return vtmp;
}


template <typename T>
const uVector< struct simpleComplex<T> > operator+(const uVector< struct simpleComplex<T> > &v1, const uVector< struct simpleComplex<T> > &v2)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp(v1.size);
    for(i=0;i<v1.size;i++)
    {
        vtmp[i] = v1[i] + v2[i];
    }

    return vtmp;
}

template <typename T>
const uVector< struct simpleComplex<T> > operator-(const uVector< struct simpleComplex<T> > &v1, const uVector< struct simpleComplex<T> > &v2)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v1.size );
    for(i=0;i<v1.size;i++)
    {
        vtmp[i] = v1[i] - v2[i];
    }

    return vtmp;
}


template <typename T>
uVector< struct simpleComplex<T> > operator*(const T &a, const uVector< struct simpleComplex<T> > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}

template <typename T>
uVector< struct simpleComplex<T> > operator*(const simpleComplex<T> &a, const uVector< struct simpleComplex<T> > &v)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );
    
	for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}


template <typename T>
uVector< struct simpleComplex<T> > operator*(const uVector< struct simpleComplex<T> > &v, const T &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );
	
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}

template <typename T>
uVector< struct simpleComplex<T> > operator*(const uVector< struct simpleComplex<T> > &v, const struct simpleComplex<T> &a)
{
    unsigned int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );
	
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = a * v[i];
    }

    return vtmp;
}


template <typename T>
const uVector< struct simpleComplex<T> > operator/(const uVector< struct simpleComplex<T> > &v, const T &a)
{
	unsigned  int i;

	uVector< struct simpleComplex<T> > vtmp( v.size );
	
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] / a;
    }

    return vtmp;
}

template <typename T>
const uVector< struct simpleComplex<T> > operator/(const uVector< struct simpleComplex<T> > &v, const struct simpleComplex<T> &a)
{
   unsigned  int i;

    uVector< struct simpleComplex<T> > vtmp( v.size );
    for(i=0;i<v.size;i++)
    {
        vtmp[i] = v[i] / a;
    }

    return vtmp;
}


template <typename T>
struct uMatrix<T> operator+(const struct uMatrix<T> &m1, const struct uMatrix<T> &m2)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp(m1.rows, m1.cols);

    for( i=0 ; i<m1.rows ; i++)
    {
        for( j=0; j < m1.cols ; j++)
        {
            vtmp(i,j) = m1(i,j) + m2(i,j);
        }
    }

    return vtmp;
}

template <typename T>
struct uMatrix<T> operator-(const struct uMatrix<T> &m1, const struct uMatrix<T> &m2)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp(m1.rows, m1.cols);

    for( i=0 ; i<m1.rows ; i++)
    {
        for( j=0; j < m1.cols ; j++)
        {
            vtmp(i,j) = m1(i,j) - m2(i,j);
        }
    }
    return vtmp;
}

// scalar and matrix difference

template <typename T>
struct uMatrix< simpleComplex<T> > operator-(const simpleComplex<T> &s, const struct uMatrix< simpleComplex<T> > &m)
{
    unsigned int i, j;
    struct uMatrix< simpleComplex<T> > vtmp(m.rows, m.cols);

    for( i=0 ; i < m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s - m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T>
struct uMatrix<simpleComplex<T> > operator-(const struct uMatrix<simpleComplex<T> > &m, const simpleComplex<T> &s)
{
    unsigned int i, j;
    struct uMatrix<simpleComplex<T> > vtmp( m.rows, m.cols);

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = m(i,j) - s;
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

// scalar and matrix mul

template <typename T>
struct uMatrix<T> operator*(const T &s, const struct uMatrix<T> &m)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp( m.rows, m.cols);

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T>
struct uMatrix<T> operator*(const struct uMatrix<T> &m, const T &s)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp( m.rows, m.cols );

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}


template <typename T>
struct uMatrix<T> operator*(const simpleComplex<T> &s, const struct uMatrix<T> &m)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp(m.rows, m.cols);

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T>
struct uMatrix<T> operator*(const struct uMatrix<T> &m, const simpleComplex<T> &s)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp(m.rows, m.cols);

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T>
struct uMatrix<T> operator*(const double &s, const struct uMatrix<T> &m)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp(m.rows, m.cols);

    for( i=0 ; i<m.rows ; i++)
    {
        for( j=0; j < m.cols ; j++)
        {
            vtmp(i,j) = s * m(i,j);
        } // for( k=0; k < sy ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T>
struct uMatrix<T> operator*(const struct uMatrix<T> &m, const double &s)
{
    unsigned int i, j;
    struct uMatrix<T> vtmp(m.rows, m.cols);

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

template <typename T>
struct uMatrix<T> operator*(const struct uMatrix<T> &m1, const struct uMatrix<T> &m2)
{
    unsigned int i, j, k;
    struct uMatrix<T> vtmp(m1.rows, m2.cols);
    T sum;

    sum = make_simpleComplex(0.0, 0.0);

    for (i = 0; i < m1.rows; i++)
    {
        for (j = 0; j < m2.cols; j++)
        {
            for (k = 0; k < m1.cols; k++)
            {
                sum = sum + m1(i,k)*m2(k,j);
            }
            vtmp(i,j) = sum;
            sum = make_simpleComplex(0.0, 0.0);
        }
    }

    return vtmp;
}

template <typename T>
const uVector< T > operator*(const struct uMatrix<T> &m, const uVector< T > &v)
{
    unsigned int i, k;
    uVector< T > vtmp( v.size ) ;

    for( i=0 ; i< m.cols ; i++)
    {
        vtmp[i] = make_simpleComplex(0.0, 0.0);
        for( k=0; k < m.rows ; k++)
        {
            vtmp[i] = vtmp[i] + m.m[ i * m.rows + k] * v[k];
        } // for( k=0; k < m.rows ; k++)
    } // for( i=0 ; i<m.cols ; i++)

    return vtmp;
}

template <typename T>
struct simpleComplex<T> expect( const unsigned int sx, const unsigned int sy, const uVector< struct simpleComplex<T> > &m, const uVector< struct simpleComplex<T> > &state)
{
        unsigned int i, j, k;

        struct simpleComplex<T> n1, n2;

        uVector < struct simpleComplex<T> > tmp( sx*sy );


        for(k=0; k < sx * sy ; k++)
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

template <typename T>
struct simpleComplex<T> expect( size_t sx, size_t sy, uMatrix< simpleComplex<T> > &m, const uVector< struct simpleComplex<T> > &state)
{
        size_t i, j, k;

        struct simpleComplex<T> n1, n2;

        uVector < struct simpleComplex<T> > tmp( sx*sy );


        for(k=0; k < sx*sy ; k++)
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

template <typename T, size_t SIZE1, size_t SIZE2 >
T expect_cnv_denmat( size_t sx, size_t sy,  const uVector< T > &m, const uVector< T > &state)
{
        size_t i, j, k;

        T n1, n2;

        //const unsigned int csize = sx*sy;

        uVector< T > tmp( sx * sy ) ;
        uVector< T > matden( sx * sy ) ;


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

template <typename T>
T expect_cnv_denmat( size_t sx, size_t sy,  const uMatrix< T > &m, const uVector< T > &state)
{
        size_t i, j, k;

        T n1, n2;

        //const unsigned int csize = sx*sy;

        uVector< T > tmp( sx * sy) ;
        uVector< T > matden( sx * sy ) ;


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

template <typename T>
T expect_cnv_denmat_ver2( size_t sx, size_t sy,  const uVector< T > &m, const uVector< T > &state)
{
        size_t i, j, k;

        T n1, n2, md, md1, md2;

        n1.re=0;
        n1.im=0;

        for( k=0; k<sx; k++)
        {
            for(i=0; i<sy; i++)
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

template <typename T >
simpleComplex<T> expect_cnv_csrdenmat(const uCSRMatrix< simpleComplex<T> > &m, const uVector< simpleComplex<T> > &state)
{
    simpleComplex<T> r, tmp;
    uVector< simpleComplex<T> > vtmp( state.size );
	
    size_t i;

    for(i=0;i<state.size;i++)
        vtmp[i] = make_simpleComplex( 0.0,  0.0 );

    vtmp = mulCSRMatByuVec(m, state);

    for(i=0 ; i<state.size ; i++)
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

template <typename T>
const uVector< struct simpleComplex<T> >  mul_mat_vec ( size_t sx, size_t sy,  const uVector< struct simpleComplex<T> > &m, const uVector< struct simpleComplex<T> > &v)
{
    size_t i, k;
	
    uVector< struct simpleComplex<T> > vtmp( sy ) ;

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

template <typename T>
void  mul_mat_vec ( const uMatrix< struct simpleComplex<T> > &m, simpleComplex<T> *vin, simpleComplex<T> *vout)
{
    size_t i, k;
	
    for( i=0 ; i< m.cols ; i++)
    {
        vout[i] = make_simpleComplex(0.0, 0.0);
        for( k=0; k < m.rows ; k++)
        {
            vout[i] = vout[i] + m.m[ i * m.rows + k] * vin[k];
        } // for( k=0; k < m.rows ; k++)
    } // for( i=0 ; i<m.cols ; i++)

}

template <typename T>
T norm( const uVector< struct simpleComplex<T> > &m )
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

template <typename T>
T normsqrt( const uVector< struct simpleComplex<T> > &m )
{
    T v = 0.0;

    v = norm( m );

    return sqrt(v);
}

template <typename T>
struct simpleComplex<T> sum( const uVector< struct simpleComplex<T> > &m )
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

template <typename T>
void normalize( uVector< struct simpleComplex<T> > &m )
{

    T v;
    v = normsqrt( m );

    m = m / v;
    return ;
}

template <typename T>
void std_base_state( uVector< struct simpleComplex<T> > &m, int state )
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
uVector< struct simpleComplex<T> > dagnotdag( const uVector< struct simpleComplex<T> > &m )
{
    size_t i;

	uVector< struct simpleComplex<T> > r( SIZE1*SIZE2 );
	uMatrix< simpleComplex<T> > d(SIZE1, SIZE2),nd(SIZE1, SIZE2);
	
	//d.cols = SIZE1;
	//d.rows = SIZE2;	
	for(i=0;i<SIZE1*SIZE2;i++)
		d.m[i] = m[i];

	//nd.cols = SIZE1;
	//nd.rows = SIZE2;
	for(i=0;i<SIZE1*SIZE2;i++)
		nd.m[i] = m[i];
	
	d = dagger(nd);
	
	d  = d * nd;
	
	//r.size=SIZE1 * SIZE2;
	for(i=0;i<SIZE1*SIZE2;i++)
		r[i] = d.m[i];
	
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
uVector< struct simpleComplex<T> > coherent( struct simpleComplex<T> alpha )
{
	uVector< struct simpleComplex<T> > x(SIZE);
	uMatrix< struct simpleComplex<T> > a(SIZE), adag(SIZE), tmpmat(SIZE), D(SIZE);

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

template <typename T, size_t SIZE>
void coherent( simpleComplex<T> *tbl, struct simpleComplex<T> alpha)
{
	size_t i;
	
	uVector< struct simpleComplex<T> > r( SIZE );
	
	r = coherent<T, SIZE>(alpha);
	
    for(i=0; i<SIZE; i++)
    {
        tbl[i] = r[i];
    }	
}


//------------------------------------------------------------------------------------------------------------
// global output operator<< for simpleComplex<T>
//------------------------------------------------------------------------------------------------------------

template <typename T>
inline ostream& operator<< (ostream & output, const simpleComplex<T> & c)
{

#define width 5

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


template <typename T>
ostream& operator<< (ostream & output, const uMatrix<T> &c)
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


template <typename T>
ostream& operator<< (ostream & output, const uVector<T> &c)
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
template <typename T>
inline void print_matrix(const uMatrix<T> &A, string name)
{
  cout << name << "[" << A.rows << "x" << A.cols << "=" << A.size << "]:\n" << A << endl;
}

template <typename T>
bool inverse_of_matrix(const uMatrix<T> &M, uMatrix<T> &A_inv)
{
	// check sizes
	unsigned int size_m = M.size;
	unsigned int size_i = A_inv.size;
	if (size_m != size_i || M.rows != A_inv.rows || M.cols != A_inv.cols)
	{
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

  uMatrix<T> A(M.rows, M.cols); //create copy of matrix M as A, because matrix A change its values

  for(i=0;i<M.size;i++)
    A.m[i]=M.m[i];

  uMatrix<T> L(M.rows, M.cols), U(M.rows, M.cols);          // matrices
  uMatrix<T> b(M.rows, M.cols), d(M.rows, M.cols), x(M.rows, M.cols);  // vectors as matrices with one column size (n,1)

  for(i=0;i<b.size;i++)
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

template <typename T>
void matprod(uMatrix<T> &A, uMatrix<T> &B, uMatrix<T> &C);

template <typename T>	
void matcopy(uMatrix<T> &A, uMatrix<T> &B);

template <typename T>
static void matexp_pade(const long int p, uMatrix< simpleComplex<T> > &A, uMatrix< simpleComplex<T> > &N);


template <typename T>
void eye_of_matrix(uMatrix<T> &m)
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

template <typename T>
void zero_matrix( uMatrix< T > &m )
{
    size_t i;

    for(i=0; i<m.size; i++)
    {
       m.m[i] = make_simpleComplex(0.0, 0.0);
    }
}

template <typename T>
void pauli_x_matrix( uMatrix<T> &m )
{
    m.m[0] = make_simpleComplex(0.0, 0.0); m.m[1] = make_simpleComplex(1.0, 0.0);
    m.m[2] = make_simpleComplex(1.0, 0.0); m.m[3] = make_simpleComplex(0.0, 0.0);	
}

template <typename T>
void pauli_x_matrix( uVector<T> &v)
{
    v[0] = make_simpleComplex(0.0, 0.0); v[1] = make_simpleComplex(1.0, 0.0);
    v[2] = make_simpleComplex(1.0, 0.0); v[3] = make_simpleComplex(0.0, 0.0);	
}

template <typename T>
void pauli_y_matrix( uVector<T> &v)
{
    v[0] = make_simpleComplex(0.0, 0.0); v[1] = make_simpleComplex(0.0,-1.0);
    v[2] = make_simpleComplex(0.0, 1.0); v[3] = make_simpleComplex(0.0, 0.0);	
}

template <typename T>
void pauli_z_matrix( uVector<T> &v)
{
    v[0] = make_simpleComplex(1.0, 0.0); v[1] = make_simpleComplex( 0.0, 0.0);
    v[2] = make_simpleComplex(0.0, 0.0); v[3] = make_simpleComplex(-1.0, 0.0);	
}

template <typename T>
void transpose_of_matrix( uMatrix<T> &m )
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

template <typename T>
uMatrix<T> dagger( const uMatrix<T> &_m )
{	
	uMatrix<T> m(_m.rows, _m.cols);
	
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


template <typename T>
void destroy_operator( uMatrix<T> &m )
{
	size_t i, y;

	y=1;
	for(i=0; i < m.rows - 1; i++)
	{
		m(i, y)= make_simpleComplex( sqrt((double)y), 0.0 );
		y++;
	}
}

template <typename T>
void create_operator( uMatrix<T> &m )
{
	size_t i, y;

	y=1;
	for(i=0; i < m.rows - 1; i++)
	{
		m(y, i) = make_simpleComplex( sqrt((double)y), 0.0 );
		y++;
	}
}


template <typename T>
void make_x_operator( uMatrix<T> &m )
{
	size_t i, y;
	
	y=m.cols - 1;
	for(i=0; i < m.rows ; i++)
	{
		m(i, y) = make_simpleComplex( 1.0, 0.0 );
		y--;
	}
}

template <typename T>
void sigma_m_matrix( uMatrix<T> &v)
{
    v(0,0) = make_simpleComplex(0.0, 0.0); v(0,1) = make_simpleComplex(0.0, 0.0);
    v(1,0) = make_simpleComplex(1.0, 0.0); v(1,1) = make_simpleComplex(0.0, 0.0);	
}


template <typename T>
T norminf(const uMatrix<T> &m)
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

template <typename T>
T norminf( const uVector<T> &v)
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

template <typename T>
void zero_vector( uVector<T> &v)
{
    size_t i;

    for(i=0; i<v.size; i++)
    {
        v[i] = make_simpleComplex(0.0, 0.0);
    }
}

template <typename T>
uMatrix<simpleComplex<T> > tensor(uMatrix<simpleComplex<T> > &m1, uMatrix<simpleComplex<T> > &m2)
{
	int x,y,i,j,ii,jj;
	uMatrix< simpleComplex<T> >  tmp(m1.rows * m2.rows, m1.cols * m2.cols);
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

template <typename T>
uMatrix< simpleComplex<T> > tensor(uMatrix< simpleComplex<T> > &m1, uMatrix< simpleComplex<T> > &m2, uMatrix< simpleComplex<T> > &m3)
{
	uMatrix< simpleComplex<T> > tmp1(m1.rows * m2.rows, m1.cols * m2.cols);
	uMatrix< simpleComplex<T> > tmp2(m1.rows * m2.rows * m3.rows, m1.cols * m2.cols * m3.cols);
	
	tmp1=tensor<T>(m1, m2);
	tmp2=tensor<T>(tmp1, m3);
	
	return tmp2;
}

template <typename T, size_t SIZE1, size_t SIZE2>
void tensor(simpleComplex<T> *m1, simpleComplex<T> *m2, simpleComplex<T>  *m3)
{

	int x,y,i,j,ii,jj;
    simpleComplex<T> num;

	for(i=0;i<SIZE1*SIZE2;i++)
    {
		m3[i].re=0.0;
		m3[i].im=0.0;
	}

    ii=0;
    jj=0;
    for(x=0;x<SIZE1;x++)
    {
        for(y=0;y<1;y++)
        {
            ii=(x*SIZE2);
            jj=(y*1);
            for(i=0;i<SIZE2;i++)
            {
                for(j=0;j<1;j++)
                {
                    num.re=0.0;
                    num.im=0.0;

					num = m1[x] * m2[i];
				 
					m3[ii] = num;
					 
                    jj++;
                }
                ii++;
                jj-=1;
            }
        }
    }
}


#endif // __complexnum_h__

