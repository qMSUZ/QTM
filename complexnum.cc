#include <math.h>

#include "complexnum.h"

template <typename T>
struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const struct simpleComplex<T> &b)
{
	struct simpleComplex<T> t;

	t.re = a.re + b.re;
	t.im = a.im + b.im;

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
struct simpleComplex<T> operator+(const struct simpleComplex<T> &a, const T &b)
{
	struct simpleComplex<T> t;

	t.re = a.re + b;
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
struct simpleComplex<T> make_simpleComplex (T r, T i)
{
	simpleComplex<T> t;

	t.re = r ;
	t.im = i ;

	return t;
}
