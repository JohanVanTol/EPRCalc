#include <iostream.h>
#include <math.h>
#include <complex.h>
#include "vector.h"

Vector::Vector(int _n)
			:n(_n)
{
	e = new double[n];
	for (int i=0; i<n; i++) e[i] = 0.0;
}

Vector::Vector(const Vector& _v)
{
	n = _v.n;
	e = new double[n];
	for (int i=0; i<n; i++) e[i] = _v.e[i];
}

void Vector::reset(int _n)
{
	if (n != _n)
	{
		delete[] e;
		e = new double[_n];
		n = _n;
	}
	for (int i=0; i<n; i++) e[i] = 0.0;
	return;
}

void Vector::add(int i, double d)
{
    if (d>1e300)
    {
        d=0.0;
        }
    if (i<n) e[i] += d;
}

void Vector::print() const
{
	if ( n < 1 ) return;

	cout << "( " ;
	for (int i=0; i<n-1; i++)
		cout << e[i] << ", ";
	cout << e[n-1] << ")" << endl;
}

double Vector::length() const
{
	double l=0;
	for (int i=0; i<n ;i++)
		l += e[i]*e[i];
	return sqrt(l);
}
/*
// Tensor Vector::matrix(const Vector &V) const
// {
	if (n != V.order()) cout << "Error, vectors not same order" ;
	Tensor T(n);
	for (int i=0; i<n;i++)
		for (int j=0; j<n; j++)
			T.set(i,j,e[i]*V.get(j));
	return T;
}
*/

Vector& Vector::operator=(const Vector& _vec)
{
	if ( this == &_vec ) return *this;

	delete [] e;
	n = _vec.n;
	e = new double[n];
	for (int i=0; i<n; i++)
		e[i] = _vec.e[i];

	return *this;
}

double Vector::operator[](int i)
{
	if ((i < 0) || (i >= n)) return 0.0;
	return e[i];
}

Vector& Vector::operator*=(double a)
{
	for (int i = 0; i<n; i++)
			e[i] *= a;
	return *this;
}

double Vector::operator*(const Vector& v) const
{
	double l=0;
	for (int i=0; i<n; i++)
			l += e[i] * v.e[i];
	return l;
}

Vector Vector::operator*(double a) const
{
	Vector temp = *this;
	temp *= a;
	return temp;
}

Vector& Vector::operator+=(const Vector& v)
{
	if ( n != v.n ) return *this;
	for (int i = 0; i<n; i++) e[i] += v.e[i];
	return *this;
}

Vector Vector::operator+(const Vector& v) const
{
	Vector temp = *this;
	temp += v;
	return temp;
}

Vector& Vector::operator-=(const Vector& v)
{
	if ( n != v.n ) return *this;
	for (int i = 0; i<n; i++) e[i] -= v.e[i];
	return *this;
}

Vector Vector::operator-(const Vector& v) const
{
	Vector temp = *this;
	temp -= v;
	return temp;
}

Vec3::Vec3()
			:Vector(3)
{}

void Vec3::setv(double x, double y, double z)
{
	set(0,x);
	set(1,y);
	set(2,z);
}

//----------------------------------------------------------
Cvector::Cvector(int _n)
			:n(_n)
{
	e = new complex[n];
	for (int i=0; i<n; i++) e[i] = complex(0.0, 0.0);
}

Cvector::Cvector(const Cvector& _v)
{
	n = _v.n;
	e = new complex[n];
	for (int i=0; i<n; i++) e[i] = _v.e[i];
}

void Cvector::reset(int _n)
{
	if (n != _n)
	{
		delete [] e;
		e = new complex[_n];
		n = _n;
	}
	for (int i=0; i<n; i++) e[i] = (0.0, 0.0);
	return;
}


void Cvector::print() const
{
	if ( n < 1 ) return;

	cout << "( " ;
	for (int i=0; i<n-1; i++)
		cout << e[i] << ", ";
	cout << e[n-1] << ")" << endl;
}

double Cvector::length() const
{
	double l=0;
	for (int i=0; i<n ;i++)
		l += real(conj(e[i])*e[i]);
	return sqrt(l);
}

complex Cvector::rephase()
{
	// first get the largest coefficient
	int i, index=0;
	double largest = norm(e[0]);
	complex factor = 0.0;
	for (i=1; i<n; i++)
		if (norm(e[i]) > largest)
		{
			largest = norm(e[i]);
			index = i;
		}
	if (largest > 0.0)
	{
		factor = conj(e[index]) / sqrt(largest);
		for (i=0; i<n; i++) e[i] *= factor;
	}
		else cout << "Smaller " << endl;

	return factor;
}


Cvector Cvector::conjugate() const
{
	Cvector temp(n);
	for (int i=0; i<n; i++) temp.set(i, conj(e[i]));
	return temp;
}

Cvector& Cvector::operator=(const Cvector& _vec)
{
	if ( this == &_vec ) return *this;

	delete [] e;
	n = _vec.n;
	e = new complex[n];
	for (int i=0; i<n; i++)
		e[i] = _vec.e[i];

	return *this;
}

complex Cvector::operator[](int i)
{
	if ((i < 0) || (i >= n)) return 0.0;
	return e[i];
}

Cvector& Cvector::operator*=(complex a)
{
	for (int i = 0; i<n; i++)
			e[i] *= a;
	return *this;
}

Cvector& Cvector::operator*=(double a)
{
	for (int i = 0; i<n; i++)
			e[i] *= a;
	return *this;
}

complex Cvector::operator*(const Cvector& v) const
{
	complex l=0;
	for (int i=0; i<n; i++)
			l += conj(e[i]) * v.e[i];
	return l;
}

Cvector Cvector::operator*(double a) const
{
	Cvector temp = *this;
	temp *= a;
	return temp;
}

Cvector Cvector::operator*(complex a) const
{
	Cvector temp = *this;
	temp *= a;
	return temp;
}

Cvector& Cvector::operator+=(const Cvector& v)
{
	if ( n != v.n ) return *this;
	for (int i = 0; i<n; i++) e[i] += v.e[i];
	return *this;
}

Cvector& Cvector::operator-=(const Cvector& v)
{
	if ( n != v.n ) return *this;
	for (int i = 0; i<n; i++) e[i] -= v.e[i];
	return *this;
}

Cvector Cvector::operator+(const Cvector& v) const
{
	Cvector temp = *this;
	temp += v;
	return temp;
}

Cvector Cvector::operator-(const Cvector& v) const
{
	Cvector temp = *this;
	temp -= v;
	return temp;
}

