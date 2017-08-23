#include <iostream.h>
#include <complex.h>
#include "vector.h"
#include "tensor.h"
#include "eigen.h"

Tensor::Tensor(int n)
					: order(n)
{
	if (order < 1) return;
	a = new double*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new double[order];
		for (int j=0; j<order; j++) a[i][j] = 0.0 ;
	}
}

Tensor::Tensor(const Tensor& T)
{
	order = T.order;
	a = new double*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new double[order];
		for (int j=0; j<order; j++) a[i][j] = T.get(i,j);
	}
}

Tensor::~Tensor()
{
	for (int i=0; i<order;i++)
				delete [] a[i];
	delete [] a;
}

void Tensor::reset(int n)
{
	for (int i=0; i<order;i++)
				delete [] a[i];
	delete [] a;

	order = n;
	a = new double*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new double[order];
		for (int j=0; j<order; j++)
			a[i][j] = 0.0;
	}
}

Tensor& Tensor::operator=(const Tensor& T)
{
	if (this == &T) return *this;
	for (int i=0; i<order;i++)
				delete [] a[i];
	delete [] a;

	order = T.order;
	a = new double*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new double[order];
		for (int j=0; j<order; j++) a[i][j] = T.a[i][j];
	}
	return *this;
}

Tensor& Tensor::operator+=(const Tensor& T)
{
	for (int i=0; i<order;i++)
		for (int j=0; j<order; j++) a[i][j] += T.a[i][j];
	return *this;
}

Tensor& Tensor::operator-=(const Tensor& T)
{
	for (int i=0; i<order;i++)
		for (int j=0; j<order; j++) a[i][j] -= T.a[i][j];
	return *this;
}


Tensor Tensor::operator*(const Tensor& T2)
{
	Tensor mult = *this;
	if (order != T2.order) return *this;
	for (int i=0; i<order; i++)
		for (int j=0; j<order; j++)
		{
			mult.a[i][j] = 0;
			for (int k=0; k<order ; k++)
				mult.a[i][j] += a[i][k] * T2.a[k][j];
		};
	return mult;
}

Vector Tensor::operator*(const Vector& V) const
{
	if (V.order() != order) return V;   // do nothing if orders are not right
	Vector New(order);
	double temp;
	for (int i=0; i<order; i++)
	{
		temp = 0.0;
		for (int j=0; j<order ; j++)
			temp += a[i][j] * V.get(j);
		New.set(i, temp);
	}
	return New;
}

Tensor Tensor::operator+(const Tensor& T2)
{
	Tensor sum = *this;
	if (order != T2.order) return *this;
	sum += T2;
	return sum;
}

Tensor& Tensor::operator*=(double c)
{
	for (int i=0; i< order; i++)
		for (int j=0; j< order; j++)
			a[i][j] *= c;
	return *this;
}

Tensor Tensor::operator*(double c)
{
	Tensor temp = *this;
	temp *= c;
	return temp;
}

Tensor Tensor::operator-(const Tensor& T2)
{
	Tensor diff = *this;
	if (order != T2.order) return *this;
	diff -= T2;
	return diff;
}

void Tensor::set(int i, int j, double r)
{
	if ((i>=order) || (j>=order))
	{
		cout << "WARNING index out of range !" << endl;
		return;
	}
	a[i][j] = r;
}

void Tensor::add(int i, int j, double r)
{
	if ((i>=order) || (j>=order)) return;
	a[i][j] += r;
}

//double Tensor::get(int i, int j) const
/*{
	if ((i>=order) || (j>=order)) return 0.0;
	return a[i][j];
}
*/

void Tensor::print() const
{
	if (order < 1) return;

	for (int i=0; i<order;i++)
	{
		cout << "( " << a[i][0];
		for (int j=1; j<order; j++)
			cout << ", " << a[i][j];
		cout << " )" << endl;
	}
}

Tensor Tensor::transpose() const
{
	int i,j;

	Tensor T(order);
	for (i=0; i<order; i++)
		for (j=0; j<order; j++)
			T.set(i,j, a[j][i]);
	return T;
}

Vector Tensor::row(int k) const
{
	Vector V(order);
	if ((k<0) || (k>=order)) return V;
	for (int i=0; i<order; i++)
		V.set(i,a[k][i]);
	return V;
}

Vector Tensor::col(int k) const
{
	Vector V(order);
	if ((k<0) || (k>=order)) return V;
	for (int i=0; i<order; i++)
		V.set(i,a[i][k]);
	return V;
}

void Tensor::GetCopy(int _order, double **M)
{
  if (order != _order) return;
  for (int i=0; i<order;i++)
    for (int j=0; j<order; j++)
      M[i][j] = a[i][j];
  return;
}

int Tensor::GetPrincipleAxes(int n, double *PV, double *Theta, double *Phi)
{
  if (order != n) return 0;
	double **A;
	double *d,*e;

	int i, order = n;
	A = new double*[order];
	for (i=0; i<order; i++)
		A[i] = new double[order];
	d = new double[order];
	e = new double[order];

  GetCopy(order,A);

	tred2(A,order, d, e);
	tqli(d,e,order, A);
  eigen_sort(d,A,n);

  double Rad = 90.0/asin(1.0);

  for (i=0; i<order;i++)
  {
     PV[i] = d[i];
     Theta[i] = acos(A[2][i]) * Rad;
     Phi[i] = 0.0;
     if (A[0][i] > 1e-12)
      Phi[i] = Rad * atan(A[1][i]/A[0][i]);
      else
        if (A[0][i] < -1e-12)
            Phi[i] = 180 + Rad * atan(A[1][i]/A[0][i]);
          else
          {
          if (A[1][i]<-1e-12) Phi[i] = 270;
            else
              if (A[1][i] > 1e-12) Phi[i] = 90;
              else Phi[i] = 0.0;
          }
   }

	for (i=0; i<order; i++)
        delete[] A[i];
    delete[] A;
    delete[] e;
    delete[] d;

   return n;
}


//------------------------------------------------------------------------
Ctensor::Ctensor(int n)
					: order(n)
{
	a = new complex*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new complex[order];
		for (int j=0; j<order; j++) a[i][j] = 0.0, 0.0 ;
	}
}

Ctensor::Ctensor(const Ctensor& T)
{
	order = T.order;
	a = new complex*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new complex[order];
		for (int j=0; j<order; j++) a[i][j] = T.a[i][j];
	}
}

Ctensor::~Ctensor()
{
	for (int i=0; i<order;i++)
				delete [] a[i];
	delete [] a;
}

void Ctensor::reset(int n)
{
	for (int i=0; i<order;i++)
				delete [] a[i];
	delete [] a;

	order = n;
	a = new complex*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new complex[order];
		for (int j=0; j<order; j++)
			a[i][j] = 0.0;
	}
}

Ctensor& Ctensor::operator=(const Ctensor& T)
{
	if (this == &T) return *this;
	for (int i=0; i<order;i++)
				delete [] a[i];
	delete [] a;

	order = T.order;
	a = new complex*[order];
	for (int i=0; i<order;i++)
	{
		a[i] = new complex[order];
		for (int j=0; j<order; j++) a[i][j] = T.a[i][j];
	}
	return *this;
}

Ctensor& Ctensor::operator+=(const Ctensor& T)
{
	for (int i=0; i<order;i++)
		for (int j=0; j<order; j++) a[i][j] += T.a[i][j];
	return *this;
}

Ctensor& Ctensor::operator-=(const Ctensor& T)
{
	for (int i=0; i<order;i++)
		for (int j=0; j<order; j++) a[i][j] -= T.a[i][j];
	return *this;
}

void Ctensor::print() const
{
	if (order<1) return;
	for (int i=0; i<order;i++)
	{
		for (int j=0; j<order; j++)
			cout << "  " << a[i][j];
		cout << endl;
	}
}

Ctensor Ctensor::operator*(const Ctensor& T2)
{
	Ctensor mult = *this;
	if (order != T2.order) return *this;
	for (int i=0; i<order; i++)
		for (int j=0; j<order; j++)
		{
			mult.a[i][j] = 0;
			for (int k=0; k<order ; k++)
				mult.a[i][j] += a[i][k] * T2.a[k][j];
		};
	return mult;
}

Cvector Ctensor::operator*(const Cvector& V) const
{
	if (V.order() != order) return V;   // do nothing if orders are not right
	Cvector New(order);
	complex temp;
	for (int i=0; i<order; i++)
	{
		temp = 0.0;
		for (int j=0; j<order ; j++)
			temp += a[i][j] * V.get(j);
		New.set(i, temp);
	}
	return New;
}

Ctensor Ctensor::operator+(const Ctensor& T2)
{
	Ctensor sum = *this;
	if (order != T2.order) return *this;
	sum += T2;
	return sum;
}

Ctensor& Ctensor::operator*=(complex c)
{
	for (int i=0; i< order; i++)
		for (int j=0; j< order; j++)
			a[i][j] *= c;
	return *this;
}

Ctensor Ctensor::operator*(complex c)
{
	Ctensor temp = *this;
	temp *= c;
	return temp;
}

Ctensor Ctensor::operator-(const Ctensor& T2)
{
	Ctensor diff = *this;
	if (order != T2.order) return *this;
	diff -= T2;
	return diff;
}

void Ctensor::set(int i, int j, complex c)
{
	if ((i>=order) || (j>=order)) return;
	a[i][j] = c;
}

void Ctensor::set(int i, int j, double re, double im)
{
	if ((i>=order) || (j>=order)) return;
	a[i][j] = complex(re,im);
}

void Ctensor::add(int i, int j, complex c)
{
	if ((i>=order) || (j>=order)) return;
	a[i][j] += c;
}

void Ctensor::add(int i, int j, double re, double im)
{
	if ((i>=order) || (j>=order)) return;
	a[i][j] += complex(re, im);
}

complex Ctensor::get(int i, int j) const
{
	if ((i>=order) || (j>=order)) return 0.0;
	return a[i][j];
}

Ctensor Ctensor::transpose() const
{
	int i,j;

	Ctensor T(order);
	for (i=0; i<order; i++)
		for (j=0; j<order; j++)
			T.set(i,j, conj(a[j][i]));
	return T;
}

Tensor Ctensor::realtensor()  const
{
	Tensor T(order);
	for (int i=0; i<order; i++)
		for (int j=0; j<order; j++)
			T.set(i,j,real(a[i][j]));
	return T;
}

Tensor Ctensor::imagtensor()  const
{
	Tensor T(order);
	for (int i=0; i<order; i++)
		for (int j=0; j<order; j++)
			T.set(i,j,imag(a[i][j]));
	return T;
}


Cvector Ctensor::row(int k) const
{
	Cvector V(order);
	if ((k<0) || (k>=order)) return V;
	for (int i=0; i<order; i++)
		V.set(i,a[k][i]);
	return V;
}

Cvector Ctensor::col(int k) const
{
	Cvector V(order);
	if ((k<0) || (k>=order)) return V;
	for (int i=0; i<order; i++)
		V.set(i,a[i][k]);
	return V;
}
