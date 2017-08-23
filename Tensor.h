
#ifndef TensorH
#define TensorH
#include <complex.h>
#include "vector.h"
class Tensor
{
	protected:
			double **a;
			int order;
	public:
			Tensor(int n=2);
			Tensor(const Tensor& _t);
			~Tensor();
			void reset(int n);

			Tensor& operator=(const Tensor& T);
			Tensor& operator*=(double r);
			Tensor& operator+=(const Tensor& T);
			Tensor& operator-=(const Tensor& T);
			Tensor operator+(const Tensor& T);
			Tensor operator-(const Tensor& T);
			Tensor operator*(const Tensor& T);
			Tensor operator*(double r);
			Vector operator*(const Vector& V) const;

			Tensor transpose() const;
			Vector row(int k) const;
			Vector col(int k) const;
			void set(int i, int j, double r);
			void add(int i, int j, double r);
			double get(int i, int j) const {return a[i][j];} 

			int get_order() const {return order;}


			void print() const;

//   New routines, which only apply to a symmetric 3x3 tensor
    void GetCopy(int _order, double **M);
    int GetPrincipleAxes(int n, double *PV, double *Theta, double *Phi);

};

class Ctensor
{
	protected:
			complex **a;
			int order;
	public:
			Ctensor(int n=3);
			Ctensor(const Ctensor& T);
			~Ctensor();
			void reset(int n);

			Ctensor& operator=(const Ctensor& T);
			Ctensor& operator*=(complex c);
			Ctensor& operator+=(const Ctensor& T);
			Ctensor& operator-=(const Ctensor& T);
			Ctensor operator+(const Ctensor& T);
			Ctensor operator-(const Ctensor& T);
			Ctensor operator*(const Ctensor& T);
			Ctensor operator*(complex c);
			Cvector operator*(const Cvector& V) const;

			Cvector row(int k) const;
			Cvector col(int k) const;

			void set(int i, int j, complex c);
			void set(int i, int j, double re, double im);
			void add(int i, int j, complex c);
			void add(int i, int j, double re, double im);
			complex get(int i, int j) const;

			Ctensor transpose() const;
			Tensor realtensor() const;
			Tensor imagtensor() const;
			int get_order() const {return order;}
			void print() const;
};
#endif
