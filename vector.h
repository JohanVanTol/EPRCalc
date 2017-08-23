#include <iostream.h>
#include <math.h>
#include <complex.h>


#ifndef VectorH
#define VectorH

#define complex complex<double>
class Vector
{
	private:
		double *e;
		int n;
	public:
		Vector(int _n=3);
		Vector(const Vector& _vec);
		~Vector() {delete [] e;}
		void reset(int n);
		int order() const {return n;}
		void set(int i, double d) {/* if (i<n) */ e[i] = d;}
		double get(int i) const {if (i<n) return e[i];
											 else return 0.0;}
		void add(int i, double d);
//		Tensor matrix(const Vector &V) const;
		void print() const;
		double length() const;
		Vector& operator=(const Vector& _v);
		Vector& operator*=(double a);
		double operator*(const Vector& _v) const;
		Vector operator*(double a) const;
		Vector& operator+=(const Vector& v);
		Vector& operator-=(const Vector& v);
		Vector operator+(const Vector& v) const;
		Vector operator-(const Vector& v) const;
		double operator[](int i);
};

class Cvector
{
	private:
		complex *e;
		int n;
	public:
		Cvector(int _n=2);
		Cvector(const Cvector& _vec);
		~Cvector() {delete [] e;}
		void reset(int n);
		int order() const {return n;}
		void set(int i, complex d) {if (i<n) e[i] = d;}
		void set(int i, double re, double im) {if (i<n) e[i] = (complex(re,im));}

		complex get(int i) const {if (i<n) return e[i];
											 else return 0.0;}

      complex rephase();

		Cvector conjugate() const;
		void print() const;
		double length() const;
		Cvector& operator=(const Cvector& _v);
		Cvector& operator*=(double a);
		Cvector& operator*=(complex a);
		complex operator*(const Cvector& _v) const;
		Cvector operator*(double a) const;
		Cvector operator*(complex a) const;
		Cvector& operator+=(const Cvector& v);
		Cvector& operator-=(const Cvector& v);
		Cvector operator+(const Cvector& v) const;
		Cvector operator-(const Cvector& v) const;
		complex operator[](int i);
};

class Vec3 : public Vector
{
	public:
		Vec3();
		void setv(double x, double y, double z);
};

#endif
