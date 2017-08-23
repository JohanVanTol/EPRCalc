#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <string.h>
#include "Vector.h"
#include "Tensor.h"
#include "my_const.h"


class Spectrum
{
	private :
		Vector *Field;
		Vector *Intens;
		int npts;
		double start;
		double stop;


	public :
		Spectrum(int _npts, double _start=0, double stop=10.0);
		Spectrum(const Spectrum& S);
		~Spectrum();
		void addGauss(double x, double w, double A);
		void addDispersGauss(double x, double w, double A);
		void addGaussArray(Vector xc, Vector w, Vector A, int n);
		void addGaussArray(double *xc, double *w, double *A, int n);
		void addGaussDispersArray(Vector xc, Vector w, Vector A, int n);
		void addGaussDispersArray(double* xc, double* w, double* A, int n);
		void addVoigtArray(double alpha, double *xc, double *w, double *A, int n);
		void addVoigtDispersArray(double alpha, double *xc, double *w, double *A, int n);
		void addLorentz(double x, double w, double A);
		void addDispersLorentz(double x, double w, double A);
		void addLorentzArray(double *xc, double *w, double *A, int n);
		void addLorentzArray(Vector xc, Vector w, Vector A, int n);
        void addLorentzDispersArray(Vector xc, Vector w, Vector A, int n);
        void addLorentzDispersArray(double* xc, double* w, double* A, int n);
//		void write(ofstream* f);
		double getField(int i) const {return Field->get(i);}
		double getIntens(int i) const {return Intens->get(i);}
		void setField(int i, double field) {if (i<npts) Field->set(i,field);}
		void setIntens(int i, double intens) {if (i<npts) Intens->set(i,intens);}
        void addIntens(int i, double intens) {if (i<npts) Intens->add(i,intens);}
		int getNpts() const {return npts;}
		void Derivative(Spectrum S, double modulation);
		int Derivative(Spectrum);
		void operator=(const Spectrum S);
        double IntPol(double field, int *status) const;
		double IntPol(double field) const;
        double Mean(int order);
};


class MultSpectrum
{
	private :
		Vector *Field;
		Vector **Intens;
		int nI;
		int npts;
		double start;
		double stop;


	public :
		MultSpectrum(int ng, int _npts, double _start=0, double stop=10.0);
		MultSpectrum(const MultSpectrum& S);
		~MultSpectrum();
		void addGauss(int i, double x, double w, double A);
		void addGaussArray(int i, Vector xc, Vector w, Vector A, int n);
		void addLorentz(int i,double x, double w, double A);
		void addLorentzArray(int i, Vector xc, Vector w, Vector A, int n);
//		void write(ofstream* f);
		double getField(int i) const {return Field->get(i);}
		double getIntens(int ig, int i) const {return Intens[ig]->get(i);}
		void setField(int i, double field) {if (i<npts) Field->set(i,field);}
		void setIntens(int ig, int i, double intens) {if (i<npts) Intens[ig]->set(i,intens);}
		int getNpts() const {return npts;}
		int getNI() const {return nI;}
		void Derivative(MultSpectrum S, double modulation);
		int Derivative(MultSpectrum);
		void operator=(const MultSpectrum S);
		double IntPol(int ig, double field) const;
};
#endif