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
#include "Spectrum.h"
#include "Gaussian2.h"
#include "Lorentz.h"


Spectrum::Spectrum(int _npts, double _start, double _stop)
		: npts(_npts), start(_start), stop(_stop)
{
	Field = new Vector(npts);
	Intens = new Vector(npts);
	for (int i=0; i<npts; i++)
		Field->set(i,start + (stop-start)*double(i)/double(npts-1));
}

Spectrum::Spectrum(const Spectrum& S)
{
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		Intens->set(i,S.getIntens(i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}

Spectrum::~Spectrum()
{
	delete Field;
	delete Intens;
}

void Spectrum::operator=(const Spectrum S)
{
	delete Field;
	delete Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		Intens->set(i,S.getIntens(i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}


//void Spectrum::addGauss(double xc, double w, double A)
//{
//	for (int i=0; i<npts; i++)
//	{
//		double x = 2.0*(Field->get(i)-xc)/w ;  // 2(x-x0)/w
//		Intens->add(i, 2.0*(A/w)*sqrt(log(2.0L)/PI) * exp(-log(2.0L)*x*x));
//	}
//}

void Spectrum::addGauss(double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
//		double x = 2.0*(Field->get(i)-xc)/w ;  // 2(x-x0)/w
		Intens->add(i, NormGauss(Field->get(i),xc,A,w));
	}
}

void Spectrum::addDispersGauss(double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
//		double x = 2.0*(Field->get(i)-xc)/w ;  // 2(x-x0)/w
		Intens->add(i, NormDispGaussian(Field->get(i),xc,A,w));
	}
}

void Spectrum::addLorentz(double xc, double w, double A)
{
	double x;
	for (int i=0; i<npts; i++)
	{
		x = 2*(Field->get(i)-xc)/w;
		Intens->add(i, (2*A)/(PI*w) * 1.0/(1.0+x*x));
	}
}

void Spectrum::addDispersLorentz(double xc, double w, double A)
{
	double x;
	for (int i=0; i<npts; i++)
	{
		x = 2*(Field->get(i)-xc)/w;
		Intens->add(i, (2*A)/(PI*w) * x/(1.0+x*x));
	}
}

void Spectrum::addGaussArray(Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double xf,x,ww;
	for (i=0; i<npts; i++)
	{
		xf = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			x = 2.0*(xf-xc.get(j))/ww;
			Intens->add(i, (2.0*A.get(j)/ww) * sqrt(log(2.0)/PI) * exp(-log(2.0)*x*x));
		}
	}
}

/*void Spectrum::addGaussArray(double *xc, double *w, double *A, int n)
{
	int i,j;
	double x;
	for (i=0; i<npts; i++)
	{
		for (j=0; j<n; j++)
		{
			x = 2.0*(Field->get(i)-xc[j])/w[j];
			if (fabs(A[j])>1e-100) Intens->add(i, (2.0*A[j]/w[j])*sqrt(log(2.0)/PI) * exp(-log(2.0)*x*x));
		}
	}
}
*/

void Spectrum::addGaussArray(double *xc, double *w, double *A, int n)
{
	for (int j=0; j<n; j++)
	{
		addGauss(xc[j],w[j],A[j]);
	}
}

void Spectrum::addGaussDispersArray(Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww,xf,x;
	for (i=0; i<npts; i++)
	{
		double xf = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			x = 2.0*(xf-xc.get(j))/ww;
			Intens->add(i, 2.0*A.get(j)/(ww*PI) * x/(1.0+x*x));
		}
	}
}

void Spectrum::addGaussDispersArray(double* xc, double* w, double* A, int n)
{
	for (int j=0; j<n; j++)
	{
		addDispersGauss(xc[j],w[j],A[j]);
	}
}

void Spectrum::addVoigtArray(double alpha, double *xc, double *w, double *A, int n)
{
	for (int i=0; i<n; i++)
	{
		addLorentz(xc[i],w[i],A[i]*alpha);
		addGauss(xc[i],w[i],A[i]*(1.0-alpha));
	}
}

void Spectrum::addVoigtDispersArray(double alpha, double *xc, double *w, double *A, int n)
{
	for (int i=0; i<n; i++)
	{
		addDispersLorentz(xc[i],w[i],A[i]*alpha);
		addDispersGauss(xc[i],w[i],A[i]*(1.0-alpha));
	}
}


void Spectrum::addLorentzArray(Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww,xf,x;
	for (i=0; i<npts; i++)
	{
		double xf = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			x = 2.0*(xf-xc.get(j))/ww;
			Intens->add(i, 2.0*A.get(j)/(ww*PI) * 1.0/(1.0+x*x));
		}
	}
}

void Spectrum::addLorentzArray(double *xc, double *w, double *A, int n)
{
	int i,j;
	double xf,x,ww;
	for (i=0; i<npts; i++)
	{
		xf = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w[j];
			x = 2.0*(xf-xc[j])/ww;
			Intens->add(i, 2.0*A[j]/(ww*PI) * 1.0/(1.0+x*x));
		}
	}
}

void Spectrum::addLorentzDispersArray(Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww,xf,x;
	for (i=0; i<npts; i++)
	{
		double xf = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			x = 2.0*(xf-xc.get(j))/ww;
			Intens->add(i, 2.0*A.get(j)/(ww*PI) * x/(1.0+x*x));
		}
	}
}

void Spectrum::addLorentzDispersArray(double* xc, double* w, double* A, int n)
{
	int i,j;
	double xf,x,ww;
	for (i=0; i<npts; i++)
	{
		xf = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w[j];
			x = 2.0*(xf-xc[j])/ww;
			Intens->add(i, 2.0*A[j]/(ww*PI) * x/(1.0+x*x));
		}
	}
}

int Spectrum::Derivative(Spectrum S)
{
	int i;
	delete Field;
	delete Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	Field->set(0,S.getField(0));
	Intens->set(0,S.getIntens(1)-S.getIntens(0));
	for (i=1; i<npts-1; i++)
	{
		Field->set(i,S.getField(i));
		Intens->set(i,(S.getIntens(i+1)-S.getIntens(i-1))/2.0);
	}
	Field->set(npts-1,S.getField(npts-1));
	Intens->set(0,S.getIntens(npts-1)-S.getIntens(npts-2));
	start = Field->get(0);
	stop = Field->get(npts-1);
	return npts;
}


void Spectrum::Derivative(Spectrum S, double mod)
{
	delete Field;
	delete Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	mod /=1000;    // from mT to Tesla
	if (fabs(mod) < 1.0e-7) mod = 1.0e-7;
	double B;
	double lval, hval;
	int lstat;
	int hstat;
	for (int i=0; i<npts; i++)
	{
		B = S.getField(i);
		Field->set(i,B);
		lval = S.IntPol(B-mod/2.0, &lstat);
		hval = S.IntPol(B+mod/2.0, &hstat);
		if (lstat < 0)
			setIntens(i,2*(hval - S.getIntens(i)));
		else
		{
			if (hstat < 0)
			   setIntens(i,2*(S.getIntens(i) - lval));
		   else
				setIntens(i,hval - lval);
		}
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
	return;
}



double Spectrum::IntPol(double field) const
{
	if ((field < Field->get(0)) || (field > Field->get(npts-1))) return -1.0;
			// Return negative value if outside range
	int i=0;
	while ((Field->get(i) < field) && (i<npts)) i++;
	if (i >= npts) return -1; //Outside of range
	double q = (field - Field->get(i-1))/(Field->get(i)-Field->get(i-1));
	return Intens->get(i-1) + q * (Intens->get(i) - Intens->get(i-1));
}

double Spectrum::IntPol(double field, int *stat) const
{
	*stat = 0;
	if (field < Field->get(0))
	{
		*stat = -1;
		return Intens->get(0);
	}
	if (field > Field->get(npts-1))
	{
		*stat = -2;
		return Intens->get(npts-1);
	}

	int i=0;
	while ((Field->get(i) < field) && (i<npts)) i++;
//	if (i >= npts) return -1; //Outside of range
	double q = (field - Field->get(i-1))/(Field->get(i)-Field->get(i-1));
	return Intens->get(i-1) + q * (Intens->get(i) - Intens->get(i-1));
}

double Spectrum::Mean(int order)
{
	double Mn = 0;
	double temp = 1;
	for (int i=0; i<npts; i++)
	{
		for (int j=0; j<order; j++)
		   temp *= Intens->get(i);
		Mn += temp;
	}
	if (npts >0)
	  Mn /= (double)npts;
	return Mn;
}



MultSpectrum::MultSpectrum(int _nI, int _npts, double _start, double _stop)
		: nI(_nI), npts(_npts), start(_start), stop(_stop)
{
	Field = new Vector(npts);
	Intens = new Vector*[nI];
	for (int ig=0; ig<nI; ig++)
		Intens[ig] = new Vector(npts);
	for (int i=0; i<npts; i++)
		Field->set(i,start + (stop-start)*double(i)/double(npts-1));
}

MultSpectrum::MultSpectrum(const MultSpectrum& S)
{
	npts = S.getNpts();
	nI = S.getNI();
	Field = new Vector(npts);
	Intens = new Vector*[nI];
	for (int ig=0; ig<nI; ig++)
		Intens[ig] = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		for (int ig=0; ig < nI; ig++)
			Intens[ig]->set(i,S.getIntens(ig,i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}

MultSpectrum::~MultSpectrum()
{
	delete Field;
	for (int ig=0; ig<nI; ig++)
		delete Intens[ig];
	delete[] Intens;
}

void MultSpectrum::operator=(const MultSpectrum S)
{
	int ig;
	delete Field;
	for (ig=0; ig<nI; ig++)
		delete Intens[ig];
	delete[] Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector*[nI];
	for (ig=0; ig<nI; ig++)
		Intens[ig] = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		for (ig=0; ig<nI; ig++)
			Intens[ig]->set(i,S.getIntens(ig, i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}


void MultSpectrum::addGauss(int ig, double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
		double x = 2.0*(Field->get(i)-xc)/w;
		Intens[ig]->add(i, (2.0*A/w) * sqrt(log(2.0)/PI)*exp(-log(2.0*x*x)));
	}
}

void MultSpectrum::addLorentz(int ig, double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
		double x = 2.0*(Field->get(i)-xc)/w;
		Intens[ig]->add(i, (2*A)/(PI*w) * 1.0/(1.0+x*x));
	}
}

void MultSpectrum::addGaussArray(int ig, Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double x;
	for (i=0; i<npts; i++)
	{
		x = Field->get(i);
		for (j=0; j<n; j++)
		{
			Intens[ig]->add(i, NormGauss(x,xc.get(j),A.get(j), w.get(j)));
		}
	}
}

void MultSpectrum::addLorentzArray(int ig, Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww;
	for (i=0; i<npts; i++)
	{
		double x = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			Intens[ig]->add(i, 2.0*A.get(j)*ww/(PI *(4*(x-xc.get(j))*(x-xc.get(j))+(ww*ww))));
		}
	}
}

int MultSpectrum::Derivative(MultSpectrum S)
{
	int i, ig;
	if ((npts!= S.getNpts()) || (nI != S.getNI()))
	{
		delete Field;
		for (ig=0; ig<nI; ig++)
			delete Intens[ig];
		delete[] Intens;
		npts = S.getNpts();
		nI = S.getNI();
		Field = new Vector(npts);
		Intens = new Vector*[nI];
		for (ig=0; ig<nI; ig++)
			Intens[ig] = new Vector(npts);
	}
	Field->set(0,S.getField(0));
	for (ig=0; ig<nI; ig++)
		Intens[ig]->set(0,S.getIntens(ig,1)-S.getIntens(ig,0));
	for (i=1; i<npts-1; i++)
	{
		Field->set(i,S.getField(i));
		for (ig=0; ig<nI; ig++)
			Intens[ig]->set(i,(S.getIntens(ig,i+1)-S.getIntens(ig,i-1))/2.0);
	}
	Field->set(npts-1,S.getField(npts-1));
	for (ig=0; ig<nI; ig++)
		Intens[ig]->set(0,S.getIntens(ig, npts-1)-S.getIntens(ig, npts-2));
	start = Field->get(0);
	stop = Field->get(npts-1);
	return npts;
}




void MultSpectrum::Derivative(MultSpectrum S, double mod)
{
	int ig;
	if ((npts!= S.getNpts()) || (nI != S.getNI()))
	{
		delete Field;
		for (ig=0; ig<nI; ig++)
			delete Intens[ig];
		delete[] Intens;
		npts = S.getNpts();
		nI = S.getNI();
		Field = new Vector(npts);
		Intens = new Vector*[nI];
		for (ig=0; ig<nI; ig++)
			Intens[ig] = new Vector(npts);
	}
	mod /=1000;
	double B;
	for (int i=0; i<npts; i++)
	{
		B = S.getField(i);
		Field->set(i,B);
		if (S.IntPol(ig, B-mod/2.0) < -0.5)
		for (ig=0; ig<nI; ig++)
			setIntens(ig,i,2*(S.IntPol(ig,B+mod/2.0) - S.getIntens(ig,i)));
		else
		{
			if (S.IntPol(ig, B+mod/2) < -0.5)
				for (ig=0; ig<nI; ig++)
					setIntens(ig,i,2*(S.getIntens(ig,i) - S.IntPol(ig,B-mod/2.0)));
				else
				for (ig=0; ig<nI; ig++)
					setIntens(ig,i,S.IntPol(ig,B+mod/2.0) - S.IntPol(ig, B-mod/2.0));
		}
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
	return;
}



double MultSpectrum::IntPol(int ig, double field) const
{
	if ((field < Field->get(0)) || (field > Field->get(npts-1))) return -1.0;
			// Return negative value if outside range
	int i=0;
	while ((Field->get(i) < field) && (i<npts)) i++;
	if (i >= npts) return -1; //Outside of range
	double q = (field - Field->get(i-1))/(Field->get(i)-Field->get(i-1));
	return Intens[ig]->get(i-1) + q * (Intens[ig]->get(i) - Intens[ig]->get(i-1));
}


	//void Spectrum::write(ofstream* f)
//{
//	for (int i =0; i<npts; i++)
//	*f << Field->get(i) << " " << Intens->get(i) << endl;
//}


