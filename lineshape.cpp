//---------------------------------------------------------------------------


#pragma hdrstop

#include "Lineshape.h"
#include "Gaussian2.h"
#include "Lorentz.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

//Lineshapes package

double PseudoVoigt(double alpha, double x, double xc, double A, double w)
{
	return (1.0-alpha)*NormGauss(x,xc,A,w)+ alpha*NormLorentz(x,xc,A,w);
}

double DispPseudoVoigt(double alpha, double x, double xc, double A, double w)
{
	return (1.0-alpha)*NormDispGaussian(x,xc,A,w)+ alpha*NormDispLorentz(x,xc,A,w);
}

void MultPhasedSplitPseudoVoigt(double x, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
	*y = a[0] + a[1]*x;
	double yL = 0.0;
	double anglederiv = 0.0;
	dyda[na-1] = 0.0;
	for (int i=2; i<na-3; i+=7)
	{
		PhasedSplitPseudoVoigt(x,a+i,&yL,dyda+i,7);
		*y += yL;
	}
	dyda[0] = 1.0;
	dyda[1] = x;
}

void MultDerivPhasedSplitPseudoVoigt(double x, double *a, double *y, double *dyda, int na)
{
	*y = a[0] + a[1]*x;
	double yL = 0.0;
	for (int i=2; i<na-3; i+=7)
	{
		DerivPhasedSplitPseudoVoigt(x,a+i,&yL,dyda+i,7);
		*y += yL;
	}
	dyda[0] = 1.0;
	dyda[1] = x;
}

void MultDerivPhasedSplitPseudoVoigt(double *x, int nx, double *a, double *y, double *dyda, int na)
{
	*y = a[0] + a[1]*x[0];
	double yL = 0.0;
	for (int i=2; i<na-6; i+=7)
	{
		DerivPhasedSplitPseudoVoigt(x[0],a+i,&yL,dyda+i,7);
		*y += yL;
	}
	dyda[0] = 1.0;
	dyda[1] = x[0];

}

void MultPhasedSplitPseudoVoigt(double *x, int nx, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
	*y = a[0] + a[1]*x[0];
	double yL = 0.0;
	dyda[na-1] = 0.0;
	for (int i=2; i<na-6; i+=7)
	{
		PhasedSplitPseudoVoigt(x[0] ,a+i,&yL,dyda+i,7);
		*y += yL;
	}
	dyda[0] = 1.0;
	dyda[1] = x[0];
}

void DerivPhasedSplitPseudoVoigt(double x, double *a, double *y, double *dyda, int na)
{
	double alpha = a[5];
	if (alpha > 1.0)  {alpha = 1.0; a[5] = 1.0;}
	if (alpha < 0.0)  {alpha = 0.0; a[5] = 0.0;}

	double Gy = 0.0;
	double Ly = 0.0;
	double* Gdyda  = new double[na];
	for (int i=0; i < na; i++) Gdyda[i] = 0.0;
	*y = 0.0;

	DerivPhasedSplitGaussian(x, a, &Gy, Gdyda, na);
	DerivPhasedSplitLorentz(x, a, &Ly, dyda, na);

	*y = alpha*Ly + (1.0 - alpha)*Gy;

	for (int i=0; i <= 5; i++)
		dyda[i] = alpha*dyda[i] + (1.0 - alpha)* Gdyda[i];

	dyda[5] = Ly - Gy;
	delete[] Gdyda;
	return;
}

void PhasedSplitPseudoVoigt(double x, double *a, double *y, double *dyda, int na)
{
	double alpha = a[5];
	if (alpha > 1.0)  {alpha = 1.0; a[5] = 1.0;}
	if (alpha < 0.0)  {alpha = 0.0; a[5] = 0.0;}
	double Gy = 0.0;
	double Ly = 0.0;
	double* Gdyda  = new double[na];
	for (int i=0; i < na; i++) Gdyda[i] = 0.0;
	*y = 0.0;

	PhasedSplitGaussian(x, a, &Gy, Gdyda, na);
	PhasedSplitLorentz(x, a, &Ly, dyda, na);

	*y = alpha*Ly + (1.0 - alpha)*Gy;

	for (int i=0; i <= 5; i++)
		dyda[i] = alpha*dyda[i] + (1.0 - alpha)* Gdyda[i];

	dyda[5] = Ly - Gy;
    delete[] Gdyda;
	return;
}

