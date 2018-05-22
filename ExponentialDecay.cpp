//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include "ExponentialDecay.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

void ExponentialDecay(double x, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;

	dyda[0] = 1.0;
	dyda[1] = exp(-x/a[2]);
	dyda[2] = 1.0*(a[1]*x/(a[2]*a[2]))*exp(-x/a[2]);
	*(yL) = a[0]+a[1]*exp(-x/a[2]);
}

void ExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;

	dyda[0] = 1.0;
	dyda[1] = exp(-x[0]/a[2]);
	dyda[2] = 1.0*(a[1]*x[0]/(a[2]*a[2]))*exp(-x[0]/a[2]);
	*(yL) = a[0]+a[1]*exp(-x[0]/a[2]);
}
void BiExponentialDecay(double x, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;
   if (a[4]<=0) a[4] = 1.0e-10;

	dyda[0] = 1.0;
	dyda[1] = exp(-x/a[2]);
	dyda[2] = 1.0*(a[1]*x/(a[2]*a[2]))*exp(-x/a[2]);
	dyda[3] = exp(-x/a[4]);
	dyda[4] = 1.0*(a[3]*x/(a[4]*a[4]))*exp(-x/a[4]);

	*(yL) = a[0]+a[1]*exp(-x/a[2])+a[3]*exp(-x/a[4]);

}

void BiExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;
   if (a[4]<=0) a[4] = 1.0e-10;

	dyda[0] = 1.0;
	dyda[1] = exp(-x[0]/a[2]);
	dyda[2] = 1.0*(a[1]*x[0]/(a[2]*a[2]))*exp(-x[0]/a[2]);
	dyda[3] = exp(-x[0]/a[4]);
	dyda[4] = 1.0*(a[3]*x[0]/(a[4]*a[4]))*exp(-x[0]/a[4]);

	*(yL) = a[0]+a[1]*exp(-x[0]/a[2])+a[3]*exp(-x[0]/a[4]);

}

void GaussExponentialDecay(double x, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;
   if (a[5]<=0) a[5] = 1.0e-10;

	dyda[0] = 1.0;
	dyda[1] = exp(-x/a[2]);
	dyda[2] = 1.0*(a[1]*x/(a[2]*a[2]))*exp(-x/a[2]);

	double Expon2 = (x-a[4])*(x-a[4])/(a[5]*a[5]) ;
	dyda[3] = exp(-Expon2);

	dyda[4] = a[3] * dyda[3] * 2.0 * (x-a[4])/(a[5]*a[5]);

	dyda[5] = a[3] * dyda[3] * 2.0 * Expon2/a[5];



	*(yL) = a[0]+a[1]*exp(-x/a[2]) + a[3]*exp(-Expon2);

}

void GaussExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;
   if (a[5]<=0) a[5] = 1.0e-10;

	dyda[0] = 1.0;
	dyda[1] = exp(-x[0]/a[2]);
	dyda[2] = 1.0*(a[1]*x[0]/(a[2]*a[2]))*exp(-x[0]/a[2]);

	double Expon2 = (x[0]-a[4])*(x[0]-a[4])/(a[5]*a[5]) ;
	dyda[3] = exp(-Expon2);

	dyda[4] = a[3]*dyda[3] * 2.0 * (x[0]-a[4])/(a[5]*a[5]);

	dyda[5] = a[3]*dyda[3] * 2.0 * Expon2/a[5];



	*(yL) = a[0]+a[1]*exp(-x[0]/a[2]) + a[3]*exp(-Expon2);

}

void StretchedExponentialDecay(double x, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;
   if (a[3]<=0.1) a[3] = 0.1;
   if (a[3]>3) a[3] = 3;

	dyda[0] = 1.0;
	dyda[1] = exp(-(pow((x/a[2]),a[3])));
	if (x<=0.0)  {
		dyda[2] = 0.0;
		dyda[3] = 0.0;
	  }
	  else {
		dyda[2] = a[1]*dyda[1]*a[3]*pow((x/a[2]),a[3]-1)*x/(a[2]*a[2]);
		dyda[3] = -a[1]*dyda[1]*pow((x/a[2]),a[3])*log(x/a[2]);
      }

	*(yL) = a[0]+a[1]*exp(-pow((x/a[2]), a[3]));

}

void StretchedExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na)
{
   if (a[2]<=0) a[2] = 1.0e-10;
   if (a[3]<=0.1) a[3] = 0.1;
   if (a[3]>3) a[3] = 3;

	dyda[0] = 1.0;
	dyda[1] = exp(-(pow((x[0]/a[2]),a[3])));
	if (x[0]<=0.0)
	{
	   dyda[2] = 0.0;
	   dyda[3] = 0.0;
	}
	  else
	  {
		dyda[2] = a[1]*dyda[1]*a[3]*pow((x[0]/a[2]),a[3]-1)*x[0]/(a[2]*a[2]);
		dyda[3] = -a[1]*dyda[1]*pow((x[0]/a[2]),a[3])*log(x[0]/a[2]);
      }

	*(yL) = a[0]+a[1]*exp(-pow((x[0]/a[2]), a[3]));

}
