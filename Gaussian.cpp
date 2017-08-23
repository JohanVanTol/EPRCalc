#include "MY_CONST.H"
#include "Gaussian.h"

extern double GaussRange;
extern double* GaussDisp;
extern double* GaussDerivDisp;
extern int GaussNpts;

double qromb(double (*func)(double, double), double a, double b);
///////////////////////////////////////////////////////////
//
//  Gaussian simply returns the function value. To make the calculations
//  simpler and faster, not the intensity A(integrated value), but the
//  maximum value is used. The intensity A is related by A = Ymax*w*PI/2
//
double Gaussian(double x, double xc, double Ymax, double w)
{
		if (w<=0.0) return 0;
		double xl = 2.0*(x-xc)/w;
		return Ymax * exp(-xl*xl);

		// Note that in this case if x-xc = w/2 the exponent is
		// unity: the width is defined as the full width at 1/e
		// of the intensity
}

double QExp(double x, double xx)
{
		return exp(x*x - xx*xx);
}
////////////////////////////////////////////////////////////
//
//    The dispersion of a Gaussian
//
double DispGaussian(double x, double xc, double Ymax, double w)
{
		if (w<=0.0) return 0;
		double xl = 2.0*(x-xc)/w;
		double dis = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl);

		return Ymax*dis;
}

//////////////////////////////////////////////////////////////
//
//     Mixture of dispersive and real parts
//
double PhasedGaussian(double x, double xc, double Ymax, double w, double angle)
{
//  this includes a mixture of absorption and dispersion through 'angle'
//  angle is aasumed to be in degrees
		if (w<=0.0) return 0;
		angle *= PI/180.0;
		double xl = 2.0*(x-xc)/w;
		return (cos(angle) + xl*sin(angle)) * Ymax/(1.0 + xl*xl);

}

/////////////////////////////////////////////////////////////////////
//
//  This implementation calculates the value of the Gaussian line at a
//  point x, and the derivatives of the parameters.
//
//  The array a contains the parameters a[0] = Ymax = 2*A/(PI*w)
//                                      a[1] = xc  center;
//                                      a[2] = w   the FWHH
//


void Gaussian(double x, double *a, double *yL, double *dyda, int na)
{
	if (a[2]<=0) a[2] = 1.0e-10;

	double xl = 2.0*(x-a[1])/a[2];
	dyda[0] = exp(-xl*xl);
	dyda[1] = dyda[0] * a[0] * 4.0 * xl/a[2];
	dyda[2] = dyda[0] * a[0] * 2.0 * xl*xl/a[2];
	*(yL) = a[0]*dyda[0];
}

void DispGaussian(double x, double *a, double *yL, double *dyda, int na)
{
	if (a[2]<=0) a[2] = 1.0e-10;
	double xl = 2*(x-a[1])/a[2];
	int j;
	double pos, fraction, change;
	if (GaussNpts > 0)  // If we have a Gaussian Line Shape file loaded
	{
		pos = (xl+GaussRange/2.0)*double(GaussNpts-1)/GaussRange;
		if (pos <= 0.0)
		{
			dyda[0] = GaussDisp[0]*GaussRange/(-2.0*xl);         // not a great approximation...
			dyda[1] = 0.0;
			dyda[2] = 0.0;
			*yL = a[0] * dyda[0];
		}
		  else if (pos >= GaussNpts) {
			dyda[0] = GaussDisp[GaussNpts-1]*GaussRange/(2.0*xl);
			dyda[1] = 0.0;
			dyda[2] = 0.0;
			*yL = a[0] * dyda[0];
			}
			  else
			  {
				j = (int)floor(pos);
				fraction = pos - j;
				change = GaussDisp[j+1]-GaussDisp[j];
				dyda[0] = GaussDisp[j] + fraction * change;
				change /= GaussRange/(GaussNpts-1);
				dyda[1] = -2.0*change*a[0]/a[2];
				dyda[2] = -xl*change*a[0]/a[2];
				*(yL) = a[0]*dyda[0];
			  }
		}
		else{

		dyda[0] = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl);
		double deltaxl = 0.01;
		double plus = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl+deltaxl);
		double min = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl-deltaxl);
		change = (plus-min)/2*deltaxl;
		dyda[1] = -2.0 * a[0]* change/a[2];
		dyda[2] = -xl * a[0] * change/a[2];
		*(yL) = a[0]*dyda[0];
		}
		return;
}

void PhasedGaussian(double x, double *a, double *yL, double *dyda, int na, double angle, double *dangle)
{
	if (a[2]<=0) a[2] = 1.0e-10;

	double xl = 2*(x-a[1])/a[2];
	angle *= PI/180.0;
	double* dispdyda = new double[na];

	double dispval;

	Gaussian(x,a,yL,dyda,na);
	DispGaussian(x,a,&dispval,dispdyda,na);

	*dangle = -1.0*sin(angle)* (*yL) + cos(angle)*dispval;

	*yL *= cos(angle);
	*yL += (sin(angle) * dispval);

	dyda[0] *= cos(angle);
	dyda[0] += sin(angle) * dispdyda[0];
	dyda[1] *= cos(angle);
	dyda[1] += sin(angle) * dispdyda[1];
	dyda[2] *= cos(angle);
	dyda[2] += sin(angle) * dispdyda[2];
	delete [] dispdyda;

}

void PhasedGaussian(double x, double *a, double *yL, double *dyda, int na)
{
//    na should be 4 ...
//  a[0] = Amplitude
//  a[1] = Center field
//  a[2] = Width
//  a[3] = Phase in radians

	if (a[2]<=0) a[2] = 1.0e-10;  // width should be positive

	double* dispdyda  = new double[na];
	double dispval;

	Gaussian(x,a,yL,dyda,na);
	DispGaussian(x,a,&dispval,dispdyda,na);

	double angle = a[3];
	dyda[3] = -1.0*sin(angle)* (*yL) + cos(angle)*dispval;

	*yL *= cos(angle);
	*yL += (sin(angle) * dispval);

	dyda[0] *= cos(angle);
	dyda[0] += sin(angle) * dispdyda[0];
	dyda[1] *= cos(angle);
	dyda[1] += sin(angle) * dispdyda[1];
	dyda[2] *= cos(angle);
	dyda[2] += sin(angle) * dispdyda[2];

	delete[] dispdyda;

}

////////////////////////////////////////////////////////////////
//
//  MultLorentz calculates a linear component plus a series of Lorentz
//  lineshapes. The na parameter contains the clue to the number of
//  lorentz lines
//

void MultGaussian(double x, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-2; i+=3)
    {
        Gaussian(x,a+i,&yL,dyda+i,3);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultPhasedGaussian(double x, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=4)
	{
        PhasedGaussian(x,a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
	dyda[1] = x;
}

void MultPhasedSplitGaussian(double x, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=7)
	{
        PhasedGaussian(x,a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultPhasedGaussian(double *x, int nx, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
        PhasedGaussian(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
	dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultPhasedSplitGaussian(double *x, int nx, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=7)
    {
        PhasedGaussian(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
	dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultGaussian(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double *yL;
    yL = new double;
    for (int i=2; i<na-2; i+=3)
    {
        Gaussian(x[0],a+i,yL,dyda+i,3);
        *y += *yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
    delete yL;
}

///////////////////////////////////////////////////////////////////////////
//
//  Derivative Gaussian
//
double DerivGaussian(double x, double xc, double Ymax, double w)
{
		if (w<=0.0) return 0;
        double xl = 2.0*(x-xc)/w;
        double noemer = 1+xl*xl;
         return Ymax*-4.0*xl/(w*noemer*noemer);
}

double DerivDispersGaussian(double x, double xc, double Ymax, double w)
{
        if (w<=0.0) return 0;
        double xl = 2.0*(x-xc)/w;
        double noemer = 1 + xl*xl;
		 return Ymax*2.0*(1-xl*xl)/(w*noemer*noemer);
}

/////////////////////////////////////////////////////////////////////
//
//  This implementation calculates the derivative value of the Gaussian line at a
//  point x, and the derivatives of the parameters.
//
//  The array a contains the parameters a[0] = Ymax = 2*A/(PI*w)
//                                      a[1] = xc  center;
//                                      a[2] = w   the FWHH
//


void DerivGaussian(double x, double *a, double *y, double *dyda, int na)
{
    if (a[2]<=0) a[2] = 1.0e-10;
    double xl = 2*(x-a[1])/a[2];

	dyda[0] = -4.0 * xl*exp(-xl*xl)/a[2];

	dyda[1] = a[0]* 8.0 * (1.0 - 2*xl*xl)*exp(-xl*xl)/(a[2]*a[2]);

	dyda[2] = a[0]* 4.0 * xl * (1.0 - 2*xl*xl)*exp(-xl*xl)/(a[2]*a[2]);

	*y = a[0] * dyda[0];
}


void DerivDispersGaussian(double x, double *a, double *yG, double *dyda, int na)
{
	if (a[2]<=0) a[2] = 1.0e-10;
	double xl = 2*(x-a[1])/a[2];
	int j;
	double pos, fraction, change;
	if (GaussNpts > 0)  // If we have a Gaussian Line Shape file loaded
	{
		pos = (xl+GaussRange/2.0)*double(GaussNpts-1)/GaussRange;
		if (pos <= 0.0)       // left of the range
		{
			dyda[0] = GaussDerivDisp[0]*GaussRange/(-2.0*xl);         // not a great approximation...
			dyda[1] = 0.0;
			dyda[2] = 0.0;
			*yG = a[0] * dyda[0];
		}
		  else if (pos >= GaussNpts) {       //right of range
			dyda[0] = GaussDerivDisp[GaussNpts-1]*GaussRange/(2.0*xl);
			dyda[1] = 0.0;
			dyda[2] = 0.0;
			*yG = a[0] * dyda[0];
			}
			  else
			  {
				j = (int)floor(pos);
				fraction = pos - j;
				change = GaussDerivDisp[j+1]-GaussDerivDisp[j];
				dyda[0] = 2*(GaussDerivDisp[j] + fraction * change)/a[2];
				change /= GaussRange/(GaussNpts-1);
				dyda[1] = -4.0*change*a[0]/(a[2]*a[2]);
				dyda[2] = -2.0*(1+xl)*change*a[0]/(a[2]*a[2]);
				*(yG) = a[0]*dyda[0];
			  }
		}
		else{

		dyda[0] = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl);
		double deltaxl = 0.01;
		double plus = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl+deltaxl);
		double min = -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl-deltaxl);
		change = (plus-min)/2*deltaxl;
		dyda[1] = -2.0 * a[0]* change/a[2];
		dyda[2] = -xl * a[0] * change/a[2];
		*(yG) = a[0]*change;
		}
		return;
}

void DerivPhasedGaussian(double x, double *a, double *yG, double *dyda, int na)
{
	if (a[2]<=0) a[2] = 1.0e-10;

	double xl = 2*(x-a[1])/a[2];
	double angle = a[3];
	double* dispdyda = new double[na];

	double dispval;

	DerivGaussian(x,a,yG,dyda,na);
	DerivDispersGaussian(x,a,&dispval,dispdyda,na);

	dyda[3] = -1.0*sin(angle)* (*yG) + cos(angle)*dispval;

	*yG *= cos(angle);
	*yG += (sin(angle) * dispval);

	dyda[0] *= cos(angle);
	dyda[0] += sin(angle) * dispdyda[0];
	dyda[1] *= cos(angle);
	dyda[1] += sin(angle) * dispdyda[1];
	dyda[2] *= cos(angle);
	dyda[2] += sin(angle) * dispdyda[2];
	delete [] dispdyda;

}




void MultDerivGaussian(double x, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-2; i+=3)
	{
        DerivGaussian(x,a+i,&yL,dyda+i,3);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultDerivPhasedGaussian(double x, double *a, double *y, double *dyda, int na)
{
	*y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
        DerivPhasedGaussian(x,a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultDerivPhasedSplitGaussian(double x, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-3; i+=7)
    {
        DerivPhasedGaussian(x,a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultDerivGaussian(double *x, int nx, double *a, double *y, double *dyda, int na)
{
	*y = a[0] + a[1]*x[0];
    double yL = 0.0;
    for (int i=2; i<na-2; i+=3)
    {
        DerivPhasedGaussian(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultDerivPhasedGaussian(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
        DerivPhasedGaussian(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultDerivPhasedSplitGaussian(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    for (int i=2; i<na-3; i+=7)
    {
        DerivPhasedGaussian(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

