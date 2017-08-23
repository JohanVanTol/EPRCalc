#include "MY_CONST.H"
#include "Gaussian2.h"
#define SQRT_LN2_BY_PI 0.4697186393498257

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

double NormGauss(double xf, double xc, double A, double w)
{
// This returns a the value at position x of a Normalized Gaussian
// centered at xc with Full Width at Half Maximum (FWHM) w and intensity A.
// The normalization makes sure that the integrated intensity of the gaussian
// will be unity if A=1
	double x = 2.0*(xf-xc)/w;
	return 2.0*(A/w)* SQRT_LN2_BY_PI *exp(-log(2.0)*x*x);
}




double Gaussian(double x, double xc, double Ymax, double w)
{
		if (w<=0.0) return 0;
		double xl = 2.0*(x-xc)/w;
		return Ymax * exp(-xl*xl*log(2));

		// Note that in this case if x-xc = w/2 the exponent is
		// -ln(2): the width is defined as the full width at 1/2
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
		double xl = sqrt(log(2))*2.0*(x-xc)/w;
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
		return (cos(angle)*Gaussian(x,xc,Ymax,w)+ sin(angle)*DispGaussian(x,xc,Ymax,w));

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
	dyda[0] = exp(-xl*xl*log(2));      // note the ln(2)
									   // makes a[2] defined as FWHM
	dyda[1] = dyda[0] * a[0] * 4.0 * log(2) * xl/a[2];
	dyda[2] = dyda[0] * a[0] * 2.0 * log(2) * xl*xl/a[2];
	*(yL) = a[0]*dyda[0];
}

double InterPolDispGaussian(double x, double xc, double Ymax, double w)
{
	if (w<=0) w = 1.0e-10;
	double xl = sqrt(log(2.0))*2.0*(x-xc)/w;
	int j;
	double pos, fraction, change;
	if (GaussNpts > 0)  // If we have a Gaussian Line Shape file loaded
	{
		pos = (xl+GaussRange/2.0)*double(GaussNpts-1)/GaussRange;
//		Note that GaussRange is expressed in xl = 2*(x-xo)/w
		if (pos <= 0.0)
		{
			return Ymax*GaussDisp[0]*GaussRange/(-2.0*xl);         // not a great approximation...
		}
		  else if (pos >= GaussNpts) {
			return Ymax * GaussDisp[GaussNpts-1]*GaussRange/(2.0*xl);
			}
			  else
			  {
				j = (int)floor(pos);
				fraction = pos - j;
				change = GaussDisp[j+1]-GaussDisp[j];
				return  Ymax*(GaussDisp[j] + fraction * change);
			  }
		}
		else{

		return Ymax * -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl);
		}
		return 0.0;
}

double NormDispGaussian(double x, double xc, double A, double w)
{
	if (w<=0) w = 1.0e-10;
	double xl = sqrt(log(2.0))*2.0*(x-xc)/w;
	int j;
	double pos, fraction, change;
	if (GaussNpts > 0)  // If we have a Gaussian Line Shape file loaded
	{
		pos = (xl+GaussRange/2.0)*double(GaussNpts-1)/GaussRange;
//		Note that GaussRange is expressed in xl = 2*(x-xo)/w
		if (pos <= 0.0)
		{
			return 2.0*(A/w)* SQRT_LN2_BY_PI*GaussDisp[0]*GaussRange/(-2.0*xl);         // not a great approximation...
		}
		  else if (pos >= GaussNpts) {
			return 2.0*(A/w)* SQRT_LN2_BY_PI * GaussDisp[GaussNpts-1]*GaussRange/(2.0*xl);
			}
			  else
			  {
				j = (int)floor(pos);
				fraction = pos - j;
				change = GaussDisp[j+1]-GaussDisp[j];
				return  2.0*(A/w)* SQRT_LN2_BY_PI*(GaussDisp[j] + fraction * change);
			  }
		}
		else{

		return 2.0*(A/w)* SQRT_LN2_BY_PI * -sqrt(2.0*3.14159265) * qromb(QExp, 0, xl);
		}
		return 0.0;
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
//		Note that GaussRange is expressed in xl = 2*(x-xo)/w
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

	Gaussian(x,a,yL,dyda,3);
	DispGaussian(x,a,&dispval,dispdyda,3);

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

void PhasedSplitGaussian(double x, double *a, double *yL, double *dyda, int na)
{
//    na should be 7 ...
//  a[0] = Amplitude
//  a[1] = Center field
//  a[2] = Width
//  a[3] = Phase in radians
//  a[4] = Splitting between the multiple lines
//  a[5] = Multiplicity (number of lines)

	int nlines = int(a[6]);
/*  NOTE  NOTE This routine is unfinished
*/
	double *cl = new double[nlines];
	double xl, dispval;
	double tempY;
	tempY = 0;
	if (a[2]<=0) a[2] = 1.0e-10;  // width should be positive
	*(yL) = 0.0;
	for (int i=0;i<5;i++) dyda[i]= 0.0;

	double* tempdyda  = new double[na];
	double* dispdyda  = new double[na];
	double* tempa = new double[na];

	for (int i=0; i < na; i++) {
		tempdyda[i] = 0.0;
		dispdyda[i] = 0.0;
		dyda[i] = 0.0;
		tempa[i] = a[i];
	}
	*yL = 0.0;

	for (int i=0; i < nlines; i++)
	{
		tempa[1] = a[1] + 0.5*a[4]*(2*i-(nlines-1));

		Gaussian(x,	tempa, &tempY,tempdyda,4);
		DispGaussian(x, tempa, &dispval, dispdyda,na);

		dyda[0] += ( cos(a[3])* tempdyda[0] + sin(a[3]) * dispdyda[0]);
		dyda[1] += ( cos(a[3])* tempdyda[1] + sin(a[3]) * dispdyda[1]);
		dyda[2] += ( cos(a[3])* tempdyda[2] + sin(a[3]) * dispdyda[2]);

		*(yL) += tempY * cos(a[3]) + dispval* sin(a[3]);

		dyda[3] += (-1.0*sin(a[3]))* (tempY) + cos(a[3]) *dispval;
		dyda[4] += 0.5*(2*i-(nlines-1)) * (cos(a[3])*tempdyda[1] + sin(a[3])*dispdyda[1]);

	}
	delete [] tempa;
	delete [] tempdyda;
	delete [] dispdyda;
	delete [] cl;
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
		PhasedSplitGaussian(x,a+i,&yL,dyda+i,7);
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
        PhasedSplitGaussian(x[0],a+i,&yL,dyda+i,7);
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
 // THIS IS NOT A GAUSSIAN !!
		if (w<=0.0) return 0;
		double xl = 2.0*(x-xc)/w;

		double noemer = 1+xl*xl;
		 return Ymax*-4.0*xl/(w*noemer*noemer);
}

double DerivDispersGaussian(double x, double xc, double Ymax, double w)
{
// THIS IS NOT A DERUV GAUSSIAN
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

	dyda[0] = -4.0 * log(2) * xl*exp(-xl*xl*log(2))/a[2];

	dyda[1] = a[0]* 8.0 * log(2)*log(2)*(1.0 - 2*xl*xl)*exp(-xl*xl)/(a[2]*a[2]);

	dyda[2] = a[0]* 4.0 * xl * log(2)*log(2)*(1.0 - 2*xl*xl)*exp(-xl*xl)/(a[2]*a[2]);

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

void DerivPhasedSplitGaussian(double x, double *a, double *yL, double *dyda, int na)
{
//    na should be 7 ...
//  a[0] = Amplitude
//  a[1] = Center field
//  a[2] = Width
//  a[3] = Phase in radians
//  a[4] = Splitting between the multiple lines
//  a[5] = Multiplicity (number of lines)

	int nlines = int(a[6]);
/*  NOTE  NOTE This routine is unfinished
*/
	double *cl = new double[nlines];
	double xl, dispval;
	double tempY;
	tempY = 0;
	if (a[2]<=0) a[2] = 1.0e-10;  // width should be positive
	*(yL) = 0.0;
	for (int i=0;i<5;i++) dyda[i]= 0.0;

	double* tempdyda  = new double[na];
	double* dispdyda  = new double[na];
	double* tempa = new double[na];

	for (int i=0; i < na; i++) {
		tempdyda[i] = 0.0;
		dispdyda[i] = 0.0;
		dyda[i] = 0.0;
		tempa[i] = a[i];
	}
	*yL = 0.0;

	for (int i=0; i < nlines; i++)
	{
		tempa[1] = a[1] + 0.5*a[4]*(2*i-(nlines-1));

		DerivGaussian(x,	tempa, &tempY,tempdyda,4);
		DerivDispersGaussian(x, tempa, &dispval, dispdyda,na);

		dyda[0] += ( cos(a[3])* tempdyda[0] + sin(a[3]) * dispdyda[0]);
		dyda[1] += ( cos(a[3])* tempdyda[1] + sin(a[3]) * dispdyda[1]);
		dyda[2] += ( cos(a[3])* tempdyda[2] + sin(a[3]) * dispdyda[2]);

		*(yL) += tempY * cos(a[3]) + dispval* sin(a[3]);

		dyda[3] += (-1.0*sin(a[3]))* (tempY) + cos(a[3]) *dispval;
		dyda[4] += 0.5*(2*i-(nlines-1)) * (cos(a[3])*tempdyda[1] + sin(a[3])*dispdyda[1]);

	}
	delete [] tempa;
	delete [] tempdyda;
	delete [] dispdyda;
	delete [] cl;
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
		DerivPhasedSplitGaussian(x,a+i,&yL,dyda+i,7);
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
	for (int i=2; i<na-6; i+=7)
	{
		DerivPhasedSplitGaussian(x[0],a+i,&yL,dyda+i,7);
		*y += yL;
	}
	dyda[0] = 1.0;
	dyda[1] = x[0];
}

