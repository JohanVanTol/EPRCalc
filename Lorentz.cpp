#include "MY_CONST.H"
#include "Lorentz.h"
///////////////////////////////////////////////////////////
//
//  Lorentz simply returns the function value. To make the calculations
//  simpler and faster, not the intensity A(integrated value), but the
//  maximum value is used. The intensity A is related by A = Ymax*w*PI/2
//
double Lorentz(double x, double xc, double Ymax, double w)
{
        if (w<=0.0) return 0;
        double xl = 2.0*(x-xc)/w;
        return Ymax/(1.0 + xl*xl);
}

double NormLorentz(double xf, double xc, double A, double w)
{
// This returns a the value at position x of a Normalized Gaussian
// centered at xc with Full Width at Half Maximum (FWHM) w and intensity A.
// The normalization makes sure that the integrated intensity of the gaussian
// will be unity if A=1
	double x = 2*(xf-xc)/w;
	return (A/w)*(2.0/PI)*1.0/(1.0+x*x);
}

double NormDispLorentz(double xf, double xc, double A, double w)
{
// This returns a the value at position x of a Normalized Gaussian
// centered at xc with Full Width at Half Maximum (FWHM) w and intensity A.
// The normalization makes sure that the integrated intensity of the gaussian
// will be unity if A=1
	double x = 2*(xf-xc)/w;
	return (A/w)*(2.0/PI)*x/(1.0+x*x);
}
////////////////////////////////////////////////////////////
//
//    The dispersion of a lorentz
//
double DispLorentz(double x, double xc, double Ymax, double w)
{
        if (w<=0.0) return 0;
        double xl = 2.0*(x-xc)/w;
        return Ymax * xl /(1.0 + xl*xl);
}

//////////////////////////////////////////////////////////////
//
//     Mixture of dispersive and real parts
//
double PhasedLorentz(double x, double xc, double Ymax, double w, double angle)
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
//  This implementation calculates the value of the Lorentz line at a
//  point x, and the derivatives of the parameters.
//
//  The array a contains the parameters a[0] = Ymax = 2*A/(PI*w)
//                                      a[1] = xc  center;
//                                      a[2] = w   the FWHH
//


void Lorentz(double x, double *a, double *yL, double *dyda, int na)
{
    if (a[2]<=0) a[2] = 1.0e-10;

    double xl = 2*(x-a[1])/a[2];
    dyda[0] = 1.0/(1.0 + xl*xl);
    dyda[1] = dyda[0]*dyda[0] * a[0] * 4.0 * xl/a[2];
    dyda[2] = dyda[0]*dyda[0] * a[0] * 2.0 * xl*xl/a[2];
    *(yL) = a[0]*dyda[0];
}

void DispLorentz(double x, double *a, double *yL, double *dyda, int na)
{
    if (a[2]<=0) a[2] = 1.0e-10;

    double xl = 2*(x-a[1])/a[2];
    dyda[0] = xl /(1.0 + xl*xl);
    dyda[1] = a[0] * dyda[0]*dyda[0]/(xl*xl) * ( xl*xl - 1.0) * 2.0/a[2];
    dyda[2] = a[0] * dyda[0]*dyda[0]/(xl*xl) * ( xl*xl - 1.0) * xl/a[2];
    *(yL) = a[0]*dyda[0];
}

void PhasedLorentz(double x, double *a, double *yL, double *dyda, int na, double angle, double *dangle)
{
    if (a[2]<=0) a[2] = 1.0e-10;

    double xl = 2*(x-a[1])/a[2];
    angle *= PI/180.0;

    dyda[0] = (cos(angle) + xl * sin(angle)) /(1.0 + xl*xl);
    dyda[1] = a[0] * (cos(angle)*2.0*xl + sin(angle)*( xl*xl - 1.0))* 2.0/a[2];
    dyda[1] /= ((1.0 + xl*xl)*(1.0 + xl*xl));
    dyda[2] = a[0] * (cos(angle)*2.0*xl + sin(angle)*( xl*xl - 1.0)) * xl/a[2];
    dyda[2] /= ((1.0 + xl*xl)*(1.0 + xl*xl));

    *(yL) = (cos(angle) + xl*sin(angle)) * a[0]/(1.0 + xl*xl);

    *dangle = (-1.0*sin(angle) + xl*cos(angle)) * a[0]/(1.0 + xl*xl);
}

void PhasedLorentz(double x, double *a, double *yL, double *dyda, int na)
{
//    na should be 4 ...
//  a[0] = Amplitude
//  a[1] = Center field
//  a[2] = Width
//  a[3] = Phase in radians

	if (a[2]<=0) a[2] = 1.0e-10;  // width should be positive

	double xl = 2*(x-a[1])/a[2];

	dyda[0] = (cos(a[3]) + xl * sin(a[3])) /(1.0 + xl*xl);
	dyda[1] = a[0] * (cos(a[3])*2.0*xl + sin(a[3])*( xl*xl - 1.0))* 2.0/a[2];
	dyda[1] /= ((1.0 + xl*xl)*(1.0 + xl*xl));
	dyda[2] = a[0] * (cos(a[3])*2.0*xl + sin(a[3])*( xl*xl - 1.0)) * xl/a[2];
	dyda[2] /= ((1.0 + xl*xl)*(1.0 + xl*xl));

	*(yL) = (cos(a[3]) + xl*sin(a[3])) * a[0]/(1.0 + xl*xl);

	dyda[3] = (-1.0*sin(a[3]) + xl*cos(a[3])) * a[0]/(1.0 + xl*xl);
}

void PhasedSplitLorentz(double x, double *a, double *yL, double *dyda, int na)
{
//    na should be 4 ...
//  a[0] = Amplitude
//  a[1] = Center field
//  a[2] = Width
//  a[3] = Phase in radians
//  a[4] = Splitting between the multiple lines
//  a[5] = Multiplicity (number of lines)

	int nlines = int(a[6]);
	double *cl = new double[nlines];
	double xl;
	if (a[2]<=0) a[2] = 1.0e-10;  // width should be positive
	*(yL) = 0.0;
	for (int i=0;i<5;i++) dyda[i]= 0.0;


	for (int i=0; i < nlines; i++)
	{
		cl[i] = a[1] + 0.5*a[4]*(2*i-(nlines-1));

		xl = 2*(x-cl[i])/a[2];

		dyda[0] += (cos(a[3]) + xl * sin(a[3])) /(1.0 + xl*xl);
		dyda[1] += (a[0] * (cos(a[3])*2.0*xl + sin(a[3])*( xl*xl - 1.0))* 2.0/a[2])
					/ ((1.0 + xl*xl)*(1.0 + xl*xl));
		dyda[2] += (a[0] * (cos(a[3])*2.0*xl + sin(a[3])*( xl*xl - 1.0)) * xl/a[2])
					/ ((1.0 + xl*xl)*(1.0 + xl*xl));

		*(yL) += (cos(a[3]) + xl*sin(a[3])) * a[0]/(1.0 + xl*xl);

		dyda[3] += (-1.0*sin(a[3]) + xl*cos(a[3])) * a[0]/(1.0 + xl*xl);
		dyda[4] += 0.5*(2*i-(nlines-1))
				* (a[0] * (cos(a[3])*2.0*xl + sin(a[3])*( xl*xl - 1.0))* 2.0/a[2])
					/ ((1.0 + xl*xl)*(1.0 + xl*xl));

	}
}
////////////////////////////////////////////////////////////////
//
//  MultLorentz calculates a linear component plus a series of Lorentz
//  lineshapes. The na parameter contains the clue to the number of
//  lorentz lines
//

void MultLorentz(double x, double *a, double *y, double *dyda, int na)
{
	*y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-2; i+=3)
    {
        Lorentz(x,a+i,&yL,dyda+i,3);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultPhasedLorentz(double x, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
        PhasedLorentz(x,a+i,&yL,dyda+i,4);
        *y += yL;
    }
	dyda[0] = 1.0;
    dyda[1] = x;
}

void MultPhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=7)
    {
		PhasedSplitLorentz(x,a+i,&yL,dyda+i,7);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultPhasedLorentz(double *x, int nx, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
        PhasedLorentz(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultPhasedSplitLorentz(double *x, int nx, double *a, double *y, double *dyda, int na)
{
//  The phase angle is the last parameter a[na-1]
	*y = a[0] + a[1]*x[0];
    double yL = 0.0;
    double anglederiv = 0.0;
    dyda[na-1] = 0.0;
    for (int i=2; i<na-3; i+=7)
    {
		PhasedSplitLorentz(x[0],a+i,&yL,dyda+i,7);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultLorentz(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double *yL;
    yL = new double;
    for (int i=2; i<na-2; i+=3)
    {
        Lorentz(x[0],a+i,yL,dyda+i,3);
        *y += *yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
    delete yL;
}

///////////////////////////////////////////////////////////////////////////
//
//  Derivative Lorentz
//
double DerivLorentz(double x, double xc, double Ymax, double w)
{
        if (w<=0.0) return 0;
		double xl = 2.0*(x-xc)/w;
        double noemer = 1+xl*xl;
         return Ymax*-4.0*xl/(w*noemer*noemer);
}

double DerivDispersLorentz(double x, double xc, double Ymax, double w)
{
        if (w<=0.0) return 0;
        double xl = 2.0*(x-xc)/w;
        double noemer = 1 + xl*xl;
         return Ymax*2.0*(1-xl*xl)/(w*noemer*noemer);
}

/////////////////////////////////////////////////////////////////////
//
//  This implementation calculates the derivative value of the Lorentz line at a
//  point x, and the derivatives of the parameters.
//
//  The array a contains the parameters a[0] = Ymax = 2*A/(PI*w)
//                                      a[1] = xc  center;
//                                      a[2] = w   the FWHH
//


void DerivLorentz(double x, double *a, double *y, double *dyda, int na)
{
    if (a[2]<=0) a[2] = 1.0e-10;
    double xl = 2*(x-a[1])/a[2];
    double nn = 1.0/(1.0 + xl*xl);

    dyda[0] = -4.0*xl*nn*nn/a[2];

    dyda[1] = a[0]*8.0*nn*nn*nn*(1-3*xl*xl)/(a[2]*a[2]);

    dyda[2] = a[0]*8.0*nn*nn*nn*xl*(1-xl*xl)/(a[2]*a[2]);

    *y = a[0] * dyda[0];
}

void DerivDispersLorentz(double x, double *a, double *y, double *dyda, int na)
{
    if (a[2]<=0) a[2] = 1.0e-10;
    double xl = 2*(x-a[1])/a[2];

    double nn = 1.0/(1.0 + xl*xl);

    dyda[0] = 2.0*(1-xl*xl)*nn*nn/a[2];

    dyda[1] = a[0]*8.0*nn*nn*nn*xl*(3.0-xl*xl)/(a[2]*a[2]);

    dyda[2] = a[0]*-2.0*nn*nn*nn*(xl*xl*xl*xl-6*xl*xl+1)/(a[2]*a[2]);

    *y = a[0] * dyda[0];
}

void DerivPhasedLorentz(double x, double *a, double *y, double *dyda, int na)
{
	if (a[2]<=0) a[2] = 1.0e-10;
	double xl = 2*(x-a[1])/a[2];
	double cc = cos(a[3]);
	double ss = sin(a[3]);
	double nn = 1.0/(1.0 + xl*xl);

	dyda[0] = (-4.0*xl*cc + 2.0*(1.0-xl*xl)*ss)*nn*nn/a[2];

	dyda[1] = a[0]*8.0*(cc*(1.0-3.0*xl*xl)+ss*xl*(3.0-xl*xl))*nn*nn*nn/(a[2]*a[2]);

	dyda[2] = a[0]*2.0*(4.0*xl*(1.0-xl*xl)*cc - ss*(xl*xl*xl*xl-6*xl*xl+1.0))*nn*nn*nn/(a[2]*a[2]);

	dyda[3] = a[0]*(4.0*xl*ss + 2.0*cc*(1.0-xl*xl))*nn*nn/a[2];

	*y = a[0] * dyda[0];
}

void DerivPhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na)
{
	int nlines = int(a[6]);
	double *cl = new double[nlines];

	if (a[2]<=0) a[2] = 1.0e-10;  // width should be positive
	*y = 0.0;
	for (int i=0;i<5;i++) dyda[i]= 0.0;
	double cc = cos(a[3]);
	double ss = sin(a[3]);
	double xl = 0.0;
	double nn = 0.0;

	for (int i=0; i < nlines; i++)
	{
		cl[i] = a[1] + 0.5*a[4]*(2*i-(nlines-1));

		xl = 2*(x-cl[i])/a[2];
		nn = 1.0/(1.0 + xl*xl);

		dyda[0] += (-4.0*xl*cc + 2.0*(1.0-xl*xl)*ss)*nn*nn/a[2];

		dyda[1] += a[0]*8.0*(cc*(1.0-3.0*xl*xl)+ss*xl*(3.0-xl*xl))*nn*nn*nn/(a[2]*a[2]);

		dyda[2] += a[0]*2.0*(4.0*xl*(1.0-xl*xl)*cc - ss*(xl*xl*xl*xl-6*xl*xl+1.0))*nn*nn*nn/(a[2]*a[2]);

		dyda[3] += a[0]*(4.0*xl*ss + 2.0*cc*(1.0-xl*xl))*nn*nn/a[2];

		dyda[4] += 0.5*(2*i-(nlines-1))* a[0]*8.0*(cc*(1.0-3.0*xl*xl)+ss*xl*(3.0-xl*xl))*nn*nn*nn/(a[2]*a[2]);
		*y += a[0] * (-4.0*xl*cc + 2.0*(1.0-xl*xl)*ss)*nn*nn/a[2];
    }
}





void MultDerivLorentz(double x, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-2; i+=3)
    {
        DerivLorentz(x,a+i,&yL,dyda+i,3);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultDerivPhasedLorentz(double x, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
        DerivPhasedLorentz(x,a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultDerivPhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x;
    double yL = 0.0;
    for (int i=2; i<na-3; i+=7)
    {
		DerivPhasedSplitLorentz(x,a+i,&yL,dyda+i,7);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x;
}

void MultDerivLorentz(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    for (int i=2; i<na-2; i+=3)
    {
        DerivPhasedLorentz(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultDerivPhasedLorentz(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    for (int i=2; i<na-3; i+=4)
    {
		DerivPhasedLorentz(x[0],a+i,&yL,dyda+i,4);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];
}

void MultDerivPhasedSplitLorentz(double *x, int nx, double *a, double *y, double *dyda, int na)
{
    *y = a[0] + a[1]*x[0];
    double yL = 0.0;
    for (int i=2; i<na-6; i+=7)
    {
		DerivPhasedSplitLorentz(x[0],a+i,&yL,dyda+i,7);
        *y += yL;
    }
    dyda[0] = 1.0;
    dyda[1] = x[0];

}
