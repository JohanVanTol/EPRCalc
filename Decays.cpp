#include <math.h>
#include <algorithm>
#include <stdlib>

double** Decay;

void rk4(double *y, double *dydx, int n, double x, double h, double *yout,
	void (*derivs)(double , double*, double*));

void derivs(double x, double* y, double *dy);

int CalcDecayedPops(double kx, double ky, double kz, double SLRx, double SLRy,
            double SLRz, double tt, double freq, double kT,
            double Px, double Py,double Pz, double* pops, double l, double m, double n)
{
    int nn=3;

    double* y = new double[nn];
    double* dydx = new double[nn];
    double* newy = new double[nn];

    double kSL = l*l*SLRx + m*m*SLRy + n*n*SLRz;
    double Bfac = exp(-freq/kT);

//  Starting Populations
    y[0] = l*l*(Py+Pz)/2.0 + m*m*(Px+Pz)/2.0 + n*n*(Px+Py)/2.0;
    y[2] = y[0];
    y[1] = l*l*Px + m*m*Py + n*n*Pz;

    //  Excited state decay
    Decay = new double*[3];
        for (int i=0; i<3;i++)
            Decay[i] = new double[3];

    Decay[0][0] = -(l*l*(ky+kz) + m*m*(kx+kz) + n*n*(kx+ky))/2.0;
    Decay[1][1] = -(l*l*kx + m*m*ky + n*n*kz);
    Decay[2][2] = Decay[0][0];

    // Spin Lattice Relaxation
    Decay[0][0] -= kSL * Bfac;
    Decay[0][1] = kSL;
    Decay[0][2] = 0.0;
    Decay[1][0] = kSL*Bfac;
    Decay[1][1] -= kSL + kSL*Bfac;
    Decay[1][2] = kSL;
    Decay[2][0] = 0.0;
    Decay[2][1] = kSL * Bfac;
    Decay[2][2] -= kSL;

    double maxim = kSL;
    if (kx > maxim) maxim = kx;
    if (ky > maxim) maxim = ky;
    if (kz > maxim) maxim = kz;

    int nsteps = ceil(3*tt*maxim);
    if (nsteps == 0) nsteps = 1;
    double stepsize = tt/nsteps;
    double x=0.0;
    for (int i=0;i<nsteps;i++)
    {
        derivs(x,y,dydx);
        rk4(y,dydx,nn ,x,stepsize,newy,derivs);
        x += stepsize;
        for (int j=0;j<nn;j++)
            y[j]= newy[j];
    }

    for (int i=0; i<3; i++) pops[i] = newy[i];

    for (int i=0; i<3;i++)
        delete[] Decay[i];
    delete[] Decay;
    delete[] newy;
    delete[] dydx;
    delete[] y;

    return 1;
}

int CalcDecayedPops_RP(double kx, double ky, double kz, double SLRx, double SLRy,
            double SLRz, double tt, double freq, double kT,
            double Px, double Py,double Pz, double* pops, double l, double m, double n)
{
    int nn=3;

    double* y = new double[nn];
    double* dydx = new double[nn];
    double* newy = new double[nn];

    double kSL = l*l*SLRx + m*m*SLRy + n*n*SLRz;
    double Bfac = exp(-freq/kT);

//  Starting Populations
    y[0] = 0.0;
    y[2] = 0.0;
    y[1] = l*l*Px + m*m*Py + n*n*Pz;

    //  Excited state decay
    Decay = new double*[3];
        for (int i=0; i<3;i++)
            Decay[i] = new double[3];

    Decay[0][0] = -(l*l*(ky+kz) + m*m*(kx+kz) + n*n*(kx+ky))/2.0;
    Decay[1][1] = -(l*l*kx + m*m*ky + n*n*kz);
    Decay[2][2] = Decay[0][0];

    // Spin Lattice Relaxation
    Decay[0][0] -= kSL * Bfac;
    Decay[0][1] = kSL;
    Decay[0][2] = 0.0;
    Decay[1][0] = kSL*Bfac;
    Decay[1][1] -= kSL + kSL*Bfac;
    Decay[1][2] = kSL;
    Decay[2][0] = 0.0;
    Decay[2][1] = kSL * Bfac;
    Decay[2][2] -= kSL;

    double maxim = kSL;
    if (kx > maxim) maxim = kx;
    if (ky > maxim) maxim = ky;
    if (kz > maxim) maxim = kz;

    int nsteps = ceil(3*tt*maxim);
    if (nsteps == 0) nsteps = 1;
    double stepsize = tt/nsteps;
    double x=0.0;
    for (int i=0;i<nsteps;i++)
    {
        derivs(x,y,dydx);
        rk4(y,dydx,nn ,x,stepsize,newy,derivs);
        x += stepsize;
        for (int j=0;j<nn;j++)
            y[j]= newy[j];
    }

    for (int i=0; i<3; i++) pops[i] = newy[i];

    for (int i=0; i<3;i++)
        delete[] Decay[i];
    delete[] Decay;
    delete[] newy;
    delete[] dydx;
    delete[] y;

    return 1;
}


void derivs(double x, double* y, double *dy)
{
    for (int i=0; i<3; i++)
    {
        dy[i] = 0.0;
        for (int j=0;j<3;j++)
           dy[i] += Decay[i][j]*y[j];
    }
    return;
}

void rk4(double *y, double *dydx, int n, double x, double h, double *yout,
	void (*derivs)(double , double*, double*))
{
	double xh,hh,h6,*dym,*dyt,*yt;

	dym = new double[n];
	dyt = new double[n];
	yt = new double[n];
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (int i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (int i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (int i=0;i<n;i++)
    {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (int i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	delete[] yt;
	delete[] dyt;
	delete[] dym;
}

