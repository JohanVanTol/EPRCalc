#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

void gaussj2(double **a, int n, double *b);
void covsrt(double **covar, int ma, int *ia, int mfit);
void mrqcof(double **x, double *y, double *sy, int ndata, int nx,
				double *a, int *ia, int ma, double **alpha, double *beta,
				double *chisq, void (*funcs)(double *x, int nx, double *a, double *y,
									double *dyda, int ma));


// Adapted from the numerical recipes
void fitcycle(double **x, double *y, double *sy, int ndata, int nx,
					double *a, int *ia, int ma, double **covar, double **alpha,
					double *chisq, void (*funcs)(double *x, int nx, double *a, double *y,
									double *dyda, int ma), double *alambda)
{
	int j,k,l,m;
	static int mfit;
	static double ochisq, *atry, *beta, *da, *oneda;

	if (*alambda < 0)   //If first call  , initialisation
	{
		atry = new double[ma];
		beta = new double[ma];
		da = new double[ma];
		for (mfit=0, j=0; j<ma; j++) if (ia[j]) mfit++;
//
		oneda = new double[mfit];
//
		*alambda=0.1;
		mrqcof(x,y,sy,ndata,nx,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=0; j<ma;j++) atry[j]=a[j];
	}

	for (j=-1, l=0; l<ma;l++)
	{
		if (ia[l])
		{
			for (j++,k=0, m=0;m<ma ;m++)
			if (ia[m])
			{
				covar[j][k] = alpha[j][k];
				k++;
			}
			covar[j][j] = alpha[j][j]*(1.0 + (*alambda));
			oneda[j] = beta[j];
		}
	}

	gaussj2(covar, mfit, oneda);             //Matrix solution

	for (j=0; j<mfit; j++) da[j] = oneda[j];
	if (*alambda == 0.0)                       //If converged
	{
		covsrt(covar,ma,ia,mfit);
		delete [] oneda;
		delete [] da;
		delete [] beta;
		delete [] atry;
		return;
	}

	for (j=0, l=0;l<ma;l++) if (ia[l]) atry[l] = a[l] + da[j++];

	mrqcof(x,y,sy,ndata,nx,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq)            //  if succesful
	{
		*alambda *= 0.1;
		ochisq = (*chisq);
		for (j=0,l=0; l<ma; l++)
		{
			if (ia[l])
			{
				for (k=0, m=0;m<ma ;m++)
				if (ia[m])
				{
					alpha[j][k] = covar[j][k];
					k++;
				}
				beta[j] = da[j];
				j++;
				a[l] = atry[l];
			}
		}
	}
	else
	 {
			*alambda *= 10;
			*chisq = ochisq;
	 }
}



