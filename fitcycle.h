#ifndef FITCYCLE_H
#define FITCYCLE_H
void fitcycle(double **x, double *y, double *sy, int ndata, int nx,
					double *a, int *ia, int ma, double **covar, double **alpha,
					double *chisq, void (*funcs)(double *x, int nx, double *a, double *y,
									double *dyda, int ma), double *alambda);
#endif