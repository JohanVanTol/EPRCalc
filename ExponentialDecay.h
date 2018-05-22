//---------------------------------------------------------------------------

#ifndef ExponentialDecayH
#define ExponentialDecayH
//---------------------------------------------------------------------------
#endif
void ExponentialDecay(double x, double *a, double *yL, double *dyda, int na);
void ExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na);

void BiExponentialDecay(double x, double *a, double *yL, double *dyda, int na);
void BiExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na);

void GaussExponentialDecay(double x, double *a, double *yL, double *dyda, int na);
void GaussExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na);

void StretchedExponentialDecay(double x, double *a, double *yL, double *dyda, int na);
void StretchedExponentialDecay(double *x, int nx, double *a, double *yL, double *dyda, int na);
