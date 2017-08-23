//---------------------------------------------------------------------------
#ifndef GaussianH
#define GaussianH

double Gaussian(double x, double xc, double Ymax, double w);
double DispGaussian(double x, double xc, double Ymax, double w);
double PhasedGaussian(double x, double xc, double Ymax, double w, double angle);

void Gaussian(double x, double *a, double *y, double *dyda, int na);
void DispGaussian(double x, double *a, double *y, double *dyda, int na);
void PhasedGaussian(double x, double *a, double *y, double *dyda, int na, double angle, double *anglederiv);

void MultGaussian(double x, double *a, double *y, double *dyda, int na);
void MultGaussian(double *x, int nx, double *a, double *y, double *dyda, int na);

void MultPhasedGaussian(double x, double *a, double *y, double *dyda, int na);
void MultPhasedGaussian(double *x, int nx, double *a, double *y, double *dyda, int na);

void MultPhasedSplitGaussian(double x, double *a, double *y, double *dyda, int na);
void MultPhasedSplitGaussian(double *x, int nx, double *a, double *y, double *dyda, int na);

//---

double DerivGaussian(double x, double xc, double Ymax, double w);
double DerivDispGaussian(double x, double xc, double Ymax, double w);
double DerivPhasedGaussian(double x, double xc, double Ymax, double w, double angle);

void DerivGaussian(double x, double *a, double *y, double *dyda, int na);
void DerivDispGaussian(double x, double *a, double *y, double *dyda, int na);
void DerivPhasedGaussian(double x, double *a, double *y, double *dyda, int na, double angle, double *anglederiv);

void MultDerivGaussian(double x, double *a, double *y, double *dyda, int na);
void MultDerivGaussian(double *x, int nx, double *a, double *y, double *dyda, int na);

void MultDerivPhasedGaussian(double x, double *a, double *y, double *dyda, int na);
void MultDerivPhasedGaussian(double *x, int nx, double *a, double *y, double *dyda, int na);

void MultDerivPhasedSplitGaussian(double x, double *a, double *y, double *dyda, int na);
void MultDerivPhasedSplitGaussian(double *x, int nx, double *a, double *y, double *dyda, int na);
//---------------------------------------------------------------------------
#endif
