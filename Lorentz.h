//---------------------------------------------------------------------------
#ifndef LorentzH
#define LorentzH

double NormLorentz(double xf, double xc, double A, double w);
double NormDispLorentz(double xf, double xc, double A, double w);
double Lorentz(double x, double xc, double Ymax, double w);
double DispLorentz(double x, double xc, double Ymax, double w);
double PhasedLorentz(double x, double xc, double Ymax, double w, double angle);

void Lorentz(double x, double *a, double *y, double *dyda, int na);
void DispLorentz(double x, double *a, double *y, double *dyda, int na);
void PhasedLorentz(double x, double *a, double *y, double *dyda, int na, double angle, double *anglederiv);
void PhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na);

void MultLorentz(double x, double *a, double *y, double *dyda, int na);
void MultLorentz(double *x, int nx, double *a, double *y, double *dyda, int na);

void MultPhasedLorentz(double x, double *a, double *y, double *dyda, int na);
void MultPhasedLorentz(double *x, int nx, double *a, double *y, double *dyda, int na);
void MultPhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na);
void MultPhasedSplitLorentz(double *x, int nx, double *a, double *y, double *dyda, int na);

//---

double DerivLorentz(double x, double xc, double Ymax, double w);
double DerivDispLorentz(double x, double xc, double Ymax, double w);
double DerivPhasedLorentz(double x, double xc, double Ymax, double w, double angle);
void DerivPhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na);

void DerivLorentz(double x, double *a, double *y, double *dyda, int na);
void DerivDispLorentz(double x, double *a, double *y, double *dyda, int na);
void DerivPhasedLorentz(double x, double *a, double *y, double *dyda, int na, double angle, double *anglederiv);

void MultDerivLorentz(double x, double *a, double *y, double *dyda, int na);
void MultDerivLorentz(double *x, int nx, double *a, double *y, double *dyda, int na);

void MultDerivPhasedLorentz(double x, double *a, double *y, double *dyda, int na);
void MultDerivPhasedSplitLorentz(double x, double *a, double *y, double *dyda, int na);
void MultDerivPhasedLorentz(double *x, int nx, double *a, double *y, double *dyda, int na);
void MultDerivPhasedSplitLorentz(double *x, int nx, double *a, double *y, double *dyda, int na);
//---------------------------------------------------------------------------
#endif
