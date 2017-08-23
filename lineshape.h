//---------------------------------------------------------------------------

#ifndef LineshapeH
#define LineshapeH
//---------------------------------------------------------------------------
double PseudoVoigt(double alpha, double x, double xc, double A, double w);
double DispPseudoVoigt(double alpha, double x,double xc,double A,double w);
void MultDerivPhasedSplitPseudoVoigt(double x, double *a, double *y, double *dyda, int na);
void MultPhasedSplitPseudoVoigt(double x, double *a, double *y, double *dyda, int na);
void MultDerivPhasedSplitPseudoVoigt(double *x, int nx, double *a, double *y, double *dyda, int na);
void MultPhasedSplitPseudoVoigt(double *x, int nx, double *a, double *y, double *dyda, int na);
void PhasedSplitPseudoVoigt(double x, double *a, double *y,double *dyda,int na);
void DerivPhasedSplitPseudoVoigt(double x, double *a, double *y,double *dyda,int na);


#endif
