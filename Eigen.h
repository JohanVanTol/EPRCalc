void jacobi(double **a, int n, double d[], double **v, int *nrot);
void tred2(double **a, int n, double d[], double e[]);
void tred2_values_only(double **a, int n, double d[], double e[]);
void tred2_c(double **a, int n, double d[], double e[]);
void tqli(double d[], double e[], int n, double **z);
void tqli_values_only(double d[], double e[], int n, double **z);
double pythag(double a, double b);
void eigen_sort(double *v, double **a, int n);
Vector eigenval(const Ctensor &H);
