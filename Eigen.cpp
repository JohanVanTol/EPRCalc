////////////////////////////////////////////////////////////////////////////
//		EIGEN.CPP
//			Last revised 9 sept 1995
//			Copyright Hans van Tol
//
// This is a series of routines that can be useful for
// eigenvalue problems. The main included programs
// are the Jacobi and Householder/QL routines from
// the numerical recipes packet.
// These however are useful only for real symmetric
// matrices of floats, if they are used as they are.
// At present complex routines are not given here,
// allthough they might be more efficient.
// The versions of NR have simply been adapted for
// doubles, and are rewritten for C++ !
// The correct declarations can be found in eigen.h
// Another change with respect to the original numerical
//	recipes routines concerns the indices, which run from
// 0 to n-1, instead of the NR convention 1..n
//
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "vector.h"
#include "tensor.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

//////////////////////////////////////////////////////////////////////////
//
//	nrerror
//
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

//////////////////////////////////////////////////////////////////////////
//
//	jacobi
//
//		Diagonalizes a real symmetric matrix
//
//	PARAMETERS
//		**a	matrix to be diagonalized     			input, NB NOT UNHARMED
//		n		order of matrix                        input
//		*d		diagonal values (eigenvalues)          output
//		**v	transformation matrix (eigenfunctions) output
//		*nrot number of iterations							output
//

void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b= new double(n);
	z= new double(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			delete [] z;
			delete [] b;
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE

//////////////////////////////////////////////////////////////////////////
//
//	tqli
//
// 	Diagonalizes a real triangular matrix created by tred2 e.g.
//
//	PARAMETERS
//		See book
//

void tqli(double d[], double e[], int n, double **z)
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) nrerror("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

//////////////////////////////////////////////////////////////////////////
//
//	tqli_values_only
//
// 	Diagonalizes a real triangular matrix created by tred2 e.g.
//		N.B. Only use if you are not interested in eigenfunctions
//		It will calculate only eigenvalues
//
//	PARAMETERS
//		See book
//

void tqli_values_only(double d[], double e[], int n, double **z)
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) nrerror("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

//////////////////////////////////////////////////////////////////////////
//
//	tred2
//
// 	Triangulates a real symmetric matrix
//
//	PARAMETERS
//		See book
//


void tred2(double **a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>=1;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[0]=0.0;
	e[0]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=0;i<n;i++) {
		l=i-1;
		if (d[i]) {
			for (j=0;j<=l;j++) {
				g=0.0;
				for (k=0;k<=l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

//////////////////////////////////////////////////////////////////////////
//
//	tred2_values_only
//
// 	Triangulates a real triangular matrix
//
//	PARAMETERS
//		See book
//

void tred2_values_only(double **a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>=1;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[0]=0.0;
	e[0]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=0;i<n;i++) {
//		l=i-1;
//		if (d[i]) {
//			for (j=0;j<=l;j++) {
//				g=0.0;
//				for (k=0;k<=l;k++)
//					g += a[i][k]*a[k][j];
//				for (k=0;k<=l;k++)
//					a[k][j] -= g*a[k][i];
//			}
//		}
		d[i]=a[i][i];
//		a[i][i]=1.0;
//		for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

//////////////////////////////////////////////////////////////////////////
//
//	eigen_sort
//
// 	Sorts the eigenvalues and functions in decsending order,
//    thus starting with the largest value
//

void eigen_sort(double *v, double **a,int n)
{
	double h;
	int i,j,k;

	for (i=0; i<n-1;i++)
		for (j=i+1; j<n;j++)
			if (v[i]<v[j])
			{
				h=v[i];
				v[i]=v[j];
				v[j]=h;
				for (k=0; k<n; k++)
				{
					h = a[k][i];
					a[k][i] = a[k][j];
					a[k][j] = h;
				}
			}
}

//////////////////////////////////////////////////////////////////////////
//
//	pythag
//
//		returns c = root(a^2 + b^2) in a complicated manner
//

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+absb*absb/(absa*absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+absa*absa/(absb*absb)));
}

/* (C) Copr. 1986-92 Numerical Recipes Software Q.36=41,. */

// Now to solve the eigenvalue problem for a hermitian
// Hamiltonian, we use the suggestion in the NR book.
// We construct a real symmetric matrix of order 2n and
// solve this one. We then do find all eigenvalues and
// functions twice, however.
// If A is the real part, B the matrix with the imaginary
// part we make
//  ( A -B )
//  ( B  A )
//  If u + iv is an eigenvector, i(u+iv) = -v+iu will also
// be an eigenvector with the same eigenvalue.
// The eigenvalues and vectors, however, are not ordered and
// one cannot simply take the first half of the eigenvalues.

//////////////////////////////////////////////////////////////////////////
//
//	eigenval
//
//		Calculates eigenvalues of complex Hermitian matrix
//

Vector eigenval(const Ctensor &H)
{
	double **M, *d, *e;
	int n = H.get_order();

	d = new double[2*n];
	e = new double[2*n];

	M = new double*[2*n];
	for (int i=0; i<2*n; i++)
		M[i] = new double[2*n];

	Tensor A = H.realtensor();
	Tensor B = H.imagtensor();

	for ( int i=0; i<n; i++)
		for ( int j=0; j<n; j++)
		{
			M[i][j] = M[i+n][j+n] = A.get(i,j);
			M[i+n][j] = B.get(i,j);
			M[i][j+n] = - M[i+n][j];
		}

	tred2_values_only(M, 2*n, d, e);
	tqli_values_only(d, e, 2*n, M);

	Vector prop(2*n);
	for (int i=0; i<2*n; i++)
		prop.set(i,d[i]);
	return prop;
}


