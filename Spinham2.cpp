#include <stdio.h>
#include <values.h>
#include "complex.h"
#include "vector.h"
#include "tensor.h"
#include "eigen.h"
#include "spin.h"
#include "Spinham2.h"
#include "my_const.h"
#include "my_util.h"

extern "C"
{
float ran1(long* seed);
}

SpinHam::SpinHam(int em, int im)
					:Emult(em), Imult(im)
{
	order = Emult*Imult;
	Ham0.reset(order);
	Ham.reset(order);
	Eivec.reset(order);
	Eival.reset(order);
	Trans.reset((order*(order-1))/2);
	Tprob.reset((order*(order-1))/2);
	g.reset(3);
	A.reset(3);
	Q.reset(3);
	for (int i=0; i<3; i++)
		g.set(i,i, 2.0023);
	B=0.0;
	theta=0.0;
	phi=0.0;
	update = 0;
	gamma = 3.0756;
	A20 = 0.0;
	A22 = 0.0;
	A40 = 0.0;
	A42 = 0.0;
	A43 = 0.0;
	A44 = 0.0;
	A60 = 0.0;
	A63 = 0.0;
	A64 = 0.0;
	A66 = 0.0;
}

void SpinHam::Reset(int em, int im)
{
	Emult = em;
	Imult = im;
	order = Emult*Imult;
	Ham0.reset(order);
	Ham.reset(order);
	Eivec.reset(order);
	Eival.reset(order);
	Trans.reset((order*(order-1))/2);
	Tprob.reset((order*(order-1))/2);
	g.reset(3);
	A.reset(3);
	Q.reset(3);
	for (int i=0; i<3; i++)
		g.set(i,i, 2.0023);
	B=0.0;
	theta=0.0;
	phi=0.0;
	update = 0;
	gamma = 3.0756;
	A20 = 0.0;
	A22 = 0.0;
	A40 = 0.0;
	A42 = 0.0;
	A43 = 0.0;
	A44 = 0.0;
	A60 = 0.0;
	A63 = 0.0;
	A64 = 0.0;
	A66 = 0.0;
}

void SpinHam::set_g_tensor(double xx, double yy, double zz,
					double xy, double xz, double yz)
{
	g.set(0,0, xx);
	g.set(0,1, xy);
	g.set(0,2, xz);
	g.set(1,0, xy);
	g.set(1,1, yy);
	g.set(1,2, yz);
	g.set(2,0, xz);
	g.set(2,1, yz);
	g.set(2,2, zz);
}


void SpinHam::set_Q_tensor(double xx, double yy, double zz,
					double xy, double xz, double yz)
{
	Q.set(0,0, xx);
	Q.set(0,1, xy);
	Q.set(0,2, xz);
	Q.set(1,0, xy);
	Q.set(1,1, yy);
	Q.set(1,2, yz);
	Q.set(2,0, xz);
	Q.set(2,1, yz);
	Q.set(2,2, zz);
}

void SpinHam::set_A_tensor(double xx, double yy, double zz,
					double xy, double xz, double yz)
{
	A.set(0,0, xx);
	A.set(0,1, xy);
	A.set(0,2, xz);
	A.set(1,0, xy);
	A.set(1,1, yy);
	A.set(1,2, yz);
	A.set(2,0, xz);
	A.set(2,1, yz);
	A.set(2,2, zz);
}

void SpinHam::set_field(const Vector& Field)
{
	B = Field.length();
	if (B<1.0e-30) {
		theta=0.0; phi=0.0; return; }
	theta = acos(Field.get(2)/B) * 180.0/PI;
	if ( fabs(Field.get(0)) < 1.0e-30 ) { phi=90.0; return;}
	phi = atan(Field.get(1)/Field.get(0)) * 180.0/PI;
	return;
}

double SpinHam::radical_resonance(double freq) const
{
	double temp, g_eff = 0.0;
	int i,j;
	Vec3 c;
	double t = theta*PI/180.0;
	double p = phi*PI/180.0;
	c.setv(sin(t)*cos(p), sin(t)*sin(p), cos(t));

	for (i=0; i<3; i++)
	{
		temp = 0.0;
		for (j=0; j<3; j++) temp += g.get(i,j)*c[j];
		g_eff += temp*temp;
	}
	g_eff = sqrt(g_eff); // Now we have the effective g-value

	return  freq * 1.0e9 * h_Planck /(mu_B * g_eff) ;
}

Ctensor SpinHam::setHF(Vector SEV)
{
	Spin I(Imult);
	if (Emult != 1) return Ham0;

	complex help;
	int u,v,k,l;
	Ham0.reset(Imult);

// The hyperfine term <S>.A.I

	for (u=0; u<3; u++)		 // sum over SEVx, SEVy, SEVz;
		for ( v=0; v<3; v++)	 // sum over Ix, Iy, Iz;

	for ( k=0; k<Imult; k++)       // sum over elements of I[v]
		for ( l=0; l<Imult; l++)
				{
					help = A.get(u,v) * SEV.get(u) * I[v].get(k,l);
					Ham0.add(k,l,help); // NB in MHz!
				}
// The Quadrupole term I.Q.I
	Ctensor helpt(Imult);

	for ( u=0; u<3; u++)		 // sum over Ix, Iy, Iz;
		for ( v=0; v<3; v++)	 // sum over Ix, Iy, Iz;
			helpt += (I[u]*I[v])*Q.get(u,v);

		for ( k=0; k<Imult; k++)       // sum over elements of I[v]
			for ( l=0; l<Imult; l++)
				Ham0.add(k, l, helpt.get(k,l) ); // NB in MHz

	return Ham0;
}

Ctensor SpinHam::setH0()
{
	int i,j,k,l,u,v;
	complex help;

	Spin S(Emult);
	Spin I(Imult);

// First the hyperfine term  S . A . I
	for (u=0; u<3; u++)		 // sum over Sx, Sy, Sz;
		for (v=0; v<3; v++)	 // sum over Ix, Iy, Iz;

	for (i=0; i<Emult; i++)       // sum over elements of S[u]
		for (j=0; j<Emult; j++)

	for (k=0; k<Imult; k++)       // sum over elements of I[v]
		for (l=0; l<Imult; l++)
				{
					help = A.get(u,v) * S[u].get(i,j) * I[v].get(k,l);
					Ham0.add(i*Imult + k, j*Imult + l, help * 1e-3); // NB GHz!
				}

// Now the Quadrupole term  I . Q . I
	Ctensor helpt(Imult);

	for (u=0; u<3; u++)		 // sum over Ix, Iy, Iz;
		for (v=0; v<3; v++)	 // sum over Ix, Iy, Iz;
			helpt += (I[u]*I[v])*Q.get(u,v);

	for (i=0; i<Emult; i++)       // sum over elements of S[u]
		for (k=0; k<Imult; k++)       // sum over elements of I[v]
			for (l=0; l<Imult; l++)
				Ham0.add(i*Imult + k, i*Imult + l, helpt.get(k,l) * 1e-3); // NB GHz!


// Now the crystal field terms    V = sum A_k^q * O_k^q
//	First of all the equivalence opearators
// O_2^0 and O_2^2 which are  zero if S < 1
// These prefactors for these terms are referred to as the 3*D and E
// parameters. For higher spins than S=1 higher order terms may also
// play a role. Here some of the fourth and sixth order may also
// be included.


	double spin = double(Emult-1)/2.0;    // spin value S
	double J2 = spin*(spin+1);            // eigenvalue of S^2

	if (spin < 0.9) return Ham0;

	if (fabs(A20) > 1.0e-30)
	{
		Ctensor O20 = S[2]*S[2]*3.0 + S[5]*(-J2);
		for (i=0; i<Emult; i++)
//			for (j=0; j<Emult; j++)       NB diagonal in Jz
				for (k=0; k<Imult; k++)
					Ham0.add( i*Imult+k, i*Imult+k, A20*O20.get(i,i) );
	}

	if (fabs(A22) > 1.0e-30)
	{
		Ctensor O22 = (S[3]*S[3] + S[4]*S[4]) * 0.5;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k, j*Imult+k, A22*O22.get(i,j) );
	}

	if (fabs(A40) > 1.0e-30)
	{
		Ctensor O40 = S[2]*S[2]*S[2]*S[2]*35.0 + S[2]*S[2]*(25.0-30.0*J2)
							- S[5]*(6.0*J2 - 3.0*J2*J2);
		for (i=0; i<Emult; i++)
//			for (j=0; j<Emult; j++)  NB Diagonal in Jz
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,i*Imult+k, A40*O40.get(i,i));
	}

	if (fabs(A42) > 1.0e-30)
	{
		Ctensor O42 = ((S[2]*S[2]*7.0 - S[5]*(J2+5)) * (S[3]*S[3]+S[4]*S[4])
					+ (S[3]*S[3]+S[4]*S[4]) * (S[2]*S[2]*7.0 - S[5]*(J2+5)))*0.25;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,j*Imult+k, A42*O42.get(i,j));
	}

	if (fabs(A43) > 1.0e-30)
	{
		Ctensor O43 = (S[2]*(S[3]*S[3]*S[3] + S[4]*S[4]*S[4]) +
								(S[3]*S[3]*S[3] + S[4]*S[4]*S[4])*S[2])*0.25;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,j*Imult+k, A43*O43.get(i,j));
	}

	if (fabs(A44) > 1.0e-30)
	{
		Ctensor O44 = (S[3]*S[3]*S[3]*S[3] + S[4]*S[4]*S[4]*S[4]) * 0.5;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,j*Imult+k, A44*O44.get(i,j));
	}

	if (fabs(A60) > 1.0e-30)
	{
		Ctensor O60 = S[2]*S[2]*S[2]*S[2]*S[2]*S[2]*231.0 +
							S[2]*S[2]*S[2]*S[2]*(735.0 - 315.0*J2) +
							 S[2]*S[2] * (294.0 - 525.0*J2 + 105.0*J2*J2) +
							  S[5]*(0.0 - 60.0*J2 + 40.0*J2*J2 - 5.0*J2*J2*J2);
		for (i=0; i<Emult; i++)
//			for (j=0; j<Emult; j++)  NB Diagonal in Jz
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,i*Imult+k, A60*O60.get(i,i));
	}

	if (fabs(A63) > 1.0e-30)
	{
		Ctensor O63 = ( (S[2]*S[2]*S[2]*11.0 - S[2]*(3*J2+59.0))
								*	(S[3]*S[3]*S[3] + S[4]*S[4]*S[4])
							+ (S[3]*S[3]*S[3] + S[4]*S[4]*S[4])
								* (S[2]*S[2]*S[2]*11.0 - S[2]*(3*J2+59.0)) )*0.25;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,j*Imult+k, A63*O63.get(i,j));
	}

	if (fabs(A64) > 1.0e-30)
	{
		Ctensor O64 = ( (S[2]*S[2]*11.0 - S[5]*(J2+38.0))
								*	(S[3]*S[3]*S[3]*S[3] + S[4]*S[4]*S[4]*S[4])
							+ (S[3]*S[3]*S[3]*S[3] + S[4]*S[4]*S[4]*S[4])
								* (S[2]*S[2]*11.0 - S[5]*(J2+38.0)) )*0.25;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,j*Imult+k, A64*O64.get(i,j));
	}

	if (fabs(A66) > 1.0e-30)
	{
		Ctensor O66 = ( S[3]*S[3]*S[3]*S[3]*S[3]*S[3]
								+ S[4]*S[4]*S[4]*S[4]*S[4]*S[4] )*0.5;
		for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Ham0.add(i*Imult+k,j*Imult+k, A66*O66.get(i,j));
	}
	return Ham0;
}



Ctensor SpinHam::addZeeman(const Spin &S, const Spin &I)
{
	int i,j,k;
	complex help;

	if (( S.get_mult() != Emult ) || ( I.get_mult() != Imult)) return Ham;
	Ctensor eZeeman(S.get_mult());
	Ctensor nZeeman(I.get_mult());
	Vec3 c;
	double t = theta*PI/180.0;
	double p = phi*PI/180.0;
	c.setv(sin(t)*cos(p), sin(t)*sin(p), cos(t));
	for (i=0; i<3; i++)
		{
			nZeeman += (I[i] * c[i]);
			for (j=0; j<3; j++)
							eZeeman += S[i] * g.get(i,j) * c[j];
		}
	eZeeman *= (B * mu_B * 1.0e-9 / h_Planck); // Zeeman splitting in GHz
	if (Emult == 1)
		nZeeman *= (B * gamma);    			// nuclear splitting in MHz
	  else
		nZeeman *= (B * gamma * 1.0e-3);    // nuclear splitting in GHz

	if (Imult == 1)
		Ham = Ham0 + eZeeman;
	  else
		{
			Ctensor Temp(order);
			for (i=0; i<Emult; i++)
				for (j=0; j<Emult; j++)
					for (k=0; k<Imult; k++)
						Temp.set(Imult*i+k, Imult*j+k, eZeeman.get(i,j));
			for (i=0; i<Imult; i++)
				for (j=0; j<Imult; j++)
					for (k=0; k<Emult; k++)
						Temp.add(i+Imult*k, j+Imult*k, nZeeman.get(i,j));

			Ham = Ham0 + Temp;
		}
	update = 0;
	return Ham;
}

void SpinHam::addZeeman2(const Spin &S, const Spin &I)
{
	int i,j,k;
	complex help;

	if (( S.get_mult() != Emult ) || ( I.get_mult() != Imult)) return;
	Ctensor eZeeman(S.get_mult());
	Ctensor nZeeman(I.get_mult());
	Vec3 c;
	double t = theta*PI/180.0;
	double p = phi*PI/180.0;
	c.setv(sin(t)*cos(p), sin(t)*sin(p), cos(t));
	for (i=0; i<3; i++)
		{
			nZeeman += (I[i] * c[i]);
			for (j=0; j<3; j++)
							eZeeman += S[i] * g.get(i,j) * c[j];
		}
	eZeeman *= (B * mu_B * 1.0e-9 / h_Planck); // Zeeman splitting in GHz
	if (Emult == 1)
		nZeeman *= (B * gamma);    			// nuclear splitting in MHz
	  else
		nZeeman *= (B * gamma * 1.0e-3);    // nuclear splitting in GHz

	if (Imult == 1)
		Ham = Ham0 + eZeeman;
	  else
		{
			Ctensor Temp(order);
			for (i=0; i<Emult; i++)
				for (j=0; j<Emult; j++)
					for (k=0; k<Imult; k++)
						Temp.set(Imult*i+k, Imult*j+k, eZeeman.get(i,j));
			for (i=0; i<Imult; i++)
				for (j=0; j<Imult; j++)
					for (k=0; k<Emult; k++)
						Temp.add(i+Imult*k, j+Imult*k, nZeeman.get(i,j));

			Ham = Ham0 + Temp;
		}
	update = 0;
	return;
}

void SpinHam::print() const
{
	cout << "Zero field Hamiltonian" << endl;
	Ham0.print();
	cout << "Total Hamiltonian" << endl;
	Ham.print();
	cout << "Eigenvalues" << endl;
	Eival.print();
	cout << "Eigenfunctions" << endl;
	Eivec.print();
	cout << "Transitions" << endl;
	Trans.print();
	Tprob.print();
}

Vector SpinHam::eigenvec()
{
//	We will calculate the eigenvalues and vectors
// using routines for a real matrix. Therefor the complex matrix
// of order n is put into a real matrix of order 2n. We obtain then
// all eigenvalues twice and all eigenvectors twice. At the end things
// are sorted.

	int j;
	double **M, *d, *e;
	int n = order;

	d = new double[2*n];
	e = new double[2*n];

// make the augmented matrix
	M = new double*[2*n];
	for (int i=0; i<2*n; i++)
		M[i] = new double[2*n];

// take the real and imaginairy parts of the hamiltonian
	Tensor A = Ham.realtensor();
	Tensor B = Ham.imagtensor();

//	construct the augmented matrix
	for ( int i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			M[i][j] = M[i+n][j+n] = A.get(i,j);
			M[i+n][j] = B.get(i,j);
			M[i][j+n] = - M[i+n][j];
		}

// diagonalize the augmented matrix
	tred2(M, 2*n, d, e);
	tqli(d, e, 2*n, M);

// sort the eigenvalues and vectors in descending order
	eigen_sort(d,M,2*n);

// put the results in Eivec and Eival
//	but invert the sorting to start from the lowest eigenvalue
	int index;
	for (int i=0; i<n; i++)
	{
		index = n-1-i;
		Eival.set(index,d[2*i]);
		for (j=0; j<n; j++)
			Eivec.set(j,index,M[j][2*i], M[j+n][2*i]);
	}

    delete[] d;
    delete[] e;
    for (int i=0;i<2*n;i++)
        delete[] M[i];
    delete[] M;

	return Eival;
}

Vector SpinHam::eigenval()
{
//	We will calculate only the eigenvalues
// using routines for a real matrix. Therefor the complex matrix
// of order n is put into a real matrix of order 2n. We obtain then
// all eigenvalues twice and all eigenvectors twice. At the end things
// are sorted.

	int j;
	double **M, *d, *e;
	int n = order;

	d = new double[2*n];
	e = new double[2*n];

// make the augmented matrix
	M = new double*[2*n];
	for (int i=0; i<2*n; i++)
		M[i] = new double[2*n];

// take the real and imaginairy parts of the hamiltonian
	Tensor A = Ham.realtensor();
	Tensor B = Ham.imagtensor();

//	construct the augmented matrix
	for ( int i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			M[i][j] = M[i+n][j+n] = A.get(i,j);
			M[i+n][j] = B.get(i,j);
			M[i][j+n] = - M[i+n][j];
		}

// diagonalize the augmented matrix
	tred2_values_only(M, 2*n, d, e);
	tqli_values_only(d, e, 2*n, M);

// sort the eigenvalues and vectors
	eigen_sort(d,M,2*n);

// put the results in  Eival
// but invert the sorting to start from the lowest eigenvalue
	for (int i=0; i<n; i++)
			Eival.set(n-1-i,d[2*i]);

    delete[] d;
    delete[] e;
    for (int i=0;i<2*n;i++)
        delete[] M[i];
    delete[] M;

	return Eival;
}

Cvector SpinHam::get_eigenvector(int i)
{
	if (i>order) i=order;
	Cvector eigenvector(order);
	for (int j=0; j<order; j++)
		eigenvector.set(j,Eivec.get(j,i));
	return eigenvector;
}

// To calculate the transition probabilities it will
// be assumed that the oscillating field is perpendicular
// to the static magnetic field. Let's (arbitrarily choose
// it such that its z-component is zero.

void SpinHam::transitions(int mode)
{
	Cvector c(3);
	switch (mode)
	{
		case 0:
			c.set(0,sin(phi*PI/180.0), 0.0);
			c.set(1,-cos(phi*PI/180.0), 0.0);
			c.set(2,0.0, 0.0);
// In this way we take the B1 field always perpendicular to the B0 field
// and always in the xy plane. (default)
			break;
		case 1:
			c.set(0,cos(theta*PI/180.0)*cos((phi+180.0)*PI/180.0), 0.0);
			c.set(1,cos(theta*PI/180.0)*sin((phi+180.0)*PI/180.0), 0.0);
			c.set(2,sin(theta), 0.0);
// Here B1 is perpendicular to B0, but in the plane spanned by B0 and z.
			break;
		case 2:
			c.set(0,sin(theta*PI/180.0)*cos((phi)*PI/180.0), 0.0);
			c.set(1,sin(theta*PI/180.0)*sin((phi)*PI/180.0), 0.0);
			c.set(2,cos(theta), 0.0);
// Here B1 is parallel to B0
			break;
		case 3:
			c.set(0, 1.0, 0.0);		//  B1 parallel to x
			c.set(1, 0.0, 0.0);
			c.set(2, 0.0, 0.0);
			break;
		case 4:
			c.set(0, 0.0, 0.0);		//  B1 parallel to y
			c.set(1, 1.0, 0.0);
			c.set(2, 0.0, 0.0);
			break;
		case 5:
			c.set(0, 0.0, 0.0);		//  B1 parallel to z
			c.set(1, 0.0, 0.0);
			c.set(2, 1.0, 0.0);
			break;
		case 6:                      // right circular
			c.set(0,sin(phi*PI/180.0), cos(theta*PI/180.0)*cos((phi+180.0)*PI/180.0));
			c.set(1,-cos(phi*PI/180.0),cos(theta*PI/180.0)*sin((phi+180.0)*PI/180.0));
			c.set(2,0.0, 0.0);
			break;
		case 7:                      // left circular
			c.set(0,sin(phi*PI/180.0), -cos(theta*PI/180.0)*cos((phi+180.0)*PI/180.0));
			c.set(1,-cos(phi*PI/180.0), -cos(theta*PI/180.0)*sin((phi+180.0)*PI/180.0));
			c.set(2,0.0, 0.0);
			break;
		default:
			c.set(0,sin(phi*PI/180.0), 0.0);
			c.set(1,-cos(phi*PI/180.0), 0.0);
			c.set(2,0.0, 0.0);
		// default equal to mode = 0
	}

// We want to have the matrix elements squared:
//	<i| mu_B . B1.g.S  + gamma .B1.I|j>^2
	int i,j,k=0;

// Now let's calculate the complete transition matrix
	Ctensor eZeeman(Emult);
	Ctensor nZeeman(Imult);
	Spin S(Emult);
	Spin I(Imult);
	for (i=0; i<3; i++)
		{
			nZeeman += (I[i] * (gamma*0.001*c[i]/13.99624));
			for (j=0; j<3; j++)
							eZeeman += S[i] * g.get(i,j) * c[j];
		}
	Ctensor Temp(order);
	for (i=0; i<Emult; i++)
		for (j=0; j<Emult; j++)
			for (k=0; k<Imult; k++)
				Temp.set(Imult*i+k, Imult*j+k, eZeeman.get(i,j));
	for (i=0; i<Imult; i++)
		for (j=0; j<Imult; j++)
			for (k=0; k<Emult; k++)
				Temp.add(i+Imult*k, j+Imult*k, nZeeman.get(i,j));
// We now have the disturbance matrix in the original basis
// But electronic and nuclear transitions are not scaled relatively!!

	Ctensor Transit = Eivec.transpose() * ( Temp * Eivec );
// First we calculate the transition energies in GHz
//	In total there are order(order-1)/2 transitions
	k=0;
	for (i=0; i<order-1; i++)
		for (j=i+1; j<order; j++)
		{
			Trans.set(k, Eival.get(j) - Eival.get(i));
			Tprob.set(k, real(conj(Transit.get(i,j))*Transit.get(i,j)));
			k++;
		}
}


int SpinHam::GetRandomState(double T, long *seed) const
{
	float randomnumber = ran1(seed);
	long double sum=0;
	long double pop=0;
	double help;
	for (int i=0; i<order; i++)
	{
		help = -Eival.get(i) * 1.0e9 * h_Planck /(k_B * T);
		if (help > 700.0) help = 700.0;
			sum += expl(help);
	}
	for (int i=0; i<order; i++)
	{
		help = -Eival.get(i) * 1.0e9 * h_Planck /(k_B * T);
		if (help > 700.0) help = 700.0;
		pop += expl(help)/sum;
		if (randomnumber<pop) return i;
	}
	return order;
}

double SpinHam::magnetization(double T) const
{
	int i,j,k;

	Spin S(Emult);
	Ctensor Zeeman(S.get_mult());
	Vec3 c;
	double t = theta*PI/180.0;
	double p = phi*PI/180.0;
	c.setv(sin(t)*cos(p), sin(t)*sin(p), cos(t));
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			Zeeman += S[i] * g.get(i,j) * c[j];

	Ctensor Temp=Zeeman;
	if (Imult != 1)
	{
			Temp.reset(order);
			for (i=0; i<Emult; i++)
			for (j=0; j<Emult; j++)
				for (k=0; k<Imult; k++)
					Temp.set(Imult*i+k, Imult*j+k, Zeeman.get(i,j));
	}
	Ctensor Test = Eivec.transpose() * (Temp * Eivec);

	double pop;
	double sum = 0.0;
	double mag = 0.0;
	for (k =0; k<order; k++)
	{
		pop = exp(- Eival.get(k)* 1.0e9 * h_Planck /(k_B * T));
		sum+=pop;
		mag += pop * real(Test.get(k,k));
	}
	return mag/sum;
}

Vector SpinHam::SpinMoment(int level) const
{
	int i;

	Spin S(Emult);
	Ctensor Mx(S.get_mult());
	Ctensor My(S.get_mult());
	Ctensor Mz(S.get_mult());

	double mx, my, mz;
	for (i=0; i<3; i++)
	{
		Mx += S[i] * g.get(i,0);
		My += S[i] * g.get(i,1);
		Mz += S[i] * g.get(i,2);
	}

	mx = real(Eivec.col(level) * (Mx * Eivec.col(level)));
	my = real(Eivec.col(level) * (My * Eivec.col(level)));
	mz = real(Eivec.col(level) * (Mz * Eivec.col(level)));

	Vector temp(3);
	temp.set(0,mx);
	temp.set(1,my);
	temp.set(2,mz);

	return temp;
}

int SpinHam::read_par()
{
	FILE *par;
	double d;
	int i,j,n=0;
	if ((par = fopen("spinham.par", "r")) == NULL ) return -1;
			// open parameter file

	if ((d = read_double(par)) != MAXDOUBLE) gamma = d;    // read gamma
	if ((d = read_double(par)) != MAXDOUBLE) frequency =d; // read frequency

// read the 6 elements of the symmetric g-tensor
// in the order g_xx, g_yy, g_zz, g_xy, g_xz, g_yz
	for (i=0; i<3; i++)
		if ((d = read_double(par)) != MAXDOUBLE)
			g.set(i,i,d);
		 else n=-1;

	for (i=0; i<3; i++)
		for (j=0; j<i; j++)
			if ((d = read_double(par)) != MAXDOUBLE)
			{
				g.set(i,j,d);
				g.set(j,i,d);
			}
			else n=-1;

// read the 6 elements of the symmetric hyperfine tensor
// in the order A_xx, A_yy, A_zz, A_xy, A_xz, A_yz
	for (i=0; i<3; i++)
		if ((d = read_double(par)) != MAXDOUBLE)
			A.set(i,i,d);
		 else n=-1;

	for (i=0; i<3; i++)
		for (j=0; j<i; j++)
			if ((d = read_double(par)) != MAXDOUBLE)
			{
				A.set(i,j,d);
				A.set(j,i,d);
			}
			else n=-1;
// read the 6 elements of the symmetric Quadrupole tensor
// in the order Q_xx, Q_yy, Q_zz, Q_xy, Q_xz, Q_yz
	for (i=0; i<3; i++)
		if ((d = read_double(par)) != MAXDOUBLE)
			Q.set(i,i,d);
		 else n=-1;

	for (i=0; i<3; i++)
		for (j=0; j<i; j++)
			if ((d = read_double(par)) != MAXDOUBLE)
			{
				Q.set(i,j,d);
				Q.set(j,i,d);
			}
			else n=-1;
	if ((d = read_double(par)) != MAXDOUBLE) A20 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A22 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A40 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A42 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A43 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A44 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A60 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A63 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A64 = d;
	if ((d = read_double(par)) != MAXDOUBLE) A66 = d;

	fclose(par);
	return n;
}

void SpinHam::write_parameters(FILE *f)
{
	if (f == NULL)
	{
		printf("The g-tensor for the effective electron spin: \n");
		printf("  x  ( %f, %f, %f )\n", g.get(0,0), g.get(0,1), g.get(0,2));
		printf("  y  ( %f, %f, %f )\n", g.get(1,0), g.get(1,1), g.get(1,2));
		printf("  z  ( %f, %f, %f )\n", g.get(2,0), g.get(2,1), g.get(2,2));

		printf("\n The hyperfine tensor A (in MHz) : \n");
		printf("  x  ( %f, %f, %f )\n", A.get(0,0), A.get(0,1), A.get(0,2));
		printf("  y  ( %f, %f, %f )\n", A.get(1,0), A.get(1,1), A.get(1,2));
		printf("  z  ( %f, %f, %f )\n", A.get(2,0), A.get(2,1), A.get(2,2));

		printf("\n The quadrupole splitting tensor Q (in MHz) : \n");
		printf("  x  ( %f, %f, %f )\n", Q.get(0,0), Q.get(0,1), Q.get(0,2));
		printf("  y  ( %f, %f, %f )\n", Q.get(1,0), Q.get(1,1), Q.get(1,2));
		printf("  z  ( %f, %f, %f )\n", Q.get(2,0), Q.get(2,1), Q.get(2,2));

		printf("\n The nuclear Zeeman splitting : %f MHz/T", gamma);

		printf("\n The crystal field parameters (in GHz): \n");
		printf("  A20 %f,  A22 %f \n", A20, A22);
		printf("  A40 %f,  A42 %f, A43 %f, A44 %f \n", A40, A42, A43, A44);
		printf("  A60 %f,  A63 %f, A64 %f, A66 %f \n", A60, A63, A64, A66);

	}
	else
	{
		fprintf(f,"The g-tensor for the effective electron spin: \n");
		fprintf(f,"  x  ( %f, %f, %f )\n", g.get(0,0), g.get(0,1), g.get(0,2));
		fprintf(f,"  y  ( %f, %f, %f )\n", g.get(1,0), g.get(1,1), g.get(1,2));
		fprintf(f,"  z  ( %f, %f, %f )\n", g.get(2,0), g.get(2,1), g.get(2,2));

		fprintf(f,"\n The hyperfine tensor A (in MHz) : \n");
		fprintf(f,"  x  ( %f, %f, %f )\n", A.get(0,0), A.get(0,1), A.get(0,2));
		fprintf(f,"  y  ( %f, %f, %f )\n", A.get(1,0), A.get(1,1), A.get(1,2));
		fprintf(f,"  z  ( %f, %f, %f )\n", A.get(2,0), A.get(2,1), A.get(2,2));

		fprintf(f,"\n The quadrupole splitting tensor Q (in MHz) : \n");
		fprintf(f,"  x  ( %f, %f, %f )\n", Q.get(0,0), Q.get(0,1), Q.get(0,2));
		fprintf(f,"  y  ( %f, %f, %f )\n", Q.get(1,0), Q.get(1,1), Q.get(1,2));
		fprintf(f,"  z  ( %f, %f, %f )\n", Q.get(2,0), Q.get(2,1), Q.get(2,2));

		fprintf(f,"\n The nuclear Zeeman splitting : %f MHz/T \n", gamma);

		fprintf(f,"\n The crystal field parameters (in GHz): \n");
		fprintf(f,"  A20 %+f,  A22 %+f \n", A20, A22);
		fprintf(f,"  A40 %+f,  A42 %+f, A43 %+f, A44 %+f \n", A40, A42, A43, A44);
		fprintf(f,"  A60 %+f,  A63 %+f, A64 %+f, A66 %+f \n\n", A60, A63, A64, A66);
	}
}

Vector SpinHam::SpinExpValue(int level)
{
	Spin S(Emult);
	double sx, sy, sz;
	Vector temp(3);

	if (Imult != 1) return temp;
	
	sx = real(Eivec.col(level) *( S[0] * Eivec.col(level)));
	sy = real(Eivec.col(level) *( S[1] * Eivec.col(level)));
	sz = real(Eivec.col(level) *( S[2] * Eivec.col(level)));

	temp.set(0,sx);
	temp.set(1,sy);
	temp.set(2,sz);

	return temp;
}
int SpinHam::SetCF(int i, double B)
{
    switch (i)
    {
        case 0: A20 = B; break;
        case 1: A22 = B; break;
        case 2: A40 = B; break;
        case 3: A42 = B; break;
        case 4: A43 = B; break;
        case 5: A44 = B; break;
        case 6: A60 = B; break;
        case 7: A63 = B; break;
        case 8: A64 = B; break;
        case 9: A66 = B; break;
    }
    return 0;
}

