#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <string.h>
#include "Vector.h"
#include "Tensor.h"
#include "my_const.h"

class Spectrum
{
	private :
		Vector *Field;
		Vector *Intens;
		int npts;
		double start;
		double stop;


	public :
		Spectrum(int _npts, double _start=0, double stop=10.0);
		Spectrum(const Spectrum& S);
		~Spectrum();
		void addGauss(double x, double w, double A);
		void addGaussArray(Vector xc, Vector w, Vector A, int n);
		void addLorentz(double x, double w, double A);
		void addLorentzArray(Vector xc, Vector w, Vector A, int n);
//		void write(ofstream* f);
		double getField(int i) const {return Field->get(i);}
		double getIntens(int i) const {return Intens->get(i);}
		void setField(int i, double field) {if (i<npts) Field->set(i,field);}
		void setIntens(int i, double intens) {if (i<npts) Intens->set(i,intens);}
		int getNpts() const {return npts;}
		void Derivative(Spectrum S, double modulation);
		int Derivative(Spectrum);
		void operator=(const Spectrum S);
		double IntPol(double field) const;
};

Spectrum::Spectrum(int _npts, double _start, double _stop)
		: npts(_npts), start(_start), stop(_stop)
{
	Field = new Vector(npts);
	Intens = new Vector(npts);
	for (int i=0; i<npts; i++)
		Field->set(i,start + (stop-start)*double(i)/double(npts-1));
}

Spectrum::Spectrum(const Spectrum& S)
{
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		Intens->set(i,S.getIntens(i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}

Spectrum::~Spectrum()
{
	delete Field;
	delete Intens;
}

void Spectrum::operator=(const Spectrum S)
{
	delete Field;
	delete Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		Intens->set(i,S.getIntens(i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}


void Spectrum::addGauss(double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
		double x = Field->get(i);
		Intens->add(i, (A/w) * exp(-2*(x-xc)*(x-xc)/(w*w)));
	}
}

void Spectrum::addLorentz(double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
		double x = Field->get(i);
		Intens->add(i, (2*A*w)/(PI * 4.0 *((x-xc)*(x-xc)+(w*w))));
	}
}

void Spectrum::addGaussArray(Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww;
	for (i=0; i<npts; i++)
	{
		double x = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			Intens->add(i, (A.get(j)/ww) * exp(-2*(x-xc.get(j))*(x-xc.get(j))/(ww*ww)));
		}
	}
}

void Spectrum::addLorentzArray(Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww;
	for (i=0; i<npts; i++)
	{
		double x = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			Intens->add(i, 2.0*A.get(j)*ww/(PI *(4*(x-xc.get(j))*(x-xc.get(j))+(ww*ww))));
		}
	}
}

int Spectrum::Derivative(Spectrum S)
{
	int i;
	delete Field;
	delete Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	Field->set(0,S.getField(0));
	Intens->set(0,S.getIntens(1)-S.getIntens(0));
	for (i=1; i<npts-1; i++)
	{
		Field->set(i,S.getField(i));
		Intens->set(i,(S.getIntens(i+1)-S.getIntens(i-1))/2.0);
	}
	Field->set(npts-1,S.getField(npts-1));
	Intens->set(0,S.getIntens(npts-1)-S.getIntens(npts-2));
	start = Field->get(0);
	stop = Field->get(npts-1);
	return npts;
}




void Spectrum::Derivative(Spectrum S, double mod)
{
	delete Field;
	delete Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector(npts);
	mod /=1000;
	double B;
	for (int i=0; i<npts; i++)
	{
		B = S.getField(i);
		Field->set(i,B);
		if (S.IntPol(B-mod/2.0) < -0.5)
			setIntens(i,2*(S.IntPol(B+mod/2.0) - S.getIntens(i)));
		  else
		  {
				if (S.IntPol(B+mod/2) < -0.5)
					setIntens(i,2*(S.getIntens(i) - S.IntPol(B-mod/2.0)));
					else
					setIntens(i,S.IntPol(B+mod/2.0) - S.IntPol(B-mod/2.0));
		  }
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
	return;
}



double Spectrum::IntPol(double field) const
{
	if ((field < Field->get(0)) || (field > Field->get(npts-1))) return -1.0;
			// Return negative value if outside range
	int i=0;
	while ((Field->get(i) < field) && (i<npts)) i++;
	if (i >= npts) return -1; //Outside of range
	double q = (field - Field->get(i-1))/(Field->get(i)-Field->get(i-1));
	return Intens->get(i-1) + q * (Intens->get(i) - Intens->get(i-1));
}

class MultSpectrum
{
	private :
		Vector *Field;
		Vector **Intens;
		int nI;
		int npts;
		double start;
		double stop;


	public :
		MultSpectrum(int ng, int _npts, double _start=0, double stop=10.0);
		MultSpectrum(const MultSpectrum& S);
		~MultSpectrum();
		void addGauss(int i, double x, double w, double A);
		void addGaussArray(int i, Vector xc, Vector w, Vector A, int n);
		void addLorentz(int i,double x, double w, double A);
		void addLorentzArray(int i, Vector xc, Vector w, Vector A, int n);
//		void write(ofstream* f);
		double getField(int i) const {return Field->get(i);}
		double getIntens(int ig, int i) const {return Intens[ig]->get(i);}
		void setField(int i, double field) {if (i<npts) Field->set(i,field);}
		void setIntens(int ig, int i, double intens) {if (i<npts) Intens[ig]->set(i,intens);}
		int getNpts() const {return npts;}
		int getNI() const {return nI;}
		void Derivative(MultSpectrum S, double modulation);
		int Derivative(MultSpectrum);
		void operator=(const MultSpectrum S);
		double IntPol(int ig, double field) const;
};

MultSpectrum::MultSpectrum(int _nI, int _npts, double _start, double _stop)
		: nI(_nI), npts(_npts), start(_start), stop(_stop)
{
	Field = new Vector(npts);
	Intens = new Vector*[nI];
	for (int ig=0; ig<nI; ig++)
		Intens[ig] = new Vector(npts);
	for (int i=0; i<npts; i++)
		Field->set(i,start + (stop-start)*double(i)/double(npts-1));
}

MultSpectrum::MultSpectrum(const MultSpectrum& S)
{
	npts = S.getNpts();
	nI = S.getNI();
	Field = new Vector(npts);
	Intens = new Vector*[nI];
	for (int ig=0; ig<nI; ig++)
		Intens[ig] = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		for (ig=0; ig < nI; ig++)
			Intens[ig]->set(i,S.getIntens(ig,i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}

MultSpectrum::~MultSpectrum()
{
	delete Field;
	for (int ig=0; ig<nI; ig++)
		delete Intens[ig];
	delete[] Intens;
}

void MultSpectrum::operator=(const MultSpectrum S)
{
	int ig;
	delete Field;
	for (ig=0; ig<nI; ig++)
		delete Intens[ig];
	delete[] Intens;
	npts = S.getNpts();
	Field = new Vector(npts);
	Intens = new Vector*[nI];
	for (ig=0; ig<nI; ig++)
		Intens[ig] = new Vector(npts);
	for (int i=0; i<npts; i++)
	{
		Field->set(i,S.getField(i));
		for (ig=0; ig<nI; ig++)
			Intens[ig]->set(i,S.getIntens(ig, i));
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
}


void MultSpectrum::addGauss(int ig, double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
		double x = Field->get(i);
		Intens[ig]->add(i, (A/w) * exp(-2*(x-xc)*(x-xc)/(w*w)));
	}
}

void MultSpectrum::addLorentz(int ig, double xc, double w, double A)
{
	for (int i=0; i<npts; i++)
	{
		double x = Field->get(i);
		Intens[ig]->add(i, (2*A*w)/(PI * 4.0 *((x-xc)*(x-xc)+(w*w))));
	}
}

void MultSpectrum::addGaussArray(int ig, Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww;
	for (i=0; i<npts; i++)
	{
		double x = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			Intens[ig]->add(i, (A.get(j)/ww) * exp(-2*(x-xc.get(j))*(x-xc.get(j))/(ww*ww)));
		}
	}
}

void MultSpectrum::addLorentzArray(int ig, Vector xc, Vector w, Vector A, int n)
{
	int i,j;
	double ww;
	for (i=0; i<npts; i++)
	{
		double x = Field->get(i);
		for (j=0; j<n; j++)
		{
			ww = w.get(j);
			Intens[ig]->add(i, 2.0*A.get(j)*ww/(PI *(4*(x-xc.get(j))*(x-xc.get(j))+(ww*ww))));
		}
	}
}

int MultSpectrum::Derivative(MultSpectrum S)
{
	int i, ig;
	if ((npts!= S.getNpts()) || (nI != S.getNI()))
	{
		delete Field;
		for (ig=0; ig<nI; ig++)
			delete Intens[ig];
		delete[] Intens;
		npts = S.getNpts();
		nI = S.getNI();
		Field = new Vector(npts);
		Intens = new Vector*[nI];
		for (ig=0; ig<nI; ig++)
			Intens[ig] = new Vector(npts);
	}
	Field->set(0,S.getField(0));
	for (ig=0; ig<nI; ig++)
		Intens[ig]->set(0,S.getIntens(ig,1)-S.getIntens(ig,0));
	for (i=1; i<npts-1; i++)
	{
		Field->set(i,S.getField(i));
		for (ig=0; ig<nI; ig++)
			Intens[ig]->set(i,(S.getIntens(ig,i+1)-S.getIntens(ig,i-1))/2.0);
	}
	Field->set(npts-1,S.getField(npts-1));
	for (ig=0; ig<nI; ig++)
		Intens[ig]->set(0,S.getIntens(ig, npts-1)-S.getIntens(ig, npts-2));
	start = Field->get(0);
	stop = Field->get(npts-1);
	return npts;
}




void MultSpectrum::Derivative(MultSpectrum S, double mod)
{
	int ig;
	if ((npts!= S.getNpts()) || (nI != S.getNI()))
	{
		delete Field;
		for (ig=0; ig<nI; ig++)
			delete Intens[ig];
		delete[] Intens;
		npts = S.getNpts();
		nI = S.getNI();
		Field = new Vector(npts);
		Intens = new Vector*[nI];
		for (ig=0; ig<nI; ig++)
			Intens[ig] = new Vector(npts);
	}
	mod /=1000;
	double B;
	for (int i=0; i<npts; i++)
	{
		B = S.getField(i);
		Field->set(i,B);
		if (S.IntPol(ig, B-mod/2.0) < -0.5)
		for (ig=0; ig<nI; ig++)
			setIntens(ig,i,2*(S.IntPol(ig,B+mod/2.0) - S.getIntens(ig,i)));
		else
		{
			if (S.IntPol(ig, B+mod/2) < -0.5)
				for (ig=0; ig<nI; ig++)
					setIntens(ig,i,2*(S.getIntens(ig,i) - S.IntPol(ig,B-mod/2.0)));
				else
				for (ig=0; ig<nI; ig++)
					setIntens(ig,i,S.IntPol(ig,B+mod/2.0) - S.IntPol(ig, B-mod/2.0));
		}
	}
	start = Field->get(0);
	stop = Field->get(npts-1);
	return;
}



double MultSpectrum::IntPol(int ig, double field) const
{
	if ((field < Field->get(0)) || (field > Field->get(npts-1))) return -1.0;
			// Return negative value if outside range
	int i=0;
	while ((Field->get(i) < field) && (i<npts)) i++;
	if (i >= npts) return -1; //Outside of range
	double q = (field - Field->get(i-1))/(Field->get(i)-Field->get(i-1));
	return Intens[ig]->get(i-1) + q * (Intens[ig]->get(i) - Intens[ig]->get(i-1));
}


	//void Spectrum::write(ofstream* f)
//{
//	for (int i =0; i<npts; i++)
//	*f << Field->get(i) << " " << Intens->get(i) << endl;
//}

int PowderEndor()
{
	int atom, natoms;
	int k,Imult, npts = 250;
	double mI;
	double Bobs, Bres = 0.0;
	double gxx, gyy, gzz, g_obs;
	double Axx, Ayy, Azz, Axy, Axz, Ayz;
	double modulation,gamma,ENDORwidth;
	double g_eff = 2.0;
	double width = 0.001;
	double freq = 244.996;
	double rad = asin(1.0)/90.0;
	double anglestep= 2.0;
	double first, last;
	double g_start, g_step;
	double aN, ANx, ANy, ANz, intens;
	double wx, wy, wz;
	int ng,ig,j;
	char filename1[13];
	char filename2[13];
	Tensor g(3);
	Tensor *A;

	cout << " Program Powder ENDOR " << endl;

	ifstream in("PIN2.dat");
	if (!in.fail())
	{
		in >> npts;                //Number of points of simulated spectrum
		in >> freq;                //EPR frequency

		in >> gxx;                 //gxx value
		in >> gyy;                 //gyy value
		in >> gzz;                 //gzz value

		in >> wx;                  // width along x, y and z
		in >> wy;
		in >> wz;

		in >> ANx;						// (Nitrogen) hyperfine splitting
		in >> ANy;                 //              along x, y, z
		in >> ANz;

		in >> anglestep;           //size of angle steps

		in >> gamma;               //nuclear Zeeman splitting (MHz/T)

		in >> natoms;
		if ((natoms < 1) || (natoms >20)) natoms=1;

		A = new Tensor[natoms];
		for (atom=0; atom<natoms; atom++)
		{
			A[atom].reset(3);
			in >> Axx;                 //hyperfine Axx in MHz
			in >> Ayy;                 //Ayy in MHz
			in >> Azz;                 //Azz
			in >> Axy;                 //Axy
			in >> Axz;                 //Axz
			in >> Ayz;                 //Ayz
			A[atom].set(0,0,Axx);
			A[atom].set(1,1,Ayy);
			A[atom].set(2,2,Azz);
			A[atom].set(0,1,Axy);
			A[atom].set(1,0,Axy);
			A[atom].set(0,2,Axz);
			A[atom].set(2,0,Axz);
			A[atom].set(1,2,Ayz);
			A[atom].set(2,1,Ayz);
		}

		in >> first;               //First frequency value
		in >> last;                //Last frequency value
		in >> Imult;               //Nuclear spin multiplicity
		in >> modulation;
		in >> g_start;
		in >> g_step;
		in >> ng;
		in >> ENDORwidth;

		in.close();
	}
	  else
	  {
		cout << " Error reading input file PIN2.DAT " << endl;
		cout << " Exit program...";
		exit(0);
	  }

	cout << "Calculating " << ng << " spectra at g-values from " <<
					g_start << " to " << (g_start + ng * g_step) << endl;

	cout << "        g-values " << gxx << ", " << gyy << ", " << gzz << endl;
	cout << "N-hyperfine (MHz)" << ANx << ", " << ANy << ", " << ANz << endl;
	cout << "   linewidths (T)" <<  wx << ", " <<  wy << ", " <<  wz << endl;

//	cout << "Spectrum of " << natoms << " ions with matrices:" << endl;
	for (atom=0; atom < natoms; atom++) A[atom].print();
	MultSpectrum Sum(ng,npts, first, last);
	MultSpectrum Deriv(ng,npts, first, last);
	g.set(0,0,gxx);
	g.set(1,1,gyy);
	g.set(2,2,gzz);

	double a1, a2, a3;
	double l,m,n;
	Vector R(3);
	int MaxVal = 1000;
	Vector *Center;
	Center = new Vector[natoms];
	Vector Intensity(MaxVal);
	Vector Width(MaxVal);
	for (int kk=0; kk< MaxVal; kk++) Width.set(kk,ENDORwidth);
	int count=0;


	cout << "Press any key to start calculation " << endl;
	while (!kbhit()) {}
	double theta = 90.0;
	double phi = 0.0;
	double maxphi = 360.0;
	// Integration over 1/8 of the surface : PI/2
	cout << "estimation number of calculated orientations : "
			<< int(PI/(2.0*anglestep*anglestep*rad*rad)) << endl;

	ofstream out, out2;

	sprintf(filename1,"Spectrum.dat");
	sprintf(filename2,"Derivat.dat");
	out.open(filename1, ios::out);
	out2.open(filename2, ios::out);


	for (ig=0; ig<ng; ig++)
	{

		g_obs = g_start + g_step*ig;
		Bobs = freq * h_Planck * 1.0e9 / (mu_B * g_obs);

		cout << "Calculating at g=" << g_obs
					<< " and " << Bobs << " Tesla" <<endl;

		while (theta > 0)        //theta integrated from 90 to 0
		{
			for (atom=0; atom<natoms; atom++) Center[atom].reset(MaxVal);
			Intensity.reset(MaxVal);
			count = 0;
			while (phi < maxphi)          //Sufficient if A has same principal axes as g
			{
				l = sin(theta*rad)*cos(phi*rad);
				m = sin(theta*rad)*sin(phi*rad);
				n = cos(theta*rad);

				R.set(0,l);
				R.set(1,m);
				R.set(2,n);

// N.B. Width varies due to hyperfine nitrogen. Say 5G in XY plane, 35 G along z
//				aN = ANx*l*l + ANy*m*m + ANz*n*n;
//				aN /= 28000;
				aN = 0;

//				width = sqrt(l*l*wx*wx + m*m*wy*wy + n*n*wz*wz);
				g_eff = sqrt(gxx*gxx*l*l + gyy*gyy*m*m + gzz*gzz*n*n);          //Calculate effective g-value
				Bres = freq * h_Planck * 1.0e9 / (mu_B * g_eff);

				if (kbhit())
				{
					if (getch() == 'S') exit(0);
				}

				if (fabs(Bres-Bobs) < (fabs(aN)+8*width))
				{
//				a1 = (A*R).length();
					for (atom=0; atom < natoms; atom++)
					{
//						a1 = (A[atom]*R)*R;
						a1 = A[atom].get(0,0) *l*l + A[atom].get(0,1) *m*l +
								A[atom].get(0,2) *l*n + A[atom].get(1,0) *m*l + A[atom].get(1,1) *m*m +
								A[atom].get(1,2) *m*n + A[atom].get(2,0) *n*l + A[atom].get(2,1) *m*n +
								A[atom].get(2,2) *n*n;
						Center[atom].set(count,0.5*a1);
					}
					intens = width/(4*(Bres-Bobs)*(Bres-Bobs)+width*width);
//					intens += width/(4*(Bres-Bobs + aN)*(Bres-Bobs +aN)+width*width);
//					intens += width/(4*(Bres-Bobs - aN)*(Bres-Bobs -aN)+width*width);
					Intensity.set(count,intens);
					count++;
				}
				phi += anglestep / sin(theta*rad);
			}
			for (atom=0; atom<natoms;atom++)	Sum.addGaussArray(ig,Center[atom],Width, Intensity, count);
			phi = 0.0;
			theta -= anglestep;
			cout << count << " ";
		}
		cout << endl;
		theta = 90.0;
		phi = 0.0;
	}
	cout << "Taking derivative... " << endl;
	Deriv.Derivative(Sum);
	cout << "Writing to file..." << endl;
	for (j=0; j<npts; j++)
	{
		out << Sum.getField(j) ;
		out2 << Deriv.getField(j);
		for (ig=0; ig < ng; ig++)
		{
			out << "   " << Sum.getIntens(ig, j);
			out2 << "   " << Deriv.getIntens(ig, j);
		}
		out << endl;
		out2 << endl;

	}
	out.close();
	out2.close();

	delete[] Center;
	delete[] A;
	cout << "finished";
	return 0;
}


