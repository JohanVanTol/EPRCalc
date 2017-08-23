//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop



//---------------------------------------------------------------------------
#pragma package(smart_init)

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include "HamilParam.h"
// This contains the definitions of the hamiltonian data type
//
//
HamilPar::HamilPar()
{
    nI=0;
    MultI = NULL;
    gamma = NULL;
    A = NULL;
    Q = NULL;
    MultS=2;
    Comment = new char[256];
    strcpy(Comment," Default");
    for (int i=0;i<10;i++) this->SetCF(i, 0.0);
    for (int i=0;i<2;i++) this->SetCFStrain(i, 0.0);
    int order = 3;
    g.reset(order);
    for (int i=0;i<order;i++)
    {
        for (int j=0;j<order;j++)
               g.set(i,j,0.0);
        g.set(i,i,2.00232);
	}
	gStrain = 0.0;
}


HamilPar::HamilPar(const HamilPar& H)
{
    MultS = H.GetMultS();
    nI=H.GetnI();
    MultI = new int[nI];
    gamma = new double[nI];
    A = new Tensor[nI];
    Q = new Tensor[nI];
    Comment = new char[256];
    strcpy(Comment, H.GetComment());
    int order = 3;
    g = H.Getg();
    B20 = H.GetCF(0);
    B22 = H.GetCF(1);
    B40 = H.GetCF(2);
    B42 = H.GetCF(3);
    B43 = H.GetCF(4);
    B44 = H.GetCF(5);
    B60 = H.GetCF(6);
    B63 = H.GetCF(7);
    B64 = H.GetCF(8);
    B66 = H.GetCF(9);
    for (int k=0; k<nI; k++)
    {
        A[k] = H.GetA(k);
        Q[k] = H.GetQ(k);
        gamma[k] = H.Getgamma(k);
        MultI[k] = H.GetMultI(k);
    }
    B20Strain = H.GetCFStrain(0);
	B22Strain = H.GetCFStrain(1);
	gStrain = H.GetgStrain();
}

HamilPar& HamilPar::operator=(const HamilPar& H)
{
    delete[] MultI;
    delete[] gamma;
    delete[] A;
    delete[] Q;
    MultS = H.GetMultS();
    nI=H.GetnI();
    MultI = new int[nI];
    gamma = new double[nI];
    A = new Tensor[nI];
    Q = new Tensor[nI];
    strcpy(Comment, H.GetComment());
    int order = 3;
    g = H.Getg();
    B20 = H.GetCF(0);
    B22 = H.GetCF(1);
    B40 = H.GetCF(2);
    B42 = H.GetCF(3);
    B43 = H.GetCF(4);
    B44 = H.GetCF(5);
    B60 = H.GetCF(6);
    B63 = H.GetCF(7);
    B64 = H.GetCF(8);
    B66 = H.GetCF(9);
    for (int k=0; k<nI; k++)
    {
        A[k] = H.GetA(k);
        Q[k] = H.GetQ(k);
        gamma[k] = H.Getgamma(k);
        MultI[k] = H.GetMultI(k);
    }
    B20Strain = H.GetCFStrain(0);
	B22Strain = H.GetCFStrain(1);
	gStrain = H.GetgStrain();
}

HamilPar::~HamilPar()
{
    delete[] MultI;
    delete[] gamma;
    delete[] A;
    delete[] Q;
    delete[] Comment;
}

int HamilPar::GetMultI(int i)
{
	if (i<nI)
	   return MultI[i];
	  else return 1;
}

double HamilPar::Getgamma(int i) const
{
	if ((i>=0) && (i<nI))
	   return gamma[i];
	  else return 0.0;
}

double HamilPar::GetCF(int i)
{
    switch (i)
    {
        case 0: return B20;
        case 1: return B22;
        case 2: return B40;
        case 3: return B42;
        case 4: return B43;
        case 5: return B44;
        case 6: return B60;
        case 7: return B63;
        case 8: return B64;
        case 9: return B66;
        default: return 0.0;
    }
}

double HamilPar::GetCFStrain(int i)
{
    switch (i)
    {
        case 0: return B20Strain;
        case 1: return B22Strain;

        default: return 0.0;
    }
}

int HamilPar::SetCF(int i, double B)
{
    switch (i)
    {
        case 0: B20 = B; break;
        case 1: B22 = B; break;
        case 2: B40 = B; break;
        case 3: B42 = B; break;
        case 4: B43 = B; break;
        case 5: B44 = B; break;
        case 6: B60 = B; break;
        case 7: B63 = B; break;
        case 8: B64 = B; break;
        case 9: B66 = B; break;
    }
    return 0;
}

int HamilPar::SetCFStrain(int i, double Strain)
{
    switch (i)
    {
        case 0: B20Strain = Strain; break;
        case 1: B22Strain = Strain; break;
        default: return -1;
    }
    return 0;
}

int HamilPar::SetA(int iI, Tensor hyp)
{
    if (iI >= nI) return -1;
    A[iI] = hyp;
    return 0;
}

int HamilPar::SetQ(int iI, Tensor Quad)
{
    if (iI >= nI) return -1;
    Q[iI] = Quad;
    return 0;
}

int HamilPar::Setgamma(int i, double gmm)
{
    if (i >= nI) return -1;
	gamma[i] = gmm;
	return 0;
}

int HamilPar::SetMultI(int i, int nn)
{
	if (i >= nI) return -1;
	MultI[i]=nn;
	return 0;
}

int HamilPar::Write(ofstream *somefile)
{
	*somefile << "Parameter File " << endl;
	*somefile << "Comment" << endl;
	*somefile << MultS << "   Spin multiplicity" << endl;
	for (int i=0; i<g.get_order();i++)
		for (int j=i; j<g.get_order();j++)
			*somefile << setprecision(8) << g.get(i,j) << "       g value " << i << j << endl;
	*somefile << B20 << "   B20 " << endl;
	*somefile << B22 << "   B22 " << endl;
	*somefile << B40 << "   B40 " << endl;
	*somefile << B42 << "   B42 " << endl;
	*somefile << B43 << "   B43 " << endl;
	*somefile << B44 << "   B44 " << endl;
	*somefile << B60 << "   B60 " << endl;
	*somefile << B63 << "   B63 " << endl;
	*somefile << B64 << "   B64 " << endl;
	*somefile << B66 << "   B66 " << endl;
	*somefile << B20Strain << "  B20 Strain " << endl;
	*somefile << B22Strain << "  B22 Strain " << endl;
	*somefile << nI << "  number of nuclear spins" << endl;
	for (int k=0; k<nI; k++)
	{
		*somefile << MultI[k] <<  "  Multiplicity nucleus " << k << endl;
		*somefile << gamma[k] <<  "  NMR frequency " << endl;
		for (int i=0; i<A[k].get_order();i++)
			for (int j=i; j<A[k].get_order();j++)
				*somefile << A[k].get(i,j) << "       A" << i << j << endl;
		if (MultI[k] > 2)
		{
			for (int i=0; i<Q[k].get_order();i++)
				for (int j=i; j<Q[k].get_order();j++)
					*somefile << Q[k].get(i,j) << "       Q" << i << j << endl;
		}
	}
	return 0;
}

int HamilPar::Write(const char* FileName)
{
	ofstream parfile(FileName);
	if (!parfile) return -1;

	parfile << "Parameter File " << endl;
	parfile << "Comment" << endl;
	parfile << MultS << "   Spin multiplicity" << endl;
	parfile.precision(8);
	for (int i=0; i<g.get_order();i++)
		for (int j=i; j<g.get_order();j++)
			parfile << g.get(i,j) << "       g value " << i << j << endl;
	parfile << B20 << "   B20 " << endl;
	parfile << B22 << "   B22 " << endl;
	parfile << B40 << "   B40 " << endl;
	parfile << B42 << "   B42 " << endl;
	parfile << B43 << "   B43 " << endl;
	parfile << B44 << "   B44 " << endl;
	parfile << B60 << "   B60 " << endl;
	parfile << B63 << "   B63 " << endl;
	parfile << B64 << "   B64 " << endl;
	parfile << B66 << "   B66 " << endl;
	parfile << B20Strain << "  B20 Strain halfwidth" << endl;
	parfile << B22Strain << "  B22 Strain halfwidth" << endl;
	parfile << nI << "  number of nuclear spins" << endl;
	for (int k=0; k<nI; k++)
	{
		parfile << MultI[k] <<  "  Multiplicity nucleus " << k << endl;
		parfile << gamma[k] <<  "  NMR frequency " << endl;
		for (int i=0; i<A[k].get_order();i++)
			for (int j=i; j<A[k].get_order();j++)
				parfile << A[k].get(i,j) << "       A" << i << j << endl;
		if (MultI[k] > 2)
		{
			for (int i=0; i<Q[k].get_order();i++)
				for (int j=i; j<Q[k].get_order();j++)
					parfile << Q[k].get(i,j) << "       Q" << i << j << endl;
		}
	}
	parfile.close();
	return 0;
}

int HamilPar::Write(AnsiString* ParText)
{
	*ParText = "Hamiltonian Parameter File \n";
	ParText->cat_printf("Spin multiplicity %d \n", MultS);
	for (int i=0; i<g.get_order();i++)
		for (int j=i; j<g.get_order();j++)
			ParText->cat_printf(" g_ %d %d   %f \n", i, j, g.get(i,j));
	ParText->cat_printf("B20   %f \n", B20);
	ParText->cat_printf("B22   %f \n", B22);
	ParText->cat_printf("B40   %f \n", B40);
	ParText->cat_printf("B42   %f \n", B42);
	ParText->cat_printf("B43   %f \n", B43);
	ParText->cat_printf("B44   %f \n", B44);
	ParText->cat_printf("B60   %f \n", B60);
	ParText->cat_printf("B63   %f \n", B63);
	ParText->cat_printf("B64   %f \n", B64);
	ParText->cat_printf("B66   %f \n", B66);
	ParText->cat_printf("B20Strain   %f \n", B20Strain);
	ParText->cat_printf("B22Strain   %f \n", B22Strain);

	ParText->cat_printf("Number of Nuclear Spins  %d \n", nI);
	for (int k=0; k<nI; k++)
	{
		ParText->cat_printf("Multiplicity nucleus %d , NMR Frequency %f \n", MultI[k], gamma[k]);
		for (int i=0; i<A[k].get_order();i++)
			for (int j=i; j<A[k].get_order();j++)
				ParText->cat_printf(" A_ %d %d   %f \n", i, j, A[k].get(i,j));
		if (MultI[k] > 2)
		{
			for (int i=0; i<Q[k].get_order();i++)
				for (int j=i; j<Q[k].get_order();j++)
					ParText->cat_printf(" Q_ %d %d   %f \n", i, j, Q[k].get(i,j));
		}
	}
	return ParText->Length();
}


int HamilPar::Read(const char* FileName)
{
	ifstream parfile(FileName);
	if (!parfile) return -1;

	int N;
	double R;

    int L = 80;
    char* line = new char[L];
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (strncmp(line,"Parameter File",14) != 0) return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    strcpy(Comment, line);
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%d", &N) == 1)  MultS = N; else return -1;
    for (int i=0; i<g.get_order();i++)
        for (int j=i; j< g.get_order();j++) //  read in order xx, xy, xz, yy, yz,zz
        {
            parfile.getline(line, L);
            if ((parfile.eof()) || (parfile.fail())) return -1;
            if (sscanf(line,"%lf", &R) == 1)
            {
                g.set(i,j,R)  ;
                g.set(j,i,R)  ;
            }
                else return -1;
        }
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B20 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B22= R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B40= R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B42 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B43 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B44 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B60 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B63 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B64 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  B66 = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
	if (sscanf(line,"%lf", &R) == 1)  B20Strain = R; else return -1;
	parfile.getline(line, L);
	if ((parfile.eof()) || (parfile.fail())) return -1;
	if (sscanf(line,"%lf", &R) == 1)  B22Strain = R; else return -1;

    if (A != NULL) delete[] A;
    if (Q != NULL) delete[] Q;
    if (MultI != NULL) delete[] MultI;
    if (gamma != NULL) delete[] gamma;
    A = NULL;
    Q = NULL;
    MultI = NULL;
    gamma = NULL;

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%d", &N) == 1)  nI = N; else return -1;
    A = new Tensor[nI];
    Q = new Tensor[nI];
    MultI = new int[nI];
	gamma = new double[nI];

    for (int k = 0; k<nI; k++)
    {
        A[k].reset(3);
        Q[k].reset(3);

        parfile.getline(line, L);
        if ((parfile.eof()) || (parfile.fail())) return -1;
        if (sscanf(line,"%d", &N) == 1)  MultI[k] = N; else return -1;

        parfile.getline(line, L);
        if ((parfile.eof()) || (parfile.fail())) return -1;
        if (sscanf(line,"%lf", &R) == 1)  gamma[k] = R; else return -1;

        for (int i=0; i<A[k].get_order();i++)
            for (int j=i; j<A[k].get_order();j++)
            {
                parfile.getline(line, L);
                if ((parfile.eof()) || (parfile.fail())) return -1;
                if (sscanf(line,"%lf", &R) == 1)
                {
                 A[k].set(i,j,R);
                 A[k].set(j,i,R);
                }
                 else return -1;
            }
        if (MultI[k] > 2)
        {
            for (int i=0; i<Q[k].get_order();i++)
                for (int j=i; j<Q[k].get_order();j++)
                {
                    parfile.getline(line, L);
                    if ((parfile.eof()) || (parfile.fail())) return -1;
                    if (sscanf(line,"%lf", &R) == 1)
                    {
                     Q[k].set(i,j,R);
                     Q[k].set(j,i,R);
                    }
                        else return -1;
                }
        }
    }
    parfile.close();
    return 0;
}

int HamilPar::AddNuc()
{
    Tensor* NewA = new Tensor[nI+1];
    Tensor* NewQ = new Tensor[nI+1];
    double* NewGamma = new double[nI+1];
    int* NewMultI = new int[nI+1];

    for (int k=0; k<nI; k++)
    {
        NewA[k] = A[k];
        NewQ[k] = Q[k];
        NewMultI[k] = MultI[k];
        NewGamma[k] = gamma[k];
    }

    delete[] A;
    delete[] Q;
    delete[] gamma;
    delete[] MultI;

    A = new Tensor[nI+1];
    Q = new Tensor[nI+1];
    gamma = new double[nI+1];
    MultI = new int[nI+1];

    for (int k=0; k<nI; k++)
    {
        A[k] = NewA[k];
        Q[k] = NewQ[k];
        MultI[k] = NewMultI[k];
        gamma[k] = NewGamma[k];
    }
    A[nI].reset(3);
    Q[nI].reset(3);
    gamma[nI] = 0.0;
    MultI[nI] = 2;

    nI++;
    delete[] NewGamma;
    delete[] NewMultI;
    delete[] NewA;
    delete[] NewQ;

    return nI;
}

int HamilPar::DeleteNuc(int Nuc)
{
    if ((Nuc < 0) || (Nuc >= nI)) return -1;

    if (Nuc < nI-1)  // if not the last nucleus
    {
        for (int k=Nuc; k<nI-1; k++)
        {
            A[k] = A[k+1];
            Q[k] = Q[k+1];
            MultI[k] = MultI[k+1];
            gamma[k] = gamma[k+1];
        }
    }
    nI -= 1;
    return nI;
}



