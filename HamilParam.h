//---------------------------------------------------------------------------
#ifndef HamilParamH
#define HamilParamH

#include "Tensor.h"

//---------------------------------------------------------------------------
class HamilPar
{
  public:
    HamilPar();
    HamilPar(const HamilPar& H);

    ~HamilPar();

	HamilPar& operator=(const HamilPar& H);

	int GetMultS() {return MultS;}
	int GetnI() {return nI;}
	int GetMultI(int i);
	Tensor Getg() {return g;}
	Tensor GetA(int i) {return A[i];}
	Tensor GetQ(int i) {return Q[i];}
	double Getgamma(int i) const;
	double GetCF(int i);
	double GetCFStrain(int i);
	double GetgStrain() const {return gStrain;}
	char* GetComment() {return Comment;}

	int SetMultS(int mult) {MultS = mult; return 0;}
	int SetnI(int n) {nI = n; return 0;}
	int SetMultI(int i, int nn);
	int Setg(Tensor gg) {g = gg; return 0;}
	int SetgStrain(double gS) {gStrain = gS;}
	int SetA(int i, Tensor Hyp);
	int SetQ(int i, Tensor Quad);
	int Setgamma(int i, double MHzT);
	int SetCF(int i, double B);
	int SetCFStrain(int i, double Strain);

	int AddNuc();
	int DeleteNuc(int Nuc);

	int Read(const char* FileName);
	int Write(const char* FileName);
	int Write(AnsiString* ParText);
	int Write(ofstream *fitfile);

  private:

	int MultS;
	int nI;
	Tensor g;
	Tensor* A;
	Tensor* Q;
	double* gamma;
	int* MultI;
	double B20;
	double B22;
	double B40;
	double B42;
	double B43;
	double B44;
	double B60;
	double B63;
	double B64;
	double B66;
	double B20Strain;
	double B22Strain;
	double gStrain;
    char* Comment;
};


#endif
