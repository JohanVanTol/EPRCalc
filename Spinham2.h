#include "spin.h"
class SpinHam {
	private:
		Ctensor Ham0;
		Ctensor Ham;
		Ctensor Eivec;
		Vector Eival;
		Vector Trans;
		Vector Tprob;
		int Emult;
		int Imult;
		int n_spins;
		int order;
		int update;
		double B;
		double theta;
		double phi;
		Tensor g;                  // electronic Zeeman splitting tensor
		Tensor A;                  // Hyperfine tensor in MHz
		Tensor Q;                  // Quadrupole tensor in MHz
		double gamma;					// Nuclear Zeeman splitting in MHz/tesla
		double frequency;				// Frequency EPR (microwaves/FIR)
		double A20;						// (=D) crystal field parameter (in GHz)
		double A22;                // (=E) crystal field parameter (in GHz)
		double A40;						// A_4^0 crystal field parameter (in GHz)
		double A42;                // A_4^2 crystal field parameter (in GHz)
		double A43;                // A_4^3 crystal field parameter (in GHz)
		double A44;                // A_4^4 crystal field parameter (in GHz)
		double A60;                // A_6^0 crystal field parameter (in GHz)
		double A63;                // A_6^3 crystal field parameter (in GHz)
		double A64;                // A_6^4 crystal field parameter (in GHz)
		double A66;                // A_6^6 crystal field parameter (in GHz)

	public:
		SpinHam(int em=2, int im=1);
		~SpinHam() {};
		void Reset(int em=2, int im=1);
		void set_g_tensor(double xx=2.0023, double yy=2.0023, double zz=2.0023,
					double xy=0.0, double xz=0.0, double yz=0.0);
		void set_A_tensor(double xx=0.0,double yy=0.0, double zz=0.0,
					double xy=0.0, double xz=0.0, double yz=0.0);
		void set_Q_tensor(double xx=0.0,double yy=0.0, double zz=0.0,
					double xy=0.0, double xz=0.0, double yz=0.0);
		void set_field(double field, double _theta=0.0, double _phi=0.0)
						{ B=field; theta=_theta; phi=_phi;}
		void set_field(const Vector& Field);
        int SetCF(int i, double Aij);
		Ctensor setH0();
		Ctensor setHF(Vector SEV);
		Ctensor addZeeman(const Spin &S, const Spin &I);
		void addZeeman2(const Spin &S, const Spin &I);
		Vector eigenvec();
		Vector eigenval();
		Cvector get_eigenvector(int i=0);
		double radical_resonance(double freq = 245.0) const;
		double magnetization(double T) const;
		Vector SpinMoment(int level) const;
		void print() const ;
		int read_par();
		void write_parameters(FILE *f = NULL);
		void transitions(int mode=0);
		Vector get_trans() const {return Trans;}
		Vector get_tprob() const {return Tprob;}
		int GetRandomState(double Temp, long* seed) const;
		Vector SpinExpValue(int level);
        int GetEmult() const {return Emult;}
		int GetOrder() const {return order;}
		void SetGamma(double MHz_p_T) {gamma = MHz_p_T;}
};
