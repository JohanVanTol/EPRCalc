//---------------------------------------------------------------------------
#ifndef ENDORsimParamH
#define ENDORsimParamH
//---------------------------------------------------------------------------
class PowderEndorSimPar
{
    public:
        PowderEndorSimPar();
        PowderEndorSimPar(const PowderEndorSimPar& P);
        PowderEndorSimPar& operator=(const PowderEndorSimPar& P);

        int GetNpts() const {return npts;}
        double GetFreq() const {return freq;}
        double GetEPRWidth(int i) const;
        double GetAN(int i) const;
        int GetMultN() const {return IMultN;}
        double GetAngleStep() const {return anglestep;}
        double GetFirst() const {return first;}
        double GetLast() const {return last;}
        double GetModulation() const {return modulation;}
        double GetStartg() const {return g_start;}
        double GetStepg() const {return g_step;}
        int GetNg() const {return ng;}
        double GetENDORwidth() const {return ENDORwidth;}

        int SetNpts(int _npts) {npts = _npts;return 0;}
        int SetFreq(double _freq) {freq = _freq;return 0;}
        int SetEPRWidth(int i, double ww);
        int SetAN(int i, double AA);
        int SetMultN(int _mult) {IMultN = _mult;return 0;}
        int SetAngleStep(double AS) {anglestep = AS;}
        int SetFirst(double Fi) {first = Fi;return 0;}
        int SetLast(double Fl) {last = Fl;return 0;}
        int SetModulation(double modul) {modulation = modul;}
        int SetStartg(double gi) {g_start = gi;return 0;}
        int SetStepg(double gs) {g_step = gs;return 0;}
        int SetNg(int _ng) {ng=_ng;return 0;}
        int SetENDORwidth(double ww) {ENDORwidth = ww; return 0;}

        int Read(const char* FileName);
        int Write(const char* FileName);

    private:
        int npts;
        double freq;
        double wx;
        double wy;
        double wz;
        double ANx;
        double ANy;
        double ANz;
        int IMultN;
        double anglestep;
        double first;
        double last;
        double modulation;
        double g_start;
        double g_step;
        int ng;
        double ENDORwidth;
};






#endif
 