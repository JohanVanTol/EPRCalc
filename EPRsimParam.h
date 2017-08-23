//---------------------------------------------------------------------------
#ifndef EPRsimParamH
#define EPRsimParamH
//---------------------------------------------------------------------------
class EPRsimParam
{
    public:
        EPRsimParam();
        EPRsimParam(const EPRsimParam& P);
        EPRsimParam& operator=(const EPRsimParam& P);

        int GetNpts() const {return npts;}
        double GetFreq() const {return freq;}
        double GetGaussWidth(int i) const;
        double GetLorentzWidth(int i) const;
        double GetAN(int i) const;
        double GetAngleStep() const {return anglestep;}
        double GetFirst() const {return first;}
        double GetLast() const {return last;}
        double GetModulation() const {return modulation;}

        int SetNpts(int _npts) {npts = _npts;return 0;}
        int SetFreq(double _freq) {freq = _freq;return 0;}
        int SetGaussWidth(int i, double ww);
        int SetLorentzWidth(int i, double ww);
        int SetAngleStep(double AS) {anglestep = AS;}
        int SetFirst(double Fi) {first = Fi;return 0;}
        int SetLast(double Fl) {last = Fl;return 0;}
        int SetModulation(double modul) {modulation = modul;}

        int Read(const char* FileName);
		int Write(const char* FileName);
		int Write(AnsiString* SimParString);

    private:
        int npts;
        double freq;
        double Gwx;
        double Gwy;
        double Gwz;
        double Lwx;
        double Lwy;
        double Lwz;
        double anglestep;
        double first;
        double last;
        double modulation;
};






#endif
