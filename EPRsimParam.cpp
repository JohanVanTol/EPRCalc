//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <iostream.h>
#include <fstream.h>
#include "EPRsimParam.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

EPRsimParam::EPRsimParam()
    :npts(512), freq(285.0),Gwx(0.01),Gwy(0.01),Gwz(0.01),
        Lwx(0.01),Lwy(0.01),Lwz(0.01),
        first(-10.0),last(10.0), anglestep(5.0),modulation(100.0)
{
}

EPRsimParam::EPRsimParam(const EPRsimParam& P)
{
    SetNpts(P.GetNpts());
    SetFreq(P.GetFreq());
    for (int i=0;i<3;i++)
    {
        SetGaussWidth(i,P.GetGaussWidth(i));
        SetLorentzWidth(i,P.GetLorentzWidth(i));
    }
    SetFirst(P.GetFirst());
    SetLast(P.GetLast());
    SetAngleStep(P.GetAngleStep());
    SetModulation(P.GetModulation());
}

EPRsimParam& EPRsimParam::operator=(const EPRsimParam& P)
{
    SetNpts(P.GetNpts());
    SetFreq(P.GetFreq());
    for (int i=0;i<3;i++)
    {
        SetGaussWidth(i,P.GetGaussWidth(i));
        SetLorentzWidth(i,P.GetLorentzWidth(i));
    }
    SetFirst(P.GetFirst());
    SetLast(P.GetLast());
    SetAngleStep(P.GetAngleStep());
    SetModulation(P.GetModulation());
    return *this;
}

double EPRsimParam::GetGaussWidth(int i) const
{
    switch(i)
    {
        case 0:return Gwx;break;
        case 1:return Gwy;break;
        case 2:return Gwz;break;
        default:return 0.0;break;
    }
}

double EPRsimParam::GetLorentzWidth(int i) const
{
    switch(i)
    {
        case 0:return Lwx;break;
        case 1:return Lwy;break;
        case 2:return Lwz;break;
        default:return 0.0;break;
    }
}
int EPRsimParam::SetGaussWidth(int i, double ww)
{
    switch(i)
    {
        case 0: Gwx = ww;break;
        case 1: Gwy = ww;break;
        case 2: Gwz = ww;break;
        default:break;
    }
    return 0;
}

int EPRsimParam::SetLorentzWidth(int i, double ww)
{
    switch(i)
    {
        case 0: Lwx = ww;break;
        case 1: Lwy = ww;break;
        case 2: Lwz = ww;break;
        default:break;
    }
    return 0;
}


int EPRsimParam::Read(const char* FileName)
{
    ifstream parfile(FileName);
    if (!parfile) return -1;

    int N;
    double R;

    int L = 80;
    char* line = new char[L];

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (strncmp(line,"Powder EPR Parameter file",18) != 0) return -1;

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%d", &N) == 1)  npts = N; else return -1;

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  freq = R; else return -1;

    for (int i=0; i<3;i++)
    {
        parfile.getline(line, L);
        if ((parfile.eof()) || (parfile.fail())) return -1;
        if (sscanf(line,"%lf", &R) == 1)
            SetGaussWidth(i,R);
          else return -1;
    }

    for (int i=0; i<3;i++)
    {
        parfile.getline(line, L);
        if ((parfile.eof()) || (parfile.fail())) return -1;
        if (sscanf(line,"%lf", &R) == 1)
            SetLorentzWidth(i,R);
          else return -1;
    }

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  anglestep = R; else return -1;

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  first= R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  last = R; else return -1;

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  modulation = R; else return -1;

    parfile.close();
    return 0;
}

int EPRsimParam::Write(const char* FileName)
{
    ofstream parfile(FileName);
    if (!parfile) return -1;

    parfile << "Powder EPR parameter file " << endl;
    parfile << "Comment" << endl;
    parfile << npts << "   number of points in simulated spectrum" << endl;
    parfile << freq << "   EPR frequency (GHz)" << endl;
    for (int i=0; i<3;i++)
            parfile << GetGaussWidth(i) << " Gaussisn EPR linewidth (T)" << i << endl;
    for (int i=0; i<3;i++)
            parfile << GetLorentzWidth(i) << " Lorentzian EPR linewidth" << i << endl;
    parfile << anglestep << "  Step of Angles in degrees" << endl;
    parfile << first << "  Starting Field(T)" << endl;
    parfile << last << "   Ending Field(T) " << endl;
    parfile << modulation << "  Modulation in kHz (for derivative) " << endl;

    parfile.close();
    return 0;
}

int EPRsimParam::Write(AnsiString* SimParString)
{
	*SimParString = "EPR Fit Parameters \n";
	return SimParString->Length();
}
