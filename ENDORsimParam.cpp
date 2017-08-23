//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <iostream.h>
#include <fstream.h>
#include "ENDORsimParam.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

PowderEndorSimPar::PowderEndorSimPar()
    :npts(512), freq(285.0),wx(0.01),wy(0.01),wz(0.01),
        ANx(0.0), ANy(0.0), ANz(0.0), IMultN(1), first(-10.0),
        last(10.0), g_start(2.0), g_step(0.01), ng(2), ENDORwidth(0.01),
        anglestep(5.0),modulation(100.0)
{
}

PowderEndorSimPar::PowderEndorSimPar(const PowderEndorSimPar& P)
{
    SetNpts(P.GetNpts());
    SetFreq(P.GetFreq());
    for (int i=0;i<3;i++)
    {
        SetEPRWidth(i,P.GetEPRWidth(i));
        SetAN(i,P.GetAN(i));
    }
    SetMultN(P.GetMultN());
    SetFirst(P.GetFirst());
    SetLast(P.GetLast());
    SetStartg(P.GetStartg());
    SetStepg(P.GetStepg());
    SetNg(P.GetNg());
    SetENDORwidth(P.GetENDORwidth());
    SetAngleStep(P.GetAngleStep());
    SetModulation(P.GetModulation());
}

PowderEndorSimPar& PowderEndorSimPar::operator=(const PowderEndorSimPar& P)
{
    SetNpts(P.GetNpts());
    SetFreq(P.GetFreq());
    for (int i=0;i<3;i++)
    {
        SetEPRWidth(i,P.GetEPRWidth(i));
        SetAN(i,P.GetAN(i));
    }
    SetMultN(P.GetMultN());
    SetFirst(P.GetFirst());
    SetLast(P.GetLast());
    SetStartg(P.GetStartg());
    SetStepg(P.GetStepg());
    SetNg(P.GetNg());
    SetENDORwidth(P.GetENDORwidth());
    SetAngleStep(P.GetAngleStep());
    SetModulation(P.GetModulation());
    return *this;
}

double PowderEndorSimPar::GetEPRWidth(int i) const
{
    switch(i)
    {
        case 0:return wx;break;
        case 1:return wy;break;
        case 2:return wz;break;
        default:return 0.0;break;
    }
}

double PowderEndorSimPar::GetAN(int i) const
{
    switch(i)
    {
        case 0:return ANx;break;
        case 1:return ANy;break;
        case 2:return ANz;break;
        default:return 0.0;break;
    }
}
int PowderEndorSimPar::SetEPRWidth(int i, double ww)
{
    switch(i)
    {
        case 0: wx = ww;break;
        case 1: wy = ww;break;
        case 2: wz = ww;break;
        default:break;
    }
    return 0;
}

int PowderEndorSimPar::SetAN(int i, double A)
{
    switch(i)
    {
        case 0: ANx = A;break;
        case 1: ANy = A;break;
        case 2: ANz = A;break;
        default: break;
    }
    return 0;
}

int PowderEndorSimPar::Read(const char* FileName)
{
    ifstream parfile(FileName);
    if (!parfile) return -1;

    int N;
    double R;

    int L = 80;
    char* line = new char[L];

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (strncmp(line,"Powder ENDOR Parameter file",20) != 0) return -1;

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
            SetEPRWidth(i,R);
          else return -1;
    }

    for (int i=0; i<3;i++)
    {
        parfile.getline(line, L);
        if ((parfile.eof()) || (parfile.fail())) return -1;
        if (sscanf(line,"%lf", &R) == 1)
            SetAN(i,R);
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

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  g_start = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  g_step = R; else return -1;
    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%d", &N) == 1)  ng = N; else return -1;

    parfile.getline(line, L);
    if ((parfile.eof()) || (parfile.fail())) return -1;
    if (sscanf(line,"%lf", &R) == 1)  ENDORwidth = R; else return -1;

    parfile.close();
    return 0;
}

int PowderEndorSimPar::Write(const char* FileName)
{
    ofstream parfile(FileName);
    if (!parfile) return -1;

    parfile << "Powder ENDOR parameter file " << endl;
    parfile << "Comment" << endl;
    parfile << npts << "   number of points in simulated spectrum" << endl;
    parfile << freq << "   EPR frequency (GHz)" << endl;
    for (int i=0; i<3;i++)
            parfile << GetEPRWidth(i) << " EPR linewidth (T)" << i << endl;
    for (int i=0; i<3;i++)
            parfile << GetAN(i) << " Nucleus contributing to EPR spectrum: A" << i << endl;
    parfile << anglestep << "  Step of Angles in degrees" << endl;
    parfile << first << "  Starting Frequency (MHz)" << endl;
    parfile << last << "   Ending frequency(MHz) " << endl;
    parfile << modulation << "  Modulation in kHz (for derivative) " << endl;
    parfile << g_start << "   Start g-value " << endl;
    parfile << g_step << "   g-value step size " << endl;
    parfile << ng << "   number of g-values for calculation " << endl;
    parfile << ENDORwidth << "  ENDOR linewidth (MHz)" << endl;

    parfile.close();
    return 0;
}
