//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "vector.h"
#include "tensor.h"
#include "EPRCalcMain.h"
#include "AboutDialog.h"
#include "TensorDlg.h"
#include "Hamilton.h"
#include "HamilParam.h"
#include "Spectrum.h"
#include "PowderEndorDlg.h"
#include "Eigen.h"
#include "Mydata20.h"
#include "Myplot20.h"
#include "EPRSimDialog.h"
#include "ErrorBx.h"
#include "EPRsimParam.h"
#include "Spinham2.h"
#include "Valid.h"
#include "OptionsForm.h"
#include "FitCycleDlg.h"
#include "StartDialog.h"
#include "Lorentz.h"
#include "Gaussian.h"
#include "FITCYCLE.h"
#include "decays.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TMainForm *MainForm;
HamilPar HPar;
HamilPar TempPar;
PowderEndorSimPar PendorPar;
DataArray *ExpData;
DataArray *SimData;
DataArray *AbsDispData;
MyPlot *SimPlot;
MyPlot *DiffPlot;
bool SimCheck = false;

double **Field = NULL;
double **Ampl = NULL;
double *WW = NULL;

int nangle;
double sig2 = 1.0;

bool HamilChanged = true;
bool WidthChanged = true;
bool PopChanged = true;
bool DisPersChanged = true;
bool NumberChanged = true;

int nTrans = 2;
int PopMode = 0;
int DerivMode = 0;
int LineFitMode = 0;



//---------------------------------------------------------------------------
__fastcall TMainForm::TMainForm(TComponent* Owner)
  : TForm(Owner)
{
}

//---------------------------------------------------------------------------
void __fastcall TMainForm::FileExitClick(TObject *Sender)
{
  Close();
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::AboutClick(TObject *Sender)
{
  AboutBox->ShowModal();
}
//
void __fastcall TMainForm::ParamHamiltonClick(TObject *Sender)
{
  HamiltonianDialog->SetParameters(HPar);
  if (HamiltonianDialog->ShowModal() == mbOK)
  {
        HPar = HamiltonianDialog->GetParameters();
        HamilChanged = true;
  }
}
//---------------------------------------------------------------------------


void __fastcall TMainForm::PowderEndorMenuItemClick(TObject *Sender)
{
    if (PowderENDORDialog->ShowModal() == mbOK)
    {
        PendorPar = PowderENDORDialog->GetPar();
        PowderEndor();
    }
}
//---------------------------------------------------------------------------


int TMainForm::PowderEndor()
{
    char line[80];
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
    double hypshift;

	OutputMemo->Lines->Add(" Program Powder ENDOR ");

	npts = PendorPar.GetNpts();                //Number of points of simulated spectrum
	freq = PendorPar.GetFreq();                //EPR frequency

    wx = PendorPar.GetEPRWidth(0);                  // width along x, y and z
    wy = PendorPar.GetEPRWidth(1);
    wz = PendorPar.GetEPRWidth(2);

    ANx = PendorPar.GetAN(0);
    ANy = PendorPar.GetAN(1);
    ANz = PendorPar.GetAN(2);
    int MultCorr = PendorPar.GetMultN();

    anglestep = PendorPar.GetAngleStep();           //size of angle steps
    natoms = HPar.GetnI();
    A = new Tensor[natoms];
	for (atom=0; atom<natoms; atom++)
	{
		A[atom].reset(3);
        A[atom] = HPar.GetA(atom);
	}
    first = PendorPar.GetFirst();               //First frequency value
	last = PendorPar.GetLast();                //Last frequency value
	Imult =PendorPar.GetMultN();               //Nuclear spin multiplicity
	modulation = PendorPar.GetModulation();
	g_start = PendorPar.GetStartg();
	g_step  = PendorPar.GetStepg();
	ng = PendorPar.GetNg();
	ENDORwidth = PendorPar.GetENDORwidth();

    sprintf(line,"Calculating %d spectra at g-values from %f to %f",
                    ng, g_start, g_start+ng*g_step);
	OutputMemo->Lines->Add(line);
	sprintf(line,"N-hyperfine (MHz) %f, %f, %f", ANx,ANy,ANz);
	OutputMemo->Lines->Add(line);
	sprintf(line,"   linewidths (T)  %f, %f, %f" ,wx,wy, wz);
	OutputMemo->Lines->Add(line);
//	cout << "N-hyperfine (MHz)" << ANx << ", " << ANy << ", " << ANz << endl;
//	cout << "   linewidths (T)" <<  wx << ", " <<  wy << ", " <<  wz << endl;

	sprintf(line,"Spectrum of %d nuclei with matrices:", natoms);
	OutputMemo->Lines->Add(line);
//	cout << "Spectrum of " << natoms << " ions with matrices:" << endl;
	for (atom=0; atom < natoms; atom++) A[atom].print();
	MultSpectrum Sum(ng,npts, first, last);
	MultSpectrum Deriv(ng,npts, first, last);
    g = HPar.Getg();

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


	double theta = 90.0;
	double phi = 0.0;
	double maxphi = 360.0;
	// Integration over 1/8 of the surface : PI/2
//	cout << "estimation number of calculated orientations : "
//			<< int(PI/(2.0*anglestep*anglestep*rad*rad)) << endl;
    double nfactor = PI/(2.0*anglestep*anglestep*rad*rad);

	ofstream out, out2;

	sprintf(filename1,"Spectrum.dat");
	sprintf(filename2,"Derivat.dat");
	out.open(filename1, ios::out);
	out2.open(filename2, ios::out);


	for (ig=0; ig<ng; ig++)
	{

		g_obs = g_start + g_step*ig;
		Bobs = freq * h_Planck * 1.0e9 / (mu_B * g_obs);

		sprintf(line,"Calculating at g=%f and %f Tesla", g_obs ,Bobs);
        OutputMemo->Lines->Add(line);
		while (theta > 0)        //theta integrated from 90 to 0
		{
            ThetaLabel->Caption = theta;
            ThetaLabel->Repaint();
 			for (atom=0; atom<natoms; atom++) Center[atom].reset(MaxVal);
			Intensity.reset(MaxVal);
			count = 0;
			while (phi < maxphi)          //Sufficient if A has same principal axes as g
			{
//                PhiLabel->Caption = phi;
//                PhiLabel->Repaint();
				l = sin(theta*rad)*cos(phi*rad);
				m = sin(theta*rad)*sin(phi*rad);
				n = cos(theta*rad);

				R.set(0,l);
				R.set(1,m);
				R.set(2,n);

// N.B. Width varies due to hyperfine nitrogen. Say 5G in XY plane, 35 G along z

//				width = sqrt(l*l*wx*wx + m*m*wy*wy + n*n*wz*wz);

                gxx = g.get(0,0);
                gyy = g.get(1,1);
                gzz = g.get(2,2);

				g_eff = sqrt(gxx*gxx*l*l + gyy*gyy*m*m + gzz*gzz*n*n);          //Calculate effective g-value
				Bres = freq * h_Planck * 1.0e9 / (mu_B * g_eff);

				aN = ANx*l*l + ANy*m*m + ANz*n*n;
				aN /= (14000*g_eff);

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
                    intens = 0.0;
                    for (int mi2 = -(MultCorr-1); mi2<(MultCorr);mi2+=2)
                    {
                        hypshift = double(mi2)*aN / 2.0;
                    	intens += width/(4*(Bres-Bobs+hypshift)*(Bres-Bobs+hypshift)+width*width);
                    }
//					intens += width/(4*(Bres-Bobs + aN)*(Bres-Bobs +aN)+width*width);
//					intens += width/(4*(Bres-Bobs - aN)*(Bres-Bobs -aN)+width*width);
                    intens /= nfactor; //correct for the amount of angles calculated
					Intensity.set(count,intens);
					count++;
				}
                phi += anglestep / sin(theta*rad);
			}
			for (atom=0; atom<natoms;atom++) Sum.addGaussArray(ig,Center[atom],Width, Intensity, count);
			phi = 0.0;
			theta -= anglestep;
//			cout << count << " ";
		}
//		cout << endl;
		theta = 90.0;
		phi = 0.0;
	}
    OutputMemo->Lines->Add(" Taking derivative... ");
	Deriv.Derivative(Sum);
    OutputMemo->Lines->Add(" Writing to file ");
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
    OutputMemo->Lines->Add(" Finished... ");
	return 0;
}





void __fastcall TMainForm::PowderEPRMenuItemClick(TObject *Sender)
{
    PopMode = 0;
    DerivMode = 1;
    LineFitMode = 0;
    NumberChanged = true;  //always do calculation
    EPRSimulationDialog->SetHamilChange(HamilChanged);
    EPRSimulationDialog->SetPopChange(PopChanged);
    EPRSimulationDialog->SetWidthChange(WidthChanged);
    EPRSimulationDialog->SetNumberChange(NumberChanged);
    EPRSimulationDialog->SetDisPersChange(DisPersChanged);

    EPRSimulationDialog->BoltzmannCheckBox->Checked =true;

    if (EPRSimulationDialog->ShowModal() == mbOK)
    {
// first check which parameters changed
        PopChanged = EPRSimulationDialog->GetPopChange();
        NumberChanged = EPRSimulationDialog->GetNumberChange();
        WidthChanged = EPRSimulationDialog->GetWidthChange();
        DisPersChanged = EPRSimulationDialog->GetDisPersChange();
        HamilChanged = EPRSimulationDialog->GetHamilChange();
    }
    else return;
    if (EPRSimulationDialog->IntegrRadioButton->Checked == true)
        DerivMode = 0; else DerivMode = 1;
// Set the parameters of the SpinHamiltonian
    Update();
    double chi2 = CalcCycle();
    OutputMemo->Lines->Add("done...chi^2:");
    OutputMemo->Lines->Add(chi2);

    Invalidate();

	// Integration over 1/8 of the surface : PI/2
    // The total number of points is 4*PI/(anglestep)^2 over the whole
    // surface.  We'll only ever need to integrate over half of that
    // but for same canonical directions of axes (both g and ZFS collinear)
    // we'll need about half of that.

    HamilChanged = false;
    NumberChanged = false;
    PopChanged = false;
    DisPersChanged = false;
    WidthChanged = false;
}

//---------------------------------------------------------------------------
//  Read Experimental Data file and store in ExpData (global)
//
int TMainForm::ReadExpData(AnsiString Filename)
{
    ifstream data(Filename.c_str());    // open file
    if (!data) return -1;
    if (data.eof()) return 0;

    int ndata = 0;
    const int maxline = 512;
    char* line = new char[maxline];
    DataPoint P;
    DataArray *TempDat = NULL;
   	double x;
	double y[MAXNY];
	int nItems;

    //  Get the first line
    data.getline(line, maxline);
    if (data.fail())
    {
        delete[] line;
        return -1;
    }
    if ((strncmp(line, "Measurement made ",17) == 0) ||
             (strncmp(line, "DataFile HansFormat",19) == 0))
    {
        //  This should be a good measurement file
        //  which includes a header. Let's read the header.
        int StopReadHeader = 0;
        while ((!data.eof()) && (StopReadHeader==0))
        {
            data.getline(line, maxline);
            if (data.fail())
            {
                delete[] line;
                data.close();
                return -1;
            }
            if (strncmp(line, ">BEGIN",6) == 0) StopReadHeader = 1;
                            //Serves as a more elegant break
        }
        if (data.eof()) return 0;

        // Now read the data
        while ((!data.fail()) && (!data.eof()))
        {
            data >> P;
            if ( (!data.eof()) && (!data.fail()))
	        {
	            if ((TempDat == NULL) || (TempDat->Getn() == 0))    // if first datapoint
                {
                    if (TempDat != NULL) delete TempDat;
                    if ((P.GetNy() > 10) && (P.GetNy() <=20))
                        TempDat = new DataArray(16396,P.GetNy());
                      else
                        if (P.GetNy() >20)
                        TempDat = new DataArray(8198,P.GetNy());
                        TempDat = new DataArray(MAXPNT,P.GetNy());
                }
	        TempDat->Add(P);
            ndata++;
    	    }
        }
        if (TempDat != NULL)
        {
            *ExpData = *TempDat;
            delete TempDat;          // cleanup
            TempDat = NULL;
        }
        delete[] line;
        data.close();
        return ndata;
    }
    // If not, let's try if we have a good ASCII data file
	nItems = sscanf(line,
		"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&x, y, y+1, y+2, y+3, y+4, y+5, y+6, y+7, y+8, y+9, y+10, y+11, y+12,
			y+13,y+14, y+15, y+16, y+17, y+18, y+19, y+20);
	if ((strlen(line) > 2) && (nItems > 1) && (nItems <= MAXNY+1))
    {
        // reset the pointer
        data >> P;
        delete ExpData;
        ExpData = new DataArray(MAXPNT,P.GetNy());
        data.seekg(0);
        while ((!data.fail()) && (!data.eof()))
        {
            data >> P;
	    	if ( (!data.eof()) && (!data.fail()))
		    {
			    ExpData->Add(P);
                ndata++;
    		}
	    }
	}
    delete[] line;
    data.close();
    return ndata;
}

void __fastcall TMainForm::Open1Click(TObject *Sender)
{
    if (OpenDialog->Execute())
    {
        int npts = ReadExpData(OpenDialog->FileName);
        if (npts < 1)
        {
            ErrorBox->ErrorMessage->Caption = "Not a valid Data File";
            ErrorBox->ShowModal();
            return;
        }
    }  else return;

//  The datafile has been read.
//  Set the plotparameters
    FileLabel->Caption = OpenDialog->FileName;
    int ncol = ExpData->GetNy();

    int curcol = ColumnComboBox->ItemIndex;

    SimPlot->ResetTraces();  // Set traces to plot to zero
    SimPlot->AddTrace(0,0);  // first column is X-axis

    if ((curcol > 0) && (curcol <= ncol))
        SimPlot->AddTrace(1,curcol);
      else     SimPlot->AddTrace(1,1);

    ColumnComboBox->Items->Clear();
    ColumnComboBox->Items->Add("0 (Sim)");
    for (int i=1; i<= ncol; i++)
        ColumnComboBox->Items->Add(AnsiString(i));
    if ((curcol>0) && ( curcol <= ncol) )
        ColumnComboBox->ItemIndex = curcol;
      else ColumnComboBox->ItemIndex = 0;

//    Use the first some points to get an idea of the noise
    int stop;

    if (!ValidInt(OptionsDialog->NoisePointsEdit->Text,&stop))
        stop = 12;
    if (stop > ExpData->Getn()) stop = ExpData->Getn();


    double S, S2, Help;
    S=S2=0.0;
    for (int i=0; i<stop; i++)
    {
        Help = ExpData->Get(i).Get(curcol);
        S  += Help;
        S2 += Help*Help;
    }
    S /= double(stop);
    S2 /= double(stop);
    sig2 = S2 - S*S;

    SimPlot->SetAutoRange(7);
    SimCheck = false;
    Invalidate();
    PlotData();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FormCreate(TObject *Sender)
{
    ExpData = new DataArray(100,10);
    SimData = new DataArray(100,10);
    SimPlot = new MyPlot();
    DiffPlot = new MyPlot();
    NumberChanged = true;
}
//---------------------------------------------------------------------------


void __fastcall TMainForm::FormDestroy(TObject *Sender)
{
    delete ExpData;
    delete SimData;
    delete DiffPlot;
    delete SimPlot;
}
//---------------------------------------------------------------------------
void TMainForm::PlotData()
{
    if (SimPlot == NULL) return;
    int TopWin = Panel2->Top + Panel2->Height;
    int PlotDivide = TopWin  + ceil((ClientHeight-TopWin)*0.67) ;

    SimPlot->SetLim(0, TopWin, ClientWidth - OutputMemo->Width, PlotDivide);
    DiffPlot->SetLim(0, PlotDivide,ClientWidth - OutputMemo->Width, ClientHeight);

    int ccol;
    int ncol = ExpData->GetNy();

    if ((ExpData->Getn() >0) && (!SimCheck))
	{
        ccol = ColumnComboBox->ItemIndex;

        SimPlot->ResetTraces();  // Set traces to plot to zero
        SimPlot->AddTrace(0,0);  // first column is X-axis

        if ((ccol > 0) && (ccol <= ncol))
            SimPlot->AddTrace(1,ccol);
//        else     SimPlot->AddTrace(1,1);

        SimPlot->SetRanges(ExpData->Minima(), ExpData->Maxima());
        SimPlot->PlotTheAxes(Canvas);
        SimPlot->PlotArray(Canvas, *ExpData);
    }

    if ((SimData->Getn() >0) && (SimCheck))
	{
        SimPlot->ResetTraces();  // Set traces to plot to zero
        SimPlot->AddTrace(0,0);  // first column is X-axis
        SimPlot->AddTrace(1,1);
        SimPlot->AddTrace(1,2);
        SimPlot->AddTrace(1,4);
        SimPlot->SetRanges(SimData->Minima(), SimData->Maxima());
        SimPlot->PlotTheAxes(Canvas);
        SimPlot->PlotArray(Canvas, *SimData);

        DiffPlot->ResetTraces();
        DiffPlot->AddTrace(0,0);  // first column is X-axis
        DiffPlot->AddTrace(1,3);
        DiffPlot->SetRanges(SimData->Minima(), SimData->Maxima());
        DiffPlot->PlotTheAxes(Canvas);
        DiffPlot->PlotArray(Canvas, *SimData);

    }
}


void __fastcall TMainForm::FormResize(TObject *Sender)
{
    Invalidate();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FormPaint(TObject *Sender)
{
    PlotData();
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::LorentzFit1Click(TObject *Sender)
{
    SpectralLineStartDialog->Show();
    LineFitMode = 1;
    return;
}

void __fastcall TMainForm::PowderTransMenuItemClick(TObject *Sender)
{

    PopMode = 1;
    DerivMode = 0;
    LineFitMode = 0;
    EPRSimulationDialog->SetHamilChange(HamilChanged);
    EPRSimulationDialog->SetPopChange(PopChanged);
    EPRSimulationDialog->SetWidthChange(WidthChanged);
    EPRSimulationDialog->SetNumberChange(NumberChanged);
    EPRSimulationDialog->SetDisPersChange(DisPersChanged);

    EPRSimulationDialog->BoltzmannCheckBox->Checked =false;
    if (EPRSimulationDialog->ShowModal() == mbOK)
    {
// first check which parameters changed
        PopChanged = EPRSimulationDialog->GetPopChange();
        NumberChanged = EPRSimulationDialog->GetNumberChange();
        WidthChanged = EPRSimulationDialog->GetWidthChange();
        DisPersChanged = EPRSimulationDialog->GetDisPersChange();

    }
    else return;

// Set the parameters of the SpinHamiltonian
    Update();
    double chi2 = CalcCycle();
    OutputMemo->Lines->Add("done...chi^2:");
    OutputMemo->Lines->Add(chi2);

    Invalidate();

	// Integration over 1/8 of the surface : PI/2
    // The total number of points is 4*PI/(anglestep)^2 over the whole
    // surface.  We'll only ever need to integrate over half of that
    // but for same canonical directions of axes (both g and ZFS collinear)
    // we'll need about half of that.

    HamilChanged = false;
    NumberChanged = false;
    PopChanged = false;
    DisPersChanged = false;
    WidthChanged = false;
}



//---------------------------------------------------------------------------


void __fastcall TMainForm::ColumnComboBoxChange(TObject *Sender)
{
    int *columns;
    columns = new int[1];
    columns[0] = ColumnComboBox->ItemIndex;
    SimPlot->SetDataColumns(1,columns,1);

    int stop;
    if (!ValidInt(OptionsDialog->NoisePointsEdit->Text,&stop))
        stop = 12;
    if (stop > ExpData->Getn()) stop = ExpData->Getn();

    double S, S2, Help;
    S = S2 = 0.0;
    for (int i=0; i<stop; i++)
    {
        Help = ExpData->Get(i).Get(columns[0]);
        S  += Help;
        S2 += Help*Help;
    }
    S /= double(stop);
    S2 /= double(stop);
    sig2 = S2 - S*S;

    SimCheck = false;
    HamilChanged = true;
    Invalidate();
}
//---------------------------------------------------------------------------





void __fastcall TMainForm::Save1Click(TObject *Sender)
{
   ofstream *fitfile;
    if (FitSaveDialog->Execute())
    {
        fitfile = new ofstream(FitSaveDialog->FileName.c_str());
        if (!fitfile)
        {
                ErrorBox->Caption = "File error";
//                ErrorBox->Label->Caption = "File could not be created";
                ErrorBox->ShowModal();
                return;
        }
        *fitfile << "DataFile HansFormat " << endl;
        *fitfile << "Fit of " << FileLabel->Caption.c_str() << " column " << ColumnComboBox->Text.c_str() << endl;
        if (LineFitMode == 0)
        {

            *fitfile << "Hamiltonian Parameters: " << endl;
            HPar.Write(fitfile);
            *fitfile << "Fit and Population Parameters: " << endl;
            EPRSimulationDialog->Write(fitfile);
//            *fitfile << endl <<"ChiSqr " << CalcCycle();
            *fitfile << endl;
        }
          else
          {
            SpectralLineStartDialog->Write(fitfile);
          }
        *fitfile << "a " << aLabel->Caption.c_str();
        *fitfile << "    b " << bLabel->Caption.c_str() << endl;
        *fitfile << ">BEGIN" << endl;
        for (int i=0; i<SimData->Getn() ; i++)
              *fitfile << SimData->Get(i) << endl;
        fitfile->close();

    }
}
//---------------------------------------------------------------------------

int TMainForm::ResetArrays(int MaxNpts, int NewNTrans)
{
     if (Field != NULL)
     {
        for (int i=0; i<nTrans; i++)
            delete[] Field[i];
        delete[] Field;
     }

     if (Ampl != NULL)
     {
        for (int i=0; i<nTrans; i++)
            delete[] Ampl[i];
        delete[] Ampl;
     }

     if (WW != NULL) delete[] WW;

     nTrans = NewNTrans;

     Field = new double*[nTrans];
     Ampl = new double*[nTrans];
     for (int i=0; i<nTrans; i++)
     {
        Field[i] = new double[MaxNpts];
        Ampl[i] = new double[MaxNpts];
     }
     WW = new double[MaxNpts];
}

double TMainForm::CalcCycle()
{
	int count=0;
	double theta = 90.0;
	double phi = 0;
    double SumPop = 0.0;
    double rad = asin(1.0)/90.0;

    char dummy2[40];

	double maxphi = 180;       // for a diagonal g-matrix
    double anglestep;
    if (!ValidReal(EPRSimulationDialog->AngleStepEdit->Text, &anglestep)) return -1.0;
    double freq;
    if (!ValidReal(EPRSimulationDialog->FreqEdit->Text, &freq)) return -1.0;
    double kT;
    if (!ValidReal(EPRSimulationDialog->TempEdit->Text, &kT)) return -1.0;
    kT *= (138/6.62);  //  (k/h * 10-9)
    double modul;
    if (!ValidReal(EPRSimulationDialog->ModulationEdit->Text, &modul)) return -1.0;
    double wx;
    if (!ValidReal(EPRSimulationDialog->LwxEdit->Text, &wx)) return -1.0;
    double wy;
    if (!ValidReal(EPRSimulationDialog->LwyEdit->Text, &wy)) return -1.0;
    double wz;
    if (!ValidReal(EPRSimulationDialog->LwzEdit->Text, &wz)) return -1.0;
    double Px;
    if (!ValidReal(EPRSimulationDialog->PxEdit->Text, &Px)) return -1.0;
    double Py;
    if (!ValidReal(EPRSimulationDialog->PyEdit->Text, &Py)) return -1.0;
    double Pz;
    if (!ValidReal(EPRSimulationDialog->PzEdit->Text, &Pz)) return -1.0;
    if (EPRSimulationDialog->BoltzmannCheckBox->Checked) PopMode = 0;
    if (EPRSimulationDialog->SOISCcheckbox->Checked) PopMode = 1;
    if (EPRSimulationDialog->RPISCcheckbox->Checked) PopMode = 2;
    double XZratio;
    if (!ValidReal(EPRSimulationDialog->XZratioEdit->Text, &XZratio)) return -1.0;
    double YZratio;
    if (!ValidReal(EPRSimulationDialog->YZratioEdit->Text, &YZratio)) return -1.0;

    double tt, kx, ky, kz, SLRx, SLRy, SLRz;
    if (!ValidReal(EPRSimulationDialog->TimeEdit->Text, &tt)) return -1.0;
    if (!ValidReal(EPRSimulationDialog->kxEdit->Text, &kx)) return -1.0;
    if (!ValidReal(EPRSimulationDialog->kyEdit->Text, &ky)) return -1.0;
    if (!ValidReal(EPRSimulationDialog->kzEdit->Text, &kz)) return -1.0;
    if (!ValidReal(EPRSimulationDialog->SLRxEdit->Text, &SLRx)) return -1.0;
    if (!ValidReal(EPRSimulationDialog->SLRyEdit->Text, &SLRy)) return -1.0;
    if (!ValidReal(EPRSimulationDialog->SLRzEdit->Text, &SLRz)) return -1.0;


    double first;
    if (!ValidReal(EPRSimulationDialog->StartFieldEdit->Text, &first)) return -1.0;
    double last;
    if (!ValidReal(EPRSimulationDialog->EndFieldEdit->Text, &last)) return -1.0;
    int npts;
    if (!ValidInt(EPRSimulationDialog->NptsEdit->Text, &npts)) return -1.0;
    double OrThet = 0.0;
    double OrPhi  = 0.0;
    double OrEnergy = 0.0;
    if (EPRSimulationDialog->OrientedRadioButton->Checked)
    {
        if (!ValidReal(EPRSimulationDialog->PartOrThetaEdit->Text, &OrThet)) return -1.0;
        if (!ValidReal(EPRSimulationDialog->PartOrPhiEdit->Text, &OrPhi)) return -1.0;
        if (!ValidReal(EPRSimulationDialog->OrderParEdit->Text, &OrEnergy)) return -1.0;

    }


    int CalcOrder = EPRSimulationDialog->CalcOrderComboBox->ItemIndex+1;

    int MaxNpts = ceil(8/(anglestep*anglestep*rad*rad));  //full half sphere

    if (fabs(HPar.Getg().get(0,1)) > 1e-8)
    {
        maxphi = 180; // if g-matrix non-diagonal in xy plane
    }
    if ((fabs(HPar.Getg().get(1,2)) > 1e-8) || (fabs(HPar.Getg().get(0,2)) > 1e-8) )
    {
         maxphi = 360;
    }

    // Now we'll need some kind of 'reasonable maximum'
    // for anglestep = 1 (0.01745 degree) we'd have 26000 points or so
    // with 64 byte real, this gets to 200 kb. I guess an anglestep of 0.25 degree
    // is the absolutely smallest value we can take. It would be of the order of
    // 4 Mb for each array.

    if (NumberChanged)
    {

    // We are neglecting all forbidden transitions and take only
    // the EPR transitions for each nuclear spin state
        int nlevels = (HPar.GetMultS() * HPar.GetMultI(0));

        if  (HamiltonianDialog->PerturButton->Checked == true)
            ResetArrays(MaxNpts,(HPar.GetMultS()-1) * HPar.GetMultI(0));
        else
            ResetArrays(MaxNpts,nlevels*(nlevels-1)/2);

        HamilChanged = true;
        WidthChanged = true;
        PopChanged = true;
        DisPersChanged = true;
        NumberChanged = false;
    }

//  if the Hamiltonian changed...  but we don't loop over this
    Spin S(HPar.GetMultS());
    Spin I(HPar.GetMultI(0));
    if  (HamiltonianDialog->PerturButton->Checked == true)
        I.SetMult(1);

    SpinHam H(HPar.GetMultS(),HPar.GetMultI(0));
    if (HamiltonianDialog->PerturButton->Checked )
	        H.Reset(HPar.GetMultS());

            //   Set g-tensor

    H.set_g_tensor(HPar.Getg().get(0,0),HPar.Getg().get(1,1),HPar.Getg().get(2,2),
         HPar.Getg().get(0,1),HPar.Getg().get(0,2),HPar.Getg().get(1,2));
    if (HPar.GetnI() >=1)
            H.set_A_tensor(HPar.GetA(0).get(0,0),HPar.GetA(0).get(1,1),HPar.GetA(0).get(2,2),
           HPar.GetA(0).get(0,1),HPar.GetA(0).get(0,2),HPar.GetA(0).get(1,2));

            // Set the crystal field / zero-field parameters

    for (int i=0; i<10; i++)
    {
        H.SetCF(i,HPar.GetCF(i));
    }


    H.setH0();
//  So far we have a very general case !

//	while (!kbhit()) {}


//  freq  = 240;

    double l,m,n;
    double field1, field2, intens1, intens2;

    double B0 = (first+last)/2.0;

    Spectrum Sum(npts, first, last);
    Spectrum DisPers(npts, first, last);
    Spectrum SumDeriv(npts, first, last);
    Spectrum DisPersDeriv(npts, first, last);


    Vector E(H.GetOrder());
    Vector E1(E.order());
    Vector Tran(E.order()*(E.order()-1)/2);
    Vector Prob(E.order()*(E.order()-1)/2);
    Vector Boltzmann(E.order()*(E.order()-1)/2);
    double width;

//    double *theta;
//    double *phi;
//    double *


//	cout << "estimation npts : " << PI/(2.0*anglestep*anglestep*rad*rad);

	count = 0;
    phi = maxphi;

    double phistep = anglestep;
    double B1;
    int ccol = ColumnComboBox->ItemIndex;
    double hyp;
    double ax, ay, az;
    double g2, Eshift;

    if ((NumberChanged) || (WidthChanged) || (HamilChanged) || (PopChanged))
    {
	while (theta > -0.00001)        //theta integrated from 90 to 0
	{
        if (theta > 0.01) phistep = anglestep / sin(theta*rad);
            else phistep = maxphi+1;
        if (HamiltonianDialog->AxialCheckBox->Checked == true)
            phistep = maxphi+1; //force a single phi value
		while (phi > 0.0)           //Sufficient if A has same principal axes as g
		{
            l = sin(theta*rad)*cos(phi*rad);
            m = sin(theta*rad)*sin(phi*rad);
            n = cos(theta*rad);

      		if ((HamilChanged) || ((PopChanged) && (HPar.GetnI()>0)))
            {
                H.set_field(B0, theta, phi);
       		    H.addZeeman2(S,I);

//			g_eff = (g*R).length(); //Calculate effective g-value
                if (HPar.GetnI() == 0)
                {
                    if (H.GetEmult() < 4)
                    {
       		            E = H.eigenval();
                        if (H.GetEmult() == 3) nTrans=2;
                        for (int i =0; i< H.GetEmult()-1; i++)
                        {
                            B1 = B0-((E[i+1]-E[i])-freq)*H.radical_resonance(freq)/freq;
                            for (int k=1; k< CalcOrder; k++)
                            {
                                H.set_field(B1, theta, phi);
               		            H.addZeeman2(S,I);
                                E1 = H.eigenval();
                                B1 = B1 -((E1[i+1]-E1[i])-freq)*H.radical_resonance(freq)/freq;
                            }
                            Field[i][count] = B1;
//                       sprintf(dummy2,"%lf,  %lf ", theta, B1);
//                        OutputMemo->Lines->Add(dummy2);
                        }
                    }
                      else // S>1 We need correct transition probabilities
                      {
                        E = H.eigenvec();
                        H.transitions();
                        Tran = H.get_trans();
                        Prob = H.get_tprob();
                        PopChanged = false;    //we'll do the population thing here

                        int jj = 0;

                        SumPop=0.0;
                        for (int i =0; i< E.order(); i++)
                                 SumPop += exp(-E[i]/kT);

                        for (int i =0; i< E.order()-1; i++)
                                    for (int j=i+1; j<E.order();j++)
                            Boltzmann.set(jj++,(exp(-E[i]/kT) - exp(-E[j]/kT))/SumPop);

                        int ii = 0;

                        for (int i=0; i<Tran.order(); i++)
                        {
                            Field[i][count] = B0;
                            Ampl[i][count] = 0.0;
                            if ((Tran[i]>0.6*freq) && (Tran[i]<1.4*freq))
                            {
                                B1 = B0-(Tran[i]-freq)*H.radical_resonance(freq)/freq;
                                for (int k=1; k< CalcOrder; k++)
                                {
                                    H.set_field(B1, theta, phi);
               		                H.addZeeman2(S,I);

                                    E = H.eigenvec();
                                    H.transitions();
                                    Tran = H.get_trans();
                                    Prob = H.get_tprob();

                                    B1 = B1 -(Tran[i]-freq)*H.radical_resonance(freq)/freq;
                                 }

                                 Ampl[i][count] = Prob[i]*Boltzmann[i];
                                 if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                        Ampl[i][count] *= sin(theta*rad);
                                 Field[i][count] = B1;

                            }

                                       //end if

                        }                 //end for i=
                      }  // end else
                }
                  else
                  if (HamiltonianDialog->PerturButton->Checked == true)
                  {
                    E = H.eigenval();
                    ax = l * HPar.Getg().get(0,0) * HPar.GetA(0).get(0,0)
                            + m * HPar.Getg().get(1,1) * HPar.GetA(0).get(0,1)
                            + n * HPar.Getg().get(2,2) * HPar.GetA(0).get(0,2);
                    ay = l * HPar.Getg().get(0,0) * HPar.GetA(0).get(1,0)
                            + m * HPar.Getg().get(1,1) * HPar.GetA(0).get(1,1)
                            + n * HPar.Getg().get(2,2) * HPar.GetA(0).get(1,2);
                    az = l * HPar.Getg().get(0,0) * HPar.GetA(0).get(2,0)
                            + m * HPar.Getg().get(1,1) * HPar.GetA(0).get(2,1)
                            + n * HPar.Getg().get(2,2) * HPar.GetA(0).get(2,2);
                    g2 = l*l * HPar.Getg().get(0,0) * HPar.Getg().get(0,0)
                         + m*m * HPar.Getg().get(1,1) * HPar.Getg().get(1,1)
                         + n*n * HPar.Getg().get(2,2) * HPar.Getg().get(2,2);

                    hyp = sqrt((ax*ax+ay*ay+az*az)/g2);

                    for (int i =0; i< H.GetEmult()-1; i++)
                       for (int nuc=0; nuc < HPar.GetMultI(0); nuc++)
                       {
                        Eshift = double((HPar.GetMultI(0)-1) - 2*nuc)* hyp /2000.0;
                        B1 = B0-((E[i+1]-E[i]+Eshift)-freq)*H.radical_resonance(freq)/freq;
                        for (int k=1; k< CalcOrder; k++)
                        {
                            H.set_field(B1, theta, phi);
               		        H.addZeeman2(S,I);
                            E1 = H.eigenval();
                            B1 = B1 -((E1[i+1]-E1[i]+Eshift)-freq)*H.radical_resonance(freq)/freq;
                        }
                        Field[HPar.GetMultI(0)*i+nuc][count] = B1;
                       }

                  }
                  else  //complete diagonalization
                  {
                    E = H.eigenvec();
                    H.transitions();
                    Tran = H.get_trans();
                    Prob = H.get_tprob();
//                    PopChanged = false;    //we'll do the population thing here

                    int jj = 0;

                    SumPop=0.0;
                    for (int i =0; i< E.order(); i++)
                             SumPop += exp(-E[i]/kT);

                    for (int i =0; i< E.order()-1; i++)
                                for (int j=i+1; j<E.order();j++)
                        Boltzmann.set(jj++,(exp(-E[i]/kT) - exp(-E[j]/kT))/SumPop);

                    int ii = 0;

                    for (int i=0; i<Tran.order(); i++)
                    {
                        if ((Tran[i]>0.6*freq) && (Tran[i]<1.4*freq))
                        {
                            B1 = B0-(Tran[i]-freq)*H.radical_resonance(freq)/freq;
                            for (int k=1; k< CalcOrder; k++)
                            {
                                H.set_field(B1, theta, phi);
               		            H.addZeeman2(S,I);

                                E = H.eigenvec();
                                H.transitions();
                                Tran = H.get_trans();
                                Prob = H.get_tprob();

                                B1 = B1 -(Tran[i]-freq)*H.radical_resonance(freq)/freq;
                            }
                            if (ii<nTrans)
                            {
                                Ampl[ii][count] = Prob[i]*Boltzmann[i];
                                if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[ii][count] *= sin(theta*rad);
                                Field[ii][count] = B1;
                            }
                            ii++;
                        }             //end if

                    }                 //end for i=
                }                 //end else nI

            }                     //end if HamilChanged

//      For intersystem crossing with SOC we have the following intensities
//            intens1 = (1-3*l*l)*Px + (1-3*m*m)*Py + (1-3*n*n)*Pz;
//            intens2 = -1.0 * intens1;
//
//  However, if we have to take Spin lattice relaxation into account,
//  the results will not be symmetrical, as the spectrum will evolve into
//  a strictly absorptive spectrum for a Boltzmann distribution
//  It is possible to calculate the time evolution.
//  For that we have to use some kind of differential equation and
//  treat the different levels individually.
//  Upper and lower level population at t=0 can be given by
//  P0(0) = P2(0) = l*l(py+Pz)/2 + m*m(Px+Pz)/2 + n*n(Px+Py)/2
//      center level population can be given by
//  Pl(0) = l*l*Px + m*m*Py + n*n*Pz
//
//  The triplet decay component can also be described very equivalently
//  dP0/dt = dP2/dt = -(l*l(ky+kz)/2 + m*m(kx+kz)/2 + n*n(kx+ky)/2)*P0 (or P2)
//  and
//  dP1/dt = -P1 * (l*l*kx + m*m*ky + n*n*kz)
//
//  The spin lattice relaxation might also be anisotropic
//  Let's just assume that we have relaxation between adjacent levels only
//  and let's say kSL = l*l*kSLx + m*m*kSLy + n*n*kSLz
//  then we have
//  dP0/dt = -kSL exp(-hv/kT) P0 + kSL P1
//  dP1/dt = -kSL P1 + kSL exp(-hv/kT) P0 - kSL exp(-hv/kT) P1 +kSL P2
//  dP2/dt = -kSL P2 + kSL exp(-hv/kT) P1
//
//  One method to use could be Runge Kutta
            double pops[3];
            if (PopChanged)
            {
                switch (PopMode)
                {
                    case 1:    // ISC through SOC   Only for triplets !!
                    if (tt>0.0001)
                      CalcDecayedPops(kx,ky,kz,SLRx,SLRy,SLRz,tt,freq,kT,Px,Py,Pz,pops,l,m,n);

                    for (int i =0; i< H.GetEmult()-1; i++)
                    {
                        if (tt<=0.0001)
                        {
                            Ampl[i][count] = (1-3*l*l)*Px + (1-3*m*m)*Py + (1-3*n*n)*Pz;
                            if (i == 1) Ampl[i][count] *= -1;
                            if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                        }
                          else
                          {
                            Ampl[i][count] = pops[i+1]-pops[i];
                            if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                          }

                    }

                    break;

                    case 2:
                    for (int i =0; i< H.GetEmult()-1; i++)
                    {
                        Ampl[i][count] = l*l*XZratio + m*m*YZratio + n*n;
                        if (i == 1) Ampl[i][count] *= -1;
                        if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                    }
                    break;

                    default:
                    {
                        SumPop = 0.0;
                        for (int i =0; i< H.GetEmult(); i++)
                            SumPop += exp(-E[i]/kT);
                        for (int i =0; i< H.GetEmult()-1; i++)
                        {
                            Ampl[i][count] = (exp(-E[i]/kT) - exp(-E[i+1]/kT))/SumPop;
                            if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                        }
                    }
                    break;
                }
                if (fabs(theta-90.0) < 0.00001)   // take only half intensity
                        for (int i =0; i< H.GetEmult()-1; i++) Ampl[i][count] *= 0.5;    //  in XY plane
                if (phi<phistep)
                        for (int i =0; i< H.GetEmult()-1; i++)
                            Ampl[i][count] *= phi/phistep;
                double CosineOrient;
                if (EPRSimulationDialog->OrientedRadioButton->Checked)
                {
                    CosineOrient = l*sin(OrThet*rad)*cos(OrPhi*rad)
                                + m*sin(OrThet*rad)*sin(OrPhi*rad)
                                + n*cos(OrThet*rad);
                    // We'll take a very simple distribution function
                    //    for which the energy is proportional to -cos*cos
                        for (int i =0; i< H.GetEmult()-1; i++)
                            Ampl[i][count] *= exp(CosineOrient*CosineOrient*OrEnergy);
                            // we don;t worry about normalization here
                }
                if ((HPar.GetnI()>0) && (HamiltonianDialog->PerturButton->Checked == true))
                {
                    for (int i = 0; i< H.GetEmult()-1; i++)
                    {
                        Ampl[i*HPar.GetMultI(0)][count] = Ampl[i][count];
                        for (int jj=1;jj<HPar.GetMultI(0);jj++)
                            Ampl[i*HPar.GetMultI(0)+jj][count] = Ampl[i*HPar.GetMultI(0)][count];
                    }
                }


            }

// Now the population intensity
//    let's say we have kx, ky, and kz. Put them in Axx, Ayy, Azz
//    The center line population will then be according to its
//    the amount of the x,y, z spin basis functions:
//      thus l*l*kx + m*m*ky + n*n*kz
//   The outer levels will be a mixture of the rest:
//       (1-l*l)kx+(1-m*m)ky+(1-n*n)kz / 2

//			g_eff = (g*R).length(); //Calculate effective g-value

            if (WidthChanged)
                WW[count] = sqrt( wx*wx*l*l + wy*wy*m*m + wz*wz*n*n) * 0.001;

			phi -= phistep;
//            PhiLabel->Caption = phi;
//            PhiLabel->Update();
            count++;

		}
//		Sum.addLorentzArray(Center,Width, Intensity, count);
//      DisPers.addLorentzDispersArray(Center, Width, Intensity, count);

		phi = maxphi;
		theta -= anglestep;
        ThetaLabel->Caption = theta;
        ThetaLabel->Update();
	}

    if (EPRSimulationDialog->LorentzCheckBox->Checked)
        for (int i=0; i<nTrans; i++)
        {
            Sum.addLorentzArray(Field[i],WW, Ampl[i], count);
            DisPers.addLorentzDispersArray(Field[i], WW, Ampl[i], count);
        }
      else
        if (EPRSimulationDialog->GaussCheckBox->Checked)
        {
            DataArray Temporary(npts,3);
            DataPoint Ptemp(3);
            for (int i=0; i<nTrans; i++)
            {
                     Sum.addGaussArray(Field[i],WW, Ampl[i], count);
            }
            for (int i = 0; i<npts; i++)
            {
                Ptemp.Set(0,Sum.getField(i));
                Ptemp.Set(1,Sum.getIntens(i));
                Temporary.Add(Ptemp);
            }

            Temporary.HilbertTransForm(0,1,2,0.0);

            for (int i = 0; i<npts; i++)
            {
                Ptemp = Temporary.Get(i);
                DisPers.setField(i,Ptemp.Get(0));
                DisPers.setIntens(i, -Ptemp.Get(2));
            }
    }

    DataArray TempData = ExpData->ReSample(first,last,npts);



    if (SimData != NULL) { delete SimData; SimData = NULL;}
    if (AbsDispData != NULL) { delete AbsDispData; AbsDispData = NULL;}

    SimData = new DataArray(Sum.getNpts(),4);
    AbsDispData = new DataArray(Sum.getNpts(),2);
    DataPoint Psim(4);
    DataPoint AbsDispPoint(2);

    if (DerivMode == 1)
    {
        SumDeriv.Derivative(Sum, modul);
        DisPersDeriv.Derivative(DisPers, modul);
        for (int i=0; i< Sum.getNpts(); i++)
        {
            Psim.Set(0,Sum.getField(i));
            AbsDispPoint.Set(0,Sum.getField(i));
                if (ccol != 0) Psim.Set(1,TempData.Get(i).Get(ccol));
                   else Psim.Set(1,0.0);
            Psim.Set(2,SumDeriv.getIntens(i));
            Psim.Set(3,0.0);
            Psim.Set(4,DisPersDeriv.getIntens(i));
            AbsDispPoint.Set(1,SumDeriv.getIntens(i));
            AbsDispPoint.Set(2,DisPersDeriv.getIntens(i));
            SimData->Add(Psim);
            AbsDispData->Add(AbsDispPoint);
        }
    }
        else
            for (int i=0; i< Sum.getNpts(); i++)
            {
                Psim.Set(0,Sum.getField(i));
                AbsDispPoint.Set(0,Sum.getField(i));
                if (ccol != 0) Psim.Set(1,TempData.Get(i).Get(ccol));
                   else Psim.Set(1,0.0);
                Psim.Set(2,Sum.getIntens(i));
                Psim.Set(3,0.0);
                Psim.Set(4,DisPers.getIntens(i));
                AbsDispPoint.Set(1,Sum.getIntens(i));
                AbsDispPoint.Set(2,DisPers.getIntens(i));
                SimData->Add(Psim);
                AbsDispData->Add(AbsDispPoint);
            }

    SimCheck = true;
    }
    double a,b;
    DataPoint Psim2(2);
    double DispAngle=0;
    if (ValidReal(EPRSimulationDialog->DispersEdit->Text, &DispAngle))
    {
        if (fabs(DispAngle) >0.001)
        {
            DispAngle *= (PI/180.0);
            for (int i=0; i<SimData->Getn(); i++)
            {
                Psim2 = SimData->Get(i);
//                AbsDispPoint = AbsDispData->Get(i);
                Psim2.Set(2, cos(DispAngle) * AbsDispData->Get(i).Get(1)
                               + sin(DispAngle) * AbsDispData->Get(i).Get(2));
                Psim2.Set(4, cos(DispAngle) * AbsDispData->Get(i).Get(2)
                               - sin(DispAngle) * AbsDispData->Get(i).Get(1));
                SimData->Set(i,Psim2);
            }
        }
    }

    double OldB=0.0;
    if (!ValidReal(bLabel->Caption, &OldB))
        EPRSimulationDialog->FixedbCheckBox->Checked = false;

    if (ccol != 0)
        SimData->LinearRegres(2,1,&a,&b);
      else {
        a = 0.0;
        b =  1.0/SimData->Maxima().Get(2);
      }
    if (EPRSimulationDialog->FixedbCheckBox->Checked)
        b=OldB;

    aLabel->Caption = a;
    bLabel->Caption = b;
    double chisqr=0.0;
    double OldValue;
    for (int i=0; i<SimData->Getn(); i++)
    {
        Psim2 = SimData->Get(i);
        OldValue = Psim2.Get(2);
        Psim2.Set(2, a + b * OldValue);
        Psim2.Set(3, Psim2.Get(1) - Psim2.Get(2));
        OldValue = Psim2.Get(4);
        Psim2.Set(4, b * OldValue);
        SimData->Set(i,Psim2);
        chisqr += Psim2.Get(3)*Psim2.Get(3);
    }

    chisqr /= SimData->Getn();
    if (sig2 > 0.0) chisqr /= sig2;
    return chisqr;

}

double TMainForm::HighSpinCalc()
{
	int count=0;
	double theta = 90.0;
	double phi = 0;
    double SumPop = 0.0;
    char dummy2[40];

	double maxphi = 180;       // for a diagonal g-matrix
    double anglestep;
    if (!ValidReal(EPRSimulationDialog->AngleStepEdit->Text, &anglestep)) return -1.0;
    double freq;
    if (!ValidReal(EPRSimulationDialog->FreqEdit->Text, &freq)) return -1.0;
    double kT;
    if (!ValidReal(EPRSimulationDialog->TempEdit->Text, &kT)) return -1.0;
    kT *= (138/6.62);  //  (k/h * 10-9)
    double modul;
    if (!ValidReal(EPRSimulationDialog->ModulationEdit->Text, &modul)) return -1.0;
    double wx;
    if (!ValidReal(EPRSimulationDialog->LwxEdit->Text, &wx)) return -1.0;
    double wy;
    if (!ValidReal(EPRSimulationDialog->LwyEdit->Text, &wy)) return -1.0;
    double wz;
    if (!ValidReal(EPRSimulationDialog->LwzEdit->Text, &wz)) return -1.0;
    double Px;
    if (!ValidReal(EPRSimulationDialog->PxEdit->Text, &Px)) return -1.0;
    double Py;
    if (!ValidReal(EPRSimulationDialog->PyEdit->Text, &Py)) return -1.0;
    double Pz;
    if (!ValidReal(EPRSimulationDialog->PzEdit->Text, &Pz)) return -1.0;
    double WStrain;
    if (!ValidReal(EPRSimulationDialog->StrainEdit->Text, &WStrain)) return -1.0;

    if (EPRSimulationDialog->BoltzmannCheckBox->Checked) PopMode = 0;
    if (EPRSimulationDialog->SOISCcheckbox->Checked) PopMode = 1;
    if (EPRSimulationDialog->RPISCcheckbox->Checked) PopMode = 2;
    double XZratio;
    if (!ValidReal(EPRSimulationDialog->XZratioEdit->Text, &XZratio)) return -1.0;
    double YZratio;
    if (!ValidReal(EPRSimulationDialog->YZratioEdit->Text, &YZratio)) return -1.0;
    double first;
    if (!ValidReal(EPRSimulationDialog->StartFieldEdit->Text, &first)) return -1.0;
    double last;
    if (!ValidReal(EPRSimulationDialog->EndFieldEdit->Text, &last)) return -1.0;
    double thet1;
    if (!ValidReal(EPRSimulationDialog->ThetaEdit->Text, &thet1)) return -1.0;
    double phi1;
    if (!ValidReal(EPRSimulationDialog->PhiEdit->Text, &phi1)) return -1.0;
    int npts;
    if (!ValidInt(EPRSimulationDialog->NptsEdit->Text, &npts)) return -1.0;
    int CalcOrder = EPRSimulationDialog->CalcOrderComboBox->ItemIndex+1;


    double rad = asin(1.0)/90.0;

    int MaxNpts = ceil(8/(anglestep*anglestep*rad*rad));  //full half sphere

    if (fabs(HPar.Getg().get(0,1)) > 1e-8)
    {
        maxphi = 180; // if g-matrix non-diagonal in xy plane
    }
    if ((fabs(HPar.Getg().get(1,2)) > 1e-8) || (fabs(HPar.Getg().get(0,2)) > 1e-8) )
    {
         maxphi = 360;
    }

    // Now we'll need some kind of 'reasonable maximum'
    // for anglestep = 1 (0.01745 degree) we'd have 26000 points or so
    // with 64 byte real, this gets to 200 kb. I guess an anglestep of 0.25 degree
    // is the absolutely smallest value we can take. It would be of the order of
    // 4 Mb for each array.

    if (NumberChanged)
    {

    // We are neglecting all forbidden transitions and take only
    // the EPR transitions for each nuclear spin state
        int nlevels = (HPar.GetMultS() * HPar.GetMultI(0));

        if  (HamiltonianDialog->PerturButton->Checked == true)
            ResetArrays(MaxNpts,(HPar.GetMultS()-1) * HPar.GetMultI(0));
        else
            ResetArrays(MaxNpts,nlevels*(nlevels-1)/2);

        HamilChanged = true;
        WidthChanged = true;
        PopChanged = true;
        DisPersChanged = true;
        NumberChanged = false;
    }

//  if the Hamiltonian changed...  but we don't loop over this
    Spin S(HPar.GetMultS());
    Spin I(HPar.GetMultI(0));
    if  (HamiltonianDialog->PerturButton->Checked == true)
        I.SetMult(1);

    SpinHam H(HPar.GetMultS(),HPar.GetMultI(0));
    if (HamiltonianDialog->PerturButton->Checked )
	        H.Reset(HPar.GetMultS());

            //   Set g-tensor

    H.set_g_tensor(HPar.Getg().get(0,0),HPar.Getg().get(1,1),HPar.Getg().get(2,2),
         HPar.Getg().get(0,1),HPar.Getg().get(0,2),HPar.Getg().get(1,2));
    if (HPar.GetnI() >=1)
            H.set_A_tensor(HPar.GetA(0).get(0,0),HPar.GetA(0).get(1,1),HPar.GetA(0).get(2,2),
           HPar.GetA(0).get(0,1),HPar.GetA(0).get(0,2),HPar.GetA(0).get(1,2));

            // Set the crystal field / zero-field parameters

    for (int i=0; i<10; i++)
    {
        H.SetCF(i,HPar.GetCF(i));
    }


    H.setH0();
//  So far we have a very general case !

//	while (!kbhit()) {}


//  freq  = 240;

    double l,m,n;
    double field1, field2, intens1, intens2;

    double B0 = (first+last)/2.0;

    Spectrum Sum(npts, first, last);
    Spectrum DisPers(npts, first, last);
    Spectrum SumDeriv(npts, first, last);
    Spectrum DisPersDeriv(npts, first, last);


    Vector E(H.GetOrder());
    Vector E1(E.order());
    Vector Tran(E.order()*(E.order()-1)/2);
    Vector Prob(E.order()*(E.order()-1)/2);
    Vector Boltzmann(E.order()*(E.order()-1)/2);
    double width;

//    double *theta;
//    double *phi;
//    double *


//	cout << "estimation npts : " << PI/(2.0*anglestep*anglestep*rad*rad);

	count = 0;
    phi = maxphi;

    double phistep = anglestep;
    double B1;
    int ccol = ColumnComboBox->ItemIndex;
    double hyp;
    double ax, ay, az;
    double g2, Eshift;

//    if ((NumberChanged) || (WidthChanged) || (HamilChanged) || (PopChanged))
//    {
	while (theta > -0.00001)        //theta integrated from 90 to 0
	{
        if (theta > 0.01) phistep = anglestep / sin(theta*rad);
            else phistep = maxphi+1;
        if (HamiltonianDialog->AxialCheckBox->Checked == true)
            phistep = maxphi+1; //force a single phi value
		while (phi > 0.0)           //Sufficient if A has same principal axes as g
		{
            if (EPRSimulationDialog->CrystalRadioButton->Checked)
            {
                theta = thet1;
                phi = phi1;
            }
            l = sin(theta*rad)*cos(phi*rad);
            m = sin(theta*rad)*sin(phi*rad);
            n = cos(theta*rad);

            double wid0 = sqrt( wx*wx*l*l + wy*wy*m*m + wz*wz*n*n) * 0.001; //(in MHz now)
                 // Note these are now GHz ...
//
//          if ((HamilChanged) || ((PopChanged) && (HPar.GetnI()>0)))
//          {
            for (int fieldi=0; fieldi<npts;fieldi++)
            {
                B0 = first + fieldi*(last-first)/(npts-1);
                H.set_field(B0, theta, phi);
       		    H.addZeeman2(S,I);

                width = wid0+ WStrain*0.01 * fabs(B0-8.56);
//			g_eff = (g*R).length(); //Calculate effective g-value
                if (HPar.GetnI() == 0)
                {
                    E = H.eigenvec();
                    H.transitions(1);
                    Tran = H.get_trans();
                    Prob = H.get_tprob();
//                        PopChanged = false;    //we'll do the population thing here

                    int jj = 0;
                    SumPop=0.0;
                    for (int i =0; i< E.order(); i++)
                              SumPop += exp(-E[i]/kT);

                    for (int i =0; i< E.order()-1; i++)
                                    for (int j=i+1; j<E.order();j++)
                            Boltzmann.set(jj++,(exp(-E[i]/kT) - exp(-E[j]/kT))/SumPop);

                    int ii = 0;

                    for (int i=0; i<Tran.order(); i++)
                    {

                        intens1 = width/(width*width+4.0*(Tran[i]-freq)*(Tran[i]-freq));
                        intens2 = 2*(Tran[i]-freq)/(width*width+4.0*(Tran[i]-freq)*(Tran[i]-freq));

                        intens1 *= (Prob[i]*Boltzmann[i]);
                        intens2 *= (Prob[i]*Boltzmann[i]);

                        Sum.addIntens(fieldi,intens1);
                        DisPers.addIntens(fieldi, intens2);
                    }
                }
                  else {}
            }
/* Not implemented for Nuclear Spin

                  if (HamiltonianDialog->PerturButton->Checked == true)
                  {
                    E = H.eigenval();
                    ax = l * HPar.Getg().get(0,0) * HPar.GetA(0).get(0,0)
                            + m * HPar.Getg().get(1,1) * HPar.GetA(0).get(0,1)
                            + n * HPar.Getg().get(2,2) * HPar.GetA(0).get(0,2);
                    ay = l * HPar.Getg().get(0,0) * HPar.GetA(0).get(1,0)
                            + m * HPar.Getg().get(1,1) * HPar.GetA(0).get(1,1)
                            + n * HPar.Getg().get(2,2) * HPar.GetA(0).get(1,2);
                    az = l * HPar.Getg().get(0,0) * HPar.GetA(0).get(2,0)
                            + m * HPar.Getg().get(1,1) * HPar.GetA(0).get(2,1)
                            + n * HPar.Getg().get(2,2) * HPar.GetA(0).get(2,2);
                    g2 = l*l * HPar.Getg().get(0,0) * HPar.Getg().get(0,0)
                         + m*m * HPar.Getg().get(1,1) * HPar.Getg().get(1,1)
                         + n*n * HPar.Getg().get(2,2) * HPar.Getg().get(2,2);

                    hyp = sqrt((ax*ax+ay*ay+az*az)/g2);

                    for (int i =0; i< H.GetEmult()-1; i++)
                       for (int nuc=0; nuc < HPar.GetMultI(0); nuc++)
                       {
                        Eshift = double((HPar.GetMultI(0)-1) - 2*nuc)* hyp /2000.0;
                        B1 = B0-((E[i+1]-E[i]+Eshift)-freq)*H.radical_resonance(freq)/freq;
                        for (int k=1; k< CalcOrder; k++)
                        {
                            H.set_field(B1, theta, phi);
               		        H.addZeeman2(S,I);
                            E1 = H.eigenval();
                            B1 = B1 -((E1[i+1]-E1[i]+Eshift)-freq)*H.radical_resonance(freq)/freq;
                        }
                        Field[HPar.GetMultI(0)*i+nuc][count] = B1;
                       }

                  }
                  else  //complete diagonalization
                  {
                    E = H.eigenvec();
                    H.transitions();
                    Tran = H.get_trans();
                    Prob = H.get_tprob();
//                    PopChanged = false;    //we'll do the population thing here

                    int jj = 0;

                    SumPop=0.0;
                    for (int i =0; i< E.order(); i++)
                             SumPop += exp(-E[i]/kT);

                    for (int i =0; i< E.order()-1; i++)
                                for (int j=i+1; j<E.order();j++)
                        Boltzmann.set(jj++,(exp(-E[i]/kT) - exp(-E[j]/kT))/SumPop);

                    int ii = 0;

                    for (int i=0; i<Tran.order(); i++)
                    {
                        if ((Tran[i]>0.6*freq) && (Tran[i]<1.4*freq))
                        {
                            B1 = B0-(Tran[i]-freq)*H.radical_resonance(freq)/freq;
                            for (int k=1; k< CalcOrder; k++)
                            {
                                H.set_field(B1, theta, phi);
               		            H.addZeeman2(S,I);

                                E = H.eigenvec();
                                H.transitions();
                                Tran = H.get_trans();
                                Prob = H.get_tprob();

                                B1 = B1 -(Tran[i]-freq)*H.radical_resonance(freq)/freq;
                            }
                            if (ii<nTrans)
                            {
                                Ampl[ii][count] = Prob[i]*Boltzmann[i];
                                if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[ii][count] *= sin(theta*rad);
                                Field[ii][count] = B1;
                            }
                            ii++;
                        }             //end if

                    }                 //end for i=
                }                 //end else nI

//          }                     //end if HamilChanged


//            intens1 = (1-3*l*l)*Px + (1-3*m*m)*Py + (1-3*n*n)*Pz;
//            intens2 = -1.0 * intens1;
            if (PopChanged)
            {
                switch (PopMode)
                {
                    case 1:
                    for (int i =0; i< H.GetEmult()-1; i++)
                    {
                        Ampl[i][count] = (1-3*l*l)*Px + (1-3*m*m)*Py + (1-3*n*n)*Pz;
                        if (i == 1) Ampl[i][count] *= -1;
                        if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                           }
                    break;

                    case 2:
                    for (int i =0; i< H.GetEmult()-1; i++)
                    {
                        Ampl[i][count] = l*l*XZratio + m*m*YZratio + n*n;
                        if (i == 1) Ampl[i][count] *= -1;
                        if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                    }
                    break;

                    default:
                    {
                        SumPop = 0.0;
                        for (int i =0; i< H.GetEmult(); i++)
                            SumPop += exp(-E[i]/kT);
                        for (int i =0; i< H.GetEmult()-1; i++)
                        {
                            Ampl[i][count] = (exp(-E[i]/kT) - exp(-E[i+1]/kT))/SumPop;
                            if (HamiltonianDialog->AxialCheckBox->Checked == true)
                                    Ampl[i][count] *= sin(theta*rad);
                        }
                    }
                    break;
                }
                if (fabs(theta-90.0) < 0.00001)   // take only half intensity
                        for (int i =0; i< H.GetEmult()-1; i++) Ampl[i][count] *= 0.5;    //  in XY plane
                if (phi<phistep)
                        for (int i =0; i< H.GetEmult()-1; i++)
                            Ampl[i][count] *= phi/phistep;
                if ((HPar.GetnI()>0) && (HamiltonianDialog->PerturButton->Checked == true))
                {
                    for (int i = 0; i< H.GetEmult()-1; i++)
                    {
                        Ampl[i*HPar.GetMultI(0)][count] = Ampl[i][count];
                        for (int jj=1;jj<HPar.GetMultI(0);jj++)
                            Ampl[i*HPar.GetMultI(0)+jj][count] = Ampl[i*HPar.GetMultI(0)][count];
                    }
                }

            }

// Now the population intensity
//    let's say we have kx, ky, and kz. Put them in Axx, Ayy, Azz
//    The center line population will then be according to its
//    the amount of the x,y, z spin basis functions:
//      thus l*l*kx + m*m*ky + n*n*kz
//   The outer levels will be a mixture of the rest:
//       (1-l*l)kx+(1-m*m)ky+(1-n*n)kz / 2

//			g_eff = (g*R).length(); //Calculate effective g-value

            if (WidthChanged)
                WW[count] = sqrt( wx*wx*l*l + wy*wy*m*m + wz*wz*n*n) * 0.001;
*/
			phi -= phistep;
//            PhiLabel->Caption = phi;
//            PhiLabel->Update();
            count++;
            if (EPRSimulationDialog->CrystalRadioButton->Checked) break;
		}
//		Sum.addLorentzArray(Center,Width, Intensity, count);
//      DisPers.addLorentzDispersArray(Center, Width, Intensity, count);

		phi = maxphi;
		theta -= anglestep;
        ThetaLabel->Caption = theta;
        ThetaLabel->Update();
        if (EPRSimulationDialog->CrystalRadioButton->Checked) break;
	}

/*    if (EPRSimulationDialog->LorentzCheckBox->Checked)
        for (int i=0; i<nTrans; i++)
        {
            Sum.addLorentzArray(Field[i],WW, Ampl[i], count);
            DisPers.addLorentzDispersArray(Field[i], WW, Ampl[i], count);
        }
      else
        if (EPRSimulationDialog->GaussCheckBox->Checked)
        {
            DataArray Temporary(npts,3);
            DataPoint Ptemp(3);
            for (int i=0; i<nTrans; i++)
                     Sum.addGaussArray(Field[i],WW, Ampl[i], count);
            for (int i = 0; i<npts; i++)
            {
                Ptemp.Set(0,Sum.getField(i));
                Ptemp.Set(1,Sum.getIntens(i));
                Temporary.Add(Ptemp);
            }

            Temporary.HilbertTransForm(0,1,2,0.0);

            for (int i = 0; i<npts; i++)
            {
                Ptemp = Temporary.Get(i);
                DisPers.setField(i,Ptemp.Get(0));
                DisPers.setIntens(i, -Ptemp.Get(2));
            }
*/
//    }

    DataArray TempData = ExpData->ReSample(first,last,npts);



    if (SimData != NULL) delete SimData;
    if (AbsDispData != NULL) delete AbsDispData;

    SimData = new DataArray(Sum.getNpts(),4);
    AbsDispData = new DataArray(Sum.getNpts(),2);
    DataPoint Psim(4);
    DataPoint AbsDispPoint(2);

    if (DerivMode == 1)
    {
        SumDeriv.Derivative(Sum, modul);
        DisPersDeriv.Derivative(DisPers, modul);
        for (int i=0; i< Sum.getNpts(); i++)
        {
            Psim.Set(0,Sum.getField(i));
            AbsDispPoint.Set(0,Sum.getField(i));
                if (ccol != 0) Psim.Set(1,TempData.Get(i).Get(ccol));
                   else Psim.Set(1,0.0);
            Psim.Set(2,SumDeriv.getIntens(i));
            Psim.Set(3,0.0);
            Psim.Set(4,DisPersDeriv.getIntens(i));
            AbsDispPoint.Set(1,SumDeriv.getIntens(i));
            AbsDispPoint.Set(2,DisPersDeriv.getIntens(i));
            SimData->Add(Psim);
            AbsDispData->Add(AbsDispPoint);
        }
    }
        else
            for (int i=0; i< Sum.getNpts(); i++)
            {
                Psim.Set(0,Sum.getField(i));
                AbsDispPoint.Set(0,Sum.getField(i));
                if (ccol != 0) Psim.Set(1,TempData.Get(i).Get(ccol));
                   else Psim.Set(1,0.0);
                Psim.Set(2,Sum.getIntens(i));
                Psim.Set(3,0.0);
                Psim.Set(4,DisPers.getIntens(i));
                AbsDispPoint.Set(1,Sum.getIntens(i));
                AbsDispPoint.Set(2,DisPers.getIntens(i));
                SimData->Add(Psim);
                AbsDispData->Add(AbsDispPoint);
            }

    SimCheck = true;
//    }
    double a,b;
    DataPoint Psim2(2);
    double DispAngle=0;
    if (ValidReal(EPRSimulationDialog->DispersEdit->Text, &DispAngle))
    {
        if (fabs(DispAngle) >0.001)
        {
            DispAngle *= (PI/180.0);
            for (int i=0; i<SimData->Getn(); i++)
            {
                Psim2 = SimData->Get(i);
//                AbsDispPoint = AbsDispData->Get(i);
                Psim2.Set(2, cos(DispAngle) * AbsDispData->Get(i).Get(1)
                               + sin(DispAngle) * AbsDispData->Get(i).Get(2));
                Psim2.Set(4, cos(DispAngle) * AbsDispData->Get(i).Get(2)
                               - sin(DispAngle) * AbsDispData->Get(i).Get(1));
                SimData->Set(i,Psim2);
            }
        }
    }

    if (ccol != 0) SimData->LinearRegres(2,1,&a,&b);
      else {
        a = 0.0;
        b =  1.0/SimData->Maxima().Get(2);
      }

    double chisqr=0.0;
    double OldValue;
    for (int i=0; i<SimData->Getn(); i++)
    {
        Psim2 = SimData->Get(i);
        OldValue = Psim2.Get(2);
        Psim2.Set(2, a + b * OldValue);
        Psim2.Set(3, Psim2.Get(1) - Psim2.Get(2));
        OldValue = Psim2.Get(4);
        Psim2.Set(4, b * OldValue);
        SimData->Set(i,Psim2);
        chisqr += Psim2.Get(3)*Psim2.Get(3);
    }

    chisqr /= SimData->Getn();
    if (sig2 > 0.0) chisqr /= sig2;
    return chisqr;

}

void __fastcall TMainForm::DisPersFitButtonClick(TObject *Sender)
{
    // dispersion change thing
    // A single cycle
    double Delta = 0.1; // degree change
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (ValidReal(EPRSimulationDialog->DispersEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += (CurrVal);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->DispersEdit->Text = (CurrVal - Delta);
        DisPersChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->DispersEdit->Text = (CurrVal + Delta);
        DisPersChanged = true;

        chih = CalcCycle();
        Line = "  high";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    EPRSimulationDialog->DispersEdit->Text = (CurrVal);

    DisPersChanged = true;

    chic = CalcCycle();
    DisPersChanged = false;
    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::PopFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.02; // relative ratio change of 2%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (ValidReal(EPRSimulationDialog->PopRatioEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->PopRatioEdit->Text = (CurrVal/(1 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->PopRatioEdit->Text = (CurrVal*(1.0 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature > 0)
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->PopRatioEdit->Text = (CurrVal);
    EPRSimulationDialog->PopRatioEditChange(this);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::WidthXFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.1;// width change. Let's take 1 G
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (!ValidReal(OptionsDialog->WChangeEdit->Text,&Delta)) return;
    if (ValidReal(EPRSimulationDialog->LwxEdit->Text, &CurrVal))
    {
        OutputMemo->Lines->Add("Optimizing Wx (1 cycle)");
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwxEdit->Text = (CurrVal - Delta);
        WidthChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwxEdit->Text = (CurrVal + Delta);
        WidthChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    EPRSimulationDialog->LwxEdit->Text = (CurrVal);

    WidthChanged = true;

    chic = CalcCycle();
    WidthChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::WidthYFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05;// width change. Let's take 0.5 G
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (!ValidReal(OptionsDialog->WChangeEdit->Text,&Delta)) return;
    if (ValidReal(EPRSimulationDialog->LwyEdit->Text, &CurrVal))
    {
        OutputMemo->Lines->Add("Optimizing Wy (1 cycle)");
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwyEdit->Text = (CurrVal - Delta);
        WidthChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwyEdit->Text = (CurrVal + Delta);
        WidthChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    EPRSimulationDialog->LwyEdit->Text = (CurrVal);

    WidthChanged = true;

    chic = CalcCycle();
    WidthChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::WidthZFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05;// width change. Let's take 0.5 G
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (!ValidReal(OptionsDialog->WChangeEdit->Text,&Delta)) return;
    if (ValidReal(EPRSimulationDialog->LwzEdit->Text, &CurrVal))
    {
        OutputMemo->Lines->Add("Optimizing Wz (1 cycle)");
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwzEdit->Text = (CurrVal - Delta);
        WidthChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwzEdit->Text = (CurrVal + Delta);
        WidthChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    EPRSimulationDialog->LwzEdit->Text = (CurrVal);

    WidthChanged = true;

    chic = CalcCycle();
    WidthChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::DFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta;
    if (!ValidReal(OptionsDialog->ZFSChangeEdit->Text,&Delta)) return;
            // D value change. Let's take 0.001 GHz
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (HPar.GetMultS() >2)
    {
        CurrVal = HPar.GetCF(0);
        OutputMemo->Lines->Add("Optimizing B20(D/3) (1 cycle)");
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        HPar.SetCF(0,(CurrVal - Delta));
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        HPar.SetCF(0, (CurrVal + Delta));
        HamilChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    HPar.SetCF(0,(CurrVal));

    HamilChanged = true;

    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::EFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.001;// E value change. Let's take 0.001 GHz
    if (!ValidReal(OptionsDialog->ZFSChangeEdit->Text,&Delta)) return;
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (HPar.GetMultS() >2)
    {
        CurrVal = HPar.GetCF(1);
        OutputMemo->Lines->Add("Optimizing B22(E) (1 cycle)");
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        HPar.SetCF(1,(CurrVal - Delta));
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        HPar.SetCF(1, (CurrVal + Delta));
        HamilChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    HPar.SetCF(1,(CurrVal));

    HamilChanged = true;

    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GxxFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing gxx (1 cycle)");
    GijFitCycle(0,0);
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GyyFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing gyy (1 cycle)");
    GijFitCycle(1,1);
}
//---------------------------------------------------------------------------


void TMainForm::GijFitCycle(int k, int l)
{
    // A single cycle
    double Delta = 0.00003;// g value change. Let's take 0.00003
    if (!ValidReal(OptionsDialog->gChangeEdit->Text,&Delta)) return;
    bool Axial = false;
    if (HamiltonianDialog->AxialCheckBox->Checked) Axial = true;
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;
    Tensor gT;
    if ((k < 0) || (k>2) || (l<0) || (l>2)) return;

        gT = HPar.Getg();
        CurrVal = gT.get(k,l);
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        gT.set(k,l,(CurrVal - Delta));
        gT.set(l,k,(CurrVal - Delta));

        if (Axial) gT.set(1,1,gT.get(0,0));
        HPar.Setg(gT);
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        gT.set(k,l,(CurrVal + Delta));
        gT.set(l,k,(CurrVal + Delta));
        if (Axial) gT.set(1,1,gT.get(0,0));
        HPar.Setg(gT);
        HamilChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);


    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2.0*Delta))
        CurrVal -= (Gradient/Curvature);
    else
        if (chih > chil) CurrVal -= (2.0*Delta);
           else CurrVal += (2.0*Delta);

    gT.set(k,l, CurrVal);
    gT.set(l,k, CurrVal);
    if (Axial) gT.set(1,1,gT.get(0,0));
    HPar.Setg(gT);

    HamilChanged = true;

    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();
}

void TMainForm::AijFitCycle(int k, int l)
{
    // A single cycle
    double Delta = 0.1;// A value change. Let's take 0.1 MHz
    if (!ValidReal(OptionsDialog->AChangeEdit->Text,&Delta)) return;
    bool Axial = false;
    if (HamiltonianDialog->AxialCheckBox->Checked) Axial = true;
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;
    Tensor gT;
    if ((k < 0) || (k>2) || (l<0) || (l>2)) return;

        gT = HPar.GetA(0);
        CurrVal = gT.get(k,l);
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        gT.set(k,l,(CurrVal - Delta));
        gT.set(l,k,(CurrVal - Delta));

        if (Axial) gT.set(1,1,gT.get(0,0));
        HPar.SetA(0,gT);
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        gT.set(k,l,(CurrVal + Delta));
        gT.set(l,k,(CurrVal + Delta));

        if (Axial) gT.set(1,1,gT.get(0,0));
        HPar.SetA(0,gT);
        HamilChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);


    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2.0*Delta))
        CurrVal -= (Gradient/Curvature);
    else
        if (chih > chil) CurrVal -= (2.0*Delta);
           else CurrVal += (2.0*Delta);

    gT.set(k,l, CurrVal);
    gT.set(l,k, CurrVal);

    if (Axial) gT.set(1,1,gT.get(0,0));
    HPar.SetA(0,gT);

    HamilChanged = true;

    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();
}

void TMainForm::GAngleFitCycle(int k)
{
    // A single cycle
    double Delta = 1;// g value change. Let's take 0.00003
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;
    Tensor gT;
    if ((k < 1) || (k>3)) return;
    Tensor Rot(3);
    Delta *= PI/180.0;

    Rot.set(0,0,cos(Delta));
    Rot.set(1,1,cos(Delta));
    Rot.set(0,1,sin(Delta));
    Rot.set(1,0,-sin(Delta));
    Rot.set(2,2,1.0);

    if (HPar.GetMultS() >2)
    {
        chic = CalcCycle();
        Line = "Zero rotation ";
//        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        gT = HPar.Getg();

        HPar.Setg(Rot.transpose()*(gT*Rot));
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   rotation ";
        Line +=  Delta;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        HPar.Setg(Rot*(gT*Rot.transpose()));
        HamilChanged = true;

        chih = CalcCycle();
        Line = " rotation ";
        Line += - Delta;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    Rot.set(0,0,cos(CurrVal));
    Rot.set(1,1,cos(CurrVal));
    Rot.set(0,1,sin(CurrVal));
    Rot.set(1,0,-sin(CurrVal));
    Rot.set(2,2,1.0);

    HPar.Setg(Rot*(gT*Rot.transpose()));

    HamilChanged = true;


    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();
}

void __fastcall TMainForm::GzzFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing gzz (1 cycle)");
    GijFitCycle(2,2);
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GxyFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing gxy (1 cycle)");
    GijFitCycle(0,1);

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GxzFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing gxz (1 cycle)");
    GijFitCycle(0,2);  //  for xy plane rotation 1

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GyzFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing gyz (1 cycle)");
    GijFitCycle(1,2);
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::OptionsMenuItemClick(TObject *Sender)
{
    OptionsDialog->ShowModal();
}
//---------------------------------------------------------------------------


void __fastcall TMainForm::GammaButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing XY plane angle");
    GAngleFitCycle(1);

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GaussianFit1Click(TObject *Sender)
{
// where are the gaussian routines ?
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::XZratioFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.02; // relative ratio change of 1%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (ValidReal(EPRSimulationDialog->XZratioEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->XZratioEdit->Text = (CurrVal/(1 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->XZratioEdit->Text = (CurrVal*(1.0 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature > 0)
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->XZratioEdit->Text = (CurrVal);
    EPRSimulationDialog->PopRatioEditChange(this);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------
void __fastcall TMainForm::YZratioFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.02; // relative ratio change of 1%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (ValidReal(EPRSimulationDialog->YZratioEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->YZratioEdit->Text = (CurrVal/(1 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->YZratioEdit->Text = (CurrVal*(1.0 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature > 0)
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->YZratioEdit->Text = (CurrVal);
    EPRSimulationDialog->PopRatioEditChange(this);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//-------------------------------------------
void __fastcall TMainForm::FitCyclesButtonClick(TObject *Sender)
{
    int icycle = 0;
    int MaxCycle;
    double MinImprovement;
    double Improvement = 0.0;
    double oldchi;
    AnsiString Line;
    if (FitCycleDialog->ShowModal() == mbOK)
    {
        if (!ValidInt(FitCycleDialog->NCycleEdit->Text, &MaxCycle)) return;
        if (!ValidReal(FitCycleDialog->ImprovementEdit->Text, &MinImprovement)) return;

        do
        {
            oldchi = CalcCycle();
            if (FitCycleDialog->DFitCheckBox->Checked)
                {
                    DFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->EFitCheckBox->Checked)
                {
                    EFitButtonClick(Sender);
                    Update();
                }
            
            if (FitCycleDialog->GisoCheckBox->Checked)
                {
                    GisoFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->GxxFitCheckBox->Checked)
                {
                    GxxFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->GyyFitCheckBox->Checked)
                {
                    GyyFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->GzzFitCheckBox->Checked)
                {
                    GzzFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->GxyFitCheckBox->Checked)
                {
                    GxyFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->GxzFitCheckBox->Checked)
                {
                    GxzFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->GyzFitCheckBox->Checked)
                {
                    GyzFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AbsDisFitCheckBox->Checked)
                {
                    DisPersFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->ISCPopRatioFitCheckBox->Checked) 
                {
                    PopFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->XZratioFitCheckBox->Checked) 
                {
                    XZratioFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->YZratioFitCheckBox->Checked) 
                {
                    YZratioFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->WxFitCheckBox->Checked)
                {
                    WidthXFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->WyFitCheckBox->Checked) 
                {
                    WidthYFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->WzFitCheckBox->Checked) 
                {
                    WidthZFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->WisoFitCheckBox->Checked)
                {
                    WisoFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AisoFitCheckBox->Checked)
                {
                    AisoFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->AxxFitCheckBox->Checked)
                {
                    AxxFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AyyFitCheckBox->Checked)
                {
                    AyyFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AzzFitCheckBox->Checked)
                {
                    AzzFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AxyFitCheckBox->Checked)
                {
                    AxyFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AxzFitCheckBox->Checked)
                {
                    AxzFitButtonClick(Sender);
                    Update();
                }
            if (FitCycleDialog->AyzFitCheckBox->Checked)
                {
                    AyzFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->KxFitCheckBox->Checked)
                {
                    KxFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->KyFitCheckBox->Checked)
                {
                    KyFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->KzFitCheckBox->Checked)
                {
                    KzFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->SLRxFitCheckBox->Checked)
                {
                    SLxFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->SLRyFitCheckBox->Checked)
                {
                    SLyFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->SLRzFitCheckBox->Checked)
                {
                    SLzFitButtonClick(Sender);
                    Update();
                }

            if (FitCycleDialog->PzFitCheckBox->Checked)
                {
                    PzFitButtonClick(Sender);
                    Update();
                }



            Improvement = 100.0*(oldchi-CalcCycle())/oldchi;     // in percent

            icycle++;

            Line = "*********************";
            OutputMemo->Lines->Add(Line);
            Line = "Cycle ";
            Line += icycle;
            OutputMemo->Lines->Add(Line);
            Line = "   Improvement ";
            Line += Improvement;
            OutputMemo->Lines->Add(Line);
            Update();

            Application->ProcessMessages();

        }

        while ((icycle < MaxCycle) && (Improvement > MinImprovement) && (!StopCheckBox->Checked));
        StopCheckBox->Checked = false;
    }
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::WisoFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05;// width change. Let's take 0.5 G
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (!ValidReal(OptionsDialog->WChangeEdit->Text,&Delta)) return;
    if (ValidReal(EPRSimulationDialog->LwxEdit->Text, &CurrVal))
    {
        OutputMemo->Lines->Add("Optimizing W (1 cycle)");
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwxEdit->Text = (CurrVal - Delta);
        EPRSimulationDialog->LwyEdit->Text = (CurrVal - Delta);
        EPRSimulationDialog->LwzEdit->Text = (CurrVal - Delta);
        WidthChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->LwxEdit->Text = (CurrVal + Delta);
        EPRSimulationDialog->LwyEdit->Text = (CurrVal + Delta);
        EPRSimulationDialog->LwzEdit->Text = (CurrVal + Delta);
        WidthChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal + Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if (Curvature != 0)
        CurrVal -= Gradient/Curvature;

    EPRSimulationDialog->LwxEdit->Text = (CurrVal);
    EPRSimulationDialog->LwyEdit->Text = (CurrVal);
    EPRSimulationDialog->LwzEdit->Text = (CurrVal);

    WidthChanged = true;

    chic = CalcCycle();
    WidthChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::GisoFitButtonClick(TObject *Sender)
{
    double Delta = 0.00003;// g value change. Let's take 0.00003
    if (!ValidReal(OptionsDialog->gChangeEdit->Text,&Delta)) return;
    double CurrXX, CurrYY, CurrZZ;
    double chil, chih, chic;
    AnsiString Line;
    Tensor gT;

        gT = HPar.Getg();
        CurrXX = gT.get(0,0);
        CurrYY = gT.get(1,1);
        CurrZZ = gT.get(2,2);

        chic = CalcCycle();
        Line = "center ";
        Line += CurrXX;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        gT.set(0,0,(CurrXX - Delta));
        gT.set(1,1,(CurrYY - Delta));
        gT.set(2,2,(CurrZZ - Delta));
        HPar.Setg(gT);
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += ( - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        gT.set(0,0,(CurrXX + Delta));
        gT.set(1,1,(CurrYY + Delta));
        gT.set(2,2,(CurrZZ + Delta));
        HPar.Setg(gT);
        HamilChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += ( Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);


    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2.0*Delta))
    {
        CurrXX -= (Gradient/Curvature);
        CurrYY -= (Gradient/Curvature);
        CurrZZ -= (Gradient/Curvature);
    }
     else
        if (chih > chil)
        {
            CurrXX -= (2.0*Delta);
            CurrYY -= (2.0*Delta);
            CurrZZ -= (2.0*Delta);
        }
           else
           {
                CurrXX += (2.0*Delta);
                CurrYY += (2.0*Delta);
                CurrZZ += (2.0*Delta);
           }

    gT.set(0,0,(CurrXX ));
    gT.set(1,1,(CurrYY ));
    gT.set(2,2,(CurrZZ ));
    HPar.Setg(gT);

    HamilChanged = true;

    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrXX+CurrYY+CurrZZ)/3;
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();

}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

void TMainForm::SimulateLorentz()
{
    double a[100];
    double da[100];
    int ai[100];
    double low, high, x, y;
    int na = SpectralLineStartDialog->GetParameters(a, ai);

    int npts = SpectralLineStartDialog->GetLimits(&low, &high);
    if (SimData != NULL) delete SimData;
    SimData = new DataArray(npts,4);
    DataArray TempData = ExpData->ReSample(low,high,npts);
    int ccol;

    DataPoint P(4);
    if ((ccol = ColumnComboBox->ItemIndex)>0)
    {}
    double chi2 = 0.0;
    for (int i=0;i<npts;i++)
    {
        x = TempData.Get(i).Get(0);
        P.Set(0,x);
        if (SpectralLineStartDialog->DerivCheckBox->Checked)
            MultDerivPhasedSplitLorentz(x, a, &y, da, na);
          else
            MultPhasedSplitLorentz(x, a, &y, da, na);
        if (ccol>0) P.Set(1,TempData.Get(i).Get(ccol));
        P.Set(2,y);
        P.Set(3,y - TempData.Get(i).Get(ccol));
        chi2 += (P.Get(1)-P.Get(2))*(P.Get(1)-P.Get(2));
        SimData->Add(P);
    }
//    MultDerivLorentz();
    OutputMemo->Lines->Add(chi2);
    SimCheck = true;
    Invalidate();
    return;
}

///////////////////////////////////////
//

void TMainForm::SimulateGauss()
{
    double a[100];
    double da[100];
    int ai[100];
    double low, high, x, y;
    int na = SpectralLineStartDialog->GetParameters(a, ai);

    int npts = SpectralLineStartDialog->GetLimits(&low, &high);
    if (SimData != NULL) delete SimData;
    SimData = new DataArray(npts,4);
    DataArray TempData = ExpData->ReSample(low,high,npts);
    int ccol;

    DataPoint P(4);
    if ((ccol = ColumnComboBox->ItemIndex)>0)
    {}
    double chi2 = 0.0;
    for (int i=0;i<npts;i++)
    {
        x = TempData.Get(i).Get(0);
        P.Set(0,x);
        if (SpectralLineStartDialog->DerivCheckBox->Checked)
            MultDerivPhasedSplitGaussian(x, a, &y, da, na);
          else
            MultPhasedSplitGaussian(x, a, &y, da, na);
        if (ccol>0) P.Set(1,TempData.Get(i).Get(ccol));
        P.Set(2,y);
        P.Set(3,y - TempData.Get(i).Get(ccol));
        chi2 += (P.Get(1)-P.Get(2))*(P.Get(1)-P.Get(2));
        SimData->Add(P);
    }
//    MultDerivLorentz();
    OutputMemo->Lines->Add(chi2);
    SimCheck = true;
    Invalidate();
    return;
}


void TMainForm::FitLorentz(int MaxCycle)
{
    double a[100];
    double da[100];
    int ai[100];

    double low, high, x, y;
    double chisqr;
    double alamb;

    double **xdat;
    double *ydat;
    double *ysig;

    int na = SpectralLineStartDialog->GetParameters(a, ai);
    double **covar;
    double **alpha;
    covar = new double*[na];
    alpha = new double*[na];
    for (int i = 0; i<na; i++)
    {
        covar[i] = new double[na];
        alpha[i] = new double[na];
    }

    int npts = SpectralLineStartDialog->GetLimits(&low, &high);
    if (SimData != NULL) delete SimData;
    SimData = new DataArray(npts,4);
    DataArray TempData = ExpData->ReSample(low,high,npts);
    xdat = new double*[npts];
    ydat = new double[npts];
    ysig = new double[npts];
    int ccol = ColumnComboBox->ItemIndex;

    DataPoint P(4);
    alamb  = -1.0;
    if (ccol >0)
    {
        for (int i=0; i<npts; i++)
        {
            xdat[i] = new double[1];
            xdat[i][0] = TempData.Get(i).Get(0);
            ydat[i] = TempData.Get(i).Get(ccol);
            ysig[i] = 1.0;
        }
        for (int i=0; i< MaxCycle; i++)
        {
        if (SpectralLineStartDialog->DerivCheckBox->Checked)
            fitcycle(xdat,ydat,ysig,npts,1,a,ai,na,covar, alpha, &chisqr,MultDerivPhasedSplitLorentz, &alamb);
          else
            fitcycle(xdat,ydat,ysig,npts,1,a,ai,na,covar, alpha, &chisqr,MultPhasedSplitLorentz, &alamb);
        }
        alamb = 0.0;  // This is needed for error-estimation...
        if (SpectralLineStartDialog->DerivCheckBox->Checked)
            fitcycle(xdat,ydat,ysig,npts,1,a,ai,na,covar, alpha, &chisqr,MultDerivPhasedSplitLorentz, &alamb);
          else
            fitcycle(xdat,ydat,ysig,npts,1,a,ai,na,covar, alpha, &chisqr,MultPhasedSplitLorentz, &alamb);

    }
    OutputMemo->Lines->Add(chisqr);
    SpectralLineStartDialog->SetParameters(na, a, ai);

    for (int i=0;i<npts;i++)
    {
        x = low + i*(high-low)/(npts-1);
        P.Set(0,x);
        if (SpectralLineStartDialog->DerivCheckBox->Checked)
            MultDerivPhasedSplitLorentz(x, a, &y, da, na);
          else
            MultPhasedSplitLorentz(x, a, &y, da, na);

        if (ccol>0) P.Set(1,TempData.Get(i).Get(ccol));
        P.Set(2,y);
        P.Set(3,y - TempData.Get(i).Get(ccol));
        SimData->Add(P);
    }
    for (int i=0; i<npts; i++) delete xdat[i];
    delete[] xdat;
    delete[] ydat;
    delete[] ysig;
    for (int i=0; i<na; i++) delete[] alpha[i];
    for (int i=0; i<na; i++) delete[] covar[i];
    delete[] covar;
    delete[] alpha;

//    MultDerivLorentz();
    SimCheck = true;
    Invalidate();
    return;
}



void __fastcall TMainForm::HighSpin1Click(TObject *Sender)
{
// Prepare Dialog
    PopMode = 0;           // Boltzman
    DerivMode = 1;
    LineFitMode = 0;
    NumberChanged = true;  //always do calculation
    EPRSimulationDialog->SetHamilChange(HamilChanged);
    EPRSimulationDialog->SetPopChange(PopChanged);
    EPRSimulationDialog->SetWidthChange(WidthChanged);
    EPRSimulationDialog->SetNumberChange(NumberChanged);
    EPRSimulationDialog->SetDisPersChange(DisPersChanged);

    EPRSimulationDialog->BoltzmannCheckBox->Checked =true;

    if (EPRSimulationDialog->ShowModal() == mbOK)
    {
// first check which parameters changed
        PopChanged = EPRSimulationDialog->GetPopChange();
        NumberChanged = EPRSimulationDialog->GetNumberChange();
        WidthChanged = EPRSimulationDialog->GetWidthChange();
        DisPersChanged = EPRSimulationDialog->GetDisPersChange();
        HamilChanged = EPRSimulationDialog->GetHamilChange();
    }
    else return;
    if (EPRSimulationDialog->IntegrRadioButton->Checked == true)
        DerivMode = 0; else DerivMode = 1;
// Set the parameters of the SpinHamiltonian
    Update();
    double chi2 = HighSpinCalc();
// Open Dialog
// Nothing implemented yet.
// 1. Take Theta, phi
// 2.   Take field
// 3.       Build hamiltonian
// 4.       Calculate transition frequencies and probabilities
// 5.       Calculate contribution to spectrum at specified frequency
// 6.   Next field value
// 7. Next theta, phi

//  Note number of transitions is S*(2S+1)
    OutputMemo->Lines->Add("done...chi^2:");
    OutputMemo->Lines->Add(chi2);

    Invalidate();

}
//---------------------------------------------------------------------------
/*
void __fastcall TDataFile::FormMouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
    if (Button == mbRight) // if right MouseButton pressed
    {
        FocusRect = new TRect();
        FocusRect->Left = X;
        FocusRect->Right = X;
        FocusRect->Top = Y;
        FocusRect->Bottom = Y;
        Canvas->DrawFocusRect(*FocusRect);
    }
}

//---------------------------------------------------------------------------

void __fastcall TDataFile::FormMouseMove(TObject *Sender,
      TShiftState Shift, int X, int Y)
{
	if (!FocusRect) return;
	Canvas->DrawFocusRect(*FocusRect);
	FocusRect->Right = X;
    FocusRect->Bottom = Y;
	Canvas->DrawFocusRect(*FocusRect);
}
//---------------------------------------------------------------------------

void __fastcall TDataFile::FormMouseUp(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
	if (!FocusRect) return;
	Canvas->DrawFocusRect(*FocusRect);
	if ((FocusRect->Right > FocusRect->Left+3) && (FocusRect->Bottom > FocusRect->Top + 3))
	{
		Plot->SetAutoRange(0);
		TPoint BottomLeft(FocusRect->Left, FocusRect->Bottom);
		TPoint TopRight(FocusRect->Right, FocusRect->Top);
		Plot->SetFixedRanges(Plot->GetUserCoordinates(BottomLeft),
										Plot->GetUserCoordinates(TopRight));
		Invalidate();
	}
	delete FocusRect;
	FocusRect = 0;
}
*/

void __fastcall TMainForm::FormMouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
    if (Button == mbRight) // if right MouseButton pressed
    {
        FocusRect = new TRect();
        FocusRect->Left = X;
        FocusRect->Right = X;
        FocusRect->Top = Y;
        FocusRect->Bottom = Y;
        Canvas->DrawFocusRect(*FocusRect);
    }
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FormMouseMove(TObject *Sender,
      TShiftState Shift, int X, int Y)
{
	if (!FocusRect) return;
	Canvas->DrawFocusRect(*FocusRect);
	FocusRect->Right = X;
    FocusRect->Bottom = Y;
	Canvas->DrawFocusRect(*FocusRect);
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FormMouseUp(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
	if (!FocusRect) return;
	Canvas->DrawFocusRect(*FocusRect);
	if ((FocusRect->Right > FocusRect->Left+3) && (FocusRect->Bottom > FocusRect->Top + 3))
	{
		SimPlot->SetAutoRange(0);
		TPoint BottomLeft(FocusRect->Left, FocusRect->Bottom);
		TPoint TopRight(FocusRect->Right, FocusRect->Top);
		SimPlot->SetFixedRanges(SimPlot->GetUserCoordinates(BottomLeft),
										SimPlot->GetUserCoordinates(TopRight));
		Invalidate();
	}
	delete FocusRect;
	FocusRect = 0;
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::SpeedButton1Click(TObject *Sender)
{
		SimPlot->SetAutoRange(7);
        Invalidate();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AisoFitButtonClick(TObject *Sender)
{
    double Delta = 0.00003;// g value change. Let's take 0.00003
    if (!ValidReal(OptionsDialog->AChangeEdit->Text,&Delta)) return;
    double CurrXX, CurrYY, CurrZZ;
    double chil, chih, chic;
    AnsiString Line;
    Tensor gT;

        gT = HPar.GetA(0);
        CurrXX = gT.get(0,0);
        CurrYY = gT.get(1,1);
        CurrZZ = gT.get(2,2);

        chic = CalcCycle();
        Line = "center ";
        Line += CurrXX;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        gT.set(0,0,(CurrXX - Delta));
        gT.set(1,1,(CurrYY - Delta));
        gT.set(2,2,(CurrZZ - Delta));
        HPar.SetA(0,gT);
        HamilChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += ( - Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        gT.set(0,0,(CurrXX + Delta));
        gT.set(1,1,(CurrYY + Delta));
        gT.set(2,2,(CurrZZ + Delta));
        HPar.SetA(0,gT);
        HamilChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += ( Delta);
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);


    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2.0*Delta))
    {
        CurrXX -= (Gradient/Curvature);
        CurrYY -= (Gradient/Curvature);
        CurrZZ -= (Gradient/Curvature);
    }
     else
        if (chih > chil)
        {
            CurrXX -= (2.0*Delta);
            CurrYY -= (2.0*Delta);
            CurrZZ -= (2.0*Delta);
        }
           else
           {
                CurrXX += (2.0*Delta);
                CurrYY += (2.0*Delta);
                CurrZZ += (2.0*Delta);
           }

    gT.set(0,0,(CurrXX ));
    gT.set(1,1,(CurrYY ));
    gT.set(2,2,(CurrZZ ));
    HPar.SetA(0,gT);

    HamilChanged = true;

    chic = CalcCycle();

    HamilChanged = false;

    Line = "New Value ";
    Line += (CurrXX+CurrYY+CurrZZ)/3;
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    HamiltonianDialog->SetParameters(HPar);
    Invalidate();

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AxxFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing Axx (1 cycle)");
    AijFitCycle(0,0);
// to do
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AyyFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing Ayy (1 cycle)");
    AijFitCycle(1,1);

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AzzFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing Azz (1 cycle)");
    AijFitCycle(2,2);

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AxyFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing Axy (1 cycle)");
    AijFitCycle(0,1);

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AxzFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing Axz (1 cycle)");
    AijFitCycle(0,2);

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::AyzFitButtonClick(TObject *Sender)
{
    OutputMemo->Lines->Add("Optimizing Ayz (1 cycle)");
    AijFitCycle(1,2);
}
//---------------------------------------------------------------------------


void __fastcall TMainForm::KxFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05; // relative ratio change of 5%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    Line = "Optimizing kx ";
    OutputMemo->Lines->Add(Line);
    if (ValidReal(EPRSimulationDialog->kxEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->kxEdit->Text = (CurrVal/(1 + Delta));
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->kxEdit->Text = (CurrVal*(1.0 + Delta));
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->kxEdit->Text = (CurrVal);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::KyFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05; // relative ratio change of 5%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    Line = "Optimizing ky ";
    OutputMemo->Lines->Add(Line);
    if (ValidReal(EPRSimulationDialog->kyEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->kyEdit->Text = (CurrVal/(1 + Delta));
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->kyEdit->Text = (CurrVal*(1.0 + Delta));
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->kyEdit->Text = (CurrVal);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::KzFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05; // relative ratio change of 5%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    Line = "Optimizing kz ";
    OutputMemo->Lines->Add(Line);
    if (ValidReal(EPRSimulationDialog->kzEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->kzEdit->Text = (CurrVal/(1 + Delta));
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->kzEdit->Text = (CurrVal*(1.0 + Delta));
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->kzEdit->Text = (CurrVal);
//    EPRSimulationDialog->PopRatioEditChange(this);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::SLxFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05; // relative ratio change of 5%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    Line = "Optimizing SLR along x ";
    OutputMemo->Lines->Add(Line);
    if (ValidReal(EPRSimulationDialog->SLRxEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->SLRxEdit->Text = (CurrVal/(1 + Delta));
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->SLRxEdit->Text = (CurrVal*(1.0 + Delta));
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->SLRxEdit->Text = (CurrVal);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::SLyFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05; // relative ratio change of 5%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    Line = "Optimizing SLR along y ";
    OutputMemo->Lines->Add(Line);
    if (ValidReal(EPRSimulationDialog->SLRyEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->SLRyEdit->Text = (CurrVal/(1 + Delta));
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->SLRyEdit->Text = (CurrVal*(1.0 + Delta));
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->SLRyEdit->Text = (CurrVal);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::SLzFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.05; // relative ratio change of 5%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    Line = "Optimizing SLR along z ";
    OutputMemo->Lines->Add(Line);
    if (ValidReal(EPRSimulationDialog->SLRzEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->SLRzEdit->Text = (CurrVal/(1 + Delta));
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->SLRzEdit->Text = (CurrVal*(1.0 + Delta));
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->SLRzEdit->Text = (CurrVal);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

void __fastcall TMainForm::PzFitButtonClick(TObject *Sender)
{
    // A single cycle
    double Delta = 0.02; // relative ratio change of 2%
    double CurrVal;
    double chil, chih, chic;
    AnsiString Line;

    if (ValidReal(EPRSimulationDialog->PzEdit->Text, &CurrVal))
    {
        chic = CalcCycle();
        Line = "center ";
        Line += CurrVal;
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chic;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->PzEdit->Text = (CurrVal/(1 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chil = CalcCycle();
        Line = "   low ";
        Line += (CurrVal/(1+Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chil;
        OutputMemo->Lines->Add(Line);

        EPRSimulationDialog->PzEdit->Text = (CurrVal*(1.0 + Delta));
        EPRSimulationDialog->PopRatioEditChange(this);
        PopChanged = true;

        chih = CalcCycle();
        Line = "  high ";
        Line += (CurrVal*(1.0 + Delta));
        OutputMemo->Lines->Add(Line);

        Line = "   X2 = ";
        Line += chih;
        OutputMemo->Lines->Add(Line);

    }

    double Gradient = (chih-chil)/(2*Delta);
    double Curvature = (chih+chil-2.0*chic)/(Delta*Delta);

    if ((Curvature > 0) && (fabs(Gradient/Curvature) < 2*Delta))
        CurrVal *= (1-Gradient/Curvature);
    else
        if (chih > chil) CurrVal /= (1.0+2.0*Delta);
           else CurrVal *= (1.0+2.0*Delta);
    EPRSimulationDialog->PzEdit->Text = (CurrVal);
    EPRSimulationDialog->PopRatioEditChange(this);

    PopChanged = true;

    chic = CalcCycle();
    PopChanged = false;

    Line = "New Value ";
    Line += (CurrVal);
    OutputMemo->Lines->Add(Line);

    Line = "   X2 = ";
    Line += chic;
    OutputMemo->Lines->Add(Line);
    Invalidate();


}
//---------------------------------------------------------------------------

