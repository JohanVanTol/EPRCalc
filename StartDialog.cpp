//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <sysutils.hpp>
#include <math.h>
#include "StartDialog.h"
#include "Valid.h"
#include "EPRCalcMain2.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
extern int MAXPAR;

TSpectralLineStartDialog *SpectralLineStartDialog;
//---------------------------------------------------------------------------
__fastcall TSpectralLineStartDialog::TSpectralLineStartDialog(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------



void __fastcall TSpectralLineStartDialog::FormCreate(TObject *Sender)
{
    NumberOfLinesListBox->ItemIndex = 0;
    EditArray = new TEdit*[63];
    CheckArray = new TCheckBox*[63];
    MultArray = new TComboBox*[10];

	EditArray[0] = A0Edit;                   //offset
	EditArray[1] = A1Edit;                   //slope

	EditArray[2] = A2Edit;                   //amplitude
	EditArray[3] = A3Edit;                   //position
	EditArray[4] = A4Edit;                   //width
	EditArray[5] = F1Edit;                   //phase
	EditArray[6] = S1Edit;                   //splitting
	EditArray[7] = Alpha1Edit;               //Voigt parameter

    EditArray[8] = A5Edit;
	EditArray[9] = A6Edit;
    EditArray[10] = A7Edit;
    EditArray[11] = F2Edit;
    EditArray[12] = S2Edit;
    EditArray[13] = Alpha2Edit;

    EditArray[14] = A8Edit;
    EditArray[15] = A9Edit;
    EditArray[16] = A10Edit;
    EditArray[17] = F3Edit;
    EditArray[18] = S3Edit;
	EditArray[19] = Alpha3Edit;

    EditArray[20] = A11Edit;
    EditArray[21] = A12Edit;
    EditArray[22] = A13Edit;
    EditArray[23] = F4Edit;
    EditArray[24] = S4Edit;
	EditArray[25] = Alpha4Edit;

    EditArray[26] = A14Edit;
    EditArray[27] = A15Edit;
    EditArray[28] = A16Edit;
    EditArray[29] = F5Edit;
    EditArray[30] = S5Edit;
	EditArray[31] = Alpha5Edit;

    EditArray[32] = A17Edit;
    EditArray[33] = A18Edit;
    EditArray[34] = A19Edit;
    EditArray[35] = F6Edit;
    EditArray[36] = S6Edit;
	EditArray[37] = Alpha6Edit;

    EditArray[38] = A20Edit;
    EditArray[39] = A21Edit;
    EditArray[40] = A22Edit;
    EditArray[41] = F7Edit;
    EditArray[42] = S7Edit;
	EditArray[43] = Alpha7Edit;

    EditArray[44] = A23Edit;
    EditArray[45] = A24Edit;
    EditArray[46] = A25Edit;
    EditArray[47] = F8Edit;
    EditArray[48] = S8Edit;
	EditArray[49] = Alpha8Edit;

	EditArray[50] = A26Edit;
	EditArray[51] = A27Edit;
	EditArray[52] = A28Edit;
	EditArray[53] = F9Edit;
	EditArray[54] = S9Edit;
	EditArray[55] = Alpha9Edit;

	EditArray[56] = A29Edit;
	EditArray[57] = A30Edit;
	EditArray[58] = A31Edit;
	EditArray[59] = F10Edit;
	EditArray[60] = S10Edit;
	EditArray[61] = Alpha10Edit;



	CheckArray[0] = A0CheckBox;
	CheckArray[1] = A1CheckBox;

	CheckArray[2] = A2CheckBox;
	CheckArray[3] = A3CheckBox;
	CheckArray[4] = A4CheckBox;
	CheckArray[5] = F1CheckBox;
	CheckArray[6] = S1CheckBox;
	CheckArray[7] = Alpha1CheckBox;

	CheckArray[8] = A5CheckBox;
	CheckArray[9] = A6CheckBox;
	CheckArray[10] = A7CheckBox;
	CheckArray[11] = F2CheckBox;
	CheckArray[12] = S2CheckBox;
	CheckArray[13] = Alpha2CheckBox;

	CheckArray[14] = A8CheckBox;
	CheckArray[15] = A9CheckBox;
	CheckArray[16] = A10CheckBox;
	CheckArray[17] = F3CheckBox;
	CheckArray[18] = S3CheckBox;
	CheckArray[19] = Alpha3CheckBox;

	CheckArray[20] = A11CheckBox;
	CheckArray[21] = A12CheckBox;
	CheckArray[22] = A13CheckBox;
	CheckArray[23] = F4CheckBox;
	CheckArray[24] = S4CheckBox;
	CheckArray[25] = Alpha4CheckBox;

	CheckArray[26] = A14CheckBox;
	CheckArray[27] = A15CheckBox;
	CheckArray[28] = A16CheckBox;
	CheckArray[29] = F5CheckBox;
	CheckArray[30] = S5CheckBox;
	CheckArray[31] = Alpha5CheckBox;

	CheckArray[32] = A17CheckBox;
	CheckArray[33] = A18CheckBox;
	CheckArray[34] = A19CheckBox;
	CheckArray[35] = F6CheckBox;
	CheckArray[36] = S6CheckBox;
	CheckArray[37] = Alpha6CheckBox;

	CheckArray[38] = A20CheckBox;
	CheckArray[39] = A21CheckBox;
	CheckArray[40] = A22CheckBox;
	CheckArray[41] = F7CheckBox;
	CheckArray[42] = S7CheckBox;
	CheckArray[43] = Alpha7CheckBox;

	CheckArray[44] = A23CheckBox;
	CheckArray[45] = A24CheckBox;
	CheckArray[46] = A25CheckBox;
	CheckArray[47] = F8CheckBox;
	CheckArray[48] = S8CheckBox;
	CheckArray[49] = Alpha8CheckBox;

	CheckArray[50] = A26CheckBox;
	CheckArray[51] = A27CheckBox;
	CheckArray[52] = A28CheckBox;
	CheckArray[53] = F9CheckBox;
	CheckArray[54] = S9CheckBox;
	CheckArray[55] = Alpha9CheckBox;

    CheckArray[56] = A29CheckBox;
    CheckArray[57] = A30CheckBox;
    CheckArray[58] = A31CheckBox;
    CheckArray[59] = F10CheckBox;
    CheckArray[60] = S10CheckBox;
	CheckArray[61] = Alpha10CheckBox;

    MultArray[0] = Mult1ListBox;
    MultArray[1] = Mult2ListBox;
    MultArray[2] = Mult3ListBox;
    MultArray[3] = Mult4ListBox;
    MultArray[4] = Mult5ListBox;
    MultArray[5] = Mult6ListBox;
    MultArray[6] = Mult7ListBox;
    MultArray[7] = Mult8ListBox;
    MultArray[8] = Mult9ListBox;
    MultArray[9] = Mult10ListBox;


    par = new double[73];
    fixed = new int[73];
}
//---------------------------------------------------------------------------
void TSpectralLineStartDialog::SetParameters(int n, double* param, int* fixedpar)
{
//  We have 2 parameters for offset and a linear term
//  We have 7 parameters per line:
//      Amplitude
//      Position
//      Width
//      Phase
//      Multiplicity
//      Splitting
//      Voigt parameter (not implemented yet)
//

    if ((n>72) || (n<2)) return;
    npar = n;

    for (int i=0; i<2;i++)
    {
        par[i] = param[i];
		EditArray[i]->Text = FloatToStrF((long double)param[i], ffGeneral, 6, 0);
		fixed[i] = fixedpar[i];
        if (fixed[i]==1) CheckArray[i]->Checked = true;
           else CheckArray[i]->Checked = false;
    }
    for (int iline = 0; iline < npar/7 ; iline++)
    {
		for (int i=0; i<6;i++)    //amplitude, center, width, phase, splitting, alpha
		{
            par[2 + 7*iline + i] = param[2 + 7*iline + i];
			if ((i==2) || (i==4))
				EditArray[2+6*iline+i]->Text = FloatToStrF((long double)param[2+7*iline+i]*1000.0, ffGeneral, 6, 0);
			 else {
				if (i==3)
					EditArray[2+6*iline+i]->Text = FloatToStrF((long double)param[2+7*iline+i]*180.0/M_PI, ffGeneral, 6, 0);
				 else
					EditArray[2+6*iline+i]->Text = FloatToStrF((long double)param[2+7*iline+i], ffGeneral, 8, 0);
			 }
			fixed[2 + 7*iline +i] = fixedpar[2 + 7*iline + i];
			if (fixed[2+7*iline+i]==1) CheckArray[2+6*iline+i]->Checked = true;
			   else CheckArray[2+6*iline+i]->Checked = false;
		}
		par[8+7*iline] = param[8+7*iline];
		MultArray[iline]->ItemIndex = int(par[8+7*iline])-1;
		if (MultArray[iline]->ItemIndex > 0)
        {
            EditArray[6+iline*6]->Visible = true;
            CheckArray[6+iline*6]->Visible = true;
        }
          else
          {
            EditArray[6+iline*6]->Visible = false;
            CheckArray[6+iline*6]->Visible = false;
          }
        fixed[8+7*iline] = 0;
    }
}

int TSpectralLineStartDialog::GetParameters(double* param, int* fixedpar)
{
    int nline = NumberOfLinesListBox->ItemIndex+1;
    npar = nline*7 + 2;
    for (int i = 0; i< 2; i++)
    {
		par[i] = EditArray[i]->Text.ToDouble();
		param[i] = par[i];
        if (CheckArray[i]->Checked) fixed[i] = 1;
          else fixed[i] = 0;
        fixedpar[i] = fixed[i];
    }

    for (int iline = 0; iline<nline; iline++)
    {
        for (int i = 0; i< 6; i++)
        {
			if ((i==2) || (i==4)) par[2+iline*7+i] = EditArray[2+iline*6+i]->Text.ToDouble()/1000.0;
			  else
			   if (i==3) par[2+iline*7+i] = EditArray[2+iline*6+i]->Text.ToDouble()*M_PI/180.0;
			   else 	par[2+iline*7+i] = EditArray[2+iline*6+i]->Text.ToDouble();

			param[2+iline*7+i] = par[2+iline*7+i];
            if (CheckArray[2+iline*6+i]->Checked) fixed[2+iline*7+i] = 1;
              else fixed[2+iline*7+i] = 0;
            fixedpar[2+iline*7+i] = fixed[2+iline*7+i];
        }
        par[8+iline*7] = MultArray[iline]->ItemIndex+1;
        fixed[8+iline*7] = 0;
        param[8+iline*7] = par[8+iline*7];
        fixedpar[8+iline*7] = fixed[8+iline*7];

	}
	return npar;
}


int TSpectralLineStartDialog::UpdateForm()
{
	NumberOfLinesListBoxChange(this);
}



void __fastcall TSpectralLineStartDialog::NumberOfLinesListBoxChange(
	  TObject *Sender)
{
	for (int i=2; i<62;i++)
	{
		EditArray[i]->Enabled = false;
		EditArray[i]->Visible = false;
		CheckArray[i]->Enabled = false;
		CheckArray[i]->Visible = false;
	}
	for (int i=0; i < 10; i++) {
		MultArray[i]->Enabled = false;
		MultArray[i]->Visible = false;
	}
	for (int k=0; k<=NumberOfLinesListBox->ItemIndex; k++)
	{
		EditArray[6*k+2]->Enabled = true;
		EditArray[6*k+3]->Enabled = true;
		EditArray[6*k+4]->Enabled = true;
		EditArray[6*k+5]->Enabled = true;
		EditArray[6*k+6]->Enabled = true;
		EditArray[6*k+7]->Enabled = true;
		CheckArray[6*k+2]->Enabled = true;
		CheckArray[6*k+3]->Enabled = true;
		CheckArray[6*k+4]->Enabled = true;
		CheckArray[6*k+5]->Enabled = true;
		CheckArray[6*k+6]->Enabled = true;
		CheckArray[6*k+7]->Enabled = true;
		MultArray[k]->Enabled = true;

		EditArray[6*k + 2]->Visible = true;
		EditArray[6*k + 3]->Visible = true;
		EditArray[6*k + 4]->Visible = true;
		EditArray[6*k + 5]->Visible = true;
		EditArray[6*k + 6]->Visible = true;
		EditArray[6*k + 7]->Visible = true;
		CheckArray[6*k + 2]->Visible = true;
		CheckArray[6*k + 3]->Visible = true;
		CheckArray[6*k + 4]->Visible = true;
		CheckArray[6*k + 5]->Visible = true;
		CheckArray[6*k + 6]->Visible = true;
		CheckArray[6*k + 7]->Visible = true;

		MultArray[k]->Visible = true;
	}
	if (ShapeComboBox->ItemIndex < 2) // If NOT Voigt
	{
		for (int k=0; k < 10; k++)
		{
			EditArray[6*k+7]->Visible = false;
			CheckArray[6*k+7]->Visible = false;
		}
	}
	StartParamGroupBox->Height = 136 + 23*NumberOfLinesListBox->ItemIndex;
	Height = 252 + 23*NumberOfLinesListBox->ItemIndex;
}
//---------------------------------------------------------------------------


void __fastcall TSpectralLineStartDialog::FormDestroy(TObject *Sender)
{
	delete[] EditArray;
	delete[] CheckArray;
	delete[] MultArray;
	delete[] par;
	delete[] fixed;
}
//---------------------------------------------------------------------------


void __fastcall TSpectralLineStartDialog::BitBtn1Click(TObject *Sender)
{
    ModalResult = mbOK;
    Close();
}
//---------------------------------------------------------------------------

int TSpectralLineStartDialog::GetLimits(double* start, double* stop)
{
    int np;

    if (!ValidInt(NptsEdit->Text,&np)) return -1;
    if (!ValidReal(StartEdit->Text,start)) return -1;
    if (!ValidReal(StopEdit->Text,stop)) return -1;
    return np;
}


void __fastcall TSpectralLineStartDialog::ChiSqrButtonClick(TObject *Sender)
{
    TMainForm *MF;
    if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
    {
        if (ShapeComboBox->ItemIndex == 0)
            MF->SimulateLorentz();
		if (ShapeComboBox->ItemIndex == 1)
			MF->SimulateGauss();
		if (ShapeComboBox->ItemIndex == 2)
			MF->SimulatePseudoVoigt();

    }
}
//---------------------------------------------------------------------------

void __fastcall TSpectralLineStartDialog::FitCycleButtonClick(
      TObject *Sender)
{
	TMainForm *MF;
	if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
	{
		if (ShapeComboBox->ItemIndex == 0)  MF->FitLorentz(1);
		if (ShapeComboBox->ItemIndex == 1)  MF->FitGauss(1);
		if (ShapeComboBox->ItemIndex == 2)  MF->FitPseudoVoigt(1);
	}

}
//---------------------------------------------------------------------------

void __fastcall TSpectralLineStartDialog::Fit10ButtonClick(TObject *Sender)
{
	TMainForm *MF;
	if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
	{
		if (ShapeComboBox->ItemIndex == 0)  MF->FitLorentz(10);
		if (ShapeComboBox->ItemIndex == 1)  MF->FitGauss(10);
		if (ShapeComboBox->ItemIndex == 2)  MF->FitPseudoVoigt(10);
	}

}
//---------------------------------------------------------------------------

int TSpectralLineStartDialog::Write(ofstream *somefile)
{
	if (ShapeComboBox->ItemIndex == 0) {
					*somefile << "Fit to Lorentzian: " << endl;
	}
	else if (ShapeComboBox->ItemIndex == 1) {
			*somefile << "Fit to Gaussian: " << endl;
		 }
		 else  {
	*somefile << "Fit to Pseudo-Voigt: " << endl;}

	*somefile << "Offset : " << A0Edit->Text.c_str() << endl;
	*somefile << "Linear term : " << A1Edit->Text.c_str() << endl;
	*somefile << "Line          Amplitude    Center     Width    PhaseAngle(rad)  Splitting  VoigtFactor" << endl;

	for (int i=0 ; i<=NumberOfLinesListBox->ItemIndex ; i++)
	{
		*somefile << "Line " << i+1 << "     ";
		*somefile << EditArray[6*i+2]->Text.ToDouble() << "    ";
		*somefile << EditArray[6*i+3]->Text.ToDouble() << "    ";
		*somefile << EditArray[6*i+4]->Text.ToDouble() << "    ";
		*somefile << EditArray[6*i+5]->Text.ToDouble() << "    ";
		*somefile << EditArray[6*i+6]->Text.ToDouble() << "    ";
		if (ShapeComboBox->ItemIndex == 2)	*somefile << EditArray[6*i+7]->Text.ToDouble() << "    ";
		*somefile << endl;
	}
	return 0;
}

int TSpectralLineStartDialog::Write(AnsiString* LinePar)
{
	if (ShapeComboBox->ItemIndex == 0)
					*LinePar = "Fit to Lorentzian: \n";
		else if (ShapeComboBox->ItemIndex == 1)
			*LinePar = "Fit to Gaussian: \n";
				 else  *LinePar = "Fit to Pseudo-Voigt: \n" ;

	LinePar->cat_printf("    Offset : %s \n", A0Edit->Text.c_str());
	LinePar->cat_printf("Linear Term: %s \n", A1Edit->Text.c_str());
//	*somefile << "Linear term : " << A1Edit->Text.c_str() << endl;
	*LinePar += "Line          Amplitude    Center     Width    PhaseAngle(rad)  Splitting  VoigtFactor \n";

	for (int i=0 ; i<=NumberOfLinesListBox->ItemIndex ; i++)
	{
		LinePar->cat_printf("Line %d       %s    %s     %s      %s        %s       %s    \n", i+1,
		   EditArray[6*i+2]->Text.c_str(),
		EditArray[6*i+3]->Text.c_str(),
		EditArray[6*i+4]->Text.c_str(),
		EditArray[6*i+5]->Text.c_str(),
		EditArray[6*i+6]->Text.c_str(),
		EditArray[6*i+7]->Text.c_str() );
	}
	return 0;
}

void __fastcall TSpectralLineStartDialog::BaseLineSubButtonClick(TObject *Sender)
{
    TMainForm *MF;
    if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
    {
        MF->SubtractBaseline();

    }

}
//---------------------------------------------------------------------------

void __fastcall TSpectralLineStartDialog::FormActivate(TObject *Sender)
{
	NumberOfLinesListBoxChange(Sender);
}
//---------------------------------------------------------------------------


void __fastcall TSpectralLineStartDialog::ShapeComboBoxChange(TObject *Sender)
{
	NumberOfLinesListBoxChange(Sender);
}
//---------------------------------------------------------------------------




