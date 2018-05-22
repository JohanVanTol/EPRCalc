//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "ExponentialDialog.h"
#include "EPRCalcMain2.h"
#include "Valid.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TExponentialDecayForm *ExponentialDecayForm;
//---------------------------------------------------------------------------
__fastcall TExponentialDecayForm::TExponentialDecayForm(TComponent* Owner)
	: TForm(Owner)
{
	mode = 0;
	ChiSqr =1;
}
//---------------------------------------------------------------------------
void __fastcall TExponentialDecayForm::SimulateButtonClick(TObject *Sender)
{
	TMainForm *MF;
	if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
	{
		MF->SimulateExponential(mode);
	}
}
//---------------------------------------------------------------------------

int TExponentialDecayForm::GetLimits(double* start, double* stop)
{
	int np;

	if (!ValidInt(NptsEdit->Text,&np)) return -1;
	if (!ValidReal(StartTimeEdit->Text,start)) return -1;
	if (!ValidReal(EndTimeEdit->Text,stop)) return -1;
	return np;
}

int TExponentialDecayForm::GetParameters(double* par, int* fixedpar)
{
	if (!ValidReal(y0Edit->Text, &par[0])) return -1;
	if (!ValidReal(Amp1Edit->Text, &par[1])) return -1;
	if (!ValidReal(Exp1Edit->Text, &par[2])) return -1;

//	par[0] = y0Edit->Text.ToDouble();
//	par[1] = Amp1Edit->Text.ToDouble();
//	par[2] = Exp1Edit->Text.ToDouble();

	if (y0CheckBox->Checked) fixedpar[0] = 1; else fixedpar[0] =0;
	if (Amp1CheckBox->Checked) fixedpar[1] = 1; else fixedpar[1] =0;
	if (Exp1CheckBox->Checked) fixedpar[2] = 1; else fixedpar[2] =0;

	int npar = 3;

	if (mode ==1)   // Biexponential
	{
		if (!ValidReal(Amp2Edit->Text, &par[3])) return -1;
		if (!ValidReal(Exp2Edit->Text, &par[4])) return -1;

		if (Amp2CheckBox->Checked) fixedpar[3] = 1; else fixedpar[3] =0;
		if (Exp2CheckBox->Checked) fixedpar[4] = 1; else fixedpar[4] =0;
		npar += 2;
	}

	if (mode == 2)   // Biexponential
	{
		if (!ValidReal(GaussAmpEdit->Text, &par[3])) return -1;
		if (!ValidReal(GaussZeroEdit->Text, &par[4])) return -1;
		if (!ValidReal(GaussSigmaEdit->Text, &par[5])) return -1;

		if (GaussAmpCheckBox->Checked) fixedpar[3] = 1; else fixedpar[3] =0;
		if (GaussZeroCheckBox->Checked) fixedpar[4] = 1; else fixedpar[4] =0;
		if (GaussSigmaCheckBox->Checked) fixedpar[5] = 1; else fixedpar[5] =0;

		npar += 3;
	}

	if (mode == 3) {
		if (!ValidReal(StretchEdit->Text, &par[3])) return -1;
		if (StretchCheckBox->Checked) fixedpar[3] = 1; else fixedpar[3] =0;

		npar+=1;
	}

	return npar;
}

int TExponentialDecayForm::SetParameters(int n, double* par, int* fixedpar)
{
	y0Edit->Text = FloatToStrF((long double)par[0], ffGeneral, 6, 0);
	Amp1Edit->Text = FloatToStrF((long double)par[1], ffGeneral, 6, 0);
	Exp1Edit->Text = FloatToStrF((long double)par[2], ffGeneral, 6, 0);

//	par[0] = y0Edit->Text.ToDouble();
//	par[1] = Amp1Edit->Text.ToDouble();
//	par[2] = Exp1Edit->Text.ToDouble();

	if (fixedpar[0]==1) y0CheckBox->Checked = true;
		   else y0CheckBox->Checked = false;
	if (fixedpar[1]==1) Amp1CheckBox->Checked = true;
		   else Amp1CheckBox->Checked = false;
	if (fixedpar[2]==1) Exp1CheckBox->Checked = true;
		   else Exp1CheckBox->Checked = false;

	if (mode == 1) {
		Amp2Edit->Text = FloatToStrF((long double)par[3], ffGeneral, 6, 0);
		Exp2Edit->Text = FloatToStrF((long double)par[4], ffGeneral, 6, 0);
	}

	if (mode == 2) {
		GaussAmpEdit->Text = FloatToStrF((long double)par[3], ffGeneral, 6, 0);
		GaussZeroEdit->Text = FloatToStrF((long double)par[4], ffGeneral, 6, 0);
		GaussSigmaEdit->Text = FloatToStrF((long double)par[5], ffGeneral, 6, 0);
	}

	if ((mode ==3)) {
		StretchEdit->Text = FloatToStrF((long double)par[3], ffGeneral, 6, 0);

	}

	return n;
}

void TExponentialDecayForm::SetErrors(int na, double* errors)
{
	y0ErrorEdit->Text = FloatToStrF((long double)errors[0], ffGeneral, 6, 0);
	Amp1ErrorEdit->Text = FloatToStrF((long double)errors[1], ffGeneral, 6, 0);
	Exp1ErrorEdit->Text = FloatToStrF((long double)errors[2], ffGeneral, 6, 0);

	if (mode == 1) {
		Amp2ErrorEdit->Text = FloatToStrF((long double)errors[3], ffGeneral, 6, 0);
		Exp2ErrorEdit->Text = FloatToStrF((long double)errors[4], ffGeneral, 6, 0);

	}

	if (mode == 2) {
		GaussAmpErrorEdit->Text = FloatToStrF((long double)errors[3], ffGeneral, 6, 0);
		GaussZeroErrorEdit->Text = FloatToStrF((long double)errors[4], ffGeneral, 6, 0);
		GaussSigmaErrorEdit->Text = FloatToStrF((long double)errors[5], ffGeneral, 6, 0);
	}

	if ((mode ==3)) {
		StretchErrorEdit->Text = FloatToStrF((long double)errors[3], ffGeneral, 6, 0);

	}

	return;
}

void __fastcall TExponentialDecayForm::CorrectBaseLineButtonClick(
	  TObject *Sender)
{
	TMainForm *MF;
	if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
	{
		MF->CorrectExponentialBaseLine();
	}
}


//---------------------------------------------------------------------------

void __fastcall TExponentialDecayForm::Fit10ButtonClick(TObject *Sender)
{
	TMainForm *MF;
	if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
	{
		switch (mode)
		{
			case 1: MF->FitBiExponentialDecay(10); break;
			case 2: MF->FitGaussExponentialDecay(10); break;
			case 3: MF->FitStretchedExponentialDecay(10); break;
			default: MF->FitExponentialDecay(10); break;
		}
	}
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

void __fastcall TExponentialDecayForm::FitButtonClick(TObject *Sender)
{
	TMainForm *MF;
	if ((MF = dynamic_cast<TMainForm *>(Application->MainForm)) != 0)
	{
		switch (mode)
		{
			case 1: MF->FitBiExponentialDecay(1); break;
			case 2: MF->FitGaussExponentialDecay(1); break;
			case 3: MF->FitStretchedExponentialDecay(1); break;
			default: MF->FitExponentialDecay(1); break;
		}
	}
}

void __fastcall TExponentialDecayForm::CancelButtonClick(TObject *Sender)
{
	Close();
}
//---------------------------------------------------------------------------

void __fastcall TExponentialDecayForm::BiRadioButtonClick(TObject *Sender)
{
	mode = 1;

	Amp2Edit->Enabled = true;
	Exp2Edit->Enabled = true;
	Amp2ErrorEdit->Enabled = true;
	Exp2ErrorEdit->Enabled = true;
	Amp2CheckBox->Enabled = true;
	Exp2CheckBox->Enabled = true;

	GaussAmpEdit->Enabled = false;
	GaussZeroEdit->Enabled = false;
	GaussSigmaEdit->Enabled = false;
	GaussAmpErrorEdit->Enabled = false;
	GaussZeroErrorEdit->Enabled = false;
	GaussSigmaErrorEdit->Enabled = false;
	GaussAmpCheckBox->Enabled = false;
	GaussZeroCheckBox->Enabled = false;
	GaussSigmaCheckBox->Enabled = false;

	StretchEdit->Enabled = false;
	StretchErrorEdit->Enabled = false;
	StretchCheckBox->Enabled = false;

}
//---------------------------------------------------------------------------



void __fastcall TExponentialDecayForm::MonoGaussRadioButtonClick(
	  TObject *Sender)
{
	mode = 2;

	Amp2Edit->Enabled = false;
	Exp2Edit->Enabled = false;
	Amp2ErrorEdit->Enabled = false;
	Exp2ErrorEdit->Enabled = false;
	Amp2CheckBox->Enabled = false;
	Exp2CheckBox->Enabled = false;

	GaussAmpEdit->Enabled = true;
	GaussZeroEdit->Enabled = true;
	GaussSigmaEdit->Enabled = true;
	GaussAmpErrorEdit->Enabled = true;
	GaussZeroErrorEdit->Enabled = true;
	GaussSigmaErrorEdit->Enabled = true;
	GaussAmpCheckBox->Enabled = true;
	GaussZeroCheckBox->Enabled = true;
	GaussSigmaCheckBox->Enabled = true;

	StretchEdit->Enabled = false;
	StretchErrorEdit->Enabled = false;
	StretchCheckBox->Enabled = false;


}
//---------------------------------------------------------------------------

void __fastcall TExponentialDecayForm::MonoRadioButtonClick(TObject *Sender)
{
	mode = 0;

	Amp2Edit->Enabled = false;
	Exp2Edit->Enabled = false;
	Amp2ErrorEdit->Enabled = false;
	Exp2ErrorEdit->Enabled = false;
	Amp2CheckBox->Enabled = false;
	Exp2CheckBox->Enabled = false;

	GaussAmpEdit->Enabled = false;
	GaussZeroEdit->Enabled = false;
	GaussSigmaEdit->Enabled = false;
	GaussAmpErrorEdit->Enabled = false;
	GaussZeroErrorEdit->Enabled = false;
	GaussSigmaErrorEdit->Enabled = false;
	GaussAmpCheckBox->Enabled = false;
	GaussZeroCheckBox->Enabled = false;
	GaussSigmaCheckBox->Enabled = false;

	StretchEdit->Enabled = false;
	StretchErrorEdit->Enabled = false;
	StretchCheckBox->Enabled = false;

}
//---------------------------------------------------------------------------

void __fastcall TExponentialDecayForm::OKButtonClick(TObject *Sender)
{
	Close();
}
//---------------------------------------------------------------------------

int TExponentialDecayForm::Write(ofstream *somefile)
{
	if (mode == 0)
	{
		*somefile << "Fit to Single Exponential Decay: " << endl;
	}
	else if (mode == 1)
		{
			*somefile << "Fit to Biexponential Decay: " << endl;
		}
		 else
		 {
			*somefile << "Fit to Gaussian + Mono exponential decay: " << endl;
		 }

	*somefile << "Baseline   : " << y0Edit->Text.c_str()
		<< "  +-  " << y0ErrorEdit->Text.c_str() << endl;
	*somefile << "Amplitude1 : " << Amp1Edit->Text.c_str()
		<< "  +-  " << Amp1ErrorEdit->Text.c_str() << endl;
	*somefile << "Exponent1  : " << Exp1Edit->Text.c_str()
		<< "  +-  " << Exp1ErrorEdit->Text.c_str() << endl;

	if (mode == 1)
	{
		*somefile << "Amplitude2 : " << Amp2Edit->Text.c_str()
				<< "  +-  " << Amp2ErrorEdit->Text.c_str() << endl;
		*somefile << "Exponent2  : " << Exp2Edit->Text.c_str()
				<< "  +-  " << Exp2ErrorEdit->Text.c_str() << endl;
	}

	if (mode == 2)
	{
		*somefile << "Gauss Amplit : " << GaussAmpEdit->Text.c_str()
				<< "  +-  " << GaussAmpErrorEdit->Text.c_str() << endl;
		*somefile << "Gauss Center : " << GaussZeroEdit->Text.c_str()
				<< "  +-  " << GaussZeroErrorEdit->Text.c_str() << endl;
		*somefile << "Gauss Width  : " << GaussSigmaEdit->Text.c_str()
				<< "  +-  " << GaussSigmaErrorEdit->Text.c_str() << endl;
	}

	if (mode == 3)
	{
		*somefile << "Stretch Parameter : " << StretchEdit->Text.c_str()
				<< "  +-  " << StretchErrorEdit->Text.c_str() << endl;
	}

	*somefile << "ChiSqr : " << ChiSqr << endl;

	return 0;
}

int TExponentialDecayForm::Write(AnsiString *ExParString)
{
	if (mode == 0)
	{
		*ExParString = "Fit to Single Exponential Decay: \n";
	}
	else if (mode == 1)
		{
			*ExParString = "Fit to Biexponential Decay: \n";
		}
		 else if (mode == 2)
		 {
			 *ExParString = "Fit to Gaussian + Mono exponential decay: ";
		 }
		   else
		 {
			*ExParString = "Fit to Stretched Exponential decay: ";
		 }

	ExParString->cat_printf("Baseline   : %s  +- %s \n", y0Edit->Text.c_str(),
		y0ErrorEdit->Text.c_str() );
	ExParString->cat_printf("Amplitude1 : %s  +- %s \n", Amp1Edit->Text.c_str(),
			Amp1ErrorEdit->Text.c_str());
	ExParString->cat_printf("Exponent1  : %s  +- %s \n", Exp1Edit->Text.c_str(),
			Exp1ErrorEdit->Text.c_str());

	if (mode == 1)
	{
		ExParString->cat_printf("Amplitude2 : %s  +- %s \n",
			Amp2Edit->Text.c_str(), Amp2ErrorEdit->Text.c_str());
		ExParString->cat_printf("Exponent2  : %s  +- %s \n",
			Exp2Edit->Text.c_str(), Exp2ErrorEdit->Text.c_str());
	}

	if (mode == 2)
	{
		ExParString->cat_printf("Gauss Amplit : %s  +- %s \n",
			GaussAmpEdit->Text.c_str(), GaussAmpErrorEdit->Text.c_str() );
		ExParString->cat_printf("Gauss Center : %s  +- %s \n",
			GaussZeroEdit->Text.c_str(), GaussZeroErrorEdit->Text.c_str() );
		ExParString->cat_printf("Gauss Width  : %s  +- %s \n",
			GaussSigmaEdit->Text.c_str(), GaussSigmaErrorEdit->Text.c_str() );
	}

	if (mode == 3)
	{
		ExParString->cat_printf("Stretch parameter : %s  +- %s \n",
			StretchEdit->Text.c_str(), StretchErrorEdit->Text.c_str());
	}

	ExParString->cat_printf("ChiSqr: %lf \n", getChiSqr());

	return ExParString->Length();
}

void __fastcall TExponentialDecayForm::StretchedRadioButtonClick(
	  TObject *Sender)
{
	mode = 3;

	Amp2Edit->Enabled = false;
	Exp2Edit->Enabled = false;
	Amp2ErrorEdit->Enabled = false;
	Exp2ErrorEdit->Enabled = false;
	Amp2CheckBox->Enabled = false;
	Exp2CheckBox->Enabled = false;

	GaussAmpEdit->Enabled = false;
	GaussZeroEdit->Enabled = false;
	GaussSigmaEdit->Enabled = false;
	GaussAmpErrorEdit->Enabled = false;
	GaussZeroErrorEdit->Enabled = false;
	GaussSigmaErrorEdit->Enabled = false;
	GaussAmpCheckBox->Enabled = false;
	GaussZeroCheckBox->Enabled = false;
	GaussSigmaCheckBox->Enabled = false;

	StretchEdit->Enabled = true;
	StretchErrorEdit->Enabled = true;
	StretchCheckBox->Enabled = true;


}
//---------------------------------------------------------------------------


