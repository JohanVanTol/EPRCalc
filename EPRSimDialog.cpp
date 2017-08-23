//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "EPRSimDialog.h"
#include "Valid.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TEPRSimulationDialog *EPRSimulationDialog;
//---------------------------------------------------------------------------
__fastcall TEPRSimulationDialog::TEPRSimulationDialog(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::PxEditChange(TObject *Sender)
{
//    if (PxEdit->Text.ToDouble()>1.0) PxEdit->Text = "1.0";
//    if (PxEdit->Text.ToDouble()<0.0) PxEdit->Text = "0.0";
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::BitBtn4Click(TObject *Sender)
{
    NumberChange = true;
    ModalResult = mbOK;
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::PopRatioEditChange(TObject *Sender)
{
    PopChange = true;
    double Ratio, Px;
    double Pz;
    if (ValidReal(PzEdit->Text, &Pz)) {}
        else Pz=1.0;
    if (ValidReal(PopRatioEdit->Text, &Ratio))
    {
        if (Ratio >=0)
        {
             PyEdit->Text = Pz+1;
             Px = Ratio + Pz;
             PxEdit->Text = Px;
        }
          else
          {
             PyEdit->Text = "0.0";
             Px = Ratio*Pz;
             PxEdit->Text = Px;
          }
    }

}
//---------------------------------------------------------------------------




void __fastcall TEPRSimulationDialog::AngleStepEditChange(TObject *Sender)
{
    NumberChange = true;    
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::TempEditChange(TObject *Sender)
{
    PopChange = true;    
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::FreqEditChange(TObject *Sender)
{
    HamilChange = true;    
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::LwxEditChange(TObject *Sender)
{
    WidthChange = true;    
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::LwyEditChange(TObject *Sender)
{
    WidthChange = true;
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::LwzEditChange(TObject *Sender)
{
    WidthChange = true;

}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::GwxEditChange(TObject *Sender)
{
    WidthChange = true;

}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::GwyEditChange(TObject *Sender)
{
    WidthChange = true;
    
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::GwzEditChange(TObject *Sender)
{
    WidthChange = true;
    
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::DispersEditChange(TObject *Sender)
{
    DisPersChange = true;
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::FormCreate(TObject *Sender)
{
    NumberChange = true;
    PopChange = true;
    DisPersChange = true;
    WidthChange = true;
    HamilChange = true;
    CalcOrderComboBox->ItemIndex = 0;
}
//---------------------------------------------------------------------------


void __fastcall TEPRSimulationDialog::CalcOrderComboBoxChange(TObject *Sender)
{
    HamilChange = true;
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::BitBtn5Click(TObject *Sender)
{
    ModalResult = mbCancel;
}
//---------------------------------------------------------------------------
int TEPRSimulationDialog::Write(ofstream *somefile)
{
    double dbl;
    int num;

    *somefile << "Sim Parameters " << endl;
    *somefile << "Field: " ;
    if (ValidReal(StartFieldEdit->Text, &dbl)) *somefile << dbl;
    *somefile << "T to " ;
    if (ValidReal(EndFieldEdit->Text, &dbl)) *somefile << dbl;
    *somefile << "T, npts:  " ;
    if (ValidInt(NptsEdit->Text, &num)) *somefile << num;

    *somefile << endl << "Anglesteps : ";
    if (ValidReal(AngleStepEdit->Text, &dbl)) *somefile << dbl;
    *somefile << "degr,  calculation to order " << CalcOrderComboBox->ItemIndex + 1;

	*somefile << endl << "Frequency: " ;
	if (ValidReal(FreqEdit->Text, &dbl)) *somefile << dbl;
	*somefile << " GHz, at " ;
	if (ValidReal(TempEdit->Text, &dbl)) *somefile << dbl;
	*somefile << " K. " << endl;

	if (ShapeButtonGroup->ItemIndex == 1)
		*somefile << "Lorentz Line width (x,y,z): " ;
	  else {
		if (ShapeButtonGroup->ItemIndex == 0)
			*somefile << "Gaussian Line width (x,y,z): " ;
		  else {
			*somefile << "PseudoVoight Line width " ;
			if (ValidReal(MixEdit->Text, &dbl)) *somefile << dbl;
			*somefile << "percent Lorentz, with widths (x,y,z): " ;
		  }
	  }

	if (ValidReal(LwxEdit->Text, &dbl)) *somefile << dbl;
	*somefile << ", " ;
	if (ValidReal(LwyEdit->Text, &dbl)) *somefile << dbl;
	*somefile << ", " ;
	if (ValidReal(LwzEdit->Text, &dbl)) *somefile << dbl;
	*somefile << endl ;

	*somefile << "Dispersion/Absorbtion (degrees): " ;
	if (ValidReal(DispersEdit->Text, &dbl)) *somefile << dbl;

	if (SOISCcheckbox->Checked)
	{
		*somefile << endl << "Population Ratio (Px-Pz)/(Py-Pz): " ;
		if (ValidReal(PopRatioEdit->Text, &dbl)) *somefile << dbl;
	}
	  else
		{
			*somefile << endl << "relative ms=0 population rates: Px/Pz " ;
			if (ValidReal(XZratioEdit->Text, &dbl)) *somefile << dbl;
			*somefile << "  Py/Pz " ;
			if (ValidReal(YZratioEdit->Text, &dbl)) *somefile << dbl;
		 }
	*somefile << endl << "Decay rates x, y, z  ";
	if (ValidReal(kxEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(kyEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(kzEdit->Text, &dbl)) *somefile << dbl << endl;

	*somefile << "SLR rates x, y, z  ";
	if (ValidReal(SLRxEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(SLRyEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(SLRzEdit->Text, &dbl)) *somefile << dbl << endl;

	return 0;
}

int TEPRSimulationDialog::Write(AnsiString* SimParString)
{
	double dbl;
	int num;

	*SimParString = "Sim Parameters \n";
 /*	*somefile << "Field: " ;
	if (ValidReal(StartFieldEdit->Text, &dbl)) *somefile << dbl;
	*somefile << "T to " ;
	if (ValidReal(EndFieldEdit->Text, &dbl)) *somefile << dbl;
	*somefile << "T, npts:  " ;
	if (ValidInt(NptsEdit->Text, &num)) *somefile << num;

	*somefile << endl << "Anglesteps : ";
	if (ValidReal(AngleStepEdit->Text, &dbl)) *somefile << dbl;
	*somefile << "degr,  calculation to order " << CalcOrderComboBox->ItemIndex + 1;

	*somefile << endl << "Frequency: " ;
	if (ValidReal(FreqEdit->Text, &dbl)) *somefile << dbl;
	*somefile << " GHz, at " ;
	if (ValidReal(TempEdit->Text, &dbl)) *somefile << dbl;
	*somefile << " K. " << endl;

	if (ShapeButtonGroup->ItemIndex == 1)
		*somefile << "Lorentz Line width (x,y,z): " ;
	  else {
		if (ShapeButtonGroup->ItemIndex == 0)
			*somefile << "Gaussian Line width (x,y,z): " ;
		  else {
			*somefile << "PseudoVoight Line width " ;
			if (ValidReal(MixEdit->Text, &dbl)) *somefile << dbl;
			*somefile << "percent Lorentz, with widths (x,y,z): " ;
		  }
	  }

	if (ValidReal(LwxEdit->Text, &dbl)) *somefile << dbl;
	*somefile << ", " ;
	if (ValidReal(LwyEdit->Text, &dbl)) *somefile << dbl;
	*somefile << ", " ;
	if (ValidReal(LwzEdit->Text, &dbl)) *somefile << dbl;
	*somefile << endl ;

	*somefile << "Dispersion/Absorbtion (degrees): " ;
	if (ValidReal(DispersEdit->Text, &dbl)) *somefile << dbl;

	if (SOISCcheckbox->Checked)
	{
		*somefile << endl << "Population Ratio (Px-Pz)/(Py-Pz): " ;
		if (ValidReal(PopRatioEdit->Text, &dbl)) *somefile << dbl;
	}
	  else
		{
			*somefile << endl << "relative ms=0 population rates: Px/Pz " ;
			if (ValidReal(XZratioEdit->Text, &dbl)) *somefile << dbl;
			*somefile << "  Py/Pz " ;
			if (ValidReal(YZratioEdit->Text, &dbl)) *somefile << dbl;
		 }
	*somefile << endl << "Decay rates x, y, z  ";
	if (ValidReal(kxEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(kyEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(kzEdit->Text, &dbl)) *somefile << dbl << endl;

	*somefile << "SLR rates x, y, z  ";
	if (ValidReal(SLRxEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(SLRyEdit->Text, &dbl)) *somefile << dbl << "  ";
	if (ValidReal(SLRzEdit->Text, &dbl)) *somefile << dbl << endl;
*/
	return 0;
}

void __fastcall TEPRSimulationDialog::XZratioEdit1Change(TObject *Sender)
{
	PopChange = true;
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::XZratioEditChange(TObject *Sender)
{
    PopChange = true;
    double Ratio, Px;
    if (ValidReal(XZratioEdit->Text, &Ratio))
    {
        if (Ratio >=0)
        {
             PzEdit->Text = "1.0";
             Px = Ratio ;
             PxEdit->Text = Px;
        }
    }
}
//---------------------------------------------------------------------------

void __fastcall TEPRSimulationDialog::YZratioEditChange(TObject *Sender)
{
    PopChange = true;
    double Ratio, Py;
    if (ValidReal(YZratioEdit->Text, &Ratio))
    {
        if (Ratio >=0)
        {
             PzEdit->Text = "1.0";
             Py = Ratio ;
             PyEdit->Text = Py;
        }
    }

}
//---------------------------------------------------------------------------





void __fastcall TEPRSimulationDialog::ShapeButtonGroupClick(TObject *Sender)
{
	if (ShapeButtonGroup->ItemIndex==2)
	{
		MixEdit->Enabled = true;
		VoigtAlphaLabel->Enabled = true;

	}
	  else
	  {
		MixEdit->Enabled = false;
		VoigtAlphaLabel->Enabled = false;
	  }
}
//---------------------------------------------------------------------------


