//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <iostream.h>
#include <fstream.h>
#include "Hamilton.h"
#include "Vector.h"
#include "Tensor.h"
#include "TensorDlg.h"
#include "HamilParam.h"
#include "Error.h"
#include "Valid.h"
#include "StrainForm.h"
#include "gStrainForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
THamiltonianDialog *HamiltonianDialog;
//---------------------------------------------------------------------------
__fastcall THamiltonianDialog::THamiltonianDialog(TComponent* Owner)
  : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::GTensorButtonClick(TObject *Sender)
{
 TensorDialog->SetTensor(&(Par.Getg()));
 if (AxialCheckBox->Checked) TensorDialog->SetAxial(1);
    else TensorDialog->SetAxial(0);
 TensorDialog->ShowModal();
 if (TensorDialog->ModalResult == mbOK)
    Par.Setg(TensorDialog->GetTensor());
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::LoadHamilButtonClick(TObject *Sender)
{
  if (HamiltonOpenDialog->Execute())
  {
    Par.Read(HamiltonOpenDialog->FileName.c_str());
    UpdateParameters();
  }
}
//---------------------------------------------------------------------------


void __fastcall THamiltonianDialog::SaveHamilButtonClick(TObject *Sender)
{
    if (HamiltonSaveDialog->Execute())
    {
     Par.Write(HamiltonSaveDialog->FileName.c_str());
    }
}
//---------------------------------------------------------------------------


void __fastcall THamiltonianDialog::FormActivate(TObject *Sender)
{
    UpdateParameters();
}
//---------------------------------------------------------------------------

void THamiltonianDialog::UpdateParameters()
{
    int Sm = Par.GetMultS();
    if (Sm < 2) Sm = 2;
    if (Sm > 2)
    {
        B20Edit->Enabled = true;
        B22Edit->Enabled=true;
    }
        else
        {
            B20Edit->Enabled = false;
            B22Edit->Enabled = false;
        }
    if (Sm > 4)
    {
        B40Edit->Enabled = true;
        B42Edit->Enabled=true;
        B43Edit->Enabled = true;
        B44Edit->Enabled=true;
    }
        else
    {
        B40Edit->Enabled = false;
        B42Edit->Enabled= false;
        B43Edit->Enabled = false;
        B44Edit->Enabled= false;
    }
    if (Sm > 6)
    {
        B60Edit->Enabled = true;
        B63Edit->Enabled=true;
        B64Edit->Enabled = true;
        B66Edit->Enabled=true;
    }
        else
    {
        B60Edit->Enabled = false;
        B63Edit->Enabled= false;
        B64Edit->Enabled = false;
        B66Edit->Enabled= false;
    }
    SpinComboBox->ItemIndex=Sm -2;
    B20Edit->Text = Par.GetCF(0);
    B22Edit->Text = Par.GetCF(1);
    B40Edit->Text = Par.GetCF(2);
    B42Edit->Text = Par.GetCF(3);
    B43Edit->Text = Par.GetCF(4);
    B44Edit->Text = Par.GetCF(5);
    B60Edit->Text = Par.GetCF(6);
    B63Edit->Text = Par.GetCF(7);
    B64Edit->Text = Par.GetCF(8);
    B66Edit->Text = Par.GetCF(9);
    int nn= Par.GetnI();
    NNucSpinEdit->Text = nn;
    if (nn > 0)
    {
        SpinNumberBox->Enabled = true;
        NucSpinComboBox->Enabled = true;
        NucZeemanEdit->Enabled = true;
        HyperfineButton->Enabled = true;
        int Nindex = SpinNumberBox->ItemIndex;
        if (Nindex < 0) Nindex = 0;
        if (Nindex >= nn) Nindex = nn-1;
        SpinNumberBox->Items->Clear();
        AnsiString AS;
        for (int i=1; i<=nn; i++)
        {
            AS = i;
            SpinNumberBox->Items->Add(AS);
        }
        SpinNumberBox->ItemIndex = Nindex;
        NucSpinComboBox->ItemIndex = Par.GetMultI(Nindex)-2;
        if (NucSpinComboBox->ItemIndex >0)
            QuadrupoleButton->Enabled = true;
          else QuadrupoleButton->Enabled = false;
        NucZeemanEdit->Text = Par.Getgamma(Nindex);
    }
        else
    {
        SpinNumberBox->Enabled = false;
        NucSpinComboBox->Enabled = false;
        NucZeemanEdit->Enabled = false;
        HyperfineButton->Enabled = false;
        QuadrupoleButton->Enabled = false;
        SpinNumberBox->Items->Clear();
        SpinNumberBox->ItemIndex = -1;
    }
    StrainDialog->B20StrainEdit->Text = Par.GetCFStrain(0);
    StrainDialog->B22StrainEdit->Text = Par.GetCFStrain(1);
}

void __fastcall THamiltonianDialog::SpinComboBoxChange(TObject *Sender)
{
	int ElecMult = SpinComboBox->ItemIndex + 2;
	if (ElecMult >= 42) ElecMult = 66;
	Par.SetMultS(ElecMult);
    UpdateParameters();
}
//---------------------------------------------------------------------------


void __fastcall THamiltonianDialog::NNucSpinEditExit(TObject *Sender)
{
  int N =(NNucSpinEdit->Text).ToInt();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  if (N > 20) N = 20;
  if (N < 0) N = 0;
  NNucSpinEdit->Text = N;
  Par.SetnI(N);
  UpdateParameters();

}
//---------------------------------------------------------------------------


void __fastcall THamiltonianDialog::AddButtonClick(TObject *Sender)
{
    Par.AddNuc();
    UpdateParameters();    
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::HyperfineButtonClick(TObject *Sender)
{
    int Nuc = SpinNumberBox->ItemIndex;
    TensorDialog->SetTensor(&(Par.GetA(Nuc)));
     if (AxialCheckBox->Checked) TensorDialog->SetAxial(1);
            else TensorDialog->SetAxial(0);
    TensorDialog->ShowModal();
    if (TensorDialog->ModalResult == mbOK)
    Par.SetA(Nuc,TensorDialog->GetTensor());

}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::QuadrupoleButtonClick(TObject *Sender)
{
    int Nuc = SpinNumberBox->ItemIndex;
    TensorDialog->SetTensor(&(Par.GetQ(Nuc)));
    TensorDialog->ShowModal();
    if (TensorDialog->ModalResult == mbOK)
    Par.SetQ(Nuc, TensorDialog->GetTensor());
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::SpinNumberBoxChange(TObject *Sender)
{
    UpdateParameters();
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::NucSpinComboBoxChange(TObject *Sender)
{
    Par.SetMultI(SpinNumberBox->ItemIndex,NucSpinComboBox->ItemIndex + 2);
    UpdateParameters();
}
//---------------------------------------------------------------------------


void __fastcall THamiltonianDialog::OKButtonClick(TObject *Sender)
{
    ModalResult = mbOK;
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::DeleteButtonClick(TObject *Sender)
{
    Par.DeleteNuc(SpinNumberBox->ItemIndex);
    UpdateParameters();
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::NucZeemanEditExit(TObject *Sender)
{
  double NewValue;
  if (ValidReal(NucZeemanEdit->Text, &NewValue))
  {
    Par.Setgamma(SpinNumberBox->ItemIndex, NewValue);
    UpdateParameters();
  }
    else
    {
      ErrorBox->Label->Caption =
       "The value you entered is not a valid number, re-enter";
      ErrorBox->ShowModal();
    }
  return;
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B20EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B20Edit->Text, &B))
        Par.SetCF(0,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B22EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B22Edit->Text, &B))
        Par.SetCF(1,B);
}
//---------------------------------------------------------------------------



void __fastcall THamiltonianDialog::B40EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B40Edit->Text, &B))
        Par.SetCF(2,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B42EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B42Edit->Text, &B))
        Par.SetCF(3,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B43EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B43Edit->Text, &B))
        Par.SetCF(4,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B44EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B44Edit->Text, &B))
        Par.SetCF(5,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B60EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B60Edit->Text, &B))
        Par.SetCF(6,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B63EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B63Edit->Text, &B))
        Par.SetCF(7,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B64EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B64Edit->Text, &B))
        Par.SetCF(8,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::B66EditExit(TObject *Sender)
{
    double B;
    if (ValidReal(B66Edit->Text, &B))
        Par.SetCF(9,B);
}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::StrainButtonClick(TObject *Sender)
{
	double Strain;
	StrainDialog->ShowModal();
	if (StrainDialog->ModalResult == mbOK)
	{
		if (!ValidReal(StrainDialog->B20StrainEdit->Text, &Strain)) Strain = 0.0;
		Par.SetCFStrain(0,Strain);
		if (!ValidReal(StrainDialog->B22StrainEdit->Text, &Strain)) Strain = 0.0;
		Par.SetCFStrain(1,Strain);
	}
	return;

}
//---------------------------------------------------------------------------

void __fastcall THamiltonianDialog::gStrainButtonClick(TObject *Sender)
{
	double gStrain;
	gStrainDialog->ShowModal();
	if (gStrainDialog->ModalResult == mbOK)
	{
		if (!ValidReal(gStrainDialog->gStrainEdit->Text, &gStrain)) gStrain = 0.0 ;
		Par.SetgStrain(gStrain);
	}
	return;

}
//---------------------------------------------------------------------------

