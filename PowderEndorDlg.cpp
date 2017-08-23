//---------------------------------------------------------------------------
#include <vcl.h>
#include <iostream.h>
#pragma hdrstop

#include "PowderEndorDlg.h"
#include "ENDORsimParam.h"
#include "Error.h"
#include "Valid.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TPowderENDORDialog *PowderENDORDialog;
//---------------------------------------------------------------------------
__fastcall TPowderENDORDialog::TPowderENDORDialog(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------








void __fastcall TPowderENDORDialog::FormActivate(TObject *Sender)
{
    UpdateParameters();
}
//---------------------------------------------------------------------------
void __fastcall TPowderENDORDialog::OkButtonClick(TObject *Sender)
{
    ModalResult = mbOK;
}
//---------------------------------------------------------------------------

int TPowderENDORDialog::UpdateParameters()
{
    NptsEdit->Text = Par.GetNpts();
    FirstEdit->Text = Par.GetFirst();
    LastEdit->Text = Par.GetLast();
    AngleStepEdit->Text = Par.GetAngleStep();

    FreqEdit->Text = Par.GetFreq();
    NgEdit->Text = Par.GetNg();
    gStartEdit->Text = Par.GetStartg();
    gStepEdit->Text = Par.GetStepg();

    ModulationEdit->Text = Par.GetModulation();
    ENDORwidthEdit->Text = Par.GetENDORwidth();

    WxEdit->Text = Par.GetEPRWidth(0);
    WyEdit->Text = Par.GetEPRWidth(1);
    WzEdit->Text = Par.GetEPRWidth(2);

    ANxEdit->Text = Par.GetAN(0);
    ANyEdit->Text = Par.GetAN(1);
    ANzEdit->Text = Par.GetAN(2);

    int mult = Par.GetMultN();
    if ((mult >0) && (mult < 8))
    SpinMultComboBox->ItemIndex = mult-1;
    if (mult>1)
    {
        ANxEdit->Enabled = true;
        ANyEdit->Enabled = true;
        ANzEdit->Enabled = true;
    }
      else
        {
         ANxEdit->Enabled = false;
         ANyEdit->Enabled = false;
         ANzEdit->Enabled = false;
        }

    return 0;

}



void __fastcall TPowderENDORDialog::SpinMultComboBoxChange(TObject *Sender)
{
    Par.SetMultN(SpinMultComboBox->ItemIndex +1);
    UpdateParameters();
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::CanvelButtonClick(TObject *Sender)
{
    ModalResult = mbCancel;
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::LoadButtonClick(TObject *Sender)
{
  if (SimOpenDialog->Execute())
  {
    Par.Read(SimOpenDialog->FileName.c_str());
    UpdateParameters();
  }

}
//---------------------------------------------------------------------------



void __fastcall TPowderENDORDialog::NptsEditExit(TObject *Sender)
{
  int MaxN = 4096;
  int N;
  if (!ValidInt(NptsEdit->Text, &N))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  if (N > MaxN)
  {
    ErrorBox->Label->Caption =
        "The value you entered is too large";
    ErrorBox->ShowModal();
    N = MaxN;
  }
  Par.SetNpts(N);
    return;
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::FirstEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(FirstEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetFirst(R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::LastEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(LastEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetLast(R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::AngleStepEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(AngleStepEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetAngleStep(R);
}
//---------------------------------------------------------------------------


void __fastcall TPowderENDORDialog::ENDORwidthEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(ENDORwidthEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetENDORwidth(R);

}
//---------------------------------------------------------------------------


void __fastcall TPowderENDORDialog::ModulationEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(ModulationEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetModulation(R);

}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::FreqEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(FreqEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetFreq(R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::NgEditExit(TObject *Sender)
{
  int N;
  if (!ValidInt(NgEdit->Text, &N))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetNg(N);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::gStartEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(gStartEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetStartg(R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::gStepEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(gStepEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetStepg(R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::WxEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(WxEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetEPRWidth(0,R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::WyEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(WyEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetEPRWidth(1,R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::WzEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(WzEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetEPRWidth(2,R);
}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::ANxEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(ANxEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetAN(0,R);

}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::ANyEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(ANyEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetAN(1,R);

}
//---------------------------------------------------------------------------

void __fastcall TPowderENDORDialog::ANzEditExit(TObject *Sender)
{
  double R;
  if (!ValidReal(ANzEdit->Text, &R))
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    UpdateParameters();
    return;
  }
  Par.SetAN(2,R);

}
//---------------------------------------------------------------------------

