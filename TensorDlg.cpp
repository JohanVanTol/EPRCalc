//---------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "Vector.h"
#include "Tensor.h"
#include "TensorDlg.h"
#include "Error.h"

//---------------------------------------------------------------------
#pragma resource "*.dfm"
TTensorDialog *TensorDialog;
//---------------------------------------------------------------------
__fastcall TTensorDialog::TTensorDialog(TComponent* AOwner)
	: TForm(AOwner)
{
  Tens = 0;
  for (int i=0; i<3; i++)
  {
    Val[i] = 0.0;
    Theta[i] = 0.0;
    Phi[i] = 0.0;
  }
}
//---------------------------------------------------------------------
void __fastcall TTensorDialog::FormActivate(TObject *Sender)
{
  char buffer[15];
  sprintf(buffer,"%10f", Tens->get(0,0));
  XX->Text = buffer;
  sprintf(buffer,"%10f", Tens->get(1,1));
  YY->Text = buffer;
  sprintf(buffer,"%10f", Tens->get(2,2));
  ZZ->Text = buffer;
  sprintf(buffer,"%10f", Tens->get(0,1));
  XY->Text = buffer;
  sprintf(buffer,"%10f", Tens->get(0,2));
  XZ->Text = buffer;
  sprintf(buffer,"%10f", Tens->get(1,2));
  YZ->Text = buffer;
}
//---------------------------------------------------------------------------
void TTensorDialog::UpdatePrinciples()
{
  char buffer[15];
  sprintf(buffer,"%10f", Val[0]);
  Value1->Text = buffer;
  sprintf(buffer,"%10f", Val[1]);
  Value2->Text = buffer;
  sprintf(buffer,"%10f", Val[2]);
  Value3->Text = buffer;
  sprintf(buffer,"%10f", Theta[0]);
  Theta1->Text = buffer;
  sprintf(buffer,"%10f", Theta[1]);
  Theta2->Text = buffer;
  sprintf(buffer,"%10f", Theta[2]);
  Theta3->Text = buffer;
  sprintf(buffer,"%10f", Phi[0]);
  Phi1->Text = buffer;
  sprintf(buffer,"%10f", Phi[1]);
  Phi2->Text = buffer;
  sprintf(buffer,"%10f", Phi[2]);
  Phi3->Text = buffer;
}

void __fastcall TTensorDialog::XXChange(TObject *Sender)
{
  double NewValue = (XX->Text).ToDouble();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  Tens->set(0,0,NewValue);
  if (Axial) Tens->set(1,1,NewValue);
  Tens->GetPrincipleAxes(3,Val, Theta, Phi);
  UpdatePrinciples();
}
//---------------------------------------------------------------------------
void TTensorDialog::SetTensor(Tensor *T)
{
  if (Tens == NULL) Tens = new Tensor(3);
  *Tens = *T;
  return;
}

void __fastcall TTensorDialog::XYChange(TObject *Sender)
{
  double NewValue = (XY->Text).ToDouble();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  Tens->set(0,1, NewValue);
  Tens->set(1,0, NewValue);
  Tens->GetPrincipleAxes(3, Val, Theta, Phi);
  UpdatePrinciples();
}
//---------------------------------------------------------------------------

void __fastcall TTensorDialog::XZChange(TObject *Sender)
{
  double NewValue = (XZ->Text).ToDouble();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  Tens->set(0,2,NewValue);
  Tens->set(2,0, NewValue);
  Tens->GetPrincipleAxes(3,Val, Theta, Phi);
  UpdatePrinciples();
}
//---------------------------------------------------------------------------

void __fastcall TTensorDialog::YYChange(TObject *Sender)
{
  double NewValue = (YY->Text).ToDouble();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  Tens->set(1,1,NewValue);
  Tens->GetPrincipleAxes(3,Val, Theta, Phi);
  UpdatePrinciples();
}
//---------------------------------------------------------------------------

void __fastcall TTensorDialog::ZZChange(TObject *Sender)
{
  double NewValue = (ZZ->Text).ToDouble();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  Tens->set(2,2,NewValue);
  Tens->GetPrincipleAxes(3,Val, Theta, Phi);
  UpdatePrinciples();
}
//---------------------------------------------------------------------------

void __fastcall TTensorDialog::YZChange(TObject *Sender)
{
  double NewValue = (YZ->Text).ToDouble();
  if (errno == ERANGE)
  {
    ErrorBox->Label->Caption =
        "The value you entered is not a valid number, re-enter";
    ErrorBox->ShowModal();
    return;
  }
  Tens->set(1,2,NewValue);
  Tens->set(2,1,NewValue);
  Tens->GetPrincipleAxes(3,Val, Theta, Phi);
  UpdatePrinciples();
}
//---------------------------------------------------------------------------


void __fastcall TTensorDialog::BitBtn1Click(TObject *Sender)
{
    ModalResult = mbOK;    
}
//---------------------------------------------------------------------------

void TTensorDialog::SetAxial(int ax)
{
    if (ax==1)
    {
        YY->Enabled = false;
        XY->Enabled = false;
        YZ->Enabled = false;
        XZ->Enabled = false;
        YY->Text = XX->Text;
        XY->Text = 0.0;
        YZ->Text = 0.0;
        XZ->Text = 0.0;
        XY->Visible = false;
        YZ->Visible = false;
        XZ->Visible = false;
        Axial = true;

    }
       else
       {
        YY->Enabled = true;
        XY->Enabled = true;
        YZ->Enabled = true;
        XZ->Enabled = true;
        XY->Visible = true;
        YZ->Visible = true;
        XZ->Visible = true;
        Axial = false;
       }
    return;
}
