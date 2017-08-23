//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FitCycleDlg.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TFitCycleDialog *FitCycleDialog;
//---------------------------------------------------------------------------
__fastcall TFitCycleDialog::TFitCycleDialog(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TFitCycleDialog::BitBtn1Click(TObject *Sender)
{
    ModalResult = mbCancel;    
}
//---------------------------------------------------------------------------
void __fastcall TFitCycleDialog::BitBtn2Click(TObject *Sender)
{
    ModalResult = mbOK;    
}
//---------------------------------------------------------------------------




