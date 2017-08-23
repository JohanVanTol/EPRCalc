//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "StrainForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TStrainDialog *StrainDialog;
//---------------------------------------------------------------------------
__fastcall TStrainDialog::TStrainDialog(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TStrainDialog::OKButtonClick(TObject *Sender)
{
    ModalResult = mbOK;    
}
//---------------------------------------------------------------------------
