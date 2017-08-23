//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "OptionsForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TOptionsDialog *OptionsDialog;
//---------------------------------------------------------------------------
__fastcall TOptionsDialog::TOptionsDialog(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------

