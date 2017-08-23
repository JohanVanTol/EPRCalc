//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "gStrainForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TgStrainDialog *gStrainDialog;
//---------------------------------------------------------------------------
__fastcall TgStrainDialog::TgStrainDialog(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TgStrainDialog::BitBtn1Click(TObject *Sender)
{
	ModalResult = mbOK;	
}
//---------------------------------------------------------------------------
