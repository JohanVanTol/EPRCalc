//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "StartLineDialog.h"
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
    EditArray = new TEdit*[33];
    CheckArray = new TCheckBox*[33];
    EditArray[0] = A0Edit;
    EditArray[1] = A1Edit;
    EditArray[2] = A2Edit;
    EditArray[3] = A3Edit;
    EditArray[4] = A4Edit;
    EditArray[5] = A5Edit;
    EditArray[6] = A6Edit;
    EditArray[7] = A7Edit;
    EditArray[8] = A8Edit;
    EditArray[9] = A9Edit;
    EditArray[10] = A10Edit;
    EditArray[11] = A11Edit;
    EditArray[12] = A12Edit;
    EditArray[13] = A13Edit;
    EditArray[14] = A14Edit;
    EditArray[15] = A15Edit;
    EditArray[16] = A16Edit;
    EditArray[17] = A17Edit;
    EditArray[18] = A18Edit;
    EditArray[19] = A19Edit;
    EditArray[20] = A20Edit;
    EditArray[21] = A21Edit;
    EditArray[22] = A22Edit;
    EditArray[23] = A23Edit;
    EditArray[24] = A24Edit;
    EditArray[25] = A25Edit;
    EditArray[26] = A26Edit;
    EditArray[27] = A27Edit;
    EditArray[28] = A28Edit;
    EditArray[29] = A29Edit;
    EditArray[30] = A30Edit;
    EditArray[31] = A31Edit;
    CheckArray[0] = A0CheckBox;
    CheckArray[1] = A1CheckBox;
    CheckArray[2] = A2CheckBox;
    CheckArray[3] = A3CheckBox;
    CheckArray[4] = A4CheckBox;
    CheckArray[5] = A5CheckBox;
    CheckArray[6] = A6CheckBox;
    CheckArray[7] = A7CheckBox;
    CheckArray[8] = A8CheckBox;
    CheckArray[9] = A9CheckBox;
    CheckArray[10] = A10CheckBox;
    CheckArray[11] = A11CheckBox;
    CheckArray[12] = A12CheckBox;
    CheckArray[13] = A13CheckBox;
    CheckArray[14] = A14CheckBox;
    CheckArray[15] = A15CheckBox;
    CheckArray[16] = A16CheckBox;
    CheckArray[17] = A17CheckBox;
    CheckArray[18] = A18CheckBox;
    CheckArray[19] = A19CheckBox;
    CheckArray[20] = A20CheckBox;
    CheckArray[21] = A21CheckBox;
    CheckArray[22] = A22CheckBox;
    CheckArray[23] = A23CheckBox;
    CheckArray[24] = A24CheckBox;
    CheckArray[25] = A25CheckBox;
    CheckArray[26] = A26CheckBox;
    CheckArray[27] = A27CheckBox;
    CheckArray[28] = A28CheckBox;
    CheckArray[29] = A29CheckBox;
    CheckArray[30] = A30CheckBox;
    CheckArray[31] = A31CheckBox;
    par = new double[33];
    fixed = new int[33];
}
//---------------------------------------------------------------------------
void TSpectralLineStartDialog::SetParameters(int n, double* param, int* fixedpar)
{
    if (n>32) return;
    npar = n;
    for (int i=0; i<npar;i++)
    {
        par[i] = param[i];
        EditArray[i]->Text = param[i];
        fixed[i] = fixedpar[i];
        if (fixed[i]==1) CheckArray[i]->Checked = true;
           else CheckArray[i]->Checked = false;
    }
}

int TSpectralLineStartDialog::GetParameters(double* param, int* fixedpar)
{
    npar = (NumberOfLinesListBox->ItemIndex+1)*3 + 2;
    for (int i = 0; i< npar; i++)
    {
        par[i] = EditArray[i]->Text.ToDouble();
        param[i] = par[i];
        if (CheckArray[i]->Checked) fixed[i] = 1;
          else fixed[i] = 0;
        fixedpar[i] = fixed[i];
    }
    return npar;
}


int TSpectralLineStartDialog::UpdateForm()
{

}



void __fastcall TSpectralLineStartDialog::NumberOfLinesListBoxChange(
      TObject *Sender)
{
    for (int i=2; i<32;i++)
        EditArray[i]->Enabled = false;
    for (int k=0; k<=NumberOfLinesListBox->ItemIndex; k++)
    {
        EditArray[3*(k+1)-1]->Enabled = true;
        EditArray[3*(k+1)]->Enabled = true;
        EditArray[3*(k+1)+1]->Enabled = true;
    }
}
//---------------------------------------------------------------------------


void __fastcall TSpectralLineStartDialog::FormDestroy(TObject *Sender)
{
    delete[] EditArray;
    delete[] CheckArray;
    delete[] par;
    delete[] fixed;
}
//---------------------------------------------------------------------------


void __fastcall TSpectralLineStartDialog::BitBtn1Click(TObject *Sender)
{
    ModalResult = mbOK;
}
//---------------------------------------------------------------------------


