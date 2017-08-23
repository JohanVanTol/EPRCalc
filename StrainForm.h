//---------------------------------------------------------------------------

#ifndef StrainFormH
#define StrainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TStrainDialog : public TForm
{
__published:	// IDE-managed Components
    TGroupBox *GroupBox1;
    TBitBtn *OKButton;
    TBitBtn *CancelButton;
    TEdit *B20StrainEdit;
    TCheckBox *B20StrainCheckBox;
    TLabel *Label1;
    TCheckBox *B22StrainCheckBox;
    TEdit *B22StrainEdit;
    TLabel *Label2;
    void __fastcall OKButtonClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
    __fastcall TStrainDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TStrainDialog *StrainDialog;
//---------------------------------------------------------------------------
#endif
