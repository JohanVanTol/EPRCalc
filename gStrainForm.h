//---------------------------------------------------------------------------

#ifndef gStrainFormH
#define gStrainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TgStrainDialog : public TForm
{
__published:	// IDE-managed Components
	TCheckBox *gStrainCheckBox;
	TLabel *label1;
	TEdit *gStrainEdit;
	TLabel *Label2;
	TBitBtn *BitBtn1;
	void __fastcall BitBtn1Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TgStrainDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TgStrainDialog *gStrainDialog;
//---------------------------------------------------------------------------
#endif
