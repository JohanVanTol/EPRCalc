//---------------------------------------------------------------------------

#ifndef OptionsFormH
#define OptionsFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TOptionsDialog : public TForm
{
__published:	// IDE-managed Components
    TBitBtn *BitBtn1;
    TGroupBox *GroupBox1;
    TEdit *NoisePointsEdit;
    TLabel *Label1;
    TGroupBox *GroupBox2;
    TEdit *gChangeEdit;
    TLabel *Label2;
    TLabel *Label3;
    TEdit *ZFSChangeEdit;
    TLabel *Label4;
    TLabel *Label5;
    TEdit *AChangeEdit;
    TLabel *Label6;
    TEdit *WChangeEdit;
    TLabel *Label7;
    TLabel *Label8;
	TComboBox *NoiseColumnComboBox;
	TLabel *Label9;
	TCheckBox *AutoRangeCheckBox;
	TRadioButton *FromFirstRadioButton;
	TRadioButton *FromLastRadioButton;
private:	// User declarations
public:		// User declarations
    __fastcall TOptionsDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TOptionsDialog *OptionsDialog;
//---------------------------------------------------------------------------
#endif
