//---------------------------------------------------------------------------

#ifndef HighSpinDialogH
#define HighSpinDialogH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
    TGroupBox *GroupBox1;
    TEdit *Edit1;
    TEdit *Edit2;
    TLabel *Label1;
    TLabel *Label2;
    TLabel *Label3;
    TLabel *Label4;
    TEdit *Edit3;
    TGroupBox *GroupBox2;
    TRadioButton *RadioButton1;
    TRadioButton *RadioButton2;
    TRadioButton *RadioButton3;
    TGroupBox *GroupBox3;
    TEdit *Edit4;
    TEdit *Edit5;
    TEdit *Edit6;
    TLabel *Label5;
    TLabel *Label6;
    TLabel *Label7;
private:	// User declarations
public:		// User declarations
    __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
