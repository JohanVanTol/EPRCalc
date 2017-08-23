//---------------------------------------------------------------------------
#ifndef PowderEndorDlgH
#define PowderEndorDlgH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include "ENDORsimParam.h"
#include <Dialogs.hpp>
//---------------------------------------------------------------------------
class TPowderENDORDialog : public TForm
{
__published:	// IDE-managed Components
    TGroupBox *GroupBox1;
    TEdit *FreqEdit;
    TLabel *Label1;
    TGroupBox *GroupBox2;
    TGroupBox *GroupBox3;
    TLabel *Label3;
    TEdit *NptsEdit;
    TLabel *Label4;
    TEdit *FirstEdit;
    TLabel *Label5;
    TLabel *Label6;
    TEdit *LastEdit;
    TLabel *Label7;
    TEdit *NgEdit;
    TLabel *Label2;
    TLabel *Label8;
    TEdit *gStartEdit;
    TEdit *gStepEdit;
    TLabel *Label9;
    TEdit *WyEdit;
    TLabel *Label10;
    TLabel *Label11;
    TLabel *Label12;
    TEdit *WzEdit;
    TEdit *WxEdit;
    TLabel *Label13;
    TLabel *Label14;
    TLabel *Label15;
    TLabel *Label16;
    TEdit *AngleStepEdit;
    TLabel *Label17;
    TGroupBox *GroupBox4;
    TEdit *ENDORwidthEdit;
    TRadioButton *RadioButton1;
    TRadioButton *RadioButton2;
    TLabel *Label18;
    TLabel *Label19;
    TLabel *Label20;
    TEdit *ModulationEdit;
    TLabel *Label21;
    TGroupBox *GroupBox5;
    TComboBox *SpinMultComboBox;
    TLabel *Label22;
    TEdit *ANxEdit;
    TEdit *ANyEdit;
    TEdit *ANzEdit;
    TLabel *Label23;
    TLabel *Label24;
    TLabel *Label25;
    TLabel *Label26;
    TLabel *Label27;
    TLabel *Label28;
    TBitBtn *OkButton;
    TBitBtn *CanvelButton;
    TBitBtn *LoadButton;
    TBitBtn *SaveButton;
    TOpenDialog *SimOpenDialog;
    TSaveDialog *SimSaveDialog;
    void __fastcall FormActivate(TObject *Sender);
    void __fastcall OkButtonClick(TObject *Sender);
    void __fastcall SpinMultComboBoxChange(TObject *Sender);
    void __fastcall CanvelButtonClick(TObject *Sender);
    void __fastcall LoadButtonClick(TObject *Sender);
    void __fastcall NptsEditExit(TObject *Sender);
    void __fastcall FirstEditExit(TObject *Sender);
    void __fastcall LastEditExit(TObject *Sender);
    void __fastcall AngleStepEditExit(TObject *Sender);
    void __fastcall ENDORwidthEditExit(TObject *Sender);
    void __fastcall ModulationEditExit(TObject *Sender);
    void __fastcall FreqEditExit(TObject *Sender);
    void __fastcall NgEditExit(TObject *Sender);
    void __fastcall gStartEditExit(TObject *Sender);
    void __fastcall gStepEditExit(TObject *Sender);
    void __fastcall WxEditExit(TObject *Sender);
    void __fastcall WyEditExit(TObject *Sender);
    void __fastcall WzEditExit(TObject *Sender);
    void __fastcall ANxEditExit(TObject *Sender);
    void __fastcall ANyEditExit(TObject *Sender);
    void __fastcall ANzEditExit(TObject *Sender);
private:	// User declarations
    PowderEndorSimPar Par;
    int UpdateParameters();
    int SetAllParameters();
public:		// User declarations
    __fastcall TPowderENDORDialog(TComponent* Owner);
    int SetPar(PowderEndorSimPar P) {Par = P;return 0;}
    PowderEndorSimPar GetPar() {return Par;}
};
//---------------------------------------------------------------------------
extern PACKAGE TPowderENDORDialog *PowderENDORDialog;
//---------------------------------------------------------------------------
#endif
