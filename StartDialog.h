//---------------------------------------------------------------------------
#ifndef StartDialogH
#define StartDialogH
//---------------------------------------------------------------------------
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TSpectralLineStartDialog : public TForm
{
__published:	// IDE-managed Components
    TComboBox *NumberOfLinesListBox;
    TGroupBox *StartParamGroupBox;
    TEdit *A0Edit;
    TEdit *A1Edit;
    TEdit *A2Edit;
    TEdit *A5Edit;
    TEdit *A8Edit;
    TEdit *A11Edit;
    TEdit *A14Edit;
    TEdit *A17Edit;
    TEdit *A20Edit;
    TEdit *A23Edit;
    TEdit *A26Edit;
    TEdit *A29Edit;
    TEdit *A3Edit;
    TEdit *A6Edit;
    TEdit *A9Edit;
    TEdit *A12Edit;
    TEdit *A15Edit;
    TEdit *A18Edit;
    TEdit *A21Edit;
    TEdit *A24Edit;
    TEdit *A27Edit;
    TEdit *A30Edit;
    TEdit *A4Edit;
    TEdit *A7Edit;
    TEdit *A10Edit;
    TEdit *A13Edit;
    TEdit *A16Edit;
    TEdit *A19Edit;
    TEdit *A22Edit;
    TEdit *A25Edit;
    TEdit *A28Edit;
    TEdit *A31Edit;
    TLabel *A0Label;
    TCheckBox *A0CheckBox;
    TCheckBox *A1CheckBox;
    TCheckBox *A5CheckBox;
    TCheckBox *A8CheckBox;
    TCheckBox *A11CheckBox;
    TCheckBox *A14CheckBox;
    TCheckBox *A17CheckBox;
    TCheckBox *A20CheckBox;
    TCheckBox *A23CheckBox;
    TCheckBox *A26CheckBox;
    TCheckBox *A29CheckBox;
    TCheckBox *A2CheckBox;
    TCheckBox *A7CheckBox;
    TCheckBox *A10CheckBox;
    TCheckBox *A13CheckBox;
    TCheckBox *A16CheckBox;
    TCheckBox *A19CheckBox;
    TCheckBox *A22CheckBox;
    TCheckBox *A25CheckBox;
    TCheckBox *A28CheckBox;
    TCheckBox *A31CheckBox;
    TCheckBox *A4CheckBox;
    TBitBtn *BitBtn1;
    TLabel *Label8;
    TLabel *A1Label;
    TEdit *F1Edit;
    TCheckBox *F1CheckBox;
    TButton *ChiSqrButton;
    TEdit *NptsEdit;
    TLabel *Label30;
    TEdit *StartEdit;
    TLabel *Label31;
    TEdit *StopEdit;
    TEdit *F2Edit;
    TEdit *F3Edit;
    TEdit *F4Edit;
    TEdit *F5Edit;
    TEdit *F6Edit;
    TEdit *F7Edit;
    TEdit *F8Edit;
    TEdit *F9Edit;
    TEdit *F10Edit;
    TCheckBox *F2CheckBox;
    TCheckBox *F3CheckBox;
    TCheckBox *F4CheckBox;
    TCheckBox *F5CheckBox;
    TCheckBox *F6CheckBox;
    TCheckBox *F7CheckBox;
    TCheckBox *F8CheckBox;
    TCheckBox *F9CheckBox;
    TCheckBox *F10CheckBox;
    TButton *FitCycleButton;
    TCheckBox *DerivCheckBox;
    TButton *Fit10Button;
	TLabel *Line1Label;
	TLabel *Line2Label;
	TLabel *Line4Label;
	TLabel *Line5Label;
	TLabel *Line6Label;
	TLabel *Line7Label;
	TLabel *Line8Label;
	TLabel *Line9Label;
	TLabel *Line10Label;
	TLabel *Line3Label;
    TBevel *Bevel1;
    TCheckBox *A3CheckBox;
    TCheckBox *A6CheckBox;
    TCheckBox *A9CheckBox;
    TCheckBox *A12CheckBox;
    TCheckBox *A15CheckBox;
    TCheckBox *A18CheckBox;
    TCheckBox *A21CheckBox;
    TCheckBox *A24CheckBox;
    TCheckBox *A27CheckBox;
    TCheckBox *A30CheckBox;
    TLabel *Label1;
    TLabel *Label2;
    TBevel *Bevel2;
    TBevel *Bevel3;
    TLabel *Label4;
    TLabel *Label7;
    TBevel *Bevel4;
    TComboBox *Mult1ListBox;
    TLabel *Label9;
    TLabel *Label13;
    TEdit *S1Edit;
    TCheckBox *S1CheckBox;
    TComboBox *Mult2ListBox;
    TComboBox *Mult3ListBox;
    TComboBox *Mult4ListBox;
    TComboBox *Mult5ListBox;
    TComboBox *Mult6ListBox;
    TComboBox *Mult7ListBox;
    TComboBox *Mult8ListBox;
    TComboBox *Mult9ListBox;
    TComboBox *Mult10ListBox;
    TEdit *S2Edit;
    TEdit *S3Edit;
    TEdit *S4Edit;
    TEdit *S5Edit;
    TEdit *S6Edit;
    TEdit *S7Edit;
    TEdit *S8Edit;
    TEdit *S9Edit;
    TEdit *S10Edit;
    TCheckBox *S2CheckBox;
    TCheckBox *S3CheckBox;
    TCheckBox *S4CheckBox;
    TCheckBox *S5CheckBox;
    TCheckBox *S6CheckBox;
    TCheckBox *S7CheckBox;
    TCheckBox *S8CheckBox;
    TCheckBox *S9CheckBox;
    TCheckBox *S10CheckBox;
    TEdit *ModulationEdit;
    TEdit *DummyEdit;
    TCheckBox *DummyCheckBox;
	TComboBox *ShapeComboBox;
        TButton *BaseLineSubButton;
	TBevel *Bevel5;
	TLabel *AlphaLabel;
	TEdit *Alpha1Edit;
	TEdit *Alpha8Edit;
	TEdit *Alpha2Edit;
	TEdit *Alpha3Edit;
	TEdit *Alpha4Edit;
	TEdit *Alpha5Edit;
	TEdit *Alpha6Edit;
	TEdit *Alpha7Edit;
	TEdit *Alpha9Edit;
	TEdit *Alpha10Edit;
	TCheckBox *Alpha1CheckBox;
	TCheckBox *Alpha2CheckBox;
	TCheckBox *Alpha3CheckBox;
	TCheckBox *Alpha4CheckBox;
	TCheckBox *Alpha5CheckBox;
	TCheckBox *Alpha6CheckBox;
	TCheckBox *Alpha7CheckBox;
	TCheckBox *Alpha8CheckBox;
	TCheckBox *Alpha9CheckBox;
	TCheckBox *Alpha10CheckBox;
    void __fastcall FormCreate(TObject *Sender);
    void __fastcall NumberOfLinesListBoxChange(TObject *Sender);
    void __fastcall FormDestroy(TObject *Sender);
    void __fastcall BitBtn1Click(TObject *Sender);
    void __fastcall ChiSqrButtonClick(TObject *Sender);
    void __fastcall FitCycleButtonClick(TObject *Sender);
    void __fastcall Fit10ButtonClick(TObject *Sender);
        void __fastcall BaseLineSubButtonClick(TObject *Sender);
	void __fastcall FormActivate(TObject *Sender);
	void __fastcall ShapeComboBoxChange(TObject *Sender);
private:	// User declarations
    int npar;
    double *par;
    int *fixed;
    TEdit **EditArray;
    TCheckBox **CheckArray;
    TComboBox **MultArray;

public:		// User declarations
    void SetParameters(int n, double* param, int* fixedpar);
    int GetParameters(double* par, int* fixed);
    int GetLimits(double* start, double* stop);
    int UpdateForm();
    __fastcall TSpectralLineStartDialog(TComponent* Owner);
	int Write(ofstream *somefile);
	int Write(AnsiString* LinePar);

};
//---------------------------------------------------------------------------
extern PACKAGE TSpectralLineStartDialog *SpectralLineStartDialog;
//---------------------------------------------------------------------------
#endif
