//---------------------------------------------------------------------------

#ifndef ExponentialDialogH
#define ExponentialDialogH
//---------------------------------------------------------------------------
#include <iostream.h>
#include <fstream.h>
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TExponentialDecayForm : public TForm
{
__published:	// IDE-managed Components
	TBitBtn *OKButton;
	TBitBtn *CancelButton;
	TButton *SimulateButton;
	TEdit *StartTimeEdit;
	TEdit *EndTimeEdit;
	TGroupBox *SingleExpGroupBox;
	TCheckBox *y0CheckBox;
	TCheckBox *Amp1CheckBox;
	TCheckBox *Exp1CheckBox;
	TEdit *y0Edit;
	TEdit *Amp1Edit;
	TEdit *Exp1Edit;
	TEdit *y0ErrorEdit;
	TEdit *Amp1ErrorEdit;
	TEdit *Exp1ErrorEdit;
	TEdit *Amp2ErrorEdit;
	TEdit *NptsEdit;
	TButton *CorrectBaseLineButton;
	TButton *FitButton;
	TCheckBox *Amp2CheckBox;
	TCheckBox *Exp2CheckBox;
	TCheckBox *GaussAmpCheckBox;
	TCheckBox *GaussZeroCheckBox;
	TCheckBox *GaussSigmaCheckBox;
	TEdit *Amp2Edit;
	TEdit *Exp2Edit;
	TEdit *GaussAmpEdit;
	TEdit *GaussZeroEdit;
	TEdit *GaussSigmaEdit;
	TEdit *Exp2ErrorEdit;
	TEdit *GaussAmpErrorEdit;
	TEdit *GaussZeroErrorEdit;
	TEdit *GaussSigmaErrorEdit;
	TRadioButton *BiRadioButton;
	TRadioButton *MonoGaussRadioButton;
	TRadioButton *MonoRadioButton;
	TButton *Button1;
	void __fastcall SimulateButtonClick(TObject *Sender);
	void __fastcall CorrectBaseLineButtonClick(TObject *Sender);
	void __fastcall Fit10ButtonClick(TObject *Sender);
	void __fastcall FitButtonClick(TObject *Sender);
	void __fastcall CancelButtonClick(TObject *Sender);
	void __fastcall BiRadioButtonClick(TObject *Sender);
	void __fastcall MonoGaussRadioButtonClick(TObject *Sender);
	void __fastcall MonoRadioButtonClick(TObject *Sender);
	void __fastcall OKButtonClick(TObject *Sender);
private:	// User declarations
	int mode;
public:		// User declarations
	__fastcall TExponentialDecayForm(TComponent* Owner);
	int GetLimits(double* start, double* stop);
	int GetParameters(double* par, int* fixed);
	int SetParameters(int na, double* par, int* fixed);
	void SetErrors(int na, double* errors);
	int Write(ofstream *somefile);
	int Write(AnsiString* ExParString);
};
//---------------------------------------------------------------------------
extern PACKAGE TExponentialDecayForm *ExponentialDecayForm;
//---------------------------------------------------------------------------
#endif
