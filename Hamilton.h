//---------------------------------------------------------------------------
#ifndef HamiltonH
#define HamiltonH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include <Mask.hpp>
#include <Dialogs.hpp>
#include "HamilParam.h"
//---------------------------------------------------------------------------
class THamiltonianDialog : public TForm
{
__published:	// IDE-managed Components
  TComboBox *SpinComboBox;
  TGroupBox *ElectronSpinGroup;
  TLabel *Label1;
  TBitBtn *OKButton;
  TBitBtn *CancelButton;
  TGroupBox *NucSpinGroup;
  TButton *GTensorButton;
  TGroupBox *CrystalFieldBox;
  TLabel *Label3;
  TLabel *Label4;
  TLabel *Label5;
  TLabel *Label6;
  TLabel *Label7;
  TLabel *Label8;
  TLabel *Label9;
  TLabel *Label10;
  TLabel *Label11;
  TLabel *Label12;
  TEdit *B20Edit;
  TEdit *B22Edit;
  TEdit *B40Edit;
  TEdit *B42Edit;
  TEdit *B43Edit;
  TEdit *B44Edit;
  TEdit *B60Edit;
  TEdit *B63Edit;
  TEdit *B64Edit;
  TEdit *B66Edit;
  TLabel *Label2;
  TLabel *Label13;
  TComboBox *SpinNumberBox;
  TLabel *Label14;
  TComboBox *NucSpinComboBox;
  TButton *HyperfineButton;
  TButton *QuadrupoleButton;
  TLabel *Label15;
  TEdit *NucZeemanEdit;
  TEdit *NNucSpinEdit;
  TBitBtn *LoadHamilButton;
  TBitBtn *SaveHamilButton;
  TBitBtn *SaveAsHamilButton;
  TOpenDialog *HamiltonOpenDialog;
  TSaveDialog *HamiltonSaveDialog;
    TButton *AddButton;
    TButton *DeleteButton;
    TCheckBox *AxialCheckBox;
    TRadioButton *DiagonButton;
    TRadioButton *PerturButton;
    TButton *StrainButton;
	TButton *gStrainButton;
  void __fastcall GTensorButtonClick(TObject *Sender);
  
  
  
  void __fastcall LoadHamilButtonClick(TObject *Sender);
    void __fastcall SaveHamilButtonClick(TObject *Sender);
    void __fastcall FormActivate(TObject *Sender);
    void __fastcall SpinComboBoxChange(TObject *Sender);
    void __fastcall NNucSpinEditExit(TObject *Sender);
    void __fastcall AddButtonClick(TObject *Sender);
    void __fastcall HyperfineButtonClick(TObject *Sender);
    void __fastcall QuadrupoleButtonClick(TObject *Sender);
    void __fastcall SpinNumberBoxChange(TObject *Sender);
    void __fastcall NucSpinComboBoxChange(TObject *Sender);
    void __fastcall OKButtonClick(TObject *Sender);
    void __fastcall DeleteButtonClick(TObject *Sender);
    void __fastcall NucZeemanEditExit(TObject *Sender);
    void __fastcall B20EditExit(TObject *Sender);
    void __fastcall B22EditExit(TObject *Sender);
    void __fastcall B40EditExit(TObject *Sender);
    void __fastcall B42EditExit(TObject *Sender);
    void __fastcall B43EditExit(TObject *Sender);
    void __fastcall B44EditExit(TObject *Sender);
    void __fastcall B60EditExit(TObject *Sender);
    void __fastcall B64EditExit(TObject *Sender);
    void __fastcall B66EditExit(TObject *Sender);
    void __fastcall B63EditExit(TObject *Sender);
    void __fastcall StrainButtonClick(TObject *Sender);
	void __fastcall gStrainButtonClick(TObject *Sender);

private:	// User declarations
    HamilPar Par;
public:		// User declarations
    void SetParameters(const HamilPar& Hparam) {Par = Hparam;}
    HamilPar GetParameters() { return Par;}
  __fastcall THamiltonianDialog(TComponent* Owner);
    void UpdateParameters();
};
//---------------------------------------------------------------------------
extern PACKAGE THamiltonianDialog *HamiltonianDialog;
//---------------------------------------------------------------------------
#endif
