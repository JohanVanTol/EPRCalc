//---------------------------------------------------------------------------

#ifndef EPRSimDialogH
#define EPRSimDialogH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include <ExtCtrls.hpp>
#include <ComCtrls.hpp>
#include <iostream.h>
#include <fstream.h>
//---------------------------------------------------------------------------
class TEPRSimulationDialog : public TForm
{
__published:	// IDE-managed Components
        TGroupBox *SimParGroupBox;
        TGroupBox *EPRLineGroupBox;
        TLabel *ModulationLabel;
    TEdit *ModulationEdit;
        TLabel *Label1;
        TGroupBox *GroupBox1;
        TLabel *Label2;
        TLabel *Label3;
    TBitBtn *BitBtn1;
    TBitBtn *BitBtn2;
    TBitBtn *BitBtn4;
    TBitBtn *BitBtn5;
    TLabel *Label5;
    TLabel *Label6;
    TLabel *Label7;
    TEdit *NptsEdit;
    TEdit *StartFieldEdit;
    TLabel *Label8;
    TLabel *Label9;
    TEdit *EndFieldEdit;
    TLabel *Label11;
    TLabel *Label13;
    TLabel *Label14;
    TEdit *LwxEdit;
    TEdit *LwyEdit;
    TEdit *LwzEdit;
    TGroupBox *PopGroupBox;
    TEdit *FreqEdit;
    TEdit *TempEdit;
    TCheckBox *BoltzmannCheckBox;
    TGroupBox *GroupBox3;
    TEdit *PxEdit;
    TEdit *PzEdit;
    TEdit *PyEdit;
    TLabel *Label15;
    TLabel *Label16;
    TLabel *Label17;
    TLabel *Label18;
    TEdit *AngleStepEdit;
    TLabel *Label19;
    TEdit *PopRatioEdit;
    TLabel *Label20;
    TLabel *Phaselabel;
    TEdit *DispersEdit;
    TUpDown *UpDown1;
    TComboBox *CalcOrderComboBox;
    TLabel *Label22;
        TLabel *Label23;
    TCheckBox *SOISCcheckbox;
    TCheckBox *RPISCcheckbox;
    TEdit *YZratioEdit;
    TLabel *Label24;
    TEdit *XZratioEdit;
    TRadioButton *DerivRadioButton;
    TRadioButton *IntegrRadioButton;
	TLabel *wxLabel;
	TLabel *wyLabel;
	TLabel *wzLabel;
    TLabel *Label27;
    TLabel *Label28;
    TEdit *MixEdit;
    TUpDown *MixUpDown;
	TLabel *VoigtAlphaLabel;
    TRadioButton *PowderRadioButton;
    TRadioButton *CrystalRadioButton;
    TEdit *ThetaEdit;
    TEdit *PhiEdit;
    TLabel *Label10;
    TLabel *Label12;
    TLabel *Label29;
    TLabel *Label30;
    TBevel *Bevel1;
    TEdit *StrainEdit;
    TGroupBox *PartialOrientationGroupBox;
    TRadioButton *Non_OrientedRadioButton;
    TRadioButton *OrientedRadioButton;
    TLabel *Label33;
    TLabel *Label34;
    TEdit *PartOrThetaEdit;
    TEdit *PartOrPhiEdit;
    TEdit *OrderParEdit;
    TGroupBox *GroupBox2;
    TEdit *kxEdit;
    TEdit *kyEdit;
    TEdit *kzEdit;
    TEdit *SLRxEdit;
    TEdit *SLRyEdit;
    TEdit *SLRzEdit;
    TLabel *Label35;
    TLabel *Label36;
    TLabel *Label37;
    TLabel *Label38;
    TLabel *Label39;
    TLabel *Label40;
    TEdit *TimeEdit;
    TCheckBox *FixedbCheckBox;
    TGroupBox *GroupBox4;
	TGroupBox *PolarizationGroupBox;
	TComboBox *PolarizationListBox;
	TRadioGroup *ShapeButtonGroup;
    void __fastcall PxEditChange(TObject *Sender);
    void __fastcall BitBtn4Click(TObject *Sender);
    void __fastcall PopRatioEditChange(TObject *Sender);
    void __fastcall AngleStepEditChange(TObject *Sender);
    void __fastcall TempEditChange(TObject *Sender);
    void __fastcall FreqEditChange(TObject *Sender);
    void __fastcall LwxEditChange(TObject *Sender);
    void __fastcall LwyEditChange(TObject *Sender);
    void __fastcall LwzEditChange(TObject *Sender);
    void __fastcall GwxEditChange(TObject *Sender);
    void __fastcall GwyEditChange(TObject *Sender);
    void __fastcall GwzEditChange(TObject *Sender);
    void __fastcall DispersEditChange(TObject *Sender);
    void __fastcall FormCreate(TObject *Sender);
    void __fastcall CalcOrderComboBoxChange(TObject *Sender);
    void __fastcall BitBtn5Click(TObject *Sender);
    void __fastcall XZratioEdit1Change(TObject *Sender);
        void __fastcall XZratioEditChange(TObject *Sender);
        void __fastcall YZratioEditChange(TObject *Sender);
	void __fastcall ShapeButtonGroupClick(TObject *Sender);
private:	// User declarations
    bool NumberChange;
    bool WidthChange;
    bool PopChange;
    bool DisPersChange;
    bool HamilChange;
public:		// User declarations
        __fastcall TEPRSimulationDialog(TComponent* Owner);
    bool GetPopChange() const {return PopChange;}
    bool GetWidthChange() const {return WidthChange;}
    bool GetNumberChange() const {return NumberChange;}
    bool GetDisPersChange() const {return DisPersChange;}
    bool GetHamilChange() const {return HamilChange;}
    void SetPopChange(bool val) { PopChange = val;}
    void SetWidthChange(bool val) { WidthChange = val;}
    void SetNumberChange(bool val) { NumberChange = val;}
    void SetDisPersChange(bool val) { DisPersChange = val;}
    void SetHamilChange(bool val) { HamilChange = val;}
	int Write(ofstream *outfile);
	int Write(AnsiString *SimParString);
};
//---------------------------------------------------------------------------
extern PACKAGE TEPRSimulationDialog *EPRSimulationDialog;
//---------------------------------------------------------------------------
#endif
