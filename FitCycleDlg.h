//---------------------------------------------------------------------------

#ifndef FitCycleDlgH
#define FitCycleDlgH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include <ComCtrls.hpp>
//---------------------------------------------------------------------------
class TFitCycleDialog : public TForm
{
__published:	// IDE-managed Components
    TGroupBox *HamiltonianFitGroupBox;
    TGroupBox *GroupBox2;
    TGroupBox *GroupBox3;
    TCheckBox *DFitCheckBox;
    TCheckBox *GxxFitCheckBox;
    TCheckBox *GyyFitCheckBox;
    TCheckBox *GzzFitCheckBox;
    TCheckBox *GxyFitCheckBox;
    TCheckBox *GxzFitCheckBox;
    TCheckBox *GyzFitCheckBox;
    TCheckBox *AbsDisFitCheckBox;
    TCheckBox *ISCPopRatioFitCheckBox;
    TCheckBox *XZratioFitCheckBox;
    TCheckBox *YZratioFitCheckBox;
    TCheckBox *WxFitCheckBox;
    TCheckBox *WyFitCheckBox;
    TCheckBox *WzFitCheckBox;
    TCheckBox *EFitCheckBox;
    TGroupBox *GroupBox1;
    TEdit *NCycleEdit;
    TLabel *Ncycles;
    TEdit *ImprovementEdit;
    TLabel *Label56;
    TLabel *Label1;
    TBitBtn *BitBtn1;
    TBitBtn *BitBtn2;
    TUpDown *UpDown1;
    TGroupBox *GroupBox4;
    TCheckBox *AxxFitCheckBox;
    TCheckBox *AyyFitCheckBox;
    TCheckBox *AzzFitCheckBox;
    TCheckBox *AxyFitCheckBox;
    TCheckBox *AxzFitCheckBox;
    TCheckBox *AyzFitCheckBox;
    TCheckBox *AisoFitCheckBox;
    TCheckBox *GisoCheckBox;
    TCheckBox *WisoFitCheckBox;
    TGroupBox *GroupBox5;
    TCheckBox *PzFitCheckBox;
    TCheckBox *KxFitCheckBox;
    TCheckBox *KyFitCheckBox;
    TCheckBox *KzFitCheckBox;
    TCheckBox *SLRxFitCheckBox;
    TCheckBox *SLRyFitCheckBox;
    TCheckBox *SLRzFitCheckBox;
    TCheckBox *DStrainFitCheckBox;
    TCheckBox *EStrainFitCheckBox;
	TCheckBox *B40FitCheckBox;
    void __fastcall BitBtn1Click(TObject *Sender);
    void __fastcall BitBtn2Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
    __fastcall TFitCycleDialog(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TFitCycleDialog *FitCycleDialog;
//---------------------------------------------------------------------------
#endif
