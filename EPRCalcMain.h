//---------------------------------------------------------------------------
#ifndef EPRCalcMainH
#define EPRCalcMainH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Dialogs.hpp>
#include <ExtCtrls.hpp>
#include <Menus.hpp>
#include <Buttons.hpp>
#include <ComCtrls.hpp>
//---------------------------------------------------------------------------
class TMainForm : public TForm
{
__published:	// IDE-managed Components
  TMainMenu *MainMenu1;
  TOpenDialog *OpenDialog;
  TSaveDialog *SaveDialog;
  TPanel *Panel1;
  TMenuItem *FileMenu;
  TMenuItem *FileExit;
  TMenuItem *N2;
  TMenuItem *SaveAs1;
  TMenuItem *Save1;
  TMenuItem *Open1;
  TMenuItem *New1;
  TMenuItem *HelpMenu;
  TMenuItem *About;
  TBevel *Bevel1;
  TMemo *OutputMemo;
  TSpeedButton *FileNewButton;
  TStatusBar *StatusBar;
  TMenuItem *Param;
  TMenuItem *ParamHamilton;
    TMenuItem *CalcMenuItem;
    TMenuItem *PowderEndorMenuItem;
        TMenuItem *PowderTransMenuItem;
        TMenuItem *PowderEPRMenuItem;
        TPanel *Panel2;
        TLabel *ThetaLabel;
        TLabel *Label4;
        TLabel *Label5;
        TLabel *PhiLabel;
    TMenuItem *Data;
    TMenuItem *SetColumnsMenuItem;
    TLabel *FileLabel;
    TComboBox *ColumnComboBox;
    TLabel *Label1;
    TSaveDialog *FitSaveDialog;
    TSpeedButton *DisPersFitButton;
    TSpeedButton *PopFitButton;
    TSpeedButton *WidthXFitButton;
    TSpeedButton *WidthYFitButton;
    TSpeedButton *WidthZFitButton;
    TSpeedButton *DFitButton;
    TSpeedButton *EFitButton;
    TSpeedButton *GxxFitButton;
    TSpeedButton *GyyFitButton;
    TSpeedButton *GzzFitButton;
    TSpeedButton *GxyFitButton;
    TSpeedButton *GxzFitButton;
    TSpeedButton *GyzFitButton;
    TMenuItem *OptionsMenuItem;
    TSpeedButton *GammaButton;
    TMenuItem *N1;
    TMenuItem *GaussianFit1;
    TMenuItem *LorentzFit1;
        TSpeedButton *XZratioFitButton;
        TSpeedButton *YZratioFitButton;
    TSpeedButton *FitCyclesButton;
    TSpeedButton *WisoFitButton;
    TSpeedButton *GisoFitButton;
    TMenuItem *HighSpin1;
    TSpeedButton *SpeedButton1;
    TSpeedButton *AisoFitButton;
    TSpeedButton *AxxFitButton;
    TSpeedButton *AyyFitButton;
    TSpeedButton *AzzFitButton;
    TSpeedButton *AxyFitButton;
    TSpeedButton *AxzFitButton;
    TSpeedButton *AyzFitButton;
    TCheckBox *StopCheckBox;
    TLabel *aLabel;
    TLabel *bLabel;
    TSpeedButton *KxFitButton;
    TSpeedButton *KyFitButton;
    TSpeedButton *KzFitButton;
    TSpeedButton *SLxFitButton;
    TSpeedButton *SLyFitButton;
    TSpeedButton *SLzFitButton;
    TSpeedButton *PzFitButton;
  void __fastcall FileExitClick(TObject *Sender);
  void __fastcall AboutClick(TObject *Sender);
  void __fastcall ParamHamiltonClick(TObject *Sender);
    void __fastcall PowderEndorMenuItemClick(TObject *Sender);
        void __fastcall PowderEPRMenuItemClick(TObject *Sender);
        void __fastcall Open1Click(TObject *Sender);
        void __fastcall FormCreate(TObject *Sender);
    void __fastcall FormDestroy(TObject *Sender);
    void __fastcall FormResize(TObject *Sender);
    void __fastcall FormPaint(TObject *Sender);
    void __fastcall PowderTransMenuItemClick(TObject *Sender);
    void __fastcall ColumnComboBoxChange(TObject *Sender);
    void __fastcall Save1Click(TObject *Sender);
    void __fastcall DisPersFitButtonClick(TObject *Sender);
    void __fastcall PopFitButtonClick(TObject *Sender);
    void __fastcall WidthXFitButtonClick(TObject *Sender);
    void __fastcall WidthYFitButtonClick(TObject *Sender);
    void __fastcall WidthZFitButtonClick(TObject *Sender);
    void __fastcall DFitButtonClick(TObject *Sender);
    void __fastcall EFitButtonClick(TObject *Sender);
    void __fastcall GxxFitButtonClick(TObject *Sender);
    void __fastcall GyyFitButtonClick(TObject *Sender);
    void __fastcall GzzFitButtonClick(TObject *Sender);
    void __fastcall GxyFitButtonClick(TObject *Sender);
    void __fastcall GxzFitButtonClick(TObject *Sender);
    void __fastcall GyzFitButtonClick(TObject *Sender);
    void __fastcall OptionsMenuItemClick(TObject *Sender);
    void __fastcall GammaButtonClick(TObject *Sender);
    void __fastcall GaussianFit1Click(TObject *Sender);
    void __fastcall XZratioFitButtonClick(TObject *Sender);
        void __fastcall YZratioFitButtonClick(TObject *Sender);
    void __fastcall FitCyclesButtonClick(TObject *Sender);
    void __fastcall WisoFitButtonClick(TObject *Sender);
    void __fastcall GisoFitButtonClick(TObject *Sender);
    void __fastcall LorentzFit1Click(TObject *Sender);
    void __fastcall HighSpin1Click(TObject *Sender);
    void __fastcall FormMouseDown(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
    void __fastcall FormMouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
    void __fastcall FormMouseUp(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
    void __fastcall SpeedButton1Click(TObject *Sender);
    void __fastcall AisoFitButtonClick(TObject *Sender);
    void __fastcall AxxFitButtonClick(TObject *Sender);
    void __fastcall AyyFitButtonClick(TObject *Sender);
    void __fastcall AzzFitButtonClick(TObject *Sender);
    void __fastcall AxyFitButtonClick(TObject *Sender);
    void __fastcall AxzFitButtonClick(TObject *Sender);
    void __fastcall AyzFitButtonClick(TObject *Sender);
    void __fastcall KxFitButtonClick(TObject *Sender);
    void __fastcall KyFitButtonClick(TObject *Sender);
    void __fastcall KzFitButtonClick(TObject *Sender);
    void __fastcall SLxFitButtonClick(TObject *Sender);
    void __fastcall SLyFitButtonClick(TObject *Sender);
    void __fastcall SLzFitButtonClick(TObject *Sender);
    void __fastcall PzFitButtonClick(TObject *Sender);


private:	// User declarations
    TRect *FocusRect;
    
public:		// User declarations
  __fastcall TMainForm(TComponent* Owner);
  int PowderEndor();
  int ReadExpData(AnsiString FileName);
  int PowderEPR();
  void PlotData();
  double CalcCycle();
  double HighSpinCalc();
  void SimulateLorentz();
  void SimulateGauss();
  void FitLorentz(int maxCycle);
  void GijFitCycle(int k, int l);
  void AijFitCycle(int k, int l);
  void GAngleFitCycle(int k);
  int ResetArrays(int npts, int NewNTrans=2);
  
};
//---------------------------------------------------------------------------
extern PACKAGE TMainForm *MainForm;
//---------------------------------------------------------------------------
#endif
