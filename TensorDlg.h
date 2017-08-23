//----------------------------------------------------------------------------
#ifndef OCBH
#define OCBH
//----------------------------------------------------------------------------
#include <System.hpp>
#include <Windows.hpp>
#include <SysUtils.hpp>
#include <Classes.hpp>
#include <Graphics.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Controls.hpp>
#include <Buttons.hpp>
#include <ExtCtrls.hpp>
#include <Mask.hpp>
#include "Tensor.h" 

//----------------------------------------------------------------------------
class TTensorDialog : public TForm
{
__published:
  TGroupBox *TensorBox;
  TLabel *x;
  TLabel *Label1;
  TLabel *Label2;
  TMaskEdit *XX;
  TMaskEdit *XY;
  TMaskEdit *XZ;
  TMaskEdit *YZ;
  TMaskEdit *ZZ;
  TMaskEdit *YY;
  TGroupBox *PrincipalBox;
  TLabel *Label3;
  TLabel *Label4;
  TLabel *Label5;
  TLabel *Label6;
  TLabel *Label7;
  TLabel *Label8;
  TEdit *Value1;
  TEdit *Value2;
  TEdit *Value3;
  TEdit *Theta1;
  TEdit *Theta2;
  TEdit *Theta3;
  TEdit *Phi1;
  TEdit *Phi2;
  TEdit *Phi3;
    TBitBtn *BitBtn1;
    TBitBtn *BitBtn2;
  void __fastcall FormActivate(TObject *Sender);
  
  void __fastcall XXChange(TObject *Sender);
  
  void __fastcall XYChange(TObject *Sender);
  void __fastcall XZChange(TObject *Sender);
  void __fastcall YYChange(TObject *Sender);
  void __fastcall ZZChange(TObject *Sender);
  void __fastcall YZChange(TObject *Sender);
    void __fastcall BitBtn1Click(TObject *Sender);

  
private:
  Tensor* Tens;
  double Val[3];
  double Theta[3];
  double Phi[3];
  bool Axial;

public:
	virtual __fastcall TTensorDialog(TComponent* AOwner);
  void UpdatePrinciples();
  void SetTensor(Tensor *T);
  Tensor GetTensor() {return *Tens;}
  void SetAxial(int ax=1);
};
//----------------------------------------------------------------------------
extern PACKAGE TTensorDialog *TensorDialog;
//----------------------------------------------------------------------------
#endif
