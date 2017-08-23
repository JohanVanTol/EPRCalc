//---------------------------------------------------------------------------
#ifndef ErrorH
#define ErrorH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
#include <ExtCtrls.hpp>
#include <Graphics.hpp>
//---------------------------------------------------------------------------
class TErrorBox : public TForm
{
__published:	// IDE-managed Components
  TBitBtn *BitBtn1;
  TImage *Image;
    TPanel *Panel1;
    TLabel *Label;
private:	// User declarations
public:		// User declarations
  __fastcall TErrorBox(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TErrorBox *ErrorBox;
//---------------------------------------------------------------------------
#endif
