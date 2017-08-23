//---------------------------------------------------------------------------

#ifndef EnergyFormH
#define EnergyFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
//---------------------------------------------------------------------------
class TEnergyForm : public TForm
{
__published:	// IDE-managed Components
	TRadioButton *FixedFieldRadioButton;
	TRadioButton *FixedOrientRadioButton;
	TEdit *Edit1;
	TEdit *Edit2;
	TLabel *Label1;
	TLabel *Label2;
private:	// User declarations
public:		// User declarations
	__fastcall TEnergyForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TEnergyForm *EnergyForm;
//---------------------------------------------------------------------------
#endif
