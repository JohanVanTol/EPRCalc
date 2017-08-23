//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
USEFORM("EPRCalcMain2.cpp", MainForm);
USEFORM("AboutDialog.cpp", AboutBox);
USEFORM("TensorDlg.cpp", TensorDialog);
USEFORM("Hamilton.cpp", HamiltonianDialog);
USEFORM("PowderEndorDlg.cpp", PowderENDORDialog);
USEFORM("EPRSimDialog.cpp", EPRSimulationDialog);
USEFORM("..\General\Dialogs\ErrorBx.cpp", ErrorBox);
USEFORM("OptionsForm.cpp", OptionsDialog);
USEFORM("FitCycleDlg.cpp", FitCycleDialog);
USEFORM("..\LSQFitProgram\DataDefineFrm.cpp", DefineDataForm);
USEFORM("StartDialog.cpp", SpectralLineStartDialog);
USEFORM("StrainForm.cpp", StrainDialog);
USEFORM("EnergyDlg.cpp", EnergyForm);
USEFORM("gStrainForm.cpp", gStrainDialog);
USEFORM("ExponentialDialog.cpp", ExponentialDecayForm);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
  try
  {
     Application->Initialize();
     Application->HelpFile = "C:\\Documents and Settings\\vantol\\My Documents\\Programs\\Calculat\\EPRCalcHelp.chm";
		Application->CreateForm(__classid(TMainForm), &MainForm);
		Application->CreateForm(__classid(TAboutBox), &AboutBox);
		Application->CreateForm(__classid(TTensorDialog), &TensorDialog);
		Application->CreateForm(__classid(THamiltonianDialog), &HamiltonianDialog);
		Application->CreateForm(__classid(TPowderENDORDialog), &PowderENDORDialog);
		Application->CreateForm(__classid(TEPRSimulationDialog), &EPRSimulationDialog);
		Application->CreateForm(__classid(TErrorBox), &ErrorBox);
		Application->CreateForm(__classid(TOptionsDialog), &OptionsDialog);
		Application->CreateForm(__classid(TFitCycleDialog), &FitCycleDialog);
		Application->CreateForm(__classid(TDefineDataForm), &DefineDataForm);
		Application->CreateForm(__classid(TSpectralLineStartDialog), &SpectralLineStartDialog);
		Application->CreateForm(__classid(TStrainDialog), &StrainDialog);
		Application->CreateForm(__classid(TEnergyForm), &EnergyForm);
		Application->CreateForm(__classid(TgStrainDialog), &gStrainDialog);
		Application->CreateForm(__classid(TExponentialDecayForm), &ExponentialDecayForm);
		Application->Run();
  }
  catch (Exception &exception)
  {
     Application->ShowException(&exception);
  }
  return 0;
}
//---------------------------------------------------------------------------


