/* RunFitPID
 *
 * Steer macro for fitting identied hadrons and extracting the flow.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

void RunFitPID()
{
	gSystem->AddIncludePath("-I/Users/vpacik/NBI/Flow/fitter/");
	gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/Flow/fitter/",gROOT->GetMacroPath()));
	gROOT->LoadMacro("Fitter.cpp+g"); // loading Fitter class
	gROOT->LoadMacro("FitPID.C+g"); // 

	FitPID();
	
	return;
}