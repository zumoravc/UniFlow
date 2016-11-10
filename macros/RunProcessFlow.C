/* RunProcessFlow
 *
 * Steer macro for procesing flow results of AliAnalysisTaskFlowPID
 * See ProcessFlow class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

void RunProcessFlow()
{
	gSystem->AddIncludePath("-I/Users/vpacik/NBI/Flow/classProcess/");
	gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/Flow/classProcess/",gROOT->GetMacroPath()));
	gROOT->LoadMacro("ProcessFlow.cpp+g"); // loading Fitter class
	//gROOT->LoadMacro("FitPID.C+g"); // 

	ProcessFlow* process = new ProcessFlow();
	process->SetInputFilePath("~/NBI/Flow/results/test-PID-GFK-2/merge");
	process->SetTag("FB768");
	process->Run();

	return;
}