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
	process->SetInputFilePath("~/NBI/Flow/results/V0s/PID-GFK-3/merge");
	process->SetOutputFilePath("~/NBI/Flow/classProcess");
	process->SetTag("FB768");
	process->SetNumBinsCentrality(9);
	process->SetNumSamples(5);
	process->Run();

	delete process;

	return;
}