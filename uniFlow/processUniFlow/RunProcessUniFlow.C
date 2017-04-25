/* RunProcessUniFlow
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunProcessUniFlow()
{
	gSystem->AddIncludePath("-I/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/",gROOT->GetMacroPath()));
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	ProcessUniFlow* process = new ProcessUniFlow();

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetInputFileName("AnalysisResults.root");
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test/test");
	process->SetTaskName("UniFlow");
	process->SetDebug();

	// setting multiplicity binning
	Double_t dMultBinning[] = {0,10,30,50,100,150};
	Int_t iSizeMult = sizeof(dMultBinning)/sizeof(dMultBinning[0]);
	process->SetMultiplicityBins(dMultBinning,iSizeMult);

	FlowTask* task1 = new FlowTask("testing",FlowTask::kRefs);
	task1->SetHarmonics(2);
	task1->SetEtaGap(-1.);
	process->AddTask(task1);

	Double_t dPtBinning[] = {0.5,1.,2.,3.};
	Int_t iSizePt = sizeof(dPtBinning)/sizeof(dPtBinning[0]);

	FlowTask* task2 = new FlowTask("K0s",FlowTask::kK0s);
	task2->SetHarmonics(2);
	task2->SetEtaGap(-1.);

	task2->SetPtBins(dPtBinning,iSizePt);
	process->AddTask(task2);

	process->Run();


	return;
}
