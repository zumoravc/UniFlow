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
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetTaskName("UniFlow");
	process->SetDebug();

	FlowTask* task1 = new FlowTask("testing",FlowTask::kRefs);
	task1->SetHarmonics(2);
	task1->SetEtaGap(-10);
	process->AddTask(task1);

	process->Run();


	return;
}
