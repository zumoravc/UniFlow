/* RunProcessUniFlow
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunProcessUniFlow(const char* sOutputFilePath = "")
{
	gSystem->AddIncludePath("-I/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/",gROOT->GetMacroPath()));
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	ProcessUniFlow* process = new ProcessUniFlow();

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_CENT_wSDD_16q/");
	// process->SetInputFileName("AnalysisResults_CENTwSDD_16q.root");
	// process->SetTaskName("UniFlow_V0A");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow_V0A");
	// process->SetOutputFilePath(sOutputFilePath);
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/testRefs");
	process->SetOutputFileName("UniFlow_test.root");
	// process->SetOutputFileMode("UPDATE");
	process->SetDebug();
	// process->SuggestMultBinning(4);

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,30,45,75,180};
	// Double_t dMultBinning[] = {0,10};
	Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80};
	Int_t iSizeMult = sizeof(dMultBinning)/sizeof(dMultBinning[0]);
	process->SetMultiplicityBins(dMultBinning,iSizeMult);

	FlowTask* task1 = new FlowTask("Refs",FlowTask::kRefs);
	task1->SetHarmonics(2);
	task1->SetEtaGap(0.8);
	process->AddTask(task1);

	Double_t dPt[] = {0.1,0.2,0.3,0.4,0.5,1.,2.};

	FlowTask* task2 = new FlowTask("PionsTest",FlowTask::kPion);
	task2->SetHarmonics(2);
	// task2->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	task2->SetEtaGap(0.8);
	// process->AddTask(task2);
	//
	// FlowTask* task3 = new FlowTask("Pions2",FlowTask::kPion);
	// task3->SetHarmonics(2);
	// // task3->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	// task3->SetEtaGap(0.8);
	// process->AddTask(task3);

	// Double_t dPtBinning[] = {0.5,1.};
	// Double_t dPtBinning[] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.};
	// Int_t iSizePt = sizeof(dPtBinning)/sizeof(dPtBinning[0]);

	// FlowTask* task7 = new FlowTask("K0s",FlowTask::kK0s);
	// task7->SetHarmonics(2);
	// task7->SetEtaGap(0.8);
	// task7->SetInvMassRebin(2);
	// task7->SetFlowMassRebin(2);
	// // task7->SetShowMultDist(kTRUE);
	// task7->SetPtBins(dPtBinning,iSizePt);
	// // task7->SuggestPtBinning(1,30000);
	// process->AddTask(task7);
	//

	// FlowTask* task8 = new FlowTask("K0s",FlowTask::kK0s);
	// task8->SetHarmonics(2);
	// task8->SetEtaGap(0.);
	// task8->SetPtBins(dPtBinning,iSizePt);
	// process->AddTask(task8);
	// //
	//
	// FlowTask* task9 = new FlowTask("K0s",FlowTask::kK0s);
	// task9->SetHarmonics(2);
	// task9->SetEtaGap(0.8);
	// // task9->SetPtBins(dPtBinning,iSizePt);
	// task9->SuggestPtBinning(1,150000);
	// process->AddTask(task9);
	// //
	// FlowTask* task1 = new FlowTask("Lambda",FlowTask::kLambda);
	// task1->SetHarmonics(2);
	// task1->SetEtaGap(0.);
	// task1->SetPtBins(dPtBinning,iSizePt);
	// process->AddTask(task1);
	//
	// FlowTask* task5 = new FlowTask("Lambda",FlowTask::kLambda);
	// task5->SetHarmonics(2);
	// task5->SetEtaGap(-1.);
	// task5->SetPtBins(dPtBinning,iSizePt);
	// process->AddTask(task5);
	//
	// FlowTask* task6 = new FlowTask("Lambda",FlowTask::kLambda);
	// task6->SetHarmonics(2);
	// task6->SetEtaGap(0.8);
	// // task6->SetPtBins(dPtBinning,iSizePt);
	// task6->SuggestPtBinning(1,150000);
	// process->AddTask(task6);
	//

	// Double_t dPtBinningPhi[] = {0.5,0.8,1.1,1.4,1.7,2.,2.3,2.6,2.9,3.2,3.5,4,4.5,5.};
	// Double_t dPtBinningPhi[] = {1.,1.5,2.05};
	// Int_t iSizePtPhi = sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]);


	// FlowTask* task2 = new FlowTask("Phi",FlowTask::kPhi);
	// task2->SetHarmonics(2);
	// task2->SetEtaGap(-1.);
	// task2->SetPtBins(dPtBinningPhi,iSizePtPhi);
	// process->AddTask(task2);

	// FlowTask* task3 = new FlowTask("Phi",FlowTask::kPhi);
	// task3->SetHarmonics(2);
	// task3->SetEtaGap(0.8);
	// task3->SetPtBins(dPtBinningPhi,iSizePtPhi);
	// task3->SuggestPtBinning(1,20000);
	// task3->SetInvMassRebin(2);
	// task3->SetFlowMassRebin(2);
	// process->AddTask(task3);
	//
	// FlowTask* task4 = new FlowTask("Phi",FlowTask::kPhi);
	// task4->SetHarmonics(2);
	// task4->SetEtaGap(0.8);
	// task4->SetPtBins(dPtBinningPhi,iSizePtPhi);
	// task4->SetShowMultDist(kTRUE);
 // 	process->AddTask(task4);

	process->Run();


	return;
}
