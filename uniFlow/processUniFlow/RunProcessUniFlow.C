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

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A");
	// process->SetInputFileName("AnalysisResults_CENTwSDD_16q.root");
	// process->SetTaskName("UniFlow_V0A");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetInputFileName("AnalysisResults_merged.root");
	process->SetTaskName("UniFlow_V0A");
	process->SetOutputFilePath(sOutputFilePath);
	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/testRefs");
	process->SetOutputFileName("UniFlow_DesampleTest.root");
	// process->SetOutputFileMode("UPDATE");
	// process->SetDebug();
	// process->SuggestMultBinning(4);

	// setting multiplicity binning
	Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	// Double_t dMultBinning[] = {0,10};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	FlowTask* task1 = new FlowTask("Refs_Rebin",FlowTask::kRefs);
	task1->SetHarmonics(2);
	task1->SetEtaGap(0.8);
	process->AddTask(task1);

	FlowTask* task2 = new FlowTask("Refs_noRebin",FlowTask::kRefs);
	task2->SetHarmonics(2);
	task2->SetEtaGap(0.8);
	task2->SetRebinning(kFALSE);
	process->AddTask(task2);

	// DIRECT

	// Double_t dPt[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.8,3.,3.2,3.4,3.6,3.8,4.,4.2,4.4,4.6,4.8,5.,5.5,6.,6.5,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};

	// FlowTask* task2 = new FlowTask("Pi",FlowTask::kPion);
	// task2->SetHarmonics(2);
	// task2->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	// task2->SetEtaGap(-1);
	// process->AddTask(task2);
	//
	// FlowTask* task3 = new FlowTask("K",FlowTask::kKaon);
	// task3->SetHarmonics(2);
	// task3->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	// task3->SetEtaGap(-1);
	// process->AddTask(task3);
	//
	// FlowTask* task4 = new FlowTask("p",FlowTask::kProton);
	// task4->SetHarmonics(2);
	// task4->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	// task4->SetEtaGap(-1);
	// process->AddTask(task4);
	//
	// FlowTask* task5 = new FlowTask("charged",FlowTask::kCharged);
	// task5->SetHarmonics(2);
	// // task5->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	// task5->SetEtaGap(0.8);
	// process->AddTask(task5);

	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningLambda[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,10.,20.};


	// Double_t dPtBinning[] = {0.5,1.};
	// Double_t dPtBinning[] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.};
	// Int_t iSizePt = sizeof(dPtBinning)/sizeof(dPtBinning[0]);

 // FlowTask* task7 = new FlowTask("K0s",FlowTask::kK0s);
 // task7->SetHarmonics(2);
 // task7->SetEtaGap(0.8);
 // task7->SetInvMassRebin(2);
 // task7->SetFlowMassRebin(2);
 // // task7->SetShowMultDist(kTRUE);
 // task7->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
 // // task7->SuggestPtBinning(1,30000);
 // process->AddTask(task7);
 //
 //
 // // FlowTask* task8 = new FlowTask("K0s",FlowTask::kK0s);
 // // task8->SetHarmonics(2);
 // // task8->SetEtaGap(0.);
 // // task8->SetPtBins(dPtBinning,iSizePt);
 // // process->AddTask(task8);
 // // //
 // //
 // // FlowTask* task9 = new FlowTask("K0s",FlowTask::kK0s);
 // // task9->SetHarmonics(2);
 // // task9->SetEtaGap(0.8);
 // // // task9->SetPtBins(dPtBinning,iSizePt);
 // // task9->SuggestPtBinning(1,150000);
 // // process->AddTask(task9);
 // // //
 // FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
 // taskLambda->SetHarmonics(2);
 // taskLambda->SetEtaGap(0.8);
 // taskLambda->SetInvMassRebin(2);
 // taskLambda->SetFlowMassRebin(2);
 // taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
 // process->AddTask(taskLambda);
 // //
 // // FlowTask* task5 = new FlowTask("Lambda",FlowTask::kLambda);
 // // task5->SetHarmonics(2);
 // // task5->SetEtaGap(-1.);
 // // task5->SetPtBins(dPtBinning,iSizePt);
 // // process->AddTask(task5);
 // //
 // // FlowTask* task6 = new FlowTask("Lambda",FlowTask::kLambda);
 // // task6->SetHarmonics(2);
 // // task6->SetEtaGap(0.8);
 // // // task6->SetPtBins(dPtBinning,iSizePt);
 // // task6->SuggestPtBinning(1,150000);
 // // process->AddTask(task6);
 // //
 //
 //
 // // Double_t dPtBinningPhi[] = {0.5,0.8,1.1,1.4,1.7,2.,2.3,2.6,2.9,3.2,3.5,4,4.5,5.};
 // // Double_t dPtBinningPhi[] = {1.,1.5,2.05};
 // // Int_t iSizePtPhi = sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]);
 //
 //
 // FlowTask* task2 = new FlowTask("Phi",FlowTask::kPhi);
 // task2->SetHarmonics(2);
 // task2->SetEtaGap(0.8);
 // task2->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
 // process->AddTask(task2);
 //
 // // FlowTask* task3 = new FlowTask("Phi",FlowTask::kPhi);
 // // task3->SetHarmonics(2);
 // // task3->SetEtaGap(0.8);
 // // task3->SetPtBins(dPtBinningPhi,iSizePtPhi);
 // // task3->SuggestPtBinning(1,20000);
 // // task3->SetInvMassRebin(2);
 // // task3->SetFlowMassRebin(2);
 // // process->AddTask(task3);
 // //
 // // FlowTask* task4 = new FlowTask("Phi",FlowTask::kPhi);
 // // task4->SetHarmonics(2);
 // // task4->SetEtaGap(0.8);
 // // task4->SetPtBins(dPtBinningPhi,iSizePtPhi);
 // // task4->SetShowMultDist(kTRUE);
 // // 	process->AddTask(task4);

	process->Run();


	return;
}
