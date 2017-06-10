/* RunProcessUniFlow
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunProcessUniFlow(const char* sOutputFilePath = "")
{
	// gSystem->AddIncludePath("-I/Users/vpacik/NBI/triggerHMstudies/TrigEff");
	// gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/triggerHMstudies/TrigEff",gROOT->GetMacroPath()));
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	ProcessUniFlow* process = new ProcessUniFlow();

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full/CENT_woSDD_16q/");
	// process->SetInputFileName("AnalysisResults_CENTwSDD_16q.root");
	// process->SetTaskName("UniFlow_V0A");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	// process->SetOutputFilePath(sOutputFilePath);
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow");
	process->SetOutputFileName("UniFlowTest.root");
	// process->SetOutputFileMode("UPDATE");
	// process->SetDebug();
	// process->SuggestMultBinning(4);

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	Double_t dMultBinning[] = {0,20};
	// Double_t dMultBinning[] = {0,10};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	FlowTask* taskRefs = new FlowTask("Refs",FlowTask::kRefs);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(0.8);
	process->AddTask(taskRefs);

	// DIRECT

	Double_t dPt[] = {1.,1.2};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	Double_t dPtKaon[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	Double_t dPtProton[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
	// Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
	// Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1

	FlowTask* taskPion = new FlowTask("Pi",FlowTask::kPion);
	taskPion->SetHarmonics(2);
	taskPion->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskPion->SetEtaGap(0.8);
	// process->AddTask(taskPion);


	FlowTask* taskKaon = new FlowTask("K",FlowTask::kKaon);
	taskKaon->SetHarmonics(2);
	taskKaon->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	taskKaon->SetEtaGap(0.8);
	// process->AddTask(taskKaon);

	FlowTask* taskProton = new FlowTask("p",FlowTask::kProton);
	taskProton->SetHarmonics(2);
	taskProton->SetPtBins(dPtProton,sizeof(dPtProton)/sizeof(dPtProton[0]));
	// taskProton->SetAlternativeProfileName("fp2Proton_<2>_harm2_gap08_Pos");
	taskProton->SetEtaGap(0.8);
	// process->AddTask(taskProton);

	FlowTask* task5 = new FlowTask("charged_pos",FlowTask::kCharged);
	task5->SetHarmonics(2);
	task5->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	task5->SetEtaGap(0.8);
	// task5->SetAlternativeProfileName("fp2Charged_<2>_harm2_gap08_Pos");
	task5->SetMergePosNeg();
	task5->SetNumSamples(1);
	process->AddTask(task5);
//
	FlowTask* taskChargedNeg = new FlowTask("charged_neg",FlowTask::kCharged);
	taskChargedNeg->SetHarmonics(2);
	taskChargedNeg->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskChargedNeg->SetEtaGap(0.8);
	taskChargedNeg->SetAlternativeProfileName("fp2Charged_<2>_harm2_gap08_Neg");
	// process->AddTask(taskChargedNeg);

	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningLambda[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,10.,20.};


	// Double_t dPtBinning[] = {0.5,1.};
	// Double_t dPtBinning[] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.};
	// Int_t iSizePt = sizeof(dPtBinning)/sizeof(dPtBinning[0]);

 FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
 taskK0s->SetHarmonics(2);
 taskK0s->SetEtaGap(0.8);
 taskK0s->SetInvMassRebin(2);
 taskK0s->SetFlowMassRebin(2);
 // task7->SetShowMultDist(kTRUE);
 taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
 // taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
 taskK0s->SetMergePosNeg();
 // task7->SuggestPtBinning(1,30000);
 process->AddTask(taskK0s);

 FlowTask* taskK0sNeg = new FlowTask("K0sNeg",FlowTask::kK0s);
 taskK0sNeg->SetHarmonics(2);
 taskK0sNeg->SetEtaGap(0.8);
 taskK0sNeg->SetInvMassRebin(2);
 taskK0sNeg->SetFlowMassRebin(2);
 // task7->SetShowMultDist(kTRUE);
 taskK0sNeg->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
 taskK0sNeg->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
 // task7->SuggestPtBinning(1,30000);
 // process->AddTask(taskK0sNeg);
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
 // process->AddTask(task9);
 // // //
 FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
 taskLambda->SetHarmonics(2);
 taskLambda->SetEtaGap(0.8);
 taskLambda->SetInvMassRebin(2);
 taskLambda->SetFlowMassRebin(2);
 taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
 // taskLambda->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Pos");
 taskLambda->SetMergePosNeg();
 process->AddTask(taskLambda);

 FlowTask* taskLambdaNeg = new FlowTask("LambdaNeg",FlowTask::kLambda);
 taskLambdaNeg->SetHarmonics(2);
 taskLambdaNeg->SetEtaGap(0.8);
 taskLambdaNeg->SetInvMassRebin(2);
 taskLambdaNeg->SetFlowMassRebin(2);
 taskLambdaNeg->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
 taskLambdaNeg->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Neg");
 // process->AddTask(taskLambdaNeg);

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
 // Double_t dPtBinningPhi[] = {0.5,0.8,1.1,1.4,1.7,2.,2.3,2.6,2.9,3.2,3.5,4,4.5,5.};
 // // Double_t dPtBinningPhi[] = {1.,1.5,2.05};
 // // Int_t iSizePtPhi = sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]);
 //
 //
 FlowTask* taskPhi = new FlowTask("Phi",FlowTask::kPhi);
 taskPhi->SetHarmonics(2);
 taskPhi->SetEtaGap(0.8);
 taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
 // taskPhi->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Pos");
 taskPhi->SetMergePosNeg();
 // process->AddTask(taskPhi);

 FlowTask* taskPhiNeg = new FlowTask("PhiNeg",FlowTask::kPhi);
 taskPhiNeg->SetHarmonics(2);
 taskPhiNeg->SetEtaGap(0.8);
 taskPhiNeg->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
 taskPhiNeg->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Neg");
 // process->AddTask(taskPhiNeg);

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
