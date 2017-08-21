/* RunProcessUniFlow
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunFitTesting(const char* sOutputFilePath = "")
{
	// gSystem->AddIncludePath("-I/Users/vpacik/NBI/triggerHMstudies/TrigEff");
	// gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/triggerHMstudies/TrigEff",gROOT->GetMacroPath()));
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	ProcessUniFlow* process = new ProcessUniFlow();

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/loosePhi_tightV0s");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit_withNUA/");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit");

	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	// process->SetTaskName("uniflow_Alex");
	// process->SetTaskName("uniflow_dcaZ");
	// process->SetTaskName("UniFlow_fb768");
	// process->SetOutputFilePath(sOutputFilePath);

	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/loosePhi_tightV0s/kaon");
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTest_final_plots");
	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/PbPblike");
	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit_withNUA/UniFlow_fb768_kaons");
	process->SetOutputFileName("UniFlow_test.root");
	// process->SetOutputFileMode("UPDATE");
	process->SetDebug(kFALSE);
	// process->SetDebug(kTRUE);
	// process->SuggestMultBinning(4);

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	// Double_t dMultBinning[] = {0,20,40,60,100};
	Double_t dMultBinning[] = {0,20};
	// Double_t dMultBinning[] = {40,60};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	// Double_t dPt[] = {1.,1.2,1.4,1.6};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPtKaon[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPtProton[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.};

	// Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
	Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
	// Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1


	// Double_t dPtBinningK0s[] = {2.,4.};
	// Double_t dPtBinningK0s[] = {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.6,3.,3.5,4.,5.,6.};
	// Double_t dPtBinningLambda[] = {0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.7,4.5,6.};
	Double_t dPtBinningK0s[] = {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.};
	// Double_t dPtBinningLambda[] = {0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.8,4.6,6.};
	// Double_t dPtBinningLambda[] = {0.,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.6.,4.3.,5.295};
	Double_t dPtBinningLambda[] = {0.8,1.1,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.6.,4.3.,5.295};
	Double_t dPtBinningPhi[] = {1.5,2.};
	// Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.};

	FlowTask* taskRefs = new FlowTask("Refs",FlowTask::kRefs);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(0.8);
	taskRefs->SetNumSamples(10);
	taskRefs->SetMergePosNeg(1);
	process->AddTask(taskRefs);

	FlowTask* taskCharged = new FlowTask("Charged",FlowTask::kCharged);
	taskCharged->SetHarmonics(2);
	taskCharged->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskCharged->SetEtaGap(0.8);
	// taskCharged->SetAlternativeProfileName("fp2Charged_<2>_harm2_gap08_Pos");
	taskCharged->SetMergePosNeg(1);
	taskCharged->SetNumSamples(10);
	// process->AddTask(taskCharged);

	FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(0.8);
	taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	taskK0s->SetFlowFitFixTerms(0);
	// task7->SetShowMultDist(kTRUE);
	// taskK0s->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	// taskK0s->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
	taskK0s->SetFittingRejectNumSigmas(7);
	taskK0s->SetFittingRange(0.44,0.56);
	// taskK0s->SetFittingRange(0.42,0.58);
	taskK0s->SetMergePosNeg(1);
	// task7->SuggestPtBinning(1,30000);
	// process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(0.8);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	// taskLambda->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	taskLambda->SetFittingRange(1.099,1.147);
	// taskLambda->SetFittingRange(1.095,1.155);
	// taskLambda->SetFittingRange(1.101,1.149);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Pos");
	taskLambda->SetFittingRejectNumSigmas(9);
	taskLambda->SetMergePosNeg();
	taskLambda->SetFlowFitFixTerms(0);
	// process->AddTask(taskLambda);


	FlowTask* taskPhi = new FlowTask("Phi",FlowTask::kPhi);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(0.8);
	taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	// taskPhi->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Pos");
	// taskPhi->SetFittingRange(0.999,1.065);
	taskPhi->SetFittingRange(1.001,1.061);
	taskPhi->SetMergePosNeg();
	taskPhi->SetFlowFitFixTerms(0);
	taskPhi->SetFlowPhiSubtLS(0);
	process->AddTask(taskPhi);

	FlowTask* taskKaon = new FlowTask("K",FlowTask::kKaon);
	taskKaon->SetHarmonics(2);
	// taskKaon->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskKaon->SetEtaGap(0.8);
	taskKaon->SetMergePosNeg(1);
	taskKaon->SetNumSamples(10);
	// process->AddTask(taskKaon);

	process->Run();


	return;
}
//
