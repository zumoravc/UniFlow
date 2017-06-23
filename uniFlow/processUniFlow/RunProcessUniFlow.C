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

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt/");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst/baseline_noNUA/merged_16q");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst/Vtx_z/");
	// process->SetInputFileName("AnalysisResults_CENTwSDD_16q.root");
	// process->SetTaskName("UniFlow_V0A");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	// process->SetOutputFilePath(sOutputFilePath);
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt/fittingNotFixed/");
	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst/Vtx_z/syst/vtx9");
	process->SetOutputFileName("UniFlow.root");
	// process->SetOutputFileMode("UPDATE");
	process->SetDebug(kFALSE);
	// process->SuggestMultBinning(4);

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	Double_t dMultBinning[] = {0,20,40,60,100};
	// Double_t dMultBinning[] = {0,10};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	// Double_t dPt[] = {1.,1.2,1.4,1.6};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPtKaon[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPtProton[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
	Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
	Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1

	Double_t dPtBinningK0s[] = {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.};
	Double_t dPtBinningLambda[] = {0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.};
	Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,8.};

	// // Double_t dPt[] = {1.,1.2};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// // Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPtKaon[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPtProton[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// // Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
	// // Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
	// // Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
	//
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningLambda[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,10.,20.};

	FlowTask* taskRefs = new FlowTask("Refs",FlowTask::kRefs);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(0.8);
	taskRefs->SetNumSamples(10);
	taskRefs->SetMergePosNeg();
	process->AddTask(taskRefs);

	FlowTask* taskCharged = new FlowTask("Charged",FlowTask::kCharged);
	taskCharged->SetHarmonics(2);
	taskCharged->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskCharged->SetEtaGap(0.8);
	// taskCharged->SetAlternativeProfileName("fp2Charged_<2>_harm2_gap08_Pos");
	taskCharged->SetMergePosNeg();
	taskCharged->SetNumSamples(10);
	// process->AddTask(taskCharged);

	FlowTask* taskPion = new FlowTask("Pi",FlowTask::kPion);
	taskPion->SetHarmonics(2);
	taskPion->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskPion->SetEtaGap(0.8);
	taskPion->SetMergePosNeg();
	taskPion->SetNumSamples(10);
	// process->AddTask(taskPion);

	FlowTask* taskKaon = new FlowTask("K",FlowTask::kKaon);
	taskKaon->SetHarmonics(2);
	taskKaon->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	taskKaon->SetEtaGap(0.8);
	taskKaon->SetMergePosNeg();
	taskKaon->SetNumSamples(10);
	// process->AddTask(taskKaon);

	FlowTask* taskProton = new FlowTask("P",FlowTask::kProton);
	taskProton->SetHarmonics(2);
	taskProton->SetPtBins(dPtProton,sizeof(dPtProton)/sizeof(dPtProton[0]));
	// taskProton->SetAlternativeProfileName("fp2Proton_<2>_harm2_gap08_Pos");
	taskProton->SetEtaGap(0.8);
	taskProton->SetMergePosNeg();
	taskProton->SetNumSamples(10);
	// process->AddTask(taskProton);

	FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(0.8);
	taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	// task7->SetShowMultDist(kTRUE);
	// taskK0s->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
	// taskK0s->SetFittingRange(0.46,0.56);
	// taskK0s->SetFittingRejectNumSigmas(5);
	taskK0s->SetMergePosNeg();
	// task7->SuggestPtBinning(1,30000);
	// process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(0.8);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	// taskLambda->SetFittingRange(1.1,1.14);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Pos");
	taskLambda->SetFittingRejectNumSigmas(4);
	taskLambda->SetMergePosNeg();
	// process->AddTask(taskLambda);


	FlowTask* taskPhi = new FlowTask("Phi",FlowTask::kPhi);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(0.8);
	taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	// taskPhi->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Pos");
	taskPhi->SetFittingRange(1.,1.05);
	taskPhi->SetMergePosNeg();
	// process->AddTask(taskPhi);

	FlowTask* K0sFitting = new FlowTask("K0s",FlowTask::kK0s);
	K0sFitting->SetHarmonics(2);
	K0sFitting->SetEtaGap(0.8);
	K0sFitting->SetInvMassRebin(2);
	K0sFitting->SetFlowMassRebin(2);
	// task7->SetShowMultDist(kTRUE);
	// K0sFitting->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	K0sFitting->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	// K0sFitting->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
	// K0sFitting->SetFittingRange(0.46,0.56);
	// K0sFitting->SetFittingRejectNumSigmas(5);
	K0sFitting->SetMergePosNeg();
	// task7->SuggestPtBinning(1,30000);
	// process->AddTask(K0sFitting);

	FlowTask* K0sFittingNoFix = new FlowTask("K0s",FlowTask::kK0s);
	K0sFittingNoFix->SetHarmonics(2);
	K0sFittingNoFix->SetEtaGap(0.8);
	K0sFittingNoFix->SetInvMassRebin(2);
	K0sFittingNoFix->SetFlowMassRebin(2);
	K0sFittingNoFix->SetfFlowFitFixTerms(kFALSE);
	// task7->SetShowMultDist(kTRUE);
	// K0sFittingNoFix->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	K0sFittingNoFix->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	// K0sFittingNoFix->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
	// K0sFittingNoFix->SetFittingRange(0.46,0.56);
	K0sFittingNoFix->SetFittingRejectNumSigmas(5);
	K0sFittingNoFix->SetMergePosNeg();
	// task7->SuggestPtBinning(1,30000);
	process->AddTask(K0sFittingNoFix);






	process->Run();


	return;
}
