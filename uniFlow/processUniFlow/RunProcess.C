#ifdef __CLING__
#include "ProcessUniFlow.cpp"
#endif

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunProcess()
{

	Int_t iNumSamples = 1;
	Double_t dEtaGap = 0.0;
	TString sEtaGap = "gap00";

	TString sInputPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/Flow_NUA/";
	TString sOutputFilePath = sInputPath+"/output/"+sEtaGap+"/";

	// // ====== Starting points

	// Double_t dMultBinning[] = {0,10,20,40,60,100};
	// Double_t dPtBinningPID[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0};
	Double_t dPtBinningProton[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0};
	// Double_t dPtBinningK0s[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0};
	// Double_t dPtBinningLambda[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0};
	// Double_t dPtBinningPhi[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0};

	// // Final binning
	// Double_t dPtBinningK0s[] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,4.0,6.0};

	// Seems rather good
	Double_t dPtBinningK0s[] = {0.4,0.8,1.2,1.6,2.0,2.4,3.2,4.0,5.0,6.0};
	// Double_t dPtBinningK0s[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.8,3.2,3.6,4.0,4.5,5.0};
	Double_t dPtBinningLambda[] = {0.8,1.2,1.6,2.0,2.4,2.8,3.2,4.0,6.0};
	Double_t dPtBinningPhi[] = {1.0,1.5,2.0,2.5,3.0,4.0,5.0};

	// // ===== Naghmeg binning of PID particles
	Double_t dMultBinning[] = {0,5,10,20,30,40,50,60};
	Double_t dPtBinningPID[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0};


	// ##### END Parameters setting ######

	#if defined (__CINT__)
		gROOT->LoadMacro("~vpacik/Codes/ALICE/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");
	#endif

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->SetSaveMult(0);
	process->SetFitCumulants(kFALSE);
	process->SetDebug(1);

	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	taskRefs->SetHarmonics(2);
	taskRefs->SetDoFourCorrelations(1);
	taskRefs->SetEtaGap(dEtaGap);
	taskRefs->SetMergePosNeg();
	process->AddTask(taskRefs);

	FlowTask* taskRefs2 = new FlowTask(FlowTask::kRefs);
	taskRefs2->SetNumSamples(iNumSamples);
	taskRefs2->SetHarmonics(2);
	taskRefs2->SetDoFourCorrelations(1);
	taskRefs2->SetEtaGap(-1.0);
	taskRefs2->SetMergePosNeg();
	process->AddTask(taskRefs2);

	FlowTask* taskRefs3 = new FlowTask(FlowTask::kRefs);
	taskRefs3->SetNumSamples(iNumSamples);
	taskRefs3->SetHarmonics(2);
	taskRefs3->SetDoFourCorrelations(1);
	taskRefs3->SetEtaGap(0.4);
	taskRefs3->SetMergePosNeg();
	process->AddTask(taskRefs3);

	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetHarmonics(2);
	taskCharged->SetDoFourCorrelations(1);
	taskCharged->SetEtaGap(dEtaGap);
	taskCharged->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskCharged->SetMergePosNeg();
	process->AddTask(taskCharged);

	FlowTask* taskCharged2 = new FlowTask(FlowTask::kCharged);
	taskCharged2->SetNumSamples(iNumSamples);
	taskCharged2->SetHarmonics(2);
	taskCharged2->SetDoFourCorrelations(1);
	taskCharged2->SetEtaGap(-1.0);
	taskCharged2->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskCharged2->SetMergePosNeg();
	process->AddTask(taskCharged2);

	FlowTask* taskCharged3 = new FlowTask(FlowTask::kCharged);
	taskCharged3->SetNumSamples(iNumSamples);
	taskCharged3->SetHarmonics(2);
	taskCharged3->SetDoFourCorrelations(1);
	taskCharged3->SetEtaGap(0.4);
	taskCharged3->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskCharged3->SetMergePosNeg();
	process->AddTask(taskCharged3);

	FlowTask* taskPion = new FlowTask(FlowTask::kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetHarmonics(2);
	taskPion->SetDoFourCorrelations(1);
	taskPion->SetEtaGap(dEtaGap);
	taskPion->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskPion->SetMergePosNeg();
	// process->AddTask(taskPion);

	FlowTask* taskKch = new FlowTask(FlowTask::kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(2);
	taskKch->SetDoFourCorrelations(1);
	taskKch->SetEtaGap(dEtaGap);
	taskKch->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	// process->AddTask(taskKch);

	FlowTask* taskProton = new FlowTask(FlowTask::kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetHarmonics(2);
	taskProton->SetDoFourCorrelations(1);
	taskProton->SetEtaGap(dEtaGap);
	taskProton->SetPtBins(dPtBinningProton,sizeof(dPtBinningProton)/sizeof(dPtBinningProton[0]));
	taskProton->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	// process->AddTask(taskProton);

 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	taskK0s->SetNumSamples(iNumSamples);
	taskK0s->SetHarmonics(2);
	taskK0s->SetDoFourCorrelations(1);
	taskK0s->SetEtaGap(dEtaGap);
	// taskK0s->SetInvMassRebin(2);
	// taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	// process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask(FlowTask::kLambda	);
	taskLambda->SetNumSamples(iNumSamples);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(dEtaGap);
	taskLambda->SetDoFourCorrelations(1);
	// taskLambda->SetInvMassRebin(2);
	// taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	taskLambda->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	// process->AddTask(taskLambda);

	FlowTask* taskPhi = new FlowTask(FlowTask::kPhi);
	taskPhi->SetNumSamples(iNumSamples);
	taskPhi->SetHarmonics(2);
	taskPhi->SetDoFourCorrelations(1);
	taskPhi->SetEtaGap(dEtaGap);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	// process->AddTask(taskPhi);


	process->Run();

	return;
}
