/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunProcess_pPb_per_Nch()
{
	Int_t iNumSamples = 1;
	// Double_t dEtaGap = 0.0;	TString sEtaGap = "gap00";
	Double_t dEtaGap = 0.4;	TString sEtaGap = "gap04";
	// Double_t dEtaGap = 0.8;	TString sEtaGap = "gap08";

	TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/etc/pPb-Nch-16q-FAST-nua/";
	TString sOutputFilePath = sInputPath+"/output_binning-3-Nch/"+sEtaGap+"/";

	TString sInputPathPP = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/etc/pPb-Nch-16q-FAST-nua/";
	TString sOutputFilePathPP = sInputPathPP+"/output_binning-3-Nch_020/"+sEtaGap+"/";

	Double_t dMultBinningPP[] = {0,10,20};
	Double_t dMultBinning[] = {0,10,20,40,60,80,100,120};

	Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,7.0};
	Double_t dPtBinningProton[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,7.0};
	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.0,5.0,7.0};
	Double_t dPtBinningLambda[] = {0.3,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.0,5.0,7.0};
	Double_t dPtBinningPhi[] = {0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0};

	// ##### END Parameters setting ######

	// gROOT->AddIncludePath("-I~/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~vpacik/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");

	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(dEtaGap);
	taskRefs->SetMergePosNeg();

	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(dEtaGap);
	taskCharged->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskCharged->SetMergePosNeg();

	FlowTask* taskPion = new FlowTask(FlowTask::kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetHarmonics(2);
	taskPion->SetEtaGap(dEtaGap);
	taskPion->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskPion->SetMergePosNeg();

	FlowTask* taskKch = new FlowTask(FlowTask::kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(2);
	taskKch->SetEtaGap(dEtaGap);
	taskKch->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskProton = new FlowTask(FlowTask::kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetHarmonics(2);
	taskProton->SetEtaGap(dEtaGap);
	taskProton->SetPtBins(dPtBinningProton,sizeof(dPtBinningProton)/sizeof(dPtBinningProton[0]));
	taskProton->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	taskK0s->SetNumSamples(iNumSamples);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(dEtaGap);
	// taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskLambda = new FlowTask(FlowTask::kLambda	);
	taskLambda->SetNumSamples(iNumSamples);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(dEtaGap);
	// taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	taskLambda->SetMergePosNeg();
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskPhi = new FlowTask(FlowTask::kPhi);
	taskPhi->SetNumSamples(iNumSamples);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(dEtaGap);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskPhi->SetMergePosNeg();
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");



	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->SetSaveMult(kTRUE);
	process->SetFitCumulants(kTRUE);
	process->SetDebug(kFALSE);
	process->AddTask(taskRefs);
	process->AddTask(taskCharged);
	// process->AddTask(taskPion);
	// process->AddTask(taskKch);
	// process->AddTask(taskProton);
	// process->AddTask(taskK0s);
	// process->AddTask(taskLambda);
	// process->AddTask(taskPhi);
	process->Run();

	ProcessUniFlow* processPP = new ProcessUniFlow();
	processPP->SetInputFilePath(sInputPathPP.Data());
	processPP->SetInputFileName("AnalysisResults.root");
	processPP->SetTaskName("UniFlow");
	processPP->SetOutputFilePath(sOutputFilePathPP.Data());
	processPP->SetOutputFileName("Processed.root");
	processPP->SetMultiplicityBins(dMultBinningPP,sizeof(dMultBinningPP)/sizeof(dMultBinningPP[0]));
	processPP->SetSaveMult(kTRUE);
	processPP->SetFitCumulants(kTRUE);
	processPP->SetDebug(kFALSE);
	processPP->AddTask(taskRefs);
	processPP->AddTask(taskCharged);
	// processPP->AddTask(taskPion);
	// processPP->AddTask(taskKch);
	// processPP->AddTask(taskProton);
	// processPP->AddTask(taskK0s);
	// processPP->AddTask(taskLambda);
	// processPP->AddTask(taskPhi);
	processPP->Run();

	return;
}
