// NEED to run .X Load.C beforehand

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunCumulants()
{

	Int_t iNumSamples = 1;
	Int_t iHarmonics = 2;
	TString sInputPath = "~/Codes/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/";
	// TString sOutputFilePath = sInputPath+"/output/";


	// // ====== Starting points
	// v2{4}
	TString sOutputFilePath = sInputPath+"/output_cum/";
	std::vector<Double_t> dMultBinning = {0,5,10,20,30,40,50};
	std::vector<Double_t> dPtBinningPID = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.4,4.8,5.2,5.6,6.0};
	std::vector<Double_t> dPtBinningProton = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0};
	std::vector<Double_t> dPtBinningK0s = {0.4,0.8,1.2,1.6,2.0,2.4,3.2,4.0,5.0,6.0};
	std::vector<Double_t> dPtBinningLambda = {0.8,1.2,1.6,2.0,2.4,2.8,3.2,4.0,5.0,6.0};
	std::vector<Double_t> dPtBinningPhi = {1.0,1.5,2.0,2.5,3.0,4.0,5.0};

	// v3{4}
	// TString sOutputFilePath = sInputPath+"/output_v3/";
	// Double_t dMultBinning[] = {10,20,30,40,50};
	// Double_t dPtBinningPID[] = {1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	// Double_t dPtBinningProton[] = {1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	// Double_t dPtBinningK0s[] = {1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	// Double_t dPtBinningLambda[] = {1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	// Double_t dPtBinningPhi[] = {1.0,1.5,2.0,2.5,3.0,4.0,6.0};


	// ##### END Parameters setting ######

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning);
	process->SetSaveMult(0);
	process->SetFitCumulants(kFALSE);
	process->SetDebug(1);

	FlowTask* taskRefs = new FlowTask(kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	taskRefs->SetHarmonics(iHarmonics);
	taskRefs->DoCumTwo(1);
	taskRefs->DoCumFour(1);
	taskRefs->SetEtaGap(-1.0);
	taskRefs->SetMergePosNeg();
	process->AddTask(taskRefs);

	FlowTask* taskRefs2 = new FlowTask(kRefs);
	taskRefs2->SetNumSamples(iNumSamples);
	taskRefs2->SetHarmonics(iHarmonics);
	taskRefs2->DoCumTwo(1);
	taskRefs2->DoCumFour(1);
	taskRefs2->SetEtaGap(0.0);
	taskRefs2->SetMergePosNeg();
	process->AddTask(taskRefs2);

	FlowTask* taskRefs3 = new FlowTask(kRefs);
	taskRefs3->SetNumSamples(iNumSamples);
	taskRefs3->SetHarmonics(iHarmonics);
	taskRefs3->DoCumTwo(1);
	taskRefs3->DoCumFour(1);
	taskRefs3->SetEtaGap(0.4);
	taskRefs3->SetMergePosNeg();
	process->AddTask(taskRefs3);

	FlowTask* taskCharged = new FlowTask(kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetHarmonics(iHarmonics);
	taskCharged->DoCumTwo(1);
	taskCharged->DoCumFour(1);
	taskCharged->SetEtaGap(-1.0);
	taskCharged->SetPtBins(dPtBinningPID);
	taskCharged->SetMergePosNeg();
	process->AddTask(taskCharged);

	FlowTask* taskCharged2 = new FlowTask(kCharged);
	taskCharged2->SetNumSamples(iNumSamples);
	taskCharged2->SetHarmonics(iHarmonics);
	taskCharged2->DoCumTwo(1);
	taskCharged2->DoCumFour(1);
	taskCharged2->SetEtaGap(0.0);
	taskCharged2->SetPtBins(dPtBinningPID);
	taskCharged2->SetMergePosNeg();
	process->AddTask(taskCharged2);

	FlowTask* taskCharged3 = new FlowTask(kCharged);
	taskCharged3->SetNumSamples(iNumSamples);
	taskCharged3->SetHarmonics(iHarmonics);
	taskCharged3->DoCumTwo(1);
	taskCharged3->DoCumFour(1);
	taskCharged3->SetEtaGap(0.4);
	taskCharged3->SetPtBins(dPtBinningPID);
	taskCharged3->SetMergePosNeg();
	process->AddTask(taskCharged3);

	FlowTask* taskPion = new FlowTask(kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetHarmonics(iHarmonics);
	taskPion->DoCumTwo(1);
	taskPion->DoCumFour(1);
	taskPion->SetEtaGap(-1.0);
	taskPion->SetPtBins(dPtBinningPID);
	taskPion->SetMergePosNeg();
	process->AddTask(taskPion);

	FlowTask* taskPion2 = new FlowTask(kPion);
	taskPion2->SetNumSamples(iNumSamples);
	taskPion2->SetHarmonics(iHarmonics);
	taskPion2->DoCumTwo(1);
	taskPion2->DoCumFour(1);
	taskPion2->SetEtaGap(0.0);
	taskPion2->SetPtBins(dPtBinningPID);
	taskPion2->SetMergePosNeg();
	process->AddTask(taskPion2);

	FlowTask* taskPion3 = new FlowTask(kPion);
	taskPion3->SetNumSamples(iNumSamples);
	taskPion3->SetHarmonics(iHarmonics);
	taskPion3->DoCumTwo(1);
	taskPion3->DoCumFour(1);
	taskPion3->SetEtaGap(0.4);
	taskPion3->SetPtBins(dPtBinningPID);
	taskPion3->SetMergePosNeg();
	process->AddTask(taskPion3);

	FlowTask* taskKch = new FlowTask(kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(iHarmonics);
	taskKch->DoCumTwo(1);
	taskKch->DoCumFour(1);
	taskKch->SetEtaGap(-1.0);
	taskKch->SetPtBins(dPtBinningPID);
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskKch);

	FlowTask* taskKch2 = new FlowTask(kKaon);
	taskKch2->SetNumSamples(iNumSamples);
	taskKch2->SetHarmonics(iHarmonics);
	taskKch2->DoCumTwo(1);
	taskKch2->DoCumFour(1);
	taskKch2->SetEtaGap(0.0);
	taskKch2->SetPtBins(dPtBinningPID);
	taskKch2->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskKch2);

	FlowTask* taskKch3 = new FlowTask(kKaon);
	taskKch3->SetNumSamples(iNumSamples);
	taskKch3->SetHarmonics(iHarmonics);
	taskKch3->DoCumTwo(1);
	taskKch3->DoCumFour(1);
	taskKch3->SetEtaGap(0.4);
	taskKch3->SetPtBins(dPtBinningPID);
	taskKch3->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskKch3);

	FlowTask* taskProton = new FlowTask(kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetHarmonics(iHarmonics);
	taskProton->DoCumTwo(1);
	taskProton->DoCumFour(1);
	taskProton->SetEtaGap(-1.0);
	taskProton->SetPtBins(dPtBinningProton);
	taskProton->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskProton);

	FlowTask* taskProton2 = new FlowTask(kProton);
	taskProton2->SetNumSamples(iNumSamples);
	taskProton2->SetHarmonics(iHarmonics);
	taskProton2->DoCumTwo(1);
	taskProton2->DoCumFour(1);
	taskProton2->SetEtaGap(0.0);
	taskProton2->SetPtBins(dPtBinningProton);
	taskProton2->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskProton2);

	FlowTask* taskProton3 = new FlowTask(kProton);
	taskProton3->SetNumSamples(iNumSamples);
	taskProton3->SetHarmonics(iHarmonics);
	taskProton3->DoCumTwo(1);
	taskProton3->DoCumFour(1);
	taskProton3->SetEtaGap(0.4);
	taskProton3->SetPtBins(dPtBinningProton);
	taskProton3->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskProton3);


 	FlowTask* taskK0s = new FlowTask(kK0s);
	taskK0s->SetNumSamples(iNumSamples);
	taskK0s->SetHarmonics(iHarmonics);
	taskK0s->DoCumTwo(1);
	taskK0s->DoCumFour(1);
	taskK0s->SetEtaGap(-1.0);
	// taskK0s->SetInvMassRebin(2);
	// taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s);
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskK0s);

	FlowTask* taskK0s2 = new FlowTask(kK0s);
	taskK0s2->SetNumSamples(iNumSamples);
	taskK0s2->SetHarmonics(iHarmonics);
	taskK0s2->DoCumTwo(1);
	taskK0s2->DoCumFour(1);
	taskK0s2->SetEtaGap(0.0);
	// taskK0s2->SetInvMassRebin(2);
	// taskK0s2->SetFlowMassRebin(2);
	taskK0s2->SetPtBins(dPtBinningK0s);
	taskK0s2->SetMergePosNeg();
	// taskK0s2->SetFittingRange(0.45,0.55);
	// taskK0s2->SetFittingRejectNumSigmas(3);
	// taskK0s2->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskK0s2);

	FlowTask* taskK0s3 = new FlowTask(kK0s);
	taskK0s3->SetNumSamples(iNumSamples);
	taskK0s3->SetHarmonics(iHarmonics);
	taskK0s3->DoCumTwo(1);
	taskK0s3->DoCumFour(1);
	taskK0s3->SetEtaGap(0.4);
	// taskK0s2->SetInvMassRebin(2);
	// taskK0s2->SetFlowMassRebin(2);
	taskK0s3->SetPtBins(dPtBinningK0s);
	taskK0s3->SetMergePosNeg();
	// taskK0s2->SetFittingRange(0.45,0.55);
	// taskK0s2->SetFittingRejectNumSigmas(3);
	// taskK0s2->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskK0s3);

	FlowTask* taskLambda = new FlowTask(kLambda);
	taskLambda->SetNumSamples(iNumSamples);
	taskLambda->SetHarmonics(iHarmonics);
	taskLambda->SetEtaGap(-1.0);
	taskLambda->DoCumTwo(1);
	taskLambda->DoCumFour(1);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda);
	taskLambda->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskLambda);

	FlowTask* taskLambda2 = new FlowTask(kLambda);
	taskLambda2->SetNumSamples(iNumSamples);
	taskLambda2->SetHarmonics(iHarmonics);
	taskLambda2->SetEtaGap(0.0);
	taskLambda2->DoCumTwo(1);
	taskLambda2->DoCumFour(1);
	taskLambda2->SetInvMassRebin(2);
	taskLambda2->SetFlowMassRebin(2);
	taskLambda2->SetPtBins(dPtBinningLambda);
	taskLambda2->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskLambda2);

	FlowTask* taskLambda3 = new FlowTask(kLambda);
	taskLambda3->SetNumSamples(iNumSamples);
	taskLambda3->SetHarmonics(iHarmonics);
	taskLambda3->SetEtaGap(0.4);
	taskLambda3->DoCumTwo(1);
	taskLambda3->DoCumFour(1);
	taskLambda3->SetInvMassRebin(2);
	taskLambda3->SetFlowMassRebin(2);
	taskLambda3->SetPtBins(dPtBinningLambda);
	taskLambda3->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskLambda3);

	FlowTask* taskPhi = new FlowTask(kPhi);
	taskPhi->SetNumSamples(iNumSamples);
	taskPhi->SetHarmonics(iHarmonics);
	taskPhi->DoCumTwo(1);
	taskPhi->DoCumFour(1);
	taskPhi->SetEtaGap(-1.0);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi);
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskPhi);

	FlowTask* taskPhi2 = new FlowTask(kPhi);
	taskPhi2->SetNumSamples(iNumSamples);
	taskPhi2->SetHarmonics(iHarmonics);
	taskPhi2->DoCumTwo(1);
	taskPhi2->DoCumFour(1);
	taskPhi2->SetEtaGap(0.0);
	// taskPhi->SetInvMassRebin(2);
	taskPhi2->SetFlowMassRebin(2);
	taskPhi2->SetPtBins(dPtBinningPhi);
	taskPhi2->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskPhi2);

	FlowTask* taskPhi3 = new FlowTask(kPhi);
	taskPhi3->SetNumSamples(iNumSamples);
	taskPhi3->SetHarmonics(iHarmonics);
	taskPhi3->DoCumTwo(1);
	taskPhi3->DoCumFour(1);
	taskPhi3->SetEtaGap(0.4);
	// taskPhi->SetInvMassRebin(2);
	taskPhi3->SetFlowMassRebin(2);
	taskPhi3->SetPtBins(dPtBinningPhi);
	taskPhi3->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskPhi3);


	process->Run();

	return;
}
