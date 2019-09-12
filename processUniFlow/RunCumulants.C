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
	Int_t iCumMax = 4;
	Double_t dGap = -1.0;
	Int_t iNumSamples = 5;
	Int_t iHarmonics = 2;
	TString sInputPath = "/mnt/CodesALICE/Flow/uniFlow/results/trains/CF_PbPb/6527_20190218-2140/merge/";
	// TString sOutputFilePath = sInputPath+"/output_cum/";
	TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/PbPb_cums/gap%02.f/",10*dGap);


	// // ====== Starting points
	// v2{4}
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
	process->SetDebug(0);

	FlowTask* taskRefs = new FlowTask(kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	taskRefs->SetHarmonics(iHarmonics);
	taskRefs->DoCumOrderMax(iCumMax);
	taskRefs->SetEtaGap(dGap);
	taskRefs->SetMergePosNeg();
	process->AddTask(taskRefs);

	FlowTask* taskCharged = new FlowTask(kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetHarmonics(iHarmonics);
	taskCharged->DoCumOrderMax(iCumMax);
	taskCharged->SetEtaGap(dGap);
	taskCharged->SetPtBins(dPtBinningPID);
	taskCharged->SetMergePosNeg();
	process->AddTask(taskCharged);

	FlowTask* taskPion = new FlowTask(kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetHarmonics(iHarmonics);
	taskPion->DoCumOrderMax(iCumMax);
	taskPion->SetEtaGap(dGap);
	taskPion->SetPtBins(dPtBinningPID);
	taskPion->SetMergePosNeg();
	process->AddTask(taskPion);

	FlowTask* taskKch = new FlowTask(kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(iHarmonics);
	taskKch->DoCumOrderMax(iCumMax);
	taskKch->SetEtaGap(dGap);
	taskKch->SetPtBins(dPtBinningPID);
	taskKch->SetMergePosNeg();
	process->AddTask(taskKch);

	FlowTask* taskProton = new FlowTask(kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetHarmonics(iHarmonics);
	taskProton->DoCumOrderMax(iCumMax);
	taskProton->SetEtaGap(dGap);
	taskProton->SetPtBins(dPtBinningProton);
	taskProton->SetMergePosNeg();
	process->AddTask(taskProton);

 	FlowTask* taskK0s = new FlowTask(kK0s);
	taskK0s->SetNumSamples(1);
	taskK0s->SetHarmonics(iHarmonics);
	taskK0s->DoCumOrderMax(iCumMax);
	taskK0s->SetEtaGap(dGap);
	// taskK0s->SetInvMassRebin(2);
	// taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s);
	taskK0s->SetMergePosNeg();
	process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask(kLambda);
	taskLambda->SetNumSamples(1);
	taskLambda->SetHarmonics(iHarmonics);
	taskLambda->SetEtaGap(dGap);
	taskLambda->DoCumOrderMax(iCumMax);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda);
	taskLambda->SetMergePosNeg();
	process->AddTask(taskLambda);

	FlowTask* taskPhi = new FlowTask(kPhi);
	taskPhi->SetNumSamples(1);
	taskPhi->SetHarmonics(iHarmonics);
	taskPhi->DoCumOrderMax(iCumMax);
	taskPhi->SetEtaGap(dGap);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi);
	taskPhi->SetMergePosNeg();
	process->AddTask(taskPhi);

	process->Run();

	return;
}
