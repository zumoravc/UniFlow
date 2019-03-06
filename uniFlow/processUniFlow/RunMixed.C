// NEED to run .X Load.C beforehand

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunMixed()
{

	Int_t iNumSamples = 1;
	Int_t iNumSamplesRefs = 5;
	Int_t iHarmonics = 2;
	// Double_t dGap = -1.0; TString sGap = "";
	Double_t dGap = 0.0; TString sGap = "_2sub(0)";
	// Double_t dGap = 0.4; TString sGap = "_2sub(0.4)";

	TString sInputPath = "/mnt/CodesALICE/Flow/uniFlow/results/trains/CF_PbPb/6527_20190218-2140/merge/";

	// // ====== Starting points
	// v2{4}
	TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/trains/CF_PbPb/6527_20190218-2140/output_mixed/gap%02g/K0s-2/",10*dGap);
	// TString sOutputFilePath = sInputPath+Form("output_mixed/gap%02g/K0s/",10*dGap);
	std::vector<Double_t> dMultBinning = {0,5,10,20,30,40,50,60};
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
	process->SetSaveInterSteps(0);
	process->SetFitCumulants(kFALSE);
	process->SetDebug(0);

	FlowTask* taskCharged = new FlowTask(kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetEtaGap(dGap);
	taskCharged->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()));
	taskCharged->SetPtBins(dPtBinningPID);
	taskCharged->SetMergePosNeg();

	FlowTask* taskCharged2 = new FlowTask(kCharged);
	taskCharged2->SetNumSamples(iNumSamples);
	taskCharged2->SetEtaGap(dGap);
	taskCharged2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()));
	taskCharged2->SetPtBins(dPtBinningPID);
	taskCharged2->SetMergePosNeg();

	FlowTask* taskCharged3 = new FlowTask(kCharged);
	taskCharged3->SetNumSamples(iNumSamples);
	taskCharged3->SetEtaGap(dGap);
	taskCharged3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()));
	taskCharged3->SetPtBins(dPtBinningPID);
	taskCharged3->SetMergePosNeg();

	FlowTask* taskPion = new FlowTask(kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetEtaGap(dGap);
	taskPion->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()));
	taskPion->SetPtBins(dPtBinningPID);
	taskPion->SetMergePosNeg();

	FlowTask* taskPion2 = new FlowTask(kPion);
	taskPion2->SetNumSamples(iNumSamples);
	taskPion2->SetEtaGap(dGap);
	taskPion2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()));
	taskPion2->SetPtBins(dPtBinningPID);
	taskPion2->SetMergePosNeg();
	FlowTask* taskPion3 = new FlowTask(kPion);
	taskPion3->SetNumSamples(iNumSamples);
	taskPion3->SetEtaGap(dGap);
	taskPion3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()));
	taskPion3->SetPtBins(dPtBinningPID);
	taskPion3->SetMergePosNeg();

	FlowTask* taskKaon = new FlowTask(kKaon);
	taskKaon->SetNumSamples(iNumSamples);
	taskKaon->SetEtaGap(dGap);
	taskKaon->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()));
	taskKaon->SetPtBins(dPtBinningPID);
	taskKaon->SetMergePosNeg();
	FlowTask* taskKaon2 = new FlowTask(kKaon);
	taskKaon2->SetNumSamples(iNumSamples);
	taskKaon2->SetEtaGap(dGap);
	taskKaon2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()));
	taskKaon2->SetPtBins(dPtBinningPID);
	taskKaon2->SetMergePosNeg();
	FlowTask* taskKaon3 = new FlowTask(kKaon);
	taskKaon3->SetNumSamples(iNumSamples);
	taskKaon3->SetEtaGap(dGap);
	taskKaon3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()));
	taskKaon3->SetPtBins(dPtBinningPID);
	taskKaon3->SetMergePosNeg();

	FlowTask* taskProton = new FlowTask(kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetEtaGap(dGap);
	taskProton->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()));
	taskProton->SetPtBins(dPtBinningPID);
	taskProton->SetMergePosNeg();
	FlowTask* taskProton2 = new FlowTask(kProton);
	taskProton2->SetNumSamples(iNumSamples);
	taskProton2->SetEtaGap(dGap);
	taskProton2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()));
	taskProton2->SetPtBins(dPtBinningPID);
	taskProton2->SetMergePosNeg();
	FlowTask* taskProton3 = new FlowTask(kProton);
	taskProton3->SetNumSamples(iNumSamples);
	taskProton3->SetEtaGap(dGap);
	taskProton3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()));
	taskProton3->SetPtBins(dPtBinningPID);
	taskProton3->SetMergePosNeg();

	FlowTask* taskK0s = new FlowTask(kK0s);
	taskK0s->SetNumSamples(1);
	taskK0s->SetEtaGap(dGap);
	taskK0s->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()),iNumSamplesRefs);
	taskK0s->SetPtBins(dPtBinningK0s);
	taskK0s->SetMergePosNeg();
	taskK0s->SetFlowMassRebin(2);

	FlowTask* taskK0s2 = new FlowTask(kK0s);
	taskK0s2->SetNumSamples(1);
	taskK0s2->SetEtaGap(dGap);
	taskK0s2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskK0s2->SetPtBins(dPtBinningK0s);
	taskK0s2->SetMergePosNeg();
	taskK0s2->SetFlowMassRebin(2);

	FlowTask* taskK0s3 = new FlowTask(kK0s);
	taskK0s3->SetNumSamples(1);
	taskK0s3->SetEtaGap(dGap);
	taskK0s3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskK0s3->SetPtBins(dPtBinningK0s);
	taskK0s3->SetMergePosNeg();
	taskK0s3->SetFlowMassRebin(2);

	FlowTask* taskLambda = new FlowTask(kLambda);
	taskLambda->SetNumSamples(1);
	taskLambda->SetEtaGap(dGap);
	taskLambda->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()),iNumSamplesRefs);
	taskLambda->SetPtBins(dPtBinningLambda);
	taskLambda->SetMergePosNeg();
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetInvMassRebin(2);

	FlowTask* taskLambda2 = new FlowTask(kLambda);
	taskLambda2->SetNumSamples(1);
	taskLambda2->SetEtaGap(dGap);
	taskLambda2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskLambda2->SetPtBins(dPtBinningLambda);
	taskLambda2->SetMergePosNeg();
	taskLambda2->SetFlowMassRebin(2);
	taskLambda2->SetInvMassRebin(2);

	FlowTask* taskLambda3 = new FlowTask(kLambda);
	taskLambda3->SetNumSamples(1);
	taskLambda3->SetEtaGap(dGap);
	taskLambda3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskLambda3->SetPtBins(dPtBinningLambda);
	taskLambda3->SetMergePosNeg();
	taskLambda3->SetFlowMassRebin(2);
	taskLambda3->SetInvMassRebin(2);

	FlowTask* taskPhi = new FlowTask(kPhi);
	taskPhi->SetNumSamples(1);
	taskPhi->SetEtaGap(dGap);
	taskPhi->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()),iNumSamplesRefs);
	taskPhi->SetPtBins(dPtBinningPhi);
	taskPhi->SetMergePosNeg();

	FlowTask* taskPhi2 = new FlowTask(kPhi);
	taskPhi2->SetNumSamples(1);
	taskPhi2->SetEtaGap(dGap);
	taskPhi2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskPhi2->SetPtBins(dPtBinningPhi);
	taskPhi2->SetMergePosNeg();

	FlowTask* taskPhi3 = new FlowTask(kPhi);
	taskPhi3->SetNumSamples(1);
	taskPhi3->SetEtaGap(dGap);
	taskPhi3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskPhi3->SetPtBins(dPtBinningPhi);
	taskPhi3->SetMergePosNeg();


	// process->AddTask(taskCharged);
	// process->AddTask(taskCharged2);
	// process->AddTask(taskCharged3);
	// process->AddTask(taskPion);
	// process->AddTask(taskPion2);
	// process->AddTask(taskPion3);
	// process->AddTask(taskKaon);
	// process->AddTask(taskKaon2);
	// process->AddTask(taskKaon3);
	// process->AddTask(taskProton);
	// process->AddTask(taskProton2);
	// process->AddTask(taskProton3);
	process->AddTask(taskK0s);
	// process->AddTask(taskK0s2);
	// process->AddTask(taskK0s3);
	// process->AddTask(taskLambda);
	// process->AddTask(taskLambda2);
	// process->AddTask(taskLambda3);
	// process->AddTask(taskPhi);
	// process->AddTask(taskPhi2);
	// process->AddTask(taskPhi3);

	process->Run();

	return;
}
