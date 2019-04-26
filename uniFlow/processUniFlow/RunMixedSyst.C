// NEED to run .X Load.C beforehand

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunMixedSyst(TString tag = "V0sPVDCA3")
{
	Int_t iNumSamplesRefs = 1;

	// TString sTaskTag = "FB768";
	// TString sTaskTag = "PID3sigma";
	// TString sTaskTag = "PVz8";
	// TString sTaskTag = "TPCcls90";
	// TString sTaskTag = "V0sCPA099";
	// TString sTaskTag = "V0sCrossFind1";
	// TString sTaskTag = "V0sDaugDCA3";
	// TString sTaskTag = "V0sDaugPt02";
	// TString sTaskTag = "V0sDecRad10";
	// TString sTaskTag = "V0sFinderOn";
	// TString sTaskTag = "V0sPVDCA3";

	TString sTaskTag = tag;

	// Double_t dGap = -1.0; TString sGap = "";
	Double_t dGap = 0.0; TString sGap = "_2sub(0)";
	// Double_t dGap = 0.4; TString sGap = "_2sub(0.4)";

	TString sInputPath = "/mnt/CodesALICE/Flow/uniFlow/results/trains/CF_PbPb/6615_20190311-1326/merged/";

	// // ====== Starting points

	// std::vector<Double_t> dMultBinning = {10,20,30,40,50};
	// std::vector<Double_t> dPtBinningK0s = {0.4,0.8,1.2,1.6,2.0,2.4,3.0,3.6,4.6,6.0};
	// std::vector<Double_t> dPtBinningLambda = {0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.4,6.0};
	// std::vector<Double_t> dPtBinningPhi = {0.5,1.0,1.5,2.0,3.0,4.0,6.0};

	std::vector<Double_t> dMultBinning = {0,5,10,20,30,40,50,60};
	std::vector<Double_t> dPtBinningK0s = {0.8,1.2,1.6,2.0,2.4,3.0,3.6,4.6,6.0};
	std::vector<Double_t> dPtBinningLambda = {0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.4,6.0};
	// std::vector<Double_t> dPtBinningPhi = {1.0,2.0,3.0,4.0,5.0,6.0};
	std::vector<Double_t> dPtBinningPhi = {1.0,1.5,2.0,3.0,4.5,6.0};
	// std::vector<Double_t> dPtBinningPhi = {1.0,2.0,3.0,4.0,6.0};


	// ##### END Parameters setting ######
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
	// taskLambda->SetInvMassRebin(2);

	FlowTask* taskLambda2 = new FlowTask(kLambda);
	taskLambda2->SetNumSamples(1);
	taskLambda2->SetEtaGap(dGap);
	taskLambda2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskLambda2->SetPtBins(dPtBinningLambda);
	taskLambda2->SetMergePosNeg();
	taskLambda2->SetFlowMassRebin(3);
	// taskLambda2->SetInvMassRebin(2);

	FlowTask* taskLambda3 = new FlowTask(kLambda);
	taskLambda3->SetNumSamples(1);
	taskLambda3->SetEtaGap(dGap);
	taskLambda3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskLambda3->SetPtBins(dPtBinningLambda);
	taskLambda3->SetMergePosNeg();
	taskLambda3->SetFlowMassRebin(2);
	// taskLambda3->SetInvMassRebin(2);

	FlowTask* taskPhi = new FlowTask(kPhi);
	taskPhi->SetNumSamples(1);
	taskPhi->SetEtaGap(dGap);
	taskPhi->DoCorrMixed(Form("<<3>>(4,-2,-2)%s",sGap.Data()),Form("<<4>>(2,2,-2,-2)%s",sGap.Data()),iNumSamplesRefs);
	taskPhi->SetPtBins(dPtBinningPhi);
	taskPhi->SetMergePosNeg();
	taskPhi->SetFlowMassRebin(2);

	FlowTask* taskPhi2 = new FlowTask(kPhi);
	taskPhi2->SetNumSamples(1);
	taskPhi2->SetEtaGap(dGap);
	taskPhi2->DoCorrMixed(Form("<<3>>(5,-3,-2)%s",sGap.Data()),Form("<<4>>(2,3,-2,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskPhi2->SetPtBins(dPtBinningPhi);
	taskPhi2->SetMergePosNeg();
	taskPhi2->SetFlowMassRebin(2);

	FlowTask* taskPhi3 = new FlowTask(kPhi);
	taskPhi3->SetNumSamples(1);
	taskPhi3->SetEtaGap(dGap);
	taskPhi3->DoCorrMixed(Form("<<3>>(6,-3,-3)%s",sGap.Data()),Form("<<4>>(3,3,-3,-3)%s",sGap.Data()),iNumSamplesRefs);
	taskPhi3->SetPtBins(dPtBinningPhi);
	taskPhi3->SetMergePosNeg();
	taskPhi3->SetFlowMassRebin(2);

	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/K0s");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/Lambda");
	TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/Phi");

	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/K0s/%s/",sTaskTag.Data());
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/Lambda/%s/",sTaskTag.Data());
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/Phi/%s/",sTaskTag.Data());



	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning);
	process->SetSaveMult(0);
	process->SetSaveInterSteps(0);
	process->SetFitCumulants(kFALSE);
	process->SetDebug(0);
	process->SetSaveInterSteps(1);
	process->SetTaskName(Form("UniFlow%s",sTaskTag.Data()));
	process->SetOutputFilePath(Form("%s/%s/",sOutputFilePath.Data(),sTaskTag.Data()));


	// process->AddTask(taskK0s);
	// process->AddTask(taskK0s2);
	// process->AddTask(taskK0s3);
	// process->AddTask(taskLambda);
	// process->AddTask(taskLambda2);
	// process->AddTask(taskLambda3);
	process->AddTask(taskPhi);
	// process->AddTask(taskPhi2);
	// process->AddTask(taskPhi3);

	process->Run();

	return;
}
// /
