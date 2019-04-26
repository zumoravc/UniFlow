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

	std::vector<Double_t> dMultBinning = {0,5,10,20,30,40,50,60};
	std::vector<Double_t> dPtBinningK0s = {0.8,1.2,1.6,2.0,2.4,3.0,3.6,4.6,6.0};
	std::vector<Double_t> dPtBinningLambda = {0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.4,6.0};
	// std::vector<Double_t> dPtBinningPhi = {1.0,2.0,3.0,4.0,5.0,6.0};
	std::vector<Double_t> dPtBinningPhi = {1.0,1.5,2.0,3.0,4.5,6.0};
	// std::vector<Double_t> dPtBinningPhi = {1.0,2.0,3.0,4.0,6.0};

	// std::vector<Double_t> dPtBinningK0sV633 = {0.8,1.4,2.0,2.6,3.2,4.4,6.0};
	// std::vector<Double_t> dPtBinningK0sV633 = {0.8,1.2,1.6,2.0,2.4,3.0,3.6,4.6,6.0};
	// std::vector<Double_t> dPtBinningLambdaV633 = {1.0,2.0,3.0,4.0,6.0};

	// std::vector<Double_t> dMultBinning_2 = {0,5,50,60};
	// std::vector<Double_t> dPtBinningK0s_2 = {0.8,1.2,1.6,2.0,2.4,3.0,4.0,5.0,6.0};
	// std::vector<Double_t> dPtBinningLambda_2 = {1.0,2.0,3.0,4.0,6.0};
	// std::vector<Double_t> dPtBinningPhi_2 = {0.5,1.0,1.5,2.0,3.0,4.0,6.0};

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


	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/K0s/");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/K0s_v422");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/K0s_v532");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/Lambda");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/Lambda_v532");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/Lambda_v633");
	TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/Phi");
	// TString sOutputFilePath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/All");


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
	process->SetSaveInterSteps(1);

	// process->AddTask(taskK0s);
	// process->AddTask(taskK0s2);
	// process->AddTask(taskK0s3);

	// process->AddTask(taskLambda);
	// process->AddTask(taskLambda2);
	// process->AddTask(taskLambda3);

	process->AddTask(taskPhi);
	process->AddTask(taskPhi2);
	process->AddTask(taskPhi3);

	process->Run();


	// process->SetMultiplicityBins(dMultBinning_2);
	// process->SetOutputFilePath(Form("%s_lowMult",sOutputFilePath.Data()));
	//
	// taskK0s->SetPtBins(dPtBinningK0s_2);
	// taskK0s->SetFlowMassRebin(4);
	// taskK0s2->SetPtBins(dPtBinningK0s_2);
	// taskK0s2->SetFlowMassRebin(4);
	// taskK0s3->SetPtBins(dPtBinningK0s_2);
	// taskK0s3->SetFlowMassRebin(4);
	//
	// taskLambda->SetPtBins(dPtBinningLambda_2);
	// taskLambda2->SetPtBins(dPtBinningLambda_2);
	// taskLambda3->SetPtBins(dPtBinningLambda_2);
	//
	// taskPhi->SetPtBins(dPtBinningPhi_2);
	// taskPhi2->SetPtBins(dPtBinningPhi_2);
	// taskPhi3->SetPtBins(dPtBinningPhi_2);

	// process->Run();

	return;
}
// /
