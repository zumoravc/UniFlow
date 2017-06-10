/* RunFitting
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunFitting(const char* sOutputFilePath = "")
{
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	ProcessUniFlow* process = new ProcessUniFlow();

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full/CENT_woSDD_16q/");
	// process->SetInputFileName("AnalysisResults_CENTwSDD_16q.root");
	// process->SetTaskName("UniFlow_V0A");
	// process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test");
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	// process->SetOutputFilePath(sOutputFilePath);
	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/fitTest");
	process->SetOutputFileName("UniFlowFitting.root");
	// process->SetOutputFileMode("UPDATE");
	// process->SetDebug();
	// process->SuggestMultBinning(4);

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	Double_t dMultBinning[] = {0,20};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	FlowTask* taskRefs = new FlowTask("Refs",FlowTask::kRefs);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(0.8);
	// process->AddTask(taskRefs);

	Double_t dPtBinningK0s[] = {1.,1.5};
	Double_t dPtBinningLambda[] = {1.,1.5};
	Double_t dPtBinningPhi[] = {1.,3.};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningLambda[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,10.,20.};


	// Double_t dPtBinning[] = {0.5,1.};
	// Double_t dPtBinning[] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.};
	// Int_t iSizePt = sizeof(dPtBinning)/sizeof(dPtBinning[0]);

 	FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(0.8);
	// taskK0s->SetInvMassRebin(2);
	// taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	taskK0s->SetFittingRange(0.45,0.55);
	taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(0.8);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Pos");
	taskLambda->SetMergePosNeg();
	taskLambda->SetFittingRange(1.1,1.14);
	taskLambda->SetFittingRejectNumSigmas(5);
	process->AddTask(taskLambda);

	FlowTask* taskPhi = new FlowTask("Phi",FlowTask::kPhi);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(0.8);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	// taskPhi->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Pos");
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRejectNumSigmas(1);
	taskPhi->SetFittingRange(1.,1.05);
	process->AddTask(taskPhi);

	process->Run();


	return;
}
