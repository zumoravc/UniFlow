/* RunFitting
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunFitAlex(const char* sOutputFilePath = "")
{
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	// Double_t dMultBinning[] = {20,40};
	// Double_t dPtBinningK0s[] = {3.,5.};
	Double_t dMultBinning[] = {0,20,40,60,100};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningK0s[] = {0.2,1.,2.,3.,4.,5.,7.,10.,20.};


	ProcessUniFlow* process = new ProcessUniFlow();

	process->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay");
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow_AN_PID");

	process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay/fitMine");
	process->SetOutputFileName("Processed.root");

 	FlowTask* taskK0s = new FlowTask("K0s_mine",FlowTask::kK0s);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(0.8);
	taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskK0s);

	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->Run();

	ProcessUniFlow* processAlex = new ProcessUniFlow();
	processAlex->SetInputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay/");
	processAlex->SetInputFileName("AnalysisResults.root");
	processAlex->SetTaskName("UniFlow_AN_PID");

	processAlex->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay/fitAlex");
	processAlex->SetOutputFileName("Processed.root");

	FlowTask* taskK0sAlex = new FlowTask("K0s_Alex",FlowTask::kK0s);
	taskK0sAlex->SetAlexFitting(kTRUE);
	taskK0sAlex->SetHarmonics(2);
	taskK0sAlex->SetEtaGap(0.8);
	taskK0sAlex->SetInvMassRebin(2);
	taskK0sAlex->SetFlowMassRebin(2);
	taskK0sAlex->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0sAlex->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	processAlex->AddTask(taskK0sAlex);

	processAlex->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	processAlex->Run();

	return;
}
