/* RunF
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunK0sComp(const char* sOutputFilePath = "")
{
	// ##### Parameters setting ######

	// Double_t dMultBinning[] = {20,40};
	// Double_t dPtBinningK0s[] = {3.,5.};
	// Double_t dMultBinning[] = {0,20,40,60,100};
	Double_t dMultBinning[] = {0,10,20,40,60,100}; // Alex multiplicity
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.2,1.,2.,3.,4.,5.,7.,10.};
	// Double_t dPtBinningK0s[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // Run1
	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.25,1.5,1.75,2.,2.5,3.,3.5,4.,5.,6.,8.}; // Alex binning
	const char* sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay";

	// ##### END Parameters setting ######

	// gROOT->AddIncludePath("-I~/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~vpacik/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");
	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath);
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow_AN_PID");

	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process->SetOutputFilePath(Form("%s/output_compAlex",sInputPath));
	process->SetOutputFileName("Processed.root");

 	FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
	taskK0s->SetAlexFitting(kTRUE);
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

	FlowTask* taskKch = new FlowTask("Kch",FlowTask::kKaon);
	// taskKch->SetAlexFitting(kTRUE);
	taskKch->SetHarmonics(2);
	taskKch->SetEtaGap(0.8);
	// taskKch->SetInvMassRebin(2);
	// taskKch->SetFlowMassRebin(2);
	taskKch->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");
	process->AddTask(taskKch);


	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->Run();

	return;
}
