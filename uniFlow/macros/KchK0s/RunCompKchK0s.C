/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunCompKchK0s()
{
	Double_t dEtaGap = 0.0;	TString sEtaGap = "gap00";
	// Double_t dEtaGap = 0.4;	TString sEtaGap = "gap04";
	// Double_t dEtaGap = 0.8;	TString sEtaGap = "gap08";

	Int_t iNumSamples = 10;

	// TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua";
	TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pp-16kl-nua";

	TString sOutputFilePath = sInputPath+"/output_compKchK0s_vn/"+sEtaGap+"/";



	// mult binning
	// Double_t dMultBinning[] = {0,100};
	// Double_t dMultBinning[] = {0,10,20,40,60,100};
	Double_t dMultBinning[] = {0,20,40,60,100};

	// pt binning
	Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.0,5.0,7.0};

	// ##### END Parameters setting ######

	// gROOT->AddIncludePath("-I~/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~vpacik/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");

	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(dEtaGap);
	taskRefs->SetMergePosNeg();


	FlowTask* taskKch = new FlowTask(FlowTask::kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(2);
	taskKch->SetEtaGap(dEtaGap);
	taskKch->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	// taskK0s->SetFittingOneGo(kTRUE);
	taskK0s->SetNumSamples(iNumSamples);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(dEtaGap);
	// taskK0s->SetInvMassRebin(2);
	// taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->SetSaveMult(kTRUE);
	process->SetFitCumulants(kFALSE);
	process->SetDebug();
	process->AddTask(taskRefs);
	process->AddTask(taskKch);
	process->AddTask(taskK0s);

	process->Run();

	return;
}
