#ifdef __CLING__
#include "ProcessUniFlow.cpp"
#endif

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunProcess_mixedtest()
{
	Int_t iNumSamples = 1;
	Double_t dEtaGap = 0.0	;
	TString sEtaGap = "gap00";

	TString sInputPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/";
	// TString sOutputFilePath = sInputPath+"/output_test/"+sEtaGap+"/";
	TString sOutputFilePath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/mixed";

	// Double_t dMultBinning[] = {0,10,20,40,60,100};

	// Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningProton[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningLambda[] = {0.3,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningPhi[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.0};

	// Binning Run1 SP
	Double_t dMultBinning[] = {0,5,10,20,30,40,50,60,70};
	// Double_t dMultBinning[] = {0,20,40,60,61,62,63,64,65,70,80,100,150};
	Double_t dPtBins[] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};
	// ##### END Parameters setting ######

	#if defined (__CINT__)
		gROOT->LoadMacro("~vpacik/Codes/ALICE/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");
	#endif

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(1);
	taskCharged->SetEtaGap(-1.0);
	taskCharged->SetPtBins(dPtBins,sizeof(dPtBins)/sizeof(dPtBins[0]));
	taskCharged->SetMergePosNeg(1);
	taskCharged->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
	process->AddTask(taskCharged);

	FlowTask* taskCharged2 = new FlowTask(FlowTask::kCharged);
	taskCharged2->SetNumSamples(1);
	taskCharged2->SetEtaGap(1.0);
	taskCharged2->SetPtBins(dPtBins,sizeof(dPtBins)/sizeof(dPtBins[0]));
	taskCharged2->SetMergePosNeg(1);
	taskCharged2->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
	process->AddTask(taskCharged2);

	FlowTask* taskCharged3 = new FlowTask(FlowTask::kCharged);
	taskCharged3->SetNumSamples(1);
	taskCharged3->SetEtaGap(1.0);
	taskCharged3->SetPtBins(dPtBins,sizeof(dPtBins)/sizeof(dPtBins[0]));
	taskCharged3->SetMergePosNeg(1);
	taskCharged3->SetProcessMixedHarmonics("Cor3p6m3m3","Cor4p3p3m3m3");
	process->AddTask(taskCharged3);

	FlowTask* taskCharged4 = new FlowTask(FlowTask::kCharged);
	taskCharged4->SetNumSamples(1);
	taskCharged4->SetEtaGap(1.0);
	taskCharged4->SetPtBins(dPtBins,sizeof(dPtBins)/sizeof(dPtBins[0]));
	taskCharged4->SetMergePosNeg(1);
	taskCharged4->SetProcessMixedHarmonics("Cor3p5m2m3","Cor4p2p3m2m3");
	process->AddTask(taskCharged4);

	process->Run();


	return;
}
