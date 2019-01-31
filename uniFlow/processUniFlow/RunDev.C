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

void RunDev()
{
	Int_t iNumSamples = 1;
	Double_t dEtaGap = 0.0	;
	TString sEtaGap = "gap00";

	TString sInputPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow_sampled_20runs/";
	// TString sOutputFilePath = sInputPath+"/output_test/"+sEtaGap+"/";
	// TString sOutputFilePath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/mixed";
	TString sOutputFilePath = "./test/";

	// Double_t dMultBinning[] = {0,10,20,40,60,100};

	// Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningProton[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningLambda[] = {0.3,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningPhi[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.0};

	// Double_t dMultBinning[] = {0,5,10,20,30,40,50,60,70};
	// Double_t dPtBins[] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};
	std::vector<Double_t> dMultBinning = {10,20};
	// Double_t dPtBins[] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};
	Double_t dPtBins[] = {1.0,2.0,3.0,4.0};
	std::vector<Double_t> vecPtBins = {1.0,2.0,3.0,4.0};
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
		process->SetMultiplicityBins(dMultBinning);
		process->SetDebug(1);

		FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
		taskRefs->SetNumSamples(1);
		taskRefs->SetEtaGap(0.0);
		taskRefs->SetPtBins(vecPtBins);
		taskRefs->SetMergePosNeg(1);
		taskRefs->SetHarmonics(2);
		taskRefs->SetDoFourCorrelations(1);
		// taskRefs->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
		process->AddTask(taskRefs);

		FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
		taskCharged->SetNumSamples(1);
		taskCharged->SetPtBins(vecPtBins);
		taskCharged->SetMergePosNeg(1);
		taskCharged->SetHarmonics(2);
		taskCharged->SetEtaGap(0.0);
		taskCharged->SetDoFourCorrelations(1);
		// taskCharged->SetProcessMixedHarmonics("<<3>>(4,-2,-2)","<<4>>(2,2,-2,-2)");
		process->AddTask(taskCharged);

		FlowTask* taskCharged2 = new FlowTask(FlowTask::kCharged);
		taskCharged2->SetNumSamples(1);
		taskCharged2->SetPtBins(vecPtBins);
		taskCharged2->SetMergePosNeg(1);
		taskCharged2->SetHarmonics(2);
		taskCharged2->SetEtaGap(0.0);
		taskCharged2->SetDoFourCorrelations(1);
		// taskCharged->SetProcessMixedHarmonics("<<3>>(4,-2,-2)","<<4>>(2,2,-2,-2)");
		process->AddTask(taskCharged2);

		FlowTask* taskCharged3 = new FlowTask(FlowTask::kCharged);
		taskCharged3->SetNumSamples(1);
		taskCharged3->SetPtBins({1.2,1.4});
		taskCharged3->SetMergePosNeg(1);
		taskCharged3->SetHarmonics(2);
		taskCharged3->SetEtaGap(0.0);
		taskCharged3->SetDoFourCorrelations(1);
		// taskCharged->SetProcessMixedHarmonics("<<3>>(4,-2,-2)","<<4>>(2,2,-2,-2)");
		process->AddTask(taskCharged3);


		// FlowTask* taskPion = new FlowTask(FlowTask::kPion);
		// taskPion->SetNumSamples(1);
		// taskPion->SetEtaGap(0.0);
		// taskPion->SetPtBins(dPtBins,sizeof(dPtBins)/sizeof(dPtBins[0]));
		// taskPion->SetMergePosNeg(1);
		// taskPion->SetHarmonics(2);
		// // taskCharged->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
		// process->AddTask(taskPion);
		//
		FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
		taskK0s->SetNumSamples(1);
		taskK0s->SetEtaGap(0.0);
		taskK0s->SetPtBins(vecPtBins);
		taskK0s->SetMergePosNeg(1);
		taskK0s->SetHarmonics(2);
		taskK0s->SetDoFourCorrelations(1);
		// // taskCharged->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
		// process->AddTask(taskK0s);

		FlowTask* taskLambda = new FlowTask(FlowTask::kLambda);
		taskLambda->SetNumSamples(1);
		taskLambda->SetEtaGap(0.0);
		taskLambda->SetPtBins(vecPtBins);
		taskLambda->SetMergePosNeg(1);
		taskLambda->SetHarmonics(2);
		taskLambda->SetDoFourCorrelations(1);
		// // taskCharged->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
		// process->AddTask(taskLambda);

		FlowTask* taskPhi = new FlowTask(FlowTask::kPhi);
		taskPhi->SetNumSamples(1);
		taskPhi->SetEtaGap(0.0);
		taskPhi->SetPtBins(vecPtBins);
		taskPhi->SetMergePosNeg(1);
		taskPhi->SetHarmonics(2);
		taskPhi->SetDoFourCorrelations(1);
		// // taskCharged->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
		// process->AddTask(taskPhi);

		FlowTask* taskK0sMix = new FlowTask(FlowTask::kK0s);
		taskK0sMix->SetNumSamples(1);
		taskK0sMix->SetEtaGap(0.0);
		taskK0sMix->SetPtBins(vecPtBins);
		taskK0sMix->SetMergePosNeg(1);
		taskK0sMix->SetHarmonics(2);
		// taskK0sMix->SetDoFourCorrelations(1);
		taskK0sMix->SetProcessMixedHarmonics("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)");
		// // taskCharged->SetProcessMixedHarmonics("Cor3p4m2m2","Cor4p2p2m2m2");
		// process->AddTask(taskK0sMix);


		process->Run();


	return;
}
