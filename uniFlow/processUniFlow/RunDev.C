// NEED to run .X Load.C beforehand

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

	// TString sInputPath = "~/Codes/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow_sampled_20runs/";
	TString sInputPath = "~/Codes/Flow/uniFlow/results/trains/CF_PbPb/6527_20190218-2140/merge/";
	// TString sOutputFilePath = sInputPath+"/output_test/"+sEtaGap+"/";
	// TString sOutputFilePath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/mixed";
	// TString sOutputFilePath = "./test/";
	TString sOutputFilePath = "/mnt/CodesALICE/Flow/uniFlow/processUniFlow/test/";

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
	Double_t dPtBins[] = {2.0,2.5};
	std::vector<Double_t> vecPtBins = {1.5,2.0};
	// ##### END Parameters setting ######

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning);
	process->SetDebug(1);
	process->SetSaveInterSteps(1);

	FlowTask* taskRefs = new FlowTask(kRefs);
	taskRefs->SetNumSamples(1);
	taskRefs->SetEtaGap(0.0);
	taskRefs->SetPtBins(vecPtBins);
	taskRefs->SetMergePosNeg(1);
	taskRefs->SetHarmonics(2);
	taskRefs->DoCumOrderMax(4);
	// process->AddTask(taskRefs);

	FlowTask* taskCharged = new FlowTask(kCharged);
	taskCharged->SetNumSamples(1);
	taskCharged->SetPtBins(vecPtBins);
	taskCharged->SetMergePosNeg(1);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(0.0);
	taskCharged->DoCumOrderMax(kFour);
	taskCharged->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)");
	// process->AddTask(taskCharged);

	FlowTask* taskCharged2 = new FlowTask(kCharged);
	taskCharged2->SetNumSamples(1);
	taskCharged2->SetPtBins(vecPtBins);
	taskCharged2->SetMergePosNeg(1);
	taskCharged2->SetHarmonics(2);
	taskCharged2->SetEtaGap(0.0);
	taskCharged2->DoCorrMixed("<<3>>(5,-3,-2)_2sub(0)","<<4>>(2,3,-2,-3)_2sub(0)");
	// process->AddTask(taskCharged2);

	FlowTask* taskCharged3 = new FlowTask(kCharged);
	taskCharged3->SetNumSamples(1);
	taskCharged3->SetPtBins(vecPtBins);
	taskCharged3->SetMergePosNeg(1);
	taskCharged3->SetHarmonics(2);
	taskCharged3->SetEtaGap(0.0);
	taskCharged3->DoCorrMixed("<<3>>(6,-3,-3)_2sub(0)","<<4>>(3,3,-3,-3)_2sub(0)");
	// process->AddTask(taskCharged3);


	// FlowTask* taskPion = new FlowTask(kPion);
	// taskPion->SetNumSamples(1);
	// taskPion->SetEtaGap(0.0);
	// taskPion->SetPtBins(dPtBins,sizeof(dPtBins)/sizeof(dPtBins[0]));
	// taskPion->SetMergePosNeg(1);
	// taskPion->SetHarmonics(2);
	// // taskCharged->DoCorrMixed("Cor3p4m2m2","Cor4p2p2m2m2");
	// process->AddTask(taskPion);
	//
	FlowTask* taskK0s = new FlowTask(kK0s);
	taskK0s->SetNumSamples(1);
	taskK0s->SetEtaGap(0.0);
	taskK0s->SetPtBins(vecPtBins);
	taskK0s->SetMergePosNeg(1);
	taskK0s->SetHarmonics(2);
	taskK0s->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)",5);
	// process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask(kLambda);
	taskLambda->SetNumSamples(1);
	taskLambda->SetEtaGap(0.0);
	taskLambda->SetPtBins(vecPtBins);
	taskLambda->SetMergePosNeg(1);
	taskLambda->SetHarmonics(2);
	taskLambda->DoCumOrderMax(4);
	// // taskCharged->DoCorrMixed("Cor3p4m2m2","Cor4p2p2m2m2");
	// process->AddTask(taskLambda);

	FlowTask* taskPhi = new FlowTask(kPhi);
	taskPhi->SetNumSamples(1);
	taskPhi->SetEtaGap(0.0);
	taskPhi->SetPtBins(vecPtBins);
	taskPhi->SetMergePosNeg(1);
	taskPhi->SetHarmonics(2);
	taskPhi->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)",5);
	process->AddTask(taskPhi);



	process->Run();


	return;
}
