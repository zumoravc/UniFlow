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
	// TString sInputPath = "/home/alidock/ana/UniFlow/uniFlow/task";
	TString sInputPath = "/home/alidock/ana/output/LHC15o/train_7346";
	// TString sOutputFilePath = sInputPath+"/output_test/"+sEtaGap+"/";
	// TString sOutputFilePath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/mixed";
	// TString sOutputFilePath = "./test/";
	// TString sOutputFilePath = "/mnt/CodesALICE/Flow/uniFlow/processUniFlow/test-K0s";
	// TString sOutputFilePath = "/mnt/CodesALICE/Flow/uniFlow/processUniFlow/test-Lambda";
	TString sOutputFilePath = sInputPath + "/processUniFlow";

	// Double_t dMultBinning[] = {0,10,20,40,60,100};

	// Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningProton[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningLambda[] = {0.3,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningPhi[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.0};

	// Double_t dMultBinning[] = {0,5,10,20,30,40,50,60,70};
	// Double_t dPtBins[] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};
	std::vector<Double_t> dMultBinning = {0,5,10,20,30,40,50,60,70,80};
	// Double_t dPtBins[] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};
	Double_t dPtBins[] = {2.0,2.5};
	std::vector<Double_t> vecPtBins = {0.0,0.5,1.0,1.5,2.0,2.5,3.0};
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

	FlowTask* taskRefs = new FlowTask(kRefs, "sixFMC");
	taskRefs->SetNumSamples(10);
	//taskRefs->SetEtaGap(1.0);
	//taskRefs->SetPtBins(vecPtBins);
	taskRefs->SetMergePosNeg(1);
	taskRefs->DoSixFMC({2,3,3,-2,-3,-3});
	//taskRefs->SetHarmonics(2);
	// taskRefs->DoCumOrderMax(kTwo);
	//taskRefs->SetRebinning(kFALSE);
	process->AddTask(taskRefs);

	FlowTask* taskRefs2 = new FlowTask(kRefs, "eightFMC");
	taskRefs2->SetNumSamples(10);
	taskRefs2->SetMergePosNeg(1);
	taskRefs2->DoSixFMC({2,2,3,-2,-2,-3});
	taskRefs2->DoFourFMC(kTRUE);
	process->AddTask(taskRefs2);

	FlowTask* taskRefs3 = new FlowTask(kRefs, "eightFMC1");
	taskRefs3->SetNumSamples(10);
	taskRefs3->SetMergePosNeg(1);
	taskRefs3->DoEightFMC({2,2,3,3,-2,-2,-3,-3});
	process->AddTask(taskRefs3);

	FlowTask* taskRefs4 = new FlowTask(kRefs, "eightFMC2");
	taskRefs4->SetNumSamples(10);
	taskRefs4->SetMergePosNeg(1);
	taskRefs4->DoEightFMC({2,2,2,3,-2,-2,-2,-3});
	process->AddTask(taskRefs4);

	FlowTask* taskRefs5 = new FlowTask(kRefs, "eightFMC3");
	taskRefs5->SetNumSamples(10);
	taskRefs5->SetMergePosNeg(1);
	taskRefs5->DoEightFMC({2,3,3,3,-2,-3,-3,-3});
	process->AddTask(taskRefs5);

	// FlowTask* taskRefs = new FlowTask(kRefs);
	// taskRefs->SetNumSamples(10);
	// taskRefs->SetEtaGap(1.0);
	// taskRefs->SetPtBins(vecPtBins);
	// taskRefs->SetMergePosNeg(1);
	// taskRefs->SetHarmonics(2);
	// taskRefs->DoCumOrderMax(kTwo);
	// //taskRefs->SetRebinning(kFALSE);
	// process->AddTask(taskRefs);

	// FlowTask* taskRefs3 = new FlowTask(kRefs);
	// taskRefs3->SetNumSamples(10);
	// // taskRefs2->SetEtaGap();
	// taskRefs3->SetPtBins(vecPtBins);
	// taskRefs3->SetMergePosNeg(1);
	// taskRefs3->SetHarmonics(2);
	// taskRefs3->DoCumOrderMax(kTwo);
	// process->AddTask(taskRefs3);
	//
	// FlowTask* taskRefs3 = new FlowTask(kRefs);
	// taskRefs3->SetNumSamples(10);
	// taskRefs3->SetEtaGap(1.0);
	// taskRefs3->SetPtBins(vecPtBins);
	// taskRefs3->SetMergePosNeg(1);
	// taskRefs3->SetHarmonics(3);
	// taskRefs3->DoCumOrderMax(kTwo);
	// //taskRefs->SetRebinning(kFALSE);
	// process->AddTask(taskRefs3);
	//
	// FlowTask* taskRefs4 = new FlowTask(kRefs);
	// taskRefs4->SetNumSamples(10);
	// taskRefs4->SetEtaGap(1.0);
	// taskRefs4->SetPtBins(vecPtBins);
	// taskRefs4->SetMergePosNeg(1);
	// taskRefs4->SetHarmonics(4);
	// taskRefs4->DoCumOrderMax(kTwo);
	// //taskRefs->SetRebinning(kFALSE);
	// process->AddTask(taskRefs4);
	//
	// FlowTask* taskRefs5 = new FlowTask(kRefs);
	// taskRefs5->SetNumSamples(10);
	// //taskRefs5->SetEtaGap(1.0);
	// taskRefs5->SetPtBins(vecPtBins);
	// taskRefs5->SetMergePosNeg(1);
	// taskRefs5->SetHarmonics(2);
	// taskRefs5->DoCumOrderMax(kFour);
	// process->AddTask(taskRefs5);

	// FlowTask* taskCharged = new FlowTask(kCharged);
	// taskCharged->SetNumSamples(1);
	// taskCharged->SetPtBins(vecPtBins);
	// taskCharged->SetMergePosNeg(1);
	// taskCharged->SetHarmonics(2);
	// taskCharged->SetEtaGap(1.0);
	// taskCharged->DoCumOrderMax(kTwo);
	// // taskCharged->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)");
	// process->AddTask(taskCharged);

	// FlowTask* taskCharged2 = new FlowTask(kCharged);
	// taskCharged2->SetNumSamples(1);
	// taskCharged2->SetPtBins(vecPtBins);
	// taskCharged2->SetMergePosNeg(1);
	// taskCharged2->SetHarmonics(2);
	// taskCharged2->SetEtaGap(0.0);
	// taskCharged2->DoCorrMixed("<<3>>(5,-3,-2)_2sub(0)","<<4>>(2,3,-2,-3)_2sub(0)");
	// process->AddTask(taskCharged2);

	// FlowTask* taskCharged3 = new FlowTask(kCharged);
	// taskCharged3->SetNumSamples(1);
	// taskCharged3->SetPtBins(vecPtBins);
	// taskCharged3->SetMergePosNeg(1);
	// taskCharged3->SetHarmonics(2);
	// taskCharged3->SetEtaGap(0.0);
	// taskCharged3->DoCorrMixed("<<3>>(6,-3,-3)_2sub(0)","<<4>>(3,3,-3,-3)_2sub(0)");
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
	// FlowTask* taskK0s = new FlowTask(kK0s);
	// taskK0s->SetNumSamples(1);
	// taskK0s->SetEtaGap(0.0);
	// taskK0s->SetPtBins(vecPtBins);
	// taskK0s->SetMergePosNeg(1);
	// taskK0s->SetHarmonics(2);
	// taskK0s->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)",5);
	// process->AddTask(taskK0s);

	// FlowTask* taskLambda = new FlowTask(kLambda);
	// taskLambda->SetNumSamples(1);
	// taskLambda->SetEtaGap(0.0);
	// taskLambda->SetPtBins(vecPtBins);
	// taskLambda->SetMergePosNeg(1);
	// taskLambda->SetHarmonics(2);
	// taskLambda->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)",5);
	// process->AddTask(taskLambda);

	// FlowTask* taskPhi = new FlowTask(kPhi);
	// taskPhi->SetNumSamples(1);
	// taskPhi->SetEtaGap(0.0);
	// taskPhi->SetPtBins(vecPtBins);
	// taskPhi->SetMergePosNeg(1);
	// taskPhi->SetHarmonics(2);
	// taskPhi->DoCorrMixed("<<3>>(4,-2,-2)_2sub(0)","<<4>>(2,2,-2,-2)_2sub(0)",5);
	// taskPhi->SetFitPhiSubtLS(1,1,0.99,1.01);
	// process->AddTask(taskPhi);



	process->Run();


	return;
}
