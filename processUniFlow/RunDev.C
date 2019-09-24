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
	TString sInputPath = "/home/alidock/ana/output/LHC15o/grid_fullSt/";
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
	process->SetInputFileName("Minus.root");
	process->SetTaskName("UniFlow");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Minus.root");
	process->SetMultiplicityBins(dMultBinning);
	process->SetDebug(1);
	process->SetSaveInterSteps(1);

	FlowTask* taskRefs = new FlowTask(kRefs, "");
	taskRefs->SetNumSamples(10);
	//taskRefs->SetEtaGap(1.0);
	//taskRefs->SetPtBins(vecPtBins);
	taskRefs->SetMergePosNeg(1);
	taskRefs->DoSixFMC({2,3,3,-2,-3,-3});
	// taskRefs->IsHijing(kTRUE);
	//taskRefs->SetHarmonics(2);
	// taskRefs->DoCumOrderMax(kTwo);
	//taskRefs->SetRebinning(kFALSE);
	process->AddTask(taskRefs);

	FlowTask* taskRefs2F = new FlowTask(kRefs, "");
	taskRefs2F->SetNumSamples(10);
	taskRefs2F->SetMergePosNeg(1);
	taskRefs2F->DoFourFMC({2,3,-2,-3});
	// taskRefs2->DoFourFMCwithSix(kTRUE);
	// taskRefs2->IsHijing(kTRUE);
	process->AddTask(taskRefs2F);

	FlowTask* taskRefs2 = new FlowTask(kRefs, "eightFMC");
	taskRefs2->SetNumSamples(10);
	taskRefs2->SetMergePosNeg(1);
	taskRefs2->DoSixFMC({2,2,3,-2,-2,-3});
	// taskRefs2->DoFourFMCwithSix(kTRUE);
	// taskRefs2->IsHijing(kTRUE);
	process->AddTask(taskRefs2);

	FlowTask* taskRefs3 = new FlowTask(kRefs, "eightFMC1");
	taskRefs3->SetNumSamples(10);
	taskRefs3->SetMergePosNeg(1);
	taskRefs3->DoEightFMC({2,2,3,3,-2,-2,-3,-3});
	// taskRefs3->IsHijing(kTRUE);
	process->AddTask(taskRefs3);

	FlowTask* taskRefs4 = new FlowTask(kRefs, "eightFMC2");
	taskRefs4->SetNumSamples(10);
	taskRefs4->SetMergePosNeg(1);
	taskRefs4->DoEightFMC({2,2,2,3,-2,-2,-2,-3});
	// taskRefs4->IsHijing(kTRUE);
	process->AddTask(taskRefs4);

	FlowTask* taskRefs5 = new FlowTask(kRefs, "eightFMC3");
	taskRefs5->SetNumSamples(10);
	taskRefs5->SetMergePosNeg(1);
	taskRefs5->DoEightFMC({2,3,3,3,-2,-3,-3,-3});
	// taskRefs5->IsHijing(kTRUE);
	process->AddTask(taskRefs5);

	FlowTask* taskRefsV2 = new FlowTask(kRefs, "");
	taskRefsV2->SetNumSamples(10);
	taskRefsV2->SetMergePosNeg(1);
	taskRefsV2->SetHarmonics(2);
	taskRefsV2->DoCumOrderMax(kTwo);
	process->AddTask(taskRefsV2);

	FlowTask* taskRefsV3 = new FlowTask(kRefs, "");
	taskRefsV3->SetNumSamples(10);
	taskRefsV3->SetMergePosNeg(1);
	taskRefsV3->SetHarmonics(3);
	taskRefsV3->DoCumOrderMax(kTwo);
	process->AddTask(taskRefsV3);

	FlowTask* taskRefsV2g = new FlowTask(kRefs, "");
	taskRefsV2g->SetNumSamples(10);
	taskRefsV2g->SetEtaGap(1.0);
	taskRefsV2g->SetMergePosNeg(1);
	taskRefsV2g->SetHarmonics(2);
	taskRefsV2g->DoCumOrderMax(kTwo);
	process->AddTask(taskRefsV2g);

	FlowTask* taskRefsV3g = new FlowTask(kRefs, "");
	taskRefsV3g->SetNumSamples(10);
	taskRefsV3g->SetEtaGap(1.0);
	taskRefsV3g->SetMergePosNeg(1);
	taskRefsV3g->SetHarmonics(3);
	taskRefsV3g->DoCumOrderMax(kTwo);
	process->AddTask(taskRefsV3g);

	FlowTask* taskRefsV4 = new FlowTask(kRefs, "");
	taskRefsV4->SetNumSamples(10);
	taskRefsV4->SetMergePosNeg(1);
	taskRefsV4->SetHarmonics(4);
	taskRefsV4->DoCumOrderMax(kTwo);
	process->AddTask(taskRefsV4);

	FlowTask* taskRefsV4g = new FlowTask(kRefs, "");
	taskRefsV4g->SetNumSamples(10);
	taskRefsV4g->SetEtaGap(1.0);
	taskRefsV4g->SetMergePosNeg(1);
	taskRefsV4g->SetHarmonics(4);
	taskRefsV4g->DoCumOrderMax(kTwo);
	process->AddTask(taskRefsV4g);



	FlowTask* taskCharged = new FlowTask(kCharged);
	taskCharged->SetNumSamples(10);
	taskCharged->SetPtBins(vecPtBins);
	taskCharged->SetMergePosNeg(1);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(1.0);
	taskCharged->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged);

	FlowTask* taskCharged2 = new FlowTask(kCharged);
	taskCharged2->SetNumSamples(10);
	taskCharged2->SetPtBins(vecPtBins);
	taskCharged2->SetMergePosNeg(1);
	taskCharged2->SetHarmonics(2);
	// taskCharged2->SetEtaGap(0.0);
	taskCharged2->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged2);

	FlowTask* taskCharged3 = new FlowTask(kCharged);
	taskCharged3->SetNumSamples(10);
	taskCharged3->SetPtBins(vecPtBins);
	taskCharged3->SetMergePosNeg(1);
	taskCharged3->SetHarmonics(3);
	// taskCharged3->SetEtaGap(0.0);
	taskCharged3->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged3);

	FlowTask* taskCharged4 = new FlowTask(kCharged);
	taskCharged4->SetNumSamples(10);
	taskCharged4->SetPtBins(vecPtBins);
	taskCharged4->SetMergePosNeg(1);
	taskCharged4->SetHarmonics(3);
	taskCharged4->SetEtaGap(1.0);
	taskCharged4->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged4);

	FlowTask* taskCharged5 = new FlowTask(kCharged);
	taskCharged5->SetNumSamples(10);
	taskCharged5->SetPtBins(vecPtBins);
	taskCharged5->SetMergePosNeg(1);
	taskCharged5->SetHarmonics(4);
	// taskCharged5->SetEtaGap(0.0);
	taskCharged5->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged5);

	FlowTask* taskCharged6 = new FlowTask(kCharged);
	taskCharged6->SetNumSamples(10);
	taskCharged6->SetPtBins(vecPtBins);
	taskCharged6->SetMergePosNeg(1);
	taskCharged6->SetHarmonics(4);
	taskCharged6->SetEtaGap(1.0);
	taskCharged6->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged6);

	FlowTask* taskRefsF2 = new FlowTask(kRefs, "");
	taskRefsF2->SetNumSamples(10);
	taskRefsF2->SetMergePosNeg(1);
	taskRefsF2->SetHarmonics(2);
	taskRefsF2->DoCumOrderMax(kFour);
	process->AddTask(taskRefsF2);

	FlowTask* taskRefsF3 = new FlowTask(kRefs, "");
	taskRefsF3->SetNumSamples(10);
	taskRefsF3->SetMergePosNeg(1);
	taskRefsF3->SetHarmonics(3);
	taskRefsF3->DoCumOrderMax(kFour);
	process->AddTask(taskRefsF3);

	FlowTask* taskCharged7 = new FlowTask(kCharged);
	taskCharged7->SetNumSamples(10);
	taskCharged7->SetPtBins(vecPtBins);
	taskCharged7->SetMergePosNeg(1);
	taskCharged7->SetHarmonics(2);
	taskCharged7->DoCumOrderMax(kFour);
	process->AddTask(taskCharged7);

	FlowTask* taskCharged8 = new FlowTask(kCharged);
	taskCharged8->SetNumSamples(10);
	taskCharged8->SetPtBins(vecPtBins);
	taskCharged8->SetMergePosNeg(1);
	taskCharged8->SetHarmonics(3);
	taskCharged8->DoCumOrderMax(kFour);
	process->AddTask(taskCharged8);



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
