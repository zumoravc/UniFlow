// NEED to run .X Load.C beforehand

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

//small systems
void RunDevSS()
{
	Int_t iNumSamples = 1;
	Double_t dEtaGap = 0.0	;
	TString sEtaGap = "gap00";

	TString sInputPath = "/home/alidock/ana/output/pPb_LHC16q/HADD/";
	TString sOutputFilePath = sInputPath + "processUniFlow";
	std::vector<Double_t> dMultBinning = {0,1};
	// Double_t dPtBins[] = {2.0,2.5};
	// std::vector<Double_t> vecPtBins = {0.0,0.5,1.0,1.5,2.0,2.5,3.0};
	// std::vector<Double_t> vecPtBins = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
	std::vector<Double_t> vecPtBins = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0};
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

	FlowTask* taskRefs = new FlowTask(kRefs, "");
	taskRefs->SetNumSamples(1);
	taskRefs->SetMergePosNeg(1);
	taskRefs->SetHarmonics(2);
	taskRefs->DoCumOrderMax(kTwo);
	process->AddTask(taskRefs);

	FlowTask* taskRefs2 = new FlowTask(kRefs, "");
	taskRefs2->SetNumSamples(1);
	taskRefs2->SetMergePosNeg(1);
	taskRefs2->SetHarmonics(2);
	taskRefs2->DoCumOrderMax(kTwo);
	taskRefs2->SetEtaGap(0.0);
	process->AddTask(taskRefs2);

	FlowTask* taskRefs3 = new FlowTask(kRefs, "");
	taskRefs3->SetNumSamples(1);
	taskRefs3->SetMergePosNeg(1);
	taskRefs3->SetHarmonics(2);
	taskRefs3->DoCumOrderMax(kTwo);
	taskRefs3->SetEtaGap(0.8);
	process->AddTask(taskRefs3);

	FlowTask* taskRefs4 = new FlowTask(kRefs, "");
	taskRefs4->SetNumSamples(1);
	taskRefs4->SetMergePosNeg(1);
	taskRefs4->SetHarmonics(2);
	taskRefs4->DoCumOrderMax(kTwo);
	taskRefs4->SetEtaGap(1.0);
	process->AddTask(taskRefs4);

	//charged

	FlowTask* taskCharged = new FlowTask(kCharged);
	taskCharged->SetNumSamples(1);
	taskCharged->SetPtBins(vecPtBins);
	taskCharged->SetMergePosNeg(1);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(1.0);
	taskCharged->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged);

	FlowTask* taskCharged2 = new FlowTask(kCharged);
	taskCharged2->SetNumSamples(1);
	taskCharged2->SetPtBins(vecPtBins);
	taskCharged2->SetMergePosNeg(1);
	taskCharged2->SetHarmonics(2);
	taskCharged2->SetEtaGap(0.0);
	taskCharged2->DoCumOrderMax(kTwo);
	// taskCharged2->DoCorrMixed("<<3>>(5,-3,-2)_2sub(0)","<<4>>(2,3,-2,-3)_2sub(0)");
	process->AddTask(taskCharged2);

	FlowTask* taskCharged3 = new FlowTask(kCharged);
	taskCharged3->SetNumSamples(1);
	taskCharged3->SetPtBins(vecPtBins);
	taskCharged3->SetMergePosNeg(1);
	taskCharged3->SetHarmonics(2);
	taskCharged3->SetEtaGap(0.8);
	taskCharged3->DoCumOrderMax(kTwo);
	// taskCharged3->DoCorrMixed("<<3>>(6,-3,-3)_2sub(0)","<<4>>(3,3,-3,-3)_2sub(0)");
	process->AddTask(taskCharged3);
	//
	FlowTask* taskCharged4 = new FlowTask(kCharged);
	taskCharged4->SetNumSamples(1);
	taskCharged4->SetPtBins(vecPtBins);
	taskCharged4->SetMergePosNeg(1);
	taskCharged4->SetHarmonics(2);
	taskCharged4->DoCumOrderMax(kTwo);
	process->AddTask(taskCharged4);

	//pions

	FlowTask* taskPion = new FlowTask(kPion);
	taskPion->SetNumSamples(1);
	taskPion->SetPtBins(vecPtBins);
	taskPion->SetMergePosNeg(1);
	taskPion->SetHarmonics(2);
	taskPion->DoCumOrderMax(kTwo);
	process->AddTask(taskPion);

	FlowTask* taskPion2 = new FlowTask(kPion);
	taskPion2->SetNumSamples(1);
	taskPion2->SetEtaGap(0.0);
	taskPion2->SetPtBins(vecPtBins);
	taskPion2->SetMergePosNeg(1);
	taskPion2->SetHarmonics(2);
	taskPion2->DoCumOrderMax(kTwo);
	process->AddTask(taskPion2);

	FlowTask* taskPion3 = new FlowTask(kPion);
	taskPion3->SetNumSamples(1);
	taskPion3->SetEtaGap(0.8);
	taskPion3->SetPtBins(vecPtBins);
	taskPion3->SetMergePosNeg(1);
	taskPion3->SetHarmonics(2);
	taskPion3->DoCumOrderMax(kTwo);
	process->AddTask(taskPion3);

	FlowTask* taskPion4 = new FlowTask(kPion);
	taskPion4->SetNumSamples(1);
	taskPion4->SetEtaGap(1.0);
	taskPion4->SetPtBins(vecPtBins);
	taskPion4->SetMergePosNeg(1);
	taskPion4->SetHarmonics(2);
	taskPion4->DoCumOrderMax(kTwo);
	process->AddTask(taskPion4);

	//kaons

	FlowTask* taskKaon = new FlowTask(kKaon);
	taskKaon->SetNumSamples(1);
	taskKaon->SetPtBins(vecPtBins);
	taskKaon->SetMergePosNeg(1);
	taskKaon->SetHarmonics(2);
	taskKaon->DoCumOrderMax(kTwo);
	process->AddTask(taskKaon);

	FlowTask* taskKaon2 = new FlowTask(kKaon);
	taskKaon2->SetNumSamples(1);
	taskKaon2->SetEtaGap(0.0);
	taskKaon2->SetPtBins(vecPtBins);
	taskKaon2->SetMergePosNeg(1);
	taskKaon2->SetHarmonics(2);
	taskKaon2->DoCumOrderMax(kTwo);
	process->AddTask(taskKaon2);

	FlowTask* taskKaon3 = new FlowTask(kKaon);
	taskKaon3->SetNumSamples(1);
	taskKaon3->SetEtaGap(0.8);
	taskKaon3->SetPtBins(vecPtBins);
	taskKaon3->SetMergePosNeg(1);
	taskKaon3->SetHarmonics(2);
	taskKaon3->DoCumOrderMax(kTwo);
	process->AddTask(taskKaon3);

	FlowTask* taskKaon4 = new FlowTask(kKaon);
	taskKaon4->SetNumSamples(1);
	taskKaon4->SetEtaGap(1.0);
	taskKaon4->SetPtBins(vecPtBins);
	taskKaon4->SetMergePosNeg(1);
	taskKaon4->SetHarmonics(2);
	taskKaon4->DoCumOrderMax(kTwo);
	process->AddTask(taskKaon4);

	//protons

	FlowTask* taskProton = new FlowTask(kProton);
	taskProton->SetNumSamples(1);
	taskProton->SetPtBins(vecPtBins);
	taskProton->SetMergePosNeg(1);
	taskProton->SetHarmonics(2);
	taskProton->DoCumOrderMax(kTwo);
	process->AddTask(taskProton);

	FlowTask* taskProton2 = new FlowTask(kProton);
	taskProton2->SetNumSamples(1);
	taskProton2->SetEtaGap(0.0);
	taskProton2->SetPtBins(vecPtBins);
	taskProton2->SetMergePosNeg(1);
	taskProton2->SetHarmonics(2);
	taskProton2->DoCumOrderMax(kTwo);
	process->AddTask(taskProton2);

	FlowTask* taskProton3 = new FlowTask(kProton);
	taskProton3->SetNumSamples(1);
	taskProton3->SetEtaGap(0.8);
	taskProton3->SetPtBins(vecPtBins);
	taskProton3->SetMergePosNeg(1);
	taskProton3->SetHarmonics(2);
	taskProton3->DoCumOrderMax(kTwo);
	process->AddTask(taskProton3);

	FlowTask* taskProton4 = new FlowTask(kProton);
	taskProton4->SetNumSamples(1);
	taskProton4->SetEtaGap(1.0);
	taskProton4->SetPtBins(vecPtBins);
	taskProton4->SetMergePosNeg(1);
	taskProton4->SetHarmonics(2);
	taskProton4->DoCumOrderMax(kTwo);
	process->AddTask(taskProton4);

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
