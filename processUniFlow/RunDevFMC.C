// NEED to run .X Load.C beforehand

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunDevFMC()
{
	Int_t iNumSamples = 1;
	Double_t dEtaGap = 0.0	;
	Bool_t hijingData = kFALSE;

	TString sInputPath = "/home/alidock/ana/output/LHC15o/train_7387";
	TString sOutputFilePath = sInputPath + "/processUniFlow";

	std::vector<Double_t> dMultBinning = {0,5,10,20,30,40,50,60,70,80};
	// ##### END Parameters setting ######

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlownogap");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed_NOGAP.root");
	process->SetMultiplicityBins(dMultBinning);
	process->SetDebug(1);
	process->SetSaveInterSteps(1);

	FlowTask* taskRefs = new FlowTask(kRefs, "");
	taskRefs->SetNumSamples(10);
	//taskRefs->SetEtaGap(1.0);
	//taskRefs->SetPtBins(vecPtBins);
	taskRefs->DoSixFMC({2,3,3,-2,-3,-3});
	// taskRefs->SetEtaGap(0.0);
	taskRefs->IsHijing(hijingData);
	//taskRefs->SetHarmonics(2);
	// taskRefs->DoCumOrderMax(kTwo);
	//taskRefs->SetRebinning(kFALSE);
	process->AddTask(taskRefs);

	FlowTask* taskSC23 = new FlowTask(kRefs, "");
	taskSC23->SetNumSamples(10);
	taskSC23->DoFourFMC({2,3,-2,-3});
	// taskSC23->SetEtaGap(0.0);
	taskSC23->IsHijing(hijingData);
	process->AddTask(taskSC23);

	FlowTask* taskSC24 = new FlowTask(kRefs, "");
	taskSC24->SetNumSamples(10);
	taskSC24->DoFourFMC({2,4,-2,-4});
	// taskSC24->SetEtaGap(0.0);
	taskSC24->IsHijing(hijingData);
	process->AddTask(taskSC24);

	FlowTask* taskSC34 = new FlowTask(kRefs, "");
	taskSC34->SetNumSamples(10);
	taskSC34->DoFourFMC({3,4,-3,-4});
	// taskSC34->SetEtaGap(0.0);
	taskSC34->IsHijing(hijingData);
	process->AddTask(taskSC34);

	FlowTask* taskRefs2 = new FlowTask(kRefs, "");
	taskRefs2->SetNumSamples(10);
	taskRefs2->SetMergePosNeg(1);
	taskRefs2->DoSixFMC({2,2,3,-2,-2,-3});
	// taskRefs2->SetEtaGap(0.0);
	taskRefs2->IsHijing(hijingData);
	process->AddTask(taskRefs2);

	FlowTask* taskRefsSixF = new FlowTask(kRefs, "");
	taskRefsSixF->SetNumSamples(10);
	taskRefsSixF->DoSixFMC({2,3,4,-2,-3,-4});
	// taskRefsSixF->SetEtaGap(0.0);
	taskRefsSixF->IsHijing(hijingData);
	process->AddTask(taskRefsSixF);

	FlowTask* taskRefs3 = new FlowTask(kRefs, "");
	taskRefs3->SetNumSamples(10);
	taskRefs3->DoEightFMC({2,2,3,3,-2,-2,-3,-3});
	// taskRefs3->SetEtaGap(0.0);
	taskRefs3->IsHijing(hijingData);
	process->AddTask(taskRefs3);

	FlowTask* taskRefs4 = new FlowTask(kRefs, "");
	taskRefs4->SetNumSamples(10);
	taskRefs4->DoEightFMC({2,2,2,3,-2,-2,-2,-3});
	// taskRefs4->SetEtaGap(0.0);
	taskRefs4->IsHijing(hijingData);
	process->AddTask(taskRefs4);

	FlowTask* taskRefs5 = new FlowTask(kRefs, "");
	taskRefs5->SetNumSamples(10);
	taskRefs5->DoEightFMC({2,3,3,3,-2,-3,-3,-3});
	// taskRefs5->SetEtaGap(0.0);
	taskRefs5->IsHijing(hijingData);
	process->AddTask(taskRefs5);

	FlowTask* taskRefsSix2 = new FlowTask(kRefs, "");
	taskRefsSix2->SetNumSamples(10);
	taskRefsSix2->DoSixFMC({3,4,5,-3,-4,-5});
	// taskRefsSix2->SetEtaGap(0.0);
	taskRefsSix2->IsHijing(hijingData);
	process->AddTask(taskRefsSix2);



	process->Run();


	return;
}
