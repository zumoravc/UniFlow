/* RunSystematics
 *
 * Macro for processing results of systematical studies of AliAnalysisTaskUniFlow
 * See RunProcessUniFlow.C for more information
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunSystematics()
{
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");

	const Short_t iNumProcess = 2;
	const char* sTaskTag[iNumProcess] = {"vtx8","vtx9"};


	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath("/Users/vpacik/NBI/Flow/results/uniFlow_syst/Vtx_z");
	process->SetInputFileName("AnalysisResults.root");

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	Double_t dMultBinning[] = {0,20};
	// Double_t dMultBinning[] = {0,10};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	// Double_t dPt[] = {1.,1.2,1.4,1.6};
	Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	Double_t dPtKaon[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	Double_t dPtProton[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
	// Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
	// Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
	// Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1

	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningLambda[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,10.,20.};

	FlowTask* taskRefs = new FlowTask("Refs",FlowTask::kRefs);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(0.8);
	process->AddTask(taskRefs);

	FlowTask* taskPion = new FlowTask("Pi",FlowTask::kPion);
	taskPion->SetHarmonics(2);
	taskPion->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskPion->SetEtaGap(0.8);
	// process->AddTask(taskPion);

	FlowTask* taskKaon = new FlowTask("K",FlowTask::kKaon);
	taskKaon->SetHarmonics(2);
	taskKaon->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
	taskKaon->SetEtaGap(0.8);
	// process->AddTask(taskKaon);

	FlowTask* taskProton = new FlowTask("P",FlowTask::kProton);
	taskProton->SetHarmonics(2);
	taskProton->SetPtBins(dPtProton,sizeof(dPtProton)/sizeof(dPtProton[0]));
	// taskProton->SetAlternativeProfileName("fp2Proton_<2>_harm2_gap08_Pos");
	taskProton->SetEtaGap(0.8);
	// process->AddTask(taskProton);

	FlowTask* taskCharged = new FlowTask("Charged",FlowTask::kCharged);
	taskCharged->SetHarmonics(2);
	taskCharged->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
	taskCharged->SetEtaGap(0.8);
	// taskCharged->SetAlternativeProfileName("fp2Charged_<2>_harm2_gap08_Pos");
	taskCharged->SetMergePosNeg();
	taskCharged->SetNumSamples(1);
	process->AddTask(taskCharged);

	FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(0.8);
	taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	// task7->SetShowMultDist(kTRUE);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
	taskK0s->SetMergePosNeg();
	// task7->SuggestPtBinning(1,30000);
	// process->AddTask(taskK0s);

	FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(0.8);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Pos");
	taskLambda->SetMergePosNeg();
	// process->AddTask(taskLambda);


	FlowTask* taskPhi = new FlowTask("Phi",FlowTask::kPhi);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(0.8);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	// taskPhi->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Pos");
	taskPhi->SetMergePosNeg();
	// process->AddTask(taskPhi);




	// task dependent

	for(Short_t iProcess(0); iProcess < iNumProcess; iProcess++)
	// for(Short_t iProcess(0); iProcess < 1; iProcess++)
	{
		process->SetTaskName(Form("UniFlow_%s",sTaskTag[iProcess]));
		process->SetOutputFileName(Form("UniFlow_%s.root",sTaskTag[iProcess]));
		process->SetOutputFilePath(Form("/Users/vpacik/NBI/Flow/results/uniFlow_syst/Vtx_z/syst/%s",sTaskTag[iProcess]));

		process->Run();
		// process->Clear();

	}

	return;
}
