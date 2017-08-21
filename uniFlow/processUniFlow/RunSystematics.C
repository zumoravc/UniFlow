/* RunSystematics
 *
 * Macro for processing results of systematical studies of AliAnalysisTaskUniFlow
 * See RunProcessUniFlow.C for more information
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */
 Bool_t bDoKOs = kFALSE;
 Bool_t bDoLambda = kFALSE;
 Bool_t bDoPhi = kFALSE;

void RunSystematics()
{
	gROOT->LoadMacro("ProcessUniFlow.cpp++g");



	const Short_t iNumProcess = 1; const char* sTaskTag[iNumProcess] = {"UniFlow_fb768"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit_withNUA/";
	// const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"uniflow_tpcCls80","uniflow_tpcCls90"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/tpcCls/";
	// const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"UniFlow_dcaz","UniFlow_dcaxy"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/dca/";
	// const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"UniFlow_vtx8","UniFlow_vtx9"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/Vtx_z/";
	// const Short_t iNumProcess = 3; const char* sTaskTag[iNumProcess] = {"UniFlow_CPA","UniFlow_DCA","UniFlow_DecayRadius"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/v0s/";
	// const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"uniflow_dcaZ","uniflow_NoArmernteros"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/";
	// const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"UniFlow_2sigma","UniFlow_bayes"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/pid/";
	// const Short_t iNumProcess = 1; const char* sTaskTag[iNumProcess] = {"UniFlow"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/loosePhi_tightV0s/";
	// const Short_t iNumProcess = 1; const char* sTaskTag[iNumProcess] = {"UniFlow"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/baseline_noNUA/";
	// const Short_t iNumProcess = 1; const char* sTaskTag[iNumProcess] = {"UniFlow"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/fit_fix/";

  // const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"UniFlow_vtx8","UniFlow_vtx9"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/Vtx_z_merged/";
  // const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"uniflow_tpcCls80","uniflow_tpcCls90"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/tpcCls_merged/";
  // const Short_t iNumProcess = 2; const char* sTaskTag[iNumProcess] = {"UniFlow_2sigma","UniFlow_bayes"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/pid_merged/";

 // 	const char* sOutputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/comparison"
 // const Short_t iNumProcess = 1; const char* sTaskTag[iNumProcess] = {"uniflow_Alex"}; const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/";

	bDoKOs = kTRUE;
	bDoLambda = kTRUE;
	bDoPhi = kTRUE;


	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath);
	process->SetInputFileName("AnalysisResults.root");

	// setting multiplicity binning
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100};
	Double_t dMultBinning[] = {0,20,40,60,100};
	// Double_t dMultBinning[] = {0,10};
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));

		// Double_t dPt[] = {1.,1.2,1.4,1.6};
		Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.};
		// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
		// Double_t dPt[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
		Double_t dPtKaon[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
		Double_t dPtProton[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,9.,10.,15.,20.};
		// Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
		// Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
		// Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1

		// Double_t dPtBinningK0s[] = {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.};
		// Double_t dPtBinningLambda[] = {0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.};
		// Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.,8.};


		// Double_t dPt[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1
		// Double_t dPtKaon[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // HEP Run /1
		// Double_t dPtProton[] = {0.5,0.75,1.,1.25,1.5,2.,2.5,3.,4.,}; // HEP Run 1


    // Double_t dPtBinningK0s[] = {2.,4.};
  	// Double_t dPtBinningK0s[] = {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.6,3.,3.5,4.,5.,6.};
  	// Double_t dPtBinningLambda[] = {0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.7,4.5,6.};
  	Double_t dPtBinningK0s[] = {0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.};
  	// Double_t dPtBinningLambda[] = {0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.8,4.6,6.};
  	// Double_t dPtBinningLambda[] = {0.,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.6.,4.3.,5.295};
  	Double_t dPtBinningLambda[] = {0.8,1.1,1.4,1.6,1.8,2.,2.2,2.4,2.7,3.2,3.6.,4.3.,5.295};
  	Double_t dPtBinningPhi[] = {0.5,1.,1.5,2.,2.5,3.,4.,6.};


		FlowTask* taskRefs = new FlowTask("Refs",FlowTask::kRefs);
		taskRefs->SetHarmonics(2);
		taskRefs->SetEtaGap(0.8);
		taskRefs->SetNumSamples(1);
		taskRefs->SetMergePosNeg(1);
		process->AddTask(taskRefs);

		FlowTask* taskK0s = new FlowTask("K0s",FlowTask::kK0s);
		taskK0s->SetHarmonics(2);
		taskK0s->SetEtaGap(0.8);
		taskK0s->SetInvMassRebin(2);
		taskK0s->SetFlowMassRebin(2);
		taskK0s->SetFlowFitFixTerms(kFALSE);
		// task7->SetShowMultDist(kTRUE);
		// taskK0s->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
		taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
		// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Pos");
		taskK0s->SetFittingRejectNumSigmas(7);
		taskK0s->SetFittingRange(0.44,0.56);
		taskK0s->SetMergePosNeg(1);
		taskK0s->SetFlowFitFixTerms(kFALSE);
		// task7->SuggestPtBinning(1,30000);
		if(bDoKOs) process->AddTask(taskK0s);

		FlowTask* taskLambda = new FlowTask("Lambda",FlowTask::kLambda);
		taskLambda->SetHarmonics(2);
		taskLambda->SetEtaGap(0.8);
		taskLambda->SetInvMassRebin(2);
		taskLambda->SetFlowMassRebin(2);
		taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
		taskLambda->SetFittingRange(1.099,1.147);
		// taskLambda->SetFittingRejectNumSigmas(3);
		// taskLambda->SetAlternativeProfileName("fp3V0sCorrLambda_<2>_harm2_gap08_Pos");
		taskLambda->SetFittingRejectNumSigmas(9);
		taskLambda->SetMergePosNeg(1);
		taskLambda->SetFlowFitFixTerms(0);
		if(bDoLambda) process->AddTask(taskLambda);

		FlowTask* taskPhi = new FlowTask("Phi",FlowTask::kPhi);
		taskPhi->SetHarmonics(2);
		taskPhi->SetEtaGap(0.8);
		taskPhi->SetInvMassRebin(2);
		taskPhi->SetFlowMassRebin(2);
		taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
		// taskPhi->SetAlternativeProfileName("fp3PhiCorr_<2>_harm2_gap08_Pos");
		taskPhi->SetFittingRange(1.001,1.059);
		taskPhi->SetMergePosNeg(1);
		taskPhi->SetFlowFitFixTerms(0);
		taskPhi->SetFlowPhiSubtLS(1);
		if(bDoPhi) process->AddTask(taskPhi);

		FlowTask* taskCharged = new FlowTask("Charged",FlowTask::kCharged);
		taskCharged->SetHarmonics(2);
		taskCharged->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
		taskCharged->SetEtaGap(0.8);
		// taskCharged->SetAlternativeProfileName("fp2Charged_<2>_harm2_gap08_Pos");
		taskCharged->SetMergePosNeg(1);
		taskCharged->SetNumSamples(1);
		process->AddTask(taskCharged);

		FlowTask* taskPion = new FlowTask("Pi",FlowTask::kPion);
		taskPion->SetHarmonics(2);
		taskPion->SetPtBins(dPt,sizeof(dPt)/sizeof(dPt[0]));
		taskPion->SetEtaGap(0.8);
		taskPion->SetMergePosNeg(1);
		taskPion->SetNumSamples(1);
		// process->AddTask(taskPion);

		FlowTask* taskKaon = new FlowTask("K",FlowTask::kKaon);
		taskKaon->SetHarmonics(2);
		taskKaon->SetPtBins(dPtKaon,sizeof(dPtKaon)/sizeof(dPtKaon[0]));
		taskKaon->SetEtaGap(0.8);
		taskKaon->SetMergePosNeg(1);
		taskKaon->SetNumSamples(1);
		// process->AddTask(taskKaon);

		FlowTask* taskProton = new FlowTask("P",FlowTask::kProton);
		taskProton->SetHarmonics(2);
		taskProton->SetPtBins(dPtProton,sizeof(dPtProton)/sizeof(dPtProton[0]));
		// taskProton->SetAlternativeProfileName("fp2Proton_<2>_harm2_gap08_Pos");
		taskProton->SetEtaGap(0.8);
		taskProton->SetMergePosNeg(1);
		taskProton->SetNumSamples(1);
		// process->AddTask(taskProton);

	// task dependent
	for(Short_t iProcess(0); iProcess < iNumProcess; iProcess++)
	// for(Short_t iProcess(0); iProcess < 1; iProcess++)
	{
		// process->SetTaskName(Form("uniflow_%s",sTaskTag[iProcess]));
		process->SetTaskName(sTaskTag[iProcess]);
		process->SetOutputFileName(Form("UniFlow_%s.root",sTaskTag[iProcess]));
		process->SetOutputFilePath(Form("%s/%s",sInputPath,sTaskTag[iProcess]));

		process->Run();
		// process->Clear();

	}

	return;
}
