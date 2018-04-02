/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunSystematics_v0s()
{
	Int_t iNumSamples = 1;
	Double_t dEtaGaps[] = {0.0,0.4,0.8};
	TString sEtaGaps[] = {"gap00","gap04","gap08"};
	Int_t iNumGaps = sizeof(dEtaGaps)/sizeof(dEtaGaps[0]);

	TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/v0s/merged-16qt";

	TString sTags[] = {"3sigma","decayRad","CPA","DCAdaughters"};
	Int_t iNumTags = 4;

	Double_t dMultBinning[] = {0,10,20,40,60,100};

	// old binning
	// Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningProton[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningLambda[] = {0.3,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.2,5.0,6.0,7.0};
	// Double_t dPtBinningPhi[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.0};

	// most recent
	Double_t dPtBinningPID[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,7.0};
	Double_t dPtBinningProton[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.0,3.5,4.2,5.0,7.0};
	Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0,2.5,3.0,3.5,4.0,5.0,7.0};
	Double_t dPtBinningLambda[] = {0.3,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.0,5.0,7.0};
	Double_t dPtBinningPhi[] = {0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0};

	// Binning Run1 SP
	// Double_t dMultBinning[] = {0,20,40,60,100};
	// Double_t dPtBinningPID[] = {0.3,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,4.0};
	// Double_t dPtBinningProton[] = {0.3,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,4.0};
	// Double_t dPtBinningK0s[] = {0.3,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,4.0};
	// Double_t dPtBinningLambda[] = {0.3,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,4.0};
	// Double_t dPtBinningPhi[] = {0.3,0.5,0.75,1.0,1.25,1.5,2.0,2.5,3.0,4.0};


	// ##### END Parameters setting ######

	// gSystem->AddIncludePath("-I/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");

	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	taskRefs->SetHarmonics(2);
	// taskRefs->SetEtaGap(dEtaGap);
	taskRefs->SetMergePosNeg();


	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetHarmonics(2);
	// taskCharged->SetEtaGap(dEtaGap);
	taskCharged->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskCharged->SetMergePosNeg();


	FlowTask* taskPion = new FlowTask(FlowTask::kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetHarmonics(2);
	// taskPion->SetEtaGap(dEtaGap);
	taskPion->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskPion->SetMergePosNeg();

	FlowTask* taskKch = new FlowTask(FlowTask::kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(2);
	// taskKch->SetEtaGap(dEtaGap);
	taskKch->SetPtBins(dPtBinningPID,sizeof(dPtBinningPID)/sizeof(dPtBinningPID[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskProton = new FlowTask(FlowTask::kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetHarmonics(2);
	// taskProton->SetEtaGap(dEtaGap);
	taskProton->SetPtBins(dPtBinningProton,sizeof(dPtBinningProton)/sizeof(dPtBinningProton[0]));
	taskProton->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	taskK0s->SetFittingOneGo(kTRUE);
	taskK0s->SetNumSamples(iNumSamples);
	taskK0s->SetHarmonics(2);
	// taskK0s->SetEtaGap(dEtaGap);
	// taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskLambda = new FlowTask(FlowTask::kLambda	);
	taskLambda->SetFittingOneGo(kTRUE);
	taskLambda->SetNumSamples(iNumSamples);
	taskLambda->SetHarmonics(2);
	// taskLambda->SetEtaGap(dEtaGap);
	// taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	taskLambda->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskPhi = new FlowTask(FlowTask::kPhi);
	taskPhi->SetFittingOneGo(kTRUE);
	taskPhi->SetNumSamples(iNumSamples);
	taskPhi->SetHarmonics(2);
	// taskPhi->SetEtaGap(dEtaGap);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");


	for(Int_t iGap(0); iGap < iNumGaps; ++iGap)
	{
		TString sEtaGap = sEtaGaps[iGap];
		Double_t dEtaGap = dEtaGaps[iGap];

		taskRefs->SetEtaGap(dEtaGap);
		taskCharged->SetEtaGap(dEtaGap);
		taskPion->SetEtaGap(dEtaGap);
		taskKch->SetEtaGap(dEtaGap);
		taskProton->SetEtaGap(dEtaGap);
		taskK0s->SetEtaGap(dEtaGap);
		taskLambda->SetEtaGap(dEtaGap);
		taskPhi->SetEtaGap(dEtaGap);

		for(Int_t iTag(0); iTag < iNumTags; ++iTag)
		{
			TString sTag = sTags[iTag];

			TString sOutputFilePath = sInputPath+"/output_binning-3/"+sEtaGap+"/"+sTag+"/";

			ProcessUniFlow* process = new ProcessUniFlow();
			process->SetInputFilePath(sInputPath.Data());
			process->SetInputFileName("AnalysisResults.root");
			process->SetTaskName(Form("UniFlow_%s",sTag.Data()));
			process->SetOutputFilePath(sOutputFilePath.Data());
			process->SetOutputFileName("Processed.root");
			process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
			process->SetSaveMult(kTRUE);
			process->SetFitCumulants(kTRUE);
			process->SetDebug();
			process->AddTask(taskRefs);
			// process->AddTask(taskCharged);
			// process->AddTask(taskPion);
			// process->AddTask(taskKch);
			// process->AddTask(taskProton);
			process->AddTask(taskK0s);
			process->AddTask(taskLambda);
			// process->AddTask(taskPhi);
			process->Run();
		}
	}





	return;
}
