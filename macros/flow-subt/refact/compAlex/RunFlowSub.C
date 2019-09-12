/* RunF
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunFlowSub()
{
	// ##### Parameters setting ######
	// Double_t dEtaGap = 0.0;	TString sEtaGap = "gap00";
	// Double_t dEtaGap = 0.2;	TString sEtaGap = "gap02";
	// Double_t dEtaGap = 0.8;	TString sEtaGap = "gap08";

	Double_t dEtaGap = 0.4;	TString sEtaGap = "gap04";
	// Double_t dEtaGap = 0.6;	TString sEtaGap = "gap06";
	// Double_t dEtaGap = 0.8;	TString sEtaGap = "gap08";
	// Double_t dEtaGap = 1.0;	TString sEtaGap = "gap10";
	// Double_t dEtaGap = 1.2;	TString sEtaGap = "gap12";

	Double_t dMultBinning[] = {0,100};
	// Double_t dMultBinning[] = {0,10,20,40,60,100};

	printf("bins %d\n",sizeof(dMultBinning)/sizeof(dMultBinning[0]));


	// Double_t dMultBinning[] = {1,2,3,4,5,6,7,8,9,10};

	// Double_t dMultBinning[] = {0,1,2,3,4,5,6}; // fixed mult. binning

	// Double_t dPtBinningK0s[] = {3.,5.};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.2,1.,2.,3.,4.,5.,7.,10.,20.};
	Double_t dPtBinningK0s[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0}; // Alex
	// Double_t dPtBinningK0s[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,4.0,5.0,6.0,8.0,10.0}; // Alex
	// Double_t dPtBinningK0s[] = {0.2,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.,7.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.0,2.0,4.0,6.0,8.0,10.0}; // large bins
	// Double_t dPtBinningK0s[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // Run1
	TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/multiplicity-fluctuations/pp-NchRFP-gap04-17p/";
	// TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pp-run3-gaps-04-06-10-12";
	// TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pPb-run3-gap08/";
	// TString sInputPath = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pPb-run3-%s",sEtaGap.Data());

	// TString sOutputFilePath = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pPb-run3-gap08/output_vn_compAlex/%s",sEtaGap.Data());
	TString sOutputFilePath = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/multiplicity-fluctuations/pp-NchRFP-gap04-17p/output_vn_compAlex/%s",sEtaGap.Data());


	TString sOutFileName = "Processed_neg.root";

	// ##### END Parameters setting ######

	// gROOT->AddIncludePath("-I~/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~vpacik/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");


	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(1);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(dEtaGap);
	// taskCharged->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskRefs->SetMergePosNeg(0);
	taskRefs->SetInputTag("Neg");


	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(1);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(dEtaGap);
	taskCharged->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskCharged->SetMergePosNeg(0);
	taskCharged->SetInputTag("Neg");

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process->SetOutputFilePath(Form("%s/GF_eventweighted",sOutputFilePath.Data()));
	process->SetOutputFileName(sOutFileName.Data());
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->AddTask(taskRefs);
	process->AddTask(taskCharged);
	// process->AddTask(taskPion);
	// process->AddTask(taskKch);
	// process->AddTask(taskProton);
	// process->AddTask(taskK0s);
	// process->AddTask(taskLambda);
	// process->AddTask(taskPhi);
	process->Run();

	ProcessUniFlow* process_noweight = new ProcessUniFlow();
	process_noweight->SetInputFilePath(sInputPath.Data());
	process_noweight->SetInputFileName("AnalysisResults.root");
	process_noweight->SetTaskName("UniFlow");
	process_noweight->SetGlobalProfNameLabel("multScaled_");
	// process_noweight->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_noweight->SetOutputFilePath(Form("%s/GF_noneventweighted",sOutputFilePath.Data()));
	process_noweight->SetOutputFileName(sOutFileName.Data());
	process_noweight->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_noweight->AddTask(taskRefs);
	process_noweight->AddTask(taskCharged);
	// process_noweight->AddTask(taskPion);
	// process_noweight->AddTask(taskKch);
	// process_noweight->AddTask(taskProton);
	// process_noweight->AddTask(taskK0s);
	// process_noweight->AddTask(taskLambda);
	// process_noweight->AddTask(taskPhi);
	process_noweight->Run();

	ProcessUniFlow* process_sub = new ProcessUniFlow();
	process_sub->SetInputFilePath(sInputPath.Data());
	process_sub->SetInputFileName("AnalysisResults.root");
	process_sub->SetTaskName("UniFlow_sub");
	// process_sub->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub->SetOutputFilePath(Form("%s/SP_nonscaled_noneventweighted",sOutputFilePath.Data()));
	process_sub->SetOutputFileName(sOutFileName.Data());
	process_sub->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub->SetSaveMult("Mult.root");
	process_sub->AddTask(taskRefs);
	process_sub->AddTask(taskCharged);
	// process_sub->AddTask(taskPion);
	// process_sub->AddTask(taskKch);
	// process_sub->AddTask(taskProton);
	// process_sub->AddTask(taskK0s);
	// process_sub->AddTask(taskLambda);
	// process_sub->AddTask(taskPhi);
	process_sub->Run();

	ProcessUniFlow* process_sub_norm = new ProcessUniFlow();
	process_sub_norm->SetInputFilePath(sInputPath.Data());
	process_sub_norm->SetInputFileName("AnalysisResults.root");
	process_sub_norm->SetTaskName("UniFlow_sub");
	process_sub_norm->SetGlobalProfNameLabel("multScaled_");
	// process_sub_norm->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub_norm->SetOutputFilePath(Form("%s/SP_scaled_noneventweighted",sOutputFilePath.Data()));
	process_sub_norm->SetOutputFileName(sOutFileName.Data());
	process_sub_norm->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub_norm->AddTask(taskRefs);
	process_sub_norm->AddTask(taskCharged);
	// process_sub_norm->AddTask(taskPion);
	// process_sub_norm->AddTask(taskKch);
	// process_sub_norm->AddTask(taskProton);
	// process_sub_norm->AddTask(taskK0s);
	// process_sub_norm->AddTask(taskLambda);
	// process_sub_norm->AddTask(taskPhi);
	process_sub_norm->Run();

	ProcessUniFlow* process_sub_norm_weighted = new ProcessUniFlow();
	process_sub_norm_weighted->SetInputFilePath(sInputPath.Data());
	process_sub_norm_weighted->SetInputFileName("AnalysisResults.root");
	process_sub_norm_weighted->SetTaskName("UniFlow_sub");
	process_sub_norm_weighted->SetGlobalProfNameLabel("multScaled_weighted_");
	// process_sub_norm_weighted->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub_norm_weighted->SetOutputFilePath(Form("%s/SP_scaled_eventweighted",sOutputFilePath.Data()));
	process_sub_norm_weighted->SetOutputFileName(sOutFileName.Data());
	process_sub_norm_weighted->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub_norm_weighted->AddTask(taskRefs);
	process_sub_norm_weighted->AddTask(taskCharged);
	// process_sub_norm_weighted->AddTask(taskPion);
	// process_sub_norm_weighted->AddTask(taskKch);
	// process_sub_norm_weighted->AddTask(taskProton);
	// process_sub_norm_weighted->AddTask(taskK0s);
	// process_sub_norm_weighted->AddTask(taskLambda);
	// process_sub_norm_weighted->AddTask(taskPhi);
	process_sub_norm_weighted->Run();

	return;
}
