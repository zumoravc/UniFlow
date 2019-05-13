/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunSystematics_fitting()
{
	Int_t iNumSamples = 10;
	// Double_t dEtaGaps[] = {0.0,0.4,0.8};
	// TString sEtaGaps[] = {"gap00","gap04","gap08"};
	// Int_t iNumGaps = sizeof(dEtaGaps)/sizeof(dEtaGaps[0]);

	Double_t dEtaGaps[] = {0.4}; TString sEtaGaps[] = {"gap04"};
	Double_t dEtaGap = dEtaGaps[0];	TString sEtaGap = sEtaGaps[0];

	TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/fitting/merged-pPb-16qt-nua";

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
	taskRefs->SetEtaGap(dEtaGap);
	taskRefs->SetMergePosNeg();

 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	// taskK0s->SetFittingOneGo(kTRUE);
	taskK0s->SetNumSamples(iNumSamples);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(dEtaGap);
	// taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskLambda = new FlowTask(FlowTask::kLambda	);
	// taskLambda->SetFittingOneGo(kTRUE);
	taskLambda->SetNumSamples(iNumSamples);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(dEtaGap);
	// taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningLambda,sizeof(dPtBinningLambda)/sizeof(dPtBinningLambda[0]));
	taskLambda->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskPhi = new FlowTask(FlowTask::kPhi);
	// taskPhi->SetFittingOneGo(kTRUE);
	taskPhi->SetNumSamples(iNumSamples);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(dEtaGap);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName(Form("UniFlow"));
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
	// process->AddTask(taskK0s);
	// process->AddTask(taskLambda);
	process->AddTask(taskPhi);


		// // mass range K0s
		// Double_t dMin = 0.4;
		// Double_t dMax = 0.6;
		// Int_t iNumIt = 7;

		// // // mass range Lambda
		// Double_t dMin = 1.096; // 1.106
		// Double_t dMax = 1.150;
		// Int_t iNumIt = 5;

		// mass range Phi
		Double_t dMin = 0.994;
		Double_t dMax= 1.07;
		Int_t iNumIt = 6;

		// TString sTag = "test";
		// TString sOutputFilePath = sInputPath+"/output_binning-3/"+sEtaGap+"/Lambda/"+sTag+"/";
		// process->SetOutputFilePath(sOutputFilePath.Data());
		// process->Run();


		// for(Int_t it(1); it < iNumIt+1; ++it)
		// {
		// 	sTag = Form("range%d",it);
		//
		// 	// taskK0s->SetFitMassRange(dMin+it*0.01, dMax-it*0.01); TString sOutputFilePath = sInputPath+"/output_binning-3/"+sEtaGap+"/K0s/"+sTag+"/";
		// 	// taskLambda->SetFitMassRange(dMin+it*0.002, dMax-it*0.004); TString sOutputFilePath = sInputPath+"/output_binning-3/"+sEtaGap+"/Lambda/"+sTag+"/";
		// 	taskPhi->SetFitMassRange(dMin+it*0.002, dMax-it*0.005); TString sOutputFilePath = sInputPath+"/output_binning-3/"+sEtaGap+"/Phi/"+sTag+"/";
		//
		// 	process->SetOutputFilePath(sOutputFilePath.Data());
		// 	process->Run();
		// }


	// default
	// TString sTag = "def";
	// taskK0s->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskK0s->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])",5);
	// taskK0s->SetFitFlowBG("[9]*x+[10]",2);
	//
	// Double_t dDefK0s[] =      {1.0,1.0,1.0,1.0,   10000,0.4976,0.003,10000,0.01, 1.0,1.0};
	// Double_t dLimLowK0s[] =   {-1,-1,-1,-1,    -1,0.48,0.003,-1,0.003, -1,-1};
	// Double_t dLimHighK0s[] =  {-1,-1,-1,-1,  -1,0.52,0.006,-1,0.01,  -1,-1};
	//
	// taskK0s->SetFitParDefaults(dDefK0s,11);
 	// taskK0s->SetFitParLimitsLow(dLimLowK0s,11);
 	// taskK0s->SetFitParLimitsHigh(dLimHighK0s,11);

	// TString sTag = "def";
	// taskLambda->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskLambda->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])",5);
	// taskLambda->SetFitFlowBG("[9]*x+[10]",2);
	//
	// Double_t dDef[] =      {1.0,1.0,1.0,1.0,   5000,1.115, 0.008,5000,0.01, 1.0,1.0};
	// Double_t dLimLow[] =   {-1,-1,-1,-1,    -1,1.10,0.001,-1,0.001, -1,-1};
	// Double_t dLimHigh[] =  {-1,-1,-1,-1,  -1,1.13,0.008,-1,0.01,  -1,-1};
	//
	// taskLambda->SetFitParDefaults(dDef,11);
 	// taskLambda->SetFitParLimitsLow(dLimLow,11);
 	// taskLambda->SetFitParLimitsHigh(dLimHigh,11);

	// TString sTag = "def";
	// taskPhi->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskPhi->SetFitMassSig("[4]*TMath::BreitWigner(x,[5],[6])",3);
	// taskPhi->SetFitFlowBG("[7]*x+[8]",2);
	//
	// Double_t dDef[] =      {1.0,1.0,1.0,1.0,   1000,1.019445,0.0046, 1.0,1.0};
	// Double_t dLimLow[] =   {-1,-1,-1,-1,    -1,1.018,0.001, -1,-1};
	// Double_t dLimHigh[] =  {-1,-1,-1,-1,  -1,1.022,0.006,  -1,-1};
	//
	//
	// taskPhi->SetFitParDefaults(dDef,9);
 	// taskPhi->SetFitParLimitsLow(dLimLow,9);
 	// taskPhi->SetFitParLimitsHigh(dLimHigh,9);

	// single gauss
	// TString sTag = "massSig";
	// taskK0s->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskK0s->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])",3);
	// taskK0s->SetFitFlowBG("[7]*x+[8]",2);
	//
	// Double_t dDefK0s[] =      {1.0,1.0,1.0,1.0,   10000,0.4976,0.003, 1.0,1.0};
	// Double_t dLimLowK0s[] =   {-1,-1,-1,-1,    -1,0.48,0.003, -1,-1};
	// Double_t dLimHighK0s[] =  {-1,-1,-1,-1,  -1,0.52,0.006,  -1,-1};
	//
	// taskK0s->SetFitParDefaults(dDefK0s,9);
 	// taskK0s->SetFitParLimitsLow(dLimLowK0s,9);
 	// taskK0s->SetFitParLimitsHigh(dLimHighK0s,9);

	// TString sTag = "massSig";
	// taskLambda->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskLambda->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])",3);
	// taskLambda->SetFitFlowBG("[7]*x+[8]",2);
	//
	//
	// Double_t dDef[] =      {1.0,1.0,1.0,1.0,   10000,1.115,0.001, 1.0,1.0};
	// Double_t dLimLow[] =   {-1,-1,-1,-1,    -1,1.10,0.001, -1,-1};
	// Double_t dLimHigh[] =  {-1,-1,-1,-1,  -1,1.13,0.008,  -1,-1};
	//
	// taskLambda->SetFitParDefaults(dDef,9);
 	// taskLambda->SetFitParLimitsLow(dLimLow,9);
 	// taskLambda->SetFitParLimitsHigh(dLimHigh,9);

	TString sTag = "massSig";
	taskPhi->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	taskPhi->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])",5);
	taskPhi->SetFitFlowBG("[9]*x+[10]",2);

	Double_t dDef[] =      {1.0,1.0,1.0,1.0,   10000,1.019445,0.0046,10000,0.01, 1.0,1.0};
	Double_t dLimLow[] =   {-1,-1,-1,-1,    -1,1.018,0.001,-1,0.001, -1,-1};
	Double_t dLimHigh[] =  {-1,-1,-1,-1,  -1,1.022,0.01,-1,0.01,  -1,-1};


	taskPhi->SetFitParDefaults(dDef,11);
	taskPhi->SetFitParLimitsLow(dLimLow,11);
	taskPhi->SetFitParLimitsHigh(dLimHigh,11);

	// massBg
	// TString sTag = "massBG";
	// taskK0s->SetFitMassBG("[0] + [1]*x + [2]*x*x",3);
	// taskK0s->SetFitMassSig("[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[4],[7])",5);
	// taskK0s->SetFitFlowBG("[8]*x+[9]",2);
	//
	// Double_t dDefK0s[] =      {1.0,1.0,1.0,   10000,0.4976,0.003,10000,0.01, 1.0,1.0};
	// Double_t dLimLowK0s[] =   {-1,-1,-1,    -1,0.48,0.003,-1,0.003, -1,-1};
	// Double_t dLimHighK0s[] =  {-1,-1,-1,  -1,0.52,0.006,-1,0.01,  -1,-1};
	//
	// taskK0s->SetFitParDefaults(dDefK0s,10);
	// taskK0s->SetFitParLimitsLow(dLimLowK0s,10);
	// taskK0s->SetFitParLimitsHigh(dLimHighK0s,10);
	//
	// TString sTag = "massBG";
	// taskLambda->SetFitMassBG("[0] + [1]*x + [2]*x*x",3);
	// taskLambda->SetFitMassSig("[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[4],[7])",5);
	// taskLambda->SetFitFlowBG("[8]*x+[9]",2);
	//
	// Double_t dDefLambda[] =      {1.0,1.0,1.0,   10000,1.115,0.001,10000,0.01, 1.0,1.0};
	// Double_t dLimLowLambda[] =   {-1,-1,-1,    -1,1.10,0.001,-1,0.001, -1,-1};
	// Double_t dLimHighLambda[] =  {-1,-1,-1,  -1,1.13,0.008,-1,0.01,  -1,-1};
	//
	// taskLambda->SetFitParDefaults(dDefLambda,10);
	// taskLambda->SetFitParLimitsLow(dLimLowLambda,10);
	// taskLambda->SetFitParLimitsHigh(dLimHighLambda,10);

	// TString sTag = "massBG";
	// taskPhi->SetFitMassBG("[0] + [1]*x + [2]*x*x",3);
	// taskPhi->SetFitMassSig("[3]*TMath::BreitWigner(x,[4],[5])",3);
	// taskPhi->SetFitFlowBG("[6]*x+[7]",2);
	//
	// Double_t dDef[] =      {1.0,1.0,1.0,   1000,1.019445,0.0046, 1.0,1.0};
	// Double_t dLimLow[] =   {-1,-1,-1,    -1,1.018,0.001, -1,-1};
	// Double_t dLimHigh[] =  {-1,-1,-1,  -1,1.022,0.006,  -1,-1};
	//
	// taskPhi->SetFitParDefaults(dDef,8);
	// taskPhi->SetFitParLimitsLow(dLimLow,8);
	// taskPhi->SetFitParLimitsHigh(dLimHigh,8);

	// // flowBG
	// TString sTag = "flowBG";
	// taskK0s->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskK0s->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])",5);
	// taskK0s->SetFitFlowBG("[9]",1);
	//
	// Double_t dDefK0s[] =      {1.0,1.0,1.0,1.0,   10000,0.4976,0.003,10000,0.01, 1.0};
	// Double_t dLimLowK0s[] =   {-1,-1,-1,-1,    -1,0.48,0.003,-1,0.003, -1};
	// Double_t dLimHighK0s[] =  {-1,-1,-1,-1,  -1,0.52,0.006,-1,0.01,  -1};
	//
	// taskK0s->SetFitParDefaults(dDefK0s,10);
 	// taskK0s->SetFitParLimitsLow(dLimLowK0s,10);
 	// taskK0s->SetFitParLimitsHigh(dLimHighK0s,10);

	// TString sTag = "flowBG";
	// taskLambda->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskLambda->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])",5);
	// taskLambda->SetFitFlowBG("[9]",1);
	//
	// Double_t dDef[] =      {1.0,1.0,1.0,1.0,   10000,1.115,0.001,10000,0.02, 0.2};
	// Double_t dLimLow[] =   {-1,-1,-1,-1,    -1,1.10,0.001,-1,0.001, -1};
	// Double_t dLimHigh[] =  {-1,-1,-1,-1,  -1,1.13,0.008,-1,0.02,  -1};
	//
	// taskLambda->SetFitParDefaults(dDef,10);
 	// taskLambda->SetFitParLimitsLow(dLimLow,10);
 	// taskLambda->SetFitParLimitsHigh(dLimHigh,10);
	//
	// TString sTag = "flowBG";
	// taskPhi->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskPhi->SetFitMassSig("[4]*TMath::BreitWigner(x,[5],[6])",3);
	// taskPhi->SetFitFlowBG("[7]",1);
	// //
	// Double_t dDef[] =      {1.0,1.0,1.0,1.0,   1000,1.019445,0.0046, 1.0};
	// Double_t dLimLow[] =   {-1,-1,-1,-1,    -1,1.018,0.001, -1};
	// Double_t dLimHigh[] =  {-1,-1,-1,-1,  -1,1.022,0.006,  -1};
	//
	// taskPhi->SetFitParDefaults(dDef,8);
 	// taskPhi->SetFitParLimitsLow(dLimLow,8);
 	// taskPhi->SetFitParLimitsHigh(dLimHigh,8);
	//

	TString sOutputFilePath = sInputPath+"/output_binning-3/"+sEtaGap+"/Phi/"+sTag+"/";
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->Run();




	// Double_t dDefLambda[] =      {1.0,1.0,1.0,1.0,   10000,0.4976,0.003,10000,0.01, 1.0,1.0};
	// Double_t dLimLowLambda[] =   {-1,-1,-1,-1,    -1,0.48,0.003,-1,0.003, -1,-1};
	// Double_t dLimHighLambda[] =  {-1,-1,-1,-1,  -1,0.52,0.006,-1,0.01,  -1,-1};
	// taskLambda->SetFitMassBG("[0] + [1]*x + [2]*x*x + [3]*x*x*x",4);
	// taskLambda->SetFitMassSig("[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])",5);
	// taskLambda->SetFitFlowBG("[9]*x+[10]",2);
	// taskLambda->SetFitParDefaults(dDefLambda,11);
	// taskLambda->SetFitParLimitsLow(dLimLowLambda,11);
	// taskLambda->SetFitParLimitsHigh(dLimHighLambda,11);


	return;
}
