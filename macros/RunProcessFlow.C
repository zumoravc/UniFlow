/* RunProcessFlow
 *
 * Steer macro for procesing flow results of AliAnalysisTaskFlowPID
 * See ProcessFlow class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

void RunProcessFlow()
{
	gSystem->AddIncludePath("-I/Users/vpacik/NBI/Flow/classProcess/");
	gROOT->SetMacroPath(Form("%s:/Users/vpacik/NBI/Flow/classProcess/",gROOT->GetMacroPath()));
	gROOT->LoadMacro("ProcessFlow.cpp+g"); // loading Fitter class
	//gROOT->LoadMacro("FitPID.C+g"); // 

	Short_t harmonics[] = {2,3,4};
	Short_t sizeHarmonics = sizeof(harmonics) / sizeof(harmonics[0]);

	//Double_t etaGaps[] = {1.};
	Double_t etaGaps[] = {-1.,0.,1.}; // only NoGap case
	Short_t sizeEtaGaps = sizeof(etaGaps) / sizeof(etaGaps[0]);

	Double_t cent[] = {0.,5.,10.,20,30,40,50,60,70,80}; //edges
	//Double_t cent[] = {20,30,40}; //edges
	Short_t sizeCent = sizeof(cent) / sizeof(cent[0]) - 1; 

	ProcessFlow* process = new ProcessFlow();
	process->SetDebug(0);
	process->SetInputFilePath("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/merge/");
	process->SetOutputFilePath("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/");
	process->SetOutputFileName("Flow_FB768_Old_PileON_PeriodON.root");
	process->SetTag("FB768_Old_PileON_PeriodON");
	process->SetHarmonicsArray(harmonics,sizeHarmonics);
	process->SetEtaGapsArray(etaGaps,sizeEtaGaps);
	process->SetBinsCentArray(cent,sizeCent);
	process->SetNumSamples(1);
	
	process->Run();
	
	//process->SetOutputFileName("ProcessFlow_FB768_GFK.root");
	//process->SetNumBinsCentrality(9);
	//process->SetNumSamples(5);
	//process->Run();
	

	delete process;

	return;
}