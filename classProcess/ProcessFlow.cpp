/* ProcessFlow
 *
 * Class implemented for processing AliAnalysisTaskFlowPID results.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "TFile.h"
#include "TList.h"
#include "TString.h"

class ProcessFlow
{
public:
		ProcessFlow();
    ~ProcessFlow();

    void	SetInputFilePath(TString path) { fsInputFilePath = path; }
    void	SetInputFileName(TString name) { fsInputFileName = name; }
    void	SetTag(TString sTag) { fsTag = sTag; }

    void	Run();
	private:
		TString fsInputFilePath; // path to directory with input file with AliAnalysisTaskFlowPID results
		TString fsInputFileName; // name of input file with AliAnalysisTaskFlowPID results
		TString fsTag; // AnalysisTask tag to process

		TFile* ffInputFile; //! input file with AliAnalysisTaskFlowPID results

		// output methods
		void Error(TString sMethod, TString sMsg); // printf the msg as error
		void Warning(TString sMethod, TString sMsg); // printf the msg as warning
		void Info(TString sMethod, TString sMsg); // printf the msg as info
};

//_____________________________________________________________________________
ProcessFlow::ProcessFlow()
{
	// default constructor
	fsInputFilePath = TString("/Users/vpacik/NBI/Flow/results/test");
	fsInputFileName = TString("AnalysisResults.root");
	ffInputFile = 0x0;
}
//_____________________________________________________________________________
ProcessFlow::~ProcessFlow()
{
	// default destructor
	if(ffInputFile->IsOpen())
		ffInputFile->Close();
}
//_____________________________________________________________________________
void ProcessFlow::Run()
{
	ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(), fsInputFileName.Data()),"READ");
	if(!ffInputFile->IsOpen())
	{
		Error("Run","Input file is not open (possibly do not exists). Terminating!");
		return;
	}

	// loading files
	ffInputFile->cd(Form("flowPID_%s",fsTag.Data()));
	
	TList* lTracks = (TList*) gDirectory->Get(Form("Tracks_flowPID_%s",fsTag.Data()));
	


	
	return;
}
//_____________________________________________________________________________
void ProcessFlow::Error(TString sMethod, TString sMsg)
{
	printf("E-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessFlow::Info(TString sMethod, TString sMsg)
{
	printf("I-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessFlow::Warning(TString sMethod, TString sMsg)
{
	printf("W-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}