/* ProcessFlow
 *
 * Class implemented for processing AliAnalysisTaskFlowPID results.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TProfile.h"
#include "TString.h"

class ProcessFlow
{
public:
		ProcessFlow();
    ~ProcessFlow();

    void	SetInputFilePath(TString path) { fsInputFilePath = path; }
    void	SetInputFileName(TString name) { fsInputFileName = name; }
    void 	SetOutputFilePath(TString path) { fsOutputFilePath = path; }
    void	SetOutputFileName(TString name) { fsOutputFileName = name; }
    void	SetTag(TString sTag) { fsTag = sTag; }
    void	SetNumSamples(Short_t numSamples) { fiNumSamples = numSamples; }
    void	SetNumBinsCentrality(Short_t numCent) { fiNumBinsCent = numCent; }
    
    void	Run();


	private:
    TH1D*	EstimateCn2(const TH1D* hCum2); // estimate cn{2} out of <<2>>
    TH1D*	EstimateCn4(const TH1D* hCum2, const TH1D* hCum4); // estimate cn{4} out of <<2>>,<<4>>
    TH1D*	EstimateDn4(const Double_t dRef2, const TH1D* hDiff2, const TH1D* hDiff4); // estimate cn{4} out of <<2>>,<<4>>
    void	DesampleList(TList* inList, const Short_t iNumSamples); // computes mean value and sigma of flow from samples

    void	Terminate();

		TString fsInputFilePath; // path to directory with input file with AliAnalysisTaskFlowPID results
		TString fsInputFileName; // name of input file with AliAnalysisTaskFlowPID results
		TString fsOutputFilePath; // path to directory for output
		TString fsOutputFileName; // name of output file
		TString fsTag; // AnalysisTask tag to process
		Short_t fiNumSamples; // number of samples
		Short_t fiNumBinsCent; // number of centrality bins

		TFile* ffInputFile; //! input file with AliAnalysisTaskFlowPID results
		TFile* ffOutputFile; //! output file 


		// output methods
		void Error(TString sMethod, TString sMsg); // printf the msg as error
		void Warning(TString sMethod, TString sMsg); // printf the msg as warning
		void Info(TString sMethod, TString sMsg); // printf the msg as info
};



//_____________________________________________________________________________
ProcessFlow::ProcessFlow() :
	fiNumSamples(0),
	fiNumBinsCent(0)
{
	// default constructor
	fsInputFilePath = TString("/Users/vpacik/NBI/Flow/results/test");
	fsInputFileName = TString("AnalysisResults.root");
	fsOutputFilePath = TString("/Users/vpacik/NBI/Flow/results/test");
	fsOutputFileName = TString("ProcessedFlow.root");
	fsTag = TString("FB768");
	


	ffInputFile = 0x0;
	ffOutputFile = 0x0;
}
//_____________________________________________________________________________
ProcessFlow::~ProcessFlow()
{
	// default destructor
	if(ffInputFile->IsOpen())
		ffInputFile->Close();

	if(ffOutputFile->IsOpen())
		ffOutputFile->Close();

	delete ffInputFile;
	delete ffOutputFile;
}
//_____________________________________________________________________________
void ProcessFlow::Terminate()
{
	// closes open files, etc.

	if(ffInputFile->IsOpen())
		ffInputFile->Close();

	if(ffOutputFile->IsOpen())
		ffOutputFile->Close();
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

	ffOutputFile = new TFile(Form("%s/%s",fsOutputFilePath.Data(), fsOutputFileName.Data()),"RECREATE");
	if(!ffOutputFile->IsOpen())
	{
		Error("Run","Output file is not open (possibly do not exists). Terminating!");
		return;
	}

	// loading files 
	ffInputFile->cd(Form("flowPID_%s",fsTag.Data()));
	
	//TList* lTracks = (TList*) gDirectory->Get(Form("Tracks_flowPID_%s",fsTag.Data()));
	TList* lPID = (TList*) gDirectory->Get(Form("PID_flowPID_%s",fsTag.Data()));
	
	TH1D* hGFKTracksRefTwo[10] = {0x0}; // reference <<2>> = cn2
	TH1D* hGFKTracksRefFour[10] = {0x0};	// reference <<4>> 
	
	TH1D* hGFKPionsTwo[10][10] = {0x0}; // pt diff <<2'>> (pions) = dn2
	TH1D* hGFKPionsFour[10][10] = {0x0}; // pt diff <<4'>> (pions)
	
	ffOutputFile->cd();

	TProfile* profTemp = 0x0;
	for(Short_t j(0); j < fiNumSamples; j++)
	{
		profTemp = (TProfile*) lPID->FindObject(Form("fc22ReTracks_number%d",j));
		hGFKTracksRefTwo[j] = (TH1D*) profTemp->ProjectionX();
		hGFKTracksRefTwo[j]->Write(Form("hGFKTrackRefTwo_sample%d",j));

		profTemp = (TProfile*) lPID->FindObject(Form("fc42ReTracks_number%d",j));
		hGFKTracksRefFour[j] = (TH1D*) profTemp->ProjectionX();
		hGFKTracksRefFour[j]->Write(Form("hGFKTrackRefFour_sample%d",j));
		//hGFKTracksCn4[j] = EstimateCn4(hGFKTracksRefTwo[j], hGFKTracksRefFour[j]);
		//hGFKTracksCn4[j]->Write();

		for(Short_t i(0); i < fiNumBinsCent; i++)
		{
			profTemp = (TProfile*) lPID->FindObject(Form("fd22RePion_cent%d_number%d",i,j));
			hGFKPionsTwo[i][j] = (TH1D*) profTemp->ProjectionX();
			hGFKPionsTwo[i][j]->Write(Form("hGFKPionsTwo_cent%d_sample%d",i,j));

			profTemp = (TProfile*) lPID->FindObject(Form("fd42RePion_cent%d_number%d",i,j));
			hGFKPionsFour[i][j] = (TH1D*) profTemp->ProjectionX();
			hGFKPionsFour[i][j]->Write(Form("hGFKPionsFour_cent%d_sample%d",i,j));
		}
	}		

	// loaded

	// estimating 4 particle cummulants dn4, cn4
	TH1D* hGFKTracksCn4[10] = {0x0};	// reference cn4 
	TH1D* hGFKPionsDn4[10][10] = {0x0};	// pt diff dn4 (pions)
	
	for(Short_t j(0); j < fiNumSamples; j++)
	{
		hGFKTracksCn4[j] = EstimateCn4(hGFKTracksRefTwo[j], hGFKTracksRefFour[j]);
		hGFKTracksCn4[j]->Write(Form("GFKTracksCn4_sample%d",j));

		for(Short_t i(0); i < fiNumBinsCent; i++)
		{
			hGFKPionsDn4[i][j] = EstimateDn4(hGFKTracksRefTwo[j]->GetBinContent(i+1),hGFKPionsTwo[i][j],hGFKPionsFour[i][j]);
			hGFKPionsDn4[i][j]->Write(Form("GFKPionsDn4_cent%d_sample%d",i,j));
		}
	}
	
	// estimating v22, v24, v'22, v'24
	TH1D* hGFKRefv22[10] = {0x0};
	TH1D* hGFKRefv24[10] = {0x0};
	TH1D* hGFKDiffv22[10][10] = {0x0};
	TH1D* hGFKDiffv24[10][10] = {0x0};

	Double_t dValue = 0;
	for(Short_t j(0); j < fiNumSamples; j++)
	{
		hGFKRefv22[j] = (TH1D*) hGFKTracksRefTwo[j]->Clone(Form("hGFKRefv22_sample%d",j));
		hGFKRefv22[j]->SetTitle(Form("Tracks: v_{2}{2} ref noGap sample %d",j));
		hGFKRefv24[j] = (TH1D*) hGFKTracksRefFour[j]->Clone(Form("hGFKRefv24_sample%d",j));
		hGFKRefv24[j]->SetTitle(Form("Tracks: v_{2}{4} ref noGap sample %d",j));

		for(Short_t i(0); i < fiNumBinsCent; i++)
		{
			dValue = hGFKTracksRefTwo[j]->GetBinContent(i+1);
			hGFKRefv22[j]->SetBinContent(i+1,TMath::Sqrt(dValue));
			hGFKRefv22[j]->SetBinError(i+1,0);

			dValue = hGFKTracksCn4[j]->GetBinContent(i+1);
			hGFKRefv24[j]->SetBinContent(i+1,TMath::Power(-dValue,0.25));
			hGFKRefv24[j]->SetBinError(i+1,0);


			hGFKDiffv22[i][j] = (TH1D*) hGFKPionsTwo[i][j]->Clone(Form("hGFKDiffv22_sample%d",j));
			hGFKDiffv22[i][j]->SetTitle(Form("#pi: v'_{2}{2} diff cent %d sample %d",i,j));
			hGFKDiffv24[i][j] = (TH1D*) hGFKPionsDn4[i][j]->Clone(Form("hGFKDiffv24_sample%d",j));
			hGFKDiffv24[i][j]->SetTitle(Form("#pi: v'_{2}{4} diff cent %d sample %d",i,j));

			for(Int_t k(1); k < (hGFKDiffv22[0][0]->GetNbinsX() + 1); k++)
			{
				dValue = hGFKPionsTwo[i][j]->GetBinContent(k);
				dValue = dValue / hGFKRefv22[j]->GetBinContent(i+1);
				hGFKDiffv22[i][j]->SetBinContent(k,dValue);
				hGFKDiffv22[i][j]->SetBinError(k,0);

				dValue = hGFKPionsDn4[i][j]->GetBinContent(k);
				dValue = dValue / TMath::Power(hGFKRefv24[j]->GetBinContent(i+1),3);
				hGFKDiffv24[i][j]->SetBinContent(k,-dValue);
				hGFKDiffv24[i][j]->SetBinError(k,0);
			}

		  hGFKDiffv22[i][j]->Write(Form("hGFKDiffv22_cent%d_sample%d",i,j));
		  hGFKDiffv24[i][j]->Write(Form("hGFKDiffv24_cent%d_sample%d",i,j));

		}

		hGFKRefv22[j]->Write(Form("hGFKRefv22_sample%d",j));
		hGFKRefv24[j]->Write(Form("hGFKRefv24_sample%d",j));
	}

	// up to this point: validated & works well

	// desampling
	TList* lRefv22 = new TList();
	TList* lRefv24 = new TList();
	TList* lDiffv22 = new TList();
	TList* lDiffv24 = new TList();

	
	ffOutputFile->cd();
	// filling the lists 
	// reference	
	for(Short_t j(0); j < fiNumSamples; j++)
	{
		lRefv22->Add(hGFKRefv22[j]);
		lRefv24->Add(hGFKRefv24[j]);
	}

	DesampleList(lRefv22,fiNumSamples);
	DesampleList(lRefv24,fiNumSamples);

	lRefv22->Write("lPionsRefv22",TObject::kSingleKey);
	lRefv24->Write("lPionsRefv24",TObject::kSingleKey);

	// pt diff
	for(Short_t i(0); i < fiNumBinsCent; i++)
	{
		lDiffv22->Clear();
		lDiffv24->Clear();

		for(Short_t j(0); j < fiNumSamples; j++)
		{
			lDiffv22->Add(hGFKDiffv22[i][j]);
			lDiffv24->Add(hGFKDiffv24[i][j]);
		}
		
		DesampleList(lDiffv22,fiNumSamples);		
		DesampleList(lDiffv24,fiNumSamples);		
		
		lDiffv22->Write(Form("lPionsDiffv22_cent%d",i),TObject::kSingleKey);
		lDiffv24->Write(Form("lPionsDiffv24_cent%d",i),TObject::kSingleKey);
	}

	delete lRefv22;
	delete lRefv24;
	delete lDiffv22;
	delete lDiffv24;

	// writing to output file

	//profTemp->Write();


	
	Terminate();

	return;
}
//_____________________________________________________________________________
TH1D* ProcessFlow::EstimateCn2(const TH1D* hCum2)
{
	TH1D* hOut = (TH1D*) hCum2->Clone();

	const Short_t iNumBinsPt = hOut->GetNbinsX();

	Double_t dValue = 0;
	for(Int_t i(1); i < iNumBinsPt+1; i++)
	{
		dValue = hCum2->GetBinContent(i);
		hOut->SetBinContent(i,dValue);
		hOut->SetBinError(i,0);
	}

	return hOut;
}
//_____________________________________________________________________________
TH1D* ProcessFlow::EstimateCn4(const TH1D* hCum2,const TH1D* hCum4)
{
	if(!hCum2 || !hCum4)
	{
		Error("EstimateRefFour","hCum2 || hCum4 does NOT exists!");
		return 0x0;
	}

	TH1D* hOut = (TH1D*) hCum4->Clone();

	const Short_t iNumBinsPt = hOut->GetNbinsX();

	Double_t dValue2 = 0, dValue4 = 0, dOutValue = 0;
	for(Int_t i(1); i < iNumBinsPt+1; i++)
	{
		dValue2 = hCum2->GetBinContent(i);
		dValue4 = hCum4->GetBinContent(i);
		dOutValue = dValue4 - 2*TMath::Power(dValue2,2);
	
		hOut->SetBinContent(i,dOutValue);
		hOut->SetBinError(i,0);
	}

	return hOut;
}
//_____________________________________________________________________________
TH1D* ProcessFlow::EstimateDn4(const Double_t dRef2, const TH1D* hDiff2, const TH1D* hDiff4)
{
	TH1D* hOut = (TH1D*) hDiff4->Clone();
	
	const Short_t iNumBinsPt = hOut->GetNbinsX();
	
	Double_t dValue2 = 0, dValue4 = 0, dOutValue = 0;
	for(Int_t i(1); i < iNumBinsPt+1; i++)
	{
		dValue2 = hDiff2->GetBinContent(i);
		dValue4 = hDiff4->GetBinContent(i);
		dOutValue = dValue4 - 2 * dValue2 * dRef2;
	
		hOut->SetBinContent(i,dOutValue);
		hOut->SetBinError(i,0);
	}

	return hOut;
}
//_____________________________________________________________________________
void ProcessFlow::DesampleList(TList* inList, const Short_t iNumSamples)
{
	// init 
	
	TString sName = inList->First()->GetName();
	TString sTitle = inList->First()->GetTitle();

	//printf("%s\n", sName.Data());
	sName.Remove(sName.Length()-8); // remove the "_sampleX_px" from end of the name
	//printf("%s\n", sName.Data());
	sTitle.Remove(sTitle.Length()-8); // remove the "_sampleX_px" from end of the name

	TH1D* histOut = (TH1D*) inList->First()->Clone(sName.Data());
	histOut->SetTitle(sTitle.Data());

	const Short_t iNumBinsX = histOut->GetNbinsX();

	TH1D* histTemp = 0x0;
	Double_t dValue = 0, dMean = 0, dSigma = 0;
	for(Short_t i(1); i < iNumBinsX+1; i++)
	{
		dValue = 0;
		dMean = 0;
		dSigma = 0;

		// estimating mean
		for(Short_t j(0); j < iNumSamples; j++)
		{
			histTemp = (TH1D*) inList->At(j);
			dValue += histTemp->GetBinContent(i);
		}

		dMean = dValue / iNumSamples;

		// estimating error
		for(Short_t j(0); j < iNumSamples; j++)
		{
			histTemp = (TH1D*) inList->At(j);
			dSigma += TMath::Power(histTemp->GetBinContent(i) - dMean, 2);
		}

		dSigma = TMath::Sqrt(dSigma / (iNumSamples-1) );

		histOut->SetBinContent(i, dMean);
		histOut->SetBinError(i, dSigma / TMath::Sqrt(iNumSamples));
	}

	inList->Add(histOut);


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