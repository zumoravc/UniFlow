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
#include "TH2.h"
#include "TProfile.h"
#include "TString.h"
#include "TCanvas.h"

class ProcessFlow
{
public:
		ProcessFlow();
    ~ProcessFlow();

    void  SetDebug(Bool_t debug) { fbDebug = debug; }
    void	SetInputFilePath(TString path) { fsInputFilePath = path; }
    void	SetInputFileName(TString name) { fsInputFileName = name; }
    void 	SetOutputFilePath(TString path) { fsOutputFilePath = path; }
    void	SetOutputFileName(TString name) { fsOutputFileName = name; }
    void	SetTag(TString sTag) { fsTag = sTag; }
    void	SetNumSamples(Short_t numSamples) { fiNumSamples = numSamples; }
    void	SetBinsCentArray(Double_t* array, Short_t iSize = 0);
    void	SetHarmonicsArray(Short_t* array, Short_t iSize = 0);
    void	SetEtaGapsArray(Double_t* array, Short_t iSize = 0);
		void	SetSamplingSpread(Bool_t bSpread) { fbSamplingSpread = bSpread; }
		void  SetDoV0s(Bool_t doV0s) { fbDoV0s = doV0s; }

    void	Run(); //


	private:
		// output methods
		void Fatal(TString sMethod, TString sMsg); // printf the msg as error
		void Error(TString sMethod, TString sMsg); // printf the msg as error
		void Warning(TString sMethod, TString sMsg); // printf the msg as warning
		void Info(TString sMethod, TString sMsg); // printf the msg as info
		void Debug(TString sMethod, TString sMsg); // printf the msg as info

    Bool_t	Initialize(); // check array size, etc.
    Bool_t DesampleList(TList* inList, const Short_t iNumSamples); // estimate average of sigma from list of samples
    Bool_t ProcessRefFlow(const TList* listIn, TList* listOut, const Short_t iHarmonics, const Double_t dEtaGap); // desample reference flow (indipendent on ProcessList / not needed to run first)
    Bool_t 	ProcessList(const TList* listIn, TList* listOut, const TList* listRef, const Short_t iHarmonics, const Double_t dEtaGap, const TString sSpecies); // made flow out of cumulant list
    //TH1D*	EstimateCn2(const TH1D* hCum2); // estimate cn{2} out of <<2>>
    //TH1D*	EstimateCn4(const TH1D* hCum2, const TH1D* hCum4); // estimate cn{4} out of <<2>>,<<4>>
    //TH1D*	EstimateDn4(const Double_t dRef2, const TH1D* hDiff2, const TH1D* hDiff4); // estimate cn{4} out of <<2>>,<<4>>
    //void	DesampleList(TList* inList, const Short_t iNumSamples); // computes mean value and sigma of flow from samples

    //void	Terminate();

    // members
		//static const Short_t fiMaxNumSamples = 10; // used for array initialization
		static const Short_t fiMaxNumBinsCent = 10; // used for array initialization
		static const Short_t fiMaxNumHarmonics = 10; // used for array initialization
		static const Short_t fiMaxNumEtaGaps = 10; // used for array initialization

		TFile* ffInputFile; //! input file with AliAnalysisTaskFlowPID results
		TFile* ffOutputFile; //! output file

		TH1D* fHistTest; //! testing

		Short_t	fiNumHarmonics; // number of eta gaps
		Short_t	fiNumEtaGaps; // number of eta gaps
		Short_t	fiHarmonics[fiMaxNumHarmonics]; // array for harmonics values
		Double_t	fdEtaGaps[fiMaxNumEtaGaps]; // array for eta gaps values
		TString	fsEtaGaps[fiMaxNumEtaGaps]; // array for eta gaps strings
		Double_t	fdBinsCent[fiMaxNumBinsCent]; // array for centrality bins values / edges

    // settings
		Bool_t fbDebug; // debugging flag
		TString fsInputFilePath; // path to directory with input file with AliAnalysisTaskFlowPID results
		TString fsInputFileName; // name of input file with AliAnalysisTaskFlowPID results
		TString fsOutputFilePath; // path to directory for output
		TString fsOutputFileName; // name of output file
		TString fsTag; // AnalysisTask tag to process
		Short_t fiNumSamples; // number of samples
		Short_t fiNumBinsCent; // number of centrality bins
		Bool_t fbSamplingSpread; // flag for plotting spread of values during de-sampling
		Bool_t fbDoV0s; // flag for analysing K0s, Lambda particles



};

//_____________________________________________________________________________
ProcessFlow::ProcessFlow()
{
	// default constructor
	fbDebug = kFALSE;
	fsInputFilePath = TString("/Users/vpacik/NBI/Flow/flow");
	fsInputFileName = TString("AnalysisResults.root");
	fsOutputFilePath = TString("/Users/vpacik/NBI/Flow/temp");
	fsOutputFileName = TString("Flow.root");
	fsTag = TString("FB768");
	fiNumSamples = 0;
	fiNumBinsCent = 0;
	fiNumEtaGaps = 0;
	fbSamplingSpread = kFALSE;
	fbDoV0s = kFALSE;

	ffInputFile = 0x0;
	ffOutputFile = 0x0;

	fHistTest = 0x0;

	for(Short_t i(0); i < fiMaxNumEtaGaps; i++)
	{
		fdEtaGaps[i] = -9999;
		fsEtaGaps[i] = TString();
	}

	for(Short_t i(0); i < fiMaxNumBinsCent; i++)
		fdBinsCent[i] = -9999;

}
//_____________________________________________________________________________
ProcessFlow::~ProcessFlow()
{
	// default destructor
	if(ffInputFile && ffInputFile->IsOpen())
		ffInputFile->Close();

	if(ffOutputFile && ffOutputFile->IsOpen())
		ffOutputFile->Close();

	delete ffInputFile;
	delete ffOutputFile;
}
//_____________________________________________________________________________
void ProcessFlow::SetBinsCentArray(Double_t* array, Short_t iSize)
{
	Debug("SetBinsCentArray","Running");
	Debug("SetBinsCentArray",Form("Size:%d",iSize));
	fiNumBinsCent = iSize;

	// checking array length
	if(iSize > fiMaxNumBinsCent)
	{
		Error("SetBinsCentArray","Number of centrality bins is greater then initialized array: change value of fiMaxNumBinsCent accordingly! ");
		return;
	}

	for(Short_t i(0); i < fiNumBinsCent; i++)
	{
		fdBinsCent[i] = array[i];
		Debug("SetBinsCentArray",Form("%d : %g",i,array[i]));
	}

	return;
}
//_____________________________________________________________________________
void ProcessFlow::SetHarmonicsArray(Short_t* array, Short_t iSize)
{
	Debug("SetHarmonicsArray","Running");
	Debug("SetHarmonicsArray",Form("Size:%d",iSize));
	fiNumHarmonics = iSize;

	// checking array length
	if(iSize > fiMaxNumHarmonics)
	{
		Error("SetHarmonicsArray","Number of harmonics is greater then initialized array: change value of fiMaxNumHarmonics accordingly! ");
		return;
	}

	for(Short_t i(0); i < fiNumHarmonics; i++)
	{
		fiHarmonics[i] = array[i];
		Debug("SetHarmonicsArray",Form("%d : %d",i,array[i]));
	}

	return;
}
//_____________________________________________________________________________
void ProcessFlow::SetEtaGapsArray(Double_t* array, Short_t iSize)
{
	Debug("SetEtaGapsArray","Running");
	Debug("SetEtaGapsArray",Form("Size:%d",iSize));
	fiNumEtaGaps = iSize;

	// checking array length
	if(iSize > fiMaxNumEtaGaps)
	{
		Error("SetEtaGapsArray","Number of eta gaps is greater then initialized array: change value of fiMaxNumEtaGaps accordingly! ");
		return;
	}

	for(Short_t i(0); i < fiNumEtaGaps; i++)
	{
		fdEtaGaps[i] = array[i];
		fsEtaGaps[i] = Form("Gap%02.2g", 10*array[i]);
		Debug("SetEtaGapsArray",Form("%d : %g / %s",i,array[i], fsEtaGaps[i].Data()));
	}

	return;
}
//_____________________________________________________________________________

// //_____________________________________________________________________________
// void ProcessFlow::Terminate()
// {
// 	// closes open files, etc.

// 	if(ffInputFile->IsOpen())
// 		ffInputFile->Close();

// 	if(ffOutputFile->IsOpen())
// 		ffOutputFile->Close();
// }
//_____________________________________________________________________________
Bool_t ProcessFlow::Initialize()
{
	Debug("Initialize","Running");

	// checking array size
	if(fiNumBinsCent > fiMaxNumBinsCent)
	{
		Fatal("Initialize","Number of is greater then initialized array: change value of fiMaxNumBinsCent! ");
		return kFALSE;
	}

	if(fiNumEtaGaps > fiMaxNumEtaGaps)
	{
		Fatal("Initialize","Number of is greater then initialized array: change value of fiMaxNumEtaGaps! ");
		return kFALSE;
	}

	// openning files
	ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(), fsInputFileName.Data()),"READ");
	if(!ffInputFile->IsOpen())
	{
		Fatal("Initialize","Input file is not open (possibly do not exists).");
		return kFALSE;
	}

	ffOutputFile = new TFile(Form("%s/%s",fsOutputFilePath.Data(), fsOutputFileName.Data()),"RECREATE");
	if(!ffOutputFile->IsOpen())
	{
		Fatal("Initialize","Output file is not open (possibly do not exists).");
		return kFALSE;
	}

	return kTRUE;
}
//_____________________________________________________________________________
void ProcessFlow::Run()
{
	Debug("Run","Running");

	if(Initialize() == kFALSE)
		return;

	// loading TLists
	ffInputFile->cd(Form("flowPID_%s",fsTag.Data()));
	TList* listCumulants = (TList*) gDirectory->Get(Form("Cumulants_flowPID_%s",fsTag.Data()));
	//listCumulants->ls();

	TList* listRef[fiMaxNumEtaGaps];
	TList* listTracks[fiMaxNumEtaGaps];
	TList* listPions[fiMaxNumEtaGaps];
	TList* listKaons[fiMaxNumEtaGaps];
	TList* listProtons[fiMaxNumEtaGaps];
	TList* listK0s[fiMaxNumEtaGaps];
	TList* listLambda[fiMaxNumEtaGaps];

	for(Short_t i(0); i < fiNumEtaGaps; i++)
	{
		listRef[i] = (TList*) listCumulants->FindObject(Form("fListRef_%s",fsEtaGaps[i].Data()) );
		if(!listRef[i])
		{
			Fatal("Run",Form("List \"fListRef_%s\" does not exists",fsEtaGaps[i].Data()));
			return;
		}

		listTracks[i] = (TList*) listCumulants->FindObject(Form("fListTracks_%s",fsEtaGaps[i].Data()) );
		if(!listTracks[i])
		{
			Fatal("Run",Form("List \"fListTracks_%s\" does not exists",fsEtaGaps[i].Data()));
			return;
		}

		listPions[i] = (TList*) listCumulants->FindObject(Form("fListPions_%s",fsEtaGaps[i].Data()) );
		if(!listPions[i])
		{
			Fatal("Run",Form("List \"fListPions_%s\" does not exists",fsEtaGaps[i].Data()));
			return;
		}

		listKaons[i] = (TList*) listCumulants->FindObject(Form("fListKaons_%s",fsEtaGaps[i].Data()) );
		if(!listKaons[i])
		{
			Fatal("Run",Form("List \"fListKaons_%s\" does not exists",fsEtaGaps[i].Data()));
			return;
		}

		listProtons[i] = (TList*) listCumulants->FindObject(Form("fListProtons_%s",fsEtaGaps[i].Data()) );
		if(!listProtons[i])
		{
			Fatal("Run",Form("List \"fListProtons_%s\" does not exists",fsEtaGaps[i].Data()));
			return;
		}

		if(fbDoV0s)
		{
			listK0s[i] = (TList*) listCumulants->FindObject(Form("fListK0s_%s",fsEtaGaps[i].Data()) );
			if(!listK0s[i])
			{
				Fatal("Run",Form("List \"fListK0s_%s\" does not exists",fsEtaGaps[i].Data()));
				return;
			}

			listLambda[i] = (TList*) listCumulants->FindObject(Form("fListLambda_%s",fsEtaGaps[i].Data()) );
			if(!listLambda[i])
			{
				Fatal("Run",Form("List \"fListLambda_%s\" does not exists",fsEtaGaps[i].Data()));
				return;
			}
		}
	}

	Debug("Run","List with cumulants loaded");
	// lists loaded

	// making empty output lists for ProcessList method
	/*
	TList* listOutRef = new TList();
	TList* listOutTracks = new TList();
	TList* listOutPions = new TList();
	TList* listOutKaons = new TList();
	TList* listOutProtons = new TList();
	TList* listOutK0s = new TList();
	TList* listOutLambda = new TList();
	*/

	TList* listOutRef[fiMaxNumHarmonics][fiMaxNumEtaGaps];
	TList* listOutTracks[fiMaxNumHarmonics][fiMaxNumEtaGaps];
	TList* listOutPions[fiMaxNumHarmonics][fiMaxNumEtaGaps];
	TList* listOutKaons[fiMaxNumHarmonics][fiMaxNumEtaGaps];
	TList* listOutProtons[fiMaxNumHarmonics][fiMaxNumEtaGaps];
	TList* listOutK0s[fiMaxNumHarmonics][fiMaxNumEtaGaps];
	TList* listOutLambda[fiMaxNumHarmonics][fiMaxNumEtaGaps];

	for(Short_t i(0); i < fiNumHarmonics; i++)
	{
		for(Short_t j(0); j < fiNumEtaGaps; j++)
		{
			listOutRef[i][j] = new TList();
			listOutTracks[i][j] = new TList();
			listOutPions[i][j] = new TList();
			listOutKaons[i][j] = new TList();
			listOutProtons[i][j] = new TList();
			listOutK0s[i][j] = new TList();
			listOutLambda[i][j] = new TList();
		}
	}


	TString sSpecies = "Pion";

	Bool_t bStatusProcess = kFALSE;

	for(Short_t iHarmonics(0); iHarmonics < fiNumHarmonics; iHarmonics++)
	{

		for(Short_t iEtaGap(0); iEtaGap < fiNumEtaGaps; iEtaGap++)
		{
			Debug("","");

			// estimate (desample) reference flow
			bStatusProcess = ProcessRefFlow(listRef[iEtaGap], listOutRef[iHarmonics][iEtaGap],fiHarmonics[iHarmonics], fdEtaGaps[iEtaGap]);
			if(!bStatusProcess)
			{
				Error("Run",Form("Processing of list: %s: Gap %g n=%d: Status %d (FAILED)!", "Reference", fdEtaGaps[iEtaGap], iHarmonics, bStatusProcess));
			}

			// estiamte PID flow
			bStatusProcess = ProcessList(listTracks[iEtaGap], listOutTracks[iHarmonics][iEtaGap], listRef[iEtaGap], fiHarmonics[iHarmonics], fdEtaGaps[iEtaGap], "Tracks");
			if(bStatusProcess == kFALSE)
			{
				Error("Run",Form("Processing of list: %s: Gap %g n=%d: Status %d (FAILED)!", "Track", fdEtaGaps[iEtaGap], iHarmonics, bStatusProcess));
			}

			bStatusProcess = ProcessList(listPions[iEtaGap], listOutPions[iHarmonics][iEtaGap], listRef[iEtaGap], fiHarmonics[iHarmonics], fdEtaGaps[iEtaGap], "Pion");
			if(bStatusProcess == kFALSE)
			{
				Error("Run",Form("Processing of list: %s: Gap %g n=%d: Status %d (FAILED)!", "Pion", fdEtaGaps[iEtaGap], iHarmonics, bStatusProcess));
			}

			bStatusProcess = ProcessList(listKaons[iEtaGap], listOutKaons[iHarmonics][iEtaGap], listRef[iEtaGap], fiHarmonics[iHarmonics], fdEtaGaps[iEtaGap], "Kaon");
			if(bStatusProcess == kFALSE)
			{
				Error("Run",Form("Processing of list: %s: Gap %g n=%d: Status %d (FAILED)!", "Kaon", fdEtaGaps[iEtaGap], iHarmonics, bStatusProcess));
			}

			bStatusProcess = ProcessList(listProtons[iEtaGap], listOutProtons[iHarmonics][iEtaGap], listRef[iEtaGap], fiHarmonics[iHarmonics], fdEtaGaps[iEtaGap], "Proton");
			if(bStatusProcess == kFALSE)
			{
				Error("Run",Form("Processing of list: %s: Gap %g n=%d: Status %d (FAILED)!", "Proton", fdEtaGaps[iEtaGap], iHarmonics, bStatusProcess));
			}

			ffOutputFile->cd();

			if(fdEtaGaps[iEtaGap] < 0.) // NoGap case = Four particle cumulats are last (not the two = before last)
			{
				listOutRef[iHarmonics][iEtaGap]->Last()->Write(Form("hRef_n%d4_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ));
				listOutRef[iHarmonics][iEtaGap]->RemoveLast(); // should remove the last one from TList => now the cn2 is last
				/*
				listOutTracks[iHarmonics][iEtaGap]->Write(Form("Tracks_n%d4_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
				listOutTracks[iHarmonics][iEtaGap]->RemoveLast();

				listOutPions[iHarmonics][iEtaGap]->Write(Form("Pions_n%d4_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
				listOutPions[iHarmonics][iEtaGap]->RemoveLast();

				listOutKaons[iHarmonics][iEtaGap]->Write(Form("Kaons_n%d4_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
				listOutKaons[iHarmonics][iEtaGap]->RemoveLast();

				listOutProtons[iHarmonics][iEtaGap]->Write(Form("Protons_n%d4_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
				listOutProtons[iHarmonics][iEtaGap]->RemoveLast();
				*/
			}

			listOutRef[iHarmonics][iEtaGap]->Last()->Write(Form("hRef_n%d2_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ));
			listOutTracks[iHarmonics][iEtaGap]->Write(Form("Tracks_n%d_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
			listOutPions[iHarmonics][iEtaGap]->Write(Form("Pions_n%d_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
			listOutKaons[iHarmonics][iEtaGap]->Write(Form("Kaons_n%d_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);
			listOutProtons[iHarmonics][iEtaGap]->Write(Form("Protons_n%d_%s",fiHarmonics[iHarmonics], fsEtaGaps[iEtaGap].Data() ),TObject::kSingleKey);

		} // end of loop over eta gaps
	} // end of loop over harmonics

	for(Short_t i(0); i < fiNumHarmonics; i++)
	{
		for(Short_t j(0); j < fiNumEtaGaps; j++)
		{
			delete listOutRef[i][j];
			delete listOutTracks[i][j];
			delete listOutPions[i][j];
			delete listOutKaons[i][j];
			delete listOutProtons[i][j];
			delete listOutK0s[i][j];
			delete listOutLambda[i][j];
		}
	}



	return;
}
//_____________________________________________________________________________
Bool_t ProcessFlow::ProcessRefFlow(const TList* listIn, TList* listOut, const Short_t iHarmonics, const Double_t dEtaGap)
{
	Debug("ProcessRefFlow",Form("Processing list with %d entries. First \"%s\"",listIn->GetEntries(),listIn->First()->GetName()));

	// checking in/out lists
	if(!listIn)
	{
		Error("ProcessRefFlow","Input list does not exists!");
		return kFALSE;
	}

	//listIn->ls();

	TProfile* profRef2 = 0x0; // two particle cumulants
	TProfile* profRef4 = 0x0; // four particle cumulants
	TH1D* histTemp2 = 0x0;
	TH1D* histTemp4 = 0x0;
	TList* listTemp2 = new TList();
	TList* listTemp4 = new TList();

	Short_t iNumBinsX = 0;
	Double_t dValue = 0, dCon = 0;
	Double_t dCum2 = 0, dCum4 = 0; // <<2>>,<<4>>

	for(Short_t iSample(0); iSample < fiNumSamples; iSample++)
	{
		profRef2 = (TProfile*) listIn->FindObject(Form("fTracksRef_n%d2_gap%02.2g_number%d",iHarmonics,10*dEtaGap,iSample));
		if(!profRef2)
		{
			Error("ProcessRefFlow","Input Ref TProfile not found!");
			return kFALSE;
		}


		histTemp2 = (TH1D*) profRef2->ProjectionX()->Clone();

		// making Cn out of <<2>>
		iNumBinsX = profRef2->GetNbinsX();
		for(Int_t iBinX(1); iBinX < iNumBinsX+1; iBinX++)
		{
			dCon = profRef2->GetBinContent(iBinX);
			if(dCon > 0.)
				dValue = TMath::Sqrt(dCon);
			else
				dValue = -9.;

			histTemp2->SetBinContent(iBinX, dValue);
			histTemp2->SetBinError(iBinX, 0);
		}

		listTemp2->Add(histTemp2);

		if(dEtaGap < 0) // if NoGap case: load 4 particle cumulants
		{
			profRef4 = (TProfile*) listIn->FindObject(Form("fTracksRef_n%d4_gap%02.2g_number%d",iHarmonics,10*dEtaGap,iSample));
			if(!profRef4)
			{
				Error("ProcessRefFlow","Input Ref TProfile not found!");
				return kFALSE;
			}
			histTemp4 = (TH1D*) profRef4->ProjectionX()->Clone();

			// making Cn out of <<4>>
			iNumBinsX = profRef4->GetNbinsX();
			for(Int_t iBinX(1); iBinX < iNumBinsX+1; iBinX++)
			{
				dCum4 = profRef4->GetBinContent(iBinX);
				dCum2 = profRef2->GetBinContent(iBinX);
				dValue = dCum4 - 2*TMath::Power(dCum2,2);
				if(dValue < 0)
					histTemp4->SetBinContent(iBinX, TMath::Power(-dValue,0.25));
			 	else
					histTemp4->SetBinContent(iBinX, -9);

				histTemp4->SetBinError(iBinX, 0);
			}

			listTemp4->Add(histTemp4);
		} // end of if NoGap

	} // end of loop over samples

	if(!DesampleList(listTemp2, fiNumSamples))
	{
		Error("ProcessRefFlow",Form("Desampling cn{2} FAILED!"));
		return kFALSE;
	}

	listOut->Add(listTemp2->Last());

	if(dEtaGap < 0)
	{
		if(!DesampleList(listTemp4, fiNumSamples))
		{
			Error("ProcessRefFlow",Form("Desampling cn{4} FAILED!"));
			return kFALSE;
		}

		listOut->Add(listTemp4->Last());
	}



	delete listTemp2;
	delete listTemp4;
	//histTemp = (TH1D*) listTemp->Last()->Clone();
	//histOut = histTemp;
	//histOut = (TH1D*) listTemp->Last()->Clone();

	return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessFlow::ProcessList(const TList* listIn, TList* listOut, const TList* listRef, const Short_t iHarmonics, const Double_t dEtaGap, const TString sSpecies)
{
	Debug("ProcessList",Form("Processing list with %d entries. First \"%s\"",listIn->GetEntries(),listIn->First()->GetName()));
	// checking in/out lists
	if(!listIn)
	{
		Error("ProcessList","Input list does not exists!");
		return kFALSE;
	}

	if(!listOut)
	{
		Error("ProcessList","Output list does not exists!");
		return kFALSE;
	}

	//listRef->ls();

	// initializint tempList for histograms based on centrality
	TList* listTemp[fiNumBinsCent];
	TList* listTemp4[fiNumBinsCent];

	for(Short_t iCent(0); iCent < fiNumBinsCent; iCent++)
	{
		listTemp[iCent] = new TList();
		listTemp4[iCent] = new TList();
	}

	//listIn->ls();


	TProfile* profRef2 = 0x0;
	TProfile* profRef4 = 0x0;
	TProfile* profTemp = 0x0;
	TProfile* profRebin = 0x0;

	TProfile* profTemp4 = 0x0;
	TH1D* histTemp = 0x0;
	TH1D* histTemp4 = 0x0;
	Short_t iNumBinsX = 0;
	Double_t dRef = 0, dValue = 0;
	Double_t dCum2 = 0, dCum4 = 0, dRefCum2 = 0, dRefCum4 = 0; // <<2'>>, <<4'>>, <<2>>, <<4>>

	TString sRapSign[2] = {"Pos","Neg"};

	for(Short_t iSample(0); iSample < fiNumSamples; iSample++)
	{
		profRef2 = (TProfile*) listRef->FindObject(Form("fTracksRef_n%d2_gap%02.2g_number%d",iHarmonics,10*dEtaGap,iSample));
		if(!profRef2)
		{
			Error("ProcessList","Input Ref <<2>> TProfile not found!");
			return kFALSE;
		}

		if(dEtaGap < 0.) // no gap case
		{
			profRef4 = (TProfile*) listRef->FindObject(Form("fTracksRef_n%d4_gap%02.2g_number%d",iHarmonics,10*dEtaGap,iSample));
			if(!profRef4)
			{
				Error("ProcessList","Input Ref <<4>> TProfile not found!");
				return kFALSE;
			}
		} // end of NoGap case

		for(Short_t iRap(0); iRap < 2; iRap++)
		{
			if(iRap == 1 && dEtaGap < 0) // for NoGap case skip Neg POIs
				continue;

			for(Short_t iCent(0); iCent < fiNumBinsCent; iCent++)
			{
				dRefCum2 = profRef2->GetBinContent(iCent+1);
				dRef = TMath::Sqrt(dRefCum2); // getting reference cummulant

				profTemp = (TProfile*) listIn->FindObject(Form("f%s_n%d2_%s_gap%02.2g_cent%d_number%d",sSpecies.Data(),iHarmonics,sRapSign[iRap].Data(),10*dEtaGap,iCent,iSample));
				if(!profTemp)
				{
					Error("ProcessList",Form("Input POIs %s <<2>> TProfile not found!",sSpecies.Data()));
					return kFALSE;
				}

				profTemp->SetName(Form("f%s_n%d2_gap%02.2g_cent%d_number%d_%d",sSpecies.Data(),iHarmonics,10*dEtaGap,iCent,iSample,iRap));

				/*
				// make projections on pT / rebinning
				Double_t dNewBinEdges[] = {0};
				Int_t iNumNewBins = dNewBinEdges/dNewBinEdges[0];
				profRebin = profTemp->Rebin(iNumNewBins,"profRebin",dNewBinEdges);
				*/

				iNumBinsX = profTemp->GetNbinsX(); // number of pT (X) bins
				histTemp = (TH1D*) profTemp->ProjectionX()->Clone();

				for(Int_t iBinX(1); iBinX < iNumBinsX+1; iBinX++)
				{
					dValue = histTemp->GetBinContent(iBinX);
					dValue = dValue / dRef;
					histTemp->SetBinContent(iBinX, dValue);
					histTemp->SetBinError(iBinX,0);
				} // end of loop over pt (X) bins

				listTemp[iCent]->Add(histTemp);

				if(dEtaGap < 0.) // no Gap case
				{
					profTemp4 = (TProfile*) listIn->FindObject(Form("f%s_n%d4_gap%02.2g_cent%d_number%d",sSpecies.Data(),iHarmonics,10*dEtaGap,iCent,iSample));
					if(!profTemp4)
					{
						Error("ProcessList",Form("Input POIs %s <<4>> TProfile not found! Name: \"f%s_n%d4_gap%02.2g_cent%d_number%d\"",sSpecies.Data(),sSpecies.Data(),iHarmonics,10*dEtaGap,iCent,iSample));
						return kFALSE;
					}

					histTemp4 = (TH1D*) profTemp4->ProjectionX()->Clone();
					iNumBinsX = histTemp4->GetNbinsX();

					for(Int_t iBinX(1); iBinX < iNumBinsX+1; iBinX++)
					{
						dCum4 = profTemp4->GetBinContent(iBinX);
						dCum2 = profTemp->GetBinContent(iBinX);
						dRefCum2 = profRef2->GetBinContent(iCent+1);
						dRefCum4 = profRef4->GetBinContent(iCent+1);

						dValue = dCum4 - 2 * dCum2 * dRefCum2; // dn{4}
						dRef = dRefCum4 - 2 * TMath::Power(dRefCum2,2);// cn{4}

						if(dRef < 0.)
							histTemp4->SetBinContent(iBinX, -1 * dValue / TMath::Power(-dRef, 0.75) ); // vn{4}
						else
							histTemp4->SetBinContent(iBinX, -9. ); // vn{4}

						histTemp4->SetBinError(iBinX,0);
					} // end of loop over pt (X) bins

					listTemp4[iCent]->Add(histTemp4);


				} // end of if NoGap case


			} // end of loop over centralities
		} // end of loop over rapidity sign (POS / NEG)
	} // end of loop over samples


	// listTemp[fiNumBinsCent] should be full of flow based on samples
	Bool_t bStatusDesample = kFALSE;
	Short_t iNumSamples = fiNumSamples;

	if(dEtaGap >= 0.) // if NoGap case include both POS & Neg POIs (2xiNumSamples)
		iNumSamples = iNumSamples*2;

	for(Short_t iCent(0); iCent < fiNumBinsCent; iCent++)
	{
		bStatusDesample = DesampleList(listTemp[iCent],iNumSamples);
		if(bStatusDesample == kFALSE)
		{
			Error("ProcessList",Form("Desampling: cent %d status %d (FAILED)!",iCent,bStatusDesample));
			return kFALSE;
		}

		listOut->Add(listTemp[iCent]->Last());


		if(dEtaGap < 0.)
		{


			bStatusDesample = DesampleList(listTemp4[iCent],iNumSamples);
			if(bStatusDesample == kFALSE)
			{
				Error("ProcessList",Form("Desampling: cent %d status %d (FAILED)!",iCent,bStatusDesample));
				return kFALSE;
			}

			//listTemp4[iCent]->ls();
			listOut->Add(listTemp4[iCent]->Last());

		}


		delete listTemp[iCent];
		delete listTemp4[iCent];
	}

	return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessFlow::DesampleList(TList* inList, const Short_t iNumSamples)
{
	Debug("DesampleList",Form("Desampling list with %d entries (first: %s)",inList->GetEntries(), inList->First()->GetName()));

	// checking in/out lists
	if(!inList)
	{
		Error("DesampleList","Input list does not exists!");
		return kFALSE;
	}

	if(fbDebug)
		inList->ls();

	if(iNumSamples != inList->GetEntries())
	{
		Error("DesampleList",Form("iNumSamples (%d) is NOT equal to number of entries (%d) in input list.",iNumSamples, inList->GetEntries()));
		return kFALSE;
	}
	TString sName = inList->First()->GetName();
	TString sTitle = inList->First()->GetTitle();

	//printf("%s\n", sName.Data());
	//printf("%s\n", sTitle.Data());
	//sName.Remove(sName.Length()-11); // remove the "_numberX_Z_px" from end of the name
	sName.Append("_desampled");
	//printf("%s\n", sName.Data());
	sTitle.Append("_desampled");
	//sTitle.Remove(sTitle.Length()-8); // remove the "_sampleX_px" from end of the name
	//printf("%s\n", sTitle.Data());

	TH1D* histOut = (TH1D*) inList->First()->Clone(sName.Data());
	histOut->SetStats(0);
	histOut->SetTitle(sTitle.Data());


	const Short_t iNumBinsX = histOut->GetNbinsX();

	TH2D* hValues = new TH2D(Form("Desampling_%s",sName.Data()),"Desampling values spread;values;axisX bin",1000,-0.5,0.5,iNumBinsX,1,iNumBinsX+1);

	TH1D* histTemp = 0x0;
	Double_t dValue = 0, dSum = 0, dMean = 0, dSigma = 0;
	Short_t iWrongSamples = 0;
	for(Short_t i(1); i < iNumBinsX+1; i++)
	{
		dValue = 0;
		dSum = 0;
		dMean = 0;
		dSigma = 0;
		iWrongSamples = 0;

		// estimating mean
		for(Short_t j(0); j < iNumSamples; j++)
		{
			histTemp = (TH1D*) inList->At(j);
			dValue = histTemp->GetBinContent(i);
			if(dValue > -1)
			{
				dSum += dValue;
			}
			else
			{
				iWrongSamples++;
			}

			hValues->Fill(dValue,i);
		}
		if((iNumSamples - iWrongSamples) > 0)
		{
			dMean = dSum / (iNumSamples - iWrongSamples);
		}
		else
		{
			dMean = dSum;
		}

		histOut->SetBinContent(i, dMean);

		// estimating error

		for(Short_t j(0); j < iNumSamples; j++)
		{
			histTemp = (TH1D*) inList->At(j);
			dValue = histTemp->GetBinContent(i);
			if(dValue > -1.)
			{
				dSigma += TMath::Power(dValue - dMean, 2);
			}
		}

		if(iNumSamples - iWrongSamples > 0)
		{
			histOut->SetBinError(i, TMath::Sqrt(dSigma) / (iNumSamples - iWrongSamples));
		}
		else
		{
			histOut->SetBinError(i, TMath::Sqrt(dSigma));
		}
	}

	inList->Add(histOut);

	if(fbSamplingSpread)
	{
		ffOutputFile->cd();
		hValues->Write();
	}

	return kTRUE;
}
//_____________________________________________________________________________
// void ProcessFlow::Run()
// {
// 	return;
// 	ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(), fsInputFileName.Data()),"READ");
// 	if(!ffInputFile->IsOpen())
// 	{
// 		Error("Run","Input file is not open (possibly do not exists). Terminating!");
// 		return;
// 	}

// 	ffOutputFile = new TFile(Form("%s/%s",fsOutputFilePath.Data(), fsOutputFileName.Data()),"RECREATE");
// 	if(!ffOutputFile->IsOpen())
// 	{
// 		Error("Run","Output file is not open (possibly do not exists). Terminating!");
// 		return;
// 	}

// 	// loading files
// 	ffInputFile->cd(Form("flowPID_%s",fsTag.Data()));

// 	//TList* lTracks = (TList*) gDirectory->Get(Form("Tracks_flowPID_%s",fsTag.Data()));
// 	TList* lPID = (TList*) gDirectory->Get(Form("PID_flowPID_%s",fsTag.Data()));

// 	TH1D* hGFKTracksRefTwo[10] = {0x0}; // reference <<2>> = cn2
// 	TH1D* hGFKTracksRefTwoGap00[10] = {0x0}; // reference <<2>> = cn2
// 	TH1D* hGFKTracksRefTwoGap10[10] = {0x0}; // reference <<2>> = cn2
// 	TH1D* hGFKTracksRefTwo3[10] = {0x0}; // reference <<2>> = cn2
// 	TH1D* hGFKTracksRefTwo4[10] = {0x0}; // reference <<2>> = cn2
// 	TH1D* hGFKTracksRefFour[10] = {0x0};	// reference <<4>>
// 	TH1D* hGFKTracksRefFour3[10] = {0x0};	// reference <<4>>
// 	TH1D* hGFKTracksRefFour4[10] = {0x0};	// reference <<4>>

// 	TH1D* hGFKPionsTwo[10][10] = {0x0}; // pt diff <<2'>> (pions) = dn2
// 	TH1D* hGFKPionsFour[10][10] = {0x0}; // pt diff <<4'>> (pions)
// 	TH1D* hGFKKaonsTwo[10][10] = {0x0}; // pt diff <<2'>> (pions) = dn2
// 	TH1D* hGFKKaonsFour[10][10] = {0x0}; // pt diff <<4'>> (pions)
// 	TH1D* hGFKProtonsTwo[10][10] = {0x0}; // pt diff <<2'>> (pions) = dn2
// 	TH1D* hGFKProtonsFour[10][10] = {0x0}; // pt diff <<4'>> (pions)

// 	ffOutputFile->cd();

// 	TProfile* profTemp = 0x0;
// 	TH1D* dummyHist = 0x0; // temporary hist


// 	for(Short_t j(0); j < fiNumSamples; j++)
// 	{
// 		profTemp = (TProfile*) lPID->FindObject(Form("fc22Tracks_gap-1_number%d",j));
// 		hGFKTracksRefTwo[j] = (TH1D*) profTemp->ProjectionX();

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc22Tracks_gap00_number%d",j));
// 		hGFKTracksRefTwoGap00[j] = (TH1D*) profTemp->ProjectionX();

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc22Tracks_gap01_number%d",j));
// 		hGFKTracksRefTwoGap10[j] = (TH1D*) profTemp->ProjectionX();
// 		//hGFKTracksRefTwo[j]->Write(Form("hGFKTrackRefTwo_sample%d",j));

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc32Tracks_gap-1_number%d",j));
// 		hGFKTracksRefTwo3[j] = (TH1D*) profTemp->ProjectionX();

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc42Tracks_gap-1_number%d",j));
// 		hGFKTracksRefTwo4[j] = (TH1D*) profTemp->ProjectionX();

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc24Tracks_gap-1_number%d",j));
// 		hGFKTracksRefFour[j] = (TH1D*) profTemp->ProjectionX();
// 		//hGFKTracksRefFour[j]->Write(Form("hGFKTrackRefFour_sample%d",j));
// 		//hGFKTracksCn4[j] = EstimateCn4(hGFKTracksRefTwo[j], hGFKTracksRefFour[j]);
// 		//hGFKTracksCn4[j]->Write();

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc34Tracks_gap-1_number%d",j));
// 		hGFKTracksRefFour3[j] = (TH1D*) profTemp->ProjectionX();

// 		profTemp = (TProfile*) lPID->FindObject(Form("fc44Tracks_gap-1_number%d",j));
// 		hGFKTracksRefFour4[j] = (TH1D*) profTemp->ProjectionX();

// 		for(Short_t i(0); i < fiNumBinsCent; i++)
// 		{
// 			profTemp = (TProfile*) lPID->FindObject(Form("fd22Pion_gap-1_cent%d_number%d",i,j));
// 			hGFKPionsTwo[i][j] = (TH1D*) profTemp->ProjectionX();
// 			//hGFKPionsTwo[i][j]->Write(Form("hGFKPionsTwo_cent%d_sample%d",i,j));

// 			profTemp = (TProfile*) lPID->FindObject(Form("fd24Pion_gap-1_cent%d_number%d",i,j));
// 			hGFKPionsFour[i][j] = (TH1D*) profTemp->ProjectionX();
// 			//hGFKPionsFour[i][j]->Write(Form("hGFKPionsFour_cent%d_sample%d",i,j));

// 			profTemp = (TProfile*) lPID->FindObject(Form("fd22Kaon_gap-1_cent%d_number%d",i,j));
// 			hGFKKaonsTwo[i][j] = (TH1D*) profTemp->ProjectionX();

// 			profTemp = (TProfile*) lPID->FindObject(Form("fd24Kaon_gap-1_cent%d_number%d",i,j));
// 			hGFKKaonsFour[i][j] = (TH1D*) profTemp->ProjectionX();

// 			profTemp = (TProfile*) lPID->FindObject(Form("fd22Proton_gap-1_cent%d_number%d",i,j));
// 			hGFKProtonsTwo[i][j] = (TH1D*) profTemp->ProjectionX();

// 			profTemp = (TProfile*) lPID->FindObject(Form("fd24Proton_gap-1_cent%d_number%d",i,j));
// 			hGFKProtonsFour[i][j] = (TH1D*) profTemp->ProjectionX();
// 		}

// 	}

// 	// loaded

// 	// estimating 4 particle cummulants dn4, cn4
// 	TH1D* hGFKTracksCn4[10] = {0x0};	// reference cn4
// 	TH1D* hGFKTracksC34[10] = {0x0};	// reference cn4
// 	TH1D* hGFKTracksC44[10] = {0x0};	// reference cn4

// 	TH1D* hGFKPionsDn4[10][10] = {0x0};	// pt diff dn4 (pions)
// 	TH1D* hGFKKaonsDn4[10][10] = {0x0};	// pt diff dn4 (pions)
// 	TH1D* hGFKProtonsDn4[10][10] = {0x0};	// pt diff dn4 (pions)

// 	for(Short_t j(0); j < fiNumSamples; j++)
// 	{
// 		hGFKTracksCn4[j] = EstimateCn4(hGFKTracksRefTwo[j], hGFKTracksRefFour[j]);
// 		hGFKTracksC34[j] = EstimateCn4(hGFKTracksRefTwo3[j], hGFKTracksRefFour3[j]);
// 		hGFKTracksC44[j] = EstimateCn4(hGFKTracksRefTwo4[j], hGFKTracksRefFour4[j]);
// 		//hGFKTracksCn4[j]->Write(Form("GFKTracksCn4_sample%d",j));

// 		for(Short_t i(0); i < fiNumBinsCent; i++)
// 		{
// 			hGFKPionsDn4[i][j] = EstimateDn4(hGFKTracksRefTwo[j]->GetBinContent(i+1),hGFKPionsTwo[i][j],hGFKPionsFour[i][j]);
// 			hGFKKaonsDn4[i][j] = EstimateDn4(hGFKTracksRefTwo[j]->GetBinContent(i+1),hGFKKaonsTwo[i][j],hGFKKaonsFour[i][j]);
// 			hGFKProtonsDn4[i][j] = EstimateDn4(hGFKTracksRefTwo[j]->GetBinContent(i+1),hGFKProtonsTwo[i][j],hGFKProtonsFour[i][j]);
// 			//hGFKPionsDn4[i][j]->Write(Form("GFKPionsDn4_cent%d_sample%d",i,j));
// 		}
// 	}

// 	// estimating v22, v24, v'22, v'24
// 	TH1D* hGFKRefv22[10] = {0x0};
// 	TH1D* hGFKRefv22Gap00[10] = {0x0};
// 	TH1D* hGFKRefv22Gap10[10] = {0x0};
// 	TH1D* hGFKRefv32[10] = {0x0};
// 	TH1D* hGFKRefv42[10] = {0x0};
// 	TH1D* hGFKRefv24[10] = {0x0};
// 	TH1D* hGFKRefv34[10] = {0x0};
// 	TH1D* hGFKRefv44[10] = {0x0};
// 	TH1D* hGFKPionsDiffv22[10][10] = {0x0};
// 	TH1D* hGFKPionsDiffv24[10][10] = {0x0};
// 	TH1D* hGFKKaonsDiffv22[10][10] = {0x0};
// 	TH1D* hGFKKaonsDiffv24[10][10] = {0x0};
// 	TH1D* hGFKProtonsDiffv22[10][10] = {0x0};
// 	TH1D* hGFKProtonsDiffv24[10][10] = {0x0};

// 	Double_t dValue = 0;
// 	for(Short_t j(0); j < fiNumSamples; j++)
// 	{
// 		hGFKRefv22[j] = (TH1D*) hGFKTracksRefTwo[j]->Clone(Form("hGFKRefv22_sample%d",j));
// 		hGFKRefv22[j]->SetTitle(Form("Tracks: v_{2}{2} ref noGap sample %d",j));
// 		hGFKRefv22Gap00[j] = (TH1D*) hGFKTracksRefTwoGap00[j]->Clone(Form("hGFKRefv22_Gap00_sample%d",j));
// 		hGFKRefv22Gap00[j]->SetTitle(Form("Tracks: v_{2}{2} ref Gap)) sample %d",j));
// 		hGFKRefv22Gap10[j] = (TH1D*) hGFKTracksRefTwoGap10[j]->Clone(Form("hGFKRefv22_Gap10_sample%d",j));
// 		hGFKRefv22Gap10[j]->SetTitle(Form("Tracks: v_{2}{2} ref Gap10 sample %d",j));
// 		hGFKRefv32[j] = (TH1D*) hGFKTracksRefTwo3[j]->Clone(Form("hGFKRefv32_sample%d",j));
// 		hGFKRefv32[j]->SetTitle(Form("Tracks: v_{3}{2} ref noGap sample %d",j));
// 		hGFKRefv42[j] = (TH1D*) hGFKTracksRefTwo4[j]->Clone(Form("hGFKRefv42_sample%d",j));
// 		hGFKRefv42[j]->SetTitle(Form("Tracks: v_{4}{2} ref noGap sample %d",j));
// 		hGFKRefv24[j] = (TH1D*) hGFKTracksRefFour[j]->Clone(Form("hGFKRefv24_sample%d",j));
// 		hGFKRefv24[j]->SetTitle(Form("Tracks: v_{2}{4} ref noGap sample %d",j));
// 		hGFKRefv34[j] = (TH1D*) hGFKTracksRefFour3[j]->Clone(Form("hGFKRefv34_sample%d",j));
// 		hGFKRefv34[j]->SetTitle(Form("Tracks: v_{3}{4} ref noGap sample %d",j));
// 		hGFKRefv44[j] = (TH1D*) hGFKTracksRefFour4[j]->Clone(Form("hGFKRefv44_sample%d",j));
// 		hGFKRefv44[j]->SetTitle(Form("Tracks: v_{4}{4} ref noGap sample %d",j));

// 		for(Short_t i(0); i < fiNumBinsCent; i++)
// 		{
// 			dValue = hGFKTracksRefTwo[j]->GetBinContent(i+1);
// 			hGFKRefv22[j]->SetBinContent(i+1,TMath::Sqrt(dValue));
// 			hGFKRefv22[j]->SetBinError(i+1,0);

// 			dValue = hGFKTracksRefTwoGap00[j]->GetBinContent(i+1);
// 			hGFKRefv22Gap00[j]->SetBinContent(i+1,TMath::Sqrt(dValue));
// 			hGFKRefv22Gap00[j]->SetBinError(i+1,0);\

// 			dValue = hGFKTracksRefTwoGap10[j]->GetBinContent(i+1);
// 			hGFKRefv22Gap10[j]->SetBinContent(i+1,TMath::Sqrt(dValue));
// 			hGFKRefv22Gap10[j]->SetBinError(i+1,0);

// 			dValue = hGFKTracksRefTwo3[j]->GetBinContent(i+1);
// 			hGFKRefv32[j]->SetBinContent(i+1,TMath::Sqrt(dValue));
// 			hGFKRefv32[j]->SetBinError(i+1,0);

// 			dValue = hGFKTracksRefTwo4[j]->GetBinContent(i+1);
// 			hGFKRefv42[j]->SetBinContent(i+1,TMath::Sqrt(dValue));
// 			hGFKRefv42[j]->SetBinError(i+1,0);

// 			dValue = hGFKTracksCn4[j]->GetBinContent(i+1);
// 			hGFKRefv24[j]->SetBinContent(i+1,TMath::Power(-dValue,0.25));
// 			hGFKRefv24[j]->SetBinError(i+1,0);

// 			dValue = hGFKTracksC34[j]->GetBinContent(i+1);
// 			if(dValue < 0)
// 			{
// 				hGFKRefv34[j]->SetBinContent(i+1,TMath::Power(-dValue,0.25));
// 			}
// 			else
// 			{
// 				hGFKRefv34[j]->SetBinContent(i+1,10);
// 			}
// 			hGFKRefv34[j]->SetBinError(i+1,0);

// 			dValue = hGFKTracksC44[j]->GetBinContent(i+1);
// 			if(dValue < 0)
// 			{
// 				hGFKRefv44[j]->SetBinContent(i+1,TMath::Power(-dValue,0.25));
// 			}
// 			else
// 			{
// 				hGFKRefv44[j]->SetBinContent(i+1,10);
// 			}
// 			hGFKRefv44[j]->SetBinError(i+1,0);

// 			hGFKPionsDiffv22[i][j] = (TH1D*) hGFKPionsTwo[i][j]->Clone(Form("hGFKPionsDiffv22_sample%d",j));
// 			hGFKPionsDiffv22[i][j]->SetTitle(Form("#pi: v'_{2}{2} diff cent %d sample %d",i,j));
// 			hGFKPionsDiffv24[i][j] = (TH1D*) hGFKPionsDn4[i][j]->Clone(Form("hGFKPionsDiffv24_sample%d",j));
// 			hGFKPionsDiffv24[i][j]->SetTitle(Form("#pi: v'_{2}{4} diff cent %d sample %d",i,j));

// 			hGFKKaonsDiffv22[i][j] = (TH1D*) hGFKKaonsTwo[i][j]->Clone(Form("hGFKKaonsDiffv22_sample%d",j));
// 			hGFKKaonsDiffv22[i][j]->SetTitle(Form("K: v'_{2}{2} diff cent %d sample %d",i,j));
// 			hGFKKaonsDiffv24[i][j] = (TH1D*) hGFKKaonsDn4[i][j]->Clone(Form("hGFKKaonsDiffv24_sample%d",j));
// 			hGFKKaonsDiffv24[i][j]->SetTitle(Form("K: v'_{2}{4} diff cent %d sample %d",i,j));

// 			hGFKProtonsDiffv22[i][j] = (TH1D*) hGFKProtonsTwo[i][j]->Clone(Form("hGFKProtonsDiffv22_sample%d",j));
// 			hGFKProtonsDiffv22[i][j]->SetTitle(Form("p: v'_{2}{2} diff cent %d sample %d",i,j));
// 			hGFKProtonsDiffv24[i][j] = (TH1D*) hGFKProtonsDn4[i][j]->Clone(Form("hGFKProtonsDiffv24_sample%d",j));
// 			hGFKProtonsDiffv24[i][j]->SetTitle(Form("p: v'_{2}{4} diff cent %d sample %d",i,j));


// 			for(Int_t k(1); k < (hGFKPionsDiffv22[0][0]->GetNbinsX() + 1); k++)
// 			{
// 				dValue = hGFKPionsTwo[i][j]->GetBinContent(k);
// 				dValue = dValue / hGFKRefv22[j]->GetBinContent(i+1);
// 				hGFKPionsDiffv22[i][j]->SetBinContent(k,dValue);
// 				hGFKPionsDiffv22[i][j]->SetBinError(k,0);

// 				dValue = hGFKPionsDn4[i][j]->GetBinContent(k);
// 				dValue = dValue / TMath::Power(hGFKRefv24[j]->GetBinContent(i+1),3);
// 				hGFKPionsDiffv24[i][j]->SetBinContent(k,-dValue);
// 				hGFKPionsDiffv24[i][j]->SetBinError(k,0);
// 			}

// 			for(Int_t k(1); k < (hGFKKaonsDiffv22[0][0]->GetNbinsX() + 1); k++)
// 			{
// 				dValue = hGFKKaonsTwo[i][j]->GetBinContent(k);
// 				dValue = dValue / hGFKRefv22[j]->GetBinContent(i+1);
// 				hGFKKaonsDiffv22[i][j]->SetBinContent(k,dValue);
// 				hGFKKaonsDiffv22[i][j]->SetBinError(k,0);

// 				dValue = hGFKKaonsDn4[i][j]->GetBinContent(k);
// 				dValue = dValue / TMath::Power(hGFKRefv24[j]->GetBinContent(i+1),3);
// 				hGFKKaonsDiffv24[i][j]->SetBinContent(k,-dValue);
// 				hGFKKaonsDiffv24[i][j]->SetBinError(k,0);
// 			}

// 			for(Int_t k(1); k < (hGFKProtonsDiffv22[0][0]->GetNbinsX() + 1); k++)
// 			{
// 				dValue = hGFKProtonsTwo[i][j]->GetBinContent(k);
// 				dValue = dValue / hGFKRefv22[j]->GetBinContent(i+1);
// 				hGFKProtonsDiffv22[i][j]->SetBinContent(k,dValue);
// 				hGFKProtonsDiffv22[i][j]->SetBinError(k,0);

// 				dValue = hGFKProtonsDn4[i][j]->GetBinContent(k);
// 				dValue = dValue / TMath::Power(hGFKRefv24[j]->GetBinContent(i+1),3);
// 				hGFKProtonsDiffv24[i][j]->SetBinContent(k,-dValue);
// 				hGFKProtonsDiffv24[i][j]->SetBinError(k,0);
// 			}

// 		  //hGFKPionsDiffv22[i][j]->Write(Form("hGFKPionsDiffv22_cent%d_sample%d",i,j));
// 		  //hGFKPionsDiffv24[i][j]->Write(Form("hGFKPionsDiffv24_cent%d_sample%d",i,j));

// 		}

// 		//hGFKRefv22[j]->Write(Form("hGFKRefv22_sample%d",j));
// 		//hGFKRefv24[j]->Write(Form("hGFKRefv24_sample%d",j));
// 	}

// 	// up to this point: validated & works well

// 	// desampling
// 	TList* lRefv22 = new TList();
// 	TList* lRefv22Gap00 = new TList();
// 	TList* lRefv22Gap10 = new TList();
// 	TList* lRefv32 = new TList();
// 	TList* lRefv42 = new TList();
// 	TList* lRefv24 = new TList();
// 	TList* lRefv34 = new TList();
// 	TList* lRefv44 = new TList();
// 	TList* lPionsDiffv22 = new TList();
// 	TList* lPionsDiffv24 = new TList();
// 	TList* lKaonsDiffv22 = new TList();
// 	TList* lKaonsDiffv24 = new TList();
// 	TList* lProtonsDiffv22 = new TList();
// 	TList* lProtonsDiffv24 = new TList();

// 	// output list by flow harmonics & particle species
// 	TList* lFinalPionsDiffv22 = new TList();
// 	TList* lFinalKaonsDiffv22 = new TList();
// 	TList* lFinalProtonsDiffv22 = new TList();
// 	TList* lFinalPionsDiffv24 = new TList();
// 	TList* lFinalKaonsDiffv24 = new TList();
// 	TList* lFinalProtonsDiffv24 = new TList();

// 	ffOutputFile->cd();
// 	// filling the lists
// 	// reference
// 	for(Short_t j(0); j < fiNumSamples; j++)
// 	{
// 		lRefv22->Add(hGFKRefv22[j]);
// 		lRefv22Gap00->Add(hGFKRefv22Gap00[j]);
// 		lRefv22Gap10->Add(hGFKRefv22Gap10[j]);
// 		lRefv32->Add(hGFKRefv32[j]);
// 		lRefv42->Add(hGFKRefv42[j]);
// 		lRefv24->Add(hGFKRefv24[j]);
// 		lRefv34->Add(hGFKRefv34[j]);
// 		lRefv44->Add(hGFKRefv44[j]);
// 	}

// 	DesampleList(lRefv22,fiNumSamples);
// 	DesampleList(lRefv22Gap00,fiNumSamples);
// 	DesampleList(lRefv22Gap10,fiNumSamples);
// 	DesampleList(lRefv32,fiNumSamples);
// 	DesampleList(lRefv42,fiNumSamples);
// 	DesampleList(lRefv24,fiNumSamples);
// 	DesampleList(lRefv34,fiNumSamples);
// 	DesampleList(lRefv44,fiNumSamples);

// 	(lRefv22->Last())->Write("Refv22_Tracks");
// 	(lRefv22Gap00->Last())->Write("Refv22_Gap00_Tracks");
// 	(lRefv22Gap10->Last())->Write("Refv22_Gap10_Tracks");
// 	(lRefv32->Last())->Write("Refv32_Tracks");
// 	(lRefv42->Last())->Write("Refv42_Tracks");
// 	(lRefv24->Last())->Write("Refv24_Tracks");
// 	(lRefv34->Last())->Write("Refv34_Tracks");
// 	(lRefv44->Last())->Write("Refv44_Tracks");

// 	// pt diff
// 	for(Short_t i(0); i < fiNumBinsCent; i++)
// 	{
// 		lPionsDiffv22->Clear();
// 		lPionsDiffv24->Clear();
// 		lKaonsDiffv22->Clear();
// 		lKaonsDiffv24->Clear();
// 		lProtonsDiffv22->Clear();
// 		lProtonsDiffv24->Clear();

// 		for(Short_t j(0); j < fiNumSamples; j++)
// 		{
// 			lPionsDiffv22->Add(hGFKPionsDiffv22[i][j]);
// 			lPionsDiffv24->Add(hGFKPionsDiffv24[i][j]);
// 			lKaonsDiffv22->Add(hGFKKaonsDiffv22[i][j]);
// 			lKaonsDiffv24->Add(hGFKKaonsDiffv24[i][j]);
// 			lProtonsDiffv22->Add(hGFKProtonsDiffv22[i][j]);
// 			lProtonsDiffv24->Add(hGFKProtonsDiffv24[i][j]);
// 		}

// 		DesampleList(lPionsDiffv22,fiNumSamples);
// 		DesampleList(lPionsDiffv24,fiNumSamples);
// 		DesampleList(lKaonsDiffv22,fiNumSamples);
// 		DesampleList(lKaonsDiffv24,fiNumSamples);
// 		DesampleList(lProtonsDiffv22,fiNumSamples);
// 		DesampleList(lProtonsDiffv24,fiNumSamples);

// 		//lPionsDiffv22->Write(Form("lPionsDiffv22_cent%d",i),TObject::kSingleKey);
// 		//lPionsDiffv24->Write(Form("lPionsDiffv24_cent%d",i),TObject::kSingleKey);
// 		//lKaonsDiffv22->Write(Form("lKaonsDiffv22_cent%d",i),TObject::kSingleKey);
// 		//lKaonsDiffv24->Write(Form("lKaonsDiffv24_cent%d",i),TObject::kSingleKey);
// 		//lProtonsDiffv22->Write(Form("lProtonsDiffv22_cent%d",i),TObject::kSingleKey);
// 		//lProtonsDiffv24->Write(Form("lProtonsDiffv24_cent%d",i),TObject::kSingleKey);
// 		/*
// 		(lPionsDiffv22->Last())->Write(Form("hPionsDiffv22_cent%d",i));
// 		(lPionsDiffv24->Last())->Write(Form("hPionsDiffv24_cent%d",i));
// 		(lKaonsDiffv22->Last())->Write(Form("hKaonsDiffv22_cent%d",i));
// 		(lKaonsDiffv24->Last())->Write(Form("hKaonsDiffv24_cent%d",i));
// 		(lProtonsDiffv22->Last())->Write(Form("hProtonsDiffv22_cent%d",i));
// 		(lProtonsDiffv24->Last())->Write(Form("hProtonsDiffv24_cent%d",i));
// 		*/


// 		// Writing output TLists & setting names for final output
// 		dummyHist = (TH1D*) lPionsDiffv22->Last();
// 		dummyHist->SetName(Form("Diffv22_Pions_cent%d",i));
// 		lFinalPionsDiffv22->Add(dummyHist);
// 		dummyHist = (TH1D*)lKaonsDiffv22->Last();
// 		dummyHist->SetName(Form("Diffv22_Kaons_cent%d",i));
// 		lFinalKaonsDiffv22->Add(dummyHist);
// 		dummyHist = (TH1D*)lProtonsDiffv22->Last();
// 		dummyHist->SetName(Form("Diffv22_Protons_cent%d",i));
// 		lFinalProtonsDiffv22->Add(dummyHist);

// 		dummyHist = (TH1D*) lPionsDiffv24->Last();
// 		dummyHist->SetName(Form("Diffv24_Pions_cent%d",i));
// 		lFinalPionsDiffv24->Add(dummyHist);
// 		dummyHist = (TH1D*)lKaonsDiffv24->Last();
// 		dummyHist->SetName(Form("Diffv24_Kaons_cent%d",i));
// 		lFinalKaonsDiffv24->Add(dummyHist);
// 		dummyHist = (TH1D*)lProtonsDiffv24->Last();
// 		dummyHist->SetName(Form("Diffv24_Protons_cent%d",i));
// 		lFinalProtonsDiffv24->Add(dummyHist);

// 	}

// 	delete lRefv22;
// 	delete lRefv22Gap00;
// 	delete lRefv22Gap10;
// 	delete lRefv32;
// 	delete lRefv42;
// 	delete lRefv24;
// 	delete lRefv34;
// 	delete lRefv44;
// 	delete lPionsDiffv22;
// 	delete lPionsDiffv24;
// 	delete lKaonsDiffv22;
// 	delete lKaonsDiffv24;
// 	delete lProtonsDiffv22;
// 	delete lProtonsDiffv24;

// 	// writing to output file
// 	ffOutputFile->cd();
// 	lFinalPionsDiffv22->Write("Diffv22_Pions",TObject::kSingleKey);
// 	lFinalKaonsDiffv22->Write("Diffv22_Kaons",TObject::kSingleKey);
// 	lFinalProtonsDiffv22->Write("Diffv22_Protons",TObject::kSingleKey);

// 	lFinalPionsDiffv24->Write("Diffv24_Pions",TObject::kSingleKey);
// 	lFinalKaonsDiffv24->Write("Diffv24_Kaons",TObject::kSingleKey);
// 	lFinalProtonsDiffv24->Write("Diffv24_Protons",TObject::kSingleKey);

// 	//lDiffv22->Write("Diffv22",TObject::kSingleKey);
// 	//lDiffv24->Write("Diffv24",TObject::kSingleKey);



// 	Terminate();

// 	return;
// }
// //_____________________________________________________________________________
// TH1D* ProcessFlow::EstimateCn2(const TH1D* hCum2)
// {
// 	TH1D* hOut = (TH1D*) hCum2->Clone();

// 	const Short_t iNumBinsPt = hOut->GetNbinsX();

// 	Double_t dValue = 0;
// 	for(Int_t i(1); i < iNumBinsPt+1; i++)
// 	{
// 		dValue = hCum2->GetBinContent(i);
// 		hOut->SetBinContent(i,dValue);
// 		hOut->SetBinError(i,0);
// 	}

// 	return hOut;
// }
// //_____________________________________________________________________________
// TH1D* ProcessFlow::EstimateCn4(const TH1D* hCum2,const TH1D* hCum4)
// {
// 	if(!hCum2 || !hCum4)
// 	{
// 		Error("EstimateRefFour","hCum2 || hCum4 does NOT exists!");
// 		return 0x0;
// 	}

// 	TH1D* hOut = (TH1D*) hCum4->Clone();

// 	const Short_t iNumBinsPt = hOut->GetNbinsX();

// 	Double_t dValue2 = 0, dValue4 = 0, dOutValue = 0;
// 	for(Int_t i(1); i < iNumBinsPt+1; i++)
// 	{
// 		dValue2 = hCum2->GetBinContent(i);
// 		dValue4 = hCum4->GetBinContent(i);
// 		dOutValue = dValue4 - 2*TMath::Power(dValue2,2);

// 		hOut->SetBinContent(i,dOutValue);
// 		hOut->SetBinError(i,0);
// 	}

// 	return hOut;
// }
// //_____________________________________________________________________________
// TH1D* ProcessFlow::EstimateDn4(const Double_t dRef2, const TH1D* hDiff2, const TH1D* hDiff4)
// {
// 	TH1D* hOut = (TH1D*) hDiff4->Clone();

// 	const Short_t iNumBinsPt = hOut->GetNbinsX();

// 	Double_t dValue2 = 0, dValue4 = 0, dOutValue = 0;
// 	for(Int_t i(1); i < iNumBinsPt+1; i++)
// 	{
// 		dValue2 = hDiff2->GetBinContent(i);
// 		dValue4 = hDiff4->GetBinContent(i);
// 		dOutValue = dValue4 - 2 * dValue2 * dRef2;

// 		hOut->SetBinContent(i,dOutValue);
// 		hOut->SetBinError(i,0);
// 	}

// 	return hOut;
// }
// //_____________________________________________________________________________
// void ProcessFlow::DesampleList(TList* inList, const Short_t iNumSamples)
// {
// 	// init

// 	TString sName = inList->First()->GetName();
// 	TString sTitle = inList->First()->GetTitle();

// 	//printf("%s\n", sName.Data());
// 	sName.Remove(sName.Length()-8); // remove the "_sampleX_px" from end of the name
// 	//printf("%s\n", sName.Data());
// 	sTitle.Remove(sTitle.Length()-8); // remove the "_sampleX_px" from end of the name

// 	TH1D* histOut = (TH1D*) inList->First()->Clone(sName.Data());
// 	histOut->SetTitle(sTitle.Data());

// 	const Short_t iNumBinsX = histOut->GetNbinsX();

// 	TH1D* histTemp = 0x0;
// 	Double_t dValue = 0, dMean = 0, dSigma = 0;
// 	for(Short_t i(1); i < iNumBinsX+1; i++)
// 	{
// 		dValue = 0;
// 		dMean = 0;
// 		dSigma = 0;

// 		// estimating mean
// 		for(Short_t j(0); j < iNumSamples; j++)
// 		{
// 			histTemp = (TH1D*) inList->At(j);
// 			dValue += histTemp->GetBinContent(i);
// 		}

// 		dMean = dValue / iNumSamples;

// 		// estimating error
// 		for(Short_t j(0); j < iNumSamples; j++)
// 		{
// 			histTemp = (TH1D*) inList->At(j);
// 			dSigma += TMath::Power(histTemp->GetBinContent(i) - dMean, 2);
// 		}

// 		dSigma = TMath::Sqrt(dSigma / (iNumSamples-1) );

// 		histOut->SetBinContent(i, dMean);
// 		histOut->SetBinError(i, dSigma / TMath::Sqrt(iNumSamples));
// 	}

// 	inList->Add(histOut);


// 	return;
// }
//_____________________________________________________________________________
void ProcessFlow::Fatal(TString sMethod, TString sMsg)
{
	printf("Fatal-ProcessFlow::%s: %s Terminating!\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessFlow::Error(TString sMethod, TString sMsg)
{
	printf("Error-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessFlow::Info(TString sMethod, TString sMsg)
{
	printf("Info-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessFlow::Warning(TString sMethod, TString sMsg)
{
	printf("Warning-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessFlow::Debug(TString sMethod, TString sMsg)
{
	if(fbDebug)
		printf("Debug-ProcessFlow::%s: %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
