#include "AliAnalysisTaskChargedFlow.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskChargedFlow)
//___________________________________________________________________________
AliAnalysisTaskChargedFlow::AliAnalysisTaskChargedFlow():
  AliAnalysisTaskSE(),
  fAOD(0),
  fAnalysisType("AOD"),
  fAnalysisMC(kFALSE),
  fVtxCut(10.0),  
  fFilterbit(128),
  fEtaCut(0.8),
  fNoClus(70),
  fEffsys(1),
  fULS(1),
  fMinPt(0.2),
  fMaxPt(5.0),
  fNHarm(2),
  fSample(kFALSE),
  fLHC10h(kTRUE),
  fPileUp(kTRUE),
  fNPtBins(14),
  fCent(0),
  fListOfObjects(0),
	fDifferential(1),
	fVtx(0),
	fVtxBeforeCuts(0),
	fVtxAfterCuts(0),
	fCentralityDis(0),
	fPhiDis(0),
	fEtaDis(0),
	fMultCorBeforeCuts(0),
	fMultCorAfterCuts(0),
	fCentSPDvsV0M(0),
	fMultvsCentr(0),
	/*
	fChc22NoGap(0),
	fChc22Gap00(0),
	fChc22Gap04(0),
	fChc22Gap08(0),
	fChc22Gap10(0),
	fChc32NoGap(0),
	fChc32Gap00(0),
	fChc32Gap04(0),
	fChc32Gap08(0),
	fChc32Gap10(0),
	fChc42NoGap(0),
	fChc42Gap00(0),
	fChc42Gap04(0),
	fChc42Gap08(0),
	fChc42Gap10(0)
	*/
	fPhiDisGapM(0),
	fPhiDisGapP(0),
	fPhiDisGap0M(0),
	fPhiDisGap0P(0),
	fPhiDisGap4M(0),
	fPhiDisGap4P(0),
	fPhiDisGap8M(0),
	fPhiDisGap8P(0),
	fPhiDispT2(0),
	fPhiDispT3(0),
	fPhiDispT5(0),
	fPhiDisVtx10(0),
	fPhiDisVtx9(0),
	fPhiDisVtx8(0),
	fPhiDisVtx7(0),
	fWeight(0)

{
 	for(int cent=0; cent<10; cent++)
	{
		fPhiDisCent[cent] = 0;
	}
  
	for (Int_t j = 0; j < fNPtBins+1; j++)
  {
    fPtBins[j] = 0;
		fPhiPt[j] = 0;
	}
 /*
 		// Qvectors init
		for(Int_t i = 0; i < 20; i++)
		{
			Qvector[i] = 0;
			QvectorM[i] = 0;
			QvectorP[i] = 0;
			Qvector0M[i] = 0;
			Qvector0P[i] = 0;
			Qvector4M[i] = 0;
			Qvector4P[i] = 0;
			Qvector8M[i] = 0;
			Qvector8P[i] = 0;
			pvector[i] = 0;
			pvectorM[i] = 0;
			pvectorP[i] = 0;
			qvector[i] = 0;
			pvector0M[i] = 0;
			pvector0P[i] = 0;
			pvector4M[i] = 0;
			pvector4P[i] = 0;
			pvector8M[i] = 0;
			pvector8P[i] = 0;
			
			for (Int_t j = 0; j < 20; j++)
			{
				Qvector[i][j] = 0;
				QvectorM[i][j] = 0;
				QvectorP[i][j] = 0;
				Qvector0M[i][j] = 0;
				Qvector0P[i][j] = 0;
				Qvector4M[i][j] = 0;
				Qvector4P[i][j] = 0;
				Qvector8M[i][j] = 0;
				Qvector8P[i][j] = 0;
				pvector[i][j] = 0;
				pvectorM[i][j] = 0;
				pvectorP[i][j] = 0;
				qvector[i][j] = 0;
				pvector0M[i][j] = 0;
				pvector0P[i][j] = 0;
				pvector4M[i][j] = 0;
				pvector4P[i][j] = 0;
				pvector8M[i][j] = 0;
				pvector8P[i][j] = 0;
			}
		}
 
*/
	for(int k=0; k<11; k++)
	{			

		fsc4242Re[k] = 0;
		fsc3232Re[k] = 0;
		fsc4343Re[k] = 0;
		//	 fsc5252Re[k] = 0;
		//	 fsc5353Re[k] = 0;

		for(Int_t h = 0; h < 5; h++)	
		{		
			fcn2Re[h][k] = 0;
			fcn2GapRe[h][k] = 0;
			fcn2Gap0Re[h][k] = 0;
			fcn2Gap4Re[h][k] = 0;
			fcn2Gap8Re[h][k] = 0;
			fcn4Re[h][k] = 0;
			fcn6Re[h][k] = 0;
			fcn8Re[h][k] = 0;

			for(int c=0; c<10; c++)
			{
				fdn2Re[h][c][k] = 0;
				fdn2GapQMRe[h][c][k] = 0;
				fdn2GapQPRe[h][c][k] = 0;
				fdn2Gap0QMRe[h][c][k] = 0;
				fdn2Gap0QPRe[h][c][k] = 0;
				fdn2Gap4QMRe[h][c][k] = 0;
				fdn2Gap4QPRe[h][c][k] = 0;
				fdn2Gap8QMRe[h][c][k] = 0;
				fdn2Gap8QPRe[h][c][k] = 0;
				fdn4Re[h][c][k] = 0;

				fcn2PtGap0Re[h][c][k] = 0;
				fcn2PtGap4Re[h][c][k] = 0;
				fcn2PtGap8Re[h][c][k] = 0;
				fcn2PtGap10Re[h][c][k] = 0;
				fcn4PtRe[h][c][k] = 0;
			}
 		}     
  }
    
}

//______________________________________________________________________________
AliAnalysisTaskChargedFlow::AliAnalysisTaskChargedFlow(const char *name):
  AliAnalysisTaskSE(name),
  fAOD(0),
  fAnalysisType("AOD"),
  fAnalysisMC(kFALSE),
  fVtxCut(10.0),  
  fFilterbit(128),
  fEtaCut(0.8),
  fNoClus(70),
  fEffsys(1),
  fULS(1),
  fMinPt(0.2),
  fMaxPt(6.0),
  fNHarm(2),
  fSample(kFALSE),
  fLHC10h(kTRUE),
  fPileUp(kTRUE),
  fNPtBins(14),
  fCent(0),
  fListOfObjects(0),
	fDifferential(1),
	fVtx(0),
	fVtxBeforeCuts(0),
	fVtxAfterCuts(0),
	fCentralityDis(0),
	fPhiDis(0),
	fEtaDis(0),
	fMultCorBeforeCuts(0),
	fMultCorAfterCuts(0),
	fCentSPDvsV0M(0),
	fMultvsCentr(0),
	/*
	fChc22NoGap(0),
	fChc22Gap00(0),
	fChc22Gap04(0),
	fChc22Gap08(0),
	fChc22Gap10(0),
	fChc32NoGap(0),
	fChc32Gap00(0),
	fChc32Gap04(0),
	fChc32Gap08(0),
	fChc32Gap10(0),
	fChc42NoGap(0),
	fChc42Gap00(0),
	fChc42Gap04(0),
	fChc42Gap08(0),
	fChc42Gap10(0)
	*/
	fPhiDisGapM(0),
	fPhiDisGapP(0),
	fPhiDisGap0M(0),
	fPhiDisGap0P(0),
	fPhiDisGap4M(0),
	fPhiDisGap4P(0),
	fPhiDisGap8M(0),
	fPhiDisGap8P(0),
	fPhiDispT2(0),
	fPhiDispT3(0),
	fPhiDispT5(0),
	fPhiDisVtx10(0),
	fPhiDisVtx9(0),
	fPhiDisVtx8(0),
	fPhiDisVtx7(0),
	fWeight(0)
{		
		for(int cent=0; cent<10; cent++)
		{
			fPhiDisCent[cent] = 0;
		}
 
    for (Int_t j = 0; j < fNPtBins+1; j++)
    {    
			fPtBins[j] = 0;
			fPhiPt[j] = 0;
		}



		/* 
		// Qvectors init
		for(Int_t i = 0; i < 20; i++)
		{
			Qvector[i] = 0;
			QvectorM[i] = 0;
			QvectorP[i] = 0;
			Qvector0M[i] = 0;
			Qvector0P[i] = 0;
			Qvector4M[i] = 0;
			Qvector4P[i] = 0;
			Qvector8M[i] = 0;
			Qvector8P[i] = 0;
			pvector[i] = 0;
			pvectorM[i] = 0;
			pvectorP[i] = 0;
			qvector[i] = 0;
			pvector0M[i] = 0;
			pvector0P[i] = 0;
			pvector4M[i] = 0;
			pvector4P[i] = 0;
			pvector8M[i] = 0;
			pvector8P[i] = 0;
			
			for (Int_t j = 0; j < 20; j++)
			{
				Qvector[i][j] = 0;
				QvectorM[i][j] = 0;
				QvectorP[i][j] = 0;
				Qvector0M[i][j] = 0;
				Qvector0P[i][j] = 0;
				Qvector4M[i][j] = 0;
				Qvector4P[i][j] = 0;
				Qvector8M[i][j] = 0;
				Qvector8P[i][j] = 0;
				pvector[i][j] = 0;
				pvectorM[i][j] = 0;
				pvectorP[i][j] = 0;
				qvector[i][j] = 0;
				pvector0M[i][j] = 0;
				pvector0P[i][j] = 0;
				pvector4M[i][j] = 0;
				pvector4P[i][j] = 0;
				pvector8M[i][j] = 0;
				pvector8P[i][j] = 0;
			}
		}
		*/
		

  
	  for (Int_t k = 1; k<11; k++)
		{

			fsc4242Re[k] = 0;
			fsc3232Re[k] = 0;
	 		fsc4343Re[k] = 0;
//	 		fsc5252Re[k] = 0;
//	 		fsc5353Re[k] = 0;

			for(Int_t h = 0; h < 5; h++)
			{
      	fcn2Re[h][k] = 0;
      	fcn2GapRe[h][k] = 0;
      	fcn2Gap0Re[h][k] = 0;
      	fcn2Gap4Re[h][k] = 0;
      	fcn2Gap8Re[h][k] = 0;
      	fcn4Re[h][k] = 0;
      	fcn6Re[h][k] = 0;
      	fcn8Re[h][k] = 0;
			
				for(int c=0; c<10; c++)
				{
		 			fdn2Re[h][c][k] = 0;
					fdn2GapQMRe[h][c][k] = 0;
					fdn2GapQPRe[h][c][k] = 0;
					fdn2Gap0QMRe[h][c][k] = 0;
					fdn2Gap0QPRe[h][c][k] = 0;
					fdn2Gap4QMRe[h][c][k] = 0;
					fdn2Gap4QPRe[h][c][k] = 0;
					fdn2Gap8QMRe[h][c][k] = 0;
					fdn2Gap8QPRe[h][c][k] = 0;
		 			fdn4Re[h][c][k] = 0;

					fcn2PtGap0Re[h][c][k] = 0;
					fcn2PtGap4Re[h][c][k] = 0;
					fcn2PtGap8Re[h][c][k] = 0;
					fcn2PtGap10Re[h][c][k] = 0;
					fcn4PtRe[h][c][k] = 0;

				}
			}
      
 		} 
		
  
	// Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
	
}

//_____________________________________________________________________________
AliAnalysisTaskChargedFlow::~AliAnalysisTaskChargedFlow()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects) 
    delete fListOfObjects;

}

//______________________________________________________________________________
void AliAnalysisTaskChargedFlow::UserCreateOutputObjects()
{ 
	
  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  // Histograms
  fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
  fListOfObjects->Add(fVtx);
  
  fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);
    
    fPhiDis = new TH1F("fPhiDis", "phi distribution; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDis);
    
		fPhiDisGapM = new TH1F("fPhiDisGapM", "phi distribution eta < -0.5; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGapM);
    
		fPhiDisGapP = new TH1F("fPhiDisGapP", "phi distribution eta > 0.5; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGapP);
    
		fPhiDisGap0M = new TH1F("fPhiDisGap0M", "phi distribution eta < -0; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGap0M);
    
		fPhiDisGap0P = new TH1F("fPhiDisGap0P", "phi distribution eta > 0; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGap0P);
    
		fPhiDisGap4M = new TH1F("fPhiDisGap4M", "phi distribution eta < -0.2; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGap4M);
    
		fPhiDisGap4P = new TH1F("fPhiDisGap4P", "phi distribution eta > 0.2; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGap4P);
    
		fPhiDisGap8M = new TH1F("fPhiDisGap8M", "phi distribution eta < -0.4; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGap8M);
    
		fPhiDisGap8P = new TH1F("fPhiDisGap8P", "phi distribution eta > 0.4; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDisGap8P);
    
		fPhiDispT2 = new TH1F("fPhiDispT1", "phi distribution 0.2 < pT < 1; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDispT2);
    
		fPhiDispT3 = new TH1F("fPhiDispT2", "phi distribution 1<pT<2; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDispT3);
    
		fPhiDispT5 = new TH1F("fPhiDispT5", "phi distribution 2<pT<5; #phi; Counts", 100, 0, 6.29);
    fListOfObjects->Add(fPhiDispT5);

		fPhiDisVtx10 = new TH1F("fPhiDisVtx10", "phi distribution |Vz|<10; #phi", 100, 0, 6.29);
		fListOfObjects->Add(fPhiDisVtx10);
    
		fPhiDisVtx9 = new TH1F("fPhiDisVtx9", "phi distribution |Vz|<9; #phi", 100, 0, 6.29);
		fListOfObjects->Add(fPhiDisVtx9);
		
		fPhiDisVtx8 = new TH1F("fPhiDisVtx8", "phi distribution |Vz|<8; #phi", 100, 0, 6.29);
		fListOfObjects->Add(fPhiDisVtx8);
		
		fPhiDisVtx7 = new TH1F("fPhiDisVtx7", "phi distribution |Vz|<7; #phi", 100, 0, 6.29);
		fListOfObjects->Add(fPhiDisVtx7);
		

		for(int cent=0; cent<10; cent++)
		{
			fPhiDisCent[cent] = new TH1F(Form("fPhiDis_cent%d", cent), Form("phi distribution cent%d; #phi; Counts", cent), 100, 0, 6.29);
    	fListOfObjects->Add(fPhiDisCent[cent]);
  	}  

		for(int p=0; p<fNPtBins; p++)
		{
			fPhiPt[p] = new TH1F(Form("fPhiDis_pt%d",p), "phi pt; #phi; Counts", 100, 0, 6.29);
		}
	
    fEtaDis = new TH1F("fEtaDis", "eta distribution; #eta; Counts", 100, -1., 1.);
    fListOfObjects->Add(fEtaDis);
    
  fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
  fListOfObjects->Add(fCentralityDis);
    
  fCentSPDvsV0M = new TH2F("fCentSPDvsV0M", "V0M-cent vs SPD-cent; V0M; SPD-cent", 100, 0, 100, 100, 0, 100);
  fListOfObjects->Add(fCentSPDvsV0M);
    
  fMultCorBeforeCuts = new TH2F("fMultCorBeforeCuts", "TPC vs Global multiplicity (Before cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
  fListOfObjects->Add(fMultCorBeforeCuts);
  
  fMultCorAfterCuts = new TH2F("fMultCorAfterCuts", "TPC vs Global multiplicity (After cuts); Global multiplicity; TPC multiplicity", 100, 0, 3000, 100, 0, 3000);
  fListOfObjects->Add(fMultCorAfterCuts);

  fMultvsCentr = new TH2F("fMultvsCentr", "Multiplicity vs centrality; centrality; Multiplicity", 22, -0.5, 21.5, 100, 0, 3000);
  fListOfObjects->Add(fMultvsCentr);
  
	fWeight = new TH1F("fWeight","corrected Phi dist. by weight; phi", 100, 0, 6.29);
	fListOfObjects->Add(fWeight);

    for (Int_t j=1; j < 11; j++) {

			fsc4242Re[j] = new TProfile(Form("fsc4242Re_number%d", j), "SC(4,2,-4,-2) Re; centrality", 10, -0.5, 9.5);
			fsc4242Re[j]->Sumw2();
			fListOfObjects->Add(fsc4242Re[j]);

			fsc3232Re[j] = new TProfile(Form("fsc3232Re_number%d", j), "SC(3,2,-3,-2) Re; centrality", 10, -0.5, 9.5);
			fsc3232Re[j]->Sumw2();
			fListOfObjects->Add(fsc3232Re[j]);
			
			fsc4343Re[j] = new TProfile(Form("fsc4343Re_number%d", j), "SC(4,3,-4,-3) Re; centrality", 10, -0.5, 9.5);
			fsc4343Re[j]->Sumw2();
			fListOfObjects->Add(fsc4343Re[j]);
		/*	
			fsc5252Re[j] = new TProfile(Form("fsc5252Re_number%d", j), "SC(5,2,-5,-2) Re; centrality", 10, -0.5, 9.5);
			fsc5252Re[j]->Sumw2();
			fListOfObjects->Add(fsc5252Re[j]);

			fsc5353Re[j] = new TProfile(Form("fsc5353Re_number%d", j), "SC(5,3,-5,-3) Re; centrality", 10, -0.5, 9.5);
			fsc5353Re[j]->Sumw2();
			fListOfObjects->Add(fsc5353Re[j]);
*/
			for(Int_t h=0; h<5; h++)
			{

        fcn2Re[h][j] = new TProfile(Form("fc%d2_number%dRe", h+2, j), "<<2>> Re; centrality", 10, -0.5, 9.5);
        fcn2Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn2Re[h][j]);
        
				fcn2GapRe[h][j] = new TProfile(Form("fc%d2_number%dGapRe", h+2, j), "<<2>> Re |eta| < 1.0; centrality", 10, -0.5, 9.5);
        fcn2GapRe[h][j]->Sumw2();
        fListOfObjects->Add(fcn2GapRe[h][j]);
        
				fcn2Gap0Re[h][j] = new TProfile(Form("fc%d2_number%dGap0Re", h+2, j), "<<2>> Re |eta| > 0.; centrality", 10, -0.5, 9.5);
        fcn2Gap0Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn2Gap0Re[h][j]);
        
				fcn2Gap4Re[h][j] = new TProfile(Form("fc%d2_number%dGap4Re", h+2, j), "<<2>> Re |eta| > 0.4; centrality", 10, -0.5, 9.5);
        fcn2Gap4Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn2Gap4Re[h][j]);
        
				fcn2Gap8Re[h][j] = new TProfile(Form("fc%d2_number%dGap8Re", h+2, j), "<<2>> Re |eta| > 0.8; centrality", 10, -0.5, 9.5);
        fcn2Gap8Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn2Gap8Re[h][j]);
        
        fcn4Re[h][j] = new TProfile(Form("fc%d4_number%dRe", h+2, j), "<<4>> Re; centrality", 10, -0.5, 9.5);
        fcn4Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn4Re[h][j]);
				
        fcn6Re[h][j] = new TProfile(Form("fc%d6_number%dRe", h+2, j), "<<6>> Re; centrality", 10, -0.5, 9.5);
        fcn6Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn6Re[h][j]);
        
        fcn8Re[h][j] = new TProfile(Form("fc%d8_number%dRe", h+2, j), "<<8>> Re; centrality", 10, -0.5, 9.5);
        fcn8Re[h][j]->Sumw2();
        fListOfObjects->Add(fcn8Re[h][j]);
			

	

			if(fDifferential == true)
			{

				for(int c=0; c<10; c++)
				{

					fdn2Re[h][c][j] = new TProfile(Form("fd%d2Re_cent%d_number%d", h+2, c, j), "<<2>> diff Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Re[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Re[h][c][j]);

					fdn2GapQMRe[h][c][j] = new TProfile(Form("fd%d2GapQMRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|<1.0; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2GapQMRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2GapQMRe[h][c][j]);

					fdn2GapQPRe[h][c][j] = new TProfile(Form("fd%d2GapQPRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|<1.0; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2GapQPRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2GapQPRe[h][c][j]);

					fdn2Gap0QMRe[h][c][j] = new TProfile(Form("fd%d2Gap0QMRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|>0.; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Gap0QMRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Gap0QMRe[h][c][j]);

					fdn2Gap0QPRe[h][c][j] = new TProfile(Form("fd%d2Gap0QPRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|>0.; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Gap0QPRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Gap0QPRe[h][c][j]);

					fdn2Gap4QMRe[h][c][j] = new TProfile(Form("fd%d2Gap4QMRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|>0.4; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Gap4QMRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Gap4QMRe[h][c][j]);

					fdn2Gap4QPRe[h][c][j] = new TProfile(Form("fd%d2Gap4QPRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|>0.4; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Gap4QPRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Gap4QPRe[h][c][j]);

					fdn2Gap8QMRe[h][c][j] = new TProfile(Form("fd%d2Gap8QMRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|>0.8; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Gap8QMRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Gap8QMRe[h][c][j]);

					fdn2Gap8QPRe[h][c][j] = new TProfile(Form("fd%d2Gap8QPRe_cent%d_number%d", h+2, c, j), "<<2>> diff Re |eta|>0.8; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn2Gap8QPRe[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn2Gap8QPRe[h][c][j]);

					fdn4Re[h][c][j] = new TProfile(Form("fd%d4Re_cent%d_number%d", h+2, c, j), "<<4>> diff Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fdn4Re[h][c][j]->Sumw2();
					fListOfObjects->Add(fdn4Re[h][c][j]);

					fcn2PtGap0Re[h][c][j] = new TProfile(Form("fc%d2PtGap0Re_cent%d_number%d", h+2, c, j), "v2[m] eta gap > 0. Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fcn2PtGap0Re[h][c][j]->Sumw2();
					fListOfObjects->Add(fcn2PtGap0Re[h][c][j]);
					
					fcn2PtGap4Re[h][c][j] = new TProfile(Form("fc%d2PtGap4Re_cent%d_number%d", h+2, c, j), "v2[m] eta gap > 0.4 Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fcn2PtGap4Re[h][c][j]->Sumw2();
					fListOfObjects->Add(fcn2PtGap4Re[h][c][j]);
					
					fcn2PtGap8Re[h][c][j] = new TProfile(Form("fc%d2PtGap8Re_cent%d_number%d", h+2, c, j), "v2[m] eta gap > 0.8 Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fcn2PtGap8Re[h][c][j]->Sumw2();
					fListOfObjects->Add(fcn2PtGap8Re[h][c][j]);
					
					fcn2PtGap10Re[h][c][j] = new TProfile(Form("fc%d2PtGap10Re_cent%d_number%d", h+2, c, j), "v2[m] eta gap > 1.0 Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);
					fcn2PtGap10Re[h][c][j]->Sumw2();
					fListOfObjects->Add(fcn2PtGap10Re[h][c][j]);
					
					fcn4PtRe[h][c][j] = new TProfile(Form("fc%d4PtRe_cent%d_number%d", h+2, c, j), "v4[m] Re; pT", fNPtBins, -0.5, fNPtBins - 0.5);


				}// centrality
			}//flag for diff. flow


			}// harmonics
     

    }// random sample

	

  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskChargedFlow::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
		
    if (fAnalysisType == "ESD") {
        Printf("Not implemented\n");
        return;
    }
    
    if (fAnalysisType == "AOD") {
        fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
        if(!fAOD){
            Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
            return;
        }

    
        Float_t zvtx = GetVertex(fAOD);
        
        if(zvtx< -990){
            fVtx->Fill(0);
        }
        else {
            fVtx->Fill(1);
            fVtxBeforeCuts->Fill(zvtx);
            if (TMath::Abs(zvtx) < fVtxCut) {
                fMultCorBeforeCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));
	
                if (fPileUp){
                    
                    Short_t isPileup = fAOD->IsPileupFromSPD(3);
                    if (isPileup != 0)
                        return;
                    
                    //if (fAOD->GetHeader()->GetRefMultiplicityComb08() < 0)
                    //    return;
                    
                    //New cut from Ruben for pileup hybrid
                    if (plpMV(fAOD))
                        return;
                }
                
                
                //fMultCorAfterCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));
                Short_t cenAOD  = GetCentrCode(fAOD);
                //Short_t cenAOD  = 1;
                //cout << cenAOD << endl;
                
                if (cenAOD >= 0){
                    fVtxAfterCuts->Fill(zvtx);
                    AnalyzeAOD(fAOD, cenAOD);
                }
            }
        }
    }
    // Post output data.
    PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskChargedFlow::AnalyzeAOD(AliAODEvent* aod, Short_t centrV0)
{  
	Int_t multTPC = GetTPCMult(aod);
  //fMultvsCentr->Fill(centrV0, multTPC);
    
    
  const Int_t nAODTracks = aod->GetNumberOfTracks();

	//double weight = 1;

	double Qcos[20][20] = {0};
	double Qsin[20][20] = {0};
	double QcosGapM[20][20] = {0};
	double QsinGapM[20][20] = {0};
	double QcosGapP[20][20] = {0};
	double QsinGapP[20][20] = {0};

	double QcosGap0M[20][20] = {0};
	double QsinGap0M[20][20] = {0};
	double QcosGap0P[20][20] = {0};
	double QsinGap0P[20][20] = {0};
  
	double QcosGap4M[20][20] = {0};
	double QsinGap4M[20][20] = {0};
	double QcosGap4P[20][20] = {0};
	double QsinGap4P[20][20] = {0};
	
	double QcosGap8M[20][20] = {0};
	double QsinGap8M[20][20] = {0};
	double QcosGap8P[20][20] = {0};
	double QsinGap8P[20][20] = {0};

	// adding 
    
  Int_t fBin = 0;
  
  TRandom3 rr(0);
  Double_t ranNum = 0.;
  ranNum = rr.Rndm();
  
  if(fSample)
  {
    //cout << rr.Rndm() << endl;
      
   	//if (ranNum > 0 && ranNum < 1) fBin = 0;
   	if (ranNum <= 0.1 && ranNum > 0. ) fBin = 1;
		else if (ranNum <= 0.2 && ranNum > 0.1) fBin = 2;
		else if (ranNum <= 0.3 && ranNum > 0.2) fBin = 3;
		else if (ranNum <= 0.4 && ranNum > 0.3) fBin = 4;
		else if (ranNum <= 0.5 && ranNum > 0.4) fBin = 5;
		else if (ranNum <= 0.6 && ranNum > 0.5) fBin = 6;
		else if (ranNum <= 0.7 && ranNum > 0.6) fBin = 7;
		else if (ranNum <= 0.8 && ranNum > 0.7) fBin = 8;
		else if (ranNum <= 0.9 && ranNum > 0.8) fBin = 9;
		else if (ranNum > 0.9 && ranNum < 1.) fBin = 10;
		else 
		{
      cout << ranNum << endl;       
      cout << "This is very strange!!!" << endl;
      fBin = 11;
	  }
  }
  
  //return; // WORKS

  //if(fSample && fBin = 0) continue;

	//..LOOP OVER TRACKS........	
	//........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++)
  {
      AliAODTrack* aodTrk =(AliAODTrack*) aod->GetTrack(nt);
      
        if (!aodTrk){
          delete aodTrk;
          continue;
      }
      
      if (!(aodTrk->TestFilterBit(fFilterbit)))
          continue;
      
      if ((TMath::Abs(aodTrk->Eta()) > fEtaCut) || (aodTrk->GetTPCNcls() < fNoClus) || (aodTrk->Pt() < fMinPt) || (aodTrk->Pt() > fMaxPt))
          continue;
		
        
        
      if ( (fULS == 1 ) || (fULS == 2 && aodTrk->Charge()>0)  ||  (fULS == 3 && aodTrk->Charge() < 0))
	  {
          
/*			TRandom3 ran(0);
	
			if(aodTrk->Phi() > 1 && aodTrk->Phi() < 2)
			{
	
				double num = ran.Rndm();	
				if(num > 0.5) continue;

			}
*/


      fPhiDis->Fill(aodTrk->Phi());
      fEtaDis->Fill(aodTrk->Eta());
      
			if(fVtxCut == 10.0) fPhiDisVtx10->Fill(aodTrk->Phi());
			if(fVtxCut == 9.0) fPhiDisVtx9->Fill(aodTrk->Phi());
			if(fVtxCut == 8.0) fPhiDisVtx8->Fill(aodTrk->Phi());
			if(fVtxCut == 7.0) fPhiDisVtx7->Fill(aodTrk->Phi());
 
			if(aodTrk->Pt() >= 0.2 && aodTrk->Pt() < 1.)	fPhiDispT2->Fill(aodTrk->Phi());
			if(aodTrk->Pt() >= 1. && aodTrk->Pt() < 2.) fPhiDispT3->Fill(aodTrk->Phi());	
			if(aodTrk->Pt() >= 2. && aodTrk->Pt() <= 5.) fPhiDispT5->Fill(aodTrk->Phi());

			fPhiDisCent[centrV0]->Fill(aodTrk->Phi());			
	
 
			double weight = GetWeight(aodTrk->Phi(), 0);

			fWeight->Fill(aodTrk->Phi()*weight);
			
			//..calculate Q-vectors
			for(int iharm=0; iharm<20; iharm++)
			{
				for(int ipow=0; ipow<20; ipow++)
				{
					Qcos[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*aodTrk->Phi());
					Qsin[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*aodTrk->Phi());
				}
			}

      //Gap = 0.0
      if (aodTrk->Eta() > 0.)
      {
         
				fPhiDisGap0P->Fill(aodTrk->Phi());

				double weight0P = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0P[iharm][ipow] += TMath::Power(weight0P, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap0P[iharm][ipow] += TMath::Power(weight0P, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
 
      }
      if (aodTrk->Eta() < -0.)
      {
				fPhiDisGap0M->Fill(aodTrk->Phi());

				double weight0M = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0M[iharm][ipow] += TMath::Power(weight0M, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap0M[iharm][ipow] += TMath::Power(weight0M, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
		
      // Gap = 0.4
      if (aodTrk->Eta() > 0.2)
      {
				fPhiDisGap4P->Fill(aodTrk->Phi());

				double weight4P = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4P[iharm][ipow] += TMath::Power(weight4P, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap4P[iharm][ipow] += TMath::Power(weight4P, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
      if (aodTrk->Eta() < -0.2)
      {
				fPhiDisGap4M->Fill(aodTrk->Phi());

				double weight4M = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4M[iharm][ipow] += TMath::Power(weight4M, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap4M[iharm][ipow] += TMath::Power(weight4M, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
	 	
      // Gap = 0.8
      if (aodTrk->Eta() > 0.4)
      {
				fPhiDisGap8P->Fill(aodTrk->Phi());

				double weight8P = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8P[iharm][ipow] += TMath::Power(weight8P, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap8P[iharm][ipow] += TMath::Power(weight8P, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
      }
          
      if (aodTrk->Eta() < -0.4)
      {
				fPhiDisGap8M->Fill(aodTrk->Phi());

				double weight8M = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8M[iharm][ipow] += TMath::Power(weight8M, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGap8M[iharm][ipow] += TMath::Power(weight8M, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}

      }


			//..if eta < -0.5
			if(aodTrk->Eta() < -0.5)
			{
				fPhiDisGapM->Fill(aodTrk->Phi());

				double weightM = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapM[iharm][ipow] += TMath::Power(weightM, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGapM[iharm][ipow] += TMath::Power(weightM, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}
			} else
			{
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapM[iharm][ipow] += 0;
						QsinGapM[iharm][ipow] += 0;
					}
				}	
			}
	
			//..if eta > 0.5
			if(aodTrk->Eta() > 0.5)
			{
				fPhiDisGapP->Fill(aodTrk->Phi());

				double weightP = GetWeight(aodTrk->Phi(), aodTrk->Eta());

				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapP[iharm][ipow] += TMath::Power(weightP, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						QsinGapP[iharm][ipow] += TMath::Power(weightP, ipow)*TMath::Sin(iharm*aodTrk->Phi());	
					}
				}
			} else
			{
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapP[iharm][ipow] += 0;
						QsinGapP[iharm][ipow] += 0;
					}
				}
			}


	
    } // end loop of charged particles
        
  } // end loop of all track

  //return; // WORKS !

	//............................
	//..GENERIC FRAMEWORK RP
	//............................
   
	//..calculate Q-vector for each harmonics n and power p
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			Qvector[iharm][ipow] = TComplex(Qcos[iharm][ipow], Qsin[iharm][ipow]);
			QvectorM[iharm][ipow] = TComplex(QcosGapM[iharm][ipow], QsinGapM[iharm][ipow]);
			QvectorP[iharm][ipow] = TComplex(QcosGapP[iharm][ipow], QsinGapP[iharm][ipow]);
			Qvector0M[iharm][ipow] = TComplex(QcosGap0M[iharm][ipow], QsinGap0M[iharm][ipow]);
			Qvector0P[iharm][ipow] = TComplex(QcosGap0P[iharm][ipow], QsinGap0P[iharm][ipow]);
			Qvector4M[iharm][ipow] = TComplex(QcosGap4M[iharm][ipow], QsinGap4M[iharm][ipow]);
			Qvector4P[iharm][ipow] = TComplex(QcosGap4P[iharm][ipow], QsinGap4P[iharm][ipow]);
			Qvector8M[iharm][ipow] = TComplex(QcosGap8M[iharm][ipow], QsinGap8M[iharm][ipow]);
			Qvector8P[iharm][ipow] = TComplex(QcosGap8P[iharm][ipow], QsinGap8P[iharm][ipow]);
		}
	}
	
	//return; // WORKS 

	//..calculate 2-particle correlations
	//..................................
	double Dn2 = Two(0, 0).Re();
	double Dn2gap = TwoGap(0, 0).Re();
	double Dn2gap0 = TwoGap0(0, 0).Re();
	double Dn2gap4 = TwoGap4(0, 0).Re();
	double Dn2gap8 = TwoGap8(0, 0).Re();

	//return; // NOT WORKING
	
	if(Dn2gap != 0)
	{
		//..v2{2} with eta gap < 1.0
		TComplex v22gap = TwoGap(2, -2);
		double v22Regap = v22gap.Re()/Dn2gap;
		fcn2GapRe[0][fBin]->Fill(centrV0, v22Regap, Dn2gap);

		//..v3{2} with eta gap < 1.0
		TComplex v32gap = TwoGap(3, -3);
		double v32Regap = v32gap.Re()/Dn2gap;
		fcn2GapRe[1][fBin]->Fill(centrV0, v32Regap, Dn2gap);

		//..v4{2} with eta gap < 1.0
		TComplex v42gap = TwoGap(4, -4);
		double v42Regap = v42gap.Re()/Dn2gap;
		fcn2GapRe[2][fBin]->Fill(centrV0, v42Regap, Dn2gap);

		//..v5{2} with eta gap < 1.0
		TComplex v52gap = TwoGap(5, -5);
		double v52Regap = v52gap.Re()/Dn2gap;
		fcn2GapRe[3][fBin]->Fill(centrV0, v52Regap, Dn2gap);

		//..v6{2} with eta gap < 1.0
		TComplex v62gap = TwoGap(6, -6);
		double v62Regap = v62gap.Re()/Dn2gap;
		fcn2GapRe[4][fBin]->Fill(centrV0, v62Regap, Dn2gap);

	}

	//return; // NOT WORKING

	if(Dn2gap0 != 0)
	{
		//..v2{2} with eta gap < 0.
		TComplex v22gap0 = TwoGap0(2, -2);
		double v22Regap0 = v22gap0.Re()/Dn2gap0;
		fcn2Gap0Re[0][fBin]->Fill(centrV0, v22Regap0, Dn2gap0);

		//..v3{2} with eta gap < 0.
		TComplex v32gap0 = TwoGap0(3, -3);
		double v32Regap0 = v32gap0.Re()/Dn2gap0;
		fcn2Gap0Re[1][fBin]->Fill(centrV0, v32Regap0, Dn2gap0);

		//..v4{2} with eta gap < 0.
		TComplex v42gap0 = TwoGap0(4, -4);
		double v42Regap0 = v42gap0.Re()/Dn2gap0;
		fcn2Gap0Re[2][fBin]->Fill(centrV0, v42Regap0, Dn2gap0);

		//..v5{2} with eta gap < 0.
		TComplex v52gap0 = TwoGap0(5, -5);
		double v52Regap0 = v52gap0.Re()/Dn2gap0;
		fcn2Gap0Re[3][fBin]->Fill(centrV0, v52Regap0, Dn2gap0);

		//..v6{2} with eta gap < 0.
		TComplex v62gap0 = TwoGap0(6, -6);
		double v62Regap0 = v62gap0.Re()/Dn2gap0;
		fcn2Gap0Re[4][fBin]->Fill(centrV0, v62Regap0, Dn2gap0);

	}
	

	//return; //NOT WORKING


	if(Dn2gap4 != 0)
	{
		//..v2{2} with eta gap > 0.4
		TComplex v22gap4 = TwoGap4(2, -2);
		double v22Regap4 = v22gap4.Re()/Dn2gap4;
		fcn2Gap4Re[0][fBin]->Fill(centrV0, v22Regap4, Dn2gap4);

		//..v3{2} with eta gap > 0.4
		TComplex v32gap4 = TwoGap4(3, -3);
		double v32Regap4 = v32gap4.Re()/Dn2gap4;
		fcn2Gap4Re[1][fBin]->Fill(centrV0, v32Regap4, Dn2gap4);

		//..v4{2} with eta gap > 0.4
		TComplex v42gap4 = TwoGap4(4, -4);
		double v42Regap4 = v42gap4.Re()/Dn2gap4;
		fcn2Gap4Re[2][fBin]->Fill(centrV0, v42Regap4, Dn2gap4);

		//..v5{2} with eta gap > 0.4
		TComplex v52gap4 = TwoGap4(5, -5);
		double v52Regap4 = v52gap4.Re()/Dn2gap4;
		fcn2Gap4Re[3][fBin]->Fill(centrV0, v52Regap4, Dn2gap4);

		//..v6{2} with eta gap > 0.4
		TComplex v62gap4 = TwoGap4(6, -6);
		double v62Regap4 = v62gap4.Re()/Dn2gap4;
		fcn2Gap4Re[4][fBin]->Fill(centrV0, v62Regap4, Dn2gap4);

	}

	//return; // NOT WORKING

	if(Dn2gap8 != 0)
	{
		//..v2{2} with eta gap > 0.8
		TComplex v22gap8 = TwoGap8(2, -2);
		double v22Regap8 = v22gap8.Re()/Dn2gap8;
		fcn2Gap8Re[0][fBin]->Fill(centrV0, v22Regap8, Dn2gap8);

		//..v3{2} with eta gap > 0.8
		TComplex v32gap8 = TwoGap8(3, -3);
		double v32Regap8 = v32gap8.Re()/Dn2gap8;
		fcn2Gap8Re[1][fBin]->Fill(centrV0, v32Regap8, Dn2gap8);

		//..v4{2} with eta gap > 0.8
		TComplex v42gap8 = TwoGap8(4, -4);
		double v42Regap8 = v42gap8.Re()/Dn2gap8;
		fcn2Gap8Re[2][fBin]->Fill(centrV0, v42Regap8, Dn2gap8);

		//..v5{2} with eta gap > 0.8
		TComplex v52gap8 = TwoGap8(5, -5);
		double v52Regap8 = v52gap8.Re()/Dn2gap8;
		fcn2Gap8Re[3][fBin]->Fill(centrV0, v52Regap8, Dn2gap8);

		//..v6{2} with eta gap > 0.8
		TComplex v62gap8 = TwoGap8(6, -6);
		double v62Regap8 = v62gap8.Re()/Dn2gap8;
		fcn2Gap8Re[4][fBin]->Fill(centrV0, v62Regap8, Dn2gap8);

	}

	//return; //NOT WORKING

	if(Dn2 != 0)
	{	
		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		fcn2Re[0][fBin]->Fill(centrV0, v22Re, Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		fcn2Re[1][fBin]->Fill(centrV0, v32Re, Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		fcn2Re[2][fBin]->Fill(centrV0, v42Re, Dn2);

		//..v5{2} = <cos5(phi1 - phi2)>
		TComplex v52 = Two(5, -5);
		double v52Re = v52.Re()/Dn2;
		fcn2Re[3][fBin]->Fill(centrV0, v52Re, Dn2);

		//..v6{2} = <cos6(phi1 - phi2)>
		TComplex v62 = Two(6, -6);
		double v62Re = v62.Re()/Dn2;
		fcn2Re[4][fBin]->Fill(centrV0, v62Re, Dn2);

	}	

	//return; // NOT WORKING

	if(fDifferential == true)
	{

	//............................
	//..PT DIFFERENTIAL FLOW
	//............................
   
	for(int ipt=0; ipt<fNPtBins; ipt++)
	{

		double pcosPt[20][20] = {0};
		double psinPt[20][20] = {0};
		double pcosPtM[20][20] = {0};
		double psinPtM[20][20] = {0};
		double pcosPtP[20][20] = {0};
		double psinPtP[20][20] = {0};
		
		double pcosPt0M[20][20] = {0};
		double psinPt0M[20][20] = {0};
		double pcosPt0P[20][20] = {0};
		double psinPt0P[20][20] = {0};
		
		double pcosPt4M[20][20] = {0};
		double psinPt4M[20][20] = {0};
		double pcosPt4P[20][20] = {0};
		double psinPt4P[20][20] = {0};
		
		double pcosPt8M[20][20] = {0};
		double psinPt8M[20][20] = {0};
		double pcosPt8P[20][20] = {0};
		double psinPt8P[20][20] = {0};
		
		double qcosPt[20][20] = {0};
		double qsinPt[20][20] = {0};
	

		for(Int_t nt = 0; nt < nAODTracks; nt++)
		{

      AliAODTrack* aodTrk =(AliAODTrack*) aod->GetTrack(nt);
      
        if (!aodTrk){
          delete aodTrk;
          continue;
      	}
      
      if (!(aodTrk->TestFilterBit(fFilterbit)))
          continue;

			if( aodTrk->Pt() <= fPtBins[ipt] || aodTrk->Pt() > fPtBins[ipt+1])
          continue;
      
      if ((TMath::Abs(aodTrk->Eta()) > fEtaCut) || (aodTrk->GetTPCNcls() < fNoClus) || (aodTrk->Pt() < fMinPt) || (aodTrk->Pt() > fMaxPt))
          continue;
		      
      if ( (fULS == 1 ) || (fULS == 2 && aodTrk->Charge()>0)  ||  (fULS == 3 && aodTrk->Charge() < 0))
	  	{
			
				double weightPt = GetWeight(aodTrk->Phi(), 0);

				fPhiPt[ipt]->Fill(aodTrk->Phi());

				//..calculate p-vectors
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						pcosPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						psinPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						qcosPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
						qsinPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
					}
				}

				if(aodTrk->Eta() < -0.5)
				{
					double weightPtM = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<20; iharm++)
					{
						for(int ipow=0; ipow<20; ipow++)
						{
							pcosPtM[iharm][ipow] += TMath::Power(weightPtM, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPtM[iharm][ipow] += TMath::Power(weightPtM, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta < -0.5	

				if(aodTrk->Eta() > 0.5)
				{
					double weightPtP = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<20; iharm++)
					{
						for(int ipow=0; ipow<20; ipow++)
						{
							pcosPtP[iharm][ipow] += TMath::Power(weightPtP, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPtP[iharm][ipow] += TMath::Power(weightPtP, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta > 0.5	
				

				//.. eta gap > 0.8
				if(aodTrk->Eta() < -0.4)
				{
					double weightPt8M = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<10; iharm++)
					{
						for(int ipow=0; ipow<10; ipow++)
						{
							pcosPt8M[iharm][ipow] += TMath::Power(weightPt8M, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPt8M[iharm][ipow] += TMath::Power(weightPt8M, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta < -0.4	

				if(aodTrk->Eta() > 0.4)
				{
					double weightPt8P = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<10; iharm++)
					{
						for(int ipow=0; ipow<10; ipow++)
						{
							pcosPt8P[iharm][ipow] += TMath::Power(weightPt8P, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPt8P[iharm][ipow] += TMath::Power(weightPt8P, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta > 0.4

				//.. eta gap > 0.4
				if(aodTrk->Eta() < -0.2)
				{
					double weightPt4M = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<10; iharm++)
					{
						for(int ipow=0; ipow<10; ipow++)
						{
							pcosPt4M[iharm][ipow] += TMath::Power(weightPt4M, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPt4M[iharm][ipow] += TMath::Power(weightPt4M, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta < -0.2	

				if(aodTrk->Eta() > 0.2)
				{
					double weightPt4P = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<10; iharm++)
					{
						for(int ipow=0; ipow<10; ipow++)
						{
							pcosPt4P[iharm][ipow] += TMath::Power(weightPt4P, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPt4P[iharm][ipow] += TMath::Power(weightPt4P, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta > 0.2	

				//.. eta gap > 0.0
				if(aodTrk->Eta() < -0.)
				{
					double weightPt0M = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<10; iharm++)
					{
						for(int ipow=0; ipow<10; ipow++)
						{
							pcosPt0M[iharm][ipow] += TMath::Power(weightPt0M, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPt0M[iharm][ipow] += TMath::Power(weightPt0M, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta < -0.	

				if(aodTrk->Eta() > 0.)
				{
					double weightPt0P = GetWeight(aodTrk->Phi(), aodTrk->Eta());

					for(int iharm=0; iharm<10; iharm++)
					{
						for(int ipow=0; ipow<10; ipow++)
						{
							pcosPt0P[iharm][ipow] += TMath::Power(weightPt0P, ipow)*TMath::Cos(iharm*aodTrk->Phi());
							psinPt0P[iharm][ipow] += TMath::Power(weightPt0P, ipow)*TMath::Sin(iharm*aodTrk->Phi());
						}
					}
				}	//..eta > 0.

			}//..charged particles
	
		}//..loop over tracks


		//return; // NOT WORKS

		//............................
		//..GENERIC FRAMEWORK POI
		//............................

		//..calculate p-vectors, q-vectors for pT-differential flow
		for(int iharm=0; iharm<20; iharm++)
		{
			for(int ipow=0; ipow<20; ipow++)
			{
				pvector[iharm][ipow] = TComplex(pcosPt[iharm][ipow], psinPt[iharm][ipow]);
				pvectorM[iharm][ipow] = TComplex(pcosPtM[iharm][ipow], psinPtM[iharm][ipow]);
				pvectorP[iharm][ipow] = TComplex(pcosPtP[iharm][ipow], psinPtP[iharm][ipow]);
				pvector0M[iharm][ipow] = TComplex(pcosPt0M[iharm][ipow], psinPt0M[iharm][ipow]);
				pvector0P[iharm][ipow] = TComplex(pcosPt0P[iharm][ipow], psinPt0P[iharm][ipow]);
				pvector4M[iharm][ipow] = TComplex(pcosPt4M[iharm][ipow], psinPt4M[iharm][ipow]);
				pvector4P[iharm][ipow] = TComplex(pcosPt4P[iharm][ipow], psinPt4P[iharm][ipow]);
				pvector8M[iharm][ipow] = TComplex(pcosPt8M[iharm][ipow], psinPt8M[iharm][ipow]);
				pvector8P[iharm][ipow] = TComplex(pcosPt8P[iharm][ipow], psinPt8P[iharm][ipow]);
				qvector[iharm][ipow] = TComplex(qcosPt[iharm][ipow], qsinPt[iharm][ipow]);
			}
		}

		//..calculate 2-particle correlations
		//...................................
		double DDn2 = TwoDiff(0, 0).Re();
		double DDn2gapQM = TwoDiffGapQM(0, 0).Re();
		double DDn2gapQP = TwoDiffGapQP(0, 0).Re();
		double DDn2gap0QM = TwoDiffGap0QM(0, 0).Re();
		double DDn2gap0QP = TwoDiffGap0QP(0, 0).Re();
		double DDn2gap4QM = TwoDiffGap4QM(0, 0).Re();
		double DDn2gap4QP = TwoDiffGap4QP(0, 0).Re();
		double DDn2gap8QM = TwoDiffGap8QM(0, 0).Re();
		double DDn2gap8QP = TwoDiffGap8QP(0, 0).Re();

		double DDn2Ptgap0 = TwoDiffPtGap0(0, 0).Re();
		double DDn2Ptgap4 = TwoDiffPtGap4(0, 0).Re();
		double DDn2Ptgap8 = TwoDiffPtGap8(0, 0).Re();
		double DDn2Ptgap10 = TwoDiffPtGap10(0, 0).Re();

		if(DDn2gapQM != 0)
		{
			//..v2{2} with eta gap < 1.0
			TComplex v22ptGapQM = TwoDiffGapQM(2, -2);
			double v22ptGapReQM = v22ptGapQM.Re()/DDn2gapQM;
			fdn2GapQMRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGapReQM, DDn2gapQM);

			//..v3{2} with eta gap < 1.0
			TComplex v32ptGapQM = TwoDiffGapQM(3, -3);
			double v32ptGapReQM = v32ptGapQM.Re()/DDn2gapQM;
			fdn2GapQMRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGapReQM, DDn2gapQM);

			//..v4{2} with eta gap < 1.0
			TComplex v42ptGapQM = TwoDiffGapQM(4, -4);
			double v42ptGapReQM = v42ptGapQM.Re()/DDn2gapQM;
			fdn2GapQMRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGapReQM, DDn2gapQM);

			//..v5{2} with eta gap < 1.0
			TComplex v52ptGapQM = TwoDiffGapQM(5, -5);
			double v52ptGapReQM = v52ptGapQM.Re()/DDn2gapQM;
			fdn2GapQMRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGapReQM, DDn2gapQM);

			//..v6{2} with eta gap < 1.0
			TComplex v62ptGapQM = TwoDiffGapQM(6, -6);
			double v62ptGapReQM = v62ptGapQM.Re()/DDn2gapQM;
			fdn2GapQMRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGapReQM, DDn2gapQM);

		}

		if(DDn2gapQP != 0)
		{
			//..v2{2} with eta gap < 1.0
			TComplex v22ptGapQP = TwoDiffGapQP(2, -2);
			double v22ptGapReQP = v22ptGapQP.Re()/DDn2gapQP;
			fdn2GapQPRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGapReQP, DDn2gapQP);
	
			//..v3{2} with eta gap < 1.0
			TComplex v32ptGapQP = TwoDiffGapQP(3, -3);
			double v32ptGapReQP = v32ptGapQP.Re()/DDn2gapQP;
			fdn2GapQPRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGapReQP, DDn2gapQP);

			//..v4{2} with eta gap < 1.0
			TComplex v42ptGapQP = TwoDiffGapQP(4, -4);
			double v42ptGapReQP = v42ptGapQP.Re()/DDn2gapQP;
			fdn2GapQPRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGapReQP, DDn2gapQP);

			//..v5{2} with eta gap < 1.0
			TComplex v52ptGapQP = TwoDiffGapQP(5, -5);
			double v52ptGapReQP = v52ptGapQP.Re()/DDn2gapQP;
			fdn2GapQPRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGapReQP, DDn2gapQP);

			//..v6{2} with eta gap < 1.0
			TComplex v62ptGapQP = TwoDiffGapQP(6, -6);
			double v62ptGapReQP = v62ptGapQP.Re()/DDn2gapQP;
			fdn2GapQPRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGapReQP, DDn2gapQP);

		}

		if(DDn2gap0QM != 0)
		{
			//..v2{2} with eta gap > 0.
			TComplex v22ptGap0QM = TwoDiffGap0QM(2, -2);
			double v22ptGap0ReQM = v22ptGap0QM.Re()/DDn2gap0QM;
			fdn2Gap0QMRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGap0ReQM, DDn2gap0QM);

			//..v3{2} with eta gap > 0.
			TComplex v32ptGap0QM = TwoDiffGap0QM(3, -3);
			double v32ptGap0ReQM = v32ptGap0QM.Re()/DDn2gap0QM;
			fdn2Gap0QMRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGap0ReQM, DDn2gap0QM);

			//..v4{2} with eta gap > 0.
			TComplex v42ptGap0QM = TwoDiffGap0QM(4, -4);
			double v42ptGap0ReQM = v42ptGap0QM.Re()/DDn2gap0QM;
			fdn2Gap0QMRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGap0ReQM, DDn2gap0QM);

			//..v5{2} with eta gap > 0.
			TComplex v52ptGap0QM = TwoDiffGap0QM(5, -5);
			double v52ptGap0ReQM = v52ptGap0QM.Re()/DDn2gap0QM;
			fdn2Gap0QMRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGap0ReQM, DDn2gap0QM);

			//..v6{2} with eta gap > 0.
			TComplex v62ptGap0QM = TwoDiffGap0QM(6, -6);
			double v62ptGap0ReQM = v62ptGap0QM.Re()/DDn2gap0QM;
			fdn2Gap0QMRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGap0ReQM, DDn2gap0QM);

		}

		if(DDn2gap0QP != 0)
		{
			//..v2{2} with eta gap > 0.
			TComplex v22ptGap0QP = TwoDiffGap0QP(2, -2);
			double v22ptGap0ReQP = v22ptGap0QP.Re()/DDn2gap0QP;
			fdn2Gap0QPRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGap0ReQP, DDn2gap0QP);
	
			//..v3{2} with eta gap > 0.
			TComplex v32ptGap0QP = TwoDiffGap0QP(3, -3);
			double v32ptGap0ReQP = v32ptGap0QP.Re()/DDn2gap0QP;
			fdn2Gap0QPRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGap0ReQP, DDn2gap0QP);

			//..v4{2} with eta gap > 0.
			TComplex v42ptGap0QP = TwoDiffGap0QP(4, -4);
			double v42ptGap0ReQP = v42ptGap0QP.Re()/DDn2gap0QP;
			fdn2Gap0QPRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGap0ReQP, DDn2gap0QP);

			//..v5{2} with eta gap > 0.
			TComplex v52ptGap0QP = TwoDiffGap0QP(5, -5);
			double v52ptGap0ReQP = v52ptGap0QP.Re()/DDn2gap0QP;
			fdn2Gap0QPRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGap0ReQP, DDn2gap0QP);

			//..v6{2} with eta gap > 0.
			TComplex v62ptGap0QP = TwoDiffGap0QP(6, -6);
			double v62ptGap0ReQP = v62ptGap0QP.Re()/DDn2gap0QP;
			fdn2Gap0QPRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGap0ReQP, DDn2gap0QP);

		}


		if(DDn2gap4QM != 0)
		{
			//..v2{2} with eta gap > 0.4
			TComplex v22ptGap4QM = TwoDiffGap4QM(2, -2);
			double v22ptGap4ReQM = v22ptGap4QM.Re()/DDn2gap4QM;
			fdn2Gap4QMRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGap4ReQM, DDn2gap4QM);

			//..v3{2} with eta gap > 0.4
			TComplex v32ptGap4QM = TwoDiffGap4QM(3, -3);
			double v32ptGap4ReQM = v32ptGap4QM.Re()/DDn2gap4QM;
			fdn2Gap4QMRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGap4ReQM, DDn2gap4QM);

			//..v4{2} with eta gap > 0.4
			TComplex v42ptGap4QM = TwoDiffGap4QM(4, -4);
			double v42ptGap4ReQM = v42ptGap4QM.Re()/DDn2gap4QM;
			fdn2Gap4QMRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGap4ReQM, DDn2gap4QM);

			//..v5{2} with eta gap > 0.4
			TComplex v52ptGap4QM = TwoDiffGap4QM(5, -5);
			double v52ptGap4ReQM = v52ptGap4QM.Re()/DDn2gap4QM;
			fdn2Gap4QMRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGap4ReQM, DDn2gap4QM);

			//..v6{2} with eta gap > 0.4
			TComplex v62ptGap4QM = TwoDiffGap4QM(6, -6);
			double v62ptGap4ReQM = v62ptGap4QM.Re()/DDn2gap4QM;
			fdn2Gap4QMRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGap4ReQM, DDn2gap4QM);

		}

		if(DDn2gap4QP != 0)
		{
			//..v2{2} with eta gap > 0.4
			TComplex v22ptGap4QP = TwoDiffGap4QP(2, -2);
			double v22ptGap4ReQP = v22ptGap4QP.Re()/DDn2gap4QP;
			fdn2Gap4QPRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGap4ReQP, DDn2gap4QP);
	
			//..v3{2} with eta gap > 0.4
			TComplex v32ptGap4QP = TwoDiffGap4QP(3, -3);
			double v32ptGap4ReQP = v32ptGap4QP.Re()/DDn2gap4QP;
			fdn2Gap4QPRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGap4ReQP, DDn2gap4QP);

			//..v4{2} with eta gap > 0.4
			TComplex v42ptGap4QP = TwoDiffGap4QP(4, -4);
			double v42ptGap4ReQP = v42ptGap4QP.Re()/DDn2gap4QP;
			fdn2Gap4QPRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGap4ReQP, DDn2gap4QP);

			//..v5{2} with eta gap > 0.4
			TComplex v52ptGap4QP = TwoDiffGap4QP(5, -5);
			double v52ptGap4ReQP = v52ptGap4QP.Re()/DDn2gap4QP;
			fdn2Gap4QPRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGap4ReQP, DDn2gap4QP);

			//..v6{2} with eta gap > 0.4
			TComplex v62ptGap4QP = TwoDiffGap4QP(6, -6);
			double v62ptGap4ReQP = v62ptGap4QP.Re()/DDn2gap4QP;
			fdn2Gap4QPRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGap4ReQP, DDn2gap4QP);

		}

		if(DDn2gap8QM != 0)
		{
			//..v2{2} with eta gap > 0.8
			TComplex v22ptGap8QM = TwoDiffGap8QM(2, -2);
			double v22ptGap8ReQM = v22ptGap8QM.Re()/DDn2gap8QM;
			fdn2Gap8QMRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGap8ReQM, DDn2gap8QM);

			//..v3{2} with eta gap > 0.8
			TComplex v32ptGap8QM = TwoDiffGap8QM(3, -3);
			double v32ptGap8ReQM = v32ptGap8QM.Re()/DDn2gap8QM;
			fdn2Gap8QMRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGap8ReQM, DDn2gap8QM);

			//..v4{2} with eta gap > 0.8
			TComplex v42ptGap8QM = TwoDiffGap8QM(4, -4);
			double v42ptGap8ReQM = v42ptGap8QM.Re()/DDn2gap8QM;
			fdn2Gap8QMRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGap8ReQM, DDn2gap8QM);

			//..v5{2} with eta gap > 0.8
			TComplex v52ptGap8QM = TwoDiffGap8QM(5, -5);
			double v52ptGap8ReQM = v52ptGap8QM.Re()/DDn2gap8QM;
			fdn2Gap8QMRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGap8ReQM, DDn2gap8QM);

			//..v6{2} with eta gap > 0.8
			TComplex v62ptGap8QM = TwoDiffGap8QM(6, -6);
			double v62ptGap8ReQM = v62ptGap8QM.Re()/DDn2gap8QM;
			fdn2Gap8QMRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGap8ReQM, DDn2gap8QM);

		}

		if(DDn2gap8QP != 0)
		{
			//..v2{2} with eta gap > 0.8
			TComplex v22ptGap8QP = TwoDiffGap8QP(2, -2);
			double v22ptGap8ReQP = v22ptGap8QP.Re()/DDn2gap8QP;
			fdn2Gap8QPRe[0][centrV0][fBin]->Fill(double(ipt), v22ptGap8ReQP, DDn2gap8QP);
	
			//..v3{2} with eta gap > 0.8
			TComplex v32ptGap8QP = TwoDiffGap8QP(3, -3);
			double v32ptGap8ReQP = v32ptGap8QP.Re()/DDn2gap8QP;
			fdn2Gap8QPRe[1][centrV0][fBin]->Fill(double(ipt), v32ptGap8ReQP, DDn2gap8QP);

			//..v4{2} with eta gap > 0.8
			TComplex v42ptGap8QP = TwoDiffGap8QP(4, -4);
			double v42ptGap8ReQP = v42ptGap8QP.Re()/DDn2gap8QP;
			fdn2Gap8QPRe[2][centrV0][fBin]->Fill(double(ipt), v42ptGap8ReQP, DDn2gap8QP);

			//..v5{2} with eta gap > 0.8
			TComplex v52ptGap8QP = TwoDiffGap8QP(5, -5);
			double v52ptGap8ReQP = v52ptGap8QP.Re()/DDn2gap8QP;
			fdn2Gap8QPRe[3][centrV0][fBin]->Fill(double(ipt), v52ptGap8ReQP, DDn2gap8QP);

			//..v6{2} with eta gap > 0.8
			TComplex v62ptGap8QP = TwoDiffGap8QP(6, -6);
			double v62ptGap8ReQP = v62ptGap8QP.Re()/DDn2gap8QP;
			fdn2Gap8QPRe[4][centrV0][fBin]->Fill(double(ipt), v62ptGap8ReQP, DDn2gap8QP);

		}

		if(DDn2 != 0)
		{
			//..v2{2}
			TComplex v22pt = TwoDiff(2, -2);
			double v22ptRe = v22pt.Re()/DDn2;
			fdn2Re[0][centrV0][fBin]->Fill(double(ipt), v22ptRe, DDn2);	

			//..v3{2}
			TComplex v32pt = TwoDiff(3, -3);
			double v32ptRe = v32pt.Re()/DDn2;
			fdn2Re[1][centrV0][fBin]->Fill(double(ipt), v32ptRe, DDn2);	
	
			//..v4{2}
			TComplex v42pt = TwoDiff(4, -4);
			double v42ptRe = v42pt.Re()/DDn2;
			fdn2Re[2][centrV0][fBin]->Fill(double(ipt), v42ptRe, DDn2);	

			//..v5{2}
			TComplex v52pt = TwoDiff(5, -5);
			double v52ptRe = v52pt.Re()/DDn2;
			fdn2Re[3][centrV0][fBin]->Fill(double(ipt), v52ptRe, DDn2);	

			//..v6{2}
			TComplex v62pt = TwoDiff(6, -6);
			double v62ptRe = v62pt.Re()/DDn2;
			fdn2Re[4][centrV0][fBin]->Fill(double(ipt), v62ptRe, DDn2);	

		}

		//..vn[2] with eta gap > 0.0
		if(DDn2Ptgap0 != 0) 
		{
			//..v2[2]
			TComplex v22DiffptGap0 = TwoDiffPtGap0(2, -2);
			double v22DiffptGap0Re = v22DiffptGap0.Re()/DDn2Ptgap0;
			fcn2PtGap0Re[0][centrV0][fBin]->Fill(double(ipt), v22DiffptGap0Re, DDn2Ptgap0);

			//..v3[2]
			TComplex v32DiffptGap0 = TwoDiffPtGap0(3, -3);
			double v32DiffptGap0Re = v32DiffptGap0.Re()/DDn2Ptgap0;
			fcn2PtGap0Re[1][centrV0][fBin]->Fill(double(ipt), v32DiffptGap0Re, DDn2Ptgap0);

			//..v4[2]
			TComplex v42DiffptGap0 = TwoDiffPtGap0(4, -4);
			double v42DiffptGap0Re = v42DiffptGap0.Re()/DDn2Ptgap0;
			fcn2PtGap0Re[2][centrV0][fBin]->Fill(double(ipt), v42DiffptGap0Re, DDn2Ptgap0);

		}

		//..vn[2] with eta gap > 0.4
		if(DDn2Ptgap4 != 0) 
		{
			//..v2[2]
			TComplex v22DiffptGap4 = TwoDiffPtGap4(2, -2);
			double v22DiffptGap4Re = v22DiffptGap4.Re()/DDn2Ptgap4;
			fcn2PtGap4Re[0][centrV0][fBin]->Fill(double(ipt), v22DiffptGap4Re, DDn2Ptgap4);

			//..v3[2]
			TComplex v32DiffptGap4 = TwoDiffPtGap4(3, -3);
			double v32DiffptGap4Re = v32DiffptGap4.Re()/DDn2Ptgap4;
			fcn2PtGap4Re[1][centrV0][fBin]->Fill(double(ipt), v32DiffptGap4Re, DDn2Ptgap4);

			//..v4[2]
			TComplex v42DiffptGap4 = TwoDiffPtGap4(4, -4);
			double v42DiffptGap4Re = v42DiffptGap4.Re()/DDn2Ptgap4;
			fcn2PtGap4Re[2][centrV0][fBin]->Fill(double(ipt), v42DiffptGap4Re, DDn2Ptgap4);

		}

		
		//..vn[2] with eta gap > 0.8
		if(DDn2Ptgap8 != 0) 
		{
			//..v2[2]
			TComplex v22DiffptGap8 = TwoDiffPtGap8(2, -2);
			double v22DiffptGap8Re = v22DiffptGap8.Re()/DDn2Ptgap8;
			fcn2PtGap8Re[0][centrV0][fBin]->Fill(double(ipt), v22DiffptGap8Re, DDn2Ptgap8);

			//..v3[2]
			TComplex v32DiffptGap8 = TwoDiffPtGap8(3, -3);
			double v32DiffptGap8Re = v32DiffptGap8.Re()/DDn2Ptgap8;
			fcn2PtGap8Re[1][centrV0][fBin]->Fill(double(ipt), v32DiffptGap8Re, DDn2Ptgap8);

			//..v4[2]
			TComplex v42DiffptGap8 = TwoDiffPtGap8(4, -4);
			double v42DiffptGap8Re = v42DiffptGap8.Re()/DDn2Ptgap8;
			fcn2PtGap8Re[2][centrV0][fBin]->Fill(double(ipt), v42DiffptGap8Re, DDn2Ptgap8);

		}
		
		//..vn[2] with eta gap > 1.0
		if(DDn2Ptgap10 != 0) 
		{
			//..v2[2]
			TComplex v22DiffptGap10 = TwoDiffPtGap10(2, -2);
			double v22DiffptGap10Re = v22DiffptGap10.Re()/DDn2Ptgap10;
			fcn2PtGap10Re[0][centrV0][fBin]->Fill(double(ipt), v22DiffptGap10Re, DDn2Ptgap10);

			//..v3[2]
			TComplex v32DiffptGap10 = TwoDiffPtGap10(3, -3);
			double v32DiffptGap10Re = v32DiffptGap10.Re()/DDn2Ptgap10;
			fcn2PtGap10Re[1][centrV0][fBin]->Fill(double(ipt), v32DiffptGap10Re, DDn2Ptgap10);

			//..v4[2]
			TComplex v42DiffptGap10 = TwoDiffPtGap10(4, -4);
			double v42DiffptGap10Re = v42DiffptGap10.Re()/DDn2Ptgap10;
			fcn2PtGap10Re[2][centrV0][fBin]->Fill(double(ipt), v42DiffptGap10Re, DDn2Ptgap10);

		}
	


	}//..loop over pT
 

	}//flag for differential flow

}
//____________________________________________________________________
//	END OF MAIN PROGRAM
//____________________________________________________________________
double AliAnalysisTaskChargedFlow::GetWeight(double phi, double eta)
{

	double weight = 0;
	double Wbin[25] = {0};
//	if(phi > 1.0 && phi < 2.0) weight = 2;
//	else weight = 1;

	double bin[25] = {0, 0.2516, 0.5032, 0.7548, 1.0064, 1.258, 1.5096, 1.7612, 
										2.0128, 2.2644, 2.516, 2.7676, 3.0192, 3.2708, 3.5224, 3.774, 
										4.0256, 4.2772, 4.5288, 4.7804, 5.032, 5.2836, 5.5352, 5.7868, 6.0384};

	if(fFilterbit == 96)
	{	//..global tracks
		
		if(eta < -0.)
		{
			double w[25] = {0.995612, 0.97966, 0.959755, 0.966339, 0.966736, 0.982394, 0.976147,
											1.02237, 1.06288, 1.05803, 0.996775, 1.00295, 0.990918, 0.977478, 
											0.987327, 1.03768, 1.02749, 1.07816, 0.961319, 0.991539, 0.989109,
											1.01855, 0.993953, 0.963159, 1.01367};
			for(int i=0; i<25; i++) Wbin[i] = w[i];
		}
		if(eta > 0.)
		{
			double w[25] = {0.966526, 1.04342, 1.07812, 0.988734, 0.998484, 0.982507, 0.983126,
											0.987193, 1.03048, 1.0203, 0.989095, 0.98032, 0.975003, 0.99697, 
											0.998521, 1.01176, 0.984994, 0.990083, 0.977393, 1.00804, 0.992916,
											1.02921, 0.994246, 0.974084, 1.01848};
			for(int i=0; i<25; i++) Wbin[i] = w[i];
		}
		if(eta == 0)
		{
			double w[25] = {0.976284, 1.06169, 1.00759, 0.976296, 0.980275, 0.984214, 0.983922,
												1.03738, 1.04017, 1.05883, 1.00825, 0.990612, 0.987512, 0.999679,
												0.983449, 1.00973, 0.991457, 1.01028, 0.960657, 0.999427, 0.985977, 
												1.0003, 0.991671, 0.965225, 1.01013};
			for(int i=0; i<25; i++) Wbin[i] = w[i];
		}
	}
	
	if(fFilterbit == 768)
	{	//..hybrid tracks
	
		if(eta < -0.)
		{
			double w[25] = {1.01152, 0.995529, 0.974893, 0.982568, 0.982356, 0.97648, 0.974624,
											1.00904, 0.977344, 1.04594, 1.01359, 1.02252, 1.00463, 0.98731,
											0.990769, 1.03482, 1.01629, 1.01836, 0.974305, 0.99233, 0.99382,
											1.01884, 0.99886, 0.978475, 1.02478};
			for(int i=0; i<25; i++) Wbin[i] = w[i];
		}
		if(eta > 0.)
		{
			double w[25] = {0.977171, 1.05721, 1.09089, 1.00019, 1.01046, 0.972998, 0.973528, 0.983498,
											1.0042, 1.02175, 0.994643, 0.98969, 0.982858, 0.995824, 0.98121, 0.993774,
											0.985768, 0.982977, 0.986135, 1.00766, 0.99467, 1.01105, 0.993567,
											0.984228, 1.02405};
			for(int i=0; i<25; i++) Wbin[i] = w[i];
		}
		if(eta == 0)
		{
			double w[25] = {0.987957, 1.07694, 1.021, 0.988812, 0.99264, 0.976635, 0.978696, 1.02326,
												0.986807, 1.05396, 1.01831, 1.0032, 0.997414, 1.00526, 0.978013, 1.00075,
												0.986894, 0.984098, 0.970481, 0.999465, 0.988175, 0.991922, 0.994515,
												0.977144, 1.01765};
			for(int i=0; i<25; i++) Wbin[i] = w[i];
		}
	}
	
	if(fFilterbit == 128)
	{	//..no weight (weight = 1)
		for(int i=0; i<25; i++) Wbin[i] = 1;
	}
/*	
	if(fNUA == 32)
	{//FB32 (Ante)
		double w[25] = {0.954382, 1.03806, 0.987299, 0.9545, 0.957934, 0.98018, 0.96671, 1.03826,
										1.28462, 1.09639, 0.99588, 0.969805, 0.965473, 0.980728, 0.969466, 1.028, 
										0.977306, 1.01541, 0.940225, 1.00162, 0.982006, 0.996996, 0.979599, 0.944637,
										0.99451};
		for(int i=0; i<25; i++) Wbin[i] = w[i];
	}
*/


	if(phi >= bin[0] && phi < bin[1]) weight = Wbin[0];
	else if(phi >= bin[1] && phi < bin[2]) weight = Wbin[1];
	else if(phi >= bin[2] && phi < bin[3]) weight = Wbin[2];
	else if(phi >= bin[3] && phi < bin[4]) weight = Wbin[3];
	else if(phi >= bin[4] && phi < bin[5]) weight = Wbin[4];
	else if(phi >= bin[5] && phi < bin[6]) weight = Wbin[5];
	else if(phi >= bin[6] && phi < bin[7]) weight = Wbin[6];
	else if(phi >= bin[7] && phi < bin[8]) weight = Wbin[7];
	else if(phi >= bin[8] && phi < bin[9]) weight = Wbin[8];
	else if(phi >= bin[9] && phi < bin[10]) weight = Wbin[9];
	else if(phi >= bin[10] && phi < bin[11]) weight = Wbin[10];
	else if(phi >= bin[11] && phi < bin[12]) weight = Wbin[11];
	else if(phi >= bin[12] && phi < bin[13]) weight = Wbin[12];
	else if(phi >= bin[13] && phi < bin[14]) weight = Wbin[13];
	else if(phi >= bin[14] && phi < bin[15]) weight = Wbin[14];
	else if(phi >= bin[15] && phi < bin[16]) weight = Wbin[15];
	else if(phi >= bin[16] && phi < bin[17]) weight = Wbin[16];
	else if(phi >= bin[17] && phi < bin[18]) weight = Wbin[17];
	else if(phi >= bin[18] && phi < bin[19]) weight = Wbin[18];
	else if(phi >= bin[19] && phi < bin[20]) weight = Wbin[19];
	else if(phi >= bin[20] && phi < bin[21]) weight = Wbin[20];
	else if(phi >= bin[21] && phi < bin[22]) weight = Wbin[21];
	else if(phi >= bin[22] && phi < bin[23]) weight = Wbin[22];
	else if(phi >= bin[23] && phi < bin[24]) weight = Wbin[23];
	else if(phi >= bin[24]) weight = Wbin[24];

	return weight;

}
//_____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGapM(int n, int p)
{

	if(n>=0) return QvectorM[n][p];
  else return TComplex::Conjugate(QvectorM[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGapP(int n, int p)
{

	if(n>=0) return QvectorP[n][p];
  else return TComplex::Conjugate(QvectorP[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap0M(int n, int p)
{

	if(n>=0) return Qvector0M[n][p];
  else return TComplex::Conjugate(Qvector0M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap0P(int n, int p)
{

	if(n>=0) return Qvector0P[n][p];
  else return TComplex::Conjugate(Qvector0P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap4M(int n, int p)
{

	if(n>=0) return Qvector4M[n][p];
  else return TComplex::Conjugate(Qvector4M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap4P(int n, int p)
{

	if(n>=0) return Qvector4P[n][p];
  else return TComplex::Conjugate(Qvector4P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap8M(int n, int p)
{

	if(n>=0) return Qvector8M[n][p];
  else return TComplex::Conjugate(Qvector8M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap8P(int n, int p)
{

	if(n>=0) return Qvector8P[n][p];
  else return TComplex::Conjugate(Qvector8P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::p(int n, int p)
{

	if(n>=0) return pvector[n][p];
	else return TComplex::Conjugate(pvector[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGapM(int n, int p)
{

	if(n>=0) return pvectorM[n][p];
	else return TComplex::Conjugate(pvectorM[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGapP(int n, int p)
{

	if(n>=0) return pvectorP[n][p];
	else return TComplex::Conjugate(pvectorP[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::q(int n, int p)
{

	if(n>=0) return qvector[n][p];
	else return TComplex::Conjugate(qvector[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGap0M(int n, int p)
{

	if(n>=0) return pvector0M[n][p];
	else return TComplex::Conjugate(pvector0M[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGap0P(int n, int p)
{

	if(n>=0) return pvector0P[n][p];
	else return TComplex::Conjugate(pvector0P[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGap4M(int n, int p)
{

	if(n>=0) return pvector4M[n][p];
	else return TComplex::Conjugate(pvector4M[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGap4P(int n, int p)
{

	if(n>=0) return pvector4P[n][p];
	else return TComplex::Conjugate(pvector4P[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGap8M(int n, int p)
{

	if(n>=0) return pvector8M[n][p];
	else return TComplex::Conjugate(pvector8M[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::pGap8P(int n, int p)
{

	if(n>=0) return pvector8P[n][p];
	else return TComplex::Conjugate(pvector8P[n][p]);

}
//____________________________________________________________________
void AliAnalysisTaskChargedFlow::ResetQ(const int nMaxHarm, const int nMaxPow)
{

	for(int i=0; i<nMaxHarm; i++)
	{
		for(int j=0; j<nMaxPow; j++)
		{
			Qvector[i][j] = TComplex(0.,0.);
		}
	}

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two(int n1, int n2)
{
	return Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
}

TComplex AliAnalysisTaskChargedFlow::TwoTest(int n1, int n2)
{
	return Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap(int n1, int n2)
{
	return QGapM(n1,1)*QGapP(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap0(int n1, int n2)
{
	return QGap0M(n1,1)*QGap0P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap4(int n1, int n2)
{
	return QGap4M(n1,1)*QGap4P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap8(int n1, int n2)
{
	return QGap8M(n1,1)*QGap8P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiff(int n1, int n2)
{
	return p(n1,1)*Q(n2,1) - q(n1+n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGapQP(int n1, int n2)
{
	return pGapM(n1,1)*QGapP(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGapQM(int n1, int n2)
{
	return pGapP(n1,1)*QGapM(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGap0QP(int n1, int n2)
{
	return pGap0M(n1,1)*QGap0P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGap0QM(int n1, int n2)
{
	return pGap0P(n1,1)*QGap0M(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGap4QP(int n1, int n2)
{
	return pGap4M(n1,1)*QGap4P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGap4QM(int n1, int n2)
{
	return pGap4P(n1,1)*QGap4M(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGap8QP(int n1, int n2)
{
	return pGap8M(n1,1)*QGap8P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffGap8QM(int n1, int n2)
{
	return pGap8P(n1,1)*QGap8M(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffPtGap0(int n1, int n2)
{
	return pGap0M(n1,1)*pGap0P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffPtGap4(int n1, int n2)
{
	return pGap4M(n1,1)*pGap4P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffPtGap8(int n1, int n2)
{
	return pGap8M(n1,1)*pGap8P(n2,1);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoDiffPtGap10(int n1, int n2)
{
	return pGapM(n1,1)*pGapP(n2,1);
}
//____________________________________________________________________
/*
TComplex AliAnalysisTaskChargedFlow::Three(int n1, int n2, int n3)
{
	return
	(
		Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
	  - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3)	
	); 
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::ThreeDiff(int n1, int n2, int n3)
{
	return 	
	(
		p(n1,1)*Q(n2,1)*Q(n3,1)-q(n1+n2,2)*Q(n3,1)-q(n1+n3,2)*Q(n2,1)
    - p(n1,1)*Q(n2+n3,2)+2.*q(n1+n2+n3,3)
	);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Four(int n1, int n2, int n3, int n4)
{
	return	
	( 
		Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
	  + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
	  + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
	  + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4)  
  );
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::FourDiff(int n1, int n2, int n3, int n4)
{
	return
	(
		p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*q(n1+n3,2)*Q(n4,1)
 		- p(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*q(n1+n4,2)
 		+ Q(n2+n3,2)*q(n1+n4,2)-p(n1,1)*Q(n3,1)*Q(n2+n4,2)+q(n1+n3,2)*Q(n2+n4,2)
 		+ 2.*Q(n3,1)*q(n1+n2+n4,3)-p(n1,1)*Q(n2,1)*Q(n3+n4,2)+q(n1+n2,2)*Q(n3+n4,2)
 		+ 2.*Q(n2,1)*q(n1+n3+n4,3)+2.*p(n1,1)*Q(n2+n3+n4,3)-6.*q(n1+n2+n3+n4,4) 
 	);
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::FourPtDiff(int n1, int n2, int n3, int n4)
{
	return
	(
		p(n1,1)*p(n2,1)*p(n3,1)*p(n4,1)-p(n1+n2,2)*p(n3,1)*p(n4,1)-p(n2,1)*p(n1+n3,2)*p(n4,1)
 		- p(n1,1)*p(n2+n3,2)*p(n4,1)+2.*p(n1+n2+n3,3)*p(n4,1)-p(n2,1)*p(n3,1)*p(n1+n4,2)
 		+ p(n2+n3,2)*p(n1+n4,2)-p(n1,1)*p(n3,1)*p(n2+n4,2)+p(n1+n3,2)*p(n2+n4,2)
 		+ 2.*p(n3,1)*p(n1+n2+n4,3)-p(n1,1)*p(n2,1)*p(n3+n4,2)+p(n1+n2,2)*p(n3+n4,2)
 		+ 2.*p(n2,1)*p(n1+n3+n4,3)+2.*p(n1,1)*p(n2+n3+n4,3)-6.*p(n1+n2+n3+n4,4)
 	);
}
//___________________________________________________________________ 
TComplex AliAnalysisTaskChargedFlow::Five(int n1, int n2, int n3, int n4, int n5)
{
	return
	(
		Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
		+ 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
		+ Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
		+ Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
		- Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
		- 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
		+ Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
		+ Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
		- Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
		+ Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
		- 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
		- 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
		+ Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
		+ Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
		+ 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
		+ 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
		- 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
		+ Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
		+ Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
		+ 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
		+ 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
		- 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
		- 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5)
  ); 
}
//___________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{
	return
	(
		Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
    - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
    - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
    + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
    + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
    - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
    - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
    - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
    - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
    + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
    + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
    - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
    - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
    - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
    - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
    + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
    - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
    - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
    - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
    - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
    + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
    - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
    - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
    - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
    + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
    - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
    + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
    - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
    - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
    - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
    + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
    + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
    - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
    - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
    - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
    - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
    + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
    + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
    - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
    - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
    - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
    - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
    + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
    - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
    - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
    - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
    + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
    + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
    + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
    + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
    + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
    + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
    + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
    - 120.*Q(n1+n2+n3+n4+n5+n6,6)
  );
}
//---------------------------------------------
TComplex AliAnalysisTaskChargedFlow::Eight_1(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	return 
		(
		Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)*Q(n7,1)*Q(n8,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n8,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)*Q(n7,1)*Q(n8,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)*Q(n7,1)*Q(n8,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n8,1)
    + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n8,1)
    + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n8,1)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n8,1)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)*Q(n7,1)*Q(n8,1)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)*Q(n7,1)*Q(n8,1)
    - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)*Q(n7,1)*Q(n8,1)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n8,1)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n8,1)
    + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n8,1)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)*Q(n7,1)*Q(n8,1)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n8,1)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)*Q(n7,1)*Q(n8,1)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)*Q(n7,1)*Q(n8,1)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)*Q(n7,1)*Q(n8,1)
    + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)*Q(n7,1)*Q(n8,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n8,1)
    + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n8,1)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n8,1)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)*Q(n7,1)*Q(n8,1)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n8,1)
    + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n8,1)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)*Q(n7,1)*Q(n8,1)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)*Q(n7,1)*Q(n8,1)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)*Q(n7,1)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n8,1)
    + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n8,1)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)*Q(n7,1)*Q(n8,1)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)*Q(n7,1)*Q(n8,1)
    + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)*Q(n7,1)*Q(n8,1)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)*Q(n7,1)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)*Q(n7,1)*Q(n8,1)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)*Q(n7,1)*Q(n8,1)
    + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)*Q(n7,1)*Q(n8,1)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)*Q(n7,1)*Q(n8,1)
    - 120.*Q(n1+n2+n3+n4+n5+n6,6)*Q(n7,1)*Q(n8,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    + Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)+Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)-2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    + Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)-Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)-Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)+Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    - Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)+6.*Q(n2+n3+n4+n5,4)*Q(n6,1)*Q(n1+n7,2)*Q(n8,1)
    + Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n8,1)
    - Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n3+n4+n5,3)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n8,1)+Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n8,1)
    - Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n8,1)+2.*Q(n2+n4+n5,3)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n1+n7,2)*Q(n8,1)+2.*Q(n4+n5,2)*Q(n2+n3+n6,3)*Q(n1+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n2+n3+n5,3)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n3+n5,2)*Q(n2+n4+n6,3)*Q(n1+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n2+n5,2)*Q(n3+n4+n6,3)*Q(n1+n7,2)*Q(n8,1)+6.*Q(n5,1)*Q(n2+n3+n4+n6,4)*Q(n1+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n8,1)-Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n2+n3+n4,3)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n3+n4,2)*Q(n2+n5+n6,3)*Q(n1+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n1+n7,2)*Q(n8,1)
    + 2.*Q(n2+n4,2)*Q(n3+n5+n6,3)*Q(n1+n7,2)*Q(n8,1)+6.*Q(n4,1)*Q(n2+n3+n5+n6,4)*Q(n1+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n1+n7,2)*Q(n8,1)+2.*Q(n2+n3,2)*Q(n4+n5+n6,3)*Q(n1+n7,2)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n2+n4+n5+n6,4)*Q(n1+n7,2)*Q(n8,1)+6.*Q(n2,1)*Q(n3+n4+n5+n6,4)*Q(n1+n7,2)*Q(n8,1)
    - 24.*Q(n2+n3+n4+n5+n6,5)*Q(n1+n7,2)*Q(n8,1)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    - Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)+6.*Q(n1+n3+n4+n5,4)*Q(n6,1)*Q(n2+n7,2)*Q(n8,1)
    + Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n8,1)
    - Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n3+n4+n5,3)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n8,1)+Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n8,1)
    - Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n8,1)+2.*Q(n1+n4+n5,3)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n2+n7,2)*Q(n8,1)+2.*Q(n4+n5,2)*Q(n1+n3+n6,3)*Q(n2+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n1,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n1+n3+n5,3)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n3+n5,2)*Q(n1+n4+n6,3)*Q(n2+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n1+n5,2)*Q(n3+n4+n6,3)*Q(n2+n7,2)*Q(n8,1)+6.*Q(n5,1)*Q(n1+n3+n4+n6,4)*Q(n2+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n8,1)-Q(n1,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n1+n3+n4,3)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n3+n4,2)*Q(n1+n5+n6,3)*Q(n2+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n2+n7,2)*Q(n8,1)
    + 2.*Q(n1+n4,2)*Q(n3+n5+n6,3)*Q(n2+n7,2)*Q(n8,1)+6.*Q(n4,1)*Q(n1+n3+n5+n6,4)*Q(n2+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n2+n7,2)*Q(n8,1)+2.*Q(n1+n3,2)*Q(n4+n5+n6,3)*Q(n2+n7,2)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n1+n4+n5+n6,4)*Q(n2+n7,2)*Q(n8,1)+6.*Q(n1,1)*Q(n3+n4+n5+n6,4)*Q(n2+n7,2)*Q(n8,1)
    - 24.*Q(n1+n3+n4+n5+n6,5)*Q(n2+n7,2)*Q(n8,1)+2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n8,1)
    - 2.*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n8,1)-2.*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n8,1)+4.*Q(n3+n4+n5,3)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n2+n7,3)*Q(n8,1)+2.*Q(n4+n5,2)*Q(n3+n6,2)*Q(n1+n2+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n2+n7,3)*Q(n8,1)+2.*Q(n3+n5,2)*Q(n4+n6,2)*Q(n1+n2+n7,3)*Q(n8,1)
    + 4.*Q(n5,1)*Q(n3+n4+n6,3)*Q(n1+n2+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n2+n7,3)*Q(n8,1)
    + 2.*Q(n3+n4,2)*Q(n5+n6,2)*Q(n1+n2+n7,3)*Q(n8,1)+4.*Q(n4,1)*Q(n3+n5+n6,3)*Q(n1+n2+n7,3)*Q(n8,1)
    + 4.*Q(n3,1)*Q(n4+n5+n6,3)*Q(n1+n2+n7,3)*Q(n8,1)-12.*Q(n3+n4+n5+n6,4)*Q(n1+n2+n7,3)*Q(n8,1)
    - Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)+Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    - 2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)+Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    - Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)+Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    - Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)-2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)-Q(n1+n2,2)*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)
    + 6.*Q(n1+n2+n4+n5,4)*Q(n6,1)*Q(n3+n7,2)*Q(n8,1)+Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n8,1)-Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n2+n4+n5,3)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n8,1)-Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n8,1)-Q(n1,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n8,1)
    + 2.*Q(n1+n4+n5,3)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n8,1)-2.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n3+n7,2)*Q(n8,1)
    + 2.*Q(n4+n5,2)*Q(n1+n2+n6,3)*Q(n3+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n1+n2,2)*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n8,1)-Q(n2,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n1+n2+n5,3)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n2+n5,2)*Q(n1+n4+n6,3)*Q(n3+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n1+n5,2)*Q(n2+n4+n6,3)*Q(n3+n7,2)*Q(n8,1)
    + 6.*Q(n5,1)*Q(n1+n2+n4+n6,4)*Q(n3+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n1+n2,2)*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n8,1)-Q(n2,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n1+n2+n4,3)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n2+n4,2)*Q(n1+n5+n6,3)*Q(n3+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n3+n7,2)*Q(n8,1)+2.*Q(n1+n4,2)*Q(n2+n5+n6,3)*Q(n3+n7,2)*Q(n8,1)
    + 6.*Q(n4,1)*Q(n1+n2+n5+n6,4)*Q(n3+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n2,1)*Q(n4+n5+n6,3)*Q(n3+n7,2)*Q(n8,1)
    + 2.*Q(n1+n2,2)*Q(n4+n5+n6,3)*Q(n3+n7,2)*Q(n8,1)+6.*Q(n2,1)*Q(n1+n4+n5+n6,4)*Q(n3+n7,2)*Q(n8,1)
    + 6.*Q(n1,1)*Q(n2+n4+n5+n6,4)*Q(n3+n7,2)*Q(n8,1)-24.*Q(n1+n2+n4+n5+n6,5)*Q(n3+n7,2)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n8,1)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n8,1)
    + 4.*Q(n2+n4+n5,3)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n8,1)-2.*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n3+n7,3)*Q(n8,1)
    + 2.*Q(n4+n5,2)*Q(n2+n6,2)*Q(n1+n3+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n3+n7,3)*Q(n8,1)
    + 2.*Q(n2+n5,2)*Q(n4+n6,2)*Q(n1+n3+n7,3)*Q(n8,1)+4.*Q(n5,1)*Q(n2+n4+n6,3)*Q(n1+n3+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n3+n7,3)*Q(n8,1)+2.*Q(n2+n4,2)*Q(n5+n6,2)*Q(n1+n3+n7,3)*Q(n8,1)
    + 4.*Q(n4,1)*Q(n2+n5+n6,3)*Q(n1+n3+n7,3)*Q(n8,1)+4.*Q(n2,1)*Q(n4+n5+n6,3)*Q(n1+n3+n7,3)*Q(n8,1)
    - 12.*Q(n2+n4+n5+n6,4)*Q(n1+n3+n7,3)*Q(n8,1)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n8,1)
    - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n8,1)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n8,1)+4.*Q(n1+n4+n5,3)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n8,1)
    - 2.*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n3+n7,3)*Q(n8,1)+2.*Q(n4+n5,2)*Q(n1+n6,2)*Q(n2+n3+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n3+n7,3)*Q(n8,1)+2.*Q(n1+n5,2)*Q(n4+n6,2)*Q(n2+n3+n7,3)*Q(n8,1)
    + 4.*Q(n5,1)*Q(n1+n4+n6,3)*Q(n2+n3+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n3+n7,3)*Q(n8,1)
    + 2.*Q(n1+n4,2)*Q(n5+n6,2)*Q(n2+n3+n7,3)*Q(n8,1)+4.*Q(n4,1)*Q(n1+n5+n6,3)*Q(n2+n3+n7,3)*Q(n8,1)
    + 4.*Q(n1,1)*Q(n4+n5+n6,3)*Q(n2+n3+n7,3)*Q(n8,1)-12.*Q(n1+n4+n5+n6,4)*Q(n2+n3+n7,3)*Q(n8,1)
    - 6.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n3+n7,4)*Q(n8,1)+6.*Q(n4+n5,2)*Q(n6,1)*Q(n1+n2+n3+n7,4)*Q(n8,1)
    + 6.*Q(n5,1)*Q(n4+n6,2)*Q(n1+n2+n3+n7,4)*Q(n8,1)+6.*Q(n4,1)*Q(n5+n6,2)*Q(n1+n2+n3+n7,4)*Q(n8,1)
    - 12.*Q(n4+n5+n6,3)*Q(n1+n2+n3+n7,4)*Q(n8,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)+Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)-2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)-Q(n2+n3,2)*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)-Q(n1+n3,2)*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    - Q(n1+n2,2)*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)+6.*Q(n1+n2+n3+n5,4)*Q(n6,1)*Q(n4+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n2+n3,2)*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n2,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n2+n3+n5,3)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n8,1)+Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n8,1)
    - Q(n1+n3,2)*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n3,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n8,1)+2.*Q(n1+n3+n5,3)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n4+n7,2)*Q(n8,1)+2.*Q(n3+n5,2)*Q(n1+n2+n6,3)*Q(n4+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n1+n2,2)*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n1,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n1+n2+n5,3)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n2+n5,2)*Q(n1+n3+n6,3)*Q(n4+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n1+n5,2)*Q(n2+n3+n6,3)*Q(n4+n7,2)*Q(n8,1)+6.*Q(n5,1)*Q(n1+n2+n3+n6,4)*Q(n4+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n1+n2,2)*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n8,1)-Q(n1,1)*Q(n2+n3,2)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n1+n2+n3,3)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n5+n6,3)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n2+n3,2)*Q(n1+n5+n6,3)*Q(n4+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n5+n6,3)*Q(n4+n7,2)*Q(n8,1)
    + 2.*Q(n1+n3,2)*Q(n2+n5+n6,3)*Q(n4+n7,2)*Q(n8,1)+6.*Q(n3,1)*Q(n1+n2+n5+n6,4)*Q(n4+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n5+n6,3)*Q(n4+n7,2)*Q(n8,1)+2.*Q(n1+n2,2)*Q(n3+n5+n6,3)*Q(n4+n7,2)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n1+n3+n5+n6,4)*Q(n4+n7,2)*Q(n8,1)+6.*Q(n1,1)*Q(n2+n3+n5+n6,4)*Q(n4+n7,2)*Q(n8,1)
    - 24.*Q(n1+n2+n3+n5+n6,5)*Q(n4+n7,2)*Q(n8,1)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n8,1)
    - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n8,1)+4.*Q(n2+n3+n5,3)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n4+n7,3)*Q(n8,1)+2.*Q(n3+n5,2)*Q(n2+n6,2)*Q(n1+n4+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n4+n7,3)*Q(n8,1)+2.*Q(n2+n5,2)*Q(n3+n6,2)*Q(n1+n4+n7,3)*Q(n8,1)
    + 4.*Q(n5,1)*Q(n2+n3+n6,3)*Q(n1+n4+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n1+n4+n7,3)*Q(n8,1)
    + 2.*Q(n2+n3,2)*Q(n5+n6,2)*Q(n1+n4+n7,3)*Q(n8,1)+4.*Q(n3,1)*Q(n2+n5+n6,3)*Q(n1+n4+n7,3)*Q(n8,1)
    + 4.*Q(n2,1)*Q(n3+n5+n6,3)*Q(n1+n4+n7,3)*Q(n8,1)-12.*Q(n2+n3+n5+n6,4)*Q(n1+n4+n7,3)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n8,1)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n8,1)
    + 4.*Q(n1+n3+n5,3)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n4+n7,3)*Q(n8,1)
    + 2.*Q(n3+n5,2)*Q(n1+n6,2)*Q(n2+n4+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n4+n7,3)*Q(n8,1)
    + 2.*Q(n1+n5,2)*Q(n3+n6,2)*Q(n2+n4+n7,3)*Q(n8,1)+4.*Q(n5,1)*Q(n1+n3+n6,3)*Q(n2+n4+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3,1)*Q(n5+n6,2)*Q(n2+n4+n7,3)*Q(n8,1)+2.*Q(n1+n3,2)*Q(n5+n6,2)*Q(n2+n4+n7,3)*Q(n8,1)
    + 4.*Q(n3,1)*Q(n1+n5+n6,3)*Q(n2+n4+n7,3)*Q(n8,1)+4.*Q(n1,1)*Q(n3+n5+n6,3)*Q(n2+n4+n7,3)*Q(n8,1)
    - 12.*Q(n1+n3+n5+n6,4)*Q(n2+n4+n7,3)*Q(n8,1)-6.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n4+n7,4)*Q(n8,1)
    + 6.*Q(n3+n5,2)*Q(n6,1)*Q(n1+n2+n4+n7,4)*Q(n8,1)+6.*Q(n5,1)*Q(n3+n6,2)*Q(n1+n2+n4+n7,4)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n5+n6,2)*Q(n1+n2+n4+n7,4)*Q(n8,1)-12.*Q(n3+n5+n6,3)*Q(n1+n2+n4+n7,4)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n8,1)-2.*Q(n1+n2,2)*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n8,1)
    + 4.*Q(n1+n2+n5,3)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n4+n7,3)*Q(n8,1)
    + 2.*Q(n2+n5,2)*Q(n1+n6,2)*Q(n3+n4+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n4+n7,3)*Q(n8,1)
    + 2.*Q(n1+n5,2)*Q(n2+n6,2)*Q(n3+n4+n7,3)*Q(n8,1)+4.*Q(n5,1)*Q(n1+n2+n6,3)*Q(n3+n4+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2,1)*Q(n5+n6,2)*Q(n3+n4+n7,3)*Q(n8,1)+2.*Q(n1+n2,2)*Q(n5+n6,2)*Q(n3+n4+n7,3)*Q(n8,1)
    + 4.*Q(n2,1)*Q(n1+n5+n6,3)*Q(n3+n4+n7,3)*Q(n8,1)+4.*Q(n1,1)*Q(n2+n5+n6,3)*Q(n3+n4+n7,3)*Q(n8,1)
    - 12.*Q(n1+n2+n5+n6,4)*Q(n3+n4+n7,3)*Q(n8,1)-6.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n4+n7,4)*Q(n8,1)
    + 6.*Q(n2+n5,2)*Q(n6,1)*Q(n1+n3+n4+n7,4)*Q(n8,1)+6.*Q(n5,1)*Q(n2+n6,2)*Q(n1+n3+n4+n7,4)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n5+n6,2)*Q(n1+n3+n4+n7,4)*Q(n8,1)-12.*Q(n2+n5+n6,3)*Q(n1+n3+n4+n7,4)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n4+n7,4)*Q(n8,1)+6.*Q(n1+n5,2)*Q(n6,1)*Q(n2+n3+n4+n7,4)*Q(n8,1)
    + 6.*Q(n5,1)*Q(n1+n6,2)*Q(n2+n3+n4+n7,4)*Q(n8,1)+6.*Q(n1,1)*Q(n5+n6,2)*Q(n2+n3+n4+n7,4)*Q(n8,1)
    - 12.*Q(n1+n5+n6,3)*Q(n2+n3+n4+n7,4)*Q(n8,1)+24.*Q(n5,1)*Q(n6,1)*Q(n1+n2+n3+n4+n7,5)*Q(n8,1)
    - 24.*Q(n5+n6,2)*Q(n1+n2+n3+n4+n7,5)*Q(n8,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    - Q(n1+n2,2)*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)+6.*Q(n1+n2+n3+n4,4)*Q(n6,1)*Q(n5+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n2+n3,2)*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n2+n4,2)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n2,1)*Q(n3+n4,2)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n2+n3+n4,3)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n8,1)+Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n8,1)
    - Q(n1+n3,2)*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n3,1)*Q(n1+n4,2)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n3+n4,2)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n8,1)+2.*Q(n1+n3+n4,3)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n6,3)*Q(n5+n7,2)*Q(n8,1)+2.*Q(n3+n4,2)*Q(n1+n2+n6,3)*Q(n5+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n1+n2,2)*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n1+n4,2)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n1,1)*Q(n2+n4,2)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n1+n2+n4,3)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n6,3)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n2+n4,2)*Q(n1+n3+n6,3)*Q(n5+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n6,3)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n1+n4,2)*Q(n2+n3+n6,3)*Q(n5+n7,2)*Q(n8,1)+6.*Q(n4,1)*Q(n1+n2+n3+n6,4)*Q(n5+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n8,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n8,1)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n6,3)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n2+n3,2)*Q(n1+n4+n6,3)*Q(n5+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n6,3)*Q(n5+n7,2)*Q(n8,1)
    + 2.*Q(n1+n3,2)*Q(n2+n4+n6,3)*Q(n5+n7,2)*Q(n8,1)+6.*Q(n3,1)*Q(n1+n2+n4+n6,4)*Q(n5+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n6,3)*Q(n5+n7,2)*Q(n8,1)+2.*Q(n1+n2,2)*Q(n3+n4+n6,3)*Q(n5+n7,2)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n1+n3+n4+n6,4)*Q(n5+n7,2)*Q(n8,1)+6.*Q(n1,1)*Q(n2+n3+n4+n6,4)*Q(n5+n7,2)*Q(n8,1)
    - 24.*Q(n1+n2+n3+n4+n6,5)*Q(n5+n7,2)*Q(n8,1)+2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n8,1)
    - 2.*Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n8,1)+4.*Q(n2+n3+n4,3)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n1+n5+n7,3)*Q(n8,1)+2.*Q(n3+n4,2)*Q(n2+n6,2)*Q(n1+n5+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n1+n5+n7,3)*Q(n8,1)+2.*Q(n2+n4,2)*Q(n3+n6,2)*Q(n1+n5+n7,3)*Q(n8,1)
    + 4.*Q(n4,1)*Q(n2+n3+n6,3)*Q(n1+n5+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n1+n5+n7,3)*Q(n8,1)
    + 2.*Q(n2+n3,2)*Q(n4+n6,2)*Q(n1+n5+n7,3)*Q(n8,1)+4.*Q(n3,1)*Q(n2+n4+n6,3)*Q(n1+n5+n7,3)*Q(n8,1)
    + 4.*Q(n2,1)*Q(n3+n4+n6,3)*Q(n1+n5+n7,3)*Q(n8,1)-12.*Q(n2+n3+n4+n6,4)*Q(n1+n5+n7,3)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n8,1)-2.*Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n3+n4,2)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n8,1)
    + 4.*Q(n1+n3+n4,3)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n2+n5+n7,3)*Q(n8,1)
    + 2.*Q(n3+n4,2)*Q(n1+n6,2)*Q(n2+n5+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n4,1)*Q(n3+n6,2)*Q(n2+n5+n7,3)*Q(n8,1)
    + 2.*Q(n1+n4,2)*Q(n3+n6,2)*Q(n2+n5+n7,3)*Q(n8,1)+4.*Q(n4,1)*Q(n1+n3+n6,3)*Q(n2+n5+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3,1)*Q(n4+n6,2)*Q(n2+n5+n7,3)*Q(n8,1)+2.*Q(n1+n3,2)*Q(n4+n6,2)*Q(n2+n5+n7,3)*Q(n8,1)
    + 4.*Q(n3,1)*Q(n1+n4+n6,3)*Q(n2+n5+n7,3)*Q(n8,1)+4.*Q(n1,1)*Q(n3+n4+n6,3)*Q(n2+n5+n7,3)*Q(n8,1)
    - 12.*Q(n1+n3+n4+n6,4)*Q(n2+n5+n7,3)*Q(n8,1)-6.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n2+n5+n7,4)*Q(n8,1)
    + 6.*Q(n3+n4,2)*Q(n6,1)*Q(n1+n2+n5+n7,4)*Q(n8,1)+6.*Q(n4,1)*Q(n3+n6,2)*Q(n1+n2+n5+n7,4)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n4+n6,2)*Q(n1+n2+n5+n7,4)*Q(n8,1)-12.*Q(n3+n4+n6,3)*Q(n1+n2+n5+n7,4)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n8,1)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n8,1)
    + 4.*Q(n1+n2+n4,3)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n4,1)*Q(n1+n6,2)*Q(n3+n5+n7,3)*Q(n8,1)
    + 2.*Q(n2+n4,2)*Q(n1+n6,2)*Q(n3+n5+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n4,1)*Q(n2+n6,2)*Q(n3+n5+n7,3)*Q(n8,1)
    + 2.*Q(n1+n4,2)*Q(n2+n6,2)*Q(n3+n5+n7,3)*Q(n8,1)+4.*Q(n4,1)*Q(n1+n2+n6,3)*Q(n3+n5+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2,1)*Q(n4+n6,2)*Q(n3+n5+n7,3)*Q(n8,1)+2.*Q(n1+n2,2)*Q(n4+n6,2)*Q(n3+n5+n7,3)*Q(n8,1)
    + 4.*Q(n2,1)*Q(n1+n4+n6,3)*Q(n3+n5+n7,3)*Q(n8,1)+4.*Q(n1,1)*Q(n2+n4+n6,3)*Q(n3+n5+n7,3)*Q(n8,1)
    - 12.*Q(n1+n2+n4+n6,4)*Q(n3+n5+n7,3)*Q(n8,1)-6.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n1+n3+n5+n7,4)*Q(n8,1)
    + 6.*Q(n2+n4,2)*Q(n6,1)*Q(n1+n3+n5+n7,4)*Q(n8,1)+6.*Q(n4,1)*Q(n2+n6,2)*Q(n1+n3+n5+n7,4)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n4+n6,2)*Q(n1+n3+n5+n7,4)*Q(n8,1)-12.*Q(n2+n4+n6,3)*Q(n1+n3+n5+n7,4)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n2+n3+n5+n7,4)*Q(n8,1)+6.*Q(n1+n4,2)*Q(n6,1)*Q(n2+n3+n5+n7,4)*Q(n8,1)
    + 6.*Q(n4,1)*Q(n1+n6,2)*Q(n2+n3+n5+n7,4)*Q(n8,1)+6.*Q(n1,1)*Q(n4+n6,2)*Q(n2+n3+n5+n7,4)*Q(n8,1)
    - 12.*Q(n1+n4+n6,3)*Q(n2+n3+n5+n7,4)*Q(n8,1)+24.*Q(n4,1)*Q(n6,1)*Q(n1+n2+n3+n5+n7,5)*Q(n8,1)
    - 24.*Q(n4+n6,2)*Q(n1+n2+n3+n5+n7,5)*Q(n8,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n8,1)
    - 2.*Q(n1+n2,2)*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n3,2)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n3,2)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n8,1)+4.*Q(n1+n2+n3,3)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3,1)*Q(n1+n6,2)*Q(n4+n5+n7,3)*Q(n8,1)+2.*Q(n2+n3,2)*Q(n1+n6,2)*Q(n4+n5+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3,1)*Q(n2+n6,2)*Q(n4+n5+n7,3)*Q(n8,1)+2.*Q(n1+n3,2)*Q(n2+n6,2)*Q(n4+n5+n7,3)*Q(n8,1)
    + 4.*Q(n3,1)*Q(n1+n2+n6,3)*Q(n4+n5+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n2,1)*Q(n3+n6,2)*Q(n4+n5+n7,3)*Q(n8,1)
    + 2.*Q(n1+n2,2)*Q(n3+n6,2)*Q(n4+n5+n7,3)*Q(n8,1)+4.*Q(n2,1)*Q(n1+n3+n6,3)*Q(n4+n5+n7,3)*Q(n8,1)
    + 4.*Q(n1,1)*Q(n2+n3+n6,3)*Q(n4+n5+n7,3)*Q(n8,1)-12.*Q(n1+n2+n3+n6,4)*Q(n4+n5+n7,3)*Q(n8,1)
    - 6.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n1+n4+n5+n7,4)*Q(n8,1)+6.*Q(n2+n3,2)*Q(n6,1)*Q(n1+n4+n5+n7,4)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n2+n6,2)*Q(n1+n4+n5+n7,4)*Q(n8,1)+6.*Q(n2,1)*Q(n3+n6,2)*Q(n1+n4+n5+n7,4)*Q(n8,1)
    - 12.*Q(n2+n3+n6,3)*Q(n1+n4+n5+n7,4)*Q(n8,1)-6.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n2+n4+n5+n7,4)*Q(n8,1)
    + 6.*Q(n1+n3,2)*Q(n6,1)*Q(n2+n4+n5+n7,4)*Q(n8,1)+6.*Q(n3,1)*Q(n1+n6,2)*Q(n2+n4+n5+n7,4)*Q(n8,1)
    + 6.*Q(n1,1)*Q(n3+n6,2)*Q(n2+n4+n5+n7,4)*Q(n8,1)-12.*Q(n1+n3+n6,3)*Q(n2+n4+n5+n7,4)*Q(n8,1)
    + 24.*Q(n3,1)*Q(n6,1)*Q(n1+n2+n4+n5+n7,5)*Q(n8,1)-24.*Q(n3+n6,2)*Q(n1+n2+n4+n5+n7,5)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n3+n4+n5+n7,4)*Q(n8,1)+6.*Q(n1+n2,2)*Q(n6,1)*Q(n3+n4+n5+n7,4)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n1+n6,2)*Q(n3+n4+n5+n7,4)*Q(n8,1)+6.*Q(n1,1)*Q(n2+n6,2)*Q(n3+n4+n5+n7,4)*Q(n8,1)
    - 12.*Q(n1+n2+n6,3)*Q(n3+n4+n5+n7,4)*Q(n8,1)+24.*Q(n2,1)*Q(n6,1)*Q(n1+n3+n4+n5+n7,5)*Q(n8,1)
    - 24.*Q(n2+n6,2)*Q(n1+n3+n4+n5+n7,5)*Q(n8,1)+24.*Q(n1,1)*Q(n6,1)*Q(n2+n3+n4+n5+n7,5)*Q(n8,1)
    - 24.*Q(n1+n6,2)*Q(n2+n3+n4+n5+n7,5)*Q(n8,1)-120.*Q(n6,1)*Q(n1+n2+n3+n4+n5+n7,6)*Q(n8,1)
    - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)+Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    + Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    - 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)+Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    - Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)+Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    - Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)-Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)
    + 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6+n7,2)*Q(n8,1)+Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n8,1)-Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n8,1)
    + Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n8,1)-Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n8,1)-Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n8,1)
    + 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n8,1)-2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6+n7,2)*Q(n8,1)
    + 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n8,1)-Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6+n7,2)*Q(n8,1)
    + 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6+n7,2)*Q(n8,1)+Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n8,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6+n7,2)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6+n7,2)*Q(n8,1)+2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6+n7,2)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6+n7,2)*Q(n8,1)-2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6+n7,2)*Q(n8,1)
    + 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6+n7,2)*Q(n8,1)+6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6+n7,2)*Q(n8,1)
    + 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6+n7,2)*Q(n8,1)-24.*Q(n1+n2+n3+n4+n5,5)*Q(n6+n7,2)*Q(n8,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n8,1)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n8,1)
    + 4.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6+n7,3)*Q(n8,1)
    + 2.*Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6+n7,3)*Q(n8,1)
    + 2.*Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6+n7,3)*Q(n8,1)+4.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6+n7,3)*Q(n8,1)+2.*Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6+n7,3)*Q(n8,1)
    + 4.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6+n7,3)*Q(n8,1)+4.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6+n7,3)*Q(n8,1)
    - 12.*Q(n2+n3+n4+n5,4)*Q(n1+n6+n7,3)*Q(n8,1)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n8,1)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n8,1)+4.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n8,1)
    - 2.*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6+n7,3)*Q(n8,1)+2.*Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6+n7,3)*Q(n8,1)+2.*Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6+n7,3)*Q(n8,1)
    + 4.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6+n7,3)*Q(n8,1)
    + 2.*Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6+n7,3)*Q(n8,1)+4.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6+n7,3)*Q(n8,1)
    + 4.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6+n7,3)*Q(n8,1)-12.*Q(n1+n3+n4+n5,4)*Q(n2+n6+n7,3)*Q(n8,1)
    - 6.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6+n7,4)*Q(n8,1)+6.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6+n7,4)*Q(n8,1)
    + 6.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6+n7,4)*Q(n8,1)+6.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6+n7,4)*Q(n8,1)
    - 12.*Q(n3+n4+n5,3)*Q(n1+n2+n6+n7,4)*Q(n8,1)+2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n8,1)+4.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6+n7,3)*Q(n8,1)+2.*Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6+n7,3)*Q(n8,1)+2.*Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6+n7,3)*Q(n8,1)
    + 4.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6+n7,3)*Q(n8,1)
    + 2.*Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6+n7,3)*Q(n8,1)+4.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6+n7,3)*Q(n8,1)
    + 4.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6+n7,3)*Q(n8,1)-12.*Q(n1+n2+n4+n5,4)*Q(n3+n6+n7,3)*Q(n8,1)
    - 6.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6+n7,4)*Q(n8,1)+6.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6+n7,4)*Q(n8,1)
    + 6.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6+n7,4)*Q(n8,1)+6.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6+n7,4)*Q(n8,1)
    - 12.*Q(n2+n4+n5,3)*Q(n1+n3+n6+n7,4)*Q(n8,1)-6.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6+n7,4)*Q(n8,1)
    + 6.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6+n7,4)*Q(n8,1)+6.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6+n7,4)*Q(n8,1)
    + 6.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6+n7,4)*Q(n8,1)-12.*Q(n1+n4+n5,3)*Q(n2+n3+n6+n7,4)*Q(n8,1)
    + 24.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6+n7,5)*Q(n8,1)-24.*Q(n4+n5,2)*Q(n1+n2+n3+n6+n7,5)*Q(n8,1)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n8,1)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n8,1)
    + 4.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6+n7,3)*Q(n8,1)
    + 2.*Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6+n7,3)*Q(n8,1)
    + 2.*Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6+n7,3)*Q(n8,1)+4.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6+n7,3)*Q(n8,1)+2.*Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6+n7,3)*Q(n8,1)
    + 4.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6+n7,3)*Q(n8,1)+4.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6+n7,3)*Q(n8,1)
    - 12.*Q(n1+n2+n3+n5,4)*Q(n4+n6+n7,3)*Q(n8,1)-6.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6+n7,4)*Q(n8,1)
    + 6.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6+n7,4)*Q(n8,1)+6.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6+n7,4)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6+n7,4)*Q(n8,1)-12.*Q(n2+n3+n5,3)*Q(n1+n4+n6+n7,4)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6+n7,4)*Q(n8,1)+6.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6+n7,4)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6+n7,4)*Q(n8,1)+6.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6+n7,4)*Q(n8,1)
    - 12.*Q(n1+n3+n5,3)*Q(n2+n4+n6+n7,4)*Q(n8,1)+24.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6+n7,5)*Q(n8,1)
    - 24.*Q(n3+n5,2)*Q(n1+n2+n4+n6+n7,5)*Q(n8,1)-6.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6+n7,4)*Q(n8,1)
    + 6.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6+n7,4)*Q(n8,1)+6.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6+n7,4)*Q(n8,1)
    + 6.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6+n7,4)*Q(n8,1)-12.*Q(n1+n2+n5,3)*Q(n3+n4+n6+n7,4)*Q(n8,1)
    + 24.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6+n7,5)*Q(n8,1)-24.*Q(n2+n5,2)*Q(n1+n3+n4+n6+n7,5)*Q(n8,1)
    + 24.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6+n7,5)*Q(n8,1)-24.*Q(n1+n5,2)*Q(n2+n3+n4+n6+n7,5)*Q(n8,1)
    - 120.*Q(n5,1)*Q(n1+n2+n3+n4+n6+n7,6)*Q(n8,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n8,1)-2.*Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n8,1)+4.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n8,1)
    - 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6+n7,3)*Q(n8,1)+2.*Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6+n7,3)*Q(n8,1)
    - 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6+n7,3)*Q(n8,1)+2.*Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6+n7,3)*Q(n8,1)
    + 4.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6+n7,3)*Q(n8,1)-2.*Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6+n7,3)*Q(n8,1)
    + 2.*Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6+n7,3)*Q(n8,1)+4.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6+n7,3)*Q(n8,1)
    + 4.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6+n7,3)*Q(n8,1)-12.*Q(n1+n2+n3+n4,4)*Q(n5+n6+n7,3)*Q(n8,1)
    - 6.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6+n7,4)*Q(n8,1)+6.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6+n7,4)*Q(n8,1)
    + 6.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6+n7,4)*Q(n8,1)+6.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6+n7,4)*Q(n8,1)
    - 12.*Q(n2+n3+n4,3)*Q(n1+n5+n6+n7,4)*Q(n8,1)-6.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6+n7,4)*Q(n8,1)
    + 6.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6+n7,4)*Q(n8,1)+6.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6+n7,4)*Q(n8,1)
    + 6.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6+n7,4)*Q(n8,1)-12.*Q(n1+n3+n4,3)*Q(n2+n5+n6+n7,4)*Q(n8,1)
    + 24.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6+n7,5)*Q(n8,1)-24.*Q(n3+n4,2)*Q(n1+n2+n5+n6+n7,5)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6+n7,4)*Q(n8,1)+6.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6+n7,4)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6+n7,4)*Q(n8,1)+6.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6+n7,4)*Q(n8,1)
    - 12.*Q(n1+n2+n4,3)*Q(n3+n5+n6+n7,4)*Q(n8,1)+24.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6+n7,5)*Q(n8,1)
    - 24.*Q(n2+n4,2)*Q(n1+n3+n5+n6+n7,5)*Q(n8,1)+24.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6+n7,5)*Q(n8,1)
    - 24.*Q(n1+n4,2)*Q(n2+n3+n5+n6+n7,5)*Q(n8,1)-120.*Q(n4,1)*Q(n1+n2+n3+n5+n6+n7,6)*Q(n8,1)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6+n7,4)*Q(n8,1)+6.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6+n7,4)*Q(n8,1)
    + 6.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6+n7,4)*Q(n8,1)+6.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6+n7,4)*Q(n8,1)
    - 12.*Q(n1+n2+n3,3)*Q(n4+n5+n6+n7,4)*Q(n8,1)+24.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6+n7,5)*Q(n8,1)
    - 24.*Q(n2+n3,2)*Q(n1+n4+n5+n6+n7,5)*Q(n8,1)+24.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6+n7,5)*Q(n8,1)
    - 24.*Q(n1+n3,2)*Q(n2+n4+n5+n6+n7,5)*Q(n8,1)-120.*Q(n3,1)*Q(n1+n2+n4+n5+n6+n7,6)*Q(n8,1)
    + 24.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6+n7,5)*Q(n8,1)-24.*Q(n1+n2,2)*Q(n3+n4+n5+n6+n7,5)*Q(n8,1)
    - 120.*Q(n2,1)*Q(n1+n3+n4+n5+n6+n7,6)*Q(n8,1)-120.*Q(n1,1)*Q(n2+n3+n4+n5+n6+n7,6)*Q(n8,1)
    + 720.*Q(n1+n2+n3+n4+n5+n6+n7,7)*Q(n8,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    + Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)+Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    + Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)-2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    + Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)-Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    + Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)-Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    - 2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)+Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    - Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)-2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    - 2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)+6.*Q(n2+n3+n4+n5,4)*Q(n6,1)*Q(n7,1)*Q(n1+n8,2)
    + Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n8,2)
    - Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n3+n4+n5,3)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n8,2)+Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n8,2)
    - Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n8,2)
    - Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n8,2)+2.*Q(n2+n4+n5,3)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n8,2)
    - 2.*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n1+n8,2)+2.*Q(n4+n5,2)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n1+n8,2)
    + Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n8,2)
    - Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n2+n3+n5,3)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n8,2)-2.*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n3+n5,2)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n1+n8,2)-2.*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n2+n5,2)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n1+n8,2)+6.*Q(n5,1)*Q(n2+n3+n4+n6,4)*Q(n7,1)*Q(n1+n8,2)
    + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n8,2)
    - Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n8,2)-Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n2+n3+n4,3)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n8,2)-2.*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n3+n4,2)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n1+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n1+n8,2)
    + 2.*Q(n2+n4,2)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n1+n8,2)+6.*Q(n4,1)*Q(n2+n3+n5+n6,4)*Q(n7,1)*Q(n1+n8,2)
    - 2.*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n1+n8,2)+2.*Q(n2+n3,2)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n1+n8,2)
    + 6.*Q(n3,1)*Q(n2+n4+n5+n6,4)*Q(n7,1)*Q(n1+n8,2)+6.*Q(n2,1)*Q(n3+n4+n5+n6,4)*Q(n7,1)*Q(n1+n8,2)
    - 24.*Q(n2+n3+n4+n5+n6,5)*Q(n7,1)*Q(n1+n8,2)+Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n8,2)
    - Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n8,2)-Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n8,2)
    - Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n8,2)+2.*Q(n3+n4+n5,3)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n8,2)
    - Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n1+n8,2)+Q(n4+n5,2)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n1+n8,2)
    - Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n1+n8,2)+Q(n3+n5,2)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n1+n8,2)
    + 2.*Q(n5,1)*Q(n3+n4+n6,3)*Q(n2+n7,2)*Q(n1+n8,2)-Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n1+n8,2)
    + Q(n3+n4,2)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n1+n8,2)+2.*Q(n4,1)*Q(n3+n5+n6,3)*Q(n2+n7,2)*Q(n1+n8,2)
    + 2.*Q(n3,1)*Q(n4+n5+n6,3)*Q(n2+n7,2)*Q(n1+n8,2)-6.*Q(n3+n4+n5+n6,4)*Q(n2+n7,2)*Q(n1+n8,2)
    + Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n8,2)-Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n8,2)
    - Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n8,2)
    + 2.*Q(n2+n4+n5,3)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n8,2)-Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n1+n8,2)
    + Q(n4+n5,2)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n1+n8,2)
    + Q(n2+n5,2)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n1+n8,2)+2.*Q(n5,1)*Q(n2+n4+n6,3)*Q(n3+n7,2)*Q(n1+n8,2)
    - Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n1+n8,2)+Q(n2+n4,2)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n1+n8,2)
    + 2.*Q(n4,1)*Q(n2+n5+n6,3)*Q(n3+n7,2)*Q(n1+n8,2)+2.*Q(n2,1)*Q(n4+n5+n6,3)*Q(n3+n7,2)*Q(n1+n8,2)
    - 6.*Q(n2+n4+n5+n6,4)*Q(n3+n7,2)*Q(n1+n8,2)-2.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n1+n8,2)
    + 2.*Q(n4+n5,2)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n1+n8,2)+2.*Q(n5,1)*Q(n4+n6,2)*Q(n2+n3+n7,3)*Q(n1+n8,2)
    + 2.*Q(n4,1)*Q(n5+n6,2)*Q(n2+n3+n7,3)*Q(n1+n8,2)-4.*Q(n4+n5+n6,3)*Q(n2+n3+n7,3)*Q(n1+n8,2)
    + Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n8,2)-Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n8,2)
    - Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n8,2)
    + 2.*Q(n2+n3+n5,3)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n8,2)-Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n1+n8,2)
    + Q(n3+n5,2)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n1+n8,2)
    + Q(n2+n5,2)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n1+n8,2)+2.*Q(n5,1)*Q(n2+n3+n6,3)*Q(n4+n7,2)*Q(n1+n8,2)
    - Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n1+n8,2)+Q(n2+n3,2)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n1+n8,2)
    + 2.*Q(n3,1)*Q(n2+n5+n6,3)*Q(n4+n7,2)*Q(n1+n8,2)+2.*Q(n2,1)*Q(n3+n5+n6,3)*Q(n4+n7,2)*Q(n1+n8,2)
    - 6.*Q(n2+n3+n5+n6,4)*Q(n4+n7,2)*Q(n1+n8,2)-2.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n1+n8,2)
    + 2.*Q(n3+n5,2)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n1+n8,2)+2.*Q(n5,1)*Q(n3+n6,2)*Q(n2+n4+n7,3)*Q(n1+n8,2)
    + 2.*Q(n3,1)*Q(n5+n6,2)*Q(n2+n4+n7,3)*Q(n1+n8,2)-4.*Q(n3+n5+n6,3)*Q(n2+n4+n7,3)*Q(n1+n8,2)
    - 2.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n1+n8,2)+2.*Q(n2+n5,2)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n1+n8,2)
    + 2.*Q(n5,1)*Q(n2+n6,2)*Q(n3+n4+n7,3)*Q(n1+n8,2)+2.*Q(n2,1)*Q(n5+n6,2)*Q(n3+n4+n7,3)*Q(n1+n8,2)
    - 4.*Q(n2+n5+n6,3)*Q(n3+n4+n7,3)*Q(n1+n8,2)+6.*Q(n5,1)*Q(n6,1)*Q(n2+n3+n4+n7,4)*Q(n1+n8,2)
    - 6.*Q(n5+n6,2)*Q(n2+n3+n4+n7,4)*Q(n1+n8,2)+Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n8,2)
    - Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n8,2)-Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n8,2)
    - Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n8,2)+2.*Q(n2+n3+n4,3)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n8,2)
    - Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n1+n8,2)+Q(n3+n4,2)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n1+n8,2)
    - Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n1+n8,2)+Q(n2+n4,2)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n1+n8,2)
  );
} 
//=======================================================
TComplex AliAnalysisTaskChargedFlow::Eight_2(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	return	
		(
			2.*Q(n4,1)*Q(n2+n3+n6,3)*Q(n5+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n1+n8,2)
      + Q(n2+n3,2)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n1+n8,2)+2.*Q(n3,1)*Q(n2+n4+n6,3)*Q(n5+n7,2)*Q(n1+n8,2)
      + 2.*Q(n2,1)*Q(n3+n4+n6,3)*Q(n5+n7,2)*Q(n1+n8,2)-6.*Q(n2+n3+n4+n6,4)*Q(n5+n7,2)*Q(n1+n8,2)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n1+n8,2)+2.*Q(n3+n4,2)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n1+n8,2)
      + 2.*Q(n4,1)*Q(n3+n6,2)*Q(n2+n5+n7,3)*Q(n1+n8,2)+2.*Q(n3,1)*Q(n4+n6,2)*Q(n2+n5+n7,3)*Q(n1+n8,2)
      - 4.*Q(n3+n4+n6,3)*Q(n2+n5+n7,3)*Q(n1+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n1+n8,2)
      + 2.*Q(n2+n4,2)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n1+n8,2)+2.*Q(n4,1)*Q(n2+n6,2)*Q(n3+n5+n7,3)*Q(n1+n8,2)
      + 2.*Q(n2,1)*Q(n4+n6,2)*Q(n3+n5+n7,3)*Q(n1+n8,2)-4.*Q(n2+n4+n6,3)*Q(n3+n5+n7,3)*Q(n1+n8,2)
      + 6.*Q(n4,1)*Q(n6,1)*Q(n2+n3+n5+n7,4)*Q(n1+n8,2)-6.*Q(n4+n6,2)*Q(n2+n3+n5+n7,4)*Q(n1+n8,2)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n1+n8,2)+2.*Q(n2+n3,2)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n1+n8,2)
      + 2.*Q(n3,1)*Q(n2+n6,2)*Q(n4+n5+n7,3)*Q(n1+n8,2)+2.*Q(n2,1)*Q(n3+n6,2)*Q(n4+n5+n7,3)*Q(n1+n8,2)
      - 4.*Q(n2+n3+n6,3)*Q(n4+n5+n7,3)*Q(n1+n8,2)+6.*Q(n3,1)*Q(n6,1)*Q(n2+n4+n5+n7,4)*Q(n1+n8,2)
      - 6.*Q(n3+n6,2)*Q(n2+n4+n5+n7,4)*Q(n1+n8,2)+6.*Q(n2,1)*Q(n6,1)*Q(n3+n4+n5+n7,4)*Q(n1+n8,2)
      - 6.*Q(n2+n6,2)*Q(n3+n4+n5+n7,4)*Q(n1+n8,2)-24.*Q(n6,1)*Q(n2+n3+n4+n5+n7,5)*Q(n1+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n8,2)-Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n8,2)
      - Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n8,2)
      + 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n8,2)-Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n1+n8,2)
      + Q(n3+n4,2)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n1+n8,2)-Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n1+n8,2)
      + Q(n2+n4,2)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n1+n8,2)+2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6+n7,2)*Q(n1+n8,2)
      - Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n1+n8,2)+Q(n2+n3,2)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n1+n8,2)
      + 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6+n7,2)*Q(n1+n8,2)+2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6+n7,2)*Q(n1+n8,2)
      - 6.*Q(n2+n3+n4+n5,4)*Q(n6+n7,2)*Q(n1+n8,2)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n1+n8,2)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6+n7,3)*Q(n1+n8,2)-4.*Q(n3+n4+n5,3)*Q(n2+n6+n7,3)*Q(n1+n8,2)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n1+n8,2)+2.*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6+n7,3)*Q(n1+n8,2)+2.*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6+n7,3)*Q(n1+n8,2)
      - 4.*Q(n2+n4+n5,3)*Q(n3+n6+n7,3)*Q(n1+n8,2)+6.*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6+n7,4)*Q(n1+n8,2)
      - 6.*Q(n4+n5,2)*Q(n2+n3+n6+n7,4)*Q(n1+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n1+n8,2)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6+n7,3)*Q(n1+n8,2)-4.*Q(n2+n3+n5,3)*Q(n4+n6+n7,3)*Q(n1+n8,2)
      + 6.*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6+n7,4)*Q(n1+n8,2)-6.*Q(n3+n5,2)*Q(n2+n4+n6+n7,4)*Q(n1+n8,2)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6+n7,4)*Q(n1+n8,2)-6.*Q(n2+n5,2)*Q(n3+n4+n6+n7,4)*Q(n1+n8,2)
      - 24.*Q(n5,1)*Q(n2+n3+n4+n6+n7,5)*Q(n1+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n1+n8,2)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6+n7,3)*Q(n1+n8,2)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6+n7,3)*Q(n1+n8,2)-4.*Q(n2+n3+n4,3)*Q(n5+n6+n7,3)*Q(n1+n8,2)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6+n7,4)*Q(n1+n8,2)-6.*Q(n3+n4,2)*Q(n2+n5+n6+n7,4)*Q(n1+n8,2)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6+n7,4)*Q(n1+n8,2)-6.*Q(n2+n4,2)*Q(n3+n5+n6+n7,4)*Q(n1+n8,2)
      - 24.*Q(n4,1)*Q(n2+n3+n5+n6+n7,5)*Q(n1+n8,2)+6.*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6+n7,4)*Q(n1+n8,2)
      - 6.*Q(n2+n3,2)*Q(n4+n5+n6+n7,4)*Q(n1+n8,2)-24.*Q(n3,1)*Q(n2+n4+n5+n6+n7,5)*Q(n1+n8,2)
      - 24.*Q(n2,1)*Q(n3+n4+n5+n6+n7,5)*Q(n1+n8,2)+120.*Q(n2+n3+n4+n5+n6+n7,6)*Q(n1+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)+Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      + Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)+Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      - 2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)+Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      - Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)+Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      - Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)-2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)-Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      - 2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)-2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)
      + 6.*Q(n1+n3+n4+n5,4)*Q(n6,1)*Q(n7,1)*Q(n2+n8,2)+Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n8,2)-Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n3+n4+n5,3)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n8,2)
      + Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n8,2)-Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n8,2)-Q(n1,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n8,2)
      + 2.*Q(n1+n4+n5,3)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n8,2)-2.*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n2+n8,2)
      + 2.*Q(n4+n5,2)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n2+n8,2)+Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n8,2)-Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n1,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n1+n3+n5,3)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - 2.*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n3+n5,2)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n2+n8,2)
      - 2.*Q(n1,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n1+n5,2)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n2+n8,2)
      + 6.*Q(n5,1)*Q(n1+n3+n4+n6,4)*Q(n7,1)*Q(n2+n8,2)+Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n8,2)-Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - Q(n1,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n1+n3+n4,3)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n8,2)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n3+n4,2)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n2+n8,2)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n2+n8,2)+2.*Q(n1+n4,2)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n2+n8,2)
      + 6.*Q(n4,1)*Q(n1+n3+n5+n6,4)*Q(n7,1)*Q(n2+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n2+n8,2)
      + 2.*Q(n1+n3,2)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n2+n8,2)+6.*Q(n3,1)*Q(n1+n4+n5+n6,4)*Q(n7,1)*Q(n2+n8,2)
      + 6.*Q(n1,1)*Q(n3+n4+n5+n6,4)*Q(n7,1)*Q(n2+n8,2)-24.*Q(n1+n3+n4+n5+n6,5)*Q(n7,1)*Q(n2+n8,2)
      + Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n8,2)-Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n8,2)
      - Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n8,2)-Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n8,2)
      + 2.*Q(n3+n4+n5,3)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n8,2)-Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n2+n8,2)
      + Q(n4+n5,2)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n2+n8,2)-Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n2+n8,2)
      + Q(n3+n5,2)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n2+n8,2)+2.*Q(n5,1)*Q(n3+n4+n6,3)*Q(n1+n7,2)*Q(n2+n8,2)
      - Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n2+n8,2)+Q(n3+n4,2)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n2+n8,2)
      + 2.*Q(n4,1)*Q(n3+n5+n6,3)*Q(n1+n7,2)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n4+n5+n6,3)*Q(n1+n7,2)*Q(n2+n8,2)
      - 6.*Q(n3+n4+n5+n6,4)*Q(n1+n7,2)*Q(n2+n8,2)+Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n8,2)
      - Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n8,2)-Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n8,2)+2.*Q(n1+n4+n5,3)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n8,2)
      - Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n2+n8,2)+Q(n4+n5,2)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n2+n8,2)+Q(n1+n5,2)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n2+n8,2)
      + 2.*Q(n5,1)*Q(n1+n4+n6,3)*Q(n3+n7,2)*Q(n2+n8,2)-Q(n1,1)*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n2+n8,2)
      + Q(n1+n4,2)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n2+n8,2)+2.*Q(n4,1)*Q(n1+n5+n6,3)*Q(n3+n7,2)*Q(n2+n8,2)
      + 2.*Q(n1,1)*Q(n4+n5+n6,3)*Q(n3+n7,2)*Q(n2+n8,2)-6.*Q(n1+n4+n5+n6,4)*Q(n3+n7,2)*Q(n2+n8,2)
      - 2.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n2+n8,2)+2.*Q(n4+n5,2)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n2+n8,2)
      + 2.*Q(n5,1)*Q(n4+n6,2)*Q(n1+n3+n7,3)*Q(n2+n8,2)+2.*Q(n4,1)*Q(n5+n6,2)*Q(n1+n3+n7,3)*Q(n2+n8,2)
      - 4.*Q(n4+n5+n6,3)*Q(n1+n3+n7,3)*Q(n2+n8,2)+Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n8,2)
      - Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n8,2)-Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n8,2)+2.*Q(n1+n3+n5,3)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n8,2)
      - Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n2+n8,2)+Q(n3+n5,2)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n2+n8,2)+Q(n1+n5,2)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n2+n8,2)
      + 2.*Q(n5,1)*Q(n1+n3+n6,3)*Q(n4+n7,2)*Q(n2+n8,2)-Q(n1,1)*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n2+n8,2)
      + Q(n1+n3,2)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n1+n5+n6,3)*Q(n4+n7,2)*Q(n2+n8,2)
      + 2.*Q(n1,1)*Q(n3+n5+n6,3)*Q(n4+n7,2)*Q(n2+n8,2)-6.*Q(n1+n3+n5+n6,4)*Q(n4+n7,2)*Q(n2+n8,2)
      - 2.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n2+n8,2)+2.*Q(n3+n5,2)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n2+n8,2)
      + 2.*Q(n5,1)*Q(n3+n6,2)*Q(n1+n4+n7,3)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n5+n6,2)*Q(n1+n4+n7,3)*Q(n2+n8,2)
      - 4.*Q(n3+n5+n6,3)*Q(n1+n4+n7,3)*Q(n2+n8,2)-2.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n2+n8,2)
      + 2.*Q(n1+n5,2)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n2+n8,2)+2.*Q(n5,1)*Q(n1+n6,2)*Q(n3+n4+n7,3)*Q(n2+n8,2)
      + 2.*Q(n1,1)*Q(n5+n6,2)*Q(n3+n4+n7,3)*Q(n2+n8,2)-4.*Q(n1+n5+n6,3)*Q(n3+n4+n7,3)*Q(n2+n8,2)
      + 6.*Q(n5,1)*Q(n6,1)*Q(n1+n3+n4+n7,4)*Q(n2+n8,2)-6.*Q(n5+n6,2)*Q(n1+n3+n4+n7,4)*Q(n2+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n8,2)-Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n8,2)
      - Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n8,2)-Q(n1,1)*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n8,2)
      + 2.*Q(n1+n3+n4,3)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n8,2)-Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n2+n8,2)
      + Q(n3+n4,2)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n2+n8,2)-Q(n1,1)*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n2+n8,2)
      + Q(n1+n4,2)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n2+n8,2)+2.*Q(n4,1)*Q(n1+n3+n6,3)*Q(n5+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n2+n8,2)+Q(n1+n3,2)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n2+n8,2)
      + 2.*Q(n3,1)*Q(n1+n4+n6,3)*Q(n5+n7,2)*Q(n2+n8,2)+2.*Q(n1,1)*Q(n3+n4+n6,3)*Q(n5+n7,2)*Q(n2+n8,2)
      - 6.*Q(n1+n3+n4+n6,4)*Q(n5+n7,2)*Q(n2+n8,2)-2.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n2+n8,2)
      + 2.*Q(n3+n4,2)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n2+n8,2)+2.*Q(n4,1)*Q(n3+n6,2)*Q(n1+n5+n7,3)*Q(n2+n8,2)
      + 2.*Q(n3,1)*Q(n4+n6,2)*Q(n1+n5+n7,3)*Q(n2+n8,2)-4.*Q(n3+n4+n6,3)*Q(n1+n5+n7,3)*Q(n2+n8,2)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n2+n8,2)+2.*Q(n1+n4,2)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n2+n8,2)
      + 2.*Q(n4,1)*Q(n1+n6,2)*Q(n3+n5+n7,3)*Q(n2+n8,2)+2.*Q(n1,1)*Q(n4+n6,2)*Q(n3+n5+n7,3)*Q(n2+n8,2)
      - 4.*Q(n1+n4+n6,3)*Q(n3+n5+n7,3)*Q(n2+n8,2)+6.*Q(n4,1)*Q(n6,1)*Q(n1+n3+n5+n7,4)*Q(n2+n8,2)
      - 6.*Q(n4+n6,2)*Q(n1+n3+n5+n7,4)*Q(n2+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n2+n8,2)
      + 2.*Q(n1+n3,2)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n1+n6,2)*Q(n4+n5+n7,3)*Q(n2+n8,2)
      + 2.*Q(n1,1)*Q(n3+n6,2)*Q(n4+n5+n7,3)*Q(n2+n8,2)-4.*Q(n1+n3+n6,3)*Q(n4+n5+n7,3)*Q(n2+n8,2)
      + 6.*Q(n3,1)*Q(n6,1)*Q(n1+n4+n5+n7,4)*Q(n2+n8,2)-6.*Q(n3+n6,2)*Q(n1+n4+n5+n7,4)*Q(n2+n8,2)
      + 6.*Q(n1,1)*Q(n6,1)*Q(n3+n4+n5+n7,4)*Q(n2+n8,2)-6.*Q(n1+n6,2)*Q(n3+n4+n5+n7,4)*Q(n2+n8,2)
      - 24.*Q(n6,1)*Q(n1+n3+n4+n5+n7,5)*Q(n2+n8,2)+Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n8,2)
      - Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n8,2)-Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n8,2)+2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n8,2)
      - Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n2+n8,2)+Q(n3+n4,2)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n2+n8,2)
      - Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n2+n8,2)+Q(n1+n4,2)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n2+n8,2)
      + 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6+n7,2)*Q(n2+n8,2)-Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n2+n8,2)
      + Q(n1+n3,2)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6+n7,2)*Q(n2+n8,2)
      + 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n6+n7,2)*Q(n2+n8,2)-6.*Q(n1+n3+n4+n5,4)*Q(n6+n7,2)*Q(n2+n8,2)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n2+n8,2)+2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n2+n8,2)
      + 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6+n7,3)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6+n7,3)*Q(n2+n8,2)
      - 4.*Q(n3+n4+n5,3)*Q(n1+n6+n7,3)*Q(n2+n8,2)-2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n2+n8,2)
      + 2.*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n2+n8,2)+2.*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6+n7,3)*Q(n2+n8,2)
      + 2.*Q(n1,1)*Q(n4+n5,2)*Q(n3+n6+n7,3)*Q(n2+n8,2)-4.*Q(n1+n4+n5,3)*Q(n3+n6+n7,3)*Q(n2+n8,2)
      + 6.*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6+n7,4)*Q(n2+n8,2)-6.*Q(n4+n5,2)*Q(n1+n3+n6+n7,4)*Q(n2+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n2+n8,2)+2.*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n2+n8,2)
      + 2.*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6+n7,3)*Q(n2+n8,2)+2.*Q(n1,1)*Q(n3+n5,2)*Q(n4+n6+n7,3)*Q(n2+n8,2)
      - 4.*Q(n1+n3+n5,3)*Q(n4+n6+n7,3)*Q(n2+n8,2)+6.*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6+n7,4)*Q(n2+n8,2)
      - 6.*Q(n3+n5,2)*Q(n1+n4+n6+n7,4)*Q(n2+n8,2)+6.*Q(n1,1)*Q(n5,1)*Q(n3+n4+n6+n7,4)*Q(n2+n8,2)
      - 6.*Q(n1+n5,2)*Q(n3+n4+n6+n7,4)*Q(n2+n8,2)-24.*Q(n5,1)*Q(n1+n3+n4+n6+n7,5)*Q(n2+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n2+n8,2)+2.*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n2+n8,2)
      + 2.*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6+n7,3)*Q(n2+n8,2)+2.*Q(n1,1)*Q(n3+n4,2)*Q(n5+n6+n7,3)*Q(n2+n8,2)
      - 4.*Q(n1+n3+n4,3)*Q(n5+n6+n7,3)*Q(n2+n8,2)+6.*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6+n7,4)*Q(n2+n8,2)
      - 6.*Q(n3+n4,2)*Q(n1+n5+n6+n7,4)*Q(n2+n8,2)+6.*Q(n1,1)*Q(n4,1)*Q(n3+n5+n6+n7,4)*Q(n2+n8,2)
      - 6.*Q(n1+n4,2)*Q(n3+n5+n6+n7,4)*Q(n2+n8,2)-24.*Q(n4,1)*Q(n1+n3+n5+n6+n7,5)*Q(n2+n8,2)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n4+n5+n6+n7,4)*Q(n2+n8,2)-6.*Q(n1+n3,2)*Q(n4+n5+n6+n7,4)*Q(n2+n8,2)
      - 24.*Q(n3,1)*Q(n1+n4+n5+n6+n7,5)*Q(n2+n8,2)-24.*Q(n1,1)*Q(n3+n4+n5+n6+n7,5)*Q(n2+n8,2)
      + 120.*Q(n1+n3+n4+n5+n6+n7,6)*Q(n2+n8,2)+2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n8,3)
      - 2.*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n8,3)-2.*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n8,3)
      - 2.*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n8,3)+4.*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n8,3)
      - 2.*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n2+n8,3)+2.*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n2+n8,3)
      - 2.*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n2+n8,3)+2.*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n2+n8,3)
      + 4.*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n1+n2+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n2+n8,3)
      + 2.*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n2+n8,3)+4.*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n1+n2+n8,3)
      + 4.*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n1+n2+n8,3)-12.*Q(n3+n4+n5+n6,4)*Q(n7,1)*Q(n1+n2+n8,3)
      - 2.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n2+n8,3)+2.*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n2+n8,3)
      + 2.*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n1+n2+n8,3)+2.*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n1+n2+n8,3)
      - 4.*Q(n4+n5+n6,3)*Q(n3+n7,2)*Q(n1+n2+n8,3)-2.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n2+n8,3)
      + 2.*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n2+n8,3)+2.*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n1+n2+n8,3)
      + 2.*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n1+n2+n8,3)-4.*Q(n3+n5+n6,3)*Q(n4+n7,2)*Q(n1+n2+n8,3)
      + 4.*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n1+n2+n8,3)-4.*Q(n5+n6,2)*Q(n3+n4+n7,3)*Q(n1+n2+n8,3)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n2+n8,3)+2.*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n2+n8,3)
      + 2.*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n1+n2+n8,3)+2.*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n1+n2+n8,3)
      - 4.*Q(n3+n4+n6,3)*Q(n5+n7,2)*Q(n1+n2+n8,3)+4.*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n1+n2+n8,3)
      - 4.*Q(n4+n6,2)*Q(n3+n5+n7,3)*Q(n1+n2+n8,3)+4.*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n1+n2+n8,3)
      - 4.*Q(n3+n6,2)*Q(n4+n5+n7,3)*Q(n1+n2+n8,3)-12.*Q(n6,1)*Q(n3+n4+n5+n7,4)*Q(n1+n2+n8,3)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n2+n8,3)+2.*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n2+n8,3)
      + 2.*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n1+n2+n8,3)+2.*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n1+n2+n8,3)
      - 4.*Q(n3+n4+n5,3)*Q(n6+n7,2)*Q(n1+n2+n8,3)+4.*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n1+n2+n8,3)
      - 4.*Q(n4+n5,2)*Q(n3+n6+n7,3)*Q(n1+n2+n8,3)+4.*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n1+n2+n8,3)
      - 4.*Q(n3+n5,2)*Q(n4+n6+n7,3)*Q(n1+n2+n8,3)-12.*Q(n5,1)*Q(n3+n4+n6+n7,4)*Q(n1+n2+n8,3)
      + 4.*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n1+n2+n8,3)-4.*Q(n3+n4,2)*Q(n5+n6+n7,3)*Q(n1+n2+n8,3)
      - 12.*Q(n4,1)*Q(n3+n5+n6+n7,4)*Q(n1+n2+n8,3)-12.*Q(n3,1)*Q(n4+n5+n6+n7,4)*Q(n1+n2+n8,3)
      + 48.*Q(n3+n4+n5+n6+n7,5)*Q(n1+n2+n8,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      - Q(n1+n2,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)+6.*Q(n1+n2+n4+n5,4)*Q(n6,1)*Q(n7,1)*Q(n3+n8,2)
      + Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n8,2)
      - Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n2,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n2+n4+n5,3)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n8,2)+Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n8,2)
      - Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n8,2)
      - Q(n1,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n8,2)+2.*Q(n1+n4+n5,3)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n8,2)
      - 2.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n3+n8,2)+2.*Q(n4+n5,2)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n3+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n1+n2,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n8,2)
      - Q(n2,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n1,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n1+n2+n5,3)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n8,2)-2.*Q(n2,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n2+n5,2)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n3+n8,2)-2.*Q(n1,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n1+n5,2)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n3+n8,2)+6.*Q(n5,1)*Q(n1+n2+n4+n6,4)*Q(n7,1)*Q(n3+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n1+n2,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n8,2)
      - Q(n2,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n8,2)-Q(n1,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n1+n2+n4,3)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n2+n4,2)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n3+n8,2)-2.*Q(n1,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n3+n8,2)
      + 2.*Q(n1+n4,2)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n3+n8,2)+6.*Q(n4,1)*Q(n1+n2+n5+n6,4)*Q(n7,1)*Q(n3+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n3+n8,2)+2.*Q(n1+n2,2)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n3+n8,2)
      + 6.*Q(n2,1)*Q(n1+n4+n5+n6,4)*Q(n7,1)*Q(n3+n8,2)+6.*Q(n1,1)*Q(n2+n4+n5+n6,4)*Q(n7,1)*Q(n3+n8,2)
      - 24.*Q(n1+n2+n4+n5+n6,5)*Q(n7,1)*Q(n3+n8,2)+Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n8,2)
      - Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n8,2)-Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n8,2)
      - Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n8,2)+2.*Q(n2+n4+n5,3)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n8,2)
      - Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n3+n8,2)+Q(n4+n5,2)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n3+n8,2)
      - Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n3+n8,2)+Q(n2+n5,2)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n3+n8,2)
      + 2.*Q(n5,1)*Q(n2+n4+n6,3)*Q(n1+n7,2)*Q(n3+n8,2)-Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n3+n8,2)
      + Q(n2+n4,2)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n3+n8,2)+2.*Q(n4,1)*Q(n2+n5+n6,3)*Q(n1+n7,2)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n4+n5+n6,3)*Q(n1+n7,2)*Q(n3+n8,2)-6.*Q(n2+n4+n5+n6,4)*Q(n1+n7,2)*Q(n3+n8,2)
      + Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n8,2)-Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n8,2)
      - Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n8,2)
      + 2.*Q(n1+n4+n5,3)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n8,2)-Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n3+n8,2)
      + Q(n4+n5,2)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n3+n8,2)
      + Q(n1+n5,2)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n3+n8,2)+2.*Q(n5,1)*Q(n1+n4+n6,3)*Q(n2+n7,2)*Q(n3+n8,2)
      - Q(n1,1)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n3+n8,2)+Q(n1+n4,2)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n3+n8,2)
      + 2.*Q(n4,1)*Q(n1+n5+n6,3)*Q(n2+n7,2)*Q(n3+n8,2)+2.*Q(n1,1)*Q(n4+n5+n6,3)*Q(n2+n7,2)*Q(n3+n8,2)
      - 6.*Q(n1+n4+n5+n6,4)*Q(n2+n7,2)*Q(n3+n8,2)-2.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n3+n8,2)
      + 2.*Q(n4+n5,2)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n3+n8,2)+2.*Q(n5,1)*Q(n4+n6,2)*Q(n1+n2+n7,3)*Q(n3+n8,2)
      + 2.*Q(n4,1)*Q(n5+n6,2)*Q(n1+n2+n7,3)*Q(n3+n8,2)-4.*Q(n4+n5+n6,3)*Q(n1+n2+n7,3)*Q(n3+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n8,2)-Q(n1+n2,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n8,2)
      - Q(n2,1)*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n8,2)
      + 2.*Q(n1+n2+n5,3)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n8,2)-Q(n2,1)*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n3+n8,2)
      + Q(n2+n5,2)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n3+n8,2)
      + Q(n1+n5,2)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n3+n8,2)+2.*Q(n5,1)*Q(n1+n2+n6,3)*Q(n4+n7,2)*Q(n3+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n3+n8,2)+Q(n1+n2,2)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n1+n5+n6,3)*Q(n4+n7,2)*Q(n3+n8,2)+2.*Q(n1,1)*Q(n2+n5+n6,3)*Q(n4+n7,2)*Q(n3+n8,2)
      - 6.*Q(n1+n2+n5+n6,4)*Q(n4+n7,2)*Q(n3+n8,2)-2.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n3+n8,2)
      + 2.*Q(n2+n5,2)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n3+n8,2)+2.*Q(n5,1)*Q(n2+n6,2)*Q(n1+n4+n7,3)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n5+n6,2)*Q(n1+n4+n7,3)*Q(n3+n8,2)-4.*Q(n2+n5+n6,3)*Q(n1+n4+n7,3)*Q(n3+n8,2)
      - 2.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n3+n8,2)+2.*Q(n1+n5,2)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n3+n8,2)
      + 2.*Q(n5,1)*Q(n1+n6,2)*Q(n2+n4+n7,3)*Q(n3+n8,2)+2.*Q(n1,1)*Q(n5+n6,2)*Q(n2+n4+n7,3)*Q(n3+n8,2)
      - 4.*Q(n1+n5+n6,3)*Q(n2+n4+n7,3)*Q(n3+n8,2)+6.*Q(n5,1)*Q(n6,1)*Q(n1+n2+n4+n7,4)*Q(n3+n8,2)
      - 6.*Q(n5+n6,2)*Q(n1+n2+n4+n7,4)*Q(n3+n8,2)+Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n8,2)
      - Q(n1+n2,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n8,2)-Q(n2,1)*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n8,2)
      - Q(n1,1)*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n8,2)+2.*Q(n1+n2+n4,3)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n8,2)
      - Q(n2,1)*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n3+n8,2)+Q(n2+n4,2)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n3+n8,2)
      - Q(n1,1)*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n3+n8,2)+Q(n1+n4,2)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n3+n8,2)
      + 2.*Q(n4,1)*Q(n1+n2+n6,3)*Q(n5+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n2,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n3+n8,2)
      + Q(n1+n2,2)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n3+n8,2)+2.*Q(n2,1)*Q(n1+n4+n6,3)*Q(n5+n7,2)*Q(n3+n8,2)
      + 2.*Q(n1,1)*Q(n2+n4+n6,3)*Q(n5+n7,2)*Q(n3+n8,2)-6.*Q(n1+n2+n4+n6,4)*Q(n5+n7,2)*Q(n3+n8,2)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n3+n8,2)+2.*Q(n2+n4,2)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n3+n8,2)
      + 2.*Q(n4,1)*Q(n2+n6,2)*Q(n1+n5+n7,3)*Q(n3+n8,2)+2.*Q(n2,1)*Q(n4+n6,2)*Q(n1+n5+n7,3)*Q(n3+n8,2)
      - 4.*Q(n2+n4+n6,3)*Q(n1+n5+n7,3)*Q(n3+n8,2)-2.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n3+n8,2)
      + 2.*Q(n1+n4,2)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n3+n8,2)+2.*Q(n4,1)*Q(n1+n6,2)*Q(n2+n5+n7,3)*Q(n3+n8,2)
      + 2.*Q(n1,1)*Q(n4+n6,2)*Q(n2+n5+n7,3)*Q(n3+n8,2)-4.*Q(n1+n4+n6,3)*Q(n2+n5+n7,3)*Q(n3+n8,2)
      + 6.*Q(n4,1)*Q(n6,1)*Q(n1+n2+n5+n7,4)*Q(n3+n8,2)-6.*Q(n4+n6,2)*Q(n1+n2+n5+n7,4)*Q(n3+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n3+n8,2)+2.*Q(n1+n2,2)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n1+n6,2)*Q(n4+n5+n7,3)*Q(n3+n8,2)+2.*Q(n1,1)*Q(n2+n6,2)*Q(n4+n5+n7,3)*Q(n3+n8,2)
      - 4.*Q(n1+n2+n6,3)*Q(n4+n5+n7,3)*Q(n3+n8,2)+6.*Q(n2,1)*Q(n6,1)*Q(n1+n4+n5+n7,4)*Q(n3+n8,2)
      - 6.*Q(n2+n6,2)*Q(n1+n4+n5+n7,4)*Q(n3+n8,2)+6.*Q(n1,1)*Q(n6,1)*Q(n2+n4+n5+n7,4)*Q(n3+n8,2)
      - 6.*Q(n1+n6,2)*Q(n2+n4+n5+n7,4)*Q(n3+n8,2)-24.*Q(n6,1)*Q(n1+n2+n4+n5+n7,5)*Q(n3+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n8,2)-Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n8,2)
      - Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n8,2)
      + 2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n8,2)-Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n3+n8,2)
      + Q(n2+n4,2)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n3+n8,2)-Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n3+n8,2)
      + Q(n1+n4,2)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n3+n8,2)+2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6+n7,2)*Q(n3+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n3+n8,2)+Q(n1+n2,2)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n6+n7,2)*Q(n3+n8,2)+2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n6+n7,2)*Q(n3+n8,2)
      - 6.*Q(n1+n2+n4+n5,4)*Q(n6+n7,2)*Q(n3+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n3+n8,2)+2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n6+n7,3)*Q(n3+n8,2)-4.*Q(n2+n4+n5,3)*Q(n1+n6+n7,3)*Q(n3+n8,2)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n3+n8,2)+2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6+n7,3)*Q(n3+n8,2)+2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n6+n7,3)*Q(n3+n8,2)
      - 4.*Q(n1+n4+n5,3)*Q(n2+n6+n7,3)*Q(n3+n8,2)+6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6+n7,4)*Q(n3+n8,2)
      - 6.*Q(n4+n5,2)*Q(n1+n2+n6+n7,4)*Q(n3+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n1+n2,2)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n3+n8,2)+2.*Q(n2,1)*Q(n1+n5,2)*Q(n4+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n1,1)*Q(n2+n5,2)*Q(n4+n6+n7,3)*Q(n3+n8,2)-4.*Q(n1+n2+n5,3)*Q(n4+n6+n7,3)*Q(n3+n8,2)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n1+n4+n6+n7,4)*Q(n3+n8,2)-6.*Q(n2+n5,2)*Q(n1+n4+n6+n7,4)*Q(n3+n8,2)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n2+n4+n6+n7,4)*Q(n3+n8,2)-6.*Q(n1+n5,2)*Q(n2+n4+n6+n7,4)*Q(n3+n8,2)
      - 24.*Q(n5,1)*Q(n1+n2+n4+n6+n7,5)*Q(n3+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n1+n2,2)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n3+n8,2)+2.*Q(n2,1)*Q(n1+n4,2)*Q(n5+n6+n7,3)*Q(n3+n8,2)
      + 2.*Q(n1,1)*Q(n2+n4,2)*Q(n5+n6+n7,3)*Q(n3+n8,2)-4.*Q(n1+n2+n4,3)*Q(n5+n6+n7,3)*Q(n3+n8,2)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n1+n5+n6+n7,4)*Q(n3+n8,2)-6.*Q(n2+n4,2)*Q(n1+n5+n6+n7,4)*Q(n3+n8,2)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n2+n5+n6+n7,4)*Q(n3+n8,2)-6.*Q(n1+n4,2)*Q(n2+n5+n6+n7,4)*Q(n3+n8,2)
      - 24.*Q(n4,1)*Q(n1+n2+n5+n6+n7,5)*Q(n3+n8,2)+6.*Q(n1,1)*Q(n2,1)*Q(n4+n5+n6+n7,4)*Q(n3+n8,2)
      - 6.*Q(n1+n2,2)*Q(n4+n5+n6+n7,4)*Q(n3+n8,2)-24.*Q(n2,1)*Q(n1+n4+n5+n6+n7,5)*Q(n3+n8,2)
      - 24.*Q(n1,1)*Q(n2+n4+n5+n6+n7,5)*Q(n3+n8,2)+120.*Q(n1+n2+n4+n5+n6+n7,6)*Q(n3+n8,2)
      + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n8,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n8,3)
      - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n8,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n8,3)
      + 4.*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n8,3)-2.*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n3+n8,3)
      + 2.*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n3+n8,3)-2.*Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n3+n8,3)
      + 2.*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n3+n8,3)+4.*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n1+n3+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n3+n8,3)+2.*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n3+n8,3)
      + 4.*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n1+n3+n8,3)+4.*Q(n2,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n1+n3+n8,3)
      - 12.*Q(n2+n4+n5+n6,4)*Q(n7,1)*Q(n1+n3+n8,3)-2.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n3+n8,3)+2.*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n1+n3+n8,3)-4.*Q(n4+n5+n6,3)*Q(n2+n7,2)*Q(n1+n3+n8,3)
      - 2.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n3+n8,3)+2.*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n1+n3+n8,3)+2.*Q(n2,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n1+n3+n8,3)
      - 4.*Q(n2+n5+n6,3)*Q(n4+n7,2)*Q(n1+n3+n8,3)+4.*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n1+n3+n8,3)
      - 4.*Q(n5+n6,2)*Q(n2+n4+n7,3)*Q(n1+n3+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n3+n8,3)+2.*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n2,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n1+n3+n8,3)-4.*Q(n2+n4+n6,3)*Q(n5+n7,2)*Q(n1+n3+n8,3)
      + 4.*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n1+n3+n8,3)-4.*Q(n4+n6,2)*Q(n2+n5+n7,3)*Q(n1+n3+n8,3)
      + 4.*Q(n2,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n1+n3+n8,3)-4.*Q(n2+n6,2)*Q(n4+n5+n7,3)*Q(n1+n3+n8,3)
      - 12.*Q(n6,1)*Q(n2+n4+n5+n7,4)*Q(n1+n3+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n3+n8,3)+2.*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n1+n3+n8,3)
      + 2.*Q(n2,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n1+n3+n8,3)-4.*Q(n2+n4+n5,3)*Q(n6+n7,2)*Q(n1+n3+n8,3)
      + 4.*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n1+n3+n8,3)-4.*Q(n4+n5,2)*Q(n2+n6+n7,3)*Q(n1+n3+n8,3)
      + 4.*Q(n2,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n1+n3+n8,3)-4.*Q(n2+n5,2)*Q(n4+n6+n7,3)*Q(n1+n3+n8,3)
      - 12.*Q(n5,1)*Q(n2+n4+n6+n7,4)*Q(n1+n3+n8,3)+4.*Q(n2,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n1+n3+n8,3)
      - 4.*Q(n2+n4,2)*Q(n5+n6+n7,3)*Q(n1+n3+n8,3)-12.*Q(n4,1)*Q(n2+n5+n6+n7,4)*Q(n1+n3+n8,3)
      - 12.*Q(n2,1)*Q(n4+n5+n6+n7,4)*Q(n1+n3+n8,3)+48.*Q(n2+n4+n5+n6+n7,5)*Q(n1+n3+n8,3)
      + 2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n8,3)-2.*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n8,3)
      - 2.*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n8,3)-2.*Q(n1,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n8,3)
      + 4.*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n8,3)-2.*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n3+n8,3)
      + 2.*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n3+n8,3)-2.*Q(n1,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n3+n8,3)
      + 2.*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n3+n8,3)+4.*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n2+n3+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n3+n8,3)+2.*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n3+n8,3)
      + 4.*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n2+n3+n8,3)+4.*Q(n1,1)*Q(n4+n5+n6,3)*Q(n7,1)*Q(n2+n3+n8,3)
      - 12.*Q(n1+n4+n5+n6,4)*Q(n7,1)*Q(n2+n3+n8,3)-2.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n3+n8,3)+2.*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n2+n3+n8,3)-4.*Q(n4+n5+n6,3)*Q(n1+n7,2)*Q(n2+n3+n8,3)
      - 2.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n3+n8,3)+2.*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n2+n3+n8,3)+2.*Q(n1,1)*Q(n5+n6,2)*Q(n4+n7,2)*Q(n2+n3+n8,3)
      - 4.*Q(n1+n5+n6,3)*Q(n4+n7,2)*Q(n2+n3+n8,3)+4.*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n2+n3+n8,3)
      - 4.*Q(n5+n6,2)*Q(n1+n4+n7,3)*Q(n2+n3+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n3+n8,3)+2.*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n1,1)*Q(n4+n6,2)*Q(n5+n7,2)*Q(n2+n3+n8,3)-4.*Q(n1+n4+n6,3)*Q(n5+n7,2)*Q(n2+n3+n8,3)
      + 4.*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n2+n3+n8,3)-4.*Q(n4+n6,2)*Q(n1+n5+n7,3)*Q(n2+n3+n8,3)
      + 4.*Q(n1,1)*Q(n6,1)*Q(n4+n5+n7,3)*Q(n2+n3+n8,3)-4.*Q(n1+n6,2)*Q(n4+n5+n7,3)*Q(n2+n3+n8,3)
      - 12.*Q(n6,1)*Q(n1+n4+n5+n7,4)*Q(n2+n3+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n3+n8,3)+2.*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n2+n3+n8,3)
      + 2.*Q(n1,1)*Q(n4+n5,2)*Q(n6+n7,2)*Q(n2+n3+n8,3)-4.*Q(n1+n4+n5,3)*Q(n6+n7,2)*Q(n2+n3+n8,3)
      + 4.*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n2+n3+n8,3)-4.*Q(n4+n5,2)*Q(n1+n6+n7,3)*Q(n2+n3+n8,3)
      + 4.*Q(n1,1)*Q(n5,1)*Q(n4+n6+n7,3)*Q(n2+n3+n8,3)-4.*Q(n1+n5,2)*Q(n4+n6+n7,3)*Q(n2+n3+n8,3)
      - 12.*Q(n5,1)*Q(n1+n4+n6+n7,4)*Q(n2+n3+n8,3)+4.*Q(n1,1)*Q(n4,1)*Q(n5+n6+n7,3)*Q(n2+n3+n8,3)
      - 4.*Q(n1+n4,2)*Q(n5+n6+n7,3)*Q(n2+n3+n8,3)-12.*Q(n4,1)*Q(n1+n5+n6+n7,4)*Q(n2+n3+n8,3)
      - 12.*Q(n1,1)*Q(n4+n5+n6+n7,4)*Q(n2+n3+n8,3)+48.*Q(n1+n4+n5+n6+n7,5)*Q(n2+n3+n8,3)
      - 6.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n3+n8,4)+6.*Q(n4+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n3+n8,4)
      + 6.*Q(n5,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n2+n3+n8,4)+6.*Q(n4,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n2+n3+n8,4)
      - 12.*Q(n4+n5+n6,3)*Q(n7,1)*Q(n1+n2+n3+n8,4)+6.*Q(n5,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n2+n3+n8,4)
      - 6.*Q(n5+n6,2)*Q(n4+n7,2)*Q(n1+n2+n3+n8,4)+6.*Q(n4,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n2+n3+n8,4)
      - 6.*Q(n4+n6,2)*Q(n5+n7,2)*Q(n1+n2+n3+n8,4)-12.*Q(n6,1)*Q(n4+n5+n7,3)*Q(n1+n2+n3+n8,4)
      + 6.*Q(n4,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n2+n3+n8,4)-6.*Q(n4+n5,2)*Q(n6+n7,2)*Q(n1+n2+n3+n8,4)
      - 12.*Q(n5,1)*Q(n4+n6+n7,3)*Q(n1+n2+n3+n8,4)-12.*Q(n4,1)*Q(n5+n6+n7,3)*Q(n1+n2+n3+n8,4)
      + 36.*Q(n4+n5+n6+n7,4)*Q(n1+n2+n3+n8,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      + Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)+Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      + Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)-2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)-Q(n2+n3,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)-Q(n1+n3,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      - 2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      - Q(n1+n2,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)-2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      - 2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)+6.*Q(n1+n2+n3+n5,4)*Q(n6,1)*Q(n7,1)*Q(n4+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n2+n3,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n8,2)
      - Q(n3,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n2,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n2+n3+n5,3)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n8,2)+Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n8,2)
      - Q(n1+n3,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n3,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n8,2)
      - Q(n1,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n8,2)+2.*Q(n1+n3+n5,3)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n8,2)
      - 2.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n4+n8,2)+2.*Q(n3+n5,2)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n4+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n1+n2,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n8,2)
      - Q(n2,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n1,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n1+n2+n5,3)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n8,2)-2.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n2+n5,2)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n4+n8,2)-2.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n1+n5,2)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n4+n8,2)+6.*Q(n5,1)*Q(n1+n2+n3+n6,4)*Q(n7,1)*Q(n4+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n5+n6,2)*Q(n7,1)*Q(n4+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n5+n6,2)*Q(n7,1)*Q(n4+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n5+n6,2)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n5+n6,2)*Q(n7,1)*Q(n4+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n2+n3,2)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n4+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n4+n8,2)
      + 2.*Q(n1+n3,2)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n4+n8,2)+6.*Q(n3,1)*Q(n1+n2+n5+n6,4)*Q(n7,1)*Q(n4+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n4+n8,2)+2.*Q(n1+n2,2)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n4+n8,2)
      + 6.*Q(n2,1)*Q(n1+n3+n5+n6,4)*Q(n7,1)*Q(n4+n8,2)+6.*Q(n1,1)*Q(n2+n3+n5+n6,4)*Q(n7,1)*Q(n4+n8,2)
      - 24.*Q(n1+n2+n3+n5+n6,5)*Q(n7,1)*Q(n4+n8,2)+Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n8,2)
      - Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n8,2)-Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n8,2)
      - Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n8,2)+2.*Q(n2+n3+n5,3)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n8,2)
      - Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n4+n8,2)+Q(n3+n5,2)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n4+n8,2)
      - Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n4+n8,2)+Q(n2+n5,2)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n4+n8,2)
      + 2.*Q(n5,1)*Q(n2+n3+n6,3)*Q(n1+n7,2)*Q(n4+n8,2)-Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n4+n8,2)
      + Q(n2+n3,2)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n4+n8,2)+2.*Q(n3,1)*Q(n2+n5+n6,3)*Q(n1+n7,2)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n3+n5+n6,3)*Q(n1+n7,2)*Q(n4+n8,2)-6.*Q(n2+n3+n5+n6,4)*Q(n1+n7,2)*Q(n4+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n8,2)-Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n8,2)
      - Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n8,2)
      + 2.*Q(n1+n3+n5,3)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n8,2)-Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n4+n8,2)
      + Q(n3+n5,2)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n4+n8,2)
      + Q(n1+n5,2)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n4+n8,2)+2.*Q(n5,1)*Q(n1+n3+n6,3)*Q(n2+n7,2)*Q(n4+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n4+n8,2)+Q(n1+n3,2)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n4+n8,2)
      + 2.*Q(n3,1)*Q(n1+n5+n6,3)*Q(n2+n7,2)*Q(n4+n8,2)+2.*Q(n1,1)*Q(n3+n5+n6,3)*Q(n2+n7,2)*Q(n4+n8,2)
      - 6.*Q(n1+n3+n5+n6,4)*Q(n2+n7,2)*Q(n4+n8,2)-2.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n4+n8,2)
      + 2.*Q(n3+n5,2)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n4+n8,2)+2.*Q(n5,1)*Q(n3+n6,2)*Q(n1+n2+n7,3)*Q(n4+n8,2)
      + 2.*Q(n3,1)*Q(n5+n6,2)*Q(n1+n2+n7,3)*Q(n4+n8,2)-4.*Q(n3+n5+n6,3)*Q(n1+n2+n7,3)*Q(n4+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n8,2)-Q(n1+n2,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n8,2)
      - Q(n2,1)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n8,2)
      + 2.*Q(n1+n2+n5,3)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n8,2)-Q(n2,1)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n4+n8,2)
      + Q(n2+n5,2)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n4+n8,2)
      + Q(n1+n5,2)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n4+n8,2)+2.*Q(n5,1)*Q(n1+n2+n6,3)*Q(n3+n7,2)*Q(n4+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n4+n8,2)+Q(n1+n2,2)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n1+n5+n6,3)*Q(n3+n7,2)*Q(n4+n8,2)+2.*Q(n1,1)*Q(n2+n5+n6,3)*Q(n3+n7,2)*Q(n4+n8,2)
      - 6.*Q(n1+n2+n5+n6,4)*Q(n3+n7,2)*Q(n4+n8,2)-2.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n4+n8,2)
      + 2.*Q(n2+n5,2)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n4+n8,2)+2.*Q(n5,1)*Q(n2+n6,2)*Q(n1+n3+n7,3)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n5+n6,2)*Q(n1+n3+n7,3)*Q(n4+n8,2)-4.*Q(n2+n5+n6,3)*Q(n1+n3+n7,3)*Q(n4+n8,2)
      - 2.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n4+n8,2)+2.*Q(n1+n5,2)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n4+n8,2)
      + 2.*Q(n5,1)*Q(n1+n6,2)*Q(n2+n3+n7,3)*Q(n4+n8,2)+2.*Q(n1,1)*Q(n5+n6,2)*Q(n2+n3+n7,3)*Q(n4+n8,2)
      - 4.*Q(n1+n5+n6,3)*Q(n2+n3+n7,3)*Q(n4+n8,2)+6.*Q(n5,1)*Q(n6,1)*Q(n1+n2+n3+n7,4)*Q(n4+n8,2)
      - 6.*Q(n5+n6,2)*Q(n1+n2+n3+n7,4)*Q(n4+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n5+n7,2)*Q(n4+n8,2)
      - Q(n1+n2,2)*Q(n3,1)*Q(n6,1)*Q(n5+n7,2)*Q(n4+n8,2)-Q(n2,1)*Q(n1+n3,2)*Q(n6,1)*Q(n5+n7,2)*Q(n4+n8,2)
      - Q(n1,1)*Q(n2+n3,2)*Q(n6,1)*Q(n5+n7,2)*Q(n4+n8,2)+2.*Q(n1+n2+n3,3)*Q(n6,1)*Q(n5+n7,2)*Q(n4+n8,2)
      - Q(n2,1)*Q(n3,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n4+n8,2)+Q(n2+n3,2)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n4+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n4+n8,2)+Q(n1+n3,2)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n4+n8,2)
      + 2.*Q(n3,1)*Q(n1+n2+n6,3)*Q(n5+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n2,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n4+n8,2)
      + Q(n1+n2,2)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n4+n8,2)+2.*Q(n2,1)*Q(n1+n3+n6,3)*Q(n5+n7,2)*Q(n4+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3+n6,3)*Q(n5+n7,2)*Q(n4+n8,2)-6.*Q(n1+n2+n3+n6,4)*Q(n5+n7,2)*Q(n4+n8,2)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n4+n8,2)+2.*Q(n2+n3,2)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n4+n8,2)
      + 2.*Q(n3,1)*Q(n2+n6,2)*Q(n1+n5+n7,3)*Q(n4+n8,2)+2.*Q(n2,1)*Q(n3+n6,2)*Q(n1+n5+n7,3)*Q(n4+n8,2)
      - 4.*Q(n2+n3+n6,3)*Q(n1+n5+n7,3)*Q(n4+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n4+n8,2)
      + 2.*Q(n1+n3,2)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n4+n8,2)+2.*Q(n3,1)*Q(n1+n6,2)*Q(n2+n5+n7,3)*Q(n4+n8,2)
      + 2.*Q(n1,1)*Q(n3+n6,2)*Q(n2+n5+n7,3)*Q(n4+n8,2)-4.*Q(n1+n3+n6,3)*Q(n2+n5+n7,3)*Q(n4+n8,2)
      + 6.*Q(n3,1)*Q(n6,1)*Q(n1+n2+n5+n7,4)*Q(n4+n8,2)-6.*Q(n3+n6,2)*Q(n1+n2+n5+n7,4)*Q(n4+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n4+n8,2)+2.*Q(n1+n2,2)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n1+n6,2)*Q(n3+n5+n7,3)*Q(n4+n8,2)+2.*Q(n1,1)*Q(n2+n6,2)*Q(n3+n5+n7,3)*Q(n4+n8,2)
      - 4.*Q(n1+n2+n6,3)*Q(n3+n5+n7,3)*Q(n4+n8,2)+6.*Q(n2,1)*Q(n6,1)*Q(n1+n3+n5+n7,4)*Q(n4+n8,2)
      - 6.*Q(n2+n6,2)*Q(n1+n3+n5+n7,4)*Q(n4+n8,2)+6.*Q(n1,1)*Q(n6,1)*Q(n2+n3+n5+n7,4)*Q(n4+n8,2)
      - 6.*Q(n1+n6,2)*Q(n2+n3+n5+n7,4)*Q(n4+n8,2)-24.*Q(n6,1)*Q(n1+n2+n3+n5+n7,5)*Q(n4+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6+n7,2)*Q(n4+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n6+n7,2)*Q(n4+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n6+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n6+n7,2)*Q(n4+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n6+n7,2)*Q(n4+n8,2)-Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n4+n8,2)
      + Q(n2+n3,2)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n4+n8,2)-Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n4+n8,2)
      + Q(n1+n3,2)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n4+n8,2)+2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n6+n7,2)*Q(n4+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n4+n8,2)+Q(n1+n2,2)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n6+n7,2)*Q(n4+n8,2)+2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n6+n7,2)*Q(n4+n8,2)
      - 6.*Q(n1+n2+n3+n5,4)*Q(n6+n7,2)*Q(n4+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n4+n8,2)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n6+n7,3)*Q(n4+n8,2)-4.*Q(n2+n3+n5,3)*Q(n1+n6+n7,3)*Q(n4+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n4+n8,2)+2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n6+n7,3)*Q(n4+n8,2)+2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n6+n7,3)*Q(n4+n8,2)
      - 4.*Q(n1+n3+n5,3)*Q(n2+n6+n7,3)*Q(n4+n8,2)+6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n6+n7,4)*Q(n4+n8,2)
      - 6.*Q(n3+n5,2)*Q(n1+n2+n6+n7,4)*Q(n4+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n4+n8,2)+2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n6+n7,3)*Q(n4+n8,2)-4.*Q(n1+n2+n5,3)*Q(n3+n6+n7,3)*Q(n4+n8,2)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n6+n7,4)*Q(n4+n8,2)-6.*Q(n2+n5,2)*Q(n1+n3+n6+n7,4)*Q(n4+n8,2)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n6+n7,4)*Q(n4+n8,2)-6.*Q(n1+n5,2)*Q(n2+n3+n6+n7,4)*Q(n4+n8,2)
      - 24.*Q(n5,1)*Q(n1+n2+n3+n6+n7,5)*Q(n4+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n5+n6+n7,3)*Q(n4+n8,2)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n5+n6+n7,3)*Q(n4+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n5+n6+n7,3)*Q(n4+n8,2)-4.*Q(n1+n2+n3,3)*Q(n5+n6+n7,3)*Q(n4+n8,2)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n1+n5+n6+n7,4)*Q(n4+n8,2)-6.*Q(n2+n3,2)*Q(n1+n5+n6+n7,4)*Q(n4+n8,2)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n2+n5+n6+n7,4)*Q(n4+n8,2)-6.*Q(n1+n3,2)*Q(n2+n5+n6+n7,4)*Q(n4+n8,2)
      - 24.*Q(n3,1)*Q(n1+n2+n5+n6+n7,5)*Q(n4+n8,2)+6.*Q(n1,1)*Q(n2,1)*Q(n3+n5+n6+n7,4)*Q(n4+n8,2)
      - 6.*Q(n1+n2,2)*Q(n3+n5+n6+n7,4)*Q(n4+n8,2)-24.*Q(n2,1)*Q(n1+n3+n5+n6+n7,5)*Q(n4+n8,2)
      - 24.*Q(n1,1)*Q(n2+n3+n5+n6+n7,5)*Q(n4+n8,2)+120.*Q(n1+n2+n3+n5+n6+n7,6)*Q(n4+n8,2)
      + 2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n8,3)-2.*Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n8,3)
      - 2.*Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n8,3)-2.*Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n8,3)
      + 4.*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n8,3)-2.*Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n4+n8,3)
      + 2.*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n4+n8,3)-2.*Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n4+n8,3)
      + 2.*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n4+n8,3)+4.*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n1+n4+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n4+n8,3)+2.*Q(n2+n3,2)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n4+n8,3)
      + 4.*Q(n3,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n1+n4+n8,3)+4.*Q(n2,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n1+n4+n8,3)
      - 12.*Q(n2+n3+n5+n6,4)*Q(n7,1)*Q(n1+n4+n8,3)-2.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n4+n8,3)+2.*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n3,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n1+n4+n8,3)-4.*Q(n3+n5+n6,3)*Q(n2+n7,2)*Q(n1+n4+n8,3)
      - 2.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n4+n8,3)+2.*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n1+n4+n8,3)+2.*Q(n2,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n1+n4+n8,3)
      - 4.*Q(n2+n5+n6,3)*Q(n3+n7,2)*Q(n1+n4+n8,3)+4.*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n1+n4+n8,3)
      - 4.*Q(n5+n6,2)*Q(n2+n3+n7,3)*Q(n1+n4+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n2+n3,2)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n4+n8,3)+2.*Q(n3,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n2,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n1+n4+n8,3)-4.*Q(n2+n3+n6,3)*Q(n5+n7,2)*Q(n1+n4+n8,3)
      + 4.*Q(n3,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n1+n4+n8,3)-4.*Q(n3+n6,2)*Q(n2+n5+n7,3)*Q(n1+n4+n8,3)
      + 4.*Q(n2,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n1+n4+n8,3)-4.*Q(n2+n6,2)*Q(n3+n5+n7,3)*Q(n1+n4+n8,3)
      - 12.*Q(n6,1)*Q(n2+n3+n5+n7,4)*Q(n1+n4+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n4+n8,3)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n1+n4+n8,3)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n1+n4+n8,3)-4.*Q(n2+n3+n5,3)*Q(n6+n7,2)*Q(n1+n4+n8,3)
      + 4.*Q(n3,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n1+n4+n8,3)-4.*Q(n3+n5,2)*Q(n2+n6+n7,3)*Q(n1+n4+n8,3)
      + 4.*Q(n2,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n1+n4+n8,3)-4.*Q(n2+n5,2)*Q(n3+n6+n7,3)*Q(n1+n4+n8,3)
      - 12.*Q(n5,1)*Q(n2+n3+n6+n7,4)*Q(n1+n4+n8,3)+4.*Q(n2,1)*Q(n3,1)*Q(n5+n6+n7,3)*Q(n1+n4+n8,3)
      - 4.*Q(n2+n3,2)*Q(n5+n6+n7,3)*Q(n1+n4+n8,3)-12.*Q(n3,1)*Q(n2+n5+n6+n7,4)*Q(n1+n4+n8,3)
      - 12.*Q(n2,1)*Q(n3+n5+n6+n7,4)*Q(n1+n4+n8,3)+48.*Q(n2+n3+n5+n6+n7,5)*Q(n1+n4+n8,3)
      + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n8,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n8,3)
      - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n8,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n8,3)
      + 4.*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n8,3)-2.*Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n4+n8,3)
      + 2.*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n4+n8,3)-2.*Q(n1,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n4+n8,3)
      + 2.*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n4+n8,3)+4.*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n2+n4+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n4+n8,3)+2.*Q(n1+n3,2)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n4+n8,3)
      + 4.*Q(n3,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n2+n4+n8,3)+4.*Q(n1,1)*Q(n3+n5+n6,3)*Q(n7,1)*Q(n2+n4+n8,3)
      - 12.*Q(n1+n3+n5+n6,4)*Q(n7,1)*Q(n2+n4+n8,3)-2.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n4+n8,3)+2.*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n3,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n2+n4+n8,3)-4.*Q(n3+n5+n6,3)*Q(n1+n7,2)*Q(n2+n4+n8,3)
      - 2.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n4+n8,3)+2.*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n2+n4+n8,3)+2.*Q(n1,1)*Q(n5+n6,2)*Q(n3+n7,2)*Q(n2+n4+n8,3)
      - 4.*Q(n1+n5+n6,3)*Q(n3+n7,2)*Q(n2+n4+n8,3)+4.*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n2+n4+n8,3)
      - 4.*Q(n5+n6,2)*Q(n1+n3+n7,3)*Q(n2+n4+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n1+n3,2)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n4+n8,3)+2.*Q(n3,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n1,1)*Q(n3+n6,2)*Q(n5+n7,2)*Q(n2+n4+n8,3)-4.*Q(n1+n3+n6,3)*Q(n5+n7,2)*Q(n2+n4+n8,3)
      + 4.*Q(n3,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n2+n4+n8,3)-4.*Q(n3+n6,2)*Q(n1+n5+n7,3)*Q(n2+n4+n8,3)
      + 4.*Q(n1,1)*Q(n6,1)*Q(n3+n5+n7,3)*Q(n2+n4+n8,3)-4.*Q(n1+n6,2)*Q(n3+n5+n7,3)*Q(n2+n4+n8,3)
      - 12.*Q(n6,1)*Q(n1+n3+n5+n7,4)*Q(n2+n4+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n1+n3,2)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n4+n8,3)+2.*Q(n3,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n2+n4+n8,3)
      + 2.*Q(n1,1)*Q(n3+n5,2)*Q(n6+n7,2)*Q(n2+n4+n8,3)-4.*Q(n1+n3+n5,3)*Q(n6+n7,2)*Q(n2+n4+n8,3)
      + 4.*Q(n3,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n2+n4+n8,3)-4.*Q(n3+n5,2)*Q(n1+n6+n7,3)*Q(n2+n4+n8,3)
      + 4.*Q(n1,1)*Q(n5,1)*Q(n3+n6+n7,3)*Q(n2+n4+n8,3)-4.*Q(n1+n5,2)*Q(n3+n6+n7,3)*Q(n2+n4+n8,3)
      - 12.*Q(n5,1)*Q(n1+n3+n6+n7,4)*Q(n2+n4+n8,3)+4.*Q(n1,1)*Q(n3,1)*Q(n5+n6+n7,3)*Q(n2+n4+n8,3)
      - 4.*Q(n1+n3,2)*Q(n5+n6+n7,3)*Q(n2+n4+n8,3)-12.*Q(n3,1)*Q(n1+n5+n6+n7,4)*Q(n2+n4+n8,3)
      - 12.*Q(n1,1)*Q(n3+n5+n6+n7,4)*Q(n2+n4+n8,3)+48.*Q(n1+n3+n5+n6+n7,5)*Q(n2+n4+n8,3)
      - 6.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n4+n8,4)+6.*Q(n3+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n4+n8,4)
      + 6.*Q(n5,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n2+n4+n8,4)+6.*Q(n3,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n2+n4+n8,4)
      - 12.*Q(n3+n5+n6,3)*Q(n7,1)*Q(n1+n2+n4+n8,4)+6.*Q(n5,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n2+n4+n8,4)
      - 6.*Q(n5+n6,2)*Q(n3+n7,2)*Q(n1+n2+n4+n8,4)+6.*Q(n3,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n2+n4+n8,4)
      - 6.*Q(n3+n6,2)*Q(n5+n7,2)*Q(n1+n2+n4+n8,4)-12.*Q(n6,1)*Q(n3+n5+n7,3)*Q(n1+n2+n4+n8,4)
      + 6.*Q(n3,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n2+n4+n8,4)-6.*Q(n3+n5,2)*Q(n6+n7,2)*Q(n1+n2+n4+n8,4)
      - 12.*Q(n5,1)*Q(n3+n6+n7,3)*Q(n1+n2+n4+n8,4)-12.*Q(n3,1)*Q(n5+n6+n7,3)*Q(n1+n2+n4+n8,4)
      + 36.*Q(n3+n5+n6+n7,4)*Q(n1+n2+n4+n8,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n8,3)
      - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n8,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n8,3)
      - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n8,3)+4.*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n8,3)
      - 2.*Q(n2,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n4+n8,3)+2.*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n4+n8,3)
      - 2.*Q(n1,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n4+n8,3)+2.*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n4+n8,3)
      + 4.*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n3+n4+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n4+n8,3)
      + 2.*Q(n1+n2,2)*Q(n5+n6,2)*Q(n7,1)*Q(n3+n4+n8,3)+4.*Q(n2,1)*Q(n1+n5+n6,3)*Q(n7,1)*Q(n3+n4+n8,3)
      + 4.*Q(n1,1)*Q(n2+n5+n6,3)*Q(n7,1)*Q(n3+n4+n8,3)-12.*Q(n1+n2+n5+n6,4)*Q(n7,1)*Q(n3+n4+n8,3)
      - 2.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n4+n8,3)+2.*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n4+n8,3)
      + 2.*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n3+n4+n8,3)+2.*Q(n2,1)*Q(n5+n6,2)*Q(n1+n7,2)*Q(n3+n4+n8,3)
      - 4.*Q(n2+n5+n6,3)*Q(n1+n7,2)*Q(n3+n4+n8,3)-2.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n4+n8,3)
    );
} 
//=======================================================
TComplex AliAnalysisTaskChargedFlow::Eight_3(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	return	
		(
			2.*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n4+n8,3)+2.*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n3+n4+n8,3)
      + 2.*Q(n1,1)*Q(n5+n6,2)*Q(n2+n7,2)*Q(n3+n4+n8,3)-4.*Q(n1+n5+n6,3)*Q(n2+n7,2)*Q(n3+n4+n8,3)
      + 4.*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n3+n4+n8,3)-4.*Q(n5+n6,2)*Q(n1+n2+n7,3)*Q(n3+n4+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n4+n8,3)+2.*Q(n1+n2,2)*Q(n6,1)*Q(n5+n7,2)*Q(n3+n4+n8,3)
      + 2.*Q(n2,1)*Q(n1+n6,2)*Q(n5+n7,2)*Q(n3+n4+n8,3)+2.*Q(n1,1)*Q(n2+n6,2)*Q(n5+n7,2)*Q(n3+n4+n8,3)
      - 4.*Q(n1+n2+n6,3)*Q(n5+n7,2)*Q(n3+n4+n8,3)+4.*Q(n2,1)*Q(n6,1)*Q(n1+n5+n7,3)*Q(n3+n4+n8,3)
      - 4.*Q(n2+n6,2)*Q(n1+n5+n7,3)*Q(n3+n4+n8,3)+4.*Q(n1,1)*Q(n6,1)*Q(n2+n5+n7,3)*Q(n3+n4+n8,3)
      - 4.*Q(n1+n6,2)*Q(n2+n5+n7,3)*Q(n3+n4+n8,3)-12.*Q(n6,1)*Q(n1+n2+n5+n7,4)*Q(n3+n4+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n4+n8,3)+2.*Q(n1+n2,2)*Q(n5,1)*Q(n6+n7,2)*Q(n3+n4+n8,3)
      + 2.*Q(n2,1)*Q(n1+n5,2)*Q(n6+n7,2)*Q(n3+n4+n8,3)+2.*Q(n1,1)*Q(n2+n5,2)*Q(n6+n7,2)*Q(n3+n4+n8,3)
      - 4.*Q(n1+n2+n5,3)*Q(n6+n7,2)*Q(n3+n4+n8,3)+4.*Q(n2,1)*Q(n5,1)*Q(n1+n6+n7,3)*Q(n3+n4+n8,3)
      - 4.*Q(n2+n5,2)*Q(n1+n6+n7,3)*Q(n3+n4+n8,3)+4.*Q(n1,1)*Q(n5,1)*Q(n2+n6+n7,3)*Q(n3+n4+n8,3)
      - 4.*Q(n1+n5,2)*Q(n2+n6+n7,3)*Q(n3+n4+n8,3)-12.*Q(n5,1)*Q(n1+n2+n6+n7,4)*Q(n3+n4+n8,3)
      + 4.*Q(n1,1)*Q(n2,1)*Q(n5+n6+n7,3)*Q(n3+n4+n8,3)-4.*Q(n1+n2,2)*Q(n5+n6+n7,3)*Q(n3+n4+n8,3)
      - 12.*Q(n2,1)*Q(n1+n5+n6+n7,4)*Q(n3+n4+n8,3)-12.*Q(n1,1)*Q(n2+n5+n6+n7,4)*Q(n3+n4+n8,3)
      + 48.*Q(n1+n2+n5+n6+n7,5)*Q(n3+n4+n8,3)-6.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n4+n8,4)
      + 6.*Q(n2+n5,2)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n4+n8,4)+6.*Q(n5,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n3+n4+n8,4)
      + 6.*Q(n2,1)*Q(n5+n6,2)*Q(n7,1)*Q(n1+n3+n4+n8,4)-12.*Q(n2+n5+n6,3)*Q(n7,1)*Q(n1+n3+n4+n8,4)
      + 6.*Q(n5,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n3+n4+n8,4)-6.*Q(n5+n6,2)*Q(n2+n7,2)*Q(n1+n3+n4+n8,4)
      + 6.*Q(n2,1)*Q(n6,1)*Q(n5+n7,2)*Q(n1+n3+n4+n8,4)-6.*Q(n2+n6,2)*Q(n5+n7,2)*Q(n1+n3+n4+n8,4)
      - 12.*Q(n6,1)*Q(n2+n5+n7,3)*Q(n1+n3+n4+n8,4)+6.*Q(n2,1)*Q(n5,1)*Q(n6+n7,2)*Q(n1+n3+n4+n8,4)
      - 6.*Q(n2+n5,2)*Q(n6+n7,2)*Q(n1+n3+n4+n8,4)-12.*Q(n5,1)*Q(n2+n6+n7,3)*Q(n1+n3+n4+n8,4)
      - 12.*Q(n2,1)*Q(n5+n6+n7,3)*Q(n1+n3+n4+n8,4)+36.*Q(n2+n5+n6+n7,4)*Q(n1+n3+n4+n8,4)
      - 6.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n4+n8,4)+6.*Q(n1+n5,2)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n4+n8,4)
      + 6.*Q(n5,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n3+n4+n8,4)+6.*Q(n1,1)*Q(n5+n6,2)*Q(n7,1)*Q(n2+n3+n4+n8,4)
      - 12.*Q(n1+n5+n6,3)*Q(n7,1)*Q(n2+n3+n4+n8,4)+6.*Q(n5,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n3+n4+n8,4)
      - 6.*Q(n5+n6,2)*Q(n1+n7,2)*Q(n2+n3+n4+n8,4)+6.*Q(n1,1)*Q(n6,1)*Q(n5+n7,2)*Q(n2+n3+n4+n8,4)
      - 6.*Q(n1+n6,2)*Q(n5+n7,2)*Q(n2+n3+n4+n8,4)-12.*Q(n6,1)*Q(n1+n5+n7,3)*Q(n2+n3+n4+n8,4)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n6+n7,2)*Q(n2+n3+n4+n8,4)-6.*Q(n1+n5,2)*Q(n6+n7,2)*Q(n2+n3+n4+n8,4)
      - 12.*Q(n5,1)*Q(n1+n6+n7,3)*Q(n2+n3+n4+n8,4)-12.*Q(n1,1)*Q(n5+n6+n7,3)*Q(n2+n3+n4+n8,4)
      + 36.*Q(n1+n5+n6+n7,4)*Q(n2+n3+n4+n8,4)+24.*Q(n5,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n3+n4+n8,5)
      - 24.*Q(n5+n6,2)*Q(n7,1)*Q(n1+n2+n3+n4+n8,5)-24.*Q(n6,1)*Q(n5+n7,2)*Q(n1+n2+n3+n4+n8,5)
      - 24.*Q(n5,1)*Q(n6+n7,2)*Q(n1+n2+n3+n4+n8,5)+48.*Q(n5+n6+n7,3)*Q(n1+n2+n3+n4+n8,5)
      - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)+Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      + Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)+Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      - 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)+Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      - Q(n2+n3,2)*Q(n1+n4,2)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)+Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      - Q(n1+n3,2)*Q(n2+n4,2)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)-2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)-Q(n1+n2,2)*Q(n3+n4,2)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      - 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)-2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)
      + 6.*Q(n1+n2+n3+n4,4)*Q(n6,1)*Q(n7,1)*Q(n5+n8,2)+Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n2+n3,2)*Q(n4,1)*Q(n1+n6,2)*Q(n7,1)*Q(n5+n8,2)-Q(n3,1)*Q(n2+n4,2)*Q(n1+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n2,1)*Q(n3+n4,2)*Q(n1+n6,2)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n2+n3+n4,3)*Q(n1+n6,2)*Q(n7,1)*Q(n5+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n7,1)*Q(n5+n8,2)-Q(n1+n3,2)*Q(n4,1)*Q(n2+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n3,1)*Q(n1+n4,2)*Q(n2+n6,2)*Q(n7,1)*Q(n5+n8,2)-Q(n1,1)*Q(n3+n4,2)*Q(n2+n6,2)*Q(n7,1)*Q(n5+n8,2)
      + 2.*Q(n1+n3+n4,3)*Q(n2+n6,2)*Q(n7,1)*Q(n5+n8,2)-2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n5+n8,2)
      + 2.*Q(n3+n4,2)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n5+n8,2)+Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n1+n2,2)*Q(n4,1)*Q(n3+n6,2)*Q(n7,1)*Q(n5+n8,2)-Q(n2,1)*Q(n1+n4,2)*Q(n3+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n1,1)*Q(n2+n4,2)*Q(n3+n6,2)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n1+n2+n4,3)*Q(n3+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n2+n4,2)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n5+n8,2)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n1+n4,2)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n5+n8,2)
      + 6.*Q(n4,1)*Q(n1+n2+n3+n6,4)*Q(n7,1)*Q(n5+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n1+n2,2)*Q(n3,1)*Q(n4+n6,2)*Q(n7,1)*Q(n5+n8,2)-Q(n2,1)*Q(n1+n3,2)*Q(n4+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - Q(n1,1)*Q(n2+n3,2)*Q(n4+n6,2)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n1+n2+n3,3)*Q(n4+n6,2)*Q(n7,1)*Q(n5+n8,2)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n2+n3,2)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n5+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n5+n8,2)+2.*Q(n1+n3,2)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n5+n8,2)
      + 6.*Q(n3,1)*Q(n1+n2+n4+n6,4)*Q(n7,1)*Q(n5+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n5+n8,2)
      + 2.*Q(n1+n2,2)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n5+n8,2)+6.*Q(n2,1)*Q(n1+n3+n4+n6,4)*Q(n7,1)*Q(n5+n8,2)
      + 6.*Q(n1,1)*Q(n2+n3+n4+n6,4)*Q(n7,1)*Q(n5+n8,2)-24.*Q(n1+n2+n3+n4+n6,5)*Q(n7,1)*Q(n5+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n7,2)*Q(n5+n8,2)-Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n1+n7,2)*Q(n5+n8,2)
      - Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n1+n7,2)*Q(n5+n8,2)-Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n1+n7,2)*Q(n5+n8,2)
      + 2.*Q(n2+n3+n4,3)*Q(n6,1)*Q(n1+n7,2)*Q(n5+n8,2)-Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n5+n8,2)
      + Q(n3+n4,2)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n5+n8,2)-Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n5+n8,2)
      + Q(n2+n4,2)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n5+n8,2)+2.*Q(n4,1)*Q(n2+n3+n6,3)*Q(n1+n7,2)*Q(n5+n8,2)
      - Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n5+n8,2)+Q(n2+n3,2)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n5+n8,2)
      + 2.*Q(n3,1)*Q(n2+n4+n6,3)*Q(n1+n7,2)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n3+n4+n6,3)*Q(n1+n7,2)*Q(n5+n8,2)
      - 6.*Q(n2+n3+n4+n6,4)*Q(n1+n7,2)*Q(n5+n8,2)+Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n2+n7,2)*Q(n5+n8,2)
      - Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n2+n7,2)*Q(n5+n8,2)-Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n2+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n3+n4,2)*Q(n6,1)*Q(n2+n7,2)*Q(n5+n8,2)+2.*Q(n1+n3+n4,3)*Q(n6,1)*Q(n2+n7,2)*Q(n5+n8,2)
      - Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n5+n8,2)+Q(n3+n4,2)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n4,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n5+n8,2)+Q(n1+n4,2)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n5+n8,2)
      + 2.*Q(n4,1)*Q(n1+n3+n6,3)*Q(n2+n7,2)*Q(n5+n8,2)-Q(n1,1)*Q(n3,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n5+n8,2)
      + Q(n1+n3,2)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n5+n8,2)+2.*Q(n3,1)*Q(n1+n4+n6,3)*Q(n2+n7,2)*Q(n5+n8,2)
      + 2.*Q(n1,1)*Q(n3+n4+n6,3)*Q(n2+n7,2)*Q(n5+n8,2)-6.*Q(n1+n3+n4+n6,4)*Q(n2+n7,2)*Q(n5+n8,2)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n5+n8,2)+2.*Q(n3+n4,2)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n5+n8,2)
      + 2.*Q(n4,1)*Q(n3+n6,2)*Q(n1+n2+n7,3)*Q(n5+n8,2)+2.*Q(n3,1)*Q(n4+n6,2)*Q(n1+n2+n7,3)*Q(n5+n8,2)
      - 4.*Q(n3+n4+n6,3)*Q(n1+n2+n7,3)*Q(n5+n8,2)+Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n3+n7,2)*Q(n5+n8,2)
      - Q(n1+n2,2)*Q(n4,1)*Q(n6,1)*Q(n3+n7,2)*Q(n5+n8,2)-Q(n2,1)*Q(n1+n4,2)*Q(n6,1)*Q(n3+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n2+n4,2)*Q(n6,1)*Q(n3+n7,2)*Q(n5+n8,2)+2.*Q(n1+n2+n4,3)*Q(n6,1)*Q(n3+n7,2)*Q(n5+n8,2)
      - Q(n2,1)*Q(n4,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n5+n8,2)+Q(n2+n4,2)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n4,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n5+n8,2)+Q(n1+n4,2)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n5+n8,2)
      + 2.*Q(n4,1)*Q(n1+n2+n6,3)*Q(n3+n7,2)*Q(n5+n8,2)-Q(n1,1)*Q(n2,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n5+n8,2)
      + Q(n1+n2,2)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n1+n4+n6,3)*Q(n3+n7,2)*Q(n5+n8,2)
      + 2.*Q(n1,1)*Q(n2+n4+n6,3)*Q(n3+n7,2)*Q(n5+n8,2)-6.*Q(n1+n2+n4+n6,4)*Q(n3+n7,2)*Q(n5+n8,2)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n5+n8,2)+2.*Q(n2+n4,2)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n5+n8,2)
      + 2.*Q(n4,1)*Q(n2+n6,2)*Q(n1+n3+n7,3)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n4+n6,2)*Q(n1+n3+n7,3)*Q(n5+n8,2)
      - 4.*Q(n2+n4+n6,3)*Q(n1+n3+n7,3)*Q(n5+n8,2)-2.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n5+n8,2)
      + 2.*Q(n1+n4,2)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n5+n8,2)+2.*Q(n4,1)*Q(n1+n6,2)*Q(n2+n3+n7,3)*Q(n5+n8,2)
      + 2.*Q(n1,1)*Q(n4+n6,2)*Q(n2+n3+n7,3)*Q(n5+n8,2)-4.*Q(n1+n4+n6,3)*Q(n2+n3+n7,3)*Q(n5+n8,2)
      + 6.*Q(n4,1)*Q(n6,1)*Q(n1+n2+n3+n7,4)*Q(n5+n8,2)-6.*Q(n4+n6,2)*Q(n1+n2+n3+n7,4)*Q(n5+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n4+n7,2)*Q(n5+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n6,1)*Q(n4+n7,2)*Q(n5+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n6,1)*Q(n4+n7,2)*Q(n5+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n6,1)*Q(n4+n7,2)*Q(n5+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n6,1)*Q(n4+n7,2)*Q(n5+n8,2)-Q(n2,1)*Q(n3,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n5+n8,2)
      + Q(n2+n3,2)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n5+n8,2)-Q(n1,1)*Q(n3,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n5+n8,2)
      + Q(n1+n3,2)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n5+n8,2)+2.*Q(n3,1)*Q(n1+n2+n6,3)*Q(n4+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n5+n8,2)+Q(n1+n2,2)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n5+n8,2)
      + 2.*Q(n2,1)*Q(n1+n3+n6,3)*Q(n4+n7,2)*Q(n5+n8,2)+2.*Q(n1,1)*Q(n2+n3+n6,3)*Q(n4+n7,2)*Q(n5+n8,2)
      - 6.*Q(n1+n2+n3+n6,4)*Q(n4+n7,2)*Q(n5+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n5+n8,2)
      + 2.*Q(n2+n3,2)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n5+n8,2)+2.*Q(n3,1)*Q(n2+n6,2)*Q(n1+n4+n7,3)*Q(n5+n8,2)
      + 2.*Q(n2,1)*Q(n3+n6,2)*Q(n1+n4+n7,3)*Q(n5+n8,2)-4.*Q(n2+n3+n6,3)*Q(n1+n4+n7,3)*Q(n5+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n5+n8,2)+2.*Q(n1+n3,2)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n5+n8,2)
      + 2.*Q(n3,1)*Q(n1+n6,2)*Q(n2+n4+n7,3)*Q(n5+n8,2)+2.*Q(n1,1)*Q(n3+n6,2)*Q(n2+n4+n7,3)*Q(n5+n8,2)
      - 4.*Q(n1+n3+n6,3)*Q(n2+n4+n7,3)*Q(n5+n8,2)+6.*Q(n3,1)*Q(n6,1)*Q(n1+n2+n4+n7,4)*Q(n5+n8,2)
      - 6.*Q(n3+n6,2)*Q(n1+n2+n4+n7,4)*Q(n5+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n5+n8,2)
      + 2.*Q(n1+n2,2)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n1+n6,2)*Q(n3+n4+n7,3)*Q(n5+n8,2)
      + 2.*Q(n1,1)*Q(n2+n6,2)*Q(n3+n4+n7,3)*Q(n5+n8,2)-4.*Q(n1+n2+n6,3)*Q(n3+n4+n7,3)*Q(n5+n8,2)
      + 6.*Q(n2,1)*Q(n6,1)*Q(n1+n3+n4+n7,4)*Q(n5+n8,2)-6.*Q(n2+n6,2)*Q(n1+n3+n4+n7,4)*Q(n5+n8,2)
      + 6.*Q(n1,1)*Q(n6,1)*Q(n2+n3+n4+n7,4)*Q(n5+n8,2)-6.*Q(n1+n6,2)*Q(n2+n3+n4+n7,4)*Q(n5+n8,2)
      - 24.*Q(n6,1)*Q(n1+n2+n3+n4+n7,5)*Q(n5+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6+n7,2)*Q(n5+n8,2)
      - Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n6+n7,2)*Q(n5+n8,2)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n6+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n6+n7,2)*Q(n5+n8,2)+2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n6+n7,2)*Q(n5+n8,2)
      - Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n6+n7,2)*Q(n5+n8,2)+Q(n2+n3,2)*Q(n1+n4,2)*Q(n6+n7,2)*Q(n5+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n6+n7,2)*Q(n5+n8,2)+Q(n1+n3,2)*Q(n2+n4,2)*Q(n6+n7,2)*Q(n5+n8,2)
      + 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n6+n7,2)*Q(n5+n8,2)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n6+n7,2)*Q(n5+n8,2)
      + Q(n1+n2,2)*Q(n3+n4,2)*Q(n6+n7,2)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n6+n7,2)*Q(n5+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n6+n7,2)*Q(n5+n8,2)-6.*Q(n1+n2+n3+n4,4)*Q(n6+n7,2)*Q(n5+n8,2)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n6+n7,3)*Q(n5+n8,2)+2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n6+n7,3)*Q(n5+n8,2)
      + 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n6+n7,3)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n6+n7,3)*Q(n5+n8,2)
      - 4.*Q(n2+n3+n4,3)*Q(n1+n6+n7,3)*Q(n5+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n6+n7,3)*Q(n5+n8,2)
      + 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n6+n7,3)*Q(n5+n8,2)+2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n6+n7,3)*Q(n5+n8,2)
      + 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n6+n7,3)*Q(n5+n8,2)-4.*Q(n1+n3+n4,3)*Q(n2+n6+n7,3)*Q(n5+n8,2)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n6+n7,4)*Q(n5+n8,2)-6.*Q(n3+n4,2)*Q(n1+n2+n6+n7,4)*Q(n5+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n6+n7,3)*Q(n5+n8,2)+2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n6+n7,3)*Q(n5+n8,2)
      + 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n6+n7,3)*Q(n5+n8,2)+2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n6+n7,3)*Q(n5+n8,2)
      - 4.*Q(n1+n2+n4,3)*Q(n3+n6+n7,3)*Q(n5+n8,2)+6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n6+n7,4)*Q(n5+n8,2)
      - 6.*Q(n2+n4,2)*Q(n1+n3+n6+n7,4)*Q(n5+n8,2)+6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n6+n7,4)*Q(n5+n8,2)
      - 6.*Q(n1+n4,2)*Q(n2+n3+n6+n7,4)*Q(n5+n8,2)-24.*Q(n4,1)*Q(n1+n2+n3+n6+n7,5)*Q(n5+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n6+n7,3)*Q(n5+n8,2)+2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n6+n7,3)*Q(n5+n8,2)
      + 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n6+n7,3)*Q(n5+n8,2)+2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n6+n7,3)*Q(n5+n8,2)
      - 4.*Q(n1+n2+n3,3)*Q(n4+n6+n7,3)*Q(n5+n8,2)+6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n6+n7,4)*Q(n5+n8,2)
      - 6.*Q(n2+n3,2)*Q(n1+n4+n6+n7,4)*Q(n5+n8,2)+6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n6+n7,4)*Q(n5+n8,2)
      - 6.*Q(n1+n3,2)*Q(n2+n4+n6+n7,4)*Q(n5+n8,2)-24.*Q(n3,1)*Q(n1+n2+n4+n6+n7,5)*Q(n5+n8,2)
      + 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n6+n7,4)*Q(n5+n8,2)-6.*Q(n1+n2,2)*Q(n3+n4+n6+n7,4)*Q(n5+n8,2)
      - 24.*Q(n2,1)*Q(n1+n3+n4+n6+n7,5)*Q(n5+n8,2)-24.*Q(n1,1)*Q(n2+n3+n4+n6+n7,5)*Q(n5+n8,2)
      + 120.*Q(n1+n2+n3+n4+n6+n7,6)*Q(n5+n8,2)+2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n1+n5+n8,3)
      - 2.*Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n1+n5+n8,3)-2.*Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n7,1)*Q(n1+n5+n8,3)
      - 2.*Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n7,1)*Q(n1+n5+n8,3)+4.*Q(n2+n3+n4,3)*Q(n6,1)*Q(n7,1)*Q(n1+n5+n8,3)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n5+n8,3)+2.*Q(n3+n4,2)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n5+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n5+n8,3)+2.*Q(n2+n4,2)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n5+n8,3)
      + 4.*Q(n4,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n1+n5+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n5+n8,3)
      + 2.*Q(n2+n3,2)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n5+n8,3)+4.*Q(n3,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n1+n5+n8,3)
      + 4.*Q(n2,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n1+n5+n8,3)-12.*Q(n2+n3+n4+n6,4)*Q(n7,1)*Q(n1+n5+n8,3)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n5+n8,3)+2.*Q(n3+n4,2)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n5+n8,3)
      + 2.*Q(n4,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n1+n5+n8,3)+2.*Q(n3,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n1+n5+n8,3)
      - 4.*Q(n3+n4+n6,3)*Q(n2+n7,2)*Q(n1+n5+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n5+n8,3)
      + 2.*Q(n2+n4,2)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n5+n8,3)+2.*Q(n4,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n1+n5+n8,3)
      + 2.*Q(n2,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n1+n5+n8,3)-4.*Q(n2+n4+n6,3)*Q(n3+n7,2)*Q(n1+n5+n8,3)
      + 4.*Q(n4,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n1+n5+n8,3)-4.*Q(n4+n6,2)*Q(n2+n3+n7,3)*Q(n1+n5+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n5+n8,3)+2.*Q(n2+n3,2)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n5+n8,3)
      + 2.*Q(n3,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n1+n5+n8,3)+2.*Q(n2,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n1+n5+n8,3)
      - 4.*Q(n2+n3+n6,3)*Q(n4+n7,2)*Q(n1+n5+n8,3)+4.*Q(n3,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n1+n5+n8,3)
      - 4.*Q(n3+n6,2)*Q(n2+n4+n7,3)*Q(n1+n5+n8,3)+4.*Q(n2,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n1+n5+n8,3)
      - 4.*Q(n2+n6,2)*Q(n3+n4+n7,3)*Q(n1+n5+n8,3)-12.*Q(n6,1)*Q(n2+n3+n4+n7,4)*Q(n1+n5+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6+n7,2)*Q(n1+n5+n8,3)+2.*Q(n2+n3,2)*Q(n4,1)*Q(n6+n7,2)*Q(n1+n5+n8,3)
      + 2.*Q(n3,1)*Q(n2+n4,2)*Q(n6+n7,2)*Q(n1+n5+n8,3)+2.*Q(n2,1)*Q(n3+n4,2)*Q(n6+n7,2)*Q(n1+n5+n8,3)
      - 4.*Q(n2+n3+n4,3)*Q(n6+n7,2)*Q(n1+n5+n8,3)+4.*Q(n3,1)*Q(n4,1)*Q(n2+n6+n7,3)*Q(n1+n5+n8,3)
      - 4.*Q(n3+n4,2)*Q(n2+n6+n7,3)*Q(n1+n5+n8,3)+4.*Q(n2,1)*Q(n4,1)*Q(n3+n6+n7,3)*Q(n1+n5+n8,3)
      - 4.*Q(n2+n4,2)*Q(n3+n6+n7,3)*Q(n1+n5+n8,3)-12.*Q(n4,1)*Q(n2+n3+n6+n7,4)*Q(n1+n5+n8,3)
      + 4.*Q(n2,1)*Q(n3,1)*Q(n4+n6+n7,3)*Q(n1+n5+n8,3)-4.*Q(n2+n3,2)*Q(n4+n6+n7,3)*Q(n1+n5+n8,3)
      - 12.*Q(n3,1)*Q(n2+n4+n6+n7,4)*Q(n1+n5+n8,3)-12.*Q(n2,1)*Q(n3+n4+n6+n7,4)*Q(n1+n5+n8,3)
      + 48.*Q(n2+n3+n4+n6+n7,5)*Q(n1+n5+n8,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n2+n5+n8,3)
      - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n2+n5+n8,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n7,1)*Q(n2+n5+n8,3)
      - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n6,1)*Q(n7,1)*Q(n2+n5+n8,3)+4.*Q(n1+n3+n4,3)*Q(n6,1)*Q(n7,1)*Q(n2+n5+n8,3)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n5+n8,3)+2.*Q(n3+n4,2)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n5+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n5+n8,3)+2.*Q(n1+n4,2)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n5+n8,3)
      + 4.*Q(n4,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n2+n5+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n5+n8,3)
      + 2.*Q(n1+n3,2)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n5+n8,3)+4.*Q(n3,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n2+n5+n8,3)
      + 4.*Q(n1,1)*Q(n3+n4+n6,3)*Q(n7,1)*Q(n2+n5+n8,3)-12.*Q(n1+n3+n4+n6,4)*Q(n7,1)*Q(n2+n5+n8,3)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n5+n8,3)+2.*Q(n3+n4,2)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n5+n8,3)
      + 2.*Q(n4,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n2+n5+n8,3)+2.*Q(n3,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n2+n5+n8,3)
      - 4.*Q(n3+n4+n6,3)*Q(n1+n7,2)*Q(n2+n5+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n5+n8,3)
      + 2.*Q(n1+n4,2)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n5+n8,3)+2.*Q(n4,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n2+n5+n8,3)
      + 2.*Q(n1,1)*Q(n4+n6,2)*Q(n3+n7,2)*Q(n2+n5+n8,3)-4.*Q(n1+n4+n6,3)*Q(n3+n7,2)*Q(n2+n5+n8,3)
      + 4.*Q(n4,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n2+n5+n8,3)-4.*Q(n4+n6,2)*Q(n1+n3+n7,3)*Q(n2+n5+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n5+n8,3)+2.*Q(n1+n3,2)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n5+n8,3)
      + 2.*Q(n3,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n2+n5+n8,3)+2.*Q(n1,1)*Q(n3+n6,2)*Q(n4+n7,2)*Q(n2+n5+n8,3)
      - 4.*Q(n1+n3+n6,3)*Q(n4+n7,2)*Q(n2+n5+n8,3)+4.*Q(n3,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n2+n5+n8,3)
      - 4.*Q(n3+n6,2)*Q(n1+n4+n7,3)*Q(n2+n5+n8,3)+4.*Q(n1,1)*Q(n6,1)*Q(n3+n4+n7,3)*Q(n2+n5+n8,3)
      - 4.*Q(n1+n6,2)*Q(n3+n4+n7,3)*Q(n2+n5+n8,3)-12.*Q(n6,1)*Q(n1+n3+n4+n7,4)*Q(n2+n5+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n6+n7,2)*Q(n2+n5+n8,3)+2.*Q(n1+n3,2)*Q(n4,1)*Q(n6+n7,2)*Q(n2+n5+n8,3)
      + 2.*Q(n3,1)*Q(n1+n4,2)*Q(n6+n7,2)*Q(n2+n5+n8,3)+2.*Q(n1,1)*Q(n3+n4,2)*Q(n6+n7,2)*Q(n2+n5+n8,3)
      - 4.*Q(n1+n3+n4,3)*Q(n6+n7,2)*Q(n2+n5+n8,3)+4.*Q(n3,1)*Q(n4,1)*Q(n1+n6+n7,3)*Q(n2+n5+n8,3)
      - 4.*Q(n3+n4,2)*Q(n1+n6+n7,3)*Q(n2+n5+n8,3)+4.*Q(n1,1)*Q(n4,1)*Q(n3+n6+n7,3)*Q(n2+n5+n8,3)
      - 4.*Q(n1+n4,2)*Q(n3+n6+n7,3)*Q(n2+n5+n8,3)-12.*Q(n4,1)*Q(n1+n3+n6+n7,4)*Q(n2+n5+n8,3)
      + 4.*Q(n1,1)*Q(n3,1)*Q(n4+n6+n7,3)*Q(n2+n5+n8,3)-4.*Q(n1+n3,2)*Q(n4+n6+n7,3)*Q(n2+n5+n8,3)
      - 12.*Q(n3,1)*Q(n1+n4+n6+n7,4)*Q(n2+n5+n8,3)-12.*Q(n1,1)*Q(n3+n4+n6+n7,4)*Q(n2+n5+n8,3)
      + 48.*Q(n1+n3+n4+n6+n7,5)*Q(n2+n5+n8,3)-6.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n5+n8,4)
      + 6.*Q(n3+n4,2)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n5+n8,4)+6.*Q(n4,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n2+n5+n8,4)
      + 6.*Q(n3,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n2+n5+n8,4)-12.*Q(n3+n4+n6,3)*Q(n7,1)*Q(n1+n2+n5+n8,4)
      + 6.*Q(n4,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n2+n5+n8,4)-6.*Q(n4+n6,2)*Q(n3+n7,2)*Q(n1+n2+n5+n8,4)
      + 6.*Q(n3,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n2+n5+n8,4)-6.*Q(n3+n6,2)*Q(n4+n7,2)*Q(n1+n2+n5+n8,4)
      - 12.*Q(n6,1)*Q(n3+n4+n7,3)*Q(n1+n2+n5+n8,4)+6.*Q(n3,1)*Q(n4,1)*Q(n6+n7,2)*Q(n1+n2+n5+n8,4)
      - 6.*Q(n3+n4,2)*Q(n6+n7,2)*Q(n1+n2+n5+n8,4)-12.*Q(n4,1)*Q(n3+n6+n7,3)*Q(n1+n2+n5+n8,4)
      - 12.*Q(n3,1)*Q(n4+n6+n7,3)*Q(n1+n2+n5+n8,4)+36.*Q(n3+n4+n6+n7,4)*Q(n1+n2+n5+n8,4)
      + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n3+n5+n8,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n3+n5+n8,3)
      - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n6,1)*Q(n7,1)*Q(n3+n5+n8,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n6,1)*Q(n7,1)*Q(n3+n5+n8,3)
      + 4.*Q(n1+n2+n4,3)*Q(n6,1)*Q(n7,1)*Q(n3+n5+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n5+n8,3)
      + 2.*Q(n2+n4,2)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n5+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n5+n8,3)
      + 2.*Q(n1+n4,2)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n5+n8,3)+4.*Q(n4,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n3+n5+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n5+n8,3)+2.*Q(n1+n2,2)*Q(n4+n6,2)*Q(n7,1)*Q(n3+n5+n8,3)
      + 4.*Q(n2,1)*Q(n1+n4+n6,3)*Q(n7,1)*Q(n3+n5+n8,3)+4.*Q(n1,1)*Q(n2+n4+n6,3)*Q(n7,1)*Q(n3+n5+n8,3)
      - 12.*Q(n1+n2+n4+n6,4)*Q(n7,1)*Q(n3+n5+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n2+n4,2)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n5+n8,3)+2.*Q(n4,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n2,1)*Q(n4+n6,2)*Q(n1+n7,2)*Q(n3+n5+n8,3)-4.*Q(n2+n4+n6,3)*Q(n1+n7,2)*Q(n3+n5+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n5+n8,3)+2.*Q(n1+n4,2)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n4,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n3+n5+n8,3)+2.*Q(n1,1)*Q(n4+n6,2)*Q(n2+n7,2)*Q(n3+n5+n8,3)
      - 4.*Q(n1+n4+n6,3)*Q(n2+n7,2)*Q(n3+n5+n8,3)+4.*Q(n4,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n3+n5+n8,3)
      - 4.*Q(n4+n6,2)*Q(n1+n2+n7,3)*Q(n3+n5+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n1+n2,2)*Q(n6,1)*Q(n4+n7,2)*Q(n3+n5+n8,3)+2.*Q(n2,1)*Q(n1+n6,2)*Q(n4+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n1,1)*Q(n2+n6,2)*Q(n4+n7,2)*Q(n3+n5+n8,3)-4.*Q(n1+n2+n6,3)*Q(n4+n7,2)*Q(n3+n5+n8,3)
      + 4.*Q(n2,1)*Q(n6,1)*Q(n1+n4+n7,3)*Q(n3+n5+n8,3)-4.*Q(n2+n6,2)*Q(n1+n4+n7,3)*Q(n3+n5+n8,3)
      + 4.*Q(n1,1)*Q(n6,1)*Q(n2+n4+n7,3)*Q(n3+n5+n8,3)-4.*Q(n1+n6,2)*Q(n2+n4+n7,3)*Q(n3+n5+n8,3)
      - 12.*Q(n6,1)*Q(n1+n2+n4+n7,4)*Q(n3+n5+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n6+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n1+n2,2)*Q(n4,1)*Q(n6+n7,2)*Q(n3+n5+n8,3)+2.*Q(n2,1)*Q(n1+n4,2)*Q(n6+n7,2)*Q(n3+n5+n8,3)
      + 2.*Q(n1,1)*Q(n2+n4,2)*Q(n6+n7,2)*Q(n3+n5+n8,3)-4.*Q(n1+n2+n4,3)*Q(n6+n7,2)*Q(n3+n5+n8,3)
      + 4.*Q(n2,1)*Q(n4,1)*Q(n1+n6+n7,3)*Q(n3+n5+n8,3)-4.*Q(n2+n4,2)*Q(n1+n6+n7,3)*Q(n3+n5+n8,3)
      + 4.*Q(n1,1)*Q(n4,1)*Q(n2+n6+n7,3)*Q(n3+n5+n8,3)-4.*Q(n1+n4,2)*Q(n2+n6+n7,3)*Q(n3+n5+n8,3)
      - 12.*Q(n4,1)*Q(n1+n2+n6+n7,4)*Q(n3+n5+n8,3)+4.*Q(n1,1)*Q(n2,1)*Q(n4+n6+n7,3)*Q(n3+n5+n8,3)
      - 4.*Q(n1+n2,2)*Q(n4+n6+n7,3)*Q(n3+n5+n8,3)-12.*Q(n2,1)*Q(n1+n4+n6+n7,4)*Q(n3+n5+n8,3)
      - 12.*Q(n1,1)*Q(n2+n4+n6+n7,4)*Q(n3+n5+n8,3)+48.*Q(n1+n2+n4+n6+n7,5)*Q(n3+n5+n8,3)
      - 6.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n5+n8,4)+6.*Q(n2+n4,2)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n5+n8,4)
      + 6.*Q(n4,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n3+n5+n8,4)+6.*Q(n2,1)*Q(n4+n6,2)*Q(n7,1)*Q(n1+n3+n5+n8,4)
      - 12.*Q(n2+n4+n6,3)*Q(n7,1)*Q(n1+n3+n5+n8,4)+6.*Q(n4,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n3+n5+n8,4)
      - 6.*Q(n4+n6,2)*Q(n2+n7,2)*Q(n1+n3+n5+n8,4)+6.*Q(n2,1)*Q(n6,1)*Q(n4+n7,2)*Q(n1+n3+n5+n8,4)
      - 6.*Q(n2+n6,2)*Q(n4+n7,2)*Q(n1+n3+n5+n8,4)-12.*Q(n6,1)*Q(n2+n4+n7,3)*Q(n1+n3+n5+n8,4)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n6+n7,2)*Q(n1+n3+n5+n8,4)-6.*Q(n2+n4,2)*Q(n6+n7,2)*Q(n1+n3+n5+n8,4)
      - 12.*Q(n4,1)*Q(n2+n6+n7,3)*Q(n1+n3+n5+n8,4)-12.*Q(n2,1)*Q(n4+n6+n7,3)*Q(n1+n3+n5+n8,4)
      + 36.*Q(n2+n4+n6+n7,4)*Q(n1+n3+n5+n8,4)-6.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n5+n8,4)
      + 6.*Q(n1+n4,2)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n5+n8,4)+6.*Q(n4,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n3+n5+n8,4)
      + 6.*Q(n1,1)*Q(n4+n6,2)*Q(n7,1)*Q(n2+n3+n5+n8,4)-12.*Q(n1+n4+n6,3)*Q(n7,1)*Q(n2+n3+n5+n8,4)
      + 6.*Q(n4,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n3+n5+n8,4)-6.*Q(n4+n6,2)*Q(n1+n7,2)*Q(n2+n3+n5+n8,4)
      + 6.*Q(n1,1)*Q(n6,1)*Q(n4+n7,2)*Q(n2+n3+n5+n8,4)-6.*Q(n1+n6,2)*Q(n4+n7,2)*Q(n2+n3+n5+n8,4)
      - 12.*Q(n6,1)*Q(n1+n4+n7,3)*Q(n2+n3+n5+n8,4)+6.*Q(n1,1)*Q(n4,1)*Q(n6+n7,2)*Q(n2+n3+n5+n8,4)
      - 6.*Q(n1+n4,2)*Q(n6+n7,2)*Q(n2+n3+n5+n8,4)-12.*Q(n4,1)*Q(n1+n6+n7,3)*Q(n2+n3+n5+n8,4)
      - 12.*Q(n1,1)*Q(n4+n6+n7,3)*Q(n2+n3+n5+n8,4)+36.*Q(n1+n4+n6+n7,4)*Q(n2+n3+n5+n8,4)
      + 24.*Q(n4,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n3+n5+n8,5)-24.*Q(n4+n6,2)*Q(n7,1)*Q(n1+n2+n3+n5+n8,5)
      - 24.*Q(n6,1)*Q(n4+n7,2)*Q(n1+n2+n3+n5+n8,5)-24.*Q(n4,1)*Q(n6+n7,2)*Q(n1+n2+n3+n5+n8,5)
      + 48.*Q(n4+n6+n7,3)*Q(n1+n2+n3+n5+n8,5)+2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n7,1)*Q(n4+n5+n8,3)
      - 2.*Q(n1+n2,2)*Q(n3,1)*Q(n6,1)*Q(n7,1)*Q(n4+n5+n8,3)-2.*Q(n2,1)*Q(n1+n3,2)*Q(n6,1)*Q(n7,1)*Q(n4+n5+n8,3)
      - 2.*Q(n1,1)*Q(n2+n3,2)*Q(n6,1)*Q(n7,1)*Q(n4+n5+n8,3)+4.*Q(n1+n2+n3,3)*Q(n6,1)*Q(n7,1)*Q(n4+n5+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n5+n8,3)+2.*Q(n2+n3,2)*Q(n1+n6,2)*Q(n7,1)*Q(n4+n5+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n5+n8,3)+2.*Q(n1+n3,2)*Q(n2+n6,2)*Q(n7,1)*Q(n4+n5+n8,3)
      + 4.*Q(n3,1)*Q(n1+n2+n6,3)*Q(n7,1)*Q(n4+n5+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n5+n8,3)
      + 2.*Q(n1+n2,2)*Q(n3+n6,2)*Q(n7,1)*Q(n4+n5+n8,3)+4.*Q(n2,1)*Q(n1+n3+n6,3)*Q(n7,1)*Q(n4+n5+n8,3)
      + 4.*Q(n1,1)*Q(n2+n3+n6,3)*Q(n7,1)*Q(n4+n5+n8,3)-12.*Q(n1+n2+n3+n6,4)*Q(n7,1)*Q(n4+n5+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n5+n8,3)+2.*Q(n2+n3,2)*Q(n6,1)*Q(n1+n7,2)*Q(n4+n5+n8,3)
      + 2.*Q(n3,1)*Q(n2+n6,2)*Q(n1+n7,2)*Q(n4+n5+n8,3)+2.*Q(n2,1)*Q(n3+n6,2)*Q(n1+n7,2)*Q(n4+n5+n8,3)
      - 4.*Q(n2+n3+n6,3)*Q(n1+n7,2)*Q(n4+n5+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n5+n8,3)
      + 2.*Q(n1+n3,2)*Q(n6,1)*Q(n2+n7,2)*Q(n4+n5+n8,3)+2.*Q(n3,1)*Q(n1+n6,2)*Q(n2+n7,2)*Q(n4+n5+n8,3)
      + 2.*Q(n1,1)*Q(n3+n6,2)*Q(n2+n7,2)*Q(n4+n5+n8,3)-4.*Q(n1+n3+n6,3)*Q(n2+n7,2)*Q(n4+n5+n8,3)
      + 4.*Q(n3,1)*Q(n6,1)*Q(n1+n2+n7,3)*Q(n4+n5+n8,3)-4.*Q(n3+n6,2)*Q(n1+n2+n7,3)*Q(n4+n5+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n5+n8,3)+2.*Q(n1+n2,2)*Q(n6,1)*Q(n3+n7,2)*Q(n4+n5+n8,3)
      + 2.*Q(n2,1)*Q(n1+n6,2)*Q(n3+n7,2)*Q(n4+n5+n8,3)+2.*Q(n1,1)*Q(n2+n6,2)*Q(n3+n7,2)*Q(n4+n5+n8,3)
      - 4.*Q(n1+n2+n6,3)*Q(n3+n7,2)*Q(n4+n5+n8,3)+4.*Q(n2,1)*Q(n6,1)*Q(n1+n3+n7,3)*Q(n4+n5+n8,3)
      - 4.*Q(n2+n6,2)*Q(n1+n3+n7,3)*Q(n4+n5+n8,3)+4.*Q(n1,1)*Q(n6,1)*Q(n2+n3+n7,3)*Q(n4+n5+n8,3)
      - 4.*Q(n1+n6,2)*Q(n2+n3+n7,3)*Q(n4+n5+n8,3)-12.*Q(n6,1)*Q(n1+n2+n3+n7,4)*Q(n4+n5+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n6+n7,2)*Q(n4+n5+n8,3)+2.*Q(n1+n2,2)*Q(n3,1)*Q(n6+n7,2)*Q(n4+n5+n8,3)
      + 2.*Q(n2,1)*Q(n1+n3,2)*Q(n6+n7,2)*Q(n4+n5+n8,3)+2.*Q(n1,1)*Q(n2+n3,2)*Q(n6+n7,2)*Q(n4+n5+n8,3)
      - 4.*Q(n1+n2+n3,3)*Q(n6+n7,2)*Q(n4+n5+n8,3)+4.*Q(n2,1)*Q(n3,1)*Q(n1+n6+n7,3)*Q(n4+n5+n8,3)
      - 4.*Q(n2+n3,2)*Q(n1+n6+n7,3)*Q(n4+n5+n8,3)+4.*Q(n1,1)*Q(n3,1)*Q(n2+n6+n7,3)*Q(n4+n5+n8,3)
      - 4.*Q(n1+n3,2)*Q(n2+n6+n7,3)*Q(n4+n5+n8,3)-12.*Q(n3,1)*Q(n1+n2+n6+n7,4)*Q(n4+n5+n8,3)
      + 4.*Q(n1,1)*Q(n2,1)*Q(n3+n6+n7,3)*Q(n4+n5+n8,3)-4.*Q(n1+n2,2)*Q(n3+n6+n7,3)*Q(n4+n5+n8,3)
      - 12.*Q(n2,1)*Q(n1+n3+n6+n7,4)*Q(n4+n5+n8,3)-12.*Q(n1,1)*Q(n2+n3+n6+n7,4)*Q(n4+n5+n8,3)
      + 48.*Q(n1+n2+n3+n6+n7,5)*Q(n4+n5+n8,3)-6.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n5+n8,4)
      + 6.*Q(n2+n3,2)*Q(n6,1)*Q(n7,1)*Q(n1+n4+n5+n8,4)+6.*Q(n3,1)*Q(n2+n6,2)*Q(n7,1)*Q(n1+n4+n5+n8,4)
      + 6.*Q(n2,1)*Q(n3+n6,2)*Q(n7,1)*Q(n1+n4+n5+n8,4)-12.*Q(n2+n3+n6,3)*Q(n7,1)*Q(n1+n4+n5+n8,4)
      + 6.*Q(n3,1)*Q(n6,1)*Q(n2+n7,2)*Q(n1+n4+n5+n8,4)-6.*Q(n3+n6,2)*Q(n2+n7,2)*Q(n1+n4+n5+n8,4)
      + 6.*Q(n2,1)*Q(n6,1)*Q(n3+n7,2)*Q(n1+n4+n5+n8,4)-6.*Q(n2+n6,2)*Q(n3+n7,2)*Q(n1+n4+n5+n8,4)
      - 12.*Q(n6,1)*Q(n2+n3+n7,3)*Q(n1+n4+n5+n8,4)+6.*Q(n2,1)*Q(n3,1)*Q(n6+n7,2)*Q(n1+n4+n5+n8,4)
      - 6.*Q(n2+n3,2)*Q(n6+n7,2)*Q(n1+n4+n5+n8,4)-12.*Q(n3,1)*Q(n2+n6+n7,3)*Q(n1+n4+n5+n8,4)
      - 12.*Q(n2,1)*Q(n3+n6+n7,3)*Q(n1+n4+n5+n8,4)+36.*Q(n2+n3+n6+n7,4)*Q(n1+n4+n5+n8,4)
      - 6.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n5+n8,4)+6.*Q(n1+n3,2)*Q(n6,1)*Q(n7,1)*Q(n2+n4+n5+n8,4)
      + 6.*Q(n3,1)*Q(n1+n6,2)*Q(n7,1)*Q(n2+n4+n5+n8,4)+6.*Q(n1,1)*Q(n3+n6,2)*Q(n7,1)*Q(n2+n4+n5+n8,4)
      - 12.*Q(n1+n3+n6,3)*Q(n7,1)*Q(n2+n4+n5+n8,4)+6.*Q(n3,1)*Q(n6,1)*Q(n1+n7,2)*Q(n2+n4+n5+n8,4)
      - 6.*Q(n3+n6,2)*Q(n1+n7,2)*Q(n2+n4+n5+n8,4)+6.*Q(n1,1)*Q(n6,1)*Q(n3+n7,2)*Q(n2+n4+n5+n8,4)
      - 6.*Q(n1+n6,2)*Q(n3+n7,2)*Q(n2+n4+n5+n8,4)-12.*Q(n6,1)*Q(n1+n3+n7,3)*Q(n2+n4+n5+n8,4)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n6+n7,2)*Q(n2+n4+n5+n8,4)-6.*Q(n1+n3,2)*Q(n6+n7,2)*Q(n2+n4+n5+n8,4)
      - 12.*Q(n3,1)*Q(n1+n6+n7,3)*Q(n2+n4+n5+n8,4)-12.*Q(n1,1)*Q(n3+n6+n7,3)*Q(n2+n4+n5+n8,4)
      + 36.*Q(n1+n3+n6+n7,4)*Q(n2+n4+n5+n8,4)+24.*Q(n3,1)*Q(n6,1)*Q(n7,1)*Q(n1+n2+n4+n5+n8,5)
      - 24.*Q(n3+n6,2)*Q(n7,1)*Q(n1+n2+n4+n5+n8,5)-24.*Q(n6,1)*Q(n3+n7,2)*Q(n1+n2+n4+n5+n8,5)
      - 24.*Q(n3,1)*Q(n6+n7,2)*Q(n1+n2+n4+n5+n8,5)+48.*Q(n3+n6+n7,3)*Q(n1+n2+n4+n5+n8,5)
      - 6.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n5+n8,4)+6.*Q(n1+n2,2)*Q(n6,1)*Q(n7,1)*Q(n3+n4+n5+n8,4)
      + 6.*Q(n2,1)*Q(n1+n6,2)*Q(n7,1)*Q(n3+n4+n5+n8,4)+6.*Q(n1,1)*Q(n2+n6,2)*Q(n7,1)*Q(n3+n4+n5+n8,4)
      - 12.*Q(n1+n2+n6,3)*Q(n7,1)*Q(n3+n4+n5+n8,4)+6.*Q(n2,1)*Q(n6,1)*Q(n1+n7,2)*Q(n3+n4+n5+n8,4)
      - 6.*Q(n2+n6,2)*Q(n1+n7,2)*Q(n3+n4+n5+n8,4)+6.*Q(n1,1)*Q(n6,1)*Q(n2+n7,2)*Q(n3+n4+n5+n8,4)
      - 6.*Q(n1+n6,2)*Q(n2+n7,2)*Q(n3+n4+n5+n8,4)-12.*Q(n6,1)*Q(n1+n2+n7,3)*Q(n3+n4+n5+n8,4)
      + 6.*Q(n1,1)*Q(n2,1)*Q(n6+n7,2)*Q(n3+n4+n5+n8,4)-6.*Q(n1+n2,2)*Q(n6+n7,2)*Q(n3+n4+n5+n8,4)
      - 12.*Q(n2,1)*Q(n1+n6+n7,3)*Q(n3+n4+n5+n8,4)-12.*Q(n1,1)*Q(n2+n6+n7,3)*Q(n3+n4+n5+n8,4)
      + 36.*Q(n1+n2+n6+n7,4)*Q(n3+n4+n5+n8,4)+24.*Q(n2,1)*Q(n6,1)*Q(n7,1)*Q(n1+n3+n4+n5+n8,5)
      - 24.*Q(n2+n6,2)*Q(n7,1)*Q(n1+n3+n4+n5+n8,5)-24.*Q(n6,1)*Q(n2+n7,2)*Q(n1+n3+n4+n5+n8,5)
      - 24.*Q(n2,1)*Q(n6+n7,2)*Q(n1+n3+n4+n5+n8,5)+48.*Q(n2+n6+n7,3)*Q(n1+n3+n4+n5+n8,5)
      + 24.*Q(n1,1)*Q(n6,1)*Q(n7,1)*Q(n2+n3+n4+n5+n8,5)-24.*Q(n1+n6,2)*Q(n7,1)*Q(n2+n3+n4+n5+n8,5)
      - 24.*Q(n6,1)*Q(n1+n7,2)*Q(n2+n3+n4+n5+n8,5)-24.*Q(n1,1)*Q(n6+n7,2)*Q(n2+n3+n4+n5+n8,5)
      + 48.*Q(n1+n6+n7,3)*Q(n2+n3+n4+n5+n8,5)-120.*Q(n6,1)*Q(n7,1)*Q(n1+n2+n3+n4+n5+n8,6)
      + 120.*Q(n6+n7,2)*Q(n1+n2+n3+n4+n5+n8,6)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n7,1)*Q(n6+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n7,1)*Q(n6+n8,2)
      - Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n7,1)*Q(n6+n8,2)+Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n7,1)*Q(n6+n8,2)
      - Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n7,1)*Q(n6+n8,2)
      - Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n7,1)*Q(n6+n8,2)+2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n7,1)*Q(n6+n8,2)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n7,1)*Q(n6+n8,2)+2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n7,1)*Q(n6+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n7,1)*Q(n6+n8,2)
      - Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n7,1)*Q(n6+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n7,1)*Q(n6+n8,2)-2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n7,1)*Q(n6+n8,2)+6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n7,1)*Q(n6+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n7,1)*Q(n6+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n7,1)*Q(n6+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n7,1)*Q(n6+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n7,1)*Q(n6+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n7,1)*Q(n6+n8,2)
      + 2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n7,1)*Q(n6+n8,2)+6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n7,1)*Q(n6+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n7,1)*Q(n6+n8,2)+2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n7,1)*Q(n6+n8,2)
      + 6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n7,1)*Q(n6+n8,2)+6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n7,1)*Q(n6+n8,2)
      - 24.*Q(n1+n2+n3+n4+n5,5)*Q(n7,1)*Q(n6+n8,2)+Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n7,2)*Q(n6+n8,2)
      - Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n7,2)*Q(n6+n8,2)-Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n7,2)*Q(n6+n8,2)
      - Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n7,2)*Q(n6+n8,2)+2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n7,2)*Q(n6+n8,2)
      - Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n7,2)*Q(n6+n8,2)+Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n7,2)*Q(n6+n8,2)
      - Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n7,2)*Q(n6+n8,2)+Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n7,2)*Q(n6+n8,2)
      + 2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n7,2)*Q(n6+n8,2)-Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n7,2)*Q(n6+n8,2)
      + Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n7,2)*Q(n6+n8,2)+2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n7,2)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n7,2)*Q(n6+n8,2)-6.*Q(n2+n3+n4+n5,4)*Q(n1+n7,2)*Q(n6+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n7,2)*Q(n6+n8,2)-Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n7,2)*Q(n6+n8,2)
      - Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n7,2)*Q(n6+n8,2)
      + 2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n7,2)*Q(n6+n8,2)-Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n7,2)*Q(n6+n8,2)
      + Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n7,2)*Q(n6+n8,2)
      + Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n7,2)*Q(n6+n8,2)+2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n7,2)*Q(n6+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n7,2)*Q(n6+n8,2)+Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n7,2)*Q(n6+n8,2)
      + 2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n7,2)*Q(n6+n8,2)+2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n7,2)*Q(n6+n8,2)
      - 6.*Q(n1+n3+n4+n5,4)*Q(n2+n7,2)*Q(n6+n8,2)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n7,3)*Q(n6+n8,2)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n7,3)*Q(n6+n8,2)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n7,3)*Q(n6+n8,2)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n7,3)*Q(n6+n8,2)-4.*Q(n3+n4+n5,3)*Q(n1+n2+n7,3)*Q(n6+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n7,2)*Q(n6+n8,2)-Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n7,2)*Q(n6+n8,2)
      - Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n7,2)*Q(n6+n8,2)
      + 2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n7,2)*Q(n6+n8,2)-Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n7,2)*Q(n6+n8,2)
      + Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n7,2)*Q(n6+n8,2)
      + Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n7,2)*Q(n6+n8,2)+2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n7,2)*Q(n6+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n7,2)*Q(n6+n8,2)+Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n7,2)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n7,2)*Q(n6+n8,2)+2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n7,2)*Q(n6+n8,2)
      - 6.*Q(n1+n2+n4+n5,4)*Q(n3+n7,2)*Q(n6+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n7,3)*Q(n6+n8,2)
      + 2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n7,3)*Q(n6+n8,2)+2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n7,3)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n7,3)*Q(n6+n8,2)-4.*Q(n2+n4+n5,3)*Q(n1+n3+n7,3)*Q(n6+n8,2)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n7,3)*Q(n6+n8,2)+2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n7,3)*Q(n6+n8,2)
      + 2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n7,3)*Q(n6+n8,2)+2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n7,3)*Q(n6+n8,2)
      - 4.*Q(n1+n4+n5,3)*Q(n2+n3+n7,3)*Q(n6+n8,2)+6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n7,4)*Q(n6+n8,2)
      - 6.*Q(n4+n5,2)*Q(n1+n2+n3+n7,4)*Q(n6+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n7,2)*Q(n6+n8,2)
      - Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n7,2)*Q(n6+n8,2)-Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n7,2)*Q(n6+n8,2)
      - Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n7,2)*Q(n6+n8,2)+2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n7,2)*Q(n6+n8,2)
      - Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n7,2)*Q(n6+n8,2)+Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n7,2)*Q(n6+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n7,2)*Q(n6+n8,2)+Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n7,2)*Q(n6+n8,2)
      + 2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n7,2)*Q(n6+n8,2)
      + Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n7,2)*Q(n6+n8,2)+2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n7,2)*Q(n6+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n7,2)*Q(n6+n8,2)-6.*Q(n1+n2+n3+n5,4)*Q(n4+n7,2)*Q(n6+n8,2)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n7,3)*Q(n6+n8,2)+2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n7,3)*Q(n6+n8,2)
      + 2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n7,3)*Q(n6+n8,2)+2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n7,3)*Q(n6+n8,2)
      - 4.*Q(n2+n3+n5,3)*Q(n1+n4+n7,3)*Q(n6+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n7,3)*Q(n6+n8,2)
      + 2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n7,3)*Q(n6+n8,2)+2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n7,3)*Q(n6+n8,2)
      + 2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n7,3)*Q(n6+n8,2)-4.*Q(n1+n3+n5,3)*Q(n2+n4+n7,3)*Q(n6+n8,2)
      + 6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n7,4)*Q(n6+n8,2)-6.*Q(n3+n5,2)*Q(n1+n2+n4+n7,4)*Q(n6+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n7,3)*Q(n6+n8,2)+2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n7,3)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n7,3)*Q(n6+n8,2)+2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n7,3)*Q(n6+n8,2)
      - 4.*Q(n1+n2+n5,3)*Q(n3+n4+n7,3)*Q(n6+n8,2)+6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n7,4)*Q(n6+n8,2)
      - 6.*Q(n2+n5,2)*Q(n1+n3+n4+n7,4)*Q(n6+n8,2)+6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n7,4)*Q(n6+n8,2)
      - 6.*Q(n1+n5,2)*Q(n2+n3+n4+n7,4)*Q(n6+n8,2)-24.*Q(n5,1)*Q(n1+n2+n3+n4+n7,5)*Q(n6+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n7,2)*Q(n6+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n7,2)*Q(n6+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n7,2)*Q(n6+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n7,2)*Q(n6+n8,2)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n7,2)*Q(n6+n8,2)
      + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n7,2)*Q(n6+n8,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n7,2)*Q(n6+n8,2)
      + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n7,2)*Q(n6+n8,2)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n7,2)*Q(n6+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n7,2)*Q(n6+n8,2)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n7,2)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n7,2)*Q(n6+n8,2)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n7,2)*Q(n6+n8,2)
      - 6.*Q(n1+n2+n3+n4,4)*Q(n5+n7,2)*Q(n6+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n7,3)*Q(n6+n8,2)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n7,3)*Q(n6+n8,2)-4.*Q(n2+n3+n4,3)*Q(n1+n5+n7,3)*Q(n6+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n7,3)*Q(n6+n8,2)+2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n7,3)*Q(n6+n8,2)+2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n7,3)*Q(n6+n8,2)
      - 4.*Q(n1+n3+n4,3)*Q(n2+n5+n7,3)*Q(n6+n8,2)+6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n7,4)*Q(n6+n8,2)
      - 6.*Q(n3+n4,2)*Q(n1+n2+n5+n7,4)*Q(n6+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n7,3)*Q(n6+n8,2)+2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n7,3)*Q(n6+n8,2)-4.*Q(n1+n2+n4,3)*Q(n3+n5+n7,3)*Q(n6+n8,2)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n7,4)*Q(n6+n8,2)-6.*Q(n2+n4,2)*Q(n1+n3+n5+n7,4)*Q(n6+n8,2)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n7,4)*Q(n6+n8,2)-6.*Q(n1+n4,2)*Q(n2+n3+n5+n7,4)*Q(n6+n8,2)
      - 24.*Q(n4,1)*Q(n1+n2+n3+n5+n7,5)*Q(n6+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n7,3)*Q(n6+n8,2)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n7,3)*Q(n6+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n7,3)*Q(n6+n8,2)-4.*Q(n1+n2+n3,3)*Q(n4+n5+n7,3)*Q(n6+n8,2)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n7,4)*Q(n6+n8,2)-6.*Q(n2+n3,2)*Q(n1+n4+n5+n7,4)*Q(n6+n8,2)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n7,4)*Q(n6+n8,2)-6.*Q(n1+n3,2)*Q(n2+n4+n5+n7,4)*Q(n6+n8,2)
      - 24.*Q(n3,1)*Q(n1+n2+n4+n5+n7,5)*Q(n6+n8,2)+6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n7,4)*Q(n6+n8,2)
      - 6.*Q(n1+n2,2)*Q(n3+n4+n5+n7,4)*Q(n6+n8,2)-24.*Q(n2,1)*Q(n1+n3+n4+n5+n7,5)*Q(n6+n8,2)
      - 24.*Q(n1,1)*Q(n2+n3+n4+n5+n7,5)*Q(n6+n8,2)+120.*Q(n1+n2+n3+n4+n5+n7,6)*Q(n6+n8,2)
      + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n1+n6+n8,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n1+n6+n8,3)
      - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n7,1)*Q(n1+n6+n8,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n7,1)*Q(n1+n6+n8,3)
      + 4.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n7,1)*Q(n1+n6+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n7,1)*Q(n1+n6+n8,3)
      + 2.*Q(n3+n4,2)*Q(n2+n5,2)*Q(n7,1)*Q(n1+n6+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n7,1)*Q(n1+n6+n8,3)
      + 2.*Q(n2+n4,2)*Q(n3+n5,2)*Q(n7,1)*Q(n1+n6+n8,3)+4.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n7,1)*Q(n1+n6+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n7,1)*Q(n1+n6+n8,3)+2.*Q(n2+n3,2)*Q(n4+n5,2)*Q(n7,1)*Q(n1+n6+n8,3)
      + 4.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n7,1)*Q(n1+n6+n8,3)+4.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n7,1)*Q(n1+n6+n8,3)
      - 12.*Q(n2+n3+n4+n5,4)*Q(n7,1)*Q(n1+n6+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n2+n7,2)*Q(n1+n6+n8,3)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n2+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n2+n7,2)*Q(n1+n6+n8,3)-4.*Q(n3+n4+n5,3)*Q(n2+n7,2)*Q(n1+n6+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n7,2)*Q(n1+n6+n8,3)+2.*Q(n2+n4,2)*Q(n5,1)*Q(n3+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n4,1)*Q(n2+n5,2)*Q(n3+n7,2)*Q(n1+n6+n8,3)+2.*Q(n2,1)*Q(n4+n5,2)*Q(n3+n7,2)*Q(n1+n6+n8,3)
      - 4.*Q(n2+n4+n5,3)*Q(n3+n7,2)*Q(n1+n6+n8,3)+4.*Q(n4,1)*Q(n5,1)*Q(n2+n3+n7,3)*Q(n1+n6+n8,3)
      - 4.*Q(n4+n5,2)*Q(n2+n3+n7,3)*Q(n1+n6+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n4+n7,2)*Q(n1+n6+n8,3)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n4+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n4+n7,2)*Q(n1+n6+n8,3)-4.*Q(n2+n3+n5,3)*Q(n4+n7,2)*Q(n1+n6+n8,3)
      + 4.*Q(n3,1)*Q(n5,1)*Q(n2+n4+n7,3)*Q(n1+n6+n8,3)-4.*Q(n3+n5,2)*Q(n2+n4+n7,3)*Q(n1+n6+n8,3)
      + 4.*Q(n2,1)*Q(n5,1)*Q(n3+n4+n7,3)*Q(n1+n6+n8,3)-4.*Q(n2+n5,2)*Q(n3+n4+n7,3)*Q(n1+n6+n8,3)
      - 12.*Q(n5,1)*Q(n2+n3+n4+n7,4)*Q(n1+n6+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n5+n7,2)*Q(n1+n6+n8,3)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n5+n7,2)*Q(n1+n6+n8,3)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n5+n7,2)*Q(n1+n6+n8,3)-4.*Q(n2+n3+n4,3)*Q(n5+n7,2)*Q(n1+n6+n8,3)
      + 4.*Q(n3,1)*Q(n4,1)*Q(n2+n5+n7,3)*Q(n1+n6+n8,3)-4.*Q(n3+n4,2)*Q(n2+n5+n7,3)*Q(n1+n6+n8,3)
      + 4.*Q(n2,1)*Q(n4,1)*Q(n3+n5+n7,3)*Q(n1+n6+n8,3)-4.*Q(n2+n4,2)*Q(n3+n5+n7,3)*Q(n1+n6+n8,3)
      - 12.*Q(n4,1)*Q(n2+n3+n5+n7,4)*Q(n1+n6+n8,3)+4.*Q(n2,1)*Q(n3,1)*Q(n4+n5+n7,3)*Q(n1+n6+n8,3)
      - 4.*Q(n2+n3,2)*Q(n4+n5+n7,3)*Q(n1+n6+n8,3)-12.*Q(n3,1)*Q(n2+n4+n5+n7,4)*Q(n1+n6+n8,3)
      - 12.*Q(n2,1)*Q(n3+n4+n5+n7,4)*Q(n1+n6+n8,3)+48.*Q(n2+n3+n4+n5+n7,5)*Q(n1+n6+n8,3)
      + 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n2+n6+n8,3)-2.*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n2+n6+n8,3)
      - 2.*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n7,1)*Q(n2+n6+n8,3)-2.*Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n7,1)*Q(n2+n6+n8,3)
      + 4.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n7,1)*Q(n2+n6+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n7,1)*Q(n2+n6+n8,3)
      + 2.*Q(n3+n4,2)*Q(n1+n5,2)*Q(n7,1)*Q(n2+n6+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n7,1)*Q(n2+n6+n8,3)
      + 2.*Q(n1+n4,2)*Q(n3+n5,2)*Q(n7,1)*Q(n2+n6+n8,3)+4.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n7,1)*Q(n2+n6+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n7,1)*Q(n2+n6+n8,3)+2.*Q(n1+n3,2)*Q(n4+n5,2)*Q(n7,1)*Q(n2+n6+n8,3)
      + 4.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n7,1)*Q(n2+n6+n8,3)+4.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n7,1)*Q(n2+n6+n8,3)
      - 12.*Q(n1+n3+n4+n5,4)*Q(n7,1)*Q(n2+n6+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n7,2)*Q(n2+n6+n8,3)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n7,2)*Q(n2+n6+n8,3)-4.*Q(n3+n4+n5,3)*Q(n1+n7,2)*Q(n2+n6+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n3+n7,2)*Q(n2+n6+n8,3)+2.*Q(n1+n4,2)*Q(n5,1)*Q(n3+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n4,1)*Q(n1+n5,2)*Q(n3+n7,2)*Q(n2+n6+n8,3)+2.*Q(n1,1)*Q(n4+n5,2)*Q(n3+n7,2)*Q(n2+n6+n8,3)
      - 4.*Q(n1+n4+n5,3)*Q(n3+n7,2)*Q(n2+n6+n8,3)+4.*Q(n4,1)*Q(n5,1)*Q(n1+n3+n7,3)*Q(n2+n6+n8,3)
      - 4.*Q(n4+n5,2)*Q(n1+n3+n7,3)*Q(n2+n6+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n4+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n1+n3,2)*Q(n5,1)*Q(n4+n7,2)*Q(n2+n6+n8,3)+2.*Q(n3,1)*Q(n1+n5,2)*Q(n4+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n1,1)*Q(n3+n5,2)*Q(n4+n7,2)*Q(n2+n6+n8,3)-4.*Q(n1+n3+n5,3)*Q(n4+n7,2)*Q(n2+n6+n8,3)
      + 4.*Q(n3,1)*Q(n5,1)*Q(n1+n4+n7,3)*Q(n2+n6+n8,3)-4.*Q(n3+n5,2)*Q(n1+n4+n7,3)*Q(n2+n6+n8,3)
      + 4.*Q(n1,1)*Q(n5,1)*Q(n3+n4+n7,3)*Q(n2+n6+n8,3)-4.*Q(n1+n5,2)*Q(n3+n4+n7,3)*Q(n2+n6+n8,3)
      - 12.*Q(n5,1)*Q(n1+n3+n4+n7,4)*Q(n2+n6+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n1+n3,2)*Q(n4,1)*Q(n5+n7,2)*Q(n2+n6+n8,3)+2.*Q(n3,1)*Q(n1+n4,2)*Q(n5+n7,2)*Q(n2+n6+n8,3)
      + 2.*Q(n1,1)*Q(n3+n4,2)*Q(n5+n7,2)*Q(n2+n6+n8,3)-4.*Q(n1+n3+n4,3)*Q(n5+n7,2)*Q(n2+n6+n8,3)
      + 4.*Q(n3,1)*Q(n4,1)*Q(n1+n5+n7,3)*Q(n2+n6+n8,3)-4.*Q(n3+n4,2)*Q(n1+n5+n7,3)*Q(n2+n6+n8,3)
      + 4.*Q(n1,1)*Q(n4,1)*Q(n3+n5+n7,3)*Q(n2+n6+n8,3)-4.*Q(n1+n4,2)*Q(n3+n5+n7,3)*Q(n2+n6+n8,3)
      - 12.*Q(n4,1)*Q(n1+n3+n5+n7,4)*Q(n2+n6+n8,3)+4.*Q(n1,1)*Q(n3,1)*Q(n4+n5+n7,3)*Q(n2+n6+n8,3)
      - 4.*Q(n1+n3,2)*Q(n4+n5+n7,3)*Q(n2+n6+n8,3)-12.*Q(n3,1)*Q(n1+n4+n5+n7,4)*Q(n2+n6+n8,3)
      - 12.*Q(n1,1)*Q(n3+n4+n5+n7,4)*Q(n2+n6+n8,3)+48.*Q(n1+n3+n4+n5+n7,5)*Q(n2+n6+n8,3)
      - 6.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n1+n2+n6+n8,4)+6.*Q(n3+n4,2)*Q(n5,1)*Q(n7,1)*Q(n1+n2+n6+n8,4)
      + 6.*Q(n4,1)*Q(n3+n5,2)*Q(n7,1)*Q(n1+n2+n6+n8,4)+6.*Q(n3,1)*Q(n4+n5,2)*Q(n7,1)*Q(n1+n2+n6+n8,4)
      - 12.*Q(n3+n4+n5,3)*Q(n7,1)*Q(n1+n2+n6+n8,4)+6.*Q(n4,1)*Q(n5,1)*Q(n3+n7,2)*Q(n1+n2+n6+n8,4)
      - 6.*Q(n4+n5,2)*Q(n3+n7,2)*Q(n1+n2+n6+n8,4)+6.*Q(n3,1)*Q(n5,1)*Q(n4+n7,2)*Q(n1+n2+n6+n8,4)
      - 6.*Q(n3+n5,2)*Q(n4+n7,2)*Q(n1+n2+n6+n8,4)-12.*Q(n5,1)*Q(n3+n4+n7,3)*Q(n1+n2+n6+n8,4)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n5+n7,2)*Q(n1+n2+n6+n8,4)-6.*Q(n3+n4,2)*Q(n5+n7,2)*Q(n1+n2+n6+n8,4)
      - 12.*Q(n4,1)*Q(n3+n5+n7,3)*Q(n1+n2+n6+n8,4)-12.*Q(n3,1)*Q(n4+n5+n7,3)*Q(n1+n2+n6+n8,4)
      + 36.*Q(n3+n4+n5+n7,4)*Q(n1+n2+n6+n8,4)+2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n3+n6+n8,3)
      - 2.*Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n3+n6+n8,3)-2.*Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n7,1)*Q(n3+n6+n8,3)
      - 2.*Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n7,1)*Q(n3+n6+n8,3)+4.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n7,1)*Q(n3+n6+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n7,1)*Q(n3+n6+n8,3)+2.*Q(n2+n4,2)*Q(n1+n5,2)*Q(n7,1)*Q(n3+n6+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n7,1)*Q(n3+n6+n8,3)+2.*Q(n1+n4,2)*Q(n2+n5,2)*Q(n7,1)*Q(n3+n6+n8,3)
      + 4.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n7,1)*Q(n3+n6+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n7,1)*Q(n3+n6+n8,3)
      + 2.*Q(n1+n2,2)*Q(n4+n5,2)*Q(n7,1)*Q(n3+n6+n8,3)+4.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n7,1)*Q(n3+n6+n8,3)
      + 4.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n7,1)*Q(n3+n6+n8,3)-12.*Q(n1+n2+n4+n5,4)*Q(n7,1)*Q(n3+n6+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n7,2)*Q(n3+n6+n8,3)+2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n7,2)*Q(n3+n6+n8,3)
      + 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n7,2)*Q(n3+n6+n8,3)+2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n7,2)*Q(n3+n6+n8,3)
      - 4.*Q(n2+n4+n5,3)*Q(n1+n7,2)*Q(n3+n6+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n7,2)*Q(n3+n6+n8,3)
      + 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n7,2)*Q(n3+n6+n8,3)+2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n7,2)*Q(n3+n6+n8,3)
      + 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n7,2)*Q(n3+n6+n8,3)-4.*Q(n1+n4+n5,3)*Q(n2+n7,2)*Q(n3+n6+n8,3)
      + 4.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n7,3)*Q(n3+n6+n8,3)-4.*Q(n4+n5,2)*Q(n1+n2+n7,3)*Q(n3+n6+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n4+n7,2)*Q(n3+n6+n8,3)+2.*Q(n1+n2,2)*Q(n5,1)*Q(n4+n7,2)*Q(n3+n6+n8,3)
      + 2.*Q(n2,1)*Q(n1+n5,2)*Q(n4+n7,2)*Q(n3+n6+n8,3)+2.*Q(n1,1)*Q(n2+n5,2)*Q(n4+n7,2)*Q(n3+n6+n8,3)
      - 4.*Q(n1+n2+n5,3)*Q(n4+n7,2)*Q(n3+n6+n8,3)+4.*Q(n2,1)*Q(n5,1)*Q(n1+n4+n7,3)*Q(n3+n6+n8,3)
      - 4.*Q(n2+n5,2)*Q(n1+n4+n7,3)*Q(n3+n6+n8,3)+4.*Q(n1,1)*Q(n5,1)*Q(n2+n4+n7,3)*Q(n3+n6+n8,3)
      - 4.*Q(n1+n5,2)*Q(n2+n4+n7,3)*Q(n3+n6+n8,3)-12.*Q(n5,1)*Q(n1+n2+n4+n7,4)*Q(n3+n6+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5+n7,2)*Q(n3+n6+n8,3)+2.*Q(n1+n2,2)*Q(n4,1)*Q(n5+n7,2)*Q(n3+n6+n8,3)
      + 2.*Q(n2,1)*Q(n1+n4,2)*Q(n5+n7,2)*Q(n3+n6+n8,3)+2.*Q(n1,1)*Q(n2+n4,2)*Q(n5+n7,2)*Q(n3+n6+n8,3)
      - 4.*Q(n1+n2+n4,3)*Q(n5+n7,2)*Q(n3+n6+n8,3)+4.*Q(n2,1)*Q(n4,1)*Q(n1+n5+n7,3)*Q(n3+n6+n8,3)
      - 4.*Q(n2+n4,2)*Q(n1+n5+n7,3)*Q(n3+n6+n8,3)+4.*Q(n1,1)*Q(n4,1)*Q(n2+n5+n7,3)*Q(n3+n6+n8,3)
      - 4.*Q(n1+n4,2)*Q(n2+n5+n7,3)*Q(n3+n6+n8,3)-12.*Q(n4,1)*Q(n1+n2+n5+n7,4)*Q(n3+n6+n8,3)
      + 4.*Q(n1,1)*Q(n2,1)*Q(n4+n5+n7,3)*Q(n3+n6+n8,3)-4.*Q(n1+n2,2)*Q(n4+n5+n7,3)*Q(n3+n6+n8,3)
      - 12.*Q(n2,1)*Q(n1+n4+n5+n7,4)*Q(n3+n6+n8,3)-12.*Q(n1,1)*Q(n2+n4+n5+n7,4)*Q(n3+n6+n8,3)
      + 48.*Q(n1+n2+n4+n5+n7,5)*Q(n3+n6+n8,3)-6.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n1+n3+n6+n8,4)
      + 6.*Q(n2+n4,2)*Q(n5,1)*Q(n7,1)*Q(n1+n3+n6+n8,4)+6.*Q(n4,1)*Q(n2+n5,2)*Q(n7,1)*Q(n1+n3+n6+n8,4)
      + 6.*Q(n2,1)*Q(n4+n5,2)*Q(n7,1)*Q(n1+n3+n6+n8,4)-12.*Q(n2+n4+n5,3)*Q(n7,1)*Q(n1+n3+n6+n8,4)
      + 6.*Q(n4,1)*Q(n5,1)*Q(n2+n7,2)*Q(n1+n3+n6+n8,4)-6.*Q(n4+n5,2)*Q(n2+n7,2)*Q(n1+n3+n6+n8,4)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n4+n7,2)*Q(n1+n3+n6+n8,4)-6.*Q(n2+n5,2)*Q(n4+n7,2)*Q(n1+n3+n6+n8,4)
      - 12.*Q(n5,1)*Q(n2+n4+n7,3)*Q(n1+n3+n6+n8,4)+6.*Q(n2,1)*Q(n4,1)*Q(n5+n7,2)*Q(n1+n3+n6+n8,4)
      - 6.*Q(n2+n4,2)*Q(n5+n7,2)*Q(n1+n3+n6+n8,4)-12.*Q(n4,1)*Q(n2+n5+n7,3)*Q(n1+n3+n6+n8,4)
      - 12.*Q(n2,1)*Q(n4+n5+n7,3)*Q(n1+n3+n6+n8,4)+36.*Q(n2+n4+n5+n7,4)*Q(n1+n3+n6+n8,4)
      - 6.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n2+n3+n6+n8,4)+6.*Q(n1+n4,2)*Q(n5,1)*Q(n7,1)*Q(n2+n3+n6+n8,4)
      + 6.*Q(n4,1)*Q(n1+n5,2)*Q(n7,1)*Q(n2+n3+n6+n8,4)+6.*Q(n1,1)*Q(n4+n5,2)*Q(n7,1)*Q(n2+n3+n6+n8,4)
      - 12.*Q(n1+n4+n5,3)*Q(n7,1)*Q(n2+n3+n6+n8,4)+6.*Q(n4,1)*Q(n5,1)*Q(n1+n7,2)*Q(n2+n3+n6+n8,4)
      - 6.*Q(n4+n5,2)*Q(n1+n7,2)*Q(n2+n3+n6+n8,4)+6.*Q(n1,1)*Q(n5,1)*Q(n4+n7,2)*Q(n2+n3+n6+n8,4)
      - 6.*Q(n1+n5,2)*Q(n4+n7,2)*Q(n2+n3+n6+n8,4)-12.*Q(n5,1)*Q(n1+n4+n7,3)*Q(n2+n3+n6+n8,4)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n5+n7,2)*Q(n2+n3+n6+n8,4)-6.*Q(n1+n4,2)*Q(n5+n7,2)*Q(n2+n3+n6+n8,4)
      - 12.*Q(n4,1)*Q(n1+n5+n7,3)*Q(n2+n3+n6+n8,4)-12.*Q(n1,1)*Q(n4+n5+n7,3)*Q(n2+n3+n6+n8,4)
      + 36.*Q(n1+n4+n5+n7,4)*Q(n2+n3+n6+n8,4)+24.*Q(n4,1)*Q(n5,1)*Q(n7,1)*Q(n1+n2+n3+n6+n8,5)
      - 24.*Q(n4+n5,2)*Q(n7,1)*Q(n1+n2+n3+n6+n8,5)-24.*Q(n5,1)*Q(n4+n7,2)*Q(n1+n2+n3+n6+n8,5)
      - 24.*Q(n4,1)*Q(n5+n7,2)*Q(n1+n2+n3+n6+n8,5)+48.*Q(n4+n5+n7,3)*Q(n1+n2+n3+n6+n8,5)
      + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n7,1)*Q(n4+n6+n8,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n7,1)*Q(n4+n6+n8,3)
    );
} 
//=======================================================
TComplex AliAnalysisTaskChargedFlow::Eight_4(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
	return
		(
			- 2.*Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n7,1)*Q(n4+n6+n8,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n7,1)*Q(n4+n6+n8,3)
      + 4.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n7,1)*Q(n4+n6+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n7,1)*Q(n4+n6+n8,3)
      + 2.*Q(n2+n3,2)*Q(n1+n5,2)*Q(n7,1)*Q(n4+n6+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n7,1)*Q(n4+n6+n8,3)
      + 2.*Q(n1+n3,2)*Q(n2+n5,2)*Q(n7,1)*Q(n4+n6+n8,3)+4.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n7,1)*Q(n4+n6+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n7,1)*Q(n4+n6+n8,3)+2.*Q(n1+n2,2)*Q(n3+n5,2)*Q(n7,1)*Q(n4+n6+n8,3)
      + 4.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n7,1)*Q(n4+n6+n8,3)+4.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n7,1)*Q(n4+n6+n8,3)
      - 12.*Q(n1+n2+n3+n5,4)*Q(n7,1)*Q(n4+n6+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n7,2)*Q(n4+n6+n8,3)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n7,2)*Q(n4+n6+n8,3)-4.*Q(n2+n3+n5,3)*Q(n1+n7,2)*Q(n4+n6+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n7,2)*Q(n4+n6+n8,3)+2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n7,2)*Q(n4+n6+n8,3)+2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n7,2)*Q(n4+n6+n8,3)
      - 4.*Q(n1+n3+n5,3)*Q(n2+n7,2)*Q(n4+n6+n8,3)+4.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n7,3)*Q(n4+n6+n8,3)
      - 4.*Q(n3+n5,2)*Q(n1+n2+n7,3)*Q(n4+n6+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n7,2)*Q(n4+n6+n8,3)+2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n7,2)*Q(n4+n6+n8,3)-4.*Q(n1+n2+n5,3)*Q(n3+n7,2)*Q(n4+n6+n8,3)
      + 4.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n7,3)*Q(n4+n6+n8,3)-4.*Q(n2+n5,2)*Q(n1+n3+n7,3)*Q(n4+n6+n8,3)
      + 4.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n7,3)*Q(n4+n6+n8,3)-4.*Q(n1+n5,2)*Q(n2+n3+n7,3)*Q(n4+n6+n8,3)
      - 12.*Q(n5,1)*Q(n1+n2+n3+n7,4)*Q(n4+n6+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n5+n7,2)*Q(n4+n6+n8,3)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n5+n7,2)*Q(n4+n6+n8,3)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n5+n7,2)*Q(n4+n6+n8,3)-4.*Q(n1+n2+n3,3)*Q(n5+n7,2)*Q(n4+n6+n8,3)
      + 4.*Q(n2,1)*Q(n3,1)*Q(n1+n5+n7,3)*Q(n4+n6+n8,3)-4.*Q(n2+n3,2)*Q(n1+n5+n7,3)*Q(n4+n6+n8,3)
      + 4.*Q(n1,1)*Q(n3,1)*Q(n2+n5+n7,3)*Q(n4+n6+n8,3)-4.*Q(n1+n3,2)*Q(n2+n5+n7,3)*Q(n4+n6+n8,3)
      - 12.*Q(n3,1)*Q(n1+n2+n5+n7,4)*Q(n4+n6+n8,3)+4.*Q(n1,1)*Q(n2,1)*Q(n3+n5+n7,3)*Q(n4+n6+n8,3)
      - 4.*Q(n1+n2,2)*Q(n3+n5+n7,3)*Q(n4+n6+n8,3)-12.*Q(n2,1)*Q(n1+n3+n5+n7,4)*Q(n4+n6+n8,3)
      - 12.*Q(n1,1)*Q(n2+n3+n5+n7,4)*Q(n4+n6+n8,3)+48.*Q(n1+n2+n3+n5+n7,5)*Q(n4+n6+n8,3)
      - 6.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n7,1)*Q(n1+n4+n6+n8,4)+6.*Q(n2+n3,2)*Q(n5,1)*Q(n7,1)*Q(n1+n4+n6+n8,4)
      + 6.*Q(n3,1)*Q(n2+n5,2)*Q(n7,1)*Q(n1+n4+n6+n8,4)+6.*Q(n2,1)*Q(n3+n5,2)*Q(n7,1)*Q(n1+n4+n6+n8,4)
      - 12.*Q(n2+n3+n5,3)*Q(n7,1)*Q(n1+n4+n6+n8,4)+6.*Q(n3,1)*Q(n5,1)*Q(n2+n7,2)*Q(n1+n4+n6+n8,4)
      - 6.*Q(n3+n5,2)*Q(n2+n7,2)*Q(n1+n4+n6+n8,4)+6.*Q(n2,1)*Q(n5,1)*Q(n3+n7,2)*Q(n1+n4+n6+n8,4)
      - 6.*Q(n2+n5,2)*Q(n3+n7,2)*Q(n1+n4+n6+n8,4)-12.*Q(n5,1)*Q(n2+n3+n7,3)*Q(n1+n4+n6+n8,4)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n5+n7,2)*Q(n1+n4+n6+n8,4)-6.*Q(n2+n3,2)*Q(n5+n7,2)*Q(n1+n4+n6+n8,4)
      - 12.*Q(n3,1)*Q(n2+n5+n7,3)*Q(n1+n4+n6+n8,4)-12.*Q(n2,1)*Q(n3+n5+n7,3)*Q(n1+n4+n6+n8,4)
      + 36.*Q(n2+n3+n5+n7,4)*Q(n1+n4+n6+n8,4)-6.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n7,1)*Q(n2+n4+n6+n8,4)
      + 6.*Q(n1+n3,2)*Q(n5,1)*Q(n7,1)*Q(n2+n4+n6+n8,4)+6.*Q(n3,1)*Q(n1+n5,2)*Q(n7,1)*Q(n2+n4+n6+n8,4)
      + 6.*Q(n1,1)*Q(n3+n5,2)*Q(n7,1)*Q(n2+n4+n6+n8,4)-12.*Q(n1+n3+n5,3)*Q(n7,1)*Q(n2+n4+n6+n8,4)
      + 6.*Q(n3,1)*Q(n5,1)*Q(n1+n7,2)*Q(n2+n4+n6+n8,4)-6.*Q(n3+n5,2)*Q(n1+n7,2)*Q(n2+n4+n6+n8,4)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n3+n7,2)*Q(n2+n4+n6+n8,4)-6.*Q(n1+n5,2)*Q(n3+n7,2)*Q(n2+n4+n6+n8,4)
      - 12.*Q(n5,1)*Q(n1+n3+n7,3)*Q(n2+n4+n6+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n5+n7,2)*Q(n2+n4+n6+n8,4)
      - 6.*Q(n1+n3,2)*Q(n5+n7,2)*Q(n2+n4+n6+n8,4)-12.*Q(n3,1)*Q(n1+n5+n7,3)*Q(n2+n4+n6+n8,4)
      - 12.*Q(n1,1)*Q(n3+n5+n7,3)*Q(n2+n4+n6+n8,4)+36.*Q(n1+n3+n5+n7,4)*Q(n2+n4+n6+n8,4)
      + 24.*Q(n3,1)*Q(n5,1)*Q(n7,1)*Q(n1+n2+n4+n6+n8,5)-24.*Q(n3+n5,2)*Q(n7,1)*Q(n1+n2+n4+n6+n8,5)
      - 24.*Q(n5,1)*Q(n3+n7,2)*Q(n1+n2+n4+n6+n8,5)-24.*Q(n3,1)*Q(n5+n7,2)*Q(n1+n2+n4+n6+n8,5)
      + 48.*Q(n3+n5+n7,3)*Q(n1+n2+n4+n6+n8,5)-6.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n7,1)*Q(n3+n4+n6+n8,4)
      + 6.*Q(n1+n2,2)*Q(n5,1)*Q(n7,1)*Q(n3+n4+n6+n8,4)+6.*Q(n2,1)*Q(n1+n5,2)*Q(n7,1)*Q(n3+n4+n6+n8,4)
      + 6.*Q(n1,1)*Q(n2+n5,2)*Q(n7,1)*Q(n3+n4+n6+n8,4)-12.*Q(n1+n2+n5,3)*Q(n7,1)*Q(n3+n4+n6+n8,4)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n1+n7,2)*Q(n3+n4+n6+n8,4)-6.*Q(n2+n5,2)*Q(n1+n7,2)*Q(n3+n4+n6+n8,4)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n2+n7,2)*Q(n3+n4+n6+n8,4)-6.*Q(n1+n5,2)*Q(n2+n7,2)*Q(n3+n4+n6+n8,4)
      - 12.*Q(n5,1)*Q(n1+n2+n7,3)*Q(n3+n4+n6+n8,4)+6.*Q(n1,1)*Q(n2,1)*Q(n5+n7,2)*Q(n3+n4+n6+n8,4)
      - 6.*Q(n1+n2,2)*Q(n5+n7,2)*Q(n3+n4+n6+n8,4)-12.*Q(n2,1)*Q(n1+n5+n7,3)*Q(n3+n4+n6+n8,4)
      - 12.*Q(n1,1)*Q(n2+n5+n7,3)*Q(n3+n4+n6+n8,4)+36.*Q(n1+n2+n5+n7,4)*Q(n3+n4+n6+n8,4)
      + 24.*Q(n2,1)*Q(n5,1)*Q(n7,1)*Q(n1+n3+n4+n6+n8,5)-24.*Q(n2+n5,2)*Q(n7,1)*Q(n1+n3+n4+n6+n8,5)
      - 24.*Q(n5,1)*Q(n2+n7,2)*Q(n1+n3+n4+n6+n8,5)-24.*Q(n2,1)*Q(n5+n7,2)*Q(n1+n3+n4+n6+n8,5)
      + 48.*Q(n2+n5+n7,3)*Q(n1+n3+n4+n6+n8,5)+24.*Q(n1,1)*Q(n5,1)*Q(n7,1)*Q(n2+n3+n4+n6+n8,5)
      - 24.*Q(n1+n5,2)*Q(n7,1)*Q(n2+n3+n4+n6+n8,5)-24.*Q(n5,1)*Q(n1+n7,2)*Q(n2+n3+n4+n6+n8,5)
      - 24.*Q(n1,1)*Q(n5+n7,2)*Q(n2+n3+n4+n6+n8,5)+48.*Q(n1+n5+n7,3)*Q(n2+n3+n4+n6+n8,5)
      - 120.*Q(n5,1)*Q(n7,1)*Q(n1+n2+n3+n4+n6+n8,6)+120.*Q(n5+n7,2)*Q(n1+n2+n3+n4+n6+n8,6)
      + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n7,1)*Q(n5+n6+n8,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n7,1)*Q(n5+n6+n8,3)
      - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n7,1)*Q(n5+n6+n8,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n7,1)*Q(n5+n6+n8,3)
      + 4.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n7,1)*Q(n5+n6+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n7,1)*Q(n5+n6+n8,3)
      + 2.*Q(n2+n3,2)*Q(n1+n4,2)*Q(n7,1)*Q(n5+n6+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n7,1)*Q(n5+n6+n8,3)
      + 2.*Q(n1+n3,2)*Q(n2+n4,2)*Q(n7,1)*Q(n5+n6+n8,3)+4.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n7,1)*Q(n5+n6+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n7,1)*Q(n5+n6+n8,3)+2.*Q(n1+n2,2)*Q(n3+n4,2)*Q(n7,1)*Q(n5+n6+n8,3)
      + 4.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n7,1)*Q(n5+n6+n8,3)+4.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n7,1)*Q(n5+n6+n8,3)
      - 12.*Q(n1+n2+n3+n4,4)*Q(n7,1)*Q(n5+n6+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n7,2)*Q(n5+n6+n8,3)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n7,2)*Q(n5+n6+n8,3)-4.*Q(n2+n3+n4,3)*Q(n1+n7,2)*Q(n5+n6+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n7,2)*Q(n5+n6+n8,3)+2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n7,2)*Q(n5+n6+n8,3)+2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n7,2)*Q(n5+n6+n8,3)
      - 4.*Q(n1+n3+n4,3)*Q(n2+n7,2)*Q(n5+n6+n8,3)+4.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n7,3)*Q(n5+n6+n8,3)
      - 4.*Q(n3+n4,2)*Q(n1+n2+n7,3)*Q(n5+n6+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n7,2)*Q(n5+n6+n8,3)+2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n7,2)*Q(n5+n6+n8,3)-4.*Q(n1+n2+n4,3)*Q(n3+n7,2)*Q(n5+n6+n8,3)
      + 4.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n7,3)*Q(n5+n6+n8,3)-4.*Q(n2+n4,2)*Q(n1+n3+n7,3)*Q(n5+n6+n8,3)
      + 4.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n7,3)*Q(n5+n6+n8,3)-4.*Q(n1+n4,2)*Q(n2+n3+n7,3)*Q(n5+n6+n8,3)
      - 12.*Q(n4,1)*Q(n1+n2+n3+n7,4)*Q(n5+n6+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n7,2)*Q(n5+n6+n8,3)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n7,2)*Q(n5+n6+n8,3)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n7,2)*Q(n5+n6+n8,3)-4.*Q(n1+n2+n3,3)*Q(n4+n7,2)*Q(n5+n6+n8,3)
      + 4.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n7,3)*Q(n5+n6+n8,3)-4.*Q(n2+n3,2)*Q(n1+n4+n7,3)*Q(n5+n6+n8,3)
      + 4.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n7,3)*Q(n5+n6+n8,3)-4.*Q(n1+n3,2)*Q(n2+n4+n7,3)*Q(n5+n6+n8,3)
      - 12.*Q(n3,1)*Q(n1+n2+n4+n7,4)*Q(n5+n6+n8,3)+4.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n7,3)*Q(n5+n6+n8,3)
      - 4.*Q(n1+n2,2)*Q(n3+n4+n7,3)*Q(n5+n6+n8,3)-12.*Q(n2,1)*Q(n1+n3+n4+n7,4)*Q(n5+n6+n8,3)
      - 12.*Q(n1,1)*Q(n2+n3+n4+n7,4)*Q(n5+n6+n8,3)+48.*Q(n1+n2+n3+n4+n7,5)*Q(n5+n6+n8,3)
      - 6.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n7,1)*Q(n1+n5+n6+n8,4)+6.*Q(n2+n3,2)*Q(n4,1)*Q(n7,1)*Q(n1+n5+n6+n8,4)
      + 6.*Q(n3,1)*Q(n2+n4,2)*Q(n7,1)*Q(n1+n5+n6+n8,4)+6.*Q(n2,1)*Q(n3+n4,2)*Q(n7,1)*Q(n1+n5+n6+n8,4)
      - 12.*Q(n2+n3+n4,3)*Q(n7,1)*Q(n1+n5+n6+n8,4)+6.*Q(n3,1)*Q(n4,1)*Q(n2+n7,2)*Q(n1+n5+n6+n8,4)
      - 6.*Q(n3+n4,2)*Q(n2+n7,2)*Q(n1+n5+n6+n8,4)+6.*Q(n2,1)*Q(n4,1)*Q(n3+n7,2)*Q(n1+n5+n6+n8,4)
      - 6.*Q(n2+n4,2)*Q(n3+n7,2)*Q(n1+n5+n6+n8,4)-12.*Q(n4,1)*Q(n2+n3+n7,3)*Q(n1+n5+n6+n8,4)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n4+n7,2)*Q(n1+n5+n6+n8,4)-6.*Q(n2+n3,2)*Q(n4+n7,2)*Q(n1+n5+n6+n8,4)
      - 12.*Q(n3,1)*Q(n2+n4+n7,3)*Q(n1+n5+n6+n8,4)-12.*Q(n2,1)*Q(n3+n4+n7,3)*Q(n1+n5+n6+n8,4)
      + 36.*Q(n2+n3+n4+n7,4)*Q(n1+n5+n6+n8,4)-6.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n7,1)*Q(n2+n5+n6+n8,4)
      + 6.*Q(n1+n3,2)*Q(n4,1)*Q(n7,1)*Q(n2+n5+n6+n8,4)+6.*Q(n3,1)*Q(n1+n4,2)*Q(n7,1)*Q(n2+n5+n6+n8,4)
      + 6.*Q(n1,1)*Q(n3+n4,2)*Q(n7,1)*Q(n2+n5+n6+n8,4)-12.*Q(n1+n3+n4,3)*Q(n7,1)*Q(n2+n5+n6+n8,4)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n1+n7,2)*Q(n2+n5+n6+n8,4)-6.*Q(n3+n4,2)*Q(n1+n7,2)*Q(n2+n5+n6+n8,4)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n3+n7,2)*Q(n2+n5+n6+n8,4)-6.*Q(n1+n4,2)*Q(n3+n7,2)*Q(n2+n5+n6+n8,4)
      - 12.*Q(n4,1)*Q(n1+n3+n7,3)*Q(n2+n5+n6+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n4+n7,2)*Q(n2+n5+n6+n8,4)
      - 6.*Q(n1+n3,2)*Q(n4+n7,2)*Q(n2+n5+n6+n8,4)-12.*Q(n3,1)*Q(n1+n4+n7,3)*Q(n2+n5+n6+n8,4)
      - 12.*Q(n1,1)*Q(n3+n4+n7,3)*Q(n2+n5+n6+n8,4)+36.*Q(n1+n3+n4+n7,4)*Q(n2+n5+n6+n8,4)
      + 24.*Q(n3,1)*Q(n4,1)*Q(n7,1)*Q(n1+n2+n5+n6+n8,5)-24.*Q(n3+n4,2)*Q(n7,1)*Q(n1+n2+n5+n6+n8,5)
      - 24.*Q(n4,1)*Q(n3+n7,2)*Q(n1+n2+n5+n6+n8,5)-24.*Q(n3,1)*Q(n4+n7,2)*Q(n1+n2+n5+n6+n8,5)
      + 48.*Q(n3+n4+n7,3)*Q(n1+n2+n5+n6+n8,5)-6.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n7,1)*Q(n3+n5+n6+n8,4)
      + 6.*Q(n1+n2,2)*Q(n4,1)*Q(n7,1)*Q(n3+n5+n6+n8,4)+6.*Q(n2,1)*Q(n1+n4,2)*Q(n7,1)*Q(n3+n5+n6+n8,4)
      + 6.*Q(n1,1)*Q(n2+n4,2)*Q(n7,1)*Q(n3+n5+n6+n8,4)-12.*Q(n1+n2+n4,3)*Q(n7,1)*Q(n3+n5+n6+n8,4)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n1+n7,2)*Q(n3+n5+n6+n8,4)-6.*Q(n2+n4,2)*Q(n1+n7,2)*Q(n3+n5+n6+n8,4)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n2+n7,2)*Q(n3+n5+n6+n8,4)-6.*Q(n1+n4,2)*Q(n2+n7,2)*Q(n3+n5+n6+n8,4)
      - 12.*Q(n4,1)*Q(n1+n2+n7,3)*Q(n3+n5+n6+n8,4)+6.*Q(n1,1)*Q(n2,1)*Q(n4+n7,2)*Q(n3+n5+n6+n8,4)
      - 6.*Q(n1+n2,2)*Q(n4+n7,2)*Q(n3+n5+n6+n8,4)-12.*Q(n2,1)*Q(n1+n4+n7,3)*Q(n3+n5+n6+n8,4)
      - 12.*Q(n1,1)*Q(n2+n4+n7,3)*Q(n3+n5+n6+n8,4)+36.*Q(n1+n2+n4+n7,4)*Q(n3+n5+n6+n8,4)
      + 24.*Q(n2,1)*Q(n4,1)*Q(n7,1)*Q(n1+n3+n5+n6+n8,5)-24.*Q(n2+n4,2)*Q(n7,1)*Q(n1+n3+n5+n6+n8,5)
      - 24.*Q(n4,1)*Q(n2+n7,2)*Q(n1+n3+n5+n6+n8,5)-24.*Q(n2,1)*Q(n4+n7,2)*Q(n1+n3+n5+n6+n8,5)
      + 48.*Q(n2+n4+n7,3)*Q(n1+n3+n5+n6+n8,5)+24.*Q(n1,1)*Q(n4,1)*Q(n7,1)*Q(n2+n3+n5+n6+n8,5)
      - 24.*Q(n1+n4,2)*Q(n7,1)*Q(n2+n3+n5+n6+n8,5)-24.*Q(n4,1)*Q(n1+n7,2)*Q(n2+n3+n5+n6+n8,5)
      - 24.*Q(n1,1)*Q(n4+n7,2)*Q(n2+n3+n5+n6+n8,5)+48.*Q(n1+n4+n7,3)*Q(n2+n3+n5+n6+n8,5)
      - 120.*Q(n4,1)*Q(n7,1)*Q(n1+n2+n3+n5+n6+n8,6)+120.*Q(n4+n7,2)*Q(n1+n2+n3+n5+n6+n8,6)
      - 6.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n7,1)*Q(n4+n5+n6+n8,4)+6.*Q(n1+n2,2)*Q(n3,1)*Q(n7,1)*Q(n4+n5+n6+n8,4)
      + 6.*Q(n2,1)*Q(n1+n3,2)*Q(n7,1)*Q(n4+n5+n6+n8,4)+6.*Q(n1,1)*Q(n2+n3,2)*Q(n7,1)*Q(n4+n5+n6+n8,4)
      - 12.*Q(n1+n2+n3,3)*Q(n7,1)*Q(n4+n5+n6+n8,4)+6.*Q(n2,1)*Q(n3,1)*Q(n1+n7,2)*Q(n4+n5+n6+n8,4)
      - 6.*Q(n2+n3,2)*Q(n1+n7,2)*Q(n4+n5+n6+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n2+n7,2)*Q(n4+n5+n6+n8,4)
      - 6.*Q(n1+n3,2)*Q(n2+n7,2)*Q(n4+n5+n6+n8,4)-12.*Q(n3,1)*Q(n1+n2+n7,3)*Q(n4+n5+n6+n8,4)
      + 6.*Q(n1,1)*Q(n2,1)*Q(n3+n7,2)*Q(n4+n5+n6+n8,4)-6.*Q(n1+n2,2)*Q(n3+n7,2)*Q(n4+n5+n6+n8,4)
      - 12.*Q(n2,1)*Q(n1+n3+n7,3)*Q(n4+n5+n6+n8,4)-12.*Q(n1,1)*Q(n2+n3+n7,3)*Q(n4+n5+n6+n8,4)
      + 36.*Q(n1+n2+n3+n7,4)*Q(n4+n5+n6+n8,4)+24.*Q(n2,1)*Q(n3,1)*Q(n7,1)*Q(n1+n4+n5+n6+n8,5)
      - 24.*Q(n2+n3,2)*Q(n7,1)*Q(n1+n4+n5+n6+n8,5)-24.*Q(n3,1)*Q(n2+n7,2)*Q(n1+n4+n5+n6+n8,5)
      - 24.*Q(n2,1)*Q(n3+n7,2)*Q(n1+n4+n5+n6+n8,5)+48.*Q(n2+n3+n7,3)*Q(n1+n4+n5+n6+n8,5)
      + 24.*Q(n1,1)*Q(n3,1)*Q(n7,1)*Q(n2+n4+n5+n6+n8,5)-24.*Q(n1+n3,2)*Q(n7,1)*Q(n2+n4+n5+n6+n8,5)
      - 24.*Q(n3,1)*Q(n1+n7,2)*Q(n2+n4+n5+n6+n8,5)-24.*Q(n1,1)*Q(n3+n7,2)*Q(n2+n4+n5+n6+n8,5)
      + 48.*Q(n1+n3+n7,3)*Q(n2+n4+n5+n6+n8,5)-120.*Q(n3,1)*Q(n7,1)*Q(n1+n2+n4+n5+n6+n8,6)
      + 120.*Q(n3+n7,2)*Q(n1+n2+n4+n5+n6+n8,6)+24.*Q(n1,1)*Q(n2,1)*Q(n7,1)*Q(n3+n4+n5+n6+n8,5)
      - 24.*Q(n1+n2,2)*Q(n7,1)*Q(n3+n4+n5+n6+n8,5)-24.*Q(n2,1)*Q(n1+n7,2)*Q(n3+n4+n5+n6+n8,5)
      - 24.*Q(n1,1)*Q(n2+n7,2)*Q(n3+n4+n5+n6+n8,5)+48.*Q(n1+n2+n7,3)*Q(n3+n4+n5+n6+n8,5)
      - 120.*Q(n2,1)*Q(n7,1)*Q(n1+n3+n4+n5+n6+n8,6)+120.*Q(n2+n7,2)*Q(n1+n3+n4+n5+n6+n8,6)
      - 120.*Q(n1,1)*Q(n7,1)*Q(n2+n3+n4+n5+n6+n8,6)+120.*Q(n1+n7,2)*Q(n2+n3+n4+n5+n6+n8,6)
      + 720.*Q(n7,1)*Q(n1+n2+n3+n4+n5+n6+n8,7)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)*Q(n7+n8,2)
      + Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n7+n8,2)
      - Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)*Q(n7+n8,2)+Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7+n8,2)
      - Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7+n8,2)
      - Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n7+n8,2)+2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)*Q(n7+n8,2)
      - 2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7+n8,2)+2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n7+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n7+n8,2)
      - Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)*Q(n7+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n7+n8,2)-2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n7+n8,2)+6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)*Q(n7+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n7+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)*Q(n7+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n7+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7+n8,2)
      + 2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n7+n8,2)+6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)*Q(n7+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7+n8,2)+2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n7+n8,2)
      + 6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)*Q(n7+n8,2)+6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)*Q(n7+n8,2)
      - 24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)*Q(n7+n8,2)+Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7+n8,2)
      - Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n7+n8,2)-Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7+n8,2)
      - Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n7+n8,2)+2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)*Q(n7+n8,2)
      - Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7+n8,2)+Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n7+n8,2)
      - Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7+n8,2)+Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n7+n8,2)
      + 2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)*Q(n7+n8,2)-Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7+n8,2)
      + Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n7+n8,2)+2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)*Q(n7+n8,2)-6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)*Q(n7+n8,2)
      + Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7+n8,2)-Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n7+n8,2)
      - Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n7+n8,2)
      + 2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)*Q(n7+n8,2)-Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7+n8,2)
      + Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7+n8,2)
      + Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n7+n8,2)+2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)*Q(n7+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7+n8,2)+Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n7+n8,2)
      + 2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)*Q(n7+n8,2)+2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)*Q(n7+n8,2)
      - 6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)*Q(n7+n8,2)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7+n8,2)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n7+n8,2)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)*Q(n7+n8,2)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)*Q(n7+n8,2)-4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)*Q(n7+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7+n8,2)-Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n7+n8,2)
      - Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n7+n8,2)
      + 2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)*Q(n7+n8,2)-Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7+n8,2)
      + Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7+n8,2)
      + Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n7+n8,2)+2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)*Q(n7+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7+n8,2)+Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)*Q(n7+n8,2)+2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)*Q(n7+n8,2)
      - 6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)*Q(n7+n8,2)-2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7+n8,2)
      + 2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n7+n8,2)+2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)*Q(n7+n8,2)-4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)*Q(n7+n8,2)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7+n8,2)+2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n7+n8,2)
      + 2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)*Q(n7+n8,2)+2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)*Q(n7+n8,2)
      - 4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)*Q(n7+n8,2)+6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)*Q(n7+n8,2)
      - 6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)*Q(n7+n8,2)+Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7+n8,2)
      - Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n7+n8,2)-Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7+n8,2)
      - Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n7+n8,2)+2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)*Q(n7+n8,2)
      - Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7+n8,2)+Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n7+n8,2)
      - Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7+n8,2)+Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n7+n8,2)
      + 2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7+n8,2)
      + Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n7+n8,2)+2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)*Q(n7+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)*Q(n7+n8,2)-6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)*Q(n7+n8,2)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7+n8,2)+2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n7+n8,2)
      + 2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)*Q(n7+n8,2)+2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)*Q(n7+n8,2)
      - 4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)*Q(n7+n8,2)-2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7+n8,2)
      + 2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n7+n8,2)+2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)*Q(n7+n8,2)
      + 2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)*Q(n7+n8,2)-4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)*Q(n7+n8,2)
      + 6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)*Q(n7+n8,2)-6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)*Q(n7+n8,2)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7+n8,2)+2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)*Q(n7+n8,2)+2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)*Q(n7+n8,2)
      - 4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)*Q(n7+n8,2)+6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)*Q(n7+n8,2)
      - 6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)*Q(n7+n8,2)+6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)*Q(n7+n8,2)
      - 6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)*Q(n7+n8,2)-24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)*Q(n7+n8,2)
      + Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7+n8,2)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n7+n8,2)
      - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n7+n8,2)
      + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)*Q(n7+n8,2)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7+n8,2)
      + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n7+n8,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7+n8,2)
      + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n7+n8,2)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)*Q(n7+n8,2)
      - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7+n8,2)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)*Q(n7+n8,2)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)*Q(n7+n8,2)
      - 6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)*Q(n7+n8,2)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n7+n8,2)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)*Q(n7+n8,2)-4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)*Q(n7+n8,2)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7+n8,2)+2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)*Q(n7+n8,2)+2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)*Q(n7+n8,2)
      - 4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)*Q(n7+n8,2)+6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)*Q(n7+n8,2)
      - 6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)*Q(n7+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n7+n8,2)+2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)*Q(n7+n8,2)-4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)*Q(n7+n8,2)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)*Q(n7+n8,2)-6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)*Q(n7+n8,2)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)*Q(n7+n8,2)-6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)*Q(n7+n8,2)
      - 24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)*Q(n7+n8,2)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n7+n8,2)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)*Q(n7+n8,2)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)*Q(n7+n8,2)-4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)*Q(n7+n8,2)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)*Q(n7+n8,2)-6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)*Q(n7+n8,2)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)*Q(n7+n8,2)-6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)*Q(n7+n8,2)
      - 24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)*Q(n7+n8,2)+6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)*Q(n7+n8,2)
      - 6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)*Q(n7+n8,2)-24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)*Q(n7+n8,2)
      - 24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)*Q(n7+n8,2)+120.*Q(n1+n2+n3+n4+n5+n6,6)*Q(n7+n8,2)
      + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7+n8,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n7+n8,3)
      - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7+n8,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n7+n8,3)
      + 4.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n1+n7+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7+n8,3)
      + 2.*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n7+n8,3)-2.*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7+n8,3)
      + 2.*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n7+n8,3)+4.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n1+n7+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7+n8,3)+2.*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n7+n8,3)
      + 4.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n1+n7+n8,3)+4.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n1+n7+n8,3)
      - 12.*Q(n2+n3+n4+n5,4)*Q(n6,1)*Q(n1+n7+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n7+n8,3)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n1+n7+n8,3)-4.*Q(n3+n4+n5,3)*Q(n2+n6,2)*Q(n1+n7+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7+n8,3)+2.*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n1+n7+n8,3)+2.*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n1+n7+n8,3)
      - 4.*Q(n2+n4+n5,3)*Q(n3+n6,2)*Q(n1+n7+n8,3)+4.*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n1+n7+n8,3)
      - 4.*Q(n4+n5,2)*Q(n2+n3+n6,3)*Q(n1+n7+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n7+n8,3)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n1+n7+n8,3)-4.*Q(n2+n3+n5,3)*Q(n4+n6,2)*Q(n1+n7+n8,3)
      + 4.*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n1+n7+n8,3)-4.*Q(n3+n5,2)*Q(n2+n4+n6,3)*Q(n1+n7+n8,3)
      + 4.*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n1+n7+n8,3)-4.*Q(n2+n5,2)*Q(n3+n4+n6,3)*Q(n1+n7+n8,3)
      - 12.*Q(n5,1)*Q(n2+n3+n4+n6,4)*Q(n1+n7+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n7+n8,3)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n1+n7+n8,3)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n1+n7+n8,3)-4.*Q(n2+n3+n4,3)*Q(n5+n6,2)*Q(n1+n7+n8,3)
      + 4.*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n1+n7+n8,3)-4.*Q(n3+n4,2)*Q(n2+n5+n6,3)*Q(n1+n7+n8,3)
      + 4.*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n1+n7+n8,3)-4.*Q(n2+n4,2)*Q(n3+n5+n6,3)*Q(n1+n7+n8,3)
      - 12.*Q(n4,1)*Q(n2+n3+n5+n6,4)*Q(n1+n7+n8,3)+4.*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n1+n7+n8,3)
      - 4.*Q(n2+n3,2)*Q(n4+n5+n6,3)*Q(n1+n7+n8,3)-12.*Q(n3,1)*Q(n2+n4+n5+n6,4)*Q(n1+n7+n8,3)
      - 12.*Q(n2,1)*Q(n3+n4+n5+n6,4)*Q(n1+n7+n8,3)+48.*Q(n2+n3+n4+n5+n6,5)*Q(n1+n7+n8,3)
      + 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7+n8,3)-2.*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n7+n8,3)
      - 2.*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7+n8,3)-2.*Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n7+n8,3)
      + 4.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)*Q(n2+n7+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7+n8,3)
      + 2.*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n7+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7+n8,3)
      + 2.*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n7+n8,3)+4.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n2+n7+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7+n8,3)+2.*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n7+n8,3)
      + 4.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n2+n7+n8,3)+4.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n6,1)*Q(n2+n7+n8,3)
      - 12.*Q(n1+n3+n4+n5,4)*Q(n6,1)*Q(n2+n7+n8,3)-2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n7+n8,3)+2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n2+n7+n8,3)-4.*Q(n3+n4+n5,3)*Q(n1+n6,2)*Q(n2+n7+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7+n8,3)+2.*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n2+n7+n8,3)+2.*Q(n1,1)*Q(n4+n5,2)*Q(n3+n6,2)*Q(n2+n7+n8,3)
      - 4.*Q(n1+n4+n5,3)*Q(n3+n6,2)*Q(n2+n7+n8,3)+4.*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n2+n7+n8,3)
      - 4.*Q(n4+n5,2)*Q(n1+n3+n6,3)*Q(n2+n7+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n7+n8,3)+2.*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n1,1)*Q(n3+n5,2)*Q(n4+n6,2)*Q(n2+n7+n8,3)-4.*Q(n1+n3+n5,3)*Q(n4+n6,2)*Q(n2+n7+n8,3)
      + 4.*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n2+n7+n8,3)-4.*Q(n3+n5,2)*Q(n1+n4+n6,3)*Q(n2+n7+n8,3)
      + 4.*Q(n1,1)*Q(n5,1)*Q(n3+n4+n6,3)*Q(n2+n7+n8,3)-4.*Q(n1+n5,2)*Q(n3+n4+n6,3)*Q(n2+n7+n8,3)
      - 12.*Q(n5,1)*Q(n1+n3+n4+n6,4)*Q(n2+n7+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n7+n8,3)+2.*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n2+n7+n8,3)
      + 2.*Q(n1,1)*Q(n3+n4,2)*Q(n5+n6,2)*Q(n2+n7+n8,3)-4.*Q(n1+n3+n4,3)*Q(n5+n6,2)*Q(n2+n7+n8,3)
      + 4.*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n2+n7+n8,3)-4.*Q(n3+n4,2)*Q(n1+n5+n6,3)*Q(n2+n7+n8,3)
      + 4.*Q(n1,1)*Q(n4,1)*Q(n3+n5+n6,3)*Q(n2+n7+n8,3)-4.*Q(n1+n4,2)*Q(n3+n5+n6,3)*Q(n2+n7+n8,3)
      - 12.*Q(n4,1)*Q(n1+n3+n5+n6,4)*Q(n2+n7+n8,3)+4.*Q(n1,1)*Q(n3,1)*Q(n4+n5+n6,3)*Q(n2+n7+n8,3)
      - 4.*Q(n1+n3,2)*Q(n4+n5+n6,3)*Q(n2+n7+n8,3)-12.*Q(n3,1)*Q(n1+n4+n5+n6,4)*Q(n2+n7+n8,3)
      - 12.*Q(n1,1)*Q(n3+n4+n5+n6,4)*Q(n2+n7+n8,3)+48.*Q(n1+n3+n4+n5+n6,5)*Q(n2+n7+n8,3)
      - 6.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7+n8,4)+6.*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n7+n8,4)
      + 6.*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n2+n7+n8,4)+6.*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n2+n7+n8,4)
      - 12.*Q(n3+n4+n5,3)*Q(n6,1)*Q(n1+n2+n7+n8,4)+6.*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n2+n7+n8,4)
      - 6.*Q(n4+n5,2)*Q(n3+n6,2)*Q(n1+n2+n7+n8,4)+6.*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n2+n7+n8,4)
      - 6.*Q(n3+n5,2)*Q(n4+n6,2)*Q(n1+n2+n7+n8,4)-12.*Q(n5,1)*Q(n3+n4+n6,3)*Q(n1+n2+n7+n8,4)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n2+n7+n8,4)-6.*Q(n3+n4,2)*Q(n5+n6,2)*Q(n1+n2+n7+n8,4)
      - 12.*Q(n4,1)*Q(n3+n5+n6,3)*Q(n1+n2+n7+n8,4)-12.*Q(n3,1)*Q(n4+n5+n6,3)*Q(n1+n2+n7+n8,4)
      + 36.*Q(n3+n4+n5+n6,4)*Q(n1+n2+n7+n8,4)+2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7+n8,3)
      - 2.*Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n3+n7+n8,3)-2.*Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n3+n7+n8,3)+4.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)*Q(n3+n7+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7+n8,3)+2.*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n7+n8,3)
      - 2.*Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7+n8,3)+2.*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n7+n8,3)
      + 4.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n3+n7+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7+n8,3)
      + 2.*Q(n1+n2,2)*Q(n4+n5,2)*Q(n6,1)*Q(n3+n7+n8,3)+4.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n6,1)*Q(n3+n7+n8,3)
      + 4.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n6,1)*Q(n3+n7+n8,3)-12.*Q(n1+n2+n4+n5,4)*Q(n6,1)*Q(n3+n7+n8,3)
      - 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7+n8,3)+2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n7+n8,3)
      + 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n3+n7+n8,3)+2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n6,2)*Q(n3+n7+n8,3)
      - 4.*Q(n2+n4+n5,3)*Q(n1+n6,2)*Q(n3+n7+n8,3)-2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7+n8,3)
      + 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n7+n8,3)+2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n3+n7+n8,3)
      + 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n6,2)*Q(n3+n7+n8,3)-4.*Q(n1+n4+n5,3)*Q(n2+n6,2)*Q(n3+n7+n8,3)
      + 4.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n3+n7+n8,3)-4.*Q(n4+n5,2)*Q(n1+n2+n6,3)*Q(n3+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7+n8,3)+2.*Q(n1+n2,2)*Q(n5,1)*Q(n4+n6,2)*Q(n3+n7+n8,3)
      + 2.*Q(n2,1)*Q(n1+n5,2)*Q(n4+n6,2)*Q(n3+n7+n8,3)+2.*Q(n1,1)*Q(n2+n5,2)*Q(n4+n6,2)*Q(n3+n7+n8,3)
      - 4.*Q(n1+n2+n5,3)*Q(n4+n6,2)*Q(n3+n7+n8,3)+4.*Q(n2,1)*Q(n5,1)*Q(n1+n4+n6,3)*Q(n3+n7+n8,3)
      - 4.*Q(n2+n5,2)*Q(n1+n4+n6,3)*Q(n3+n7+n8,3)+4.*Q(n1,1)*Q(n5,1)*Q(n2+n4+n6,3)*Q(n3+n7+n8,3)
      - 4.*Q(n1+n5,2)*Q(n2+n4+n6,3)*Q(n3+n7+n8,3)-12.*Q(n5,1)*Q(n1+n2+n4+n6,4)*Q(n3+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7+n8,3)+2.*Q(n1+n2,2)*Q(n4,1)*Q(n5+n6,2)*Q(n3+n7+n8,3)
      + 2.*Q(n2,1)*Q(n1+n4,2)*Q(n5+n6,2)*Q(n3+n7+n8,3)+2.*Q(n1,1)*Q(n2+n4,2)*Q(n5+n6,2)*Q(n3+n7+n8,3)
      - 4.*Q(n1+n2+n4,3)*Q(n5+n6,2)*Q(n3+n7+n8,3)+4.*Q(n2,1)*Q(n4,1)*Q(n1+n5+n6,3)*Q(n3+n7+n8,3)
      - 4.*Q(n2+n4,2)*Q(n1+n5+n6,3)*Q(n3+n7+n8,3)+4.*Q(n1,1)*Q(n4,1)*Q(n2+n5+n6,3)*Q(n3+n7+n8,3)
      - 4.*Q(n1+n4,2)*Q(n2+n5+n6,3)*Q(n3+n7+n8,3)-12.*Q(n4,1)*Q(n1+n2+n5+n6,4)*Q(n3+n7+n8,3)
      + 4.*Q(n1,1)*Q(n2,1)*Q(n4+n5+n6,3)*Q(n3+n7+n8,3)-4.*Q(n1+n2,2)*Q(n4+n5+n6,3)*Q(n3+n7+n8,3)
      - 12.*Q(n2,1)*Q(n1+n4+n5+n6,4)*Q(n3+n7+n8,3)-12.*Q(n1,1)*Q(n2+n4+n5+n6,4)*Q(n3+n7+n8,3)
      + 48.*Q(n1+n2+n4+n5+n6,5)*Q(n3+n7+n8,3)-6.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7+n8,4)
      + 6.*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n7+n8,4)+6.*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n3+n7+n8,4)
      + 6.*Q(n2,1)*Q(n4+n5,2)*Q(n6,1)*Q(n1+n3+n7+n8,4)-12.*Q(n2+n4+n5,3)*Q(n6,1)*Q(n1+n3+n7+n8,4)
      + 6.*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n3+n7+n8,4)-6.*Q(n4+n5,2)*Q(n2+n6,2)*Q(n1+n3+n7+n8,4)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n4+n6,2)*Q(n1+n3+n7+n8,4)-6.*Q(n2+n5,2)*Q(n4+n6,2)*Q(n1+n3+n7+n8,4)
      - 12.*Q(n5,1)*Q(n2+n4+n6,3)*Q(n1+n3+n7+n8,4)+6.*Q(n2,1)*Q(n4,1)*Q(n5+n6,2)*Q(n1+n3+n7+n8,4)
      - 6.*Q(n2+n4,2)*Q(n5+n6,2)*Q(n1+n3+n7+n8,4)-12.*Q(n4,1)*Q(n2+n5+n6,3)*Q(n1+n3+n7+n8,4)
      - 12.*Q(n2,1)*Q(n4+n5+n6,3)*Q(n1+n3+n7+n8,4)+36.*Q(n2+n4+n5+n6,4)*Q(n1+n3+n7+n8,4)
      - 6.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7+n8,4)+6.*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n7+n8,4)
      + 6.*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n3+n7+n8,4)+6.*Q(n1,1)*Q(n4+n5,2)*Q(n6,1)*Q(n2+n3+n7+n8,4)
      - 12.*Q(n1+n4+n5,3)*Q(n6,1)*Q(n2+n3+n7+n8,4)+6.*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n3+n7+n8,4)
      - 6.*Q(n4+n5,2)*Q(n1+n6,2)*Q(n2+n3+n7+n8,4)+6.*Q(n1,1)*Q(n5,1)*Q(n4+n6,2)*Q(n2+n3+n7+n8,4)
      - 6.*Q(n1+n5,2)*Q(n4+n6,2)*Q(n2+n3+n7+n8,4)-12.*Q(n5,1)*Q(n1+n4+n6,3)*Q(n2+n3+n7+n8,4)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n5+n6,2)*Q(n2+n3+n7+n8,4)-6.*Q(n1+n4,2)*Q(n5+n6,2)*Q(n2+n3+n7+n8,4)
      - 12.*Q(n4,1)*Q(n1+n5+n6,3)*Q(n2+n3+n7+n8,4)-12.*Q(n1,1)*Q(n4+n5+n6,3)*Q(n2+n3+n7+n8,4)
      + 36.*Q(n1+n4+n5+n6,4)*Q(n2+n3+n7+n8,4)+24.*Q(n4,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n3+n7+n8,5)
      - 24.*Q(n4+n5,2)*Q(n6,1)*Q(n1+n2+n3+n7+n8,5)-24.*Q(n5,1)*Q(n4+n6,2)*Q(n1+n2+n3+n7+n8,5)
      - 24.*Q(n4,1)*Q(n5+n6,2)*Q(n1+n2+n3+n7+n8,5)+48.*Q(n4+n5+n6,3)*Q(n1+n2+n3+n7+n8,5)
      + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7+n8,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n4+n7+n8,3)
      - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7+n8,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n4+n7+n8,3)
      + 4.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n6,1)*Q(n4+n7+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7+n8,3)
      + 2.*Q(n2+n3,2)*Q(n1+n5,2)*Q(n6,1)*Q(n4+n7+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7+n8,3)
      + 2.*Q(n1+n3,2)*Q(n2+n5,2)*Q(n6,1)*Q(n4+n7+n8,3)+4.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n6,1)*Q(n4+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7+n8,3)+2.*Q(n1+n2,2)*Q(n3+n5,2)*Q(n6,1)*Q(n4+n7+n8,3)
      + 4.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n6,1)*Q(n4+n7+n8,3)+4.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n6,1)*Q(n4+n7+n8,3)
      - 12.*Q(n1+n2+n3+n5,4)*Q(n6,1)*Q(n4+n7+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n6,2)*Q(n4+n7+n8,3)+2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n6,2)*Q(n4+n7+n8,3)-4.*Q(n2+n3+n5,3)*Q(n1+n6,2)*Q(n4+n7+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7+n8,3)+2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n6,2)*Q(n4+n7+n8,3)+2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n6,2)*Q(n4+n7+n8,3)
      - 4.*Q(n1+n3+n5,3)*Q(n2+n6,2)*Q(n4+n7+n8,3)+4.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n6,3)*Q(n4+n7+n8,3)
      - 4.*Q(n3+n5,2)*Q(n1+n2+n6,3)*Q(n4+n7+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n6,2)*Q(n4+n7+n8,3)+2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n6,2)*Q(n4+n7+n8,3)-4.*Q(n1+n2+n5,3)*Q(n3+n6,2)*Q(n4+n7+n8,3)
      + 4.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n6,3)*Q(n4+n7+n8,3)-4.*Q(n2+n5,2)*Q(n1+n3+n6,3)*Q(n4+n7+n8,3)
      + 4.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n6,3)*Q(n4+n7+n8,3)-4.*Q(n1+n5,2)*Q(n2+n3+n6,3)*Q(n4+n7+n8,3)
      - 12.*Q(n5,1)*Q(n1+n2+n3+n6,4)*Q(n4+n7+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n5+n6,2)*Q(n4+n7+n8,3)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n5+n6,2)*Q(n4+n7+n8,3)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n5+n6,2)*Q(n4+n7+n8,3)-4.*Q(n1+n2+n3,3)*Q(n5+n6,2)*Q(n4+n7+n8,3)
      + 4.*Q(n2,1)*Q(n3,1)*Q(n1+n5+n6,3)*Q(n4+n7+n8,3)-4.*Q(n2+n3,2)*Q(n1+n5+n6,3)*Q(n4+n7+n8,3)
      + 4.*Q(n1,1)*Q(n3,1)*Q(n2+n5+n6,3)*Q(n4+n7+n8,3)-4.*Q(n1+n3,2)*Q(n2+n5+n6,3)*Q(n4+n7+n8,3)
      - 12.*Q(n3,1)*Q(n1+n2+n5+n6,4)*Q(n4+n7+n8,3)+4.*Q(n1,1)*Q(n2,1)*Q(n3+n5+n6,3)*Q(n4+n7+n8,3)
      - 4.*Q(n1+n2,2)*Q(n3+n5+n6,3)*Q(n4+n7+n8,3)-12.*Q(n2,1)*Q(n1+n3+n5+n6,4)*Q(n4+n7+n8,3)
      - 12.*Q(n1,1)*Q(n2+n3+n5+n6,4)*Q(n4+n7+n8,3)+48.*Q(n1+n2+n3+n5+n6,5)*Q(n4+n7+n8,3)
      - 6.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7+n8,4)+6.*Q(n2+n3,2)*Q(n5,1)*Q(n6,1)*Q(n1+n4+n7+n8,4)
      + 6.*Q(n3,1)*Q(n2+n5,2)*Q(n6,1)*Q(n1+n4+n7+n8,4)+6.*Q(n2,1)*Q(n3+n5,2)*Q(n6,1)*Q(n1+n4+n7+n8,4)
      - 12.*Q(n2+n3+n5,3)*Q(n6,1)*Q(n1+n4+n7+n8,4)+6.*Q(n3,1)*Q(n5,1)*Q(n2+n6,2)*Q(n1+n4+n7+n8,4)
      - 6.*Q(n3+n5,2)*Q(n2+n6,2)*Q(n1+n4+n7+n8,4)+6.*Q(n2,1)*Q(n5,1)*Q(n3+n6,2)*Q(n1+n4+n7+n8,4)
      - 6.*Q(n2+n5,2)*Q(n3+n6,2)*Q(n1+n4+n7+n8,4)-12.*Q(n5,1)*Q(n2+n3+n6,3)*Q(n1+n4+n7+n8,4)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n5+n6,2)*Q(n1+n4+n7+n8,4)-6.*Q(n2+n3,2)*Q(n5+n6,2)*Q(n1+n4+n7+n8,4)
      - 12.*Q(n3,1)*Q(n2+n5+n6,3)*Q(n1+n4+n7+n8,4)-12.*Q(n2,1)*Q(n3+n5+n6,3)*Q(n1+n4+n7+n8,4)
      + 36.*Q(n2+n3+n5+n6,4)*Q(n1+n4+n7+n8,4)-6.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7+n8,4)
      + 6.*Q(n1+n3,2)*Q(n5,1)*Q(n6,1)*Q(n2+n4+n7+n8,4)+6.*Q(n3,1)*Q(n1+n5,2)*Q(n6,1)*Q(n2+n4+n7+n8,4)
      + 6.*Q(n1,1)*Q(n3+n5,2)*Q(n6,1)*Q(n2+n4+n7+n8,4)-12.*Q(n1+n3+n5,3)*Q(n6,1)*Q(n2+n4+n7+n8,4)
      + 6.*Q(n3,1)*Q(n5,1)*Q(n1+n6,2)*Q(n2+n4+n7+n8,4)-6.*Q(n3+n5,2)*Q(n1+n6,2)*Q(n2+n4+n7+n8,4)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n3+n6,2)*Q(n2+n4+n7+n8,4)-6.*Q(n1+n5,2)*Q(n3+n6,2)*Q(n2+n4+n7+n8,4)
      - 12.*Q(n5,1)*Q(n1+n3+n6,3)*Q(n2+n4+n7+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n5+n6,2)*Q(n2+n4+n7+n8,4)
      - 6.*Q(n1+n3,2)*Q(n5+n6,2)*Q(n2+n4+n7+n8,4)-12.*Q(n3,1)*Q(n1+n5+n6,3)*Q(n2+n4+n7+n8,4)
      - 12.*Q(n1,1)*Q(n3+n5+n6,3)*Q(n2+n4+n7+n8,4)+36.*Q(n1+n3+n5+n6,4)*Q(n2+n4+n7+n8,4)
      + 24.*Q(n3,1)*Q(n5,1)*Q(n6,1)*Q(n1+n2+n4+n7+n8,5)-24.*Q(n3+n5,2)*Q(n6,1)*Q(n1+n2+n4+n7+n8,5)
      - 24.*Q(n5,1)*Q(n3+n6,2)*Q(n1+n2+n4+n7+n8,5)-24.*Q(n3,1)*Q(n5+n6,2)*Q(n1+n2+n4+n7+n8,5)
      + 48.*Q(n3+n5+n6,3)*Q(n1+n2+n4+n7+n8,5)-6.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7+n8,4)
      + 6.*Q(n1+n2,2)*Q(n5,1)*Q(n6,1)*Q(n3+n4+n7+n8,4)+6.*Q(n2,1)*Q(n1+n5,2)*Q(n6,1)*Q(n3+n4+n7+n8,4)
      + 6.*Q(n1,1)*Q(n2+n5,2)*Q(n6,1)*Q(n3+n4+n7+n8,4)-12.*Q(n1+n2+n5,3)*Q(n6,1)*Q(n3+n4+n7+n8,4)
      + 6.*Q(n2,1)*Q(n5,1)*Q(n1+n6,2)*Q(n3+n4+n7+n8,4)-6.*Q(n2+n5,2)*Q(n1+n6,2)*Q(n3+n4+n7+n8,4)
      + 6.*Q(n1,1)*Q(n5,1)*Q(n2+n6,2)*Q(n3+n4+n7+n8,4)-6.*Q(n1+n5,2)*Q(n2+n6,2)*Q(n3+n4+n7+n8,4)
      - 12.*Q(n5,1)*Q(n1+n2+n6,3)*Q(n3+n4+n7+n8,4)+6.*Q(n1,1)*Q(n2,1)*Q(n5+n6,2)*Q(n3+n4+n7+n8,4)
      - 6.*Q(n1+n2,2)*Q(n5+n6,2)*Q(n3+n4+n7+n8,4)-12.*Q(n2,1)*Q(n1+n5+n6,3)*Q(n3+n4+n7+n8,4)
      - 12.*Q(n1,1)*Q(n2+n5+n6,3)*Q(n3+n4+n7+n8,4)+36.*Q(n1+n2+n5+n6,4)*Q(n3+n4+n7+n8,4)
      + 24.*Q(n2,1)*Q(n5,1)*Q(n6,1)*Q(n1+n3+n4+n7+n8,5)-24.*Q(n2+n5,2)*Q(n6,1)*Q(n1+n3+n4+n7+n8,5)
      - 24.*Q(n5,1)*Q(n2+n6,2)*Q(n1+n3+n4+n7+n8,5)-24.*Q(n2,1)*Q(n5+n6,2)*Q(n1+n3+n4+n7+n8,5)
      + 48.*Q(n2+n5+n6,3)*Q(n1+n3+n4+n7+n8,5)+24.*Q(n1,1)*Q(n5,1)*Q(n6,1)*Q(n2+n3+n4+n7+n8,5)
      - 24.*Q(n1+n5,2)*Q(n6,1)*Q(n2+n3+n4+n7+n8,5)-24.*Q(n5,1)*Q(n1+n6,2)*Q(n2+n3+n4+n7+n8,5)
      - 24.*Q(n1,1)*Q(n5+n6,2)*Q(n2+n3+n4+n7+n8,5)+48.*Q(n1+n5+n6,3)*Q(n2+n3+n4+n7+n8,5)
      - 120.*Q(n5,1)*Q(n6,1)*Q(n1+n2+n3+n4+n7+n8,6)+120.*Q(n5+n6,2)*Q(n1+n2+n3+n4+n7+n8,6)
      + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7+n8,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n5+n7+n8,3)
      - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7+n8,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n5+n7+n8,3)
      + 4.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n6,1)*Q(n5+n7+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7+n8,3)
      + 2.*Q(n2+n3,2)*Q(n1+n4,2)*Q(n6,1)*Q(n5+n7+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7+n8,3)
      + 2.*Q(n1+n3,2)*Q(n2+n4,2)*Q(n6,1)*Q(n5+n7+n8,3)+4.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n6,1)*Q(n5+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7+n8,3)+2.*Q(n1+n2,2)*Q(n3+n4,2)*Q(n6,1)*Q(n5+n7+n8,3)
      + 4.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n6,1)*Q(n5+n7+n8,3)+4.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n6,1)*Q(n5+n7+n8,3)
      - 12.*Q(n1+n2+n3+n4,4)*Q(n6,1)*Q(n5+n7+n8,3)-2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n6,2)*Q(n5+n7+n8,3)+2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n6,2)*Q(n5+n7+n8,3)-4.*Q(n2+n3+n4,3)*Q(n1+n6,2)*Q(n5+n7+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7+n8,3)+2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n6,2)*Q(n5+n7+n8,3)+2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n6,2)*Q(n5+n7+n8,3)
      - 4.*Q(n1+n3+n4,3)*Q(n2+n6,2)*Q(n5+n7+n8,3)+4.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n6,3)*Q(n5+n7+n8,3)
      - 4.*Q(n3+n4,2)*Q(n1+n2+n6,3)*Q(n5+n7+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n6,2)*Q(n5+n7+n8,3)+2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n6,2)*Q(n5+n7+n8,3)-4.*Q(n1+n2+n4,3)*Q(n3+n6,2)*Q(n5+n7+n8,3)
      + 4.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n6,3)*Q(n5+n7+n8,3)-4.*Q(n2+n4,2)*Q(n1+n3+n6,3)*Q(n5+n7+n8,3)
      + 4.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n6,3)*Q(n5+n7+n8,3)-4.*Q(n1+n4,2)*Q(n2+n3+n6,3)*Q(n5+n7+n8,3)
      - 12.*Q(n4,1)*Q(n1+n2+n3+n6,4)*Q(n5+n7+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n6,2)*Q(n5+n7+n8,3)+2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n6,2)*Q(n5+n7+n8,3)
      + 2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n6,2)*Q(n5+n7+n8,3)-4.*Q(n1+n2+n3,3)*Q(n4+n6,2)*Q(n5+n7+n8,3)
      + 4.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n6,3)*Q(n5+n7+n8,3)-4.*Q(n2+n3,2)*Q(n1+n4+n6,3)*Q(n5+n7+n8,3)
      + 4.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n6,3)*Q(n5+n7+n8,3)-4.*Q(n1+n3,2)*Q(n2+n4+n6,3)*Q(n5+n7+n8,3)
      - 12.*Q(n3,1)*Q(n1+n2+n4+n6,4)*Q(n5+n7+n8,3)+4.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n6,3)*Q(n5+n7+n8,3)
      - 4.*Q(n1+n2,2)*Q(n3+n4+n6,3)*Q(n5+n7+n8,3)-12.*Q(n2,1)*Q(n1+n3+n4+n6,4)*Q(n5+n7+n8,3)
      - 12.*Q(n1,1)*Q(n2+n3+n4+n6,4)*Q(n5+n7+n8,3)+48.*Q(n1+n2+n3+n4+n6,5)*Q(n5+n7+n8,3)
      - 6.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7+n8,4)+6.*Q(n2+n3,2)*Q(n4,1)*Q(n6,1)*Q(n1+n5+n7+n8,4)
      + 6.*Q(n3,1)*Q(n2+n4,2)*Q(n6,1)*Q(n1+n5+n7+n8,4)+6.*Q(n2,1)*Q(n3+n4,2)*Q(n6,1)*Q(n1+n5+n7+n8,4)
      - 12.*Q(n2+n3+n4,3)*Q(n6,1)*Q(n1+n5+n7+n8,4)+6.*Q(n3,1)*Q(n4,1)*Q(n2+n6,2)*Q(n1+n5+n7+n8,4)
      - 6.*Q(n3+n4,2)*Q(n2+n6,2)*Q(n1+n5+n7+n8,4)+6.*Q(n2,1)*Q(n4,1)*Q(n3+n6,2)*Q(n1+n5+n7+n8,4)
      - 6.*Q(n2+n4,2)*Q(n3+n6,2)*Q(n1+n5+n7+n8,4)-12.*Q(n4,1)*Q(n2+n3+n6,3)*Q(n1+n5+n7+n8,4)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n4+n6,2)*Q(n1+n5+n7+n8,4)-6.*Q(n2+n3,2)*Q(n4+n6,2)*Q(n1+n5+n7+n8,4)
      - 12.*Q(n3,1)*Q(n2+n4+n6,3)*Q(n1+n5+n7+n8,4)-12.*Q(n2,1)*Q(n3+n4+n6,3)*Q(n1+n5+n7+n8,4)
      + 36.*Q(n2+n3+n4+n6,4)*Q(n1+n5+n7+n8,4)-6.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7+n8,4)
      + 6.*Q(n1+n3,2)*Q(n4,1)*Q(n6,1)*Q(n2+n5+n7+n8,4)+6.*Q(n3,1)*Q(n1+n4,2)*Q(n6,1)*Q(n2+n5+n7+n8,4)
      + 6.*Q(n1,1)*Q(n3+n4,2)*Q(n6,1)*Q(n2+n5+n7+n8,4)-12.*Q(n1+n3+n4,3)*Q(n6,1)*Q(n2+n5+n7+n8,4)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n1+n6,2)*Q(n2+n5+n7+n8,4)-6.*Q(n3+n4,2)*Q(n1+n6,2)*Q(n2+n5+n7+n8,4)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n3+n6,2)*Q(n2+n5+n7+n8,4)-6.*Q(n1+n4,2)*Q(n3+n6,2)*Q(n2+n5+n7+n8,4)
      - 12.*Q(n4,1)*Q(n1+n3+n6,3)*Q(n2+n5+n7+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n4+n6,2)*Q(n2+n5+n7+n8,4)
      - 6.*Q(n1+n3,2)*Q(n4+n6,2)*Q(n2+n5+n7+n8,4)-12.*Q(n3,1)*Q(n1+n4+n6,3)*Q(n2+n5+n7+n8,4)
      - 12.*Q(n1,1)*Q(n3+n4+n6,3)*Q(n2+n5+n7+n8,4)+36.*Q(n1+n3+n4+n6,4)*Q(n2+n5+n7+n8,4)
      + 24.*Q(n3,1)*Q(n4,1)*Q(n6,1)*Q(n1+n2+n5+n7+n8,5)-24.*Q(n3+n4,2)*Q(n6,1)*Q(n1+n2+n5+n7+n8,5)
      - 24.*Q(n4,1)*Q(n3+n6,2)*Q(n1+n2+n5+n7+n8,5)-24.*Q(n3,1)*Q(n4+n6,2)*Q(n1+n2+n5+n7+n8,5)
      + 48.*Q(n3+n4+n6,3)*Q(n1+n2+n5+n7+n8,5)-6.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7+n8,4)
      + 6.*Q(n1+n2,2)*Q(n4,1)*Q(n6,1)*Q(n3+n5+n7+n8,4)+6.*Q(n2,1)*Q(n1+n4,2)*Q(n6,1)*Q(n3+n5+n7+n8,4)
      + 6.*Q(n1,1)*Q(n2+n4,2)*Q(n6,1)*Q(n3+n5+n7+n8,4)-12.*Q(n1+n2+n4,3)*Q(n6,1)*Q(n3+n5+n7+n8,4)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n1+n6,2)*Q(n3+n5+n7+n8,4)-6.*Q(n2+n4,2)*Q(n1+n6,2)*Q(n3+n5+n7+n8,4)
      + 6.*Q(n1,1)*Q(n4,1)*Q(n2+n6,2)*Q(n3+n5+n7+n8,4)-6.*Q(n1+n4,2)*Q(n2+n6,2)*Q(n3+n5+n7+n8,4)
      - 12.*Q(n4,1)*Q(n1+n2+n6,3)*Q(n3+n5+n7+n8,4)+6.*Q(n1,1)*Q(n2,1)*Q(n4+n6,2)*Q(n3+n5+n7+n8,4)
      - 6.*Q(n1+n2,2)*Q(n4+n6,2)*Q(n3+n5+n7+n8,4)-12.*Q(n2,1)*Q(n1+n4+n6,3)*Q(n3+n5+n7+n8,4)
      - 12.*Q(n1,1)*Q(n2+n4+n6,3)*Q(n3+n5+n7+n8,4)+36.*Q(n1+n2+n4+n6,4)*Q(n3+n5+n7+n8,4)
      + 24.*Q(n2,1)*Q(n4,1)*Q(n6,1)*Q(n1+n3+n5+n7+n8,5)-24.*Q(n2+n4,2)*Q(n6,1)*Q(n1+n3+n5+n7+n8,5)
      - 24.*Q(n4,1)*Q(n2+n6,2)*Q(n1+n3+n5+n7+n8,5)-24.*Q(n2,1)*Q(n4+n6,2)*Q(n1+n3+n5+n7+n8,5)
      + 48.*Q(n2+n4+n6,3)*Q(n1+n3+n5+n7+n8,5)+24.*Q(n1,1)*Q(n4,1)*Q(n6,1)*Q(n2+n3+n5+n7+n8,5)
      - 24.*Q(n1+n4,2)*Q(n6,1)*Q(n2+n3+n5+n7+n8,5)-24.*Q(n4,1)*Q(n1+n6,2)*Q(n2+n3+n5+n7+n8,5)
      - 24.*Q(n1,1)*Q(n4+n6,2)*Q(n2+n3+n5+n7+n8,5)+48.*Q(n1+n4+n6,3)*Q(n2+n3+n5+n7+n8,5)
      - 120.*Q(n4,1)*Q(n6,1)*Q(n1+n2+n3+n5+n7+n8,6)+120.*Q(n4+n6,2)*Q(n1+n2+n3+n5+n7+n8,6)
      - 6.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7+n8,4)+6.*Q(n1+n2,2)*Q(n3,1)*Q(n6,1)*Q(n4+n5+n7+n8,4)
      + 6.*Q(n2,1)*Q(n1+n3,2)*Q(n6,1)*Q(n4+n5+n7+n8,4)+6.*Q(n1,1)*Q(n2+n3,2)*Q(n6,1)*Q(n4+n5+n7+n8,4)
      - 12.*Q(n1+n2+n3,3)*Q(n6,1)*Q(n4+n5+n7+n8,4)+6.*Q(n2,1)*Q(n3,1)*Q(n1+n6,2)*Q(n4+n5+n7+n8,4)
      - 6.*Q(n2+n3,2)*Q(n1+n6,2)*Q(n4+n5+n7+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n2+n6,2)*Q(n4+n5+n7+n8,4)
      - 6.*Q(n1+n3,2)*Q(n2+n6,2)*Q(n4+n5+n7+n8,4)-12.*Q(n3,1)*Q(n1+n2+n6,3)*Q(n4+n5+n7+n8,4)
      + 6.*Q(n1,1)*Q(n2,1)*Q(n3+n6,2)*Q(n4+n5+n7+n8,4)-6.*Q(n1+n2,2)*Q(n3+n6,2)*Q(n4+n5+n7+n8,4)
      - 12.*Q(n2,1)*Q(n1+n3+n6,3)*Q(n4+n5+n7+n8,4)-12.*Q(n1,1)*Q(n2+n3+n6,3)*Q(n4+n5+n7+n8,4)
      + 36.*Q(n1+n2+n3+n6,4)*Q(n4+n5+n7+n8,4)+24.*Q(n2,1)*Q(n3,1)*Q(n6,1)*Q(n1+n4+n5+n7+n8,5)
      - 24.*Q(n2+n3,2)*Q(n6,1)*Q(n1+n4+n5+n7+n8,5)-24.*Q(n3,1)*Q(n2+n6,2)*Q(n1+n4+n5+n7+n8,5)
      - 24.*Q(n2,1)*Q(n3+n6,2)*Q(n1+n4+n5+n7+n8,5)+48.*Q(n2+n3+n6,3)*Q(n1+n4+n5+n7+n8,5)
      + 24.*Q(n1,1)*Q(n3,1)*Q(n6,1)*Q(n2+n4+n5+n7+n8,5)-24.*Q(n1+n3,2)*Q(n6,1)*Q(n2+n4+n5+n7+n8,5)
      - 24.*Q(n3,1)*Q(n1+n6,2)*Q(n2+n4+n5+n7+n8,5)-24.*Q(n1,1)*Q(n3+n6,2)*Q(n2+n4+n5+n7+n8,5)
      + 48.*Q(n1+n3+n6,3)*Q(n2+n4+n5+n7+n8,5)-120.*Q(n3,1)*Q(n6,1)*Q(n1+n2+n4+n5+n7+n8,6)
      + 120.*Q(n3+n6,2)*Q(n1+n2+n4+n5+n7+n8,6)+24.*Q(n1,1)*Q(n2,1)*Q(n6,1)*Q(n3+n4+n5+n7+n8,5)
      - 24.*Q(n1+n2,2)*Q(n6,1)*Q(n3+n4+n5+n7+n8,5)-24.*Q(n2,1)*Q(n1+n6,2)*Q(n3+n4+n5+n7+n8,5)
      - 24.*Q(n1,1)*Q(n2+n6,2)*Q(n3+n4+n5+n7+n8,5)+48.*Q(n1+n2+n6,3)*Q(n3+n4+n5+n7+n8,5)
      - 120.*Q(n2,1)*Q(n6,1)*Q(n1+n3+n4+n5+n7+n8,6)+120.*Q(n2+n6,2)*Q(n1+n3+n4+n5+n7+n8,6)
      - 120.*Q(n1,1)*Q(n6,1)*Q(n2+n3+n4+n5+n7+n8,6)+120.*Q(n1+n6,2)*Q(n2+n3+n4+n5+n7+n8,6)
      + 720.*Q(n6,1)*Q(n1+n2+n3+n4+n5+n7+n8,7)+2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7+n8,3)
      - 2.*Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6+n7+n8,3)-2.*Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6+n7+n8,3)+4.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6+n7+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7+n8,3)+2.*Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6+n7+n8,3)
      - 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7+n8,3)+2.*Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6+n7+n8,3)
      + 4.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6+n7+n8,3)-2.*Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7+n8,3)
      + 2.*Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6+n7+n8,3)+4.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6+n7+n8,3)
      + 4.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6+n7+n8,3)-12.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6+n7+n8,3)
      - 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7+n8,3)+2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6+n7+n8,3)
      + 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6+n7+n8,3)+2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6+n7+n8,3)
      - 4.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6+n7+n8,3)-2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7+n8,3)
      + 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6+n7+n8,3)+2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6+n7+n8,3)
      + 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6+n7+n8,3)-4.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6+n7+n8,3)
      + 4.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6+n7+n8,3)-4.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7+n8,3)+2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6+n7+n8,3)
      + 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6+n7+n8,3)+2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6+n7+n8,3)
      - 4.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6+n7+n8,3)+4.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6+n7+n8,3)
      - 4.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6+n7+n8,3)+4.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6+n7+n8,3)
      - 4.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6+n7+n8,3)-12.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6+n7+n8,3)
      - 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7+n8,3)+2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6+n7+n8,3)
      + 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6+n7+n8,3)+2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6+n7+n8,3)
      - 4.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6+n7+n8,3)+4.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6+n7+n8,3)
      - 4.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6+n7+n8,3)+4.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6+n7+n8,3)
      - 4.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6+n7+n8,3)-12.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6+n7+n8,3)
      + 4.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6+n7+n8,3)-4.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6+n7+n8,3)
      - 12.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6+n7+n8,3)-12.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6+n7+n8,3)
      + 48.*Q(n1+n2+n3+n4+n5,5)*Q(n6+n7+n8,3)-6.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7+n8,4)
      + 6.*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6+n7+n8,4)+6.*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6+n7+n8,4)
      + 6.*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6+n7+n8,4)-12.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6+n7+n8,4)
      + 6.*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6+n7+n8,4)-6.*Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6+n7+n8,4)
      + 6.*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6+n7+n8,4)-6.*Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6+n7+n8,4)
      - 12.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6+n7+n8,4)+6.*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6+n7+n8,4)
      - 6.*Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6+n7+n8,4)-12.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6+n7+n8,4)
      - 12.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6+n7+n8,4)+36.*Q(n2+n3+n4+n5,4)*Q(n1+n6+n7+n8,4)
      - 6.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7+n8,4)+6.*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6+n7+n8,4)
      + 6.*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6+n7+n8,4)+6.*Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6+n7+n8,4)
      - 12.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6+n7+n8,4)+6.*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6+n7+n8,4)
      - 6.*Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6+n7+n8,4)+6.*Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6+n7+n8,4)
      - 6.*Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6+n7+n8,4)-12.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6+n7+n8,4)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6+n7+n8,4)-6.*Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6+n7+n8,4)
      - 12.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6+n7+n8,4)-12.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6+n7+n8,4)
      + 36.*Q(n1+n3+n4+n5,4)*Q(n2+n6+n7+n8,4)+24.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6+n7+n8,5)
      - 24.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6+n7+n8,5)-24.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6+n7+n8,5)
      - 24.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6+n7+n8,5)+48.*Q(n3+n4+n5,3)*Q(n1+n2+n6+n7+n8,5)
      - 6.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7+n8,4)+6.*Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6+n7+n8,4)
      + 6.*Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6+n7+n8,4)+6.*Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6+n7+n8,4)
      - 12.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6+n7+n8,4)+6.*Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6+n7+n8,4)
      - 6.*Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6+n7+n8,4)+6.*Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6+n7+n8,4)
      - 6.*Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6+n7+n8,4)-12.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6+n7+n8,4)
      + 6.*Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6+n7+n8,4)-6.*Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6+n7+n8,4)
      - 12.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6+n7+n8,4)-12.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6+n7+n8,4)
      + 36.*Q(n1+n2+n4+n5,4)*Q(n3+n6+n7+n8,4)+24.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6+n7+n8,5)
      - 24.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6+n7+n8,5)-24.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6+n7+n8,5)
      - 24.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6+n7+n8,5)+48.*Q(n2+n4+n5,3)*Q(n1+n3+n6+n7+n8,5)
      + 24.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6+n7+n8,5)-24.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6+n7+n8,5)
      - 24.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6+n7+n8,5)-24.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6+n7+n8,5)
      + 48.*Q(n1+n4+n5,3)*Q(n2+n3+n6+n7+n8,5)-120.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6+n7+n8,6)
      + 120.*Q(n4+n5,2)*Q(n1+n2+n3+n6+n7+n8,6)-6.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7+n8,4)
      + 6.*Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6+n7+n8,4)+6.*Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6+n7+n8,4)
      + 6.*Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6+n7+n8,4)-12.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6+n7+n8,4)
      + 6.*Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6+n7+n8,4)-6.*Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6+n7+n8,4)
      + 6.*Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6+n7+n8,4)-6.*Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6+n7+n8,4)
      - 12.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6+n7+n8,4)+6.*Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6+n7+n8,4)
      - 6.*Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6+n7+n8,4)-12.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6+n7+n8,4)
      - 12.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6+n7+n8,4)+36.*Q(n1+n2+n3+n5,4)*Q(n4+n6+n7+n8,4)
      + 24.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6+n7+n8,5)-24.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6+n7+n8,5)
      - 24.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6+n7+n8,5)-24.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6+n7+n8,5)
      + 48.*Q(n2+n3+n5,3)*Q(n1+n4+n6+n7+n8,5)+24.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6+n7+n8,5)
      - 24.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6+n7+n8,5)-24.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6+n7+n8,5)
      - 24.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6+n7+n8,5)+48.*Q(n1+n3+n5,3)*Q(n2+n4+n6+n7+n8,5)
      - 120.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6+n7+n8,6)+120.*Q(n3+n5,2)*Q(n1+n2+n4+n6+n7+n8,6)
      + 24.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6+n7+n8,5)-24.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6+n7+n8,5)
      - 24.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6+n7+n8,5)-24.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6+n7+n8,5)
      + 48.*Q(n1+n2+n5,3)*Q(n3+n4+n6+n7+n8,5)-120.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6+n7+n8,6)
      + 120.*Q(n2+n5,2)*Q(n1+n3+n4+n6+n7+n8,6)-120.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6+n7+n8,6)
      + 120.*Q(n1+n5,2)*Q(n2+n3+n4+n6+n7+n8,6)+720.*Q(n5,1)*Q(n1+n2+n3+n4+n6+n7+n8,7)
      - 6.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7+n8,4)+6.*Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6+n7+n8,4)
      + 6.*Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6+n7+n8,4)+6.*Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6+n7+n8,4)
      - 12.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6+n7+n8,4)+6.*Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6+n7+n8,4)
      - 6.*Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6+n7+n8,4)+6.*Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6+n7+n8,4)
      - 6.*Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6+n7+n8,4)-12.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6+n7+n8,4)
      + 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6+n7+n8,4)-6.*Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6+n7+n8,4)
      - 12.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6+n7+n8,4)-12.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6+n7+n8,4)
      + 36.*Q(n1+n2+n3+n4,4)*Q(n5+n6+n7+n8,4)+24.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6+n7+n8,5)
      - 24.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6+n7+n8,5)-24.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6+n7+n8,5)
      - 24.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6+n7+n8,5)+48.*Q(n2+n3+n4,3)*Q(n1+n5+n6+n7+n8,5)
      + 24.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6+n7+n8,5)-24.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6+n7+n8,5)
      - 24.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6+n7+n8,5)-24.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6+n7+n8,5)
      + 48.*Q(n1+n3+n4,3)*Q(n2+n5+n6+n7+n8,5)-120.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6+n7+n8,6)
      + 120.*Q(n3+n4,2)*Q(n1+n2+n5+n6+n7+n8,6)+24.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6+n7+n8,5)
      - 24.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6+n7+n8,5)-24.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6+n7+n8,5)
      - 24.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6+n7+n8,5)+48.*Q(n1+n2+n4,3)*Q(n3+n5+n6+n7+n8,5)
      - 120.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6+n7+n8,6)+120.*Q(n2+n4,2)*Q(n1+n3+n5+n6+n7+n8,6)
      - 120.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6+n7+n8,6)+120.*Q(n1+n4,2)*Q(n2+n3+n5+n6+n7+n8,6)
      + 720.*Q(n4,1)*Q(n1+n2+n3+n5+n6+n7+n8,7)+24.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6+n7+n8,5)
      - 24.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6+n7+n8,5)-24.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6+n7+n8,5)
      - 24.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6+n7+n8,5)+48.*Q(n1+n2+n3,3)*Q(n4+n5+n6+n7+n8,5)
      - 120.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6+n7+n8,6)+120.*Q(n2+n3,2)*Q(n1+n4+n5+n6+n7+n8,6)
      - 120.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6+n7+n8,6)+120.*Q(n1+n3,2)*Q(n2+n4+n5+n6+n7+n8,6)
      + 720.*Q(n3,1)*Q(n1+n2+n4+n5+n6+n7+n8,7)-120.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6+n7+n8,6)
      + 120.*Q(n1+n2,2)*Q(n3+n4+n5+n6+n7+n8,6)+720.*Q(n2,1)*Q(n1+n3+n4+n5+n6+n7+n8,7)
      + 720.*Q(n1,1)*Q(n2+n3+n4+n5+n6+n7+n8,7)-5040.*Q(n1+n2+n3+n4+n5+n6+n7+n8,8)
    );
} 
//=======================================================
TComplex AliAnalysisTaskChargedFlow::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{
 	return	
 		(
			Eight_1(n1,n2,n3,n4,n5,n6,n7,n8)
	    + Eight_2(n1,n2,n3,n4,n5,n6,n7,n8)  
	    + Eight_3(n1,n2,n3,n4,n5,n6,n7,n8)
      + Eight_4(n1,n2,n3,n4,n5,n6,n7,n8)
	  );
} 
*/
//=======================================================
Int_t AliAnalysisTaskChargedFlow::GetTPCMult(AliVEvent* ev) const
{
	
  // function to return the tpc only multiplicity as in the flow package; it is used to cut on the event multiplcity to remove outliers in the centrality
	
  Int_t multTPC = 0;
    
  if (fAnalysisType == "AOD") {
		
    AliAODEvent* aod = (AliAODEvent*)ev;
    const Int_t nGoodTracks = aod->GetNumberOfTracks();
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {     
      AliAODTrack* trackAOD =(AliAODTrack*) aod->GetTrack(iTracks);
			
      if (!trackAOD)
	continue;
			
      if (!(trackAOD->TestFilterBit(1)))
	continue;
			
      if ((trackAOD->Pt() < fMinPt) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > fEtaCut) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2))
	continue;
			
      multTPC++;
    }//track loop
		
  }//analysis type
	
  return multTPC;
}

//____________________________________________________________________
Int_t AliAnalysisTaskChargedFlow::GetGlobalMult(AliVEvent* ev) const
{
	
  // function to return the global multiplicity as in the flow package; it is used to cut on the event multiplcity to remove outliers in the centrality
	
  Int_t multGlobal = 0;
  if (fAnalysisType == "AOD") {
		
    AliAODEvent* aod = (AliAODEvent*)ev;
    const Int_t nGoodTracks = aod->GetNumberOfTracks();
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {     
      AliAODTrack* trackAOD = (AliAODTrack*)aod->GetTrack(iTracks);
			
      if (!trackAOD)
	continue;
			
      if (!(trackAOD->TestFilterBit(16)))
	continue;
			
      if ((trackAOD->Pt() < fMinPt) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > fEtaCut) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1))
	continue;
			
      
        // clone for constraining
        Double_t b[2] = { -99., -99.};
        Double_t bCov[3] = { -99., -99., -99.};
        
        AliAODTrack* trackAODC = new AliAODTrack(*trackAOD);
        if (!trackAODC) {
            AliWarning("Clone of AOD track failed.");
            delete trackAODC;
            continue;
        }
        
        if (!trackAODC->PropagateToDCA(aod->GetPrimaryVertex(), aod->GetMagneticField(), 100., b, bCov)){
            delete trackAODC;
            continue;
        } else {
            delete trackAODC;
        }

        
      if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3))
          continue;
			
      multGlobal++;
    }//track loop
		
  }//analysis type
	
  return multGlobal;
}


//_____________________________________________________________________________
Float_t AliAnalysisTaskChargedFlow::GetVertex(AliVEvent* ev) const
{

  Float_t vtxz = -999.;

  if (fAnalysisType == "AOD"){

    AliAODEvent* aod = (AliAODEvent*)ev;
      
      
      const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
      if (!trkVtx || trkVtx->GetNContributors()<=0)
          return vtxz;
      TString vtxTtl = trkVtx->GetTitle();
      if (!vtxTtl.Contains("VertexerTracks"))
          return vtxz;
      
      // comment out on Nov 30,2015 for a test of vertex cut
      const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
      if (!spdVtx || spdVtx->GetNContributors()<=0)
          return vtxz;
      if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5)
          return vtxz;
      
      vtxz = trkVtx->GetZ();

      
      
  }
  
  return vtxz;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskChargedFlow::GetCentrCode(AliVEvent* ev)
{
    Short_t centrCode = -1;
    Float_t lPercentile = 0;
    Float_t V0M_Cent = 0, SPD_Cent = 0;
    
    if (fAnalysisType == "AOD"){
        AliAODEvent* aod = (AliAODEvent*)ev;
        AliMultSelection *MultSelection = 0;
        MultSelection = (AliMultSelection * ) aod->FindListObject("MultSelection");
        if(!MultSelection){
            lPercentile = -100;
        }
        else{
            if (fCent == 0)
                lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
            
            if (fCent == 1)
                lPercentile = MultSelection->GetMultiplicityPercentile("CL0");
            
            if (fCent == 2)
                lPercentile = MultSelection->GetMultiplicityPercentile("CL1");
            
            V0M_Cent = MultSelection->GetMultiplicityPercentile("V0M");
            SPD_Cent = MultSelection->GetMultiplicityPercentile("CL1");
            
        }
    }

    //cout << "lPercentile=" << lPercentile << endl;
    
    fCentralityDis->Fill(lPercentile);
    fCentSPDvsV0M->Fill(V0M_Cent, V0M_Cent);
    

  if (fLHC10h) {
	
      if (lPercentile <= 80 && lPercentile > 0 && TMath::Abs(V0M_Cent - SPD_Cent) < 5){
        
      if ((lPercentile > 0) && (lPercentile <= 5.0) && (TMath::Abs(V0M_Cent - SPD_Cent) < 1))
          centrCode = 0;
      else if ((lPercentile > 5.0) && (lPercentile <= 10.0) && (TMath::Abs(V0M_Cent - SPD_Cent) < 1))
          centrCode = 1;
      else if ((lPercentile > 10.0) && (lPercentile <= 20.0))
          centrCode = 2;
      else if ((lPercentile > 20.0) && (lPercentile <= 30.0))
          centrCode = 3;
      else if ((lPercentile > 30.0) && (lPercentile <= 40.0))
          centrCode = 4;
      else if ((lPercentile > 40.0) && (lPercentile <= 50.0))
          centrCode = 5;
      else if ((lPercentile > 50.0) && (lPercentile <= 60.0))
          centrCode = 6;
      else if ((lPercentile > 60.0) && (lPercentile <= 70.0))
          centrCode = 7;
      else if ((lPercentile > 70.0) && (lPercentile <= 80.0))
          centrCode = 8;
      else if ((lPercentile > 80.0) && (lPercentile <= 90.0))
          centrCode = 9;
 
    }
      
  }
	
	
  return centrCode;

}
    
    
//_____________________________________________________________________________
Double_t AliAnalysisTaskChargedFlow::GetWDist(const AliVVertex* v0, const AliVVertex* v1)    {
        
        // calculate sqrt of weighted distance to other vertex
        if (!v0 || !v1) {
            printf("One of vertices is not valid\n");
            return 0;
        }
        static TMatrixDSym vVb(3);
        double dist = -1;
        double dx = v0->GetX()-v1->GetX();
        double dy = v0->GetY()-v1->GetY();
        double dz = v0->GetZ()-v1->GetZ();
        double cov0[6],cov1[6];
        v0->GetCovarianceMatrix(cov0);
        v1->GetCovarianceMatrix(cov1);
        vVb(0,0) = cov0[0]+cov1[0];
        vVb(1,1) = cov0[2]+cov1[2];
        vVb(2,2) = cov0[5]+cov1[5];
        vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
        vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
        vVb.InvertFast();
        if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
        dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
        +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
        return dist>0 ? TMath::Sqrt(dist) : -1;
        
}
    
//_____________________________________________________________________________
Bool_t AliAnalysisTaskChargedFlow::plpMV(const AliVEvent *event)
{
        // check for multi-vertexer pile-up
        const AliAODEvent *aod = (const AliAODEvent*)event;
        const AliESDEvent *esd = (const AliESDEvent*)event;
        //
        const int    kMinPlpContrib = 5;
        const double kMaxPlpChi2 = 5.0;
        const double kMinWDist = 15;
        //
        if (!aod && !esd) {
            printf("Event is neither of AOD nor ESD\n");
            exit(1);
        }
        //
        const AliVVertex* vtPrm = 0;
        const AliVVertex* vtPlp = 0;
        int nPlp = 0;
        //
        if (aod) {
            if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
            vtPrm = aod->GetPrimaryVertex();
            if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
        }
        else {
            if ( !(nPlp=esd->GetNumberOfPileupVerticesTracks())) return kFALSE;
            vtPrm = esd->GetPrimaryVertexTracks();
            if (((AliESDVertex*)vtPrm)->GetStatus()!=1) return kTRUE; // there are pile-up vertices but no primary
        }
        
        //int bcPrim = vtPrm->GetBC();
        //
        for (int ipl=0;ipl<nPlp;ipl++) {
            vtPlp = aod ? (const AliVVertex*)aod->GetPileupVertexTracks(ipl) : (const AliVVertex*)esd->GetPileupVertexTracks(ipl);
            //
            if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
            if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
            //  int bcPlp = vtPlp->GetBC();
            //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
            //
            double wDst = GetWDist(vtPrm,vtPlp);
            if (wDst<kMinWDist) continue;
            //
            return kTRUE; // pile-up: well separated vertices
        }
        //
        return kFALSE;
        //
}


//_____________________________________________________________________________
//void AliAnalysisTaskChargedFlow::ProcessMCTruth(AliAODEvent* aod, Short_t centrV0)




//_____________________________________________________________________________
void AliAnalysisTaskChargedFlow::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
