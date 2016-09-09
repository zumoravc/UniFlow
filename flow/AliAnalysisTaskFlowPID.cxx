/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskFlowPID
 *
 * analysis task for flow study of PID particles
 * Note: So far implemented only for AOD analysis!
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "TChain.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TList.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliVTrack.h"
#include "TComplex.h"
#include "AliAnalysisTaskFlowPID.h"

#include "AliLog.h" 

class AliAnalysisTaskFlowPID;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskFlowPID) // classimp: necessary for root

AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID() : AliAnalysisTaskSE(), 
  fAOD(0),
  fTrack(0),
  fQvec(0),
  fQvecTest(0),
  fPOIvec(0),
  fRFPvec(0),
  fArrTracksSelected("AliAODTrack",5000),
  fLocalEventCounter(0),
  fCent(0),
  fOutputList(0),
  fOutputListQA(0),
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fPVtxCutZ(10),
  fCentEdgeLow(0),
  fCentEdgeUp(0),
  fTrackEtaMax(0.8),
  fTrackPtMax(10),
  fTrackPtMin(0.2),
  fNumTPCclsMin(70),
  fTrackFilterBit(128),
  fEventCounter(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentSPDvsV0M(0),
  fMultTracksSelected(0),
  fTracksPt(0),
  fTracksEta(0),
  fTracksPhi(0),
  fRefCor2(0),
  fRefCor2Test(0),
  fDiffCor2(0),
  fQAPVz(0),
  fQANumTracks(0),
  fQATrackPt(0),
  fQATrackEta(0),
  fQATrackPhi(0),
  fQATrackFilterMap(0)
{
	  // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),
  fTrack(0),
  fQvec(0),
  fQvecTest(0),
  fPOIvec(0),
  fRFPvec(0),
  fArrTracksSelected("AliAODTrack",5000),
  fLocalEventCounter(0),
  fCent(0),
  fOutputList(0),
  fOutputListQA(0),
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fPVtxCutZ(10),
  fCentEdgeLow(0),
  fCentEdgeUp(0),
  fTrackEtaMax(0.8),
  fTrackPtMax(10),
  fTrackPtMin(0.2),
  fNumTPCclsMin(70),
  fTrackFilterBit(128),
  fEventCounter(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentSPDvsV0M(0),
  fMultTracksSelected(0),
  fTracksPt(0),
  fTracksEta(0),
  fTracksPhi(0),
  fRefCor2(0),
  fRefCor2Test(0),
  fDiffCor2(0),
  fQAPVz(0),
  fQANumTracks(0),
  fQATrackPt(0),
  fQATrackEta(0),
  fQATrackPhi(0),
  fQATrackFilterMap(0)
{
  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                      // this chain is created by the analysis manager, so no need to worry about it, 
                                      // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                      // you can add more output objects by calling DefineOutput(2, classname::Class())
                                      // if you add more output objects, make sure to call PostData for all of them, and to
                                      // make changes to your AddTask macro!
  DefineOutput(2, TList::Class());	
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::~AliAnalysisTaskFlowPID()
{
  // destructor
  if(fOutputList) 
  {
      delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserCreateOutputObjects()
{
	// create output objects
	// this function is called ONCE at the start of your analysis (RUNTIME)
	// here you create the histograms that you want to use 
	// the histograms are in this case added to a TList, this list is in the end saved to an output file

	fOutputList = new TList();          // this is a list which will contain all of your histograms at the end of the analysis, the contents of this list are written to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested (dont worry about this now)

	fOutputListQA = new TList();
	fOutputListQA->SetOwner(kTRUE);

	// main output
	fEventMult = new TH1D("fEventMult","Event multiplicity (selected)",100,0,5000);
	fOutputList->Add(fEventMult);
	fCentralityDis = new TH1D("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
  fOutputList->Add(fCentralityDis);
	fCentSPDvsV0M = new TH2D("fCentSPDvsV0M", "V0M-cent vs SPD-cent; V0M; SPD-cent", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fCentSPDvsV0M);
	fMultTracksSelected = new TH1D("fMultTracksSelected","Multiplicity of selected tracks",100,0,5000);
	fOutputList->Add(fMultTracksSelected);
	fTracksPt = new TH1D("fTracksPt", "Tracks #it{p}_{T} (selected)", 100, 0, 10);    
	fOutputList->Add(fTracksPt);          
	fTracksEta = new TH1D("fTracksEta", "Tracks #it{#eta} (selected)", 300, -1.5, 1.5);    
	fOutputList->Add(fTracksEta);          
	fTracksPhi = new TH1D("fTracksPhi", "Tracks #it{#varphi} (selected)", 360, 0., TMath::TwoPi());    
	fOutputList->Add(fTracksPhi);          
	fRefCor2 = new TProfile("fRefCor2","#LT#LT2#GT#GT (ref. flow)",10,-0.5,9.5);
	fRefCor2->Sumw2();
	fOutputList->Add(fRefCor2);
	fRefCor2Test = new TProfile("fRefCor2Test","#LT#LT2#GT#GT (ref. flow)",10,-0.5,9.5);
	fRefCor2Test->Sumw2();
	fOutputList->Add(fRefCor2Test);
	fDiffCor2 = new TProfile("fDiffCor2","#LT#LT2'#GT#GT (diff. flow)",1,0,1);
	fDiffCor2->Sumw2();
	fOutputList->Add(fDiffCor2);

	// QA output
	Int_t iNEventCounterBins = 3;
	TString sEventCounterLabel[] = {"Input","AOD OK","Selected"};
	fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
	for(Int_t i = 0; i < iNEventCounterBins; i++)
		fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
	fOutputListQA->Add(fEventCounter);
	fQAPVz = new TH1D("fQAPVz","QA: PV #it{z}",100,-50,50);
	fOutputListQA->Add(fQAPVz);
	fQANumTracks = new TH1D("fQANumTracks","QA: Number of AOD tracks",1000,0,10000);
	fOutputListQA->Add(fQANumTracks);
	fQATrackPt = new TH1D("fQATrackPt","QA: Track #it{p}_{T} (all)",100,0,10);
	fOutputListQA->Add(fQATrackPt);
	fQATrackEta = new TH1D("fQATrackEta","QA: Track #it{#eta} (all)",300,-1.5,1.5);
	fOutputListQA->Add(fQATrackEta);
	fQATrackPhi = new TH1D("fQATrackPhi","QA: Track #it{#varphi} (all)",300,0,TMath::TwoPi());
	fOutputListQA->Add(fQATrackPhi);
	fQATrackFilterMap = new TH1D("fQATrackFilterMap","QA: Tracks filter map (all)",1000,0,1000);
	fOutputListQA->Add(fQATrackFilterMap);

	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
	PostData(2, fOutputListQA);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserExec(Option_t *)
{
  // this function is called once for each event
	//if(fLocalEventCounter > 0) return; // just for debugging purposes: stops after first selected event

	fEventCounter->Fill(0); // input event

  if( !fAODAnalysis || dynamic_cast<AliESDEvent*>(InputEvent()) ) // ESD analysis
  {
  	::Warning("UserExec","ESD event: not implemented. Terminating!");
  	return;	
  }

  if(fAODAnalysis) // AOD analysis 
  {
  	fAOD = dynamic_cast<AliAODEvent*>(InputEvent()); 
	  if(!fAOD) return;
  }

  fEventCounter->Fill(1); // event AOD ok
  
  // basic event QA
  EventQA(fAOD);

  // event selection
  if(!IsEventSelected(fAOD))
	{
		return;
	}
	fCent = GetCentrCode(fAOD);

  fEventCounter->Fill(2); // event selected
	// only events passing selection criteria defined @ IsEventSelected()

  fArrTracksSelected.Clear("C");
  //printf("Array:pre-Enties: %d\n",fArrTracksSelected.GetEntries());

	const Int_t iTracks(fAOD->GetNumberOfTracks());           
  fEventMult->Fill(iTracks);
  
  Int_t iNumTracksSelected = 0;
  // loop over all tracks
  for(Int_t i(0); i < iTracks; i++) 
  {                 
      fTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!fTrack) continue;
      
      if(!IsTrackSelected(fTrack)) continue;
      // only selected tracks

      //cloning the track to TClonesArray of selected tracks
      new(fArrTracksSelected[iNumTracksSelected]) AliAODTrack(*fTrack);
      iNumTracksSelected++;
      
      fTracksPt->Fill(fTrack->Pt());                     
      fTracksEta->Fill(fTrack->Eta());   
      fTracksPhi->Fill(fTrack->Phi());   
  }  	// end of loop over all tracks

  // input tracks fileter and stored in fArrTracksSelected array
  if(iNumTracksSelected == 0) return; // 0 tracks selected in this event
	fLocalEventCounter++;

  fMultTracksSelected->Fill(iNumTracksSelected);
	
	// estimating flow vectors for 2-part correlations
  fQvec = TComplex(0,0,kFALSE); 
  fQvecTest = TComplex(0,0,kFALSE); 
  fPOIvec = TComplex(0,0,kFALSE);
  fRFPvec = TComplex(0,0,kFALSE);

  Int_t iCounterPOI = 0;

  AliAODTrack* track1 = NULL;
  AliAODTrack* track2 = NULL;
  Double_t dPhiTrack1 = 0;
  Double_t dPhiTrack2 = 0;

  const Int_t iHarmonic = 2; // harmonic component (=n)

  // loop over selected (filtered) tracks in given event
  for(Int_t i = 0; i < iNumTracksSelected; i++)
  {
  	track1 = static_cast<AliAODTrack*>(fArrTracksSelected.At(i));
  	dPhiTrack1 = track1->Phi();
		fQvecTest += TComplex(TMath::Cos(iHarmonic*(dPhiTrack1)),TMath::Sin(iHarmonic*(dPhiTrack1)),kFALSE);
  
/*
  	if( (track1->Pt() > 1) && (track1->Pt() < 2) )
  	{
  		iCounterPOI++;
  		fPOIvec += TComplex(TMath::Cos(iHarmonic*dPhiTrack1), TMath::Sin(iHarmonic*dPhiTrack1));
  	}
*/
  	for(Int_t j = 0; j < iNumTracksSelected; j++)
  	{
  		track2 = static_cast<AliAODTrack*>(fArrTracksSelected.At(j));
  		dPhiTrack2 = track2->Phi();

  		fQvec += TComplex(TMath::Cos(iHarmonic*(dPhiTrack1-dPhiTrack2)),TMath::Sin(iHarmonic*(dPhiTrack1-dPhiTrack2)),kFALSE);
  		


  		//fRFPvec += TComplex(TMath::Cos(iHarmonic*dPhiTrack2), TMath::Sin(iHarmonic*dPhiTrack2));
  	}
  	//printf("Re(Q): %f // Phi1: %f // Phi 2:%f \n",fQvec.Re(), dPhiTrack1, dPhiTrack2 );
  }

  // Reference Flow
  Double_t dWeight = iNumTracksSelected*(iNumTracksSelected-1);
  Double_t dNom = (fQvec.Re() - iNumTracksSelected)/dWeight;
  fRefCor2->Fill(fCent, dNom, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
	//printf("Ref: %f / weight: %f \n", dNom, dWeight);


  //printf("Orig: %f(%f) //",dNom, dWeight);

  dNom = (fQvecTest*(TComplex::Conjugate(fQvecTest))).Re();
  dNom = (dNom - iNumTracksSelected) / dWeight;
  fRefCor2Test->Fill(fCent, dNom, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  //printf("Test: %f(%f)\n",dNom, dWeight);


/*
  // differential flow 
	Double_t dotProduct = (fPOIvec*(TComplex::Conjugate(fRFPvec))).Re(); // Re() of Tcom*Tcom <=> dot product of TComplex
	dWeight = iCounterPOI*(iNumTracksSelected-1);
	dNom = (dotProduct - iCounterPOI) / dWeight;
	
	printf("Diff: %f / product %f / weight: %f \n", dNom,dotProduct, dWeight);
	fDiffCor2->Fill(0.5,dNom,dWeight);
*/


  PostData(1, fOutputList);	// stream the results the analysis of this event to the output manager which will take care of writing it to a file
  PostData(2, fOutputListQA);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsEventSelected(const AliAODEvent* event)
{
  // event selection criteria

  // Pileup rejection
  if( event->IsPileupFromSPD(3) ) //min contributors ???  	
	{
		if(fDebug) ::Info("IsEventSelected","Event rejected: SPD pile up");
		return kFALSE;
	}

  	// consider including plpMV pileup code from Katarina?
 
  // Primary vertex criteria
  const AliAODVertex* aodVtx = event->GetPrimaryVertex();
  if(!aodVtx)
  {
  	if(fDebug > 0) 
  		::Info("IsEventSelected","Event rejected: Primary vertex not found");
  	return kFALSE;
  }

 Int_t iNContrib = aodVtx->GetNContributors();
  if(iNContrib <= 0) 
  {	
  	if(fDebug > 0)
  	{
  		::Info("IsEventSelected","Event rejected: Low number of PV contributors");
  	} 
  	return kFALSE;
  }

  TString sVtxTitle = aodVtx->GetTitle();
  if( !sVtxTitle.Contains("VertexerTracks") )
  {
		if(fDebug > 1)
			::Info("IsEventSelected","Event rejected: Title does not contain 'VertexerTracks'");
  	return kFALSE;
  }

  const AliAODVertex* spdVtx = event->GetPrimaryVertexSPD();
  if(!spdVtx)
  {
  	if(fDebug) 
  		::Info("IsEventSelected","Event rejected: SPD vertex not found");
  	return kFALSE;
  }

  Int_t iNContribSPD = spdVtx->GetNContributors();
  if(iNContribSPD <= 0)
  {
  	if(fDebug) 
  		::Info("IsEventSelected","Event rejected: Low number of SPD contributors");
  	return kFALSE;
  }

  const Double_t aodVtxZ = aodVtx->GetZ();
  const Double_t spdVtxZ = spdVtx->GetZ();
  if(TMath::Abs(aodVtxZ - spdVtxZ) > 0.5)
  {
  	if(fDebug) 
  		::Info("IsEventSelected","Event rejected: High z-distance diference between AOD and SPD vertex");
  	return kFALSE;
  }	

  if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
	{
		if(fDebug) 
			::Info("IsEventSelected","Event rejected: PV z-distance cut");
  	return kFALSE;
	}

	// centrality rejection
	
	if(fPbPb)
	{
		// not implemented yet
		// centrality solving issue?
	}

 	return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsTrackSelected(const AliAODTrack* track)
{
	// track selection procedure

	if(!track)
	{
		return kFALSE;
	}

	if( !track->TestFilterBit(fTrackFilterBit) )
	{
		return kFALSE;	
	}

	if(track->GetTPCNcls() < fNumTPCclsMin)
	{
		return kFALSE;
	}

	if(TMath::Abs(track->Eta()) > fTrackEtaMax)
	{
		return kFALSE;
	}

	if( (track->Pt() > fTrackPtMax) || (track->Pt() < fTrackPtMin) )
	{
		return kFALSE;
	}

	return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::EventQA(const AliAODEvent* event)
{
  // event QA procedure
	const AliAODVertex* aodVtx = event->GetPrimaryVertex();
	const Double_t dVtxZ = aodVtx->GetZ();
	//const Double_t dNumContrinutors = aodVtx->GetNContributors();

	fQAPVz->Fill(dVtxZ);


	// tracks QA procedure
	const Int_t iNumAODTracks = event->GetNumberOfTracks();
	fQANumTracks->Fill(iNumAODTracks);


	const AliAODTrack* track = NULL;
	for(Int_t i = 0; i < iNumAODTracks; i++)
	{
		track = static_cast<AliAODTrack*>(event->GetTrack(i));
		if(!track) continue;

		fQATrackPt->Fill(track->Pt());
		fQATrackEta->Fill(track->Eta());
		fQATrackPhi->Fill(track->Phi());
		fQATrackFilterMap->Fill(track->GetFilterMap());
	}

	return; 
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowPID::GetCentrCode(AliVEvent* ev)
{
    Short_t centrCode = -1;
    Float_t lPercentile = 0;
    Float_t V0M_Cent = 0, SPD_Cent = 0;
    
    if (fAODAnalysis){
        AliAODEvent* aod = (AliAODEvent*)ev;
        AliMultSelection* MultSelection = 0;
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
