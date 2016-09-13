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
  fQvec2(0),
  fQvec3(0),
  fQvec4(0),
  fQvec5(0),
  fQvec2Gap00P(0),
  fQvec2Gap04P(0),
  fQvec2Gap08P(0),
  fQvec2Gap10P(0),
  fQvec2Gap00N(0),
  fQvec2Gap04N(0),
  fQvec2Gap08N(0),
  fQvec2Gap10N(0),
  fPOIvec(0),
  fRFPvec(0),
  fLocalEventCounter(0),
  fCent(0),
  fOutputList(0),
  fOutputListQA(0),
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fCentFlag(0),
  fPVtxCutZ(10),
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
  fRefCorTwo2(0),
  fRefCorTwo3(0),
  fRefCorTwo4(0),
  fRefCorTwo5(0),
  fRefCorTwo2Gap00(0),
  fRefCorTwo2Gap04(0),
  fRefCorTwo2Gap08(0),
  fRefCorTwo2Gap10(0),
  fDiffCorTwo(0),
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
  fQvec2(0),
  fQvec3(0),
  fQvec4(0),
  fQvec5(0),
  fQvec2Gap00P(0),
  fQvec2Gap04P(0),
  fQvec2Gap08P(0),
  fQvec2Gap10P(0),
  fQvec2Gap00N(0),
  fQvec2Gap04N(0),
  fQvec2Gap08N(0),
  fQvec2Gap10N(0),
  fPOIvec(0),
  fRFPvec(0),
  fLocalEventCounter(0),
  fCent(0),
  fOutputList(0),
  fOutputListQA(0),
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fCentFlag(0),
  fPVtxCutZ(10),
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
  fRefCorTwo2(0),
  fRefCorTwo3(0),
  fRefCorTwo4(0),
  fRefCorTwo5(0),
  fRefCorTwo2Gap00(0),
  fRefCorTwo2Gap04(0),
  fRefCorTwo2Gap08(0),
  fRefCorTwo2Gap10(0),
  fDiffCorTwo(0),
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
  fRefCorTwo2 = new TProfile("fRefCorTwo2","#LT#LT2#GT#GT (ref. flow) v2",10,-0.5,9.5);
  fRefCorTwo2->Sumw2();
	fOutputList->Add(fRefCorTwo2);
	fRefCorTwo3 = new TProfile("fRefCorTwo3","#LT#LT2#GT#GT (ref. flow) v3",10,-0.5,9.5);
	fRefCorTwo3->Sumw2();
	fOutputList->Add(fRefCorTwo3);
	fRefCorTwo4 = new TProfile("fRefCorTwo4","#LT#LT2#GT#GT (ref. flow) v4",10,-0.5,9.5);
	fRefCorTwo4->Sumw2();
	fOutputList->Add(fRefCorTwo4);
	fRefCorTwo5 = new TProfile("fRefCorTwo5","#LT#LT2#GT#GT (ref. flow) v5",10,-0.5,9.5);
	fRefCorTwo5->Sumw2();
	fOutputList->Add(fRefCorTwo5);

  fRefCorTwo2Gap00 = new TProfile("fRefCorTwo2Gap00","#LT#LT2#GT#GT (ref. flow) v2 Gap00",10,-0.5,9.5);
  fRefCorTwo2Gap00->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap00);
  fRefCorTwo2Gap04 = new TProfile("fRefCorTwo2Gap04","#LT#LT2#GT#GT (ref. flow) v2 Gap04",10,-0.5,9.5);
  fRefCorTwo2Gap04->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap04);
  fRefCorTwo2Gap08 = new TProfile("fRefCorTwo2Gap08","#LT#LT2#GT#GT (ref. flow) v2 Gap08",10,-0.5,9.5);
  fRefCorTwo2Gap08->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap08);
  fRefCorTwo2Gap10 = new TProfile("fRefCorTwo2Gap10","#LT#LT2#GT#GT (ref. flow) v2 Gap10",10,-0.5,9.5);
  fRefCorTwo2Gap10->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap10);
  

  fDiffCorTwo = new TProfile("fDiffCorTwo","#LT#LT2'#GT#GT (diff. flow)",1,0,1);
	fDiffCorTwo->Sumw2();
	fOutputList->Add(fDiffCorTwo);

	// QA output
	Int_t iNEventCounterBins = 4;
	TString sEventCounterLabel[] = {"Input","AOD OK","Selected (me)","Selected (Katka)"};
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
  fEventCounter->Fill(2); // event selected

  // only events passing selection criteria defined @ IsEventSelected()

  const Int_t iTracks(fAOD->GetNumberOfTracks());           
  fEventMult->Fill(iTracks);
  
  Int_t iNumTracksSelected = 0;
  Int_t iNumGap00P = 0;
  Int_t iNumGap00N = 0;
  Int_t iNumGap04P = 0;
  Int_t iNumGap04N = 0;
  Int_t iNumGap08P = 0;
  Int_t iNumGap08N = 0;
  Int_t iNumGap10P = 0;
  Int_t iNumGap10N = 0;
  
  
	// estimating flow vectors for 2-part correlations
  //fQvec = TComplex(0,0,kFALSE); 
  fQvec2 = TComplex(0,0,kFALSE); 
  fQvec3 = TComplex(0,0,kFALSE); 
  fQvec4 = TComplex(0,0,kFALSE); 
  fQvec5 = TComplex(0,0,kFALSE); 

  fQvec2Gap00P = TComplex(0,0,kFALSE);
  fQvec2Gap04P = TComplex(0,0,kFALSE);
  fQvec2Gap08P = TComplex(0,0,kFALSE);
  fQvec2Gap10P = TComplex(0,0,kFALSE);
  fQvec2Gap00N = TComplex(0,0,kFALSE);
  fQvec2Gap04N = TComplex(0,0,kFALSE);
  fQvec2Gap08N = TComplex(0,0,kFALSE);
  fQvec2Gap10N = TComplex(0,0,kFALSE);

  fPOIvec = TComplex(0,0,kFALSE);
  fRFPvec = TComplex(0,0,kFALSE);

  Int_t iCounterPOI = 0;

  // loop over all tracks
  for(Int_t i(0); i < iTracks; i++) 
  {                 
    fTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!fTrack) continue;
    
    if(!IsTrackSelected(fTrack)) continue;
    // only selected tracks
    iNumTracksSelected++;

    fTracksPt->Fill(fTrack->Pt());                     
    fTracksEta->Fill(fTrack->Eta());   
    fTracksPhi->Fill(fTrack->Phi());   

    fQvec2 += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
    fQvec3 += TComplex(TMath::Cos(3*(fTrack->Phi())),TMath::Sin(3*(fTrack->Phi())),kFALSE);
    fQvec4 += TComplex(TMath::Cos(4*(fTrack->Phi())),TMath::Sin(4*(fTrack->Phi())),kFALSE);
	  fQvec5 += TComplex(TMath::Cos(5*(fTrack->Phi())),TMath::Sin(5*(fTrack->Phi())),kFALSE);

    // eta gap
    Double_t dTrackEta = fTrack->Eta();

    if(dTrackEta > 0.)
    {
      fQvec2Gap00P += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap00P++;
    }

    if(dTrackEta < 0.)
    {
      fQvec2Gap00N += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap00N++;
    }
    
    if(dTrackEta > 0.2)
    {
      fQvec2Gap04P += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap04P++;
    }

    if(dTrackEta < -0.2)
    {
      fQvec2Gap04N += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap04N++;
    }
    
    if(dTrackEta > 0.4)
    {
      fQvec2Gap08P += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap08P++;
    }

    if(dTrackEta < -0.4)
    {
      fQvec2Gap08N += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap08N++;
    }
    
    if(dTrackEta > 0.5)
    {
      fQvec2Gap10P += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap10P++;
    }

    if(dTrackEta < -0.5)
    {
      fQvec2Gap10N += TComplex(TMath::Cos(2*(fTrack->Phi())),TMath::Sin(2*(fTrack->Phi())),kFALSE);
      iNumGap10N++;
    }    
  }   // end of loop over all tracks

  if(iNumTracksSelected == 0) return; // 0 tracks selected in this event
  fLocalEventCounter++;

  fMultTracksSelected->Fill(iNumTracksSelected);
  
  // Reference Flow
  Double_t dAmp = 0;
  Double_t dVal = 0;
  Double_t dWeight = 0;

  dWeight = iNumTracksSelected*(iNumTracksSelected-1);

  dAmp = (fQvec2*(TComplex::Conjugate(fQvec2))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCent != -1) )
    fRefCorTwo2->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dAmp = (fQvec3*(TComplex::Conjugate(fQvec3))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCent != -1) )
    fRefCorTwo3->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile

	dAmp = (fQvec4*(TComplex::Conjugate(fQvec4))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCent != -1) )
    fRefCorTwo4->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dAmp = (fQvec5*(TComplex::Conjugate(fQvec5))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCent != -1) )
    fRefCorTwo5->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile

  // ref flow with eta gap

  dWeight = iNumGap00P*iNumGap00N;
  dAmp = (fQvec2Gap00P*(TComplex::Conjugate(fQvec2Gap00N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCent != -1) )
    fRefCorTwo2Gap00->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile

  dWeight = iNumGap04P*iNumGap04N;
  dAmp = (fQvec2Gap04P*(TComplex::Conjugate(fQvec2Gap04N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCent != -1) )
    fRefCorTwo2Gap04->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dWeight = iNumGap08P*iNumGap08N;
  dAmp = (fQvec2Gap08P*(TComplex::Conjugate(fQvec2Gap08N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCent != -1) )
    fRefCorTwo2Gap08->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dWeight = iNumGap10P*iNumGap10N;
  dAmp = (fQvec2Gap10P*(TComplex::Conjugate(fQvec2Gap10N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCent != -1) )
    fRefCorTwo2Gap10->Fill(fCent, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
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
		fCent = GetCentrCode(fAOD);
    if (fCent < 0)
      return kFALSE;
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
  
  if (fAODAnalysis)
  {
    AliAODEvent* aod = (AliAODEvent*) ev;
    AliMultSelection* MultSelection = 0;
    MultSelection = (AliMultSelection*) aod->FindListObject("MultSelection");
    
    if(!MultSelection)
    {
      lPercentile = -100;
    }
    else
    {
      if(fCentFlag == 0)
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
      
      if(fCentFlag == 1)
        lPercentile = MultSelection->GetMultiplicityPercentile("CL0");
      
      if(fCentFlag == 2)
        lPercentile = MultSelection->GetMultiplicityPercentile("CL1");
      
      V0M_Cent = MultSelection->GetMultiplicityPercentile("V0M");
      SPD_Cent = MultSelection->GetMultiplicityPercentile("CL1");   
    }
  }

  //cout << "lPercentile=" << lPercentile << endl;
  
  fCentralityDis->Fill(lPercentile);
  fCentSPDvsV0M->Fill(V0M_Cent, SPD_Cent);

  if (fLHC10h) {
      if (lPercentile <= 80 && lPercentile > 0 && TMath::Abs(V0M_Cent - SPD_Cent) < 5)
      {  
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

Bool_t AliAnalysisTaskFlowPID::IsEventSelectedKatarina(const AliAODEvent* event)
{
  Float_t zvtx = GetVertex(fAOD);
  
  if(zvtx< -990)
  {
    //fVtx->Fill(0);
    return kFALSE;
  }
  else 
  {
    //fVtx->Fill(1);
    //fVtxBeforeCuts->Fill(zvtx);
    if (TMath::Abs(zvtx) < fPVtxCutZ) 
    {
      //fMultCorBeforeCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));

    
      Short_t isPileup = fAOD->IsPileupFromSPD(3);
      if (isPileup != 0)
        return kFALSE;
                  
      //if (fAOD->GetHeader()->GetRefMultiplicityComb08() < 0)
      //    return;
      
      //New cut from Ruben for pileup hybrid
      if (plpMV(fAOD))
          return kFALSE;
    
                
                
      //fMultCorAfterCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));
      Short_t cenAOD  = GetCentrCode(fAOD);
      //Short_t cenAOD  = 1;
      //cout << cenAOD << endl;
            
      if (cenAOD >= 0){
        //fVtxAfterCuts->Fill(zvtx);
        return kTRUE;
      }
    }
  }
  return kFALSE;
}

Float_t AliAnalysisTaskFlowPID::GetVertex(AliVEvent* ev) const
{

  Float_t vtxz = -999.;

  if (1){

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

Bool_t AliAnalysisTaskFlowPID::plpMV(const AliVEvent *event)
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
Double_t AliAnalysisTaskFlowPID::GetWDist(const AliVVertex* v0, const AliVVertex* v1)    {
        
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
    