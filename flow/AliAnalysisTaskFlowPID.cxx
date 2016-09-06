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

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskFlowPID.h"

#include "AliLog.h" 

class AliAnalysisTaskFlowPID;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskFlowPID) // classimp: necessary for root

AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID() : AliAnalysisTaskSE(), 
  fAOD(0),
  fOutputList(0),
  fOutputListQA(0),
  fHistPt(0),
  fQAPVz(0),
  fQANumTracks(0),
  fEventCounter(0),
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fPVtxCutZ(0),
  fCentEdgeLow(0),
  fCentEdgeUp(0)
{
	  // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),
  fOutputList(0),
  fOutputListQA(0),
  fHistPt(0),
  fEventCounter(0),
  fQAPVz(0),
  fQANumTracks(0),
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fPVtxCutZ(0),
  fCentEdgeLow(0),
  fCentEdgeUp(0)
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
	Int_t iNEventCounterBins = 3;
	TString sEventCounterLabel[] = {"Input","AOD OK","Selected"};
	fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
	for(Int_t i = 0; i < iNEventCounterBins; i++)
		fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
	fOutputList->Add(fEventCounter);

	fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);    
	fOutputList->Add(fHistPt);          
	
	// QA output
	fQAPVz = new TH1D("fQAPVz","QA: PV #it{z}",100,-50,50);
	fOutputListQA->Add(fQAPVz);
	fQANumTracks = new TH1D("fQANumTracks","QA: Number of AOD tracks",100,0,1000);
	fOutputListQA->Add(fQANumTracks);

	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
	PostData(2, fOutputListQA);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserExec(Option_t *)
{
  // this function is called once for each event

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

  Int_t iTracks(fAOD->GetNumberOfTracks());           
  for(Int_t i(0); i < iTracks; i++) 
  {                 
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!track) continue;                            
      fHistPt->Fill(track->Pt());                     
  }                                                   
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

	const Double_t dNumAODTracks = event->GetNumberOfTracks();
	fQANumTracks->Fill(dNumAODTracks);


	return; 
}

