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

class AliAnalysisTaskFlowPID;    

ClassImp(AliAnalysisTaskFlowPID) // classimp: necessary for root

Double_t AliAnalysisTaskFlowPID::fPtBinEdges[] = {0.,1.,2.,3.,4.};
Double_t AliAnalysisTaskFlowPID::fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID() : AliAnalysisTaskSE(), 
  fAOD(0),
  fTrack(0),
  fTrackPt(0),
  fTrackEta(0),
  fTrackPhi(0),
  fCentBinIndex(0),
  fCentPercentile(0),
  fPtBinIndex(0),
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
  fDiffFlow(kTRUE),
  fEventCounter(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentSPDvsV0M(0),
  fMultTracksSelected(0),
  fTracksPtCent(0),
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
  fTrackPt(0),
  fTrackEta(0),
  fTrackPhi(0),
  fCentBinIndex(0),
  fCentPercentile(0),
  fPtBinIndex(0),
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
  
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fCentFlag(0),
  fPVtxCutZ(10),
  fTrackEtaMax(0.8),
  fTrackPtMax(5),
  fTrackPtMin(0.2),
  fNumTPCclsMin(70),
  fTrackFilterBit(128),
  fDiffFlow(kTRUE),

  fOutputList(0),
  fOutputListQA(0),
  fEventCounter(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentSPDvsV0M(0),
  fMultTracksSelected(0),
  fTracksPtCent(0),
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
  fQAPVz(0),
  fQANumTracks(0),
  fQATrackPt(0),
  fQATrackEta(0),
  fQATrackPhi(0),
  fQATrackFilterMap(0)
{
  // constructor

  for(Int_t i = 0; i < fNumPtBins; i++)
  {
     fPvec2[i] = 0;
     fPvec2Gap00P[i] = 0;
     fPvec2Gap04P[i] = 0;
     fPvec2Gap08P[i] = 0;
     fPvec2Gap10P[i] = 0;
     fPvec3[i] = 0;
  }

  for(Int_t i(0); i < fNumCentBins; i++)
  {
    fDiffCorTwo2[i] = 0;
    fDiffCorTwo2Gap00[i] = 0;
    fDiffCorTwo2Gap04[i] = 0;
    fDiffCorTwo2Gap08[i] = 0;
    fDiffCorTwo2Gap10[i] = 0;
    fDiffCorTwo3[i] = 0;
  }
  

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
  fEventMult = new TH1D("fEventMult","Track multiplicity (all tracks in selected events); tracks;",100,0,10000);
  fOutputList->Add(fEventMult);
  fCentralityDis = new TH1D("fCentralityDis", "centrality distribution; centrality;", fNumCentBins,fCentBinEdges);
  fOutputList->Add(fCentralityDis);
  fMultTracksSelected = new TH1D("fMultTracksSelected","Track multiplicity (selected tracks in selected events); tracks;",100,0,5000);
  fOutputList->Add(fMultTracksSelected);
  fTracksPtCent = new TH2D("fTracksPtCent", "Tracks #it{p}_{T} vs. centrality (selected); #it{p}^{track}_{T} (GeV/#it{c}); centrality;", fNumPtBins,fPtBinEdges,fNumCentBins,fCentBinEdges);    
  fOutputList->Add(fTracksPtCent);          
  fTracksPt = new TH1D("fTracksPt", "Tracks #it{p}_{T} (selected); #it{p}^{track}_{T} (GeV/#it{c});", 100, 0, 10);    
  fOutputList->Add(fTracksPt);          
  fTracksEta = new TH1D("fTracksEta", "Tracks #it{#eta} (selected); #it{#eta}^{track};", 300, -1.5, 1.5);    
  fOutputList->Add(fTracksEta);          
  fTracksPhi = new TH1D("fTracksPhi", "Tracks #it{#varphi} (selected); #it{#varphi}^{track};", 360, 0., TMath::TwoPi());    
  fOutputList->Add(fTracksPhi);          
  fRefCorTwo2 = new TProfile("fRefCorTwo2","#LT#LT2#GT#GT_{2} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2->Sumw2();
  fOutputList->Add(fRefCorTwo2);
  fRefCorTwo3 = new TProfile("fRefCorTwo3","#LT#LT2#GT#GT_{3} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo3->Sumw2();
  fOutputList->Add(fRefCorTwo3);
  fRefCorTwo4 = new TProfile("fRefCorTwo4","#LT#LT2#GT#GT_{4} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo4->Sumw2();
  fOutputList->Add(fRefCorTwo4);
  fRefCorTwo5 = new TProfile("fRefCorTwo5","#LT#LT2#GT#GT_{5} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo5->Sumw2();
  fOutputList->Add(fRefCorTwo5);

  fRefCorTwo2Gap00 = new TProfile("fRefCorTwo2_Gap00","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap00->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap00);
  fRefCorTwo2Gap04 = new TProfile("fRefCorTwo2_Gap04","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.2} (ref. flow); centrality",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap04->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap04);
  fRefCorTwo2Gap08 = new TProfile("fRefCorTwo2_Gap08","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.4} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap08->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap08);
  fRefCorTwo2Gap10 = new TProfile("fRefCorTwo2_Gap10","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.5} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap10->Sumw2();
  fOutputList->Add(fRefCorTwo2Gap10);
  
  if(fDiffFlow) // do differential flow switch
  {
    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2[i] = new TProfile(Form("fDiffCorTwo2_Cent%d",i),Form("#LT#LT2'#GT#GT_{2} Cent %g-%g%% (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2[i]->Sumw2();
      fOutputList->Add(fDiffCorTwo2[i]);
    }
    
    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap00[i] = new TProfile(Form("fDiffCorTwo2_Gap00_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0} Cent %g-%g%% |#it{#eta}^{POI}|>0 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap00[i]->Sumw2();
      fOutputList->Add(fDiffCorTwo2Gap00[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap04[i] = new TProfile(Form("fDiffCorTwo2_Gap04_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0.2} Cent %g-%g%% #it{#eta}^{POI}>0.2 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap04[i]->Sumw2();
      fOutputList->Add(fDiffCorTwo2Gap04[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap08[i] = new TProfile(Form("fDiffCorTwo2_Gap08_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0.4} Cent %g-%g%% #it{#eta}^{POI}>0.4 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap08[i]->Sumw2();
      fOutputList->Add(fDiffCorTwo2Gap08[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap10[i] = new TProfile(Form("fDiffCorTwo2_Gap10_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0.5} Cent %g-%g%% #it{#eta}^{POI}>0.5 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap10[i]->Sumw2();
      fOutputList->Add(fDiffCorTwo2Gap10[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo3[i] = new TProfile(Form("fDiffCorTwo3_Cent%d",i),Form("#LT#LT2'#GT#GT_{3} Cent %g-%g%% (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo3[i]->Sumw2();
      fOutputList->Add(fDiffCorTwo3[i]);
    }
    
  }
  

  // QA output
  Int_t iNEventCounterBins = 8;
  TString sEventCounterLabel[] = {"Input","AOD OK","Pile-up OK","PV OK","SPD Vtx OK","PV #it{z} OK","Centrality OK","At least 2 selected tracks"};
  fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
  for(Int_t i = 0; i < iNEventCounterBins; i++)
    fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
  fOutputListQA->Add(fEventCounter);
  fQAPVz = new TH1D("fQAPVz","QA: PV #it{z}; #it{z} (cm);",100,-50,50);
  fOutputListQA->Add(fQAPVz);
  fCentSPDvsV0M = new TH2D("fCentSPDvsV0M", "V0M-cent vs SPD-cent; V0M; SPD-cent", 100, 0, 100, 100, 0, 100);
  fOutputListQA->Add(fCentSPDvsV0M);
  fQANumTracks = new TH1D("fQANumTracks","QA: Number of AOD tracks; tracks;",100,0,10000);
  fOutputListQA->Add(fQANumTracks);
  fQATrackPt = new TH1D("fQATrackPt","QA: Track #it{p}_{T} (all); #it{p}^{track}_{T} (GeV/#it{c});",100,0,10);
  fOutputListQA->Add(fQATrackPt);
  fQATrackEta = new TH1D("fQATrackEta","QA: Track #it{#eta} (all); #it{#eta}^{track};",300,-1.5,1.5);
  fOutputListQA->Add(fQATrackEta);
  fQATrackPhi = new TH1D("fQATrackPhi","QA: Track #it{#varphi} (all); #it{#varphi}^{track};",300,0,TMath::TwoPi());
  fOutputListQA->Add(fQATrackPhi);
	fQATrackFilterMap = new TH1D("fQATrackFilterMap","QA: Tracks filter map (all); filter bit;",1000,0,1000);
	fOutputListQA->Add(fQATrackFilterMap);

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
  // only events passing selection criteria defined @ IsEventSelected()
  
  const Int_t iTracks(fAOD->GetNumberOfTracks());           
  fEventMult->Fill(iTracks);
  
  // track counters in different regions
  Int_t iNumTracksSelected = 0;
  Int_t iNumGap00P = 0;
  Int_t iNumGap00N = 0;
  Int_t iNumGap04P = 0;
  Int_t iNumGap04N = 0;
  Int_t iNumGap08P = 0;
  Int_t iNumGap08N = 0;
  Int_t iNumGap10P = 0;
  Int_t iNumGap10N = 0;
  
	// estimating flow vectors for <2> with no eta gap (diff harmonics)
  fQvec2 = TComplex(0,0,kFALSE); 
  fQvec3 = TComplex(0,0,kFALSE); 
  fQvec4 = TComplex(0,0,kFALSE); 
  fQvec5 = TComplex(0,0,kFALSE); 

  // <2> eta gap and n = 2
  fQvec2Gap00P = TComplex(0,0,kFALSE);
  fQvec2Gap04P = TComplex(0,0,kFALSE);
  fQvec2Gap08P = TComplex(0,0,kFALSE);
  fQvec2Gap10P = TComplex(0,0,kFALSE);
  fQvec2Gap00N = TComplex(0,0,kFALSE);
  fQvec2Gap04N = TComplex(0,0,kFALSE);
  fQvec2Gap08N = TComplex(0,0,kFALSE);
  fQvec2Gap10N = TComplex(0,0,kFALSE);

  // diff flow
  Int_t iNumP2[fNumPtBins] = {0};
  Int_t iNumP2Gap00P[fNumPtBins] = {0};
  Int_t iNumP2Gap04P[fNumPtBins] = {0};
  Int_t iNumP2Gap08P[fNumPtBins] = {0};
  Int_t iNumP2Gap10P[fNumPtBins] = {0};

  for(Int_t i(0); i < fNumPtBins; i++)
  {
    fPvec2[i] = TComplex(0,0,kFALSE);
    fPvec2Gap00P[i] = TComplex(0,0,kFALSE);
    fPvec2Gap04P[i] = TComplex(0,0,kFALSE);
    fPvec2Gap08P[i] = TComplex(0,0,kFALSE);
    fPvec2Gap10P[i] = TComplex(0,0,kFALSE);
    fPvec3[i] = TComplex(0,0,kFALSE);
  }

  // loop over all tracks
  for(Int_t i(0); i < iTracks; i++) 
  {                 
    fTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!fTrack) continue;
    
    if(!IsTrackSelected(fTrack)) continue;
    
    // only selected tracks
    iNumTracksSelected++;

    fTrackEta = fTrack->Eta();
    fTrackPhi = fTrack->Phi();
    fTrackPt = fTrack->Pt();
    fPtBinIndex = GetPtBinIndex(fTrackPt);

    fTracksPtCent->Fill(fTrackPt,fCentPercentile);
    fTracksPt->Fill(fTrackPt);                     
    fTracksEta->Fill(fTrackEta);   
    fTracksPhi->Fill(fTrackPhi);   

    fQvec2 += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
    fQvec3 += TComplex(TMath::Cos(3*fTrackPhi),TMath::Sin(3*fTrackPhi),kFALSE);
    fQvec4 += TComplex(TMath::Cos(4*fTrackPhi),TMath::Sin(4*fTrackPhi),kFALSE);
	  fQvec5 += TComplex(TMath::Cos(5*fTrackPhi),TMath::Sin(5*fTrackPhi),kFALSE);

    if(fDiffFlow) // do differential flow switch
    {
      if(fPtBinIndex != -1)
      {
        iNumP2[fPtBinIndex]++;

        fPvec2[fPtBinIndex] += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
        fPvec3[fPtBinIndex] += TComplex(TMath::Cos(3*fTrackPhi),TMath::Sin(3*fTrackPhi),kFALSE);

        if(fTrackEta > 0.)
        {
          fPvec2Gap00P[fPtBinIndex] += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
          iNumP2Gap00P[fPtBinIndex]++;
        }
        
        if(fTrackEta > 0.2)
        {
          fPvec2Gap04P[fPtBinIndex] += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
          iNumP2Gap04P[fPtBinIndex]++;
        }
        
        if(fTrackEta > 0.4)
        {
          fPvec2Gap08P[fPtBinIndex] += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
          iNumP2Gap08P[fPtBinIndex]++;
        }

        if(fTrackEta > 0.5)
        {
          fPvec2Gap10P[fPtBinIndex] += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
          iNumP2Gap10P[fPtBinIndex]++;
        }

      }
    }

    // eta gap
    if(fTrackEta > 0.)
    {
      fQvec2Gap00P += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
      iNumGap00P++;
    }

    if(fTrackEta < 0.)
    {
      fQvec2Gap00N += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
      iNumGap00N++;
    }
    
    if(fTrackEta > 0.2)
    {
      fQvec2Gap04P += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);;
      iNumGap04P++;
    }

    if(fTrackEta < -0.2)
    {
      fQvec2Gap04N += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
      iNumGap04N++;
    }
    
    if(fTrackEta > 0.4)
    {
      fQvec2Gap08P += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
      iNumGap08P++;
    }

    if(fTrackEta < -0.4)
    {
      fQvec2Gap08N += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);;
      iNumGap08N++;
    }
    
    if(fTrackEta > 0.5)
    {
      fQvec2Gap10P += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);;
      iNumGap10P++;
    }

    if(fTrackEta < -0.5)
    {
      fQvec2Gap10N += TComplex(TMath::Cos(2*fTrackPhi),TMath::Sin(2*fTrackPhi),kFALSE);
      iNumGap10N++;
    }    
  }   // end of loop over all tracks

  if(iNumTracksSelected < 2) return; // 0 tracks selected in this event

  // events with at least two selected track
  fEventCounter->Fill(7); 
  fCentralityDis->Fill(fCentPercentile);
  fMultTracksSelected->Fill(iNumTracksSelected);

  // Reference Flow
  Double_t dAmp = 0;
  Double_t dVal = 0;
  Double_t dWeight = 0;

  dWeight = iNumTracksSelected*(iNumTracksSelected-1);

  dAmp = (fQvec2*(TComplex::Conjugate(fQvec2))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCentPercentile > 0) )
    fRefCorTwo2->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dAmp = (fQvec3*(TComplex::Conjugate(fQvec3))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCentPercentile > 0) )
    fRefCorTwo3->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile

	dAmp = (fQvec4*(TComplex::Conjugate(fQvec4))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCentPercentile > 0) )
    fRefCorTwo4->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dAmp = (fQvec5*(TComplex::Conjugate(fQvec5))).Re();
  dVal = (dAmp - iNumTracksSelected) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCentPercentile > 0) )
    fRefCorTwo5->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile

  // ref flow with eta gap

  dWeight = iNumGap00P*iNumGap00N;
  dAmp = (fQvec2Gap00P*(TComplex::Conjugate(fQvec2Gap00N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCentPercentile > 0) )
    fRefCorTwo2Gap00->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile

  dWeight = iNumGap04P*iNumGap04N;
  dAmp = (fQvec2Gap04P*(TComplex::Conjugate(fQvec2Gap04N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCentPercentile > 0) )
    fRefCorTwo2Gap04->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dWeight = iNumGap08P*iNumGap08N;
  dAmp = (fQvec2Gap08P*(TComplex::Conjugate(fQvec2Gap08N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCentPercentile > 0) )
    fRefCorTwo2Gap08->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  
  dWeight = iNumGap10P*iNumGap10N;
  dAmp = (fQvec2Gap10P*(TComplex::Conjugate(fQvec2Gap10N))).Re();
  dVal = (dAmp) / dWeight;
  if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCentPercentile > 0) )
    fRefCorTwo2Gap10->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  

  // differential flow 
  if(fDiffFlow) // do a differential flow switch
  {
    for(Int_t i(0); i < fNumPtBins; i++)
    {
      dWeight = iNumP2[i]*(iNumTracksSelected-1);
      dAmp = (fPvec2[i]*(TComplex::Conjugate(fQvec2))).Re();
      dVal = (dAmp - iNumP2[i]) / dWeight;
      if( TMath::Abs(dVal < 1) && (dWeight > 1))
        fDiffCorTwo2[fCentBinIndex]->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2 ,dVal, dWeight);
      //else // test for not filling
        //printf("DiffFlow[%d] not filled: pT: %g dVal: %g dWeight: %g NumOfPart: %d dAmp: %g \n",fCentBinIndex,(fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal,dWeight, iNumP2[i],dAmp);  

      dWeight = iNumP2[i]*(iNumTracksSelected-1);
      dAmp = (fPvec3[i]*(TComplex::Conjugate(fQvec3))).Re();
      dVal = (dAmp - iNumP2[i]) / dWeight;
      if( TMath::Abs(dVal < 1) && (dWeight > 1))
        fDiffCorTwo3[fCentBinIndex]->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2 ,dVal, dWeight);

      dWeight = iNumP2Gap00P[i]*(iNumGap00N);
      dAmp = (fPvec2Gap00P[i]*(TComplex::Conjugate(fQvec2Gap00N))).Re();
      dVal = dAmp / dWeight;
      if(TMath::Abs(dVal< 1) && (dWeight > 0))
        fDiffCorTwo2Gap00[fCentBinIndex]->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal,dWeight);

      dWeight = iNumP2Gap04P[i]*(iNumGap04N);
      dAmp = (fPvec2Gap04P[i]*(TComplex::Conjugate(fQvec2Gap04N))).Re();
      dVal = dAmp / dWeight;
      if(TMath::Abs(dVal< 1) && (dWeight > 0))
        fDiffCorTwo2Gap04[fCentBinIndex]->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal,dWeight);     

      dWeight = iNumP2Gap08P[i]*(iNumGap08N);
      dAmp = (fPvec2Gap08P[i]*(TComplex::Conjugate(fQvec2Gap08N))).Re();
      dVal = dAmp / dWeight;
      if(TMath::Abs(dVal< 1) && (dWeight > 0))
        fDiffCorTwo2Gap08[fCentBinIndex]->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal,dWeight);     

      dWeight = iNumP2Gap10P[i]*(iNumGap10N);
      dAmp = (fPvec2Gap10P[i]*(TComplex::Conjugate(fQvec2Gap10N))).Re();
      dVal = dAmp / dWeight;
      if(TMath::Abs(dVal< 1) && (dWeight > 0))
        fDiffCorTwo2Gap10[fCentBinIndex]->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal,dWeight);     
    }
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
  
  fEventCounter->Fill(2);

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

  fEventCounter->Fill(3);

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

  fEventCounter->Fill(4);

  if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
	{
		if(fDebug) 
			::Info("IsEventSelected","Event rejected: PV z-distance cut");
  	return kFALSE;
	}

  fEventCounter->Fill(5);

	// centrality rejection
	
	if(fPbPb)
	{
		// not implemented yet
		EstimateCentrality(fAOD);
    if (fCentBinIndex < 0 || fCentBinIndex > fNumCentBins || fCentPercentile < 0 || fCentPercentile > 80)
      return kFALSE;
	}
  
  fEventCounter->Fill(6);

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
void AliAnalysisTaskFlowPID::EstimateCentrality(AliVEvent* ev)
{
  Short_t centrCode = -1;
  Float_t lPercentile = -100;
  Float_t V0M_Cent = -100, SPD_Cent = -100;
  
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

  fCentSPDvsV0M->Fill(V0M_Cent, SPD_Cent);
  
  for(Int_t i(0); i < fNumCentBins; i++)
  {
    if( (lPercentile > fCentBinEdges[i]) && (lPercentile <= fCentBinEdges[i+1]) )
      centrCode = i;
  }
  
  fCentPercentile = lPercentile;
  fCentBinIndex = centrCode;

  if(fLHC10h)
  {
    if( !( (fCentPercentile <= 80) && (fCentPercentile > 0) ) || !(TMath::Abs(V0M_Cent - SPD_Cent) < 5) )
    {
      fCentPercentile = -100.;
      fCentBinIndex = -1;
      return;      
    }

    if( (fCentBinIndex == 0) && !(TMath::Abs(V0M_Cent - SPD_Cent) < 1) ) // 10h check for 0-5% centrality
    {
      fCentPercentile = -100.;
      fCentBinIndex = -1;
      return;
    }

    if( (fCentBinIndex == 1) && !(TMath::Abs(V0M_Cent - SPD_Cent) < 1) ) // 10h check for 5-10% centrality
    {
      fCentPercentile = -100.;
      fCentBinIndex = -1;
      return;
    }
  }
  
  return;

  /*
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
  */
}
//_____________________________________________________________________________
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
//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowPID::GetPtBinIndex(const Double_t dPt)
{
  for(Int_t i(0); i < fNumPtBins; i++)
  {
    if( (dPt >= fPtBinEdges[i]) && (dPt < fPtBinEdges[i+1]) )
      return i;
  }

  return -1;
}