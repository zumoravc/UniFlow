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

/* AliAnalysisTaskUniFlow - ALICE Unified Flow framework
*
* ALICE analysis task for universal study of flow.
* Note: So far implemented only for AOD analysis!
*
* Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
*/

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TChain.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TList.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVTrack.h"
#include "TComplex.h"
#include "TRandom3.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAnalysisTaskUniFlow.h"
#include "AliLog.h"

class AliAnalysisTaskUniFlow;

ClassImp(AliAnalysisTaskUniFlow) // classimp: necessary for root

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fEventAOD(0x0),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(0),
  fAnalType(kAOD),
  fColSystem(kPP),
  fPeriod(kNon),
  fSampling(kFALSE),
  fNumSamples(10),
  fProcessCharged(kFALSE),
  fProcessPID(kFALSE),
  fProcessV0s(kFALSE),

  fArrChargedRPF(0x0),
  fArrChargedPOI(0x0),
  fArrPion(0x0),
  fArrKaon(0x0),
  fArrProton(0x0),
  fArrK0s(0x0),
  fArrLambda(0x0),
  fArrALambda(0x0),

  fPVtxCutZ(0.),
  fTrackEtaMax(0),
  fTrackPtMax(0),
  fTrackPtMin(0),
  fTracksDCAzMax(0),
  fNumTPCclsMin(0),
  fTrackFilterBit(0),
  fCutPionNumSigmaMax(0),
  fCutKaonNumSigmaMax(0),
  fCutProtonNumSigmaMax(0),
  fCutBayesPIDPionMin(0),
  fCutBayesPIDKaonMin(0),
  fCutBayesPIDProtonMin(0),

  fCutV0onFly(kFALSE),
  fCutV0rejectKinks(kFALSE),
  fCutV0refitTPC(kTRUE),
  fCutV0MinCPALambda(0.),
  fCutV0MinCPAK0s(0.),
  fCutV0MinDCAtoPV(0.),
  fCutV0MaxDCAtoPV(0.),
  fCutV0MaxDCADaughters(0.),
  fCutV0MinDecayRadius(0.),
  fCutV0MaxDecayRadius(0.),
  fCutV0DaughterPtMin(0.),
  fCutV0DaughterEtaMax(0.),
  fCutV0MotherEtaMax(0.),
  fCutV0MotherRapMax(0.),
  fCutV0MotherPtMin(0.),
  fCutV0MotherPtMax(0.),
  fCutV0NumTauK0sMax(0.),
  fCutV0K0sArmenterosAlphaMin(0.),
  fCutV0NumTauLambdaMax(0.),
  fCutV0ProtonNumSigmaMax(0),
  fCutV0ProtonPIDPtMax(0.),

  fhEventSampling(0x0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow(const char* name) : AliAnalysisTaskSE(name),
  fEventAOD(0x0),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(0),

  fAnalType(kAOD),
  fColSystem(kPP),
  fPeriod(kNon),
  fSampling(kFALSE),
  fNumSamples(10),
  fProcessCharged(kFALSE),
  fProcessPID(kFALSE),
  fProcessV0s(kFALSE),

  fArrChargedRPF(0x0),
  fArrChargedPOI(0x0),
  fArrPion(0x0),
  fArrKaon(0x0),
  fArrProton(0x0),
  fArrK0s(0x0),
  fArrLambda(0x0),
  fArrALambda(0x0),

  fPVtxCutZ(0.),
  fTrackEtaMax(0),
  fTrackPtMax(0),
  fTrackPtMin(0),
  fTracksDCAzMax(0),
  fNumTPCclsMin(0),
  fTrackFilterBit(0),
  fCutPionNumSigmaMax(0),
  fCutKaonNumSigmaMax(0),
  fCutProtonNumSigmaMax(0),
  fCutBayesPIDPionMin(0),
  fCutBayesPIDKaonMin(0),
  fCutBayesPIDProtonMin(0),

  fCutV0onFly(kFALSE),
  fCutV0rejectKinks(kFALSE),
  fCutV0refitTPC(kTRUE),
  fCutV0MinCPALambda(0.),
  fCutV0MinCPAK0s(0.),
  fCutV0MinDCAtoPV(0.),
  fCutV0MaxDCAtoPV(0.),
  fCutV0MaxDCADaughters(0.),
  fCutV0MinDecayRadius(0.),
  fCutV0MaxDecayRadius(0.),
  fCutV0DaughterPtMin(0.),
  fCutV0DaughterEtaMax(0.),
  fCutV0MotherEtaMax(0.),
  fCutV0MotherRapMax(0.),
  fCutV0MotherPtMin(0.),
  fCutV0MotherPtMax(0.),
  fCutV0NumTauK0sMax(0.),
  fCutV0K0sArmenterosAlphaMin(0.),
  fCutV0NumTauLambdaMax(0.),
  fCutV0ProtonNumSigmaMax(0),
  fCutV0ProtonPIDPtMax(0.),

  fhEventSampling(0x0)
{
  // defining input/output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::~AliAnalysisTaskUniFlow()
{
  // // destructor
  // if(fPIDCombined)
  // {
  //   delete fPIDCombined;
  // }
  //
  // if(fOutListCumulants)
  // {
  //   delete fOutListCumulants;
  // }
  //
  // if(fOutListEvents)
  // {
  //   delete fOutListEvents;
  // }
  //
  // if(fOutListTracks)
  // {
  //   delete fOutListTracks;
  // }
  //
  // if(fOutListPID)
  // {
  //   delete fOutListPID;
  // }
  //
  // if(fOutListV0s)
  // {
  //   delete fOutListV0s;
  // }
  //
  // /*
  // if(fOutListQA)
  // {
  //   delete fOutListQA;
  // }
  // */
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)

  // list all parameters used in this analysis
  ListParameters();

  // task initialization
  fInit = InitializeTask();
  if(!fInit) return;

  // Creating arrays for particles
  fArrChargedRPF = new TClonesArray("AliAODTrack",10000);
  fArrChargedPOI = new TClonesArray("AliAODTrack",10000);
  fArrPion = new TClonesArray("AliAODTrack",5000);
  fArrKaon = new TClonesArray("AliAODTrack",5000);
  fArrProton = new TClonesArray("AliAODTrack",5000);
  fArrK0s = new TClonesArray("AliAODv0",5000);
  fArrLambda = new TClonesArray("AliAODv0",5000);
  fArrALambda = new TClonesArray("AliAODv0",5000);

  // creating output lists
  fOutListEvents = new TList();
  fOutListEvents->SetOwner(kTRUE);

  // creating histograms
  fhEventSampling = new TH1D("fhEventSampling","Event sampling",fNumSamples, 0, fNumSamples);
  fOutListEvents->Add(fhEventSampling);

  // posting data (mandatory)
  PostData(1, fOutListEvents);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListParameters()
{
  // lists all task parameters

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::InitializeTask()
{
  // called once on beginning of task (within UserCreateOutputObjects method)
  // check if task parameters are specified and valid
  // returns kTRUE if succesfull

  if(fAnalType != kESD && fAnalType != kAOD)
  {
    ::Error("InitializeTask","Analysis type not specified! Terminating!");
    return kFALSE;
  }

  if(fAnalType == kESD)
  {
    ::Error("InitializeTask","Analysis type: ESD not implemented! Terminating!");
    return kFALSE;
  }

  if(fColSystem != kPP && fColSystem != kPP && fColSystem != kPbPb)
  {
    ::Error("InitializeTask","Collisional system not specified! Terminating!");
    return kFALSE;
  }

  if(fPeriod == kNon)
  {
    ::Error("InitializeTask","Period of data sample not selected! Terminating!");
    return kFALSE;
  }

  // TODO check if period corresponds to selected collisional system

  if(fSampling && fNumSamples == 0)
  {
    ::Error("InitializeTask","Sampling used, but number of samples is 0! Terminating!");
    return kFALSE;
  }

  if(!fSampling)
    fNumSamples = 1;

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::UserExec(Option_t *)
{
  // main method called for each event (event loop)

  if(!fInit) return; // check if initialization succesfull

  // Fill event QA BEFORE cuts
  // TODO

  // event selection
  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!EventSelection()) return;

  // Fill event QA AFTER cuts
  // TODO

  // filter all particle of interest in given (selected) events
  // FilterTracks();

  // estimate centrality & assign indexes (centrality/percentile, sampling, ...)
  fIndexSampling = GetSamplingIndex();
  fhEventSampling->Fill(fIndexSampling);

  fIndexCentrality = GetCentralityIndex();
  // TODO implement

  // processing of selected event
  if(!ProcessEvent()) return;

  // posting data (mandatory)
  PostData(1, fOutListEvents);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::EventSelection()
{
  // main (envelope) method for event selection
  // Specific event selection methods are called from here
  // returns kTRUE if event pass all selection criteria

  if(!fEventAOD) return kFALSE;


  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::Filtering()
{
  // main (envelope) method for filtering all particles of interests (POIs) in selected events
  // All POIs passing selection criteria are saved to relevant TClonesArray for further processing
  // return kTRUE if succesfull (no errors in process)

  if(!fProcessCharged && !fProcessPID && !fProcessV0s) // if neither is ON, filtering is skipped
    return kFALSE;

  fArrChargedRPF->Clear("C");

  if(fProcessCharged)
  {
    fArrChargedPOI->Clear("C");
    // TODO
  }

  if(fProcessPID)
  {
    fArrPion->Clear("C");
    fArrKaon->Clear("C");
    fArrProton->Clear("C");

    if(!FilterPID())
      return kFALSE;

  }

  if(fProcessV0s)
  {
    fArrK0s->Clear("C");
    fArrLambda->Clear("C");
    fArrALambda->Clear("C");

    if(!FilterV0s())
      return kFALSE;
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::FilterV0s()
{
  Short_t iNumLambdaSelected = 0;

  const Short_t iNumV0s = fEventAOD->GetNumberOfV0s();
  if(iNumV0s < 1)
    return kFALSE;

  AliAODv0* v0 = 0x0;
  for(Short_t iV0(0); iV0 < iNumV0s; iV0++)
  {
    // the minimalistic dynamic allocation of the TClonesArray*
    v0 = static_cast<AliAODv0*>(fEventAOD->GetV0(iV0));
    printf("%d",iNumLambdaSelected);
    new((*fArrLambda)[iNumLambdaSelected++]) AliAODv0(*v0);
  }
  printf("Lambda selected: %d\n", iNumLambdaSelected);


  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::FilterPID()
{

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::ProcessEvent()
{
  // main method for processing of (selected) event
  // returns kTRUE if succesfull

  if(!Filtering()) return kFALSE;


  return kTRUE;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetSamplingIndex()
{
  // Assing sampling index based on generated random number
  // returns centrality index
  Short_t index = 0;

  if(fSampling && fNumSamples > 0)
  {
    TRandom3 rr(0);
    Double_t ranNum = rr.Rndm(); // getting random number in (0,1)
    Double_t generated = ranNum * fNumSamples; // getting random number in range (0, fNumSamples)

    // finding right index for sampling based on generated number and total number of samples
    for(Short_t i(0); i < fNumSamples; i++)
    {
      if(generated < (i+1) )
      {
        index = i;
        break;
      }
    }
  }

  return index;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetCentralityIndex()
{
  Short_t index = 0;
  // Assing centrality index based on provided method (centrality estimator / number of selected tracks)
  // If number of selected track method is selected, track filtering must be done first
  // returns centrality index

  return index;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::Terminate(Option_t* option)
{
  // called on end of task, after all events are processed

  return;
}
