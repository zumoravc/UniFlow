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

//Double_t AliAnalysisTaskUniFlow::fPtBinEdges[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0}; // You, Katarina binning
//Double_t AliAnalysisTaskUniFlow::fPtBinEdges[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.,5.5,6}; // PID flow v2 JHEP paper
Double_t AliAnalysisTaskUniFlow::fPtBinEdges[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.,4.2,4.4,4.6,4.8,5.};

Double_t AliAnalysisTaskUniFlow::fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; //# PbPb
//Double_t AliAnalysisTaskUniFlow::fCentBinEdges[] = {0.,50.,100.,150};
//Double_t AliAnalysisTaskUniFlow::fMinvFlowBinEdgesK0s[] = {0.4,0.42,0.44,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.54,0.56,0.58,0.6};
Double_t AliAnalysisTaskUniFlow::fMinvFlowBinEdgesK0s[] = {0.4,0.425,0.45,0.47,0.49,0.495,0.5,0.505,0.51,0.53,0.55,0.575,0.6};
Double_t AliAnalysisTaskUniFlow::fMinvFlowBinEdgesLambda[] = {1.08,1.09,1.10,1.105,1.11,1.115,1.12,1.125,1.13,1.14,1.15,1.16};
Int_t AliAnalysisTaskUniFlow::fHarmonics[AliAnalysisTaskUniFlow::fNumHarmonics] = {2};
Double_t AliAnalysisTaskUniFlow::fEtaGap[AliAnalysisTaskUniFlow::fNumEtaGap] = {-1.,0.,0.8};
Short_t AliAnalysisTaskUniFlow::fTracksScanFB[AliAnalysisTaskUniFlow::fNumScanFB] = {1,2,4,5,8,16,32,64,96,128,256,512,768,1024};

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(0),
  fAnalType(kAOD),
  fColSystem(kPP),
  fPeriod(kNon),
  fSampling(kFALSE),
  fNumSamples(10),
  fFilterCharged(kFALSE),
  fFilterPID(kFALSE),
  fFilterV0s(kFALSE),

  fArrTrackRPF(0x0),
  fArrTrackPOI(0x0),
  fArrPion(0x0),
  fArrKaon(0x0),
  fArrProton(0x0),
  fArrK0s(0x0),
  fArrLambda(0x0),
  fArrALambda(0x0),


  fLHC10h(kTRUE),
  fTrigger(kFALSE),
  fRejectPileFromSPD(kFALSE),
  fUseIsPileUpFromSPD(kFALSE),
  fUsePlpMV(kFALSE),
  fRejectOutOfBunchPU(kFALSE),
  fCentFlag(0),
  fDoFlow(0),
  fDiffFlow(0),
  fTracksScan(0),
  fPID(0),
  fUseBayesPID(0),
  fDoV0s(0),
  fOldFlow(0),
  fUseOldCent(0),
  fDoGenFramKat(0),

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
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(0),

  fAnalType(kAOD),
  fColSystem(kPP),
  fPeriod(kNon),
  fSampling(kFALSE),
  fNumSamples(10),
  fFilterCharged(kFALSE),
  fFilterPID(kFALSE),
  fFilterV0s(kFALSE),

  fArrTrackRPF(0x0),
  fArrTrackPOI(0x0),
  fArrPion(0x0),
  fArrKaon(0x0),
  fArrProton(0x0),
  fArrK0s(0x0),
  fArrLambda(0x0),
  fArrALambda(0x0),

  fLHC10h(kTRUE),
  fTrigger(kFALSE),
  fRejectPileFromSPD(kFALSE),
  fUseIsPileUpFromSPD(kFALSE),
  fRejectOutOfBunchPU(kFALSE),
  fUsePlpMV(kFALSE),
  fCentFlag(0),
  fDoFlow(0),
  fDiffFlow(0),
  fTracksScan(0),
  fPID(0),
  fUseBayesPID(0),
  fDoV0s(0),
  fOldFlow(0),
  fDoGenFramKat(0),
  fUseOldCent(0),

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
  AliAODEvent* event = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!EventSelection(event)) return;

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
  if(!ProcessEvent(event)) return;

  // posting data (mandatory)
  PostData(1, fOutListEvents);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::EventSelection(const AliAODEvent* event)
{
  // main (envelope) method for event selection
  // Specific event selection methods are called from here
  // returns kTRUE if event pass all selection criteria

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::Filtering()
{
  // main (envelope) method for filtering all particles of interests (POIs) in selected events
  // All POIs passing selection criteria are saved to relevant TClonesArray for further processing
  // return kTRUE if succesfull (no errors in process)

  if(fFilterCharged)
  {
        // TODO :: implement§
  }

  if(fFilterPID)
  {
    // TODO:: implelemt
  }

  if(fFilterV0s)
  {
    // TODO implement
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::ProcessEvent(const AliAODEvent* event)
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
