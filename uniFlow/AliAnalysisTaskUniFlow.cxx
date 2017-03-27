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
  fPIDResponse(0x0),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(0),
  fEventCounter(0),
  fNumEventsAnalyse(50),

  fRunMode(kFull),
  fAnalType(kAOD),
  fColSystem(kPP),
  fPeriod(kNon),
  fTrigger(0),
  fSampling(kFALSE),
  fNumSamples(10),
  fProcessCharged(kFALSE),
  fProcessPID(kFALSE),
  fProcessV0s(kFALSE),

  fArrCharged(0x0),
  fArrChargedRPF(0x0),
  fArrChargedPOI(0x0),
  fArrPion(0x0),
  fArrKaon(0x0),
  fArrProton(0x0),
  fArrK0s(0x0),
  fArrLambda(0x0),
  fArrALambda(0x0),

  fPVtxCutZ(0.),
  fCutChargedEtaMax(0),
  fCutChargedPtMax(0),
  fCutChargedPtMin(0),
  fCutChargedDCAzMax(0),
  fCutChargedDCAxyMax(0),
  fCutChargedNumTPCclsMin(0),
  fCutChargedTrackFilterBit(0),

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
  fCutV0MinDCADaughters(0.),
  fCutV0MaxDCADaughters(0.),
  fCutV0MinDecayRadius(0.),
  fCutV0MaxDecayRadius(0.),
  fCutV0DaughterPtMin(0.),
  fCutV0DaughterPtMax(0.),
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

  // output lists
  fOutListEvents(0x0),
  fOutListCharged(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCounter(0x0),

  // charged histogram
  fhChargedCounter(0x0),

  // V0s histogram
  fhV0sCounter(0x0)

{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow(const char* name) : AliAnalysisTaskSE(name),
  fEventAOD(0x0),
  fPIDResponse(0x0),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(0),
  fEventCounter(0),
  fNumEventsAnalyse(50),

  fRunMode(kFull),
  fAnalType(kAOD),
  fColSystem(kPP),
  fPeriod(kNon),
  fTrigger(0),
  fSampling(kFALSE),
  fNumSamples(10),
  fProcessCharged(kFALSE),
  fProcessPID(kFALSE),
  fProcessV0s(kFALSE),

  fArrCharged(0x0),
  fArrChargedRPF(0x0),
  fArrChargedPOI(0x0),
  fArrPion(0x0),
  fArrKaon(0x0),
  fArrProton(0x0),
  fArrK0s(0x0),
  fArrLambda(0x0),
  fArrALambda(0x0),

  fPVtxCutZ(0.),
  fCutChargedEtaMax(0),
  fCutChargedPtMax(0),
  fCutChargedPtMin(0),
  fCutChargedDCAzMax(0),
  fCutChargedDCAxyMax(0),
  fCutChargedTrackFilterBit(0),
  fCutChargedNumTPCclsMin(0),

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
  fCutV0MinDCADaughters(0.),
  fCutV0MaxDCADaughters(0.),
  fCutV0MinDecayRadius(0.),
  fCutV0MaxDecayRadius(0.),
  fCutV0DaughterPtMin(0.),
  fCutV0DaughterPtMax(0.),
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

  // output lists
  fOutListEvents(0x0),
  fOutListCharged(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCounter(0x0),

  // charged histogram
  fhChargedCounter(0x0),

  // V0s histogram
  fhV0sCounter(0x0)

{
  // QA histograms
  for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
  {
    // Event histograms
    fhQAEventsPVz[iQA] = 0x0;
    fhQAEventsNumContrPV[iQA] = 0x0;
    fhQAEventsNumSPDContrPV[iQA] = 0x0;
    fhQAEventsDistPVSPD[iQA] = 0x0;
    fhQAEventsSPDresol[iQA] = 0x0;

    // charged
    fhQAChargedMult[iQA] = 0x0;
    fhQAChargedPt[iQA] = 0x0;
    fhQAChargedEta[iQA] = 0x0;
    fhQAChargedPhi[iQA] = 0x0;
    fhQAChargedCharge[iQA] = 0x0;
    fhQAChargedFilterBit[iQA] = 0x0;
    fhQAChargedNumTPCcls[iQA] = 0x0;
    fhQAChargedDCAxy[iQA] = 0x0;
    fhQAChargedDCAz[iQA] = 0x0;
    fhQAChargedTPCstatus[iQA] = 0x0;
    fhQAChargedTOFstatus[iQA] = 0x0;
    fhQAChargedTPCdEdx[iQA] = 0x0;
    fhQAChargedTOFbeta[iQA] = 0x0;
  }

  // defining input/output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::~AliAnalysisTaskUniFlow()
{
  // destructor
  // if(fPIDCombined)
  // {
  //   delete fPIDCombined;
  // }

  // deleting output lists
  if(fOutListEvents) delete fOutListEvents;
  if(fOutListCharged) delete fOutListCharged;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // *************************************************************

  // list all parameters used in this analysis
  ListParameters();

  // task initialization
  fInit = InitializeTask();
  if(!fInit) return;

  // Creating arrays for particles
  fArrCharged = new TClonesArray("AliAODTrack",10000);
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
  fOutListCharged = new TList();
  fOutListCharged->SetOwner(kTRUE);

  // creating histograms
    // event histogram
    fhEventSampling = new TH1D("fhEventSampling","Event sampling",fNumSamples, 0, fNumSamples);
    fOutListEvents->Add(fhEventSampling);

    const Short_t iEventCounterBins = 7;
    TString sEventCounterLabel[iEventCounterBins] = {"Input","Physics selection OK","PV OK","SPD Vtx OK","Pileup MV OK","PV #it{z} OK","Selected"};
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    for(Short_t i(0); i < iEventCounterBins; i++) fhEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    fOutListEvents->Add(fhEventCounter);

    // charged (tracks) histograms
    TString sChargedCounterLabel[] = {"Input","FB","#TPC-Cls","DCA-z","DCA-xy","Eta","Selected"};
    const Short_t iNBinsChargedCounter = sizeof(sChargedCounterLabel)/sizeof(sChargedCounterLabel[0]);
    fhChargedCounter = new TH1D("fhChargedCounter","Charged tracks: Counter",iNBinsChargedCounter,0,iNBinsChargedCounter);
    for(Short_t i(0); i < iNBinsChargedCounter; i++) fhChargedCounter->GetXaxis()->SetBinLabel(i+1, sChargedCounterLabel[i].Data() );
    fOutListCharged->Add(fhChargedCounter);

    // QA histograms
    TString sQAindex[fiNumIndexQA] = {"Before", "After"};
    for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
    {
      // EVENTs QA histograms
      fhQAEventsPVz[iQA] = new TH1D(Form("fhQAEventsPVz_%s",sQAindex[iQA].Data()), "QA Events: PV-#it{z}", 101,-50,50);
      fOutListEvents->Add(fhQAEventsPVz[iQA]);
      fhQAEventsNumContrPV[iQA] = new TH1D(Form("fhQAEventsNumContrPV_%s",sQAindex[iQA].Data()), "QA Events: Number of contributors to AOD PV", 20,0,20);
      fOutListEvents->Add(fhQAEventsNumContrPV[iQA]);
      fhQAEventsNumSPDContrPV[iQA] = new TH1D(Form("fhQAEventsNumSPDContrPV_%s",sQAindex[iQA].Data()), "QA Events: SPD contributors to PV", 20,0,20);
      fOutListEvents->Add(fhQAEventsNumSPDContrPV[iQA]);
      fhQAEventsDistPVSPD[iQA] = new TH1D(Form("fhQAEventsDistPVSPD_%s",sQAindex[iQA].Data()), "QA Events: PV SPD vertex", 50,0,5);
      fOutListEvents->Add(fhQAEventsDistPVSPD[iQA]);
      fhQAEventsSPDresol[iQA] = new TH1D(Form("fhQAEventsSPDresol_%s",sQAindex[iQA].Data()), "QA Events: SPD resolution", 150,0,15);
      fOutListEvents->Add(fhQAEventsSPDresol[iQA]);

      // Charged tracks QA
      const Short_t iNBinsPIDstatus = 4;
      TString sPIDstatus[iNBinsPIDstatus] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};

      fhQAChargedMult[iQA] = new TH1D(Form("fhQAChargedMult_%s",sQAindex[iQA].Data()),"QA Charged: Number of Charged in selected events; #it{N}^{Charged}", 1500,0,1500);
      fOutListCharged->Add(fhQAChargedMult[iQA]);

      fhQAChargedCharge[iQA] = new TH1D(Form("fhQAChargedCharge_%s",sQAindex[iQA].Data()),"QA Charged: Track charge; charge;", 3,-1.5,1.5);
      fOutListCharged->Add(fhQAChargedCharge[iQA]);
      fhQAChargedPt[iQA] = new TH1D(Form("fhQAChargedPt_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 300,0.,30.);
      fOutListCharged->Add(fhQAChargedPt[iQA]);
      fhQAChargedEta[iQA] = new TH1D(Form("fhQAChargedEta_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#eta}; #it{#eta}", 301,-3,3);
      fOutListCharged->Add(fhQAChargedEta[iQA]);
      fhQAChargedPhi[iQA] = new TH1D(Form("fhQAChargedPhi_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#varphi}; #it{#varphi}", 100,0.,TMath::TwoPi());
      fOutListCharged->Add(fhQAChargedPhi[iQA]);

      const Short_t iNFilterMapBinBins = 32;
      fhQAChargedFilterBit[iQA] = new TH1D(Form("fhQAChargedFilterBit_%s",sQAindex[iQA].Data()), "QA Charged: Filter bit",iNFilterMapBinBins,0,iNFilterMapBinBins);
      fOutListCharged->Add(fhQAChargedFilterBit[iQA]);
      for(Int_t j = 0; j < iNFilterMapBinBins; j++) fhQAChargedFilterBit[iQA]->GetXaxis()->SetBinLabel(j+1, Form("%g",TMath::Power(2,j)));
      fOutListCharged->Add(fhQAChargedFilterBit[iQA]);

      fhQAChargedNumTPCcls[iQA] = new TH1D(Form("fhQAChargedNumTPCcls_%s",sQAindex[iQA].Data()),"QA Charged: Track number of TPC clusters; #it{N}^{TPC clusters}", 160,0,160);
      fOutListCharged->Add(fhQAChargedNumTPCcls[iQA]);
      fhQAChargedDCAxy[iQA] = new TH1D(Form("fhQAChargedDCAxy_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-xy; DCA_{#it{xy}} (cm)", 100,0.,10);
      fOutListCharged->Add(fhQAChargedDCAxy[iQA]);
      fhQAChargedDCAz[iQA] = new TH1D(Form("fhQAChargedDCAz_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-z; DCA_{#it{z}} (cm)", 200,-10.,10.);
      fOutListCharged->Add(fhQAChargedDCAz[iQA]);
      fhQAChargedTPCstatus[iQA] = new TH1D(Form("fhQAChargedTPCstatus_%s",sQAindex[iQA].Data()),"QA Charged: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
      fOutListCharged->Add(fhQAChargedTPCstatus[iQA]);
      fhQAChargedTPCdEdx[iQA] = new TH2D(Form("fhQAChargedTPCdEdx_%s",sQAindex[iQA].Data()),"QA Charged: TPC PID information; #it{p} (GeV/#it{c}); TPC dEdx (au)", 50,0,10, 101,-10,1000);
      fOutListCharged->Add(fhQAChargedTPCdEdx[iQA]);
      fhQAChargedTOFstatus[iQA] = new TH1D(Form("fhQAChargedTOFstatus_%s",sQAindex[iQA].Data()),"QA Charged: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
      fOutListCharged->Add(fhQAChargedTOFstatus[iQA]);
      fhQAChargedTOFbeta[iQA] = new TH2D(Form("fhQAChargedTOFbeta_%s",sQAindex[iQA].Data()),"QA Charged: TOF #beta information; #it{p} (GeV/#it{c}); TOF #beta", 50,0,10, 80,-0.1,1.5);
      fOutListCharged->Add(fhQAChargedTOFbeta[iQA]);

      for(Int_t j = 0; j < iNBinsPIDstatus; j++)
      {
        fhQAChargedTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
        fhQAChargedTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
      }
    }

  // posting data (mandatory)
  PostData(1, fOutListEvents);
  PostData(2, fOutListCharged);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListParameters()
{
  // lists all task parameters
  // *************************************************************
  printf("\n======= List of parameters ==================\n");
  printf("   -------- Analysis task ------------------\n");
  printf("   -------- Events -------------------------\n");
  printf("   -------- Charge tracks ------------------\n");
  printf("   fCutChargedTrackFilterBit: (UInt) %d\n", fCutChargedTrackFilterBit);
  printf("   fCutChargedNumTPCclsMin: (UShort_t) %d\n", fCutChargedNumTPCclsMin);
  printf("   fCutChargedEtaMax: (Float_t) %g\n", fCutChargedEtaMax);
  printf("   fCutChargedPtMin: (Float_t) %g (GeV/c)\n", fCutChargedPtMin);
  printf("   fCutChargedPtMax: (Float_t) %g (GeV/c)\n", fCutChargedPtMax);
  printf("   fCutChargedDCAzMax: (Float_t) %g (cm)\n", fCutChargedDCAzMax);
  printf("   fCutChargedDCAxyMax: (Float_t) %g (cm)\n", fCutChargedDCAxyMax);
  printf("=============================================\n\n");

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::InitializeTask()
{
  // called once on beginning of task (within UserCreateOutputObjects method)
  // check if task parameters are specified and valid
  // returns kTRUE if succesfull
  // *************************************************************

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

  // checking PID response
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse)
  {
    ::Error("InitializeTask","AliPIDResponse object not found! Terminating!");
    return kFALSE;
  }

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
  // *************************************************************

  if(!fInit) return; // check if initialization succesfull

  // local event counter check: if running in test mode, it runs until the 50 events are succesfully processed
  if(fRunMode == kTest && fEventCounter >= fNumEventsAnalyse) return;

  // event selection
  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!EventSelection()) return;

  // estimate centrality & assign indexes (centrality/percentile, sampling, ...)
  fIndexSampling = GetSamplingIndex();
  fhEventSampling->Fill(fIndexSampling);

  fIndexCentrality = GetCentralityIndex();
  // TODO implement

  // processing of selected event
  if(!ProcessEvent()) return;

  // posting data (mandatory)
  PostData(1, fOutListEvents);
  PostData(2, fOutListCharged);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::EventSelection()
{
  // main (envelope) method for event selection
  // Specific event selection methods are called from here
  // returns kTRUE if event pass all selection criteria
  // *************************************************************

  Bool_t eventSelected = kFALSE;

  if(!fEventAOD) return kFALSE;

  // Fill event QA BEFORE cuts
  FillEventsQA(0);

  // event selection for small systems pp, pPb in Run2 (2016)
  if( (fColSystem == kPP || fColSystem == kPPb)
      && (fPeriod == k16k || fPeriod == k16l || fPeriod == k16q || fPeriod == k16r || fPeriod == k16s || fPeriod == k16t)
    ) eventSelected = IsEventSelected_2016();

  // Fill event QA AFTER cuts
  FillEventsQA(1);

  return eventSelected;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsEventSelected_2016()
{
  // Event selection for small system collision recorder in Run 2 year 2016
  // pp (LHC16kl), pPb (LHC16rqts)
  // return kTRUE if event passes all criteria, kFALSE otherwise
  // *************************************************************

  fhEventCounter->Fill("Input",1);

  // Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();

  Bool_t isTriggerSelected = kFALSE;
  switch(fTrigger) // check for high multiplicity trigger
  {
    case 0:
      isTriggerSelected = fSelectMask& AliVEvent::kINT7;
      break;

    case 1:
      isTriggerSelected = fSelectMask& AliVEvent::kHighMultV0;
      break;

    case 2:
      isTriggerSelected = fSelectMask& AliVEvent::kHighMultSPD;
      break;

    default: isTriggerSelected = kFALSE;
  }

  if(!isTriggerSelected)
    return kFALSE;

  // events passing physics selection
  fhEventCounter->Fill("Physics selection OK",1);

  // primary vertex selection
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
  fhEventCounter->Fill("PV OK",1);

  // SPD vertex selection
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertexSPD());

  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
  fhEventCounter->Fill("SPD Vtx OK",1);

  // PileUp rejection included in Physics selection
  // but with values for high mult pp (> 5 contrib) => for low ones: do manually (> 3 contrib)

  /*
  if(fTrigger == 0 && fAOD->IsPileupFromSPD(3,0.8) )
  {
    return kFALSE;
  }
  */

  //fhEventCounter->Fill("Pileup SPD OK",1);

  // pileup rejection from multivertexer
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5);
  utils.SetMaxPlpChi2MV(5);
  utils.SetMinWDistMV(15);
  utils.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t isPileupFromMV = utils.IsPileUpMV(fEventAOD);

  if(isPileupFromMV) return kFALSE;
  fhEventCounter->Fill("Pileup MV OK",1);

  // if(fRejectOutOfBunchPU) // out-of-bunch rejection (provided by Christian)
  // {
  //   //out-of-bunch 11 BC
  //   if (utils.IsOutOfBunchPileUp(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //   fhEventCounter->Fill("OOBPU OK",1);
  //
  //   if (utils.IsSPDClusterVsTrackletBG(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //
  //   fhEventCounter->Fill("SPDClTrBG OK",1);
  //
  //   // SPD pileup
  //   if (utils.IsPileUpSPD(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //
  //   fhEventCounter->Fill("SPDPU OK",1);
  // }

  //fhEventCounter->Fill("Utils OK",1);

  // cutting on PV z-distance
  const Double_t aodVtxZ = vtx->GetZ();
  if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
  {
    return kFALSE;
  }
  fhEventCounter->Fill("PV #it{z} OK",1);

  fhEventCounter->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillEventsQA(const Short_t iQAindex)
{
  // Filling various QA plots related with event selection
  // *************************************************************

  const AliAODVertex* aodVtx = fEventAOD->GetPrimaryVertex();
  const Double_t dVtxZ = aodVtx->GetZ();
  const Int_t iNumContr = aodVtx->GetNContributors();
  const AliAODVertex* spdVtx = fEventAOD->GetPrimaryVertexSPD();
  const Int_t iNumContrSPD = spdVtx->GetNContributors();
  const Double_t spdVtxZ = spdVtx->GetZ();

  fhQAEventsPVz[iQAindex]->Fill(dVtxZ);
  fhQAEventsNumContrPV[iQAindex]->Fill(iNumContr);
  fhQAEventsNumSPDContrPV[iQAindex]->Fill(iNumContrSPD);
  fhQAEventsDistPVSPD[iQAindex]->Fill(TMath::Abs(dVtxZ - spdVtxZ));

  // // event / physics selection criteria
  // AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  // UInt_t fSelectMask = inputHandler->IsEventSelected();
  //
  // if( fSelectMask& AliVEvent::kINT7 ) { fQAEventsTriggerSelection[iQAindex]->Fill("kINT7",1); }
  // else if (fSelectMask& AliVEvent::kHighMultV0) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultV0",1); }
  // else if (fSelectMask& AliVEvent::kHighMultSPD) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultSPD",1); }
  // else { fQAEventsTriggerSelection[iQAindex]->Fill("Other",1); }

  // SPD vertexer resolution
  Double_t cov[6] = {0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  fhQAEventsSPDresol[iQAindex]->Fill(zRes);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::Filtering()
{
  // main (envelope) method for filtering all particles of interests (POIs) in selected events
  // All POIs passing selection criteria are saved to relevant TClonesArray for further processing
  // return kTRUE if succesfull (no errors in process)
  // *************************************************************

  if(!fProcessCharged && !fProcessPID && !fProcessV0s) // if neither is ON, filtering is skipped
    return kFALSE;

  if(fProcessCharged)
  {
    fArrCharged->Clear("C");
    //fArrChargedRPF->Clear("C");
    //fArrChargedPOI->Clear("C");
    if(!FilterCharged()) return kFALSE;
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
Bool_t AliAnalysisTaskUniFlow::FilterCharged()
{
  Short_t iNumTrackSelected = 0; // selected charged track counter

  const Short_t iNumTracks = fEventAOD->GetNumberOfTracks();
  if(iNumTracks < 1) return kFALSE;

  AliAODTrack* track = 0x0;
  for(Short_t iTrack(0); iTrack < iNumTracks; iTrack++)
  {
    track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
    FillQACharged(0,track); // QA before selection

    if(IsChargedSelected(track))
    {
      new((*fArrCharged)[iNumTrackSelected++]) AliAODTrack(*track);
      FillQACharged(1,track); // QA after selection
    }
  }

  // fill QA charged multiplicity
  fhQAChargedMult[0]->Fill(fEventAOD->GetNumberOfTracks());
  fhQAChargedMult[1]->Fill(fArrCharged->GetEntriesFast());

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsChargedSelected(const AliAODTrack* track)
{
  // Selection of charged track
  // returns kTRUE if track pass all requirements, kFALSE otherwise
  // *************************************************************
  if(!track) return kFALSE;
  fhChargedCounter->Fill("Input",1);

  // filter bit
  if( !track->TestFilterBit(fCutChargedTrackFilterBit) ) return kFALSE;
  fhChargedCounter->Fill("FB",1);

  // number of TPC clusters (additional check for not ITS-standalone tracks)
  if( track->GetTPCNcls() < fCutChargedNumTPCclsMin && fCutChargedTrackFilterBit != 2) return kFALSE;
  fhChargedCounter->Fill("#TPC-Cls",1);

  // track DCA coordinates
  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  Double_t dTrackXYZ[3] = {0};
  Double_t dVertexXYZ[3] = {0.};
  Double_t dDCAXYZ[3] = {0.};
  if( fCutChargedDCAzMax > 0. || fCutChargedDCAxyMax > 0.)
  {
    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) return kFALSE; // event does not have a PV

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++)
      dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
  }

  if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fCutChargedDCAzMax) return kFALSE;
  fhChargedCounter->Fill("DCA-z",1);

  if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) return kFALSE;
  fhChargedCounter->Fill("DCA-xy",1);

  // pseudorapidity (eta)
  if(fCutChargedEtaMax > 0. && TMath::Abs(track->Eta()) > fCutChargedEtaMax) return kFALSE;
  fhChargedCounter->Fill("Eta",1);

  // track passing all criteria
  fhChargedCounter->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQACharged(const Short_t iQAindex, const AliAODTrack* track)
{
  // Filling various QA plots related to charged track selection
  // *************************************************************
  if(!track) return;

  // filter bit testing
  for(Short_t i(0); i < 32; i++)
  {
    if(track->TestFilterBit(TMath::Power(2.,i)))
      fhQAChargedFilterBit[iQAindex]->Fill(i);
  }

  // track charge
  fhQAChargedCharge[iQAindex]->Fill(track->Charge());

  // number of TPC clusters
  fhQAChargedNumTPCcls[iQAindex]->Fill(track->GetTPCNcls());

  // track DCA
  Double_t dDCAXYZ[3] = {-999., -999., -999.};
  const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
  if(vertex)
  {
    Double_t dTrackXYZ[3] = {-999., -999., -999.};
    Double_t dVertexXYZ[3] = {-999., -999., -999.};

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++)
      dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
  }
  fhQAChargedDCAxy[iQAindex]->Fill(TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]));
  fhQAChargedDCAz[iQAindex]->Fill(dDCAXYZ[2]);

  // TPC & TOF statuses & measures
  if(fPIDResponse)
  {
    AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
    AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);

    fhQAChargedTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
    fhQAChargedTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );

    if(pidStatusTPC == AliPIDResponse::kDetPidOk)
    fhQAChargedTPCdEdx[iQAindex]->Fill(track->P(), track->GetTPCsignal());
    else // TPC status not OK
    fhQAChargedTPCdEdx[iQAindex]->Fill(track->P(), -5.);

    if(pidStatusTOF == AliPIDResponse::kDetPidOk && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME))
    {
      Double_t dTOF[5];
      track->GetIntegratedTimes(dTOF);
      Double_t dBetaTOF = dTOF[0] / track->GetTOFsignal();
      fhQAChargedTOFbeta[iQAindex]->Fill(track->P(),dBetaTOF);
    }
    else // TOF status not OK
      fhQAChargedTOFbeta[iQAindex]->Fill(track->P(),-0.05);
  }

  // kinematics
  fhQAChargedPt[iQAindex]->Fill(track->Pt());
  fhQAChargedPhi[iQAindex]->Fill(track->Phi());
  fhQAChargedEta[iQAindex]->Fill(track->Eta());

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::FilterV0s()
{
  // Filtering of V0s candidates
  // returns kFALSE if any complications occurs
  // *************************************************************

  Short_t iNumLambdaSelected = 0;
  Short_t iNumK0sSelected = 0;

  const Short_t iNumV0s = fEventAOD->GetNumberOfV0s();
  if(iNumV0s < 1)
    return kFALSE;

  AliAODv0* v0 = 0x0;
  for(Short_t iV0(0); iV0 < iNumV0s; iV0++)
  {
    // the minimalistic dynamic allocation of the TClonesArray*
    v0 = static_cast<AliAODv0*>(fEventAOD->GetV0(iV0));
    if(!v0) continue;

    if(IsV0aK0s(v0))
      new((*fArrK0s)[iNumK0sSelected++]) AliAODv0(*v0);

    if(IsV0aLambda(v0))
      new((*fArrLambda)[iNumLambdaSelected++]) AliAODv0(*v0);

    if(IsV0aK0s(v0) && IsV0aLambda(v0))
      ::Warning("FilterV0s","V0 passing selection for both K0s and (A)Lambda");
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0aK0s(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for K0s candidates
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;

  if(!IsV0Selected(v0)) return kFALSE; // not passing common selection

  // K0s specific criteria


  // passing all criteria
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0aLambda(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for Lambda candidates
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;

  if(!IsV0Selected(v0)) return kFALSE; // not passing common selection

  // (Anti-)Lambda specific criteria


  // passing all criteria
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0Selected(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // common for both K0s and (Anti)-Lambdas
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;

  const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  // daughter track check
  if(!daughterPos || !daughterNeg) return kFALSE;

  // daughters & mother charge checks
  if( (TMath::Abs(daughterPos->Charge()) != 1) || (TMath::Abs(daughterNeg->Charge()) != 1) ) return kFALSE;
  if(daughterPos->Charge() == daughterNeg->Charge()) return kFALSE;
  if(v0->Charge() != 0) return kFALSE;

  // reconstruction method: online (on-the-fly) OR offline
  if(v0->GetOnFlyStatus() != fCutV0onFly) return kFALSE;

  // TPC refit
  if(fCutV0refitTPC && ( !daughterPos->IsOn(AliAODTrack::kTPCrefit) || !daughterNeg->IsOn(AliAODTrack::kTPCrefit) ) ) return kFALSE;

  // Kinks
  const AliAODVertex* prodVtxDaughterPos = (AliAODVertex*) daughterPos->GetProdVertex(); // production vertex of the positive daughter track
  const AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*) daughterNeg->GetProdVertex(); // production vertex of the negative daughter track
  if(fCutV0rejectKinks && ( (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) || (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) ) ) return kFALSE;

  // Daughters DCA to PV
  const Float_t dDCAPosToPV = TMath::Abs(v0->DcaPosToPrimVertex());
  const Float_t dDCANegToPV = TMath::Abs(v0->DcaNegToPrimVertex());
  if(fCutV0MinDCAtoPV > 0. && ( dDCAPosToPV < fCutV0MinDCAtoPV || dDCANegToPV < fCutV0MinDCAtoPV ) ) return kFALSE;
  if(fCutV0MaxDCAtoPV > 0. && ( dDCAPosToPV > fCutV0MaxDCAtoPV || dDCANegToPV > fCutV0MaxDCAtoPV ) ) return kFALSE;

  // Daughter DCA among themselves
  if(fCutV0MinDCADaughters > 0. && TMath::Abs(v0->DcaV0Daughters()) < fCutV0MinDCADaughters) return kFALSE;
  if(fCutV0MaxDCADaughters > 0. && TMath::Abs(v0->DcaV0Daughters()) > fCutV0MaxDCADaughters) return kFALSE;

  // radius of decay vertex in transverse plane
  Double_t dSecVtxCoor[3] = {0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  if( fCutV0MinDecayRadius > 0. && (dDecayRadius < fCutV0MinDecayRadius) ) return kFALSE;
  if( fCutV0MaxDecayRadius > 0. && (dDecayRadius > fCutV0MaxDecayRadius) ) return kFALSE;

  // acceptance checks
  if(fCutV0DaughterEtaMax > 0. && ( (TMath::Abs(daughterNeg->Eta()) > fCutV0DaughterEtaMax) || (TMath::Abs(daughterPos->Eta()) > fCutV0DaughterEtaMax) ) ) return kFALSE;
  if(fCutV0DaughterPtMin > 0. && (daughterPos->Pt() < fCutV0DaughterPtMin  || daughterNeg->Pt() < fCutV0DaughterPtMin) ) return kFALSE;
  if(fCutV0DaughterPtMax > 0. && (daughterPos->Pt() > fCutV0DaughterPtMax  || daughterNeg->Pt() > fCutV0DaughterPtMax) ) return kFALSE;

  if(fCutV0MotherEtaMax > 0. && TMath::Abs(v0->Eta()) > fCutV0DaughterEtaMax ) return kFALSE;
  if(fCutV0MotherPtMin > 0. && v0->Pt() < fCutV0MotherPtMin) return kFALSE;
  if(fCutV0MotherPtMax > 0. && v0->Pt() > fCutV0MotherPtMax) return kFALSE;

  // passing all common criteria
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
  // *************************************************************

  // filtering particles
  if(!Filtering()) return kFALSE;

  fEventCounter++; // counter of processed events
  return kTRUE;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetSamplingIndex()
{
  // Assing sampling index based on generated random number
  // returns centrality index
  // *************************************************************

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
  // Assing centrality index based on provided method (centrality estimator / number of selected tracks)
  // If number of selected track method is selected, track filtering must be done first
  // returns centrality index
  // *************************************************************

  Short_t index = 0;

  // centrality estimation for pPb analysis in Run2
  // TODO : just moved from flow Event selection for 2016
  if(fColSystem == kPPb)
  {
    Float_t lPercentile = 900;
    AliMultSelection* MultSelection = 0x0;
    MultSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");

    if(!MultSelection)
    {
      //If you get this warning (and lPercentiles 900) please check that the AliMultSelectionTask actually ran (before your task)
      ::Warning("GetCentralityIndex","AliMultSelection object not found!");
    }
    else
    {
      lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    }

    // TODO implement centrality selection based on lPercentile
    ::Warning("GetCentralityIndex","Centrality selection not implemented");
  }



  return index;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::Terminate(Option_t* option)
{
  // called on end of task, after all events are processed
  // *************************************************************

  return;
}
