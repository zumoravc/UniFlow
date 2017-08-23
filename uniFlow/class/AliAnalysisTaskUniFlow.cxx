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

// AliAnalysisTaskUniFlow - ALICE Unified Flow framework
//
// ALICE analysis task for universal study of flow.
// Note: So far implemented only for AOD analysis!
//
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
//

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TList.h"
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
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVTrack.h"

class AliAnalysisTaskUniFlow;

ClassImp(AliAnalysisTaskUniFlow); // classimp: necessary for root

Int_t AliAnalysisTaskUniFlow::fHarmonics[] = {2,3};
Double_t AliAnalysisTaskUniFlow::fEtaGap[] = {-1.,0.,0.8};

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fEventAOD(0x0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fFlowWeightsFile(0x0),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(-1),
  fEventCounter(0),
  fNumEventsAnalyse(50),
  fRunNumber(-1),
  fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
  fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
  fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),
  fPDGMassPhi(TDatabasePDG::Instance()->GetParticle(333)->Mass()),
  fPDGMassK0s(TDatabasePDG::Instance()->GetParticle(310)->Mass()),
  fPDGMassLambda(TDatabasePDG::Instance()->GetParticle(3122)->Mass()),

  // FlowPart containers
  fVectorCharged(0x0),
  fVectorPion(0x0),
  fVectorKaon(0x0),
  fVectorProton(0x0),
  fVectorK0s(0x0),
  fVectorLambda(0x0),
  fVectorPhi(0x0),

  // analysis selection
  fRunMode(kFull),
  fAnalType(kAOD),
  fSampling(kFALSE),
  fFillQA(kTRUE),
  //fNumSamples(10),
  fProcessCharged(kFALSE),
  fProcessPID(kFALSE),
  fProcessV0s(kFALSE),
  fProcessPhi(kFALSE),

  // flow related
  fCutFlowRFPsPtMin(0),
  fCutFlowRFPsPtMax(0),
  fCutFlowDoFourCorrelations(kTRUE),
  fFlowFillWeights(kFALSE),
  fFlowPOIsPtMin(0),
  fFlowPOIsPtMax(20.),
  fFlowCentMin(0),
  fFlowCentMax(150),
  fFlowCentNumBins(150),
  fFlowWeightsPath(),
  fFlowUseWeights(kFALSE),

  // events selection
  fPVtxCutZ(0.),
  fColSystem(kPP),
  fPeriod(kNon),
  fMultEstimator(),
  fTrigger(0),

  // charged tracks selection
  fCutChargedEtaMax(0),
  fCutChargedPtMax(0),
  fCutChargedPtMin(0),
  fCutChargedDCAzMax(0),
  fCutChargedDCAxyMax(0),
  fCutChargedTrackFilterBit(0),
  fCutChargedNumTPCclsMin(0),

  // PID tracks selection
  fCutPIDUseAntiProtonOnly(kFALSE),
  fCutPIDnSigmaPionMax(3),
  fCutPIDnSigmaKaonMax(3),
  fCutPIDnSigmaProtonMax(3),
  fCutPIDnSigmaTPCRejectElectron(3),
  fCutPIDnSigmaCombinedNoTOFrejection(kFALSE),
  fCutUseBayesPID(kFALSE),
  fCutPIDBayesPionMin(0.9),
  fCutPIDBayesKaonMin(0.9),
  fCutPIDBayesProtonMin(0.9),
  fCutPIDBayesRejectElectron(0.5),
  fCutPIDBayesRejectMuon(0.5),

  // V0s selection
  fCutV0sOnFly(kFALSE),
  fCutV0srejectKinks(kFALSE),
  fCutV0srefitTPC(kFALSE),
  fCutV0sCrossMassRejection(kFALSE),
  fCutV0sCrossMassCutK0s(0.005),
  fCutV0sCrossMassCutLambda(0.010),
  fCutV0sCPAK0sMin(0.),
  fCutV0sCPALambdaMin(0.),
  fCutV0sDCAtoPVMin(0.),
  fCutV0sDCAtoPVMax(0.),
  fCutV0sDCAtoPVzMax(0.),
  fCutV0sDCADaughtersMin(0.),
  fCutV0sDCADaughtersMax(0.),
  fCutV0sDecayRadiusMin(0.),
  fCutV0sDecayRadiusMax(0.),
  fCutV0sDaughterFilterBit(0),
  fCutV0sDaughterPtMin(0.),
  fCutV0sDaughterPtMax(0.),
  fCutV0sDaughterEtaMax(0.),
  fCutV0sMotherEtaMax(0.),
  fCutV0sMotherRapMax(0.),
  fCutV0sArmenterosAlphaK0sMin(0.),
  fCutV0sArmenterosAlphaLambdaMax(0.),
  fCutV0sInvMassK0sMin(0.4),
  fCutV0sInvMassK0sMax(0.6),
  fCutV0sInvMassLambdaMin(1.08),
  fCutV0sInvMassLambdaMax(1.16),
  fCutV0sNumTauK0sMax(0.),
  fCutV0sNumTauLambdaMax(0.),
  fCutV0sProtonNumSigmaMax(0),
  fCutV0sProtonPIDPtMin(0.),
  fCutV0sProtonPIDPtMax(0.),

  // phi selection
  fCutPhiMotherEtaMax(0),
  fCutPhiInvMassMin(0.99),
  fCutPhiInvMassMax(1.07),

  // output lists
  fQAEvents(0x0),
  fQACharged(0x0),
  fQAPID(0x0),
  fQAV0s(0x0),
  fQAPhi(0x0),
  fFlowWeights(0x0),
  fFlowRefs(0x0),
  fFlowCharged(0x0),
  fFlowPID(0x0),
  fFlowPhi(0x0),
  fFlowK0s(0x0),
  fFlowLambda(0x0),

  // flow histograms & profiles
  fh3WeightsRefs(0x0),
  fh3WeightsCharged(0x0),
  fh3WeightsPion(0x0),
  fh3WeightsKaon(0x0),
  fh3WeightsProton(0x0),
  fh3WeightsPhi(0x0),
  fh3WeightsK0s(0x0),
  fh3WeightsLambda(0x0),
  fh3AfterWeightsRefs(0x0),
  fh3AfterWeightsCharged(0x0),
  fh3AfterWeightsPion(0x0),
  fh3AfterWeightsKaon(0x0),
  fh3AfterWeightsProton(0x0),
  fh3AfterWeightsPhi(0x0),
  fh3AfterWeightsK0s(0x0),
  fh3AfterWeightsLambda(0x0),
  fh2WeightRefs(0x0),
  fh2WeightCharged(0x0),
  fh2WeightPion(0x0),
  fh2WeightKaon(0x0),
  fh2WeightProton(0x0),
  fh2WeightK0s(0x0),
  fh2WeightLambda(0x0),
  fh2WeightPhi(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCentrality(0x0),
  fh2EventCentralityNumSelCharged(0x0),
  fhEventCounter(0x0),

  // charged histogram
  fhRefsMult(0x0),
  fhRefsPt(0x0),
  fhRefsEta(0x0),
  fhRefsPhi(0x0),
  fhChargedCounter(0x0),

  // PID histogram
  fhPIDPionMult(0x0),
  fhPIDPionPt(0x0),
  fhPIDPionPhi(0x0),
  fhPIDPionEta(0x0),
  fhPIDPionCharge(0x0),
  fhPIDKaonMult(0x0),
  fhPIDKaonPt(0x0),
  fhPIDKaonPhi(0x0),
  fhPIDKaonEta(0x0),
  fhPIDKaonCharge(0x0),
  fhPIDProtonMult(0x0),
  fhPIDProtonPt(0x0),
  fhPIDProtonPhi(0x0),
  fhPIDProtonEta(0x0),
  fhPIDProtonCharge(0x0),
  fh2PIDPionTPCdEdx(0x0),
  fh2PIDPionTOFbeta(0x0),
  fh2PIDKaonTPCdEdx(0x0),
  fh2PIDKaonTOFbeta(0x0),
  fh2PIDProtonTPCdEdx(0x0),
  fh2PIDProtonTOFbeta(0x0),
  fh2PIDPionTPCnSigmaPion(0x0),
  fh2PIDPionTOFnSigmaPion(0x0),
  fh2PIDPionTPCnSigmaKaon(0x0),
  fh2PIDPionTOFnSigmaKaon(0x0),
  fh2PIDPionTPCnSigmaProton(0x0),
  fh2PIDPionTOFnSigmaProton(0x0),
  fh2PIDPionBayesPion(0x0),
  fh2PIDPionBayesKaon(0x0),
  fh2PIDPionBayesProton(0x0),
  fh2PIDKaonTPCnSigmaPion(0x0),
  fh2PIDKaonTOFnSigmaPion(0x0),
  fh2PIDKaonTPCnSigmaKaon(0x0),
  fh2PIDKaonTOFnSigmaKaon(0x0),
  fh2PIDKaonTPCnSigmaProton(0x0),
  fh2PIDKaonTOFnSigmaProton(0x0),
  fh2PIDKaonBayesPion(0x0),
  fh2PIDKaonBayesKaon(0x0),
  fh2PIDKaonBayesProton(0x0),
  fh2PIDProtonTPCnSigmaPion(0x0),
  fh2PIDProtonTOFnSigmaPion(0x0),
  fh2PIDProtonTPCnSigmaKaon(0x0),
  fh2PIDProtonTOFnSigmaKaon(0x0),
  fh2PIDProtonTPCnSigmaProton(0x0),
  fh2PIDProtonTOFnSigmaProton(0x0),
  fh2PIDProtonBayesPion(0x0),
  fh2PIDProtonBayesKaon(0x0),
  fh2PIDProtonBayesProton(0x0),

  // phi histograms
  fhPhiCounter(0x0),
  fhPhiMult(0x0),
  fhPhiBGMult(0x0),
  fhPhiInvMass(0x0),
  fhPhiBGInvMass(0x0),
  fhPhiCharge(0x0),
  fhPhiBGCharge(0x0),
  fhPhiPt(0x0),
  fhPhiEta(0x0),
  fhPhiPhi(0x0),

  // V0s histogram
  fhV0sCounter(0x0),
  fhV0sCounterK0s(0x0),
  fhV0sCounterLambda(0x0),
  fhV0sCompetingInvMassK0s(0x0),
  fhV0sCompetingInvMassLambda(0x0),

  //Alex
  fEtaCut(0.8),
  fNoClus(70),
  fMinPt(0.2),
  fMaxPt(20.0),
  fNsigCut(3.),
  fNoClusPid(70),
  fNcrFind(0.8),
  fDCADghtPV(0.1),
  fMaxDCADght(0.5),
  fCosPA(0.998),
  fMinRad(5.),
  fMaxRad(100.),
  fArmPodCut(kFALSE),
  fMinPtDght(kFALSE),
  fDoAlexK0sSelection(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow(const char* name) : AliAnalysisTaskSE(name),
  fEventAOD(0x0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fFlowWeightsFile(0x0),
  fInit(kFALSE),
  fIndexSampling(0),
  fIndexCentrality(-1),
  fEventCounter(0),
  fNumEventsAnalyse(50),
  fRunNumber(-1),
  fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
  fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
  fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),
  fPDGMassPhi(TDatabasePDG::Instance()->GetParticle(333)->Mass()),
  fPDGMassK0s(TDatabasePDG::Instance()->GetParticle(310)->Mass()),
  fPDGMassLambda(TDatabasePDG::Instance()->GetParticle(3122)->Mass()),

  // FlowPart containers
  fVectorCharged(0x0),
  fVectorPion(0x0),
  fVectorKaon(0x0),
  fVectorProton(0x0),
  fVectorK0s(0x0),
  fVectorLambda(0x0),
  fVectorPhi(0x0),

  // analysis selection
  fRunMode(kFull),
  fAnalType(kAOD),
  fSampling(kFALSE),
  fFillQA(kTRUE),
  // fNumSamples(10),
  fProcessCharged(kFALSE),
  fProcessPID(kFALSE),
  fProcessV0s(kFALSE),
  fProcessPhi(kFALSE),

  // flow related
  fCutFlowRFPsPtMin(0),
  fCutFlowRFPsPtMax(0),
  fFlowPOIsPtMin(0),
  fFlowPOIsPtMax(20.),
  fCutFlowDoFourCorrelations(kTRUE),
  fFlowFillWeights(kFALSE),
  fFlowCentMin(0),
  fFlowCentMax(150),
  fFlowCentNumBins(150),
  fFlowWeightsPath(),
  fFlowUseWeights(kFALSE),

  // events selection
  fPVtxCutZ(0.),
  fColSystem(kPP),
  fPeriod(kNon),
  fMultEstimator(),
  fTrigger(0),

  // charged tracks selection
  fCutChargedEtaMax(0),
  fCutChargedPtMax(0),
  fCutChargedPtMin(0),
  fCutChargedDCAzMax(0),
  fCutChargedDCAxyMax(0),
  fCutChargedTrackFilterBit(0),
  fCutChargedNumTPCclsMin(0),

  // PID tracks selection
  fCutPIDUseAntiProtonOnly(kFALSE),
  fCutPIDnSigmaPionMax(3),
  fCutPIDnSigmaKaonMax(3),
  fCutPIDnSigmaProtonMax(3),
  fCutPIDnSigmaTPCRejectElectron(3),
  fCutPIDnSigmaCombinedNoTOFrejection(kFALSE),
  fCutUseBayesPID(kFALSE),
  fCutPIDBayesPionMin(0.9),
  fCutPIDBayesKaonMin(0.9),
  fCutPIDBayesProtonMin(0.9),
  fCutPIDBayesRejectElectron(0.5),
  fCutPIDBayesRejectMuon(0.5),

  // V0s selection
  fCutV0sOnFly(kFALSE),
  fCutV0srejectKinks(kFALSE),
  fCutV0srefitTPC(kFALSE),
  fCutV0sCrossMassRejection(kFALSE),
  fCutV0sCrossMassCutK0s(0.005),
  fCutV0sCrossMassCutLambda(0.010),
  fCutV0sCPAK0sMin(0.),
  fCutV0sCPALambdaMin(0.),
  fCutV0sDCAtoPVMin(0.),
  fCutV0sDCAtoPVMax(0.),
  fCutV0sDCAtoPVzMax(0.),
  fCutV0sDCADaughtersMin(0.),
  fCutV0sDCADaughtersMax(0.),
  fCutV0sDecayRadiusMin(0.),
  fCutV0sDecayRadiusMax(0.),
  fCutV0sDaughterFilterBit(0),
  fCutV0sDaughterPtMin(0.),
  fCutV0sDaughterPtMax(0.),
  fCutV0sDaughterEtaMax(0.),
  fCutV0sMotherEtaMax(0.),
  fCutV0sMotherRapMax(0.),
  fCutV0sArmenterosAlphaK0sMin(0.),
  fCutV0sArmenterosAlphaLambdaMax(0.),
  fCutV0sInvMassK0sMin(0.4),
  fCutV0sInvMassK0sMax(0.6),
  fCutV0sInvMassLambdaMin(1.08),
  fCutV0sInvMassLambdaMax(1.16),
  fCutV0sNumTauK0sMax(0.),
  fCutV0sNumTauLambdaMax(0.),
  fCutV0sProtonNumSigmaMax(0),
  fCutV0sProtonPIDPtMin(0.),
  fCutV0sProtonPIDPtMax(0.),

  // phi selection
  fCutPhiMotherEtaMax(0),
  fCutPhiInvMassMin(0.99),
  fCutPhiInvMassMax(1.07),

  // output lists
  fQAEvents(0x0),
  fQACharged(0x0),
  fQAPID(0x0),
  fQAV0s(0x0),
  fQAPhi(0x0),
  fFlowWeights(0x0),
  fFlowRefs(0x0),
  fFlowCharged(0x0),
  fFlowPID(0x0),
  fFlowPhi(0x0),
  fFlowK0s(0x0),
  fFlowLambda(0x0),

  // flow histograms & profiles
  fh3WeightsRefs(0x0),
  fh3WeightsCharged(0x0),
  fh3WeightsPion(0x0),
  fh3WeightsKaon(0x0),
  fh3WeightsProton(0x0),
  fh3WeightsPhi(0x0),
  fh3WeightsK0s(0x0),
  fh3WeightsLambda(0x0),
  fh3AfterWeightsRefs(0x0),
  fh3AfterWeightsCharged(0x0),
  fh3AfterWeightsPion(0x0),
  fh3AfterWeightsKaon(0x0),
  fh3AfterWeightsProton(0x0),
  fh3AfterWeightsPhi(0x0),
  fh3AfterWeightsK0s(0x0),
  fh3AfterWeightsLambda(0x0),
  fh2WeightRefs(0x0),
  fh2WeightCharged(0x0),
  fh2WeightPion(0x0),
  fh2WeightKaon(0x0),
  fh2WeightProton(0x0),
  fh2WeightK0s(0x0),
  fh2WeightLambda(0x0),
  fh2WeightPhi(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCentrality(0x0),
  fh2EventCentralityNumSelCharged(0x0),
  fhEventCounter(0x0),

  // charged histogram
  fhRefsMult(0x0),
  fhRefsPt(0x0),
  fhRefsEta(0x0),
  fhRefsPhi(0x0),
  fhChargedCounter(0x0),

  // PID histogram
  fhPIDPionMult(0x0),
  fhPIDPionPt(0x0),
  fhPIDPionPhi(0x0),
  fhPIDPionEta(0x0),
  fhPIDPionCharge(0x0),
  fhPIDKaonMult(0x0),
  fhPIDKaonPt(0x0),
  fhPIDKaonPhi(0x0),
  fhPIDKaonEta(0x0),
  fhPIDKaonCharge(0x0),
  fhPIDProtonMult(0x0),
  fhPIDProtonPt(0x0),
  fhPIDProtonPhi(0x0),
  fhPIDProtonEta(0x0),
  fhPIDProtonCharge(0x0),
  fh2PIDPionTPCdEdx(0x0),
  fh2PIDPionTOFbeta(0x0),
  fh2PIDKaonTPCdEdx(0x0),
  fh2PIDKaonTOFbeta(0x0),
  fh2PIDProtonTPCdEdx(0x0),
  fh2PIDProtonTOFbeta(0x0),
  fh2PIDPionTPCnSigmaPion(0x0),
  fh2PIDPionTOFnSigmaPion(0x0),
  fh2PIDPionTPCnSigmaKaon(0x0),
  fh2PIDPionTOFnSigmaKaon(0x0),
  fh2PIDPionTPCnSigmaProton(0x0),
  fh2PIDPionTOFnSigmaProton(0x0),
  fh2PIDPionBayesPion(0x0),
  fh2PIDPionBayesKaon(0x0),
  fh2PIDPionBayesProton(0x0),
  fh2PIDKaonTPCnSigmaPion(0x0),
  fh2PIDKaonTOFnSigmaPion(0x0),
  fh2PIDKaonTPCnSigmaKaon(0x0),
  fh2PIDKaonTOFnSigmaKaon(0x0),
  fh2PIDKaonTPCnSigmaProton(0x0),
  fh2PIDKaonTOFnSigmaProton(0x0),
  fh2PIDKaonBayesPion(0x0),
  fh2PIDKaonBayesKaon(0x0),
  fh2PIDKaonBayesProton(0x0),
  fh2PIDProtonTPCnSigmaPion(0x0),
  fh2PIDProtonTOFnSigmaPion(0x0),
  fh2PIDProtonTPCnSigmaKaon(0x0),
  fh2PIDProtonTOFnSigmaKaon(0x0),
  fh2PIDProtonTPCnSigmaProton(0x0),
  fh2PIDProtonTOFnSigmaProton(0x0),
  fh2PIDProtonBayesPion(0x0),
  fh2PIDProtonBayesKaon(0x0),
  fh2PIDProtonBayesProton(0x0),

  // phi histograms
  fhPhiCounter(0x0),
  fhPhiMult(0x0),
  fhPhiBGMult(0x0),
  fhPhiInvMass(0x0),
  fhPhiBGInvMass(0x0),
  fhPhiCharge(0x0),
  fhPhiBGCharge(0x0),
  fhPhiPt(0x0),
  fhPhiEta(0x0),
  fhPhiPhi(0x0),

  // V0s histogram
  fhV0sCounter(0x0),
  fhV0sCounterK0s(0x0),
  fhV0sCounterLambda(0x0),
  fhV0sCompetingInvMassK0s(0x0),
  fhV0sCompetingInvMassLambda(0x0),

  //Alex
  fEtaCut(0.8),
  fNoClus(70),
  fMinPt(0.2),
  fMaxPt(20.0),
  fNsigCut(3.),
  fNoClusPid(70),
  fNcrFind(0.8),
  fDCADghtPV(0.1),
  fMaxDCADght(0.5),
  fCosPA(0.998),
  fMinRad(5.),
  fMaxRad(100.),
  fArmPodCut(kFALSE),
  fMinPtDght(kFALSE),
  fDoAlexK0sSelection(kFALSE)
{
  // Flow vectors
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
      fFlowVecQpos[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecQneg[iHarm][iPower] = TComplex(0,0,kFALSE);

      for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
      {
        fFlowVecPpos[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
        fFlowVecPneg[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
        fFlowVecS[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
      }
    }
  }

  // Flow profiles & histograms
  for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
  {
    for(Short_t iSample(0); iSample < fNumSamples; iSample++)
    {
      fpRefsCor4[iSample][iHarm] = 0x0;
      fp2ChargedCor4[iSample][iHarm] = 0x0;
      fp2PionCor4[iSample][iHarm] = 0x0;
      fp2KaonCor4[iSample][iHarm] = 0x0;
      fp2ProtonCor4[iSample][iHarm] = 0x0;
    }

    fp3V0sCorrK0sCor4[iHarm] = 0x0;
    fp3V0sCorrLambdaCor4[iHarm] = 0x0;
    fp3PhiCorrCor4[iHarm] = 0x0;

    for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
    {
      if(iHarm == 0)
      {
        fh3V0sEntriesK0s[iGap] = 0x0;
        fh3V0sEntriesLambda[iGap] = 0x0;
        fh3PhiEntriesSignal[iGap] = 0x0;
        fh3PhiEntriesBG[iGap] = 0x0;
      }

      // mean Qx,Qy
      fpMeanQxRefsPos[iGap][iHarm] = 0x0;
      fpMeanQxRefsNeg[iGap][iHarm] = 0x0;
      fpMeanQyRefsPos[iGap][iHarm] = 0x0;
      fpMeanQyRefsNeg[iGap][iHarm] = 0x0;

      for(Short_t iSample(0); iSample < fNumSamples; iSample++)
      {
        fpRefsCor2[iSample][iGap][iHarm] = 0x0;
        fp2ChargedCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2ChargedCor2Neg[iSample][iGap][iHarm] = 0x0;
        fp2PionCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2PionCor2Neg[iSample][iGap][iHarm] = 0x0;
        fp2KaonCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2KaonCor2Neg[iSample][iGap][iHarm] = 0x0;
        fp2ProtonCor2Pos[iSample][iGap][iHarm] = 0x0;
        fp2ProtonCor2Neg[iSample][iGap][iHarm] = 0x0;
      }

      fp3V0sCorrK0sCor2Pos[iGap][iHarm] = 0x0;
      fp3V0sCorrK0sCor2Neg[iGap][iHarm] = 0x0;
      fp3V0sCorrLambdaCor2Pos[iGap][iHarm] = 0x0;
      fp3V0sCorrLambdaCor2Neg[iGap][iHarm] = 0x0;
      fp3PhiCorrCor2Pos[iGap][iHarm] = 0x0;
      fp3PhiCorrCor2Neg[iGap][iHarm] = 0x0;
    }
  }

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

    // PID
    fhQAPIDTPCstatus[iQA] = 0x0;
    fhQAPIDTOFstatus[iQA] = 0x0;
    fhQAPIDTPCdEdx[iQA] = 0x0;
    fhQAPIDTOFbeta[iQA] = 0x0;
    fh3QAPIDnSigmaBayesElectron[iQA] = 0x0;
    fh3QAPIDnSigmaBayesMuon[iQA] = 0x0;
    fh3QAPIDnSigmaBayesPion[iQA] = 0x0;
    fh3QAPIDnSigmaBayesKaon[iQA] = 0x0;
    fh3QAPIDnSigmaBayesProton[iQA] = 0x0;

    // V0s
    fhQAV0sMultK0s[iQA] = 0x0;
    fhQAV0sMultLambda[iQA] = 0x0;
  	fhQAV0sRecoMethod[iQA] = 0x0;
		fhQAV0sDCAtoPV[iQA] = 0x0;
		fhQAV0sDCADaughters[iQA] = 0x0;
		fhQAV0sDecayRadius[iQA] = 0x0;
    fhQAV0sDaughterTPCRefit[iQA] = 0x0;
    fhQAV0sDaughterKinks[iQA] = 0x0;
    fhQAV0sDaughterPt[iQA] = 0x0;
		fhQAV0sDaughterPhi[iQA] = 0x0;
		fhQAV0sDaughterEta[iQA] = 0x0;
    fhQAV0sDaughterCharge[iQA] = 0x0;
    fhQAV0sDaughterTPCdEdxK0s[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionK0s[iQA] = 0x0;
    fhQAV0sDaughterTPCstatus[iQA] = 0x0;
    fhQAV0sDaughterTOFstatus[iQA] = 0x0;
    fhQAV0sDaughterTPCdEdxLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaProtonLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionALambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaProtonALambda[iQA] = 0x0;
    fhQAV0sMotherPt[iQA] = 0x0;
		fhQAV0sMotherPhi[iQA] = 0x0;
		fhQAV0sMotherEta[iQA] = 0x0;
    fhQAV0sMotherCharge[iQA] = 0x0;
		fhQAV0sMotherRapK0s[iQA] = 0x0;
		fhQAV0sMotherRapLambda[iQA] = 0x0;
    fhQAV0sInvMassK0s[iQA] = 0x0;
    fhQAV0sInvMassLambda[iQA] = 0x0;
		fhQAV0sCPAK0s[iQA] = 0x0;
		fhQAV0sCPALambda[iQA] = 0x0;
		fhQAV0sNumTauK0s[iQA] = 0x0;
		fhQAV0sNumTauLambda[iQA] = 0x0;
		fhQAV0sArmenterosK0s[iQA] = 0x0;
		fhQAV0sArmenterosLambda[iQA] = 0x0;
		fhQAV0sArmenterosALambda[iQA] = 0x0;
  }

  // defining input/output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  DefineOutput(8, TList::Class());
  DefineOutput(9, TList::Class());
  DefineOutput(10, TList::Class());
  DefineOutput(11, TList::Class());
  DefineOutput(12, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::~AliAnalysisTaskUniFlow()
{
  // destructor
  // if(fPIDCombined)
  // {
  //   delete fPIDCombined;
  // }

  // deleting FlowPart vectors (containers)
  if(fVectorCharged) delete fVectorCharged;
  if(fVectorPion) delete fVectorPion;
  if(fVectorKaon) delete fVectorKaon;
  if(fVectorProton) delete fVectorProton;
  if(fVectorK0s) delete fVectorK0s;
  if(fVectorLambda) delete fVectorLambda;
  if(fVectorPhi) delete fVectorPhi;

  // deleting output lists
  if(fFlowWeights) delete fFlowWeights;
  if(fFlowRefs) delete fFlowRefs;
  if(fFlowCharged) delete fFlowCharged;
  if(fFlowPID) delete fFlowPID;
  if(fFlowPhi) delete fFlowPhi;
  if(fFlowK0s) delete fFlowK0s;
  if(fFlowLambda) delete fFlowLambda;

  if(fQAEvents) delete fQAEvents;
  if(fQACharged) delete fQACharged;
  if(fQAPID) delete fQAPID;
  if(fQAPhi) delete fQAPhi;
  if(fQAV0s) delete fQAV0s;
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

  // creating output lists
  fFlowRefs = new TList();
  fFlowRefs->SetOwner(kTRUE);
  fFlowRefs->SetName("fFlowRefs");
  fFlowCharged = new TList();
  fFlowCharged->SetOwner(kTRUE);
  fFlowCharged->SetName("fFlowCharged");
  fFlowPID = new TList();
  fFlowPID->SetOwner(kTRUE);
  fFlowPID->SetName("fFlowPID");
  fFlowPhi = new TList();
  fFlowPhi->SetOwner(kTRUE);
  fFlowPhi->SetName("fFlowPhi");
  fFlowK0s = new TList();
  fFlowK0s->SetOwner(kTRUE);
  fFlowK0s->SetName("fFlowK0s");
  fFlowLambda = new TList();
  fFlowLambda->SetOwner(kTRUE);
  fFlowLambda->SetName("fFlowLambda");
  fFlowWeights = new TList();
  fFlowWeights->SetOwner(kTRUE);
  fFlowWeights->SetName("fFlowWeights");

  fQAEvents = new TList();
  fQAEvents->SetOwner(kTRUE);
  fQACharged = new TList();
  fQACharged->SetOwner(kTRUE);
  fQAPID = new TList();
  fQAPID->SetOwner(kTRUE);
  fQAPhi = new TList();
  fQAPhi->SetOwner(kTRUE);
  fQAV0s = new TList();
  fQAV0s->SetOwner(kTRUE);

  // creating particle vectors
  fVectorCharged = new std::vector<FlowPart>;
  fVectorPion = new std::vector<FlowPart>;
  fVectorKaon = new std::vector<FlowPart>;
  fVectorProton = new std::vector<FlowPart>;
  fVectorPhi = new std::vector<FlowPart>;
  fVectorK0s = new std::vector<FlowPart>;
  fVectorLambda = new std::vector<FlowPart>;

  fVectorCharged->reserve(300);
  if(fProcessPID) { fVectorPion->reserve(200); fVectorKaon->reserve(100); fVectorProton->reserve(100); }
  if(fProcessPhi) { fVectorPhi->reserve(200); }
  if(fProcessV0s) { fVectorK0s->reserve(100); fVectorLambda->reserve(100); }

  // creating histograms
    // event histogram
    fhEventSampling = new TH2D("fhEventSampling","Event sampling; centrality/multiplicity; sample index", fFlowCentNumBins,0,fFlowCentNumBins, fNumSamples,0,fNumSamples);
    fQAEvents->Add(fhEventSampling);
    fhEventCentrality = new TH1D("fhEventCentrality",Form("Event centrality (%s); centrality/multiplicity",fMultEstimator.Data()), fFlowCentNumBins,0,fFlowCentNumBins);
    fQAEvents->Add(fhEventCentrality);
    fh2EventCentralityNumSelCharged = new TH2D("fh2EventCentralityNumSelCharged",Form("Event centrality (%s) vs. N^{sel}_{ch}; N^{sel}_{ch}; centrality/multiplicity",fMultEstimator.Data()), 150,0,150, fFlowCentNumBins,0,fFlowCentNumBins);
    fQAEvents->Add(fh2EventCentralityNumSelCharged);

    const Short_t iEventCounterBins = 7;
    TString sEventCounterLabel[iEventCounterBins] = {"Input","Physics selection OK","PV OK","SPD Vtx OK","Pileup MV OK","PV #it{z} OK","Selected"};
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    for(Short_t i(0); i < iEventCounterBins; i++) fhEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    fQAEvents->Add(fhEventCounter);

    // flow histograms & profiles
    // weights
    if(fFlowFillWeights || fRunMode == kFillWeights)
    {
      // for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
      //   for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
      //   {
      //     fpMeanQxRefsPos[iGap][iHarm] = new TProfile(Form("fpMeanQxRefs_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]),Form("<<Qx>>: Refs | Gap %g | n=%d | POS; multiplicity/centrality; <Qx>",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,0,fFlowCentNumBins);
      //     fpMeanQxRefsPos[iGap][iHarm]->Sumw2();
      //     fFlowWeights->Add(fpMeanQxRefsPos[iGap][iHarm]);
      //
      //     fpMeanQxRefsNeg[iGap][iHarm] = new TProfile(Form("fpMeanQxRefs_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]),Form("<<Qx>>: Refs | Gap %g | n=%d | NEG; multiplicity/centrality; <Qx>",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,0,fFlowCentNumBins);
      //     fpMeanQxRefsNeg[iGap][iHarm]->Sumw2();
      //     fFlowWeights->Add(fpMeanQxRefsNeg[iGap][iHarm]);
      //
      //     fpMeanQyRefsPos[iGap][iHarm] = new TProfile(Form("fpMeanQyRefs_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]),Form("<<Qy>>: Refs | Gap %g | n=%d | POS; multiplicity/centrality; <Qy>",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,0,fFlowCentNumBins);
      //     fpMeanQyRefsPos[iGap][iHarm]->Sumw2();
      //     fFlowWeights->Add(fpMeanQyRefsPos[iGap][iHarm]);
      //
      //     fpMeanQyRefsNeg[iGap][iHarm] = new TProfile(Form("fpMeanQyRefs_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]),Form("<<Qy>>: Refs | Gap %g | n=%d | NEG; multiplicity/centrality; <Qy>",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,0,fFlowCentNumBins);
      //     fpMeanQyRefsNeg[iGap][iHarm]->Sumw2();
      //     fFlowWeights->Add(fpMeanQyRefsNeg[iGap][iHarm]);
      //   }

      fh3WeightsRefs = new TH3D("fh3WeightsRefs","Weights: Refs; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsRefs->Sumw2();
      fFlowWeights->Add(fh3WeightsRefs);
      fh3WeightsCharged = new TH3D("fh3WeightsCharged","Weights: Charged; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsCharged->Sumw2();
      fFlowWeights->Add(fh3WeightsCharged);
      fh3WeightsPion = new TH3D("fh3WeightsPion","Weights: #pi; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsPion->Sumw2();
      fFlowWeights->Add(fh3WeightsPion);
      fh3WeightsKaon = new TH3D("fh3WeightsKaon","Weights: K; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsKaon->Sumw2();
      fFlowWeights->Add(fh3WeightsKaon);
      fh3WeightsProton = new TH3D("fh3WeightsProton","Weights: p; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsProton->Sumw2();
      fFlowWeights->Add(fh3WeightsProton);
      fh3WeightsPhi = new TH3D("fh3WeightsPhi","Weights: #phi; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,-TMath::Pi(),TMath::Pi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsPhi->Sumw2();
      fFlowWeights->Add(fh3WeightsPhi);
      fh3WeightsK0s = new TH3D("fh3WeightsK0s","Weights: K^{0}_{S}; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsK0s->Sumw2();
      fFlowWeights->Add(fh3WeightsK0s);
      fh3WeightsLambda = new TH3D("fh3WeightsLambda","Weights: #Lambda; #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3WeightsLambda->Sumw2();
      fFlowWeights->Add(fh3WeightsLambda);
    }

    if(fFlowUseWeights)
    {
      fh3AfterWeightsRefs = new TH3D("fh3AfterWeightsRefs","Weights: Refs; #varphi (After); #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsRefs->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsRefs);
      fh3AfterWeightsCharged = new TH3D("fh3AfterWeightsCharged","Weights: Charged (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsCharged->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsCharged);
      fh3AfterWeightsPion = new TH3D("fh3AfterWeightsPion","Weights: #pi (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsPion->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsPion);
      fh3AfterWeightsKaon = new TH3D("fh3AfterWeightsKaon","Weights: K (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsKaon->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsKaon);
      fh3AfterWeightsProton = new TH3D("fh3AfterWeightsProton","Weights: p (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsProton->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsProton);
      fh3AfterWeightsPhi = new TH3D("fh3AfterWeightsPhi","Weights: #phi (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,-TMath::Pi(),TMath::Pi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsPhi->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsPhi);
      fh3AfterWeightsK0s = new TH3D("fh3AfterWeightsK0s","Weights: K^{0}_{S} (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsK0s->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsK0s);
      fh3AfterWeightsLambda = new TH3D("fh3AfterWeightsLambda","Weights: #Lambda (After); #varphi; #eta; #it{p}_{T} (GeV/#it{c})", 100,0,TMath::TwoPi(), 151,-1.5,1.5, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
      fh3AfterWeightsLambda->Sumw2();
      fFlowWeights->Add(fh3AfterWeightsLambda);
    }

    // candidate distribution for flow-mass method
    for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
    {
      if(fProcessPhi)
      {
        fh3PhiEntriesSignal[iGap] = new TH3D(Form("fh3PhiEntriesSignal_gap%02.2g",10*fEtaGap[iGap]), Form("#phi: Distribution (Gap %g); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
        fh3PhiEntriesSignal[iGap]->Sumw2();
        fFlowPhi->Add(fh3PhiEntriesSignal[iGap]);
        fh3PhiEntriesBG[iGap] = new TH3D(Form("fh3PhiEntriesBG_gap%02.2g",10*fEtaGap[iGap]), Form("#phi (BG): Distribution (Gap %g); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
        fh3PhiEntriesBG[iGap]->Sumw2();
        fFlowPhi->Add(fh3PhiEntriesBG[iGap]);
      }

      if(fProcessV0s)
      {
        fh3V0sEntriesK0s[iGap] = new TH3D(Form("fh3V0sEntriesK0s_gap%02.2g",10*fEtaGap[iGap]), Form("K_{S}^{0}: Distribution (Gap %g); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
        fh3V0sEntriesK0s[iGap]->Sumw2();
        fFlowK0s->Add(fh3V0sEntriesK0s[iGap]);
        fh3V0sEntriesLambda[iGap] = new TH3D(Form("fh3V0sEntriesLambda_gap%02.2g",10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: Distribution (Gap %g); centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
        fh3V0sEntriesLambda[iGap]->Sumw2();
        fFlowLambda->Add(fh3V0sEntriesLambda[iGap]);
      }
    }

    // correlations
    for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
    {
      for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
      {
        for(Short_t iSample(0); iSample < fNumSamples; iSample++)
        {
          fpRefsCor2[iSample][iGap][iHarm] = new TProfile(Form("fpRefs_<2>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<2>> | Gap %g | n=%d | sample %d ; centrality/multiplicity;",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax);
          fpRefsCor2[iSample][iGap][iHarm]->Sumw2(kTRUE);
          fFlowRefs->Add(fpRefsCor2[iSample][iGap][iHarm]);

          if(fCutFlowDoFourCorrelations && iGap == 0)
          {
            fpRefsCor4[iSample][iHarm] = new TProfile(Form("fpRefs_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<4>> | Gap %g | n=%d | sample %d ; centrality/multiplicity;",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax);
            fpRefsCor4[iSample][iHarm]->Sumw2(kTRUE);
            fFlowRefs->Add(fpRefsCor4[iSample][iHarm]);
          }

          if(fProcessCharged)
          {
            fp2ChargedCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Charged_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<2'>> | Gap %g | n=%d | sample %d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
            fp2ChargedCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
            fFlowCharged->Add(fp2ChargedCor2Pos[iSample][iGap][iHarm]);

            if(fEtaGap[iGap] != -1.)
            {
              fp2ChargedCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Charged_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ChargedCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowCharged->Add(fp2ChargedCor2Neg[iSample][iGap][iHarm]);
            }

            if(fCutFlowDoFourCorrelations && iGap == 0)
            {
              fp2ChargedCor4[iSample][iHarm] = new TProfile2D(Form("fp2Charged_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ChargedCor4[iSample][iHarm]->Sumw2(kTRUE);
              fFlowCharged->Add(fp2ChargedCor4[iSample][iHarm]);
            }
          }

          if(fProcessPID)
          {
            fp2PionCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Pion_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<2'>> | Gap %g | n=%d | sample %d  | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
            fp2PionCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
            fFlowPID->Add(fp2PionCor2Pos[iSample][iGap][iHarm]);

            fp2KaonCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Kaon_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<2'>> | Gap %g | n=%d | sample %d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
            fp2KaonCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
            fFlowPID->Add(fp2KaonCor2Pos[iSample][iGap][iHarm]);

            fp2ProtonCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Proton_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<2'>> | Gap %g | n=%d | sample %d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
            fp2ProtonCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
            fFlowPID->Add(fp2ProtonCor2Pos[iSample][iGap][iHarm]);

            if(fEtaGap[iGap] != -1.)
            {
              fp2PionCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Pion_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2PionCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2PionCor2Neg[iSample][iGap][iHarm]);

              fp2KaonCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Kaon_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2KaonCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2KaonCor2Neg[iSample][iGap][iHarm]);

              fp2ProtonCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Proton_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<2'>> | Gap %g | n=%d | sample %d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ProtonCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2ProtonCor2Neg[iSample][iGap][iHarm]);
            }

            if(fCutFlowDoFourCorrelations && iGap == 0)
            {
              fp2PionCor4[iSample][iHarm] = new TProfile2D(Form("fp2Pion_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2PionCor4[iSample][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2PionCor4[iSample][iHarm]);
              fp2KaonCor4[iSample][iHarm] = new TProfile2D(Form("fp2Kaon_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2KaonCor4[iSample][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2KaonCor4[iSample][iHarm]);
              fp2ProtonCor4[iSample][iHarm] = new TProfile2D(Form("fp2Proton_<4>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<4'>> | Gap %g | n=%d | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
              fp2ProtonCor4[iSample][iHarm]->Sumw2(kTRUE);
              fFlowPID->Add(fp2ProtonCor4[iSample][iHarm]);
            }
          }
        }

        if(fProcessPhi)
        {
          fp3PhiCorrCor2Pos[iGap][iHarm] = new TProfile3D(Form("fp3PhiCorr_<2>_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#phi: <<2'>> | Gap %g | n=%d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
          fp3PhiCorrCor2Pos[iGap][iHarm]->Sumw2();
          fFlowPhi->Add(fp3PhiCorrCor2Pos[iGap][iHarm]);

          if(fEtaGap[iGap] != -1.)
          {
            fp3PhiCorrCor2Neg[iGap][iHarm] = new TProfile3D(Form("fp3PhiCorr_<2>_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#phi: <<2'>> | Gap %g | n=%d  | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
            fp3PhiCorrCor2Neg[iGap][iHarm]->Sumw2();
            fFlowPhi->Add(fp3PhiCorrCor2Neg[iGap][iHarm]);
          }

          if(fCutFlowDoFourCorrelations && iGap == 0)
          {
            fp3PhiCorrCor4[iHarm] = new TProfile3D(Form("fp3PhiCorr_<4>_harm%d_gap%02.2g",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#phi: <<4'>> | Gap %g | n=%d; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
            fp3PhiCorrCor4[iHarm]->Sumw2();
            fFlowPhi->Add(fp3PhiCorrCor4[iHarm]);
          }
        }

        if(fProcessV0s)
        {
          fp3V0sCorrK0sCor2Pos[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrK0s_<2>_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("K_{S}^{0}: <<2'>> | Gap %g | n=%d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
          fp3V0sCorrK0sCor2Pos[iGap][iHarm]->Sumw2();
          fFlowK0s->Add(fp3V0sCorrK0sCor2Pos[iGap][iHarm]);
          fp3V0sCorrLambdaCor2Pos[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrLambda_<2>_harm%d_gap%02.2g_Pos",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: <<2'>> | Gap %g | n=%d | POIs pos; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
          fp3V0sCorrLambdaCor2Pos[iGap][iHarm]->Sumw2();
          fFlowLambda->Add(fp3V0sCorrLambdaCor2Pos[iGap][iHarm]);

          if(fEtaGap[iGap] != -1.)
          {
            fp3V0sCorrK0sCor2Neg[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrK0s_<2>_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("K_{S}^{0}: <<2'>> | Gap %g | n=%d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
            fp3V0sCorrK0sCor2Neg[iGap][iHarm]->Sumw2();
            fFlowK0s->Add(fp3V0sCorrK0sCor2Neg[iGap][iHarm]);

            fp3V0sCorrLambdaCor2Neg[iGap][iHarm] = new TProfile3D(Form("fp3V0sCorrLambda_<2>_harm%d_gap%02.2g_Neg",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: <<2'>> | Gap %g | n=%d | POIs neg; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
            fp3V0sCorrLambdaCor2Neg[iGap][iHarm]->Sumw2();
            fFlowLambda->Add(fp3V0sCorrLambdaCor2Neg[iGap][iHarm]);
          }

          if(fCutFlowDoFourCorrelations && iGap == 0)
          {
            fp3V0sCorrK0sCor4[iHarm] = new TProfile3D(Form("fp3V0sCorrK0s_<4>_harm%d_gap%02.2g",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("K_{S}^{0}: <<4'>> | Gap %g | n=%d; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
            fp3V0sCorrK0sCor4[iHarm]->Sumw2();
            fFlowK0s->Add(fp3V0sCorrK0sCor4[iHarm]);
            fp3V0sCorrLambdaCor4[iHarm] = new TProfile3D(Form("fp3V0sCorrLambda_<4>_harm%d_gap%02.2g",fHarmonics[iHarm],10*fEtaGap[iGap]), Form("#Lambda/#bar{#Lambda}: <<4'>> | Gap %g | n=%d; centrality/multiplicity; #it{p}_{T} (GeV/c); #it{m}_{inv} (GeV/#it{c}^{2})",fEtaGap[iGap],fHarmonics[iHarm]), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
            fp3V0sCorrLambdaCor4[iHarm]->Sumw2();
            fFlowLambda->Add(fp3V0sCorrLambdaCor4[iHarm]);
          }
        }
      }
    }

    // charged (tracks) histograms
    fhRefsMult = new TH1D("fhRefsMult","RFPs: Multiplicity; multiplicity", 1000,0,1000);
    fQACharged->Add(fhRefsMult);
    fhRefsPt = new TH1D("fhRefsPt","RFPs: #it{p}_{T};  #it{p}_{T} (GeV/#it{c})", 300,0,30);
    fQACharged->Add(fhRefsPt);
    fhRefsEta = new TH1D("fhRefsEta","RFPs: #eta; #eta", 151,-1.5,1.5);
    fQACharged->Add(fhRefsEta);
    fhRefsPhi = new TH1D("fhRefsPhi","RFPs: #varphi; #varphi", 100,0,TMath::TwoPi());
    fQACharged->Add(fhRefsPhi);

    if(fProcessCharged)
    {
      TString sChargedCounterLabel[] = {"Input","FB","#TPC-Cls","DCA-z","DCA-xy","Eta","Selected"};
      const Short_t iNBinsChargedCounter = sizeof(sChargedCounterLabel)/sizeof(sChargedCounterLabel[0]);
      fhChargedCounter = new TH1D("fhChargedCounter","Charged tracks: Counter",iNBinsChargedCounter,0,iNBinsChargedCounter);
      for(Short_t i(0); i < iNBinsChargedCounter; i++) fhChargedCounter->GetXaxis()->SetBinLabel(i+1, sChargedCounterLabel[i].Data() );
      fQACharged->Add(fhChargedCounter);
    } // endif {fProcessCharged}

    // PID tracks histograms
    if(fProcessPID || fProcessPhi)
    {
      fhPIDPionMult = new TH1D("fhPIDPionMult","PID: #pi: Multiplicity; multiplicity", 200,0,200);
      fQAPID->Add(fhPIDPionMult);
      fhPIDKaonMult = new TH1D("fhPIDKaonMult","PID: K: Multiplicity; multiplicity", 100,0,100);
      fQAPID->Add(fhPIDKaonMult);
      fhPIDProtonMult = new TH1D("fhPIDProtonMult","PID: p: Multiplicity; multiplicity", 100,0,100);
      fQAPID->Add(fhPIDProtonMult);

      if(fFillQA)
      {
        fhPIDPionPt = new TH1D("fhPIDPionPt","PID: #pi: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDPionPt);
        fhPIDPionPhi = new TH1D("fhPIDPionPhi","PID: #pi: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDPionPhi);
        fhPIDPionEta = new TH1D("fhPIDPionEta","PID: #pi: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDPionEta);
        fhPIDPionCharge = new TH1D("fhPIDPionCharge","PID: #pi: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDPionCharge);
        fhPIDKaonPt = new TH1D("fhPIDKaonPt","PID: K: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDKaonPt);
        fhPIDKaonPhi = new TH1D("fhPIDKaonPhi","PID: K: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDKaonPhi);
        fhPIDKaonEta = new TH1D("fhPIDKaonEta","PID: K: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDKaonEta);
        fhPIDKaonCharge = new TH1D("fhPIDKaonCharge","PID: K: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDKaonCharge);
        fhPIDProtonPt = new TH1D("fhPIDProtonPt","PID: p: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDProtonPt);
        fhPIDProtonPhi = new TH1D("fhPIDProtonPhi","PID: p: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDProtonPhi);
        fhPIDProtonEta = new TH1D("fhPIDProtonEta","PID: p: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDProtonEta);
        fhPIDProtonCharge = new TH1D("fhPIDProtonCharge","PID: p: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDProtonCharge);
        fh2PIDPionTPCdEdx = new TH2D("fh2PIDPionTPCdEdx","PID: #pi: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDPionTPCdEdx);
        fh2PIDPionTOFbeta = new TH2D("fh2PIDPionTOFbeta","PID: #pi: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDPionTOFbeta);
        fh2PIDKaonTPCdEdx = new TH2D("fh2PIDKaonTPCdEdx","PID: K: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDKaonTPCdEdx);
        fh2PIDKaonTOFbeta = new TH2D("fh2PIDKaonTOFbeta","PID: K: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDKaonTOFbeta);
        fh2PIDProtonTPCdEdx = new TH2D("fh2PIDProtonTPCdEdx","PID: p: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDProtonTPCdEdx);
        fh2PIDProtonTOFbeta = new TH2D("fh2PIDProtonTOFbeta","PID: p: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDProtonTOFbeta);
        fh2PIDPionTPCnSigmaPion = new TH2D("fh2PIDPionTPCnSigmaPion","PID: #pi: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaPion);
        fh2PIDPionTOFnSigmaPion = new TH2D("fh2PIDPionTOFnSigmaPion","PID: #pi: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaPion);
        fh2PIDPionTPCnSigmaKaon = new TH2D("fh2PIDPionTPCnSigmaKaon","PID: #pi: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaKaon);
        fh2PIDPionTOFnSigmaKaon = new TH2D("fh2PIDPionTOFnSigmaKaon","PID: #pi: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaKaon);
        fh2PIDPionTPCnSigmaProton = new TH2D("fh2PIDPionTPCnSigmaProton","PID: #pi: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaProton);
        fh2PIDPionTOFnSigmaProton = new TH2D("fh2PIDPionTOFnSigmaProton","PID: #pi: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaProton);
        fh2PIDPionBayesPion = new TH2D("fh2PIDPionBayesPion","PID: #pi: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesPion);
        fh2PIDPionBayesKaon = new TH2D("fh2PIDPionBayesKaon","PID: #pi: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesKaon);
        fh2PIDPionBayesProton = new TH2D("fh2PIDPionBayesProton","PID: #pi: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesProton);
        fh2PIDKaonTPCnSigmaPion = new TH2D("fh2PIDKaonTPCnSigmaPion","PID: K: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaPion);
        fh2PIDKaonTOFnSigmaPion = new TH2D("fh2PIDKaonTOFnSigmaPion","PID: K: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaPion);
        fh2PIDKaonTPCnSigmaKaon = new TH2D("fh2PIDKaonTPCnSigmaKaon","PID: K: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaKaon);
        fh2PIDKaonTOFnSigmaKaon = new TH2D("fh2PIDKaonTOFnSigmaKaon","PID: K: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaKaon);
        fh2PIDKaonTPCnSigmaProton = new TH2D("fh2PIDKaonTPCnSigmaProton","PID: K: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaProton);
        fh2PIDKaonTOFnSigmaProton = new TH2D("fh2PIDKaonTOFnSigmaProton","PID: K: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaProton);
        fh2PIDKaonBayesPion = new TH2D("fh2PIDKaonBayesPion","PID: K: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesPion);
        fh2PIDKaonBayesKaon = new TH2D("fh2PIDKaonBayesKaon","PID: K: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesKaon);
        fh2PIDKaonBayesProton = new TH2D("fh2PIDKaonBayesProton","PID: K: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesProton);
        fh2PIDProtonTPCnSigmaPion = new TH2D("fh2PIDProtonTPCnSigmaPion","PID: p: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaPion);
        fh2PIDProtonTOFnSigmaPion = new TH2D("fh2PIDProtonTOFnSigmaPion","PID: p: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaPion);
        fh2PIDProtonTPCnSigmaKaon = new TH2D("fh2PIDProtonTPCnSigmaKaon","PID: p: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaKaon);
        fh2PIDProtonTOFnSigmaKaon = new TH2D("fh2PIDProtonTOFnSigmaKaon","PID: p: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaKaon);
        fh2PIDProtonTPCnSigmaProton = new TH2D("fh2PIDProtonTPCnSigmaProton","PID: p: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaProton);
        fh2PIDProtonTOFnSigmaProton = new TH2D("fh2PIDProtonTOFnSigmaProton","PID: p: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaProton);
        fh2PIDProtonBayesPion = new TH2D("fh2PIDProtonBayesPion","PID: p: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesPion);
        fh2PIDProtonBayesKaon = new TH2D("fh2PIDProtonBayesKaon","PID: p: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesKaon);
        fh2PIDProtonBayesProton = new TH2D("fh2PIDProtonBayesProton","PID: p: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesProton);

      }

    } //endif {fProcessPID}

    if(fProcessPhi)
    {
      TString sPhiCounterLabel[] = {"Input","InvMass","Pt","Eta","Selected","Unlike-sign","BG"};
      const Short_t iNBinsPhiCounter = sizeof(sPhiCounterLabel)/sizeof(sPhiCounterLabel[0]);
      fhPhiCounter = new TH1D("fhPhiCounter","#phi: Counter",iNBinsPhiCounter,0,iNBinsPhiCounter);
      for(Short_t i(0); i < iNBinsPhiCounter; i++) fhPhiCounter->GetXaxis()->SetBinLabel(i+1, sPhiCounterLabel[i].Data() );
      fQAPhi->Add(fhPhiCounter);

      fhPhiMult = new TH1D("fhPhiMult","#phi: Multiplicity; Multiplicity", 150,0,150);
      fQAPhi->Add(fhPhiMult);
      fhPhiBGMult = new TH1D("fhPhiBGMult","#phi (BG): Multiplicity; Multiplicity", 150,0,150);
      fQAPhi->Add(fhPhiBGMult);
      fhPhiInvMass = new TH1D("fhPhiInvMass","#phi: InvMass; #it{m}_{inv} (GeV/#it{c}^{2})", 90,fCutPhiInvMassMin,fCutPhiInvMassMax);
      fQAPhi->Add(fhPhiInvMass);
      fhPhiBGInvMass = new TH1D("fhPhiBGInvMass","#phi (BG): InvMass; #it{m}_{inv} (GeV/#it{c}^{2})", 90,fCutPhiInvMassMin,fCutPhiInvMassMax);
      fQAPhi->Add(fhPhiBGInvMass);
      fhPhiCharge = new TH1D("fhPhiCharge","#phi: charge; charge", 5,-2.5,2.5);
      fQAPhi->Add(fhPhiCharge);
      fhPhiBGCharge = new TH1D("fhPhiBGCharge","#phi (BG): charge; charge", 5,-2.5,2.5);
      fQAPhi->Add(fhPhiBGCharge);
      fhPhiPt = new TH1D("fhPhiPt","#phi: #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 150,0.,30.);
      fQAPhi->Add(fhPhiPt);
      fhPhiEta = new TH1D("fhPhiEta","#phi: #eta; #eta", 151,-1.5,1.5);
      fQAPhi->Add(fhPhiEta);
      fhPhiPhi = new TH1D("fhPhiPhi","#phi: #varphi; #varphi", 100,-TMath::Pi(),TMath::Pi());
      fQAPhi->Add(fhPhiPhi);
    } //endif {fProcessPhi}

    // V0 candidates histograms
    if(fProcessV0s)
    {
      TString sV0sCounterLabel[] = {"Input","Daughters OK","Charge","Reconstruction method","TPC refit","Kinks","DCA to PV","Daughters DCA","Decay radius","Acceptance","Common passed","K^{0}_{S}","#Lambda/#bar{#Lambda}","K^{0}_{S} && #Lambda/#bar{#Lambda}"};
      const Short_t iNBinsV0sCounter = sizeof(sV0sCounterLabel)/sizeof(sV0sCounterLabel[0]);
      fhV0sCounter = new TH1D("fhV0sCounter","V^{0}: Counter",iNBinsV0sCounter,0,iNBinsV0sCounter);
      for(Short_t i(0); i < iNBinsV0sCounter; i++) fhV0sCounter->GetXaxis()->SetBinLabel(i+1, sV0sCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounter);

      TString sV0sK0sCounterLabel[] = {"Input","CPA","c#tau","Armenteros-Podolanski","#it{y}","InvMass","Competing InvMass","Selected"};
      const Short_t iNBinsV0sK0sCounter = sizeof(sV0sK0sCounterLabel)/sizeof(sV0sK0sCounterLabel[0]);
      fhV0sCounterK0s = new TH1D("fhV0sCounterK0s","V^{0}: K^{0}_{S} Counter",iNBinsV0sK0sCounter,0,iNBinsV0sK0sCounter);
      for(Short_t i(0); i < iNBinsV0sK0sCounter; i++) fhV0sCounterK0s->GetXaxis()->SetBinLabel(i+1, sV0sK0sCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounterK0s);

      TString sV0sLambdaCounterLabel[] = {"Input","CPA","c#tau","Armenteros-Podolanski","#it{y}","InvMass","Competing InvMass","Proton PID","Selected","only #Lambda","only #bar{#Lambda}","#Lambda && #bar{#Lambda}"};
      const Short_t iNBinsV0sLambdaCounter = sizeof(sV0sLambdaCounterLabel)/sizeof(sV0sLambdaCounterLabel[0]);
      fhV0sCounterLambda = new TH1D("fhV0sCounterLambda","V^{0}: #Lambda/#bar{#Lambda} Counter",iNBinsV0sLambdaCounter,0,iNBinsV0sLambdaCounter);
      for(Short_t i(0); i < iNBinsV0sLambdaCounter; i++) fhV0sCounterLambda->GetXaxis()->SetBinLabel(i+1, sV0sLambdaCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounterLambda);

      fhV0sCompetingInvMassK0s = new TH2D("fhV0sCompetingInvMassK0s","V^{0}: K^{0}_{S}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassK0s);
      fhV0sCompetingInvMassLambda = new TH2D("fhV0sCompetingInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: Competing InvMass rejection; #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2}); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2})", 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax, 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
      fQAV0s->Add(fhV0sCompetingInvMassLambda);
    } // endif {fProcessV0s}

    const Short_t iNBinsPIDstatus = 4;
    TString sPIDstatus[iNBinsPIDstatus] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};
    const Short_t iNFilterMapBinBins = 32;

    // QA histograms
    if(fFillQA)
    {
      TString sQAindex[fiNumIndexQA] = {"Before", "After"};
      for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
      {
        // EVENTs QA histograms
        fhQAEventsPVz[iQA] = new TH1D(Form("fhQAEventsPVz_%s",sQAindex[iQA].Data()), "QA Events: PV-#it{z}", 101,-50,50);
        fQAEvents->Add(fhQAEventsPVz[iQA]);
        fhQAEventsNumContrPV[iQA] = new TH1D(Form("fhQAEventsNumContrPV_%s",sQAindex[iQA].Data()), "QA Events: Number of contributors to AOD PV", 20,0,20);
        fQAEvents->Add(fhQAEventsNumContrPV[iQA]);
        fhQAEventsNumSPDContrPV[iQA] = new TH1D(Form("fhQAEventsNumSPDContrPV_%s",sQAindex[iQA].Data()), "QA Events: SPD contributors to PV", 20,0,20);
        fQAEvents->Add(fhQAEventsNumSPDContrPV[iQA]);
        fhQAEventsDistPVSPD[iQA] = new TH1D(Form("fhQAEventsDistPVSPD_%s",sQAindex[iQA].Data()), "QA Events: PV SPD vertex", 50,0,5);
        fQAEvents->Add(fhQAEventsDistPVSPD[iQA]);
        fhQAEventsSPDresol[iQA] = new TH1D(Form("fhQAEventsSPDresol_%s",sQAindex[iQA].Data()), "QA Events: SPD resolution", 150,0,15);
        fQAEvents->Add(fhQAEventsSPDresol[iQA]);

        // Charged tracks QA
        if(fProcessCharged)
        {
          fhQAChargedMult[iQA] = new TH1D(Form("fhQAChargedMult_%s",sQAindex[iQA].Data()),"QA Charged: Number of Charged in selected events; #it{N}^{Charged}", 1500,0,1500);
          fQACharged->Add(fhQAChargedMult[iQA]);
          fhQAChargedCharge[iQA] = new TH1D(Form("fhQAChargedCharge_%s",sQAindex[iQA].Data()),"QA Charged: Track charge; charge;", 3,-1.5,1.5);
          fQACharged->Add(fhQAChargedCharge[iQA]);
          fhQAChargedPt[iQA] = new TH1D(Form("fhQAChargedPt_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 300,0.,30.);
          fQACharged->Add(fhQAChargedPt[iQA]);
          fhQAChargedEta[iQA] = new TH1D(Form("fhQAChargedEta_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#eta}; #it{#eta}", 151,-1.5,1.5);
          fQACharged->Add(fhQAChargedEta[iQA]);
          fhQAChargedPhi[iQA] = new TH1D(Form("fhQAChargedPhi_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#varphi}; #it{#varphi}", 100,0.,TMath::TwoPi());
          fQACharged->Add(fhQAChargedPhi[iQA]);
          fhQAChargedFilterBit[iQA] = new TH1D(Form("fhQAChargedFilterBit_%s",sQAindex[iQA].Data()), "QA Charged: Filter bit",iNFilterMapBinBins,0,iNFilterMapBinBins);
          for(Int_t j = 0x0; j < iNFilterMapBinBins; j++) fhQAChargedFilterBit[iQA]->GetXaxis()->SetBinLabel(j+1, Form("%g",TMath::Power(2,j)));
          fQACharged->Add(fhQAChargedFilterBit[iQA]);
          fhQAChargedNumTPCcls[iQA] = new TH1D(Form("fhQAChargedNumTPCcls_%s",sQAindex[iQA].Data()),"QA Charged: Track number of TPC clusters; #it{N}^{TPC clusters}", 160,0,160);
          fQACharged->Add(fhQAChargedNumTPCcls[iQA]);
          fhQAChargedDCAxy[iQA] = new TH1D(Form("fhQAChargedDCAxy_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-xy; DCA_{#it{xy}} (cm)", 100,0.,10);
          fQACharged->Add(fhQAChargedDCAxy[iQA]);
          fhQAChargedDCAz[iQA] = new TH1D(Form("fhQAChargedDCAz_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-z; DCA_{#it{z}} (cm)", 200,-10.,10.);
          fQACharged->Add(fhQAChargedDCAz[iQA]);
        } // endif {fProcessCharged}

        // PID tracks QA
        if(fProcessPID || fProcessPhi)
        {
          fhQAPIDTPCstatus[iQA] = new TH1D(Form("fhQAPIDTPCstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAPID->Add(fhQAPIDTPCstatus[iQA]);
          fhQAPIDTPCdEdx[iQA] = new TH2D(Form("fhQAPIDTPCdEdx_%s",sQAindex[iQA].Data()),"QA PID: TPC PID information; #it{p} (GeV/#it{c}); TPC dEdx (au)", 100,0,10, 131,-10,1000);
          fQAPID->Add(fhQAPIDTPCdEdx[iQA]);
          fhQAPIDTOFstatus[iQA] = new TH1D(Form("fhQAPIDTOFstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAPID->Add(fhQAPIDTOFstatus[iQA]);
          fhQAPIDTOFbeta[iQA] = new TH2D(Form("fhQAPIDTOFbeta_%s",sQAindex[iQA].Data()),"QA PID: TOF #beta information; #it{p} (GeV/#it{c}); TOF #beta", 100,0,10, 101,-0.1,1.5);
          fQAPID->Add(fhQAPIDTOFbeta[iQA]);
          fh3QAPIDnSigmaBayesElectron[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesElectron_%s",sQAindex[iQA].Data()),"QA PID: e; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 22,-11,10, 22,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesElectron[iQA]);
          fh3QAPIDnSigmaBayesMuon[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesMuon_%s",sQAindex[iQA].Data()),"QA PID: #mu; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 22,-11,10, 22,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesMuon[iQA]);
          fh3QAPIDnSigmaBayesPion[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesPion_%s",sQAindex[iQA].Data()),"QA PID: #pi; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 22,-11,10, 22,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesPion[iQA]);
          fh3QAPIDnSigmaBayesKaon[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesKaon_%s",sQAindex[iQA].Data()),"QA PID: K; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 22,-11,10, 22,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesKaon[iQA]);
          fh3QAPIDnSigmaBayesProton[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesProton_%s",sQAindex[iQA].Data()),"QA PID: p; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 22,-11,10, 22,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesProton[iQA]);

          for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
          {
            fhQAPIDTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
            fhQAPIDTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          }
        } // endif {fProcessPID}

        // V0s QA
        if(fProcessV0s)
        {
          fhQAV0sMultK0s[iQA] = new TH1D(Form("fhQAV0sMultK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of K^{0}_{S} candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultK0s[iQA]);
          fhQAV0sMultLambda[iQA] = new TH1D(Form("fhQAV0sMultLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #Lambda candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultLambda[iQA]);
          fhQAV0sMultALambda[iQA] = new TH1D(Form("fhQAV0sMultALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #bar{#Lambda} candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultALambda[iQA]);
          fhQAV0sRecoMethod[iQA] = new TH1D(Form("fhQAV0sRecoMethod_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Reconstruction method", 2,-0.5,1.5);
          fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(1, "offline");
          fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(2, "online (on-the-fly)");
          fQAV0s->Add(fhQAV0sRecoMethod[iQA]);
          fhQAV0sDCAtoPV[iQA] = new TH1D(Form("fhQAV0sDCAtoPV_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter DCA to PV; daughter DCA^{PV} (cm)", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDCAtoPV[iQA]);
          fhQAV0sDCADaughters[iQA] = new TH1D(Form("fhQAV0sDCADaughters_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: DCA among daughters; DCA^{daughters} (cm)", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDCADaughters[iQA]);
          fhQAV0sDecayRadius[iQA] = new TH1D(Form("fhQAV0sDecayRadius_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Decay radius; #it{r_{xy}}^{decay} (cm)", 300,0.,300.);
          fQAV0s->Add(fhQAV0sDecayRadius[iQA]);
          fhQAV0sCPAK0s[iQA] = new TH1D(Form("fhQAV0sCPAK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: CPA; CPA^{K0s}", 100,0.9,1.);
          fQAV0s->Add(fhQAV0sCPAK0s[iQA]);
          fhQAV0sCPALambda[iQA] = new TH1D(Form("fhQAV0sCPALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: CPA; CPA^{#Lambda}", 100, 0.9,1.);
          fQAV0s->Add(fhQAV0sCPALambda[iQA]);
          fhQAV0sNumTauK0s[iQA] = new TH1D(Form("fhQAV0sNumTauK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{K0s} (cm)", 100, 0.,20.);
          fQAV0s->Add(fhQAV0sNumTauK0s[iQA]);
          fhQAV0sNumTauLambda[iQA] = new TH1D(Form("fhQAV0sNumTauLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{#Lambda} (cm)", 100, 0.,60);
          fQAV0s->Add(fhQAV0sNumTauLambda[iQA]);
          fhQAV0sArmenterosK0s[iQA] = new TH2D(Form("fhQAV0sArmenterosK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosK0s[iQA]);
          fhQAV0sArmenterosLambda[iQA] = new TH2D(Form("fhQAV0sArmenterosLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosLambda[iQA]);
          fhQAV0sArmenterosALambda[iQA] = new TH2D(Form("fhQAV0sArmenterosALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosALambda[iQA]);
          fhQAV0sInvMassK0s[iQA] = new TH1D(Form("fhQAV0sInvMassK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: InvMass; #it{m}_{inv} (GeV/#it{c}^{2});", 200,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
          fQAV0s->Add(fhQAV0sInvMassK0s[iQA]);
          fhQAV0sInvMassLambda[iQA] = new TH1D(Form("fhQAV0sInvMassLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: InvMass; #it{m}_{inv} (GeV/#it{c}^{2});", 80,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
          fQAV0s->Add(fhQAV0sInvMassLambda[iQA]);
          fhQAV0sMotherPt[iQA] = new TH1D(Form("fhQAV0sMotherPt_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{p}_{T}; #it{p}_{T}^{V0} (GeV/#it{c})", 200,0.,20.);
          fQAV0s->Add(fhQAV0sMotherPt[iQA]);
          fhQAV0sMotherPhi[iQA] = new TH1D(Form("fhQAV0sMotherPhi_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{#varphi}; #it{#varphi}^{V0} (GeV/#it{c})", 100,0.,TMath::TwoPi());
          fQAV0s->Add(fhQAV0sMotherPhi[iQA]);
          fhQAV0sMotherEta[iQA] = new TH1D(Form("fhQAV0sMotherEta_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{#eta}; #it{#eta}^{V0}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherEta[iQA]);
          fhQAV0sMotherCharge[iQA] = new TH1D(Form("fhQAV0sMotherCharge_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother charge; V^{0} charge", 3,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherCharge[iQA]);
          fhQAV0sMotherRapK0s[iQA] = new TH1D(Form("fhQAV0sMotherRapK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (K^{0}_{S} hypo); #it{y}^{V0,K0s}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherRapK0s[iQA]);
          fhQAV0sMotherRapLambda[iQA] = new TH1D(Form("fhQAV0sMotherRapLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (Lambda/#bar{#Lambda} hypo); #it{y}^{V0,#Lambda}", 301,-3.,3.);
          fQAV0s->Add(fhQAV0sMotherRapLambda[iQA]);
          fhQAV0sDaughterTPCRefit[iQA] = new TH1D(Form("fhQAV0sDaughterTPCRefit_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter TPC refit", 2,-0.5,1.5);
          fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODTrack::kTPCrefit");
          fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(2, "AliAODTrack::kTPCrefit");
          fQAV0s->Add(fhQAV0sDaughterTPCRefit[iQA]);
          fhQAV0sDaughterKinks[iQA] = new TH1D(Form("fhQAV0sDaughterKinks_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter Kinks", 2,-0.5,1.5);
          fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODVertex::kKink");
          fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(2, "AliAODVertex:kKink");
          fQAV0s->Add(fhQAV0sDaughterKinks[iQA]);
          fhQAV0sDaughterPt[iQA] = new TH1D(Form("fhQAV0sDaughterPt_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{p}_{T}; #it{p}_{T}^{daughter} (GeV/#it{c})", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDaughterPt[iQA]);
          fhQAV0sDaughterPhi[iQA] = new TH1D(Form("fhQAV0sDaughterPhi_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{#varphi}; #it{#varphi}^{daughter} (GeV/#it{c})", 100,0.,TMath::TwoPi());
          fQAV0s->Add(fhQAV0sDaughterPhi[iQA]);
          fhQAV0sDaughterEta[iQA] = new TH1D(Form("fhQAV0sDaughterEta_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{#eta}; #it{#eta}^{daugter}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sDaughterEta[iQA]);
          fhQAV0sDaughterCharge[iQA] = new TH1D(Form("fhQAV0sDaughterCharge_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter charge; daughter charge", 3,-1.5,1.5);
          fQAV0s->Add(fhQAV0sDaughterCharge[iQA]);
          fhQAV0sDaughterTPCstatus[iQA] = new TH1D(Form("fhQAV0sDaughterTPCstatus_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAV0s->Add(fhQAV0sDaughterTPCstatus[iQA]);
          fhQAV0sDaughterTOFstatus[iQA] = new TH1D(Form("fhQAV0sDaughterTOFstatus_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAV0s->Add(fhQAV0sDaughterTOFstatus[iQA]);
          fhQAV0sDaughterTPCdEdxK0s[iQA] = new TH2D(Form("fhQAV0sDaughterTPCdEdxK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: TPC dEdx daughters; #it{p}^{daughter} (GeV/#it{c}); TPC dEdx (au);", 100,0.,20, 101,-10,1000);
          fQAV0s->Add(fhQAV0sDaughterTPCdEdxK0s[iQA]);
          fhQAV0sDaughterNumSigmaPionK0s[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: Daughter PID (#pi); #it{p}_{T}^{daughter} (GeV/#it{c}); #pi PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionK0s[iQA]);
          fhQAV0sDaughterTPCdEdxLambda[iQA] = new TH2D(Form("fhQAV0sDaughterTPCdEdxLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: TPC dEdx daughters; #it{p}^{daughter} (GeV/#it{c}); TPC dEdx (au);", 100,0.,20, 101,-10,1000);
          fQAV0s->Add(fhQAV0sDaughterTPCdEdxLambda[iQA]);
          fhQAV0sDaughterNumSigmaPionLambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Daughter PID (#pi); #it{p}_{T}^{pion} (GeV/#it{c}); pion PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionLambda[iQA]);
          fhQAV0sDaughterNumSigmaProtonLambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaProtonLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Daughter PID (p); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaProtonLambda[iQA]);
          fhQAV0sDaughterNumSigmaPionALambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Daughter PID (#pi); #it{p}_{T}^{pion} (GeV/#it{c}); pion PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionALambda[iQA]);
          fhQAV0sDaughterNumSigmaProtonALambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaProtonALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Daughter PID (p); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaProtonALambda[iQA]);

          for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
          {
            fhQAV0sDaughterTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
            fhQAV0sDaughterTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          }
        } // endif {fProcessV0s}
      }
    }

  // posting data (mandatory)
  PostData(1, fFlowRefs);
  PostData(2, fFlowCharged);
  PostData(3, fFlowPID);
  PostData(4, fFlowPhi);
  PostData(5, fFlowK0s);
  PostData(6, fFlowLambda);
  PostData(7, fQAEvents);
  PostData(8, fQACharged);
  PostData(9, fQAPID);
  PostData(10, fQAPhi);
  PostData(11, fQAV0s);
  PostData(12, fFlowWeights);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListParameters()
{
  // lists all task parameters
  // *************************************************************
  printf("\n======= List of parameters ========================================\n");
  printf("   -------- Analysis task ---------------------------------------\n");
  printf("      fRunMode: (RunMode) %d\n",    fRunMode);
  printf("      fAnalType: (AnalType) %d\n",    fAnalType);
  printf("      fSampling: (Bool_t) %s\n",    fSampling ? "kTRUE" : "kFALSE");
  printf("      fFillQA: (Bool_t) %s\n",    fFillQA ? "kTRUE" : "kFALSE");
  printf("      fProcessCharged: (Bool_t) %s\n",    fProcessCharged ? "kTRUE" : "kFALSE");
  printf("      fProcessPID: (Bool_t) %s\n",    fProcessPID ? "kTRUE" : "kFALSE");
  printf("      fProcessPhi: (Bool_t) %s\n",    fProcessPhi ? "kTRUE" : "kFALSE");
  printf("      fProcessV0s: (Bool_t) %s\n",    fProcessV0s ? "kTRUE" : "kFALSE");
  printf("   -------- Flow related ----------------------------------------\n");
  printf("      fCutFlowDoFourCorrelations: (Bool_t) %s\n",    fCutFlowDoFourCorrelations ? "kTRUE" : "kFALSE");
  printf("      fCutFlowRFPsPtMin: (Float_t) %g (GeV/c)\n",    fCutFlowRFPsPtMin);
  printf("      fCutFlowRFPsPtMax: (Float_t) %g (GeV/c)\n",    fCutFlowRFPsPtMax);
  printf("      fFlowPOIsPtMin: (Float_t) %g (GeV/c)\n",    fFlowPOIsPtMin);
  printf("      fFlowPOIsPtMax: (Float_t) %g (GeV/c)\n",    fFlowPOIsPtMax);
  printf("      fFlowCentNumBins: (Int_t) %d (GeV/c)\n",    fFlowCentNumBins);
  printf("      fFlowCentMin: (Int_t) %d (GeV/c)\n",    fFlowCentMin);
  printf("      fFlowCentMax: (Int_t) %d (GeV/c)\n",    fFlowCentMax);
  printf("      fFlowUseWeights: (Bool_t) %s\n",    fFlowUseWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowWeightsPath: (TString) '%s' \n",    fFlowWeightsPath.Data());
  printf("   -------- Events ----------------------------------------------\n");
  printf("      fColSystem: (ColSystem) %d\n",    fColSystem);
  printf("      fPeriod: (DataPeriod) %d\n",    fPeriod);
  printf("      fTrigger: (Short_t) %d\n",    fTrigger);
  printf("      fMultEstimator: (TString) '%s'\n",    fMultEstimator.Data());
  printf("      fPVtxCutZ: (Double_t) %g (cm)\n",    fPVtxCutZ);
  printf("   -------- Charge tracks ---------------------------------------\n");
  printf("      fCutChargedTrackFilterBit: (UInt) %d\n",    fCutChargedTrackFilterBit);
  printf("      fCutChargedNumTPCclsMin: (UShort_t) %d\n",    fCutChargedNumTPCclsMin);
  printf("      fCutChargedEtaMax: (Float_t) %g\n",    fCutChargedEtaMax);
  printf("      fCutChargedPtMin: (Float_t) %g (GeV/c)\n",    fCutChargedPtMin);
  printf("      fCutChargedPtMax: (Float_t) %g (GeV/c)\n",    fCutChargedPtMax);
  printf("      fCutChargedDCAzMax: (Float_t) %g (cm)\n",    fCutChargedDCAzMax);
  printf("      fCutChargedDCAxyMax: (Float_t) %g (cm)\n",    fCutChargedDCAxyMax);
  printf("   -------- PID (pi,K,p) tracks ---------------------------------\n");
  printf("      fCutPIDUseAntiProtonOnly: (Bool_t) %s\n",  fCutPIDUseAntiProtonOnly ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaCombinedNoTOFrejection: (Bool_t) %s\n",  fCutPIDnSigmaCombinedNoTOFrejection ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaPionMax: (Double_t) %g\n",    fCutPIDnSigmaPionMax);
  printf("      fCutPIDnSigmaKaonMax: (Double_t) %g\n",    fCutPIDnSigmaKaonMax);
  printf("      fCutPIDnSigmaProtonMax: (Double_t) %g\n",    fCutPIDnSigmaProtonMax);
  printf("      fCutUseBayesPID: (Bool_t) %s\n",    fCutUseBayesPID ? "kTRUE" : "kFALSE");
  printf("      fCutPIDBayesPionMin: (Double_t) %g\n",    fCutPIDBayesPionMin);
  printf("      fCutPIDBayesKaonMin: (Double_t) %g\n",    fCutPIDBayesKaonMin);
  printf("      fCutPIDBayesProtonMin: (Double_t) %g\n",    fCutPIDBayesProtonMin);
  printf("      SetPIDBayesRejectElectron: (Double_t) %g\n",    fCutPIDBayesRejectElectron);
  printf("      SetPIDBayesRejectMuon: (Double_t) %g\n",    fCutPIDBayesRejectMuon);
  printf("   -------- Phi candidates --------------------------------------\n");
  printf("      fCutPhiMotherEtaMax: (Double_t) %g\n",    fCutPhiMotherEtaMax);
  printf("      fCutPhiInvMassMin: (Double_t) %g\n",    fCutPhiInvMassMin);
  printf("      fCutPhiInvMassMax: (Double_t) %g\n",    fCutPhiInvMassMax);
  printf("   -------- V0s candidates --------------------------------------\n");
  printf("      fCutV0sOnFly: (Bool_t) %s\n",    fCutV0sOnFly ? "kTRUE" : "kFALSE");
  printf("      fCutV0srefitTPC: (Bool_t) %s\n",     fCutV0srefitTPC ? "kTRUE" : "kFALSE");
  printf("      fCutV0srejectKinks: (Bool_t) %s\n",     fCutV0srejectKinks ? "kTRUE" : "kFALSE");
  printf("      fCutV0sCrossMassRejection: (Bool_t) %s\n",     fCutV0sCrossMassRejection ? "kTRUE" : "kFALSE");
  printf("      fCutV0sCrossMassCutK0s: (Double_t) %g (GeV/c2)\n",     fCutV0sCrossMassCutK0s);
  printf("      fCutV0sCrossMassCutLambda: (Double_t) %g (GeV/c2)\n",     fCutV0sCrossMassCutLambda);
  printf("      fCutV0sDCAtoPVMin: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVMin);
  printf("      fCutV0sDCAtoPVMax: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVMax);
  printf("      fCutV0sDCAtoPVzMax: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVzMax);
  printf("      fCutV0sDaughterFilterBit: (UInt) %d\n",    fCutV0sDaughterFilterBit);
  printf("      fCutV0sDCADaughtersMin: (Double_t) %g (cm)\n",    fCutV0sDCADaughtersMin);
  printf("      fCutV0sDCADaughtersMax: (Double_t) %g (cm)\n",    fCutV0sDCADaughtersMax);
  printf("      fCutV0sDecayRadiusMin: (Double_t) %g (cm)\n",    fCutV0sDecayRadiusMin);
  printf("      fCutV0sDecayRadiusMax: (Double_t) %g (cm)\n",    fCutV0sDecayRadiusMax);
  printf("      fCutV0sDaughterPtMin: (Double_t) %g (GeV/c)\n",    fCutV0sDaughterPtMin);
  printf("      fCutV0sDaughterPtMax: (Double_t) %g (GeV/c)\n",    fCutV0sDaughterPtMax);
  printf("      fCutV0sDaughterEtaMax: (Double_t) %g ()\n",    fCutV0sDaughterEtaMax);
  printf("      fCutV0sMotherEtaMax: (Double_t) %g ()\n",    fCutV0sMotherEtaMax);
  printf("      fCutV0sMotherRapMax: (Double_t) %g ()\n",    fCutV0sMotherRapMax);
  printf("      fCutV0sInvMassK0sMin: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassK0sMin);
  printf("      fCutV0sInvMassK0sMax: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassK0sMax);
  printf("      fCutV0sInvMassLambdaMin: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassLambdaMin);
  printf("      fCutV0sInvMassLambdaMax: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassLambdaMax);
  printf("      fCutV0sCPAK0sMin: (Double_t) %g ()\n",    fCutV0sCPAK0sMin);
  printf("      fCutV0sCPALambdaMin: (Double_t) %g ()\n",    fCutV0sCPALambdaMin);
  printf("      fCutV0sNumTauK0sMax: (Double_t) %g (c*tau)\n",    fCutV0sNumTauK0sMax);
  printf("      fCutV0sNumTauLambdaMax: (Double_t) %g (c*tau)\n",    fCutV0sNumTauLambdaMax);
  printf("      fCutV0sArmenterosAlphaK0sMin: (Double_t) %g (alpha)\n",    fCutV0sArmenterosAlphaK0sMin);
  printf("      fCutV0sArmenterosAlphaLambdaMax: (Double_t) %g (alpha)\n",    fCutV0sArmenterosAlphaLambdaMax);
  printf("      fCutV0sProtonNumSigmaMax: (Double_t) %g (n*sigma)\n",    fCutV0sProtonNumSigmaMax);
  printf("      fCutV0sProtonPIDPtMin: (Double_t) %g (GeV/c)\n",    fCutV0sProtonPIDPtMin);
  printf("      fCutV0sProtonPIDPtMax: (Double_t) %g (GeV/c)\n",    fCutV0sProtonPIDPtMax);
  printf("=====================================================================\n\n");

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::InitializeTask()
{
  // called once on beginning of task (within UserCreateOutputObjects method)
  // check if task parameters are specified and valid
  // returns kTRUE if succesfull
  // *************************************************************

  printf("====== InitializeTask AliAnalysisTaskUniFlow =========================\n");

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

  if(fColSystem != kPP && fColSystem != kPPb && fColSystem != kPbPb)
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

  fPIDCombined = new AliPIDCombined();
  if(!fPIDCombined)
  {
    ::Error("InitializeTask","AliPIDCombined object not found! Terminating!");
    return kFALSE;
  }
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetSelectedSpecies(5); // all particle species
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask


  if(fSampling && fNumSamples == 0)
  {
    ::Error("InitializeTask","Sampling used, but number of samples is 0! Terminating!");
    return kFALSE;
  }

  if(!fSampling)
  //fNumSamples = 1;

  // checking cut setting
  ::Info("InitializeTask","Checking task parameters setting conflicts (ranges, etc)");
  if(fCutFlowRFPsPtMin > 0. && fCutFlowRFPsPtMax > 0. && fCutFlowRFPsPtMin > fCutFlowRFPsPtMax)
  {
    ::Error("InitializeTask","Cut: RFPs Pt range wrong!");
    return kFALSE;
  }
  if(fCutV0sInvMassK0sMin > fCutV0sInvMassK0sMax || fCutV0sInvMassK0sMin < 0. || fCutV0sInvMassK0sMax < 0.)
  {
    ::Error("InitializeTask","Cut: InvMass (K0s) range wrong!");
    return kFALSE;
  }
  if(fCutV0sInvMassLambdaMin > fCutV0sInvMassLambdaMax || fCutV0sInvMassLambdaMin < 0. || fCutV0sInvMassLambdaMax < 0.)
  {
    ::Error("InitializeTask","Cut: InvMass (Lambda) range wrong!");
    return kFALSE;
  }

  // upper-case for multiplicity estimator
  fMultEstimator.ToUpper();

  // checking for weights source file
  if(fFlowUseWeights && !fFlowWeightsPath.EqualTo(""))
  {
    fFlowWeightsFile = TFile::Open(Form("alien:///%s",fFlowWeightsPath.Data()));
    if(!fFlowWeightsFile)
    {
      ::Error("InitializeTask","Flow weights file not found");
      return kFALSE;
    }
  }

  ::Info("InitializeTask","Initialization succesfull!");
  printf("======================================================================\n\n");
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

  // fIndexCentrality = GetCentralityIndex();

  // processing of selected event
  if(!ProcessEvent()) return;

  // posting data (mandatory)
  PostData(1, fFlowRefs);
  PostData(2, fFlowCharged);
  PostData(3, fFlowPID);
  PostData(4, fFlowPhi);
  PostData(5, fFlowK0s);
  PostData(6, fFlowLambda);
  PostData(7, fQAEvents);
  PostData(8, fQACharged);
  PostData(9, fQAPID);
  PostData(10, fQAPhi);
  PostData(11, fQAV0s);
  PostData(12, fFlowWeights);

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
  if(fFillQA) FillEventsQA(0);

  // event selection for small systems pp, pPb in Run2 (2016)
  if( (fColSystem == kPP || fColSystem == kPPb)
      && (fPeriod == k16k || fPeriod == k16l || fPeriod == k16q || fPeriod == k16r || fPeriod == k16s || fPeriod == k16t)
    ) eventSelected = IsEventSelected_2016();

  // Fill event QA AFTER cuts
  if(fFillQA && eventSelected) FillEventsQA(1);

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
void AliAnalysisTaskUniFlow::Filtering()
{
  // main (envelope) method for filtering all particles of interests (POIs) in selected events
  // All POIs passing selection criteria are saved to relevant TClonesArray for further processing
  // return kTRUE if succesfull (no errors in process)
  // *************************************************************

  if(!fProcessCharged && !fProcessPID && !fProcessV0s && !fProcessPhi) // if neither is ON, filtering is skipped
    return;

  fVectorCharged->clear();
  FilterCharged();

  // estimate centrality & assign indexes (centrality/percentile, sampling, ...)
  fIndexCentrality = GetCentralityIndex();
  if(fIndexCentrality < 0) return; // not succesfull estimation
  fhEventCentrality->Fill(fIndexCentrality);
  fh2EventCentralityNumSelCharged->Fill(fVectorCharged->size(),fIndexCentrality);

  fIndexSampling = GetSamplingIndex();
  fhEventSampling->Fill(fIndexCentrality,fIndexSampling);

  if(fProcessPID || fProcessPhi)
  {
    fVectorPion->clear();
    fVectorKaon->clear();
    fVectorProton->clear();
    FilterPID();
  }

  if(fProcessPhi)
  {
    fVectorPhi->clear();
    FilterPhi();
    // printf("Num Phi: %d | Num BG: %d\n",fVectorPhi->size(),fVectorPhiBG->size());
  }

  if(fProcessV0s)
  {
    fVectorK0s->clear();
    fVectorLambda->clear();
    FilterV0s();
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterCharged()
{
  // Filtering input charged tracks
  // If track passes all requirements as defined in IsChargedSelected(),
  // the relevant properties (pT, eta, phi) are stored in FlowPart struct
  // and pushed to relevant vector container.
  // return kFALSE if any complications occurs
  // *************************************************************

  const Short_t iNumTracks = fEventAOD->GetNumberOfTracks();
  if(iNumTracks < 1) return;

  AliAODTrack* track = 0x0;
  Int_t iNumRefs = 0;
  Double_t weight = 0;

  for(Short_t iTrack(0); iTrack < iNumTracks; iTrack++)
  {
    track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
    if(!track) continue;

    if(fFillQA) FillQACharged(0,track); // QA before selection

    if(IsChargedSelected(track))
    {
      fVectorCharged->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kCharged) );
      if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsCharged->Fill(track->Phi(),track->Eta(),track->Pt());
      if(fFlowUseWeights)
      {
        weight = fh2WeightCharged->GetBinContent( fh2WeightCharged->FindBin(track->Eta(),track->Phi()) );
        fh3AfterWeightsCharged->Fill(track->Phi(),track->Eta(),track->Pt(),weight);
      }
      if(fFillQA) FillQACharged(1,track); // QA after selection
      //printf("pt %g | phi %g | eta %g\n",track->Pt(),track->Phi(),track->Eta());

      // Filling refs QA plots
      if(fCutFlowRFPsPtMin > 0. && track->Pt() >= fCutFlowRFPsPtMin && fCutFlowRFPsPtMax > 0. && track->Pt() <= fCutFlowRFPsPtMax)
      {
        iNumRefs++;
        if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsRefs->Fill(track->Phi(),track->Eta(),track->Pt());
        if(fFlowUseWeights)
        {
          weight = fh2WeightRefs->GetBinContent( fh2WeightRefs->FindBin(track->Eta(),track->Phi()) );
          fh3AfterWeightsRefs->Fill(track->Phi(),track->Eta(),track->Pt(),weight);
        }
        FillQARefs(1,track);
      }
    }
  }

  // fill QA charged multiplicity
  fhRefsMult->Fill(iNumRefs);
  if(fFillQA)
  {
    fhQAChargedMult[0]->Fill(fEventAOD->GetNumberOfTracks());
    fhQAChargedMult[1]->Fill(fVectorCharged->size());
  }

  return;
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
void AliAnalysisTaskUniFlow::FillQARefs(const Short_t iQAindex, const AliAODTrack* track)
{
  // Filling various QA plots related to RFPs subset of charged track selection
  // *************************************************************

  if(!track) return;
  if(iQAindex == 0) return; // NOTE implemented only for selected RFPs

  fhRefsPt->Fill(track->Pt());
  fhRefsEta->Fill(track->Eta());
  fhRefsPhi->Fill(track->Phi());

  return;
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

  // kinematics
  fhQAChargedPt[iQAindex]->Fill(track->Pt());
  fhQAChargedPhi[iQAindex]->Fill(track->Phi());
  fhQAChargedEta[iQAindex]->Fill(track->Eta());

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterV0s()
{
  // Filtering input V0s candidates (K0s, (Anti)Lambda)
  // If track passes all requirements as defined in IsV0sSelected() (and species dependent one)
  // the relevant properties (pT, eta, phi,mass,species) are stored in FlowPart struct
  // and pushed to relevant vector container.
  // return kFALSE if any complications occurs
  // *************************************************************

  Short_t iNumK0sSelected = 0;  // counter for selected K0s candidates
  Short_t iNumLambdaSelected = 0; // counter for selected Lambda candidates
  Short_t iNumALambdaSelected = 0; // counter for selected Anti-Lambda candidates

  const Short_t iNumV0s = fEventAOD->GetNumberOfV0s();
  if(iNumV0s < 1) return;

  Bool_t bIsK0s = kFALSE;
  Short_t iIsLambda = 0;
  Double_t weight = 0;

  AliAODv0* v0 = 0x0;
  for(Short_t iV0(0); iV0 < iNumV0s; iV0++)
  {
    // the minimalistic dynamic allocation of the TClonesArray*
    v0 = static_cast<AliAODv0*>(fEventAOD->GetV0(iV0));
    if(!v0) continue;

    if(fFillQA) FillQAV0s(0,v0); // QA BEFORE selection

    // Alex selections
    if(fDoAlexK0sSelection)
    {
      if(IsV0SelectedK0sAlex(v0))
      {
        // selected K0s according to Alex
        iNumK0sSelected++;
        fhV0sCounter->Fill("K^{0}_{S}",1);
        fVectorK0s->emplace_back( FlowPart(v0->Pt(),v0->Phi(),v0->Eta(), 0, kK0s, v0->MassK0Short()) );
        if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsK0s->Fill(v0->Phi(),v0->Eta(),v0->Pt());
        if(fFlowUseWeights)
        {
          weight = fh2WeightK0s->GetBinContent( fh2WeightK0s->FindBin(v0->Eta(),v0->Phi()) );
          fh3AfterWeightsK0s->Fill(v0->Phi(),v0->Eta(),v0->Pt(),weight);
        }
      }
    }
    else
    {
      if(IsV0Selected(v0))
      {
        bIsK0s = IsV0aK0s(v0);
        iIsLambda = IsV0aLambda(v0);

        if(fFillQA && (bIsK0s || iIsLambda != 0))
        FillQAV0s(1,v0,bIsK0s,iIsLambda); // QA AFTER selection

        if(bIsK0s)
        {
          iNumK0sSelected++;
          fhV0sCounter->Fill("K^{0}_{S}",1);
          fVectorK0s->emplace_back( FlowPart(v0->Pt(),v0->Phi(),v0->Eta(), 0, kK0s, v0->MassK0Short()) );
          if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsK0s->Fill(v0->Phi(),v0->Eta(),v0->Pt());
          if(fFlowUseWeights)
          {
            weight = fh2WeightK0s->GetBinContent( fh2WeightK0s->FindBin(v0->Eta(),v0->Phi()) );
            fh3AfterWeightsK0s->Fill(v0->Phi(),v0->Eta(),v0->Pt(),weight);
          }
        }

        if(iIsLambda == 1) // lambda
        {
          iNumLambdaSelected++;
          fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
          fVectorLambda->emplace_back( FlowPart(v0->Pt(),v0->Phi(),v0->Eta(), 0, kLambda, v0->MassLambda()) );
          if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsLambda->Fill(v0->Phi(),v0->Eta(),v0->Pt());
          if(fFlowUseWeights)
          {
            weight = fh2WeightLambda->GetBinContent( fh2WeightLambda->FindBin(v0->Eta(),v0->Phi()) );
            fh3AfterWeightsLambda->Fill(v0->Phi(),v0->Eta(),v0->Pt(),weight);
          }
        }

        if(iIsLambda == -1) // anti-lambda
        {
          iNumALambdaSelected++;
          fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
          fVectorLambda->emplace_back( FlowPart(v0->Pt(),v0->Phi(),v0->Eta(), 0, kLambda, v0->MassAntiLambda()) );
          if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsLambda->Fill(v0->Phi(),v0->Eta(),v0->Pt());
          if(fFlowUseWeights)
          {
            weight = fh2WeightLambda->GetBinContent( fh2WeightLambda->FindBin(v0->Eta(),v0->Phi()) );
            fh3AfterWeightsLambda->Fill(v0->Phi(),v0->Eta(),v0->Pt(),weight);
          }
        }

        if(bIsK0s && iIsLambda != 0)
        fhV0sCounter->Fill("K^{0}_{S} && #Lambda/#bar{#Lambda}",1);
      }

    }

  }

  // fill QA charged multiplicity
  if(fFillQA)
  {
    fhQAV0sMultK0s[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultLambda[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultALambda[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultK0s[1]->Fill(iNumK0sSelected);
    fhQAV0sMultLambda[1]->Fill(iNumLambdaSelected);
    fhQAV0sMultALambda[1]->Fill(iNumALambdaSelected);
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0aK0s(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for K0s candidates
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;
  fhV0sCounterK0s->Fill("Input",1);

  // cosine of pointing angle (CPA)
  if(fCutV0sCPAK0sMin > 0.)
  {
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    Double_t dCPA = v0->CosPointingAngle(primVtx);
    if( dCPA < fCutV0sCPAK0sMin ) return kFALSE;
  }
  fhV0sCounterK0s->Fill("CPA",1);

  // proper life-time
  if( fCutV0sNumTauK0sMax > 0. )
  {
    AliAODVertex* primVtx2 = fEventAOD->GetPrimaryVertex();
    Double_t dPrimVtxCoor[3] = {0}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0}; // decay vector coor {xyz}
    primVtx2->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 2; i++)
      dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];

    Double_t dPropLife = ( (fPDGMassK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    if( dPropLife > (fCutV0sNumTauK0sMax * 2.68) ) return kFALSE;
  }
  fhV0sCounterK0s->Fill("c#tau",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaK0sMin > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm < (fCutV0sArmenterosAlphaK0sMin * TMath::Abs(dAlpha))) return kFALSE;
  }
  fhV0sCounterK0s->Fill("Armenteros-Podolanski",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0. && ( TMath::Abs(v0->RapK0Short()) > fCutV0sMotherRapMax ) ) return kFALSE;
  fhV0sCounterK0s->Fill("#it{y}",1);

  // inv. mass window
  Double_t dMass = v0->MassK0Short();
  if( dMass < fCutV0sInvMassK0sMin || dMass > fCutV0sInvMassK0sMax ) return kFALSE;
  fhV0sCounterK0s->Fill("InvMass",1);

  // competing V0 rejection based on InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassLambda = v0->MassLambda();
    Double_t dMassALambda = v0->MassAntiLambda();

    // K0s candidate is within 10 MeV of (Anti)Lambda InvMass physSelTask
    if(TMath::Abs(dMassLambda - fPDGMassLambda) < fCutV0sCrossMassCutK0s)
    {
      // in Lambda peak
      fhV0sCompetingInvMassK0s->Fill(dMass,dMassLambda);
      return kFALSE;
    }

    if(TMath::Abs(dMassALambda - fPDGMassLambda) < fCutV0sCrossMassCutK0s)
    {
      // in Anti-Lambda peak
      fhV0sCompetingInvMassK0s->Fill(dMass,dMassALambda);
      return kFALSE;
    }
  }
  fhV0sCounterK0s->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterK0s->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::IsV0aLambda(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for Lambda candidates
  // return 0 if candidate does not fullfill any Lambda or Anti-Lambda requirements;
  // return 1 if a candidate fulfill all Lambda requirements;
  // return -1 if a candidate fullfill all Anti-Lambda requirements;
  // return 2 if a candidate fulfill all both Lambda & Anti-Lambda requirements
  // *************************************************************
  if(!v0) return 0;
  fhV0sCounterLambda->Fill("Input",1);

  // cosine of pointing angle (CPA)
  if( fCutV0sCPALambdaMin > 0. )
  {
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    Double_t dCPA = v0->CosPointingAngle(primVtx);
    if( dCPA < fCutV0sCPALambdaMin ) return 0;
  }
  fhV0sCounterLambda->Fill("CPA",1);

  // proper life-time
  if( fCutV0sNumTauLambdaMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0}; // decay vector coor {xyz}
    AliAODVertex* primVtx2 = fEventAOD->GetPrimaryVertex();
    primVtx2->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 2; i++)
      dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];

    Double_t dPropLife = ( (fPDGMassLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    if( dPropLife > (fCutV0sNumTauLambdaMax * 7.89) ) return 0;
  }
  fhV0sCounterLambda->Fill("c#tau",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaLambdaMax > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm > (fCutV0sArmenterosAlphaLambdaMax * TMath::Abs(dAlpha))) return 0;
  }
  fhV0sCounterLambda->Fill("Armenteros-Podolanski",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0. && ( TMath::Abs(v0->RapLambda()) > fCutV0sMotherRapMax) ) return 0;
  fhV0sCounterLambda->Fill("#it{y}",1);

  // particle species dependent
  Bool_t bIsLambda = kFALSE;
  Bool_t bIsALambda = kFALSE;

  // inv. mass window
  Double_t dMassLambda = v0->MassLambda();
  Double_t dMassALambda = v0->MassAntiLambda();
  if( dMassLambda > fCutV0sInvMassLambdaMin && dMassLambda < fCutV0sInvMassLambdaMax)
    bIsLambda = kTRUE;

  if( dMassALambda > fCutV0sInvMassLambdaMin && dMassALambda < fCutV0sInvMassLambdaMax)
    bIsALambda = kTRUE;

  if(!bIsLambda && !bIsALambda) // if neither
    return 0;
  fhV0sCounterLambda->Fill("InvMass",1);

  // Armenteros-Podolanski for candidates fullfilling both Lambda and Anti-Lambda selection
  if(bIsLambda && bIsALambda)
  {
    if(v0->AlphaV0() < 0.) bIsLambda = kFALSE;
    if(v0->AlphaV0() > 0.) bIsALambda = kFALSE;
  }

  // competing V0 rejection based on InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassK0s = v0->MassK0Short();
    if( TMath::Abs(dMassK0s - fPDGMassK0s) < fCutV0sCrossMassCutLambda)
    {
      // Lambda candidate is within 5 MeV of K0s InvMass physSelTask

      Double_t dMass = 0;
      if(bIsLambda) dMass = dMassLambda;
      if(bIsALambda) dMass = dMassALambda;

      fhV0sCompetingInvMassLambda->Fill(dMass,dMassK0s);
      return kFALSE;
    }
  }
  fhV0sCounterLambda->Fill("Competing InvMass",1);


  // proton PID of Lambda Candidates
  if( (fCutV0sProtonNumSigmaMax > 0.) && (fCutV0sProtonPIDPtMax > 0. || fCutV0sProtonPIDPtMin > 0.) && fPIDResponse)
  {
    const AliAODTrack* trackDaughterPos = (AliAODTrack*) v0->GetDaughter(0); // positive charge
    const AliAODTrack* trackDaughterNeg = (AliAODTrack*) v0->GetDaughter(1); // negative charge

    // checking detector status
    AliPIDResponse::EDetPidStatus pidStatusTPCpos = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughterPos);
    AliPIDResponse::EDetPidStatus pidStatusTPCneg = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughterNeg);
    if(pidStatusTPCpos != AliPIDResponse::kDetPidOk || pidStatusTPCneg != AliPIDResponse::kDetPidOk) return 0;

    Double_t dSigmaProtPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kProton));
    Double_t dSigmaProtNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kProton));

    // positive daughter is not a proton (within TPC sigmas)
    if(bIsLambda && (trackDaughterPos->Pt() > fCutV0sProtonPIDPtMin) && (trackDaughterPos->Pt() < fCutV0sProtonPIDPtMax) && (dSigmaProtPos > fCutV0sProtonNumSigmaMax) )
      return 0;

    // negative daughter is not a proton (within TPC sigmas)
    if(bIsALambda && (trackDaughterNeg->Pt() > fCutV0sProtonPIDPtMin) && (trackDaughterNeg->Pt() < fCutV0sProtonPIDPtMax) && (dSigmaProtNeg > fCutV0sProtonNumSigmaMax) )
      return 0;
  }
  fhV0sCounterLambda->Fill("Proton PID",1);

  // passing all criteria
  fhV0sCounterLambda->Fill("Selected",1);

  if(bIsLambda && bIsALambda) // both Lambda & Anti-Lambda
  {
    fhV0sCounterLambda->Fill("#Lambda && #bar{#Lambda}",1);
    return 2;
  }

  if(bIsLambda) // only Lambda
  {
    fhV0sCounterLambda->Fill("only #Lambda",1);
    return 1;
  }

  if(bIsALambda) // only Anti-Lambda
  {
    fhV0sCounterLambda->Fill("only #bar{#Lambda}",1);
    return -1;
  }

  // default return 0: should return earlier if candidates is either Lambda or ALambda
  return 0;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0SelectedK0sAlex(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // Copy of Alexandru's K0s selection in AliAnalysisTaskPiKpK0Lamba class.
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************

  if (!v0) return kFALSE;

  // V0 cuts
  Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
  Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
  Double_t lV0CosineOfPointingAngle = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
  Double_t lDcaV0Daughters = v0->DcaV0Daughters();
  Double_t lV0Radius = v0->RadiusV0();

  Double_t  lV0Position[3];
  v0->GetXYZ(lV0Position);
//        Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
//                TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
//                TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2));

  if ( (lDcaPosToPrimVertex < fDCADghtPV) || (lDcaNegToPrimVertex < fDCADghtPV) || (lDcaV0Daughters > fMaxDCADght) || (lV0CosineOfPointingAngle < fCosPA) || (lV0Radius < fMinRad) || (lV0Radius > fMaxRad) || (v0->Pt() < fMinPt) || (v0->Pt() >= fMaxPt) || (TMath::Abs(v0->Eta()) >= fEtaCut) )
       return kFALSE;



  if (fArmPodCut){

      Double_t ptArm = v0->PtArmV0();
      Double_t angAlpha = v0->AlphaV0();

      if (ptArm <= 0.2*TMath::Abs(angAlpha))
           return kFALSE;

  }



  // Tracks quality cuts
  const AliAODTrack* negTrack = (AliAODTrack *)v0->GetDaughter(1);
  if (!negTrack){
      delete negTrack;
       return kFALSE;
  }

  const AliAODTrack* posTrack =(AliAODTrack *)v0->GetDaughter(0);
  if (!posTrack){
      delete posTrack;
       return kFALSE;
  }


  if (!posTrack->IsOn(AliAODTrack::kTPCrefit))
       return kFALSE;

  if (!negTrack->IsOn(AliAODTrack::kTPCrefit))
       return kFALSE;



  if (posTrack->Charge() == negTrack->Charge()){
      //cout<< "like sign, continue"<< endl;
       return kFALSE;
  }


  if ((posTrack->GetTPCNcls() < fNoClus) || (negTrack->GetTPCNcls() < fNoClus))
       return kFALSE;


  Float_t nCrossedRowsTPCPos = posTrack->GetTPCNCrossedRows();
  if (nCrossedRowsTPCPos < fNoClus)
       return kFALSE;
  Int_t findablePos = posTrack->GetTPCNclsF();
  if (findablePos <= 0)
       return kFALSE;
  if (nCrossedRowsTPCPos/findablePos < fNcrFind)
       return kFALSE;


  Float_t nCrossedRowsTPCNeg = negTrack->GetTPCNCrossedRows();
  if (nCrossedRowsTPCNeg < fNoClus)
       return kFALSE;
  Int_t findableNeg = negTrack->GetTPCNclsF();
  if (findableNeg <= 0)
       return kFALSE;
  if (nCrossedRowsTPCNeg/findableNeg < fNcrFind)
       return kFALSE;


  if (TMath::Abs(posTrack->Eta()) > fEtaCut || TMath::Abs(negTrack->Eta()) > fEtaCut)
       return kFALSE;


  if (fMinPtDght){
      if (posTrack->Pt() < fMinPt || negTrack->Pt() < fMinPt)
           return kFALSE;
  }



  if (negTrack->GetTPCsignalN() < fNoClusPid || posTrack->GetTPCsignalN() < fNoClusPid)
       return kFALSE;

  Double_t nSigmaPipos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion));
  Double_t nSigmaPineg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion));

  Double_t nSigmaPpos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kProton));
  Double_t nSigmaPneg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kProton));

  Double_t lInvMassK0s = v0->MassK0Short();
  Double_t lInvMassL = v0->MassLambda();

  Short_t flagV0 = -1;

  //K0
  if (nSigmaPipos < fNsigCut && nSigmaPineg < fNsigCut && lInvMassK0s > 0.4 && lInvMassK0s < 0.6)
  {

      Double_t rapK0 = GetRapidity(0.497648, v0->Pt(), v0->Eta());

      if (TMath::Abs(rapK0) < 0.5){ flagV0 = 0; }
    }

    //Lambda
    if ( ((nSigmaPpos < fNsigCut && nSigmaPineg < fNsigCut) || (nSigmaPneg < fNsigCut || nSigmaPipos < fNsigCut)) && lInvMassL > 1.07 && lInvMassL < 1.17){

        Double_t rapL = GetRapidity(1.115683, v0->Pt(), v0->Eta());

        if (TMath::Abs(rapL) < 0.5){ flagV0 = 1; }
    }

    if(flagV0 == 0) return kTRUE;
    else return kFALSE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUniFlow::GetRapidity(Double_t mass, Double_t Pt, Double_t Eta)
{
    Double_t rapid = TMath::Log( (TMath::Sqrt(mass*mass + Pt*Pt*TMath::CosH(Eta)*TMath::CosH(Eta)) + Pt*TMath::SinH(Eta)) / TMath::Sqrt(mass*mass + Pt*Pt) );
    return rapid;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0Selected(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // common for both K0s and (Anti)-Lambdas
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;
  fhV0sCounter->Fill("Input",1);

  const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);


  // daughter track check
  if(!daughterPos || !daughterNeg) return kFALSE;
  fhV0sCounter->Fill("Daughters OK",1);

  // filter bit
  // printf("FB: %d\n",daughterPos->TestFilterBit(fCutV0sDaughterFilterBit));
  if( fCutV0sDaughterFilterBit > 0 && (!daughterPos->TestFilterBit(fCutV0sDaughterFilterBit) || !daughterNeg->TestFilterBit(fCutV0sDaughterFilterBit) ) ) return kFALSE;

  // daughters & mother charge checks
  if( (TMath::Abs(daughterPos->Charge()) != 1) || (TMath::Abs(daughterNeg->Charge()) != 1) ) return kFALSE;
  if(daughterPos->Charge() == daughterNeg->Charge()) return kFALSE;
  if(v0->Charge() != 0) return kFALSE;
  fhV0sCounter->Fill("Charge",1);

  // reconstruction method: online (on-the-fly) OR offline
  if(v0->GetOnFlyStatus() != fCutV0sOnFly) return kFALSE;
  fhV0sCounter->Fill("Reconstruction method",1);

  // TPC refit
  if(fCutV0srefitTPC && ( !daughterPos->IsOn(AliAODTrack::kTPCrefit) || !daughterNeg->IsOn(AliAODTrack::kTPCrefit) ) ) return kFALSE;
  fhV0sCounter->Fill("TPC refit",1);

  // Kinks
  const AliAODVertex* prodVtxDaughterPos = (AliAODVertex*) daughterPos->GetProdVertex(); // production vertex of the positive daughter track
  const AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*) daughterNeg->GetProdVertex(); // production vertex of the negative daughter track
  if(fCutV0srejectKinks && ( (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) || (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) ) ) return kFALSE;
  fhV0sCounter->Fill("Kinks",1);


  // Daughters DCA to PV
  const Float_t dDCAPosToPV = TMath::Abs(v0->DcaPosToPrimVertex());
  const Float_t dDCANegToPV = TMath::Abs(v0->DcaNegToPrimVertex());

  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  Double_t dVertexXYZ[3] = {0.};
  Double_t dTrackXYZpos[3] = {0};
  Double_t dTrackXYZneg[3] = {0};
  Double_t dDCAXYZpos[3] = {0.};
  Double_t dDCAXYZneg[3] = {0.};
  if( fCutV0sDCAtoPVzMax > 0. )
  {
    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) return kFALSE; // event does not have a PV

    daughterPos->GetXYZ(dTrackXYZpos);
    daughterNeg->GetXYZ(dTrackXYZneg);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++)
    {
      dDCAXYZpos[i] = dTrackXYZpos[i] - dVertexXYZ[i];
      dDCAXYZneg[i] = dTrackXYZneg[i] - dVertexXYZ[i];
    }
  }

  if(fCutV0sDCAtoPVMin > 0. && ( dDCAPosToPV < fCutV0sDCAtoPVMin || dDCANegToPV < fCutV0sDCAtoPVMin ) ) return kFALSE;
  if(fCutV0sDCAtoPVMax > 0. && ( dDCAPosToPV > fCutV0sDCAtoPVMax || dDCANegToPV > fCutV0sDCAtoPVMax ) ) return kFALSE;
  if(fCutV0sDCAtoPVzMax > 0. && ( TMath::Abs(dDCAXYZpos[2]) > fCutV0sDCAtoPVzMax || TMath::Abs(dDCAXYZneg[2]) > fCutV0sDCAtoPVzMax ) ) return kFALSE;
  fhV0sCounter->Fill("DCA to PV",1);


  // Daughter DCA among themselves
  if(fCutV0sDCADaughtersMin > 0. && TMath::Abs(v0->DcaV0Daughters()) < fCutV0sDCADaughtersMin) return kFALSE;
  if(fCutV0sDCADaughtersMax > 0. && TMath::Abs(v0->DcaV0Daughters()) > fCutV0sDCADaughtersMax) return kFALSE;
  fhV0sCounter->Fill("Daughters DCA",1);

  // radius of decay vertex in transverse plane
  Double_t dSecVtxCoor[3] = {0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  // can be substituted by AliAODv0::RadiusV0() - equivalent
  if( fCutV0sDecayRadiusMin > 0. && (dDecayRadius < fCutV0sDecayRadiusMin) ) return kFALSE;
  if( fCutV0sDecayRadiusMax > 0. && (dDecayRadius > fCutV0sDecayRadiusMax) ) return kFALSE;
  fhV0sCounter->Fill("Decay radius",1);

  // acceptance checks
  if(fCutV0sDaughterEtaMax > 0. && ( (TMath::Abs(daughterNeg->Eta()) > fCutV0sDaughterEtaMax) || (TMath::Abs(daughterPos->Eta()) > fCutV0sDaughterEtaMax) ) ) return kFALSE;
  if(fCutV0sDaughterPtMin > 0. && (daughterPos->Pt() < fCutV0sDaughterPtMin  || daughterNeg->Pt() < fCutV0sDaughterPtMin) ) return kFALSE;
  if(fCutV0sDaughterPtMax > 0. && (daughterPos->Pt() > fCutV0sDaughterPtMax  || daughterNeg->Pt() > fCutV0sDaughterPtMax) ) return kFALSE;

  if(fCutV0sMotherEtaMax > 0. && TMath::Abs(v0->Eta()) > fCutV0sDaughterEtaMax ) return kFALSE;
  fhV0sCounter->Fill("Acceptance",1);

  // passing all common criteria
  fhV0sCounter->Fill("Common passed",1);
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAV0s(const Short_t iQAindex, const AliAODv0* v0, const Bool_t bIsK0s, const Short_t bIsLambda)
{
  // Filling various QA plots related to V0 candidate selection
  // *************************************************************
  // checking mother & daughters
  if(!v0) return;
  AliAODTrack* trackDaughter[2] = {(AliAODTrack*) v0->GetDaughter(0), (AliAODTrack*) v0->GetDaughter(1)};
  if(!trackDaughter[0] || !trackDaughter[1]) return;

  // setting internal flags for Lambdas and Anti-Lambdas
  Bool_t bCandLambda;
  Bool_t bCandAntiLambda;

  switch (bIsLambda)
  {
    case 1:
    {
      bCandLambda = kTRUE;
      bCandAntiLambda = kFALSE;
      break;
    }
    case -1:
    {
      bCandLambda = kFALSE;
      bCandAntiLambda = kTRUE;
      break;
    }
    case 2:
    {
      bCandLambda = kTRUE;
      bCandAntiLambda = kTRUE;
      break;
    }
    default:
    {
      bCandLambda = kFALSE;
      bCandAntiLambda = kFALSE;
    }
  }

  // reconstruction method
  fhQAV0sRecoMethod[iQAindex]->Fill(v0->GetOnFlyStatus());

  // DCA between daughters and PV
  fhQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaPosToPrimVertex());
  fhQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaNegToPrimVertex());

  // Daughter DCA among themselves
  fhQAV0sDCADaughters[iQAindex]->Fill(v0->DcaV0Daughters());

  // charge
  fhQAV0sMotherCharge[iQAindex]->Fill(v0->Charge());

  // radius of decay vertex in transverse plane
  Double_t dSecVtxCoor[3] = {0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  fhQAV0sDecayRadius[iQAindex]->Fill(dDecayRadius);

  // mother kinematics
  fhQAV0sMotherPt[iQAindex]->Fill(v0->Pt());
  fhQAV0sMotherPhi[iQAindex]->Fill(v0->Phi());
  fhQAV0sMotherEta[iQAindex]->Fill(v0->Eta());

  // proper lifetime preparation (to be filled in particle dependent if scope)
  Double_t dPrimVtxCoor[3] = {0};
  Double_t dDecayCoor[3] = {0};
  AliAODVertex* primVtx2 = fEventAOD->GetPrimaryVertex();
  primVtx2->GetXYZ(dPrimVtxCoor);
  for(Int_t i(0); i < 2; i++)
    dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];

  // particle dependent
  if(bIsK0s)
  {
    // K0s
    fhQAV0sMotherRapK0s[iQAindex]->Fill(v0->RapK0Short());
    fhQAV0sInvMassK0s[iQAindex]->Fill(v0->MassK0Short());

    // CPA
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    fhQAV0sCPAK0s[iQAindex]->Fill(v0->CosPointingAngle(primVtx));

    // Armenteros-Podolanski
    fhQAV0sArmenterosK0s[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    // proper lifetime
    Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    Double_t dPropLifeK0s = ( (dMassPDGK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    fhQAV0sNumTauK0s[iQAindex]->Fill(dPropLifeK0s);
  }
  if(bCandLambda || bCandAntiLambda)
  {
    // (Anti)Lambda
    fhQAV0sMotherRapLambda[iQAindex]->Fill(v0->RapLambda());
    fhQAV0sInvMassLambda[iQAindex]->Fill(v0->MassLambda());
    fhQAV0sInvMassLambda[iQAindex]->Fill(v0->MassAntiLambda());

    // CPA
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    fhQAV0sCPALambda[iQAindex]->Fill(v0->CosPointingAngle(primVtx));

    // Armenteros-Podolanski
    if(bCandLambda)
      fhQAV0sArmenterosLambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    if(bCandAntiLambda)
      fhQAV0sArmenterosALambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    // proper lifetime
    Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
    Double_t dPropLifeLambda = ( (dMassPDGLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    fhQAV0sNumTauLambda[iQAindex]->Fill(dPropLifeLambda);
  }

  AliPIDResponse::EDetPidStatus pidStatusTPC;
  AliPIDResponse::EDetPidStatus pidStatusTOF;

  // daughters properties
  AliAODVertex* prodVtxDaughter = 0x0;
  for(Short_t i(0); i < 2; i++)
  {
    // TPC refit
    fhQAV0sDaughterTPCRefit[iQAindex]->Fill(trackDaughter[i]->IsOn(AliAODTrack::kTPCrefit));

    // kinks
    prodVtxDaughter = (AliAODVertex*) trackDaughter[i]->GetProdVertex();
    fhQAV0sDaughterKinks[iQAindex]->Fill(prodVtxDaughter->GetType() == AliAODVertex::kKink);

    // detector status
    pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[i]);
    pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trackDaughter[i]);
    fhQAV0sDaughterTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );
    fhQAV0sDaughterTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );

    // daughter kinematics
    fhQAV0sDaughterPt[iQAindex]->Fill(trackDaughter[i]->Pt());
    fhQAV0sDaughterPhi[iQAindex]->Fill(trackDaughter[i]->Phi());
    fhQAV0sDaughterEta[iQAindex]->Fill(trackDaughter[i]->Eta());

    // daughter charge
    fhQAV0sDaughterCharge[iQAindex]->Fill(trackDaughter[i]->Charge());
  }

  AliPIDResponse::EDetPidStatus pidStatusTPCpos;
  AliPIDResponse::EDetPidStatus pidStatusTPCneg;

  // PID checks
  if(fPIDResponse)
  {
    // checking the detector status
    pidStatusTPCpos = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[0]);
    pidStatusTPCneg = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[1]);

    if(pidStatusTPCpos == AliPIDResponse::kDetPidOk && pidStatusTPCneg == AliPIDResponse::kDetPidOk)
    {
      if(bIsK0s)
      {
        // daughter PID
        fhQAV0sDaughterTPCdEdxK0s[iQAindex]->Fill(trackDaughter[0]->P(), trackDaughter[0]->GetTPCsignal());
        fhQAV0sDaughterTPCdEdxK0s[iQAindex]->Fill(trackDaughter[1]->P(), trackDaughter[1]->GetTPCsignal());

        // Pion PID for daughters
        fhQAV0sDaughterNumSigmaPionK0s[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kPion));
        fhQAV0sDaughterNumSigmaPionK0s[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kPion));
      }

      if(bCandLambda || bCandAntiLambda)
      {
        // daughter PID
        fhQAV0sDaughterTPCdEdxLambda[iQAindex]->Fill(trackDaughter[0]->P(), trackDaughter[0]->GetTPCsignal());
        fhQAV0sDaughterTPCdEdxLambda[iQAindex]->Fill(trackDaughter[1]->P(), trackDaughter[1]->GetTPCsignal());

        if(bCandLambda)
        {
          fhQAV0sDaughterNumSigmaProtonLambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kProton));
          fhQAV0sDaughterNumSigmaPionLambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kPion));
        }

        if(bCandAntiLambda)
        {
          fhQAV0sDaughterNumSigmaProtonALambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kProton));
          fhQAV0sDaughterNumSigmaPionALambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kPion));
        }
      }
    }
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterPhi()
{
  // Reconstruction and filtering of Phi meson candidates out of selected Kaon sample
  // If track passes all requirements, the relevant properties (pT, eta, phi) are stored
  // in FlowPart struct  and pushed to relevant vector container.
  // *************************************************************

  // printf("== enter FilterPhi ==\n");

  const Int_t iNumKaons = fVectorKaon->size();
  // printf("Num Kaons: %d\n", iNumKaons);
  // check if there are at least 2 selected kaons in event (minimum for phi reconstruction)
  if(iNumKaons < 2) return;

  // start Phi reconstruction
  FlowPart* kaon1 = 0x0;
  FlowPart* kaon2 = 0x0;
  FlowPart mother = FlowPart();

  Int_t iNumBG = 0;
  for(Short_t iKaon1(0); iKaon1 < iNumKaons; iKaon1++)
  {
    kaon1 = &(fVectorKaon->at(iKaon1));
    if(!kaon1) continue;

    for(Short_t iKaon2(iKaon1+1); iKaon2 < iNumKaons; iKaon2++)
    {
      kaon2 = &(fVectorKaon->at(iKaon2));
      if(!kaon2) continue;

      fhPhiCounter->Fill("Input",1);
      mother = FlowPart::MakeMother(kaon1,kaon2,kPhi);

      // filling QA BEFORE selection
      if(fFillQA) FillQAPhi(0,&mother);

      if(mother.mass < fCutPhiInvMassMin || mother.mass > fCutPhiInvMassMax) continue;
      fhPhiCounter->Fill("InvMass",1);

      if(mother.pt < fFlowPOIsPtMin || mother.pt > fFlowPOIsPtMax) continue;
      fhPhiCounter->Fill("Pt",1);

      if(fCutPhiMotherEtaMax > 0 && TMath::Abs(mother.eta) > fCutPhiMotherEtaMax)  continue;
      fhPhiCounter->Fill("Eta",1);

      // mother (phi) candidate passing all criteria
      fhPhiCounter->Fill("Selected",1);

      // filling QA AFTER selection
      if(fFillQA) FillQAPhi(1,&mother);

      // filling weights
      if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsPhi->Fill(mother.phi, mother.eta, mother.pt);
      if(fFlowUseWeights)
      {
        Double_t weight = fh2WeightPhi->GetBinContent( fh2WeightPhi->FindBin(mother.eta,mother.phi) );
        fh3AfterWeightsPhi->Fill(mother.phi,mother.eta,mother.pt,weight);
      }

      if(mother.charge == 0)
      {
        // opposite-sign combination (signal+background)
        fhPhiCounter->Fill("Unlike-sign",1);
        fVectorPhi->emplace_back(mother);
      }

      if(TMath::Abs(mother.charge) == 2)
      {
        // like-sign combination (background)
        fhPhiCounter->Fill("BG",1);
        iNumBG++;

        // filing background entries for Phi candidates
        for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
        {
          if(TMath::Abs(mother.eta) > fEtaGap[iGap]/2 ) fh3PhiEntriesBG[iGap]->Fill(fIndexCentrality,mother.pt,mother.mass);
        }
      }
    } // endfor {iKaon2} : second kaon
  } // endfor {iKaon1} : first Kaon

  // filling multiplicity distribution
  fhPhiMult->Fill(fVectorPhi->size());
  fhPhiBGMult->Fill(iNumBG);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAPhi(const Short_t iQAindex, const FlowPart* part)
{
  if(!part) return;
  if(part->species != kPhi) return;

  if(iQAindex == 0) return; // TODO not implemented (do not know what)

  if(part->charge == 0)
  {
    // phi candidate (unlike-sign pair)
    fhPhiInvMass->Fill(part->mass);
    fhPhiCharge->Fill(part->charge);
    fhPhiPt->Fill(part->pt);
    fhPhiEta->Fill(part->eta);
    fhPhiPhi->Fill(part->phi);
  }

  if(TMath::Abs(part->charge) == 2)
  {
    // phi candidate (unlike-sign pair)
    fhPhiBGInvMass->Fill(part->mass);
    fhPhiBGCharge->Fill(part->charge);
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterPID()
{
  // Filtering input PID tracks (pi,K,p)
  // If track passes all requirements as defined in IsPIDSelected() (and species dependent),
  // the relevant properties (pT, eta, phi) are stored in FlowPart struct
  // and pushed to relevant vector container.
  // return kFALSE if any complications occurs
  // *************************************************************

  const Short_t iNumTracks = fEventAOD->GetNumberOfTracks();
  if(iNumTracks < 1) return;

  PartSpecies species = kUnknown;
  AliAODTrack* track = 0x0;
  Double_t weight = 0;

  for(Short_t iTrack(0); iTrack < iNumTracks; iTrack++)
  {
    track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
    if(!track) continue;

    // PID tracks are subset of selected charged tracks (same quality requirements)
    if(!IsChargedSelected(track)) continue;

    if(fFillQA) FillPIDQA(0,track,kUnknown);   // filling QA for tracks before selection (but after charged criteria applied)

    // PID track selection (return most favourable species)
    PartSpecies species = IsPIDSelected(track);
    // check if only protons should be used
    if(fCutPIDUseAntiProtonOnly && species == kProton && track->Charge() == 1) species = kUnknown;

    // selection of PID tracks
    switch (species)
    {
      case kPion:
        fVectorPion->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kPion, fPDGMassPion, track->Px(), track->Py(), track->Pz()) );
        if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsPion->Fill(track->Phi(), track->Eta(), track->Pt());
        if(fFlowUseWeights)
        {
          weight = fh2WeightPion->GetBinContent( fh2WeightPion->FindBin(track->Eta(),track->Phi()) );
          fh3AfterWeightsPion->Fill(track->Phi(),track->Eta(),track->Pt(),weight);
        }
        break;
      case kKaon:
        fVectorKaon->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kKaon, fPDGMassKaon, track->Px(), track->Py(), track->Pz()) );
        if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsKaon->Fill(track->Phi(), track->Eta(), track->Pt());
        if(fFlowUseWeights)
        {
          weight = fh2WeightKaon->GetBinContent( fh2WeightKaon->FindBin(track->Eta(),track->Phi()) );
          fh3AfterWeightsKaon->Fill(track->Phi(),track->Eta(),track->Pt(),weight);
        }
        break;
      case kProton:
        fVectorProton->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kProton, fPDGMassProton, track->Px(), track->Py(), track->Pz()) );
        if(fRunMode == kFillWeights || fFlowFillWeights) fh3WeightsProton->Fill(track->Phi(), track->Eta(), track->Pt());
        if(fFlowUseWeights)
        {
          weight = fh2WeightProton->GetBinContent( fh2WeightProton->FindBin(track->Eta(),track->Phi()) );
          fh3AfterWeightsProton->Fill(track->Phi(),track->Eta(),track->Pt(),weight);
        }
        break;
      default:
        continue;
    }

    if(fFillQA) FillPIDQA(1,track,species); // filling QA for tracks AFTER selection
  }

  fhPIDPionMult->Fill(fVectorPion->size());
  fhPIDKaonMult->Fill(fVectorKaon->size());
  fhPIDProtonMult->Fill(fVectorProton->size());

  return;
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::PartSpecies AliAnalysisTaskUniFlow::IsPIDSelected(const AliAODTrack* track)
{
  // Selection of PID tracks (pi,K,p) - track identification
  // Based on fCutUseBayesPID flag, either Bayes PID or nSigma cutting is used
  // returns AliAnalysisTaskUniFlow::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
  // *************************************************************

  // checking detector statuses
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);

  Bool_t bIsTPCok = (pidStatusTPC == AliPIDResponse::kDetPidOk);
  Bool_t bIsTOFok = ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME)); // checking TOF

  if(!bIsTPCok) return kUnknown;
  // TODO: TOF check???

  const Double_t dPt = track->Pt();

  if(fCutUseBayesPID)
  {
    // use Bayesian PID
    Double_t dProbPID[5] = {0}; // array for Bayes PID probabilities:  0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dProbPID); // filling probabilities to dPropPID array
    Double_t dMaxProb = TMath::MaxElement(5,dProbPID);

    // printf("PID Prob: e %g | mu %g | pi %g | K %g | p %g ||| MAX %g \n",dProbPID[0],dProbPID[1],dProbPID[2],dProbPID[3],dProbPID[4],dMaxProb);

    // electron and mion rejection
    if(dProbPID[0] >= fCutPIDBayesRejectElectron || dProbPID[1] >= fCutPIDBayesRejectMuon) return kUnknown;

    // checking the PID probability
    // TODO: think about: if Pion has maximum probibility < fCutBayesPIDPion, track is rejected -> is it good?
    if(dMaxProb == dProbPID[2] && dProbPID[2] >= fCutPIDBayesPionMin) return kPion;
    if(dMaxProb == dProbPID[3] && dProbPID[3] >= fCutPIDBayesKaonMin) return kKaon;
    if(dMaxProb == dProbPID[4] && dProbPID[4] >= fCutPIDBayesProtonMin) return kProton;
  }
  else
  {
    // use nSigma cuts (based on combination of TPC / TOF nSigma cuts)

    Double_t dNumSigmaTPC[5] = {-99,-99,-99,-99,-99}; // TPC nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    Double_t dNumSigmaTOF[5] = {-99,-99,-99,-99,-99}; // TOF nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton

    // filling nSigma arrays
    if(bIsTPCok) // should be anyway
    {
      dNumSigmaTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));
      dNumSigmaTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon));
      dNumSigmaTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
      dNumSigmaTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
      dNumSigmaTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    }

    if(bIsTOFok) // should be anyway
    {
      dNumSigmaTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron));
      dNumSigmaTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon));
      dNumSigmaTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
      dNumSigmaTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
      dNumSigmaTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton));
    }

    // TPC nSigma cuts
    if(dPt <= 0.4)
    {
      Double_t dMinSigmasTPC = TMath::MinElement(5,dNumSigmaTPC);

      // electron rejection
      if(dMinSigmasTPC == dNumSigmaTPC[0] && dNumSigmaTPC[0] <= fCutPIDnSigmaTPCRejectElectron) return kUnknown;
      if(dMinSigmasTPC == dNumSigmaTPC[2] && dNumSigmaTPC[2] <= fCutPIDnSigmaPionMax) return kPion;
      if(dMinSigmasTPC == dNumSigmaTPC[3] && dNumSigmaTPC[3] <= fCutPIDnSigmaKaonMax) return kKaon;
      if(dMinSigmasTPC == dNumSigmaTPC[4] && dNumSigmaTPC[4] <= fCutPIDnSigmaProtonMax) return kProton;
    }

    // combined TPC + TOF nSigma cuts
    if(dPt > 0.4) // && < 4 GeV TODO once TPC dEdx parametrisation is available
    {
      Double_t dNumSigmaCombined[5] = {-99,-99,-99,-99,-99};

      // discard candidates if no TOF is available if cut is on
      if(fCutPIDnSigmaCombinedNoTOFrejection && !bIsTOFok) return kUnknown;

      // calculating combined nSigmas
      for(Short_t i(0); i < 5; i++)
      {
        if(bIsTOFok) { dNumSigmaCombined[i] = TMath::Sqrt(dNumSigmaTPC[i]*dNumSigmaTPC[i] + dNumSigmaTOF[i]*dNumSigmaTOF[i]); }
        else { dNumSigmaCombined[i] = dNumSigmaTPC[i]; }
      }

      Double_t dMinSigmasCombined = TMath::MinElement(5,dNumSigmaCombined);

      // electron rejection
      if(dMinSigmasCombined == dNumSigmaCombined[0] && dNumSigmaCombined[0] <= fCutPIDnSigmaPionMax) return kUnknown;
      if(dMinSigmasCombined == dNumSigmaCombined[2] && dNumSigmaCombined[2] <= fCutPIDnSigmaPionMax) return kPion;
      if(dMinSigmasCombined == dNumSigmaCombined[3] && dNumSigmaCombined[3] <= fCutPIDnSigmaKaonMax) return kKaon;
      if(dMinSigmasCombined == dNumSigmaCombined[4] && dNumSigmaCombined[4] <= fCutPIDnSigmaProtonMax) return kProton;
    }

    // TPC dEdx parametrisation (dEdx - <dEdx>)
    // TODO: TPC dEdx parametrisation cuts
    // if(dPt > 3.)
    // {
    //
    // }
  }

  return kUnknown;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillPIDQA(const Short_t iQAindex, const AliAODTrack* track, const PartSpecies species)
{
  // Filling various QA plots related to PID (pi,K,p) track selection
  // *************************************************************
  if(!track) return;

  if(!fPIDResponse || !fPIDCombined)
  {
    ::Error("FillPIDQA","AliPIDResponse or AliPIDCombined object not found!");
    return;
  }

  // TPC & TOF statuses & measures
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);

  fhQAPIDTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
  fhQAPIDTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );

  Bool_t bIsTPCok = (pidStatusTPC == AliPIDResponse::kDetPidOk);
  Bool_t bIsTOFok = ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));

  Double_t dNumSigmaTPC[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Double_t dNumSigmaTOF[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Double_t dBayesProb[5] = {-0.1}; // Bayesian probability | array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton

  Double_t dTPCdEdx = -5; // TPC dEdx for selected particle
  Double_t dTOFbeta = -0.05; //TOF beta for selected particle

  // filling Bayesian PID probabilities to dBayesProb array
  UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dBayesProb);

  Double_t dP = track->P();
  Double_t dPt = track->Pt();

  // detector status dependent
  if(bIsTPCok)
  {
    dNumSigmaTPC[0] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    dNumSigmaTPC[1] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon);
    dNumSigmaTPC[2] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    dNumSigmaTPC[3] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    dNumSigmaTPC[4] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    dTPCdEdx = track->GetTPCsignal();
    fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), dTPCdEdx);
  }
  else // TPC status not OK
  {
    dNumSigmaTPC[0] = -11.;
    dNumSigmaTPC[1] = -11.;
    dNumSigmaTPC[2] = -11.;
    dNumSigmaTPC[3] = -11.;
    dNumSigmaTPC[4] = -11.;

    fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), -5.);
  }

  if(bIsTOFok)
  {
    dNumSigmaTOF[0] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    dNumSigmaTOF[1] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon);
    dNumSigmaTOF[2] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    dNumSigmaTOF[3] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    dNumSigmaTOF[4] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    Double_t dTOF[5];
    track->GetIntegratedTimes(dTOF);
    dTOFbeta = dTOF[0] / track->GetTOFsignal();
    fhQAPIDTOFbeta[iQAindex]->Fill(dP,dTOFbeta);
  }
  else // TOF status not OK
  {
    dNumSigmaTOF[0] = -11.;
    dNumSigmaTOF[1] = -11.;
    dNumSigmaTOF[2] = -11.;
    dNumSigmaTOF[3] = -11.;
    dNumSigmaTOF[4] = -11.;

    fhQAPIDTOFbeta[iQAindex]->Fill(track->P(),-0.05);
  }

  fh3QAPIDnSigmaBayesElectron[iQAindex]->Fill(dNumSigmaTPC[0],dNumSigmaTOF[0],dBayesProb[0]);
  fh3QAPIDnSigmaBayesMuon[iQAindex]->Fill(dNumSigmaTPC[1],dNumSigmaTOF[1],dBayesProb[1]);
  fh3QAPIDnSigmaBayesPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],dBayesProb[2]);
  fh3QAPIDnSigmaBayesKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],dBayesProb[3]);
  fh3QAPIDnSigmaBayesProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],dBayesProb[4]);

  // species dependent QA
  switch (species)
  {
    case kPion:
      fhPIDPionPt->Fill(track->Pt());
      fhPIDPionPhi->Fill(track->Phi());
      fhPIDPionEta->Fill(track->Eta());
      fhPIDPionCharge->Fill(track->Charge());
      fh2PIDPionTPCdEdx->Fill(dPt,dTPCdEdx);
      fh2PIDPionTOFbeta->Fill(dPt,dTOFbeta);
      fh2PIDPionTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDPionTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDPionTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDPionTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDPionTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDPionTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDPionBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDPionBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDPionBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    case kKaon:
      fhPIDKaonPt->Fill(track->Pt());
      fhPIDKaonPhi->Fill(track->Phi());
      fhPIDKaonEta->Fill(track->Eta());
      fhPIDKaonCharge->Fill(track->Charge());
      fh2PIDKaonTPCdEdx->Fill(dP,dTPCdEdx);
      fh2PIDKaonTOFbeta->Fill(dP,dTOFbeta);
      fh2PIDKaonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDKaonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDKaonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDKaonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDKaonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDKaonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDKaonBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDKaonBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDKaonBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    case kProton:
      fhPIDProtonPt->Fill(track->Pt());
      fhPIDProtonPhi->Fill(track->Phi());
      fhPIDProtonEta->Fill(track->Eta());
      fhPIDProtonCharge->Fill(track->Charge());
      fh2PIDProtonTPCdEdx->Fill(dP,dTPCdEdx);
      fh2PIDProtonTOFbeta->Fill(dP,dTOFbeta);
      fh2PIDProtonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDProtonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDProtonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDProtonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDProtonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDProtonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDProtonBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDProtonBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDProtonBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    default:
      break;
  }


  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::ProcessEvent()
{
  // main method for processing of (selected) event:
  // - Filtering of tracks / particles for flow calculations
  // - Phi,eta,pt weights for generic framework are calculated if specified
  // - Flow calculations
  // returns kTRUE if succesfull
  // *************************************************************

  // printf("======= EVENT ================\n");

  // checking the run number for aplying weights & loading TList with weights
  // if(fFlowUseWeights && fRunNumber < 0)
  if(fFlowUseWeights && (fRunNumber < 0 || fRunNumber != fEventAOD->GetRunNumber()) )
  {
    fRunNumber = fEventAOD->GetRunNumber();
    if(fFlowWeightsFile)
    {
      TList* listFlowWeights = (TList*) fFlowWeightsFile->Get(Form("%d",fRunNumber));
      if(!listFlowWeights) {::Error("ProcessEvent","TList from flow weights not found."); return kFALSE; }
      fh2WeightRefs = (TH2D*) listFlowWeights->FindObject("Refs"); if(!fh2WeightRefs) { ::Error("ProcessEvent","Refs weights not found"); return kFALSE; }
      fh2WeightCharged = (TH2D*) listFlowWeights->FindObject("Charged"); if(!fh2WeightCharged) { ::Error("ProcessEvent","Charged weights not found"); return kFALSE; }
      fh2WeightPion = (TH2D*) listFlowWeights->FindObject("Pion"); if(!fh2WeightPion) { ::Error("ProcessEvent","Pion weights not found"); return kFALSE; }
      fh2WeightKaon = (TH2D*) listFlowWeights->FindObject("Kaon"); if(!fh2WeightKaon) { ::Error("ProcessEvent","Kaon weights not found"); return kFALSE; }
      fh2WeightProton = (TH2D*) listFlowWeights->FindObject("Proton"); if(!fh2WeightProton) { ::Error("ProcessEvent","Proton weights not found"); return kFALSE; }
      fh2WeightK0s = (TH2D*) listFlowWeights->FindObject("K0s"); if(!fh2WeightK0s) { ::Error("ProcessEvent","K0s weights not found"); return kFALSE; }
      fh2WeightLambda = (TH2D*) listFlowWeights->FindObject("Lambda"); if(!fh2WeightLambda) { ::Error("ProcessEvent","Phi weights not found"); return kFALSE; }
      fh2WeightPhi = (TH2D*) listFlowWeights->FindObject("Phi"); if(!fh2WeightPhi) { ::Error("ProcessEvent","Phi weights not found"); return kFALSE; }
    }
  }

  // filtering particles
  Filtering();
  // at this point, centrality index (percentile) should be properly estimated, if not, skip event
  if(fIndexCentrality < 0) return kFALSE;

  // if running in kFillWeights mode, skip the remaining part
  if(fRunMode == kFillWeights) { fEventCounter++; return kTRUE; }

  // checking if there is at least one charged track selected;
  // if not, event is skipped: unable to compute Reference flow (and thus any differential flow)
  if(fVectorCharged->size() < 1)
    return kFALSE;

  // at this point, all particles fullfiling relevant POIs (and REFs) criteria are filled in TClonesArrays

  // >>>> flow starts here <<<<
  // >>>> Flow a la General Framework <<<<
  for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
  {
    // Reference (pT integrated) flow
    DoFlowRefs(iGap);

    // pT differential
    if(fProcessCharged)
    {
      // charged track flow
      DoFlowCharged(iGap);
    }

    if(fProcessPID)
    {
      const Int_t iSizePion = fVectorPion->size();
      const Int_t iSizeKaon = fVectorKaon->size();
      const Int_t iSizeProton = fVectorProton->size();

      // pi,K,p flow
      if(iSizePion > 0) DoFlowPID(iGap,kPion);
      if(iSizeKaon > 0) DoFlowPID(iGap,kKaon);
      if(iSizeProton > 0) DoFlowPID(iGap,kProton);
    }

    if(fProcessPhi)
    {
      const Int_t iSizePhi = fVectorPhi->size();
      if(iSizePhi > 0)
      {
        for(Short_t iMass(0); iMass < fPhiNumBinsMass; iMass++) DoFlowPhi(iGap,iMass);
      }
    }

    if(fProcessV0s)
    {
      const Int_t iSizeK0s = fVectorK0s->size();
      const Int_t iSizeLambda = fVectorLambda->size();

      if(iSizeK0s > 0)
      {
        for(Short_t iMass(0); iMass < fV0sNumBinsMass; iMass++)
        {
          // V0s (K0s, Lambda/ALambda) flow
          DoFlowV0s(iGap,iMass,kK0s);
        }
      }
      if(iSizeLambda > 0)
      {
        for(Short_t iMass(0); iMass < fV0sNumBinsMass; iMass++)
        {
          // V0s (K0s, Lambda/ALambda) flow
          DoFlowV0s(iGap,iMass,kLambda);
        }
      }
    }
  } // endfor {iGap} eta gaps

  fEventCounter++; // counter of processed events
  //printf("event %d\n",fEventCounter);

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowRefs(const Short_t iEtaGapIndex)
{
  // Estimate <2> for reference flow for all harmonics based on relevant flow vectors
  // *************************************************************

  Float_t dEtaGap = fEtaGap[iEtaGapIndex];
  Short_t iHarmonics = 0;
  Double_t Cn2 = 0;
  TComplex vector = TComplex(0,0,kFALSE);
  Double_t dValue = 999;

  FillRefsVectors(iEtaGapIndex); // filling RFPs (Q) flow vectors

  // // filling mean Qx, Qy
  // if(fFlowFillWeights)
  // {
  //   for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
  //   {
  //     fpMeanQxRefsPos[iEtaGapIndex][iHarm]->Fill(fIndexCentrality, fFlowVecQpos[iHarm][0].Re());
  //     fpMeanQyRefsPos[iEtaGapIndex][iHarm]->Fill(fIndexCentrality, fFlowVecQpos[iHarm][0].Im());
  //     fpMeanQxRefsNeg[iEtaGapIndex][iHarm]->Fill(fIndexCentrality, fFlowVecQneg[iHarm][0].Re());
  //     fpMeanQyRefsNeg[iEtaGapIndex][iHarm]->Fill(fIndexCentrality, fFlowVecQneg[iHarm][0].Im());
  //   }
  // }

  if(dEtaGap == -1) // no gap
  {
    // estimating <2>
    Cn2 = Two(0,0).Re();
    if(Cn2 != 0)
    {
      for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
      {
        iHarmonics = fHarmonics[iHarm];
        vector = Two(iHarmonics,-iHarmonics);
        dValue = vector.Re()/Cn2;
        // printf("Gap (RFPs): %g Harm %d | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,Cn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fIndexCentrality);
        if( TMath::Abs(dValue < 1) )
          fpRefsCor2[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, dValue, Cn2);
      }
    }

    // estimating <4>
    if(fCutFlowDoFourCorrelations)
    {
      Cn2 = Four(0,0,0,0).Re();
      if(Cn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = Four(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics);
          dValue = vector.Re()/Cn2;
          // printf("Gap (RFPs): %g Harm %d | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,Cn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          fpRefsCor4[fIndexSampling][iHarm]->Fill(fIndexCentrality, dValue, Cn2);
        }
      }
    }
  }
  else // with gap
  {
    // estimating <2>
    Cn2 = TwoGap(0,0).Re();
    if(Cn2 != 0)
    {
      for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
      {
        iHarmonics = fHarmonics[iHarm];
        vector = TwoGap(iHarmonics,-iHarmonics);
        dValue = vector.Re()/Cn2;
        // printf("Gap (RFPs): %g Harm %d | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,Cn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fIndexCentrality);
        if( TMath::Abs(dValue < 1) )
          fpRefsCor2[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, dValue, Cn2);

      }
    }
  } // endif {dEtaGap}
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowCharged(const Short_t iEtaGapIndex)
{
  // Estimate <2> for pT diff flow of charged tracks for all harmonics based on relevant flow vectors
  // *************************************************************

  FillPOIsVectors(iEtaGapIndex,kCharged);  // filling POIs (P,S) flow vectors

  const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;

  Float_t dEtaGap = fEtaGap[iEtaGapIndex];
  Short_t iHarmonics = 0;
  Double_t Dn2 = 0;
  TComplex vector = TComplex(0,0,kFALSE);
  Double_t dValue = 999;


  for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
  {
    if(dEtaGap == -1) // no eta gap
    {
      // estimating <2'>
      Dn2 = TwoDiff(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiff(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if(TMath::Abs(dValue < 1))
            fp2ChargedCor2Pos[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
        }
      }

      // estimating <4'>
      if(fCutFlowDoFourCorrelations)
      {
        Dn2 = FourDiff(0,0,0,0,iPt).Re();
        if(Dn2 != 0)
        {
          for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
          {
            iHarmonics = fHarmonics[iHarm];
            vector = FourDiff(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics,iPt);
            dValue = vector.Re()/Dn2;
            // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
            if(TMath::Abs(dValue < 1))
            fp2ChargedCor4[fIndexSampling][iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
          }
        }
      }
    }
    else // eta gap
    {
      // estimating <2'>
      // POIs in positive eta
      Dn2 = TwoDiffGapPos(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapPos(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Pos): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
            fp2ChargedCor2Pos[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
        }
      }

      // POIs in negative eta
      Dn2 = TwoDiffGapNeg(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapNeg(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Neg): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
            fp2ChargedCor2Neg[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
        }
      }
    } // endif {dEtaGap}
  } // endfor {iPt}
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowPID(const Short_t iEtaGapIndex, const PartSpecies species)
{
  // Estimate <2> for pT diff flow of pi/K/p tracks for all harmonics based on relevant flow vectors
  // *************************************************************

  TProfile2D** profile2Pos = 0x0;
  TProfile2D** profile2Neg = 0x0;
  TProfile2D** profile4 = 0x0;

  switch (species)
  {
    case kPion:
      profile2Pos = fp2PionCor2Pos[fIndexSampling][iEtaGapIndex];
      profile2Neg = fp2PionCor2Neg[fIndexSampling][iEtaGapIndex];
      profile4 = fp2PionCor4[fIndexSampling];
      break;

    case kKaon:
      profile2Pos = fp2KaonCor2Pos[fIndexSampling][iEtaGapIndex];
      profile2Neg = fp2KaonCor2Neg[fIndexSampling][iEtaGapIndex];
      profile4 = fp2KaonCor4[fIndexSampling];
      break;

    case kProton:
      profile2Pos = fp2ProtonCor2Pos[fIndexSampling][iEtaGapIndex];
      profile2Neg = fp2ProtonCor2Neg[fIndexSampling][iEtaGapIndex];
      profile4 = fp2ProtonCor4[fIndexSampling];
      break;

    default:
      ::Error("DoFlowPID","Unexpected species! Terminating!");
      return;
  }

  FillPOIsVectors(iEtaGapIndex,species); // Filling POIs vectors

  const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;

  Float_t dEtaGap = fEtaGap[iEtaGapIndex];
  Short_t iHarmonics = 0;
  Double_t Dn2 = 0;
  TComplex vector = TComplex(0,0,kFALSE);
  Double_t dValue = 999;

  // filling POIs (P,S) flow vectors

  for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
  {
    if(dEtaGap == -1) // no eta gap
    {
      // estimating <2'>
      Dn2 = TwoDiff(0,0,iPt).Re();
      if(Dn2 != 0)
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiff(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if(TMath::Abs(dValue < 1))
            profile2Pos[iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
        }

      // estimating <4'>
      if(fCutFlowDoFourCorrelations)
      {
        Dn2 = FourDiff(0,0,0,0,iPt).Re();
        if(Dn2 != 0)
        {
          for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
          {
            iHarmonics = fHarmonics[iHarm];
            vector = FourDiff(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics,iPt);
            dValue = vector.Re()/Dn2;
            // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
            if(TMath::Abs(dValue < 1))
            profile4[iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
          }
        }
      }
    }
    else // eta gap
    {
      // estimating <2'>
      // POIs in positive eta
      Dn2 = TwoDiffGapPos(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapPos(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Pos): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
            profile2Pos[iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
        }
      }
      // POIs in negative eta
      Dn2 = TwoDiffGapNeg(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapNeg(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Neg): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
            profile2Neg[iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
        }
      }
    } // endif {dEtaGap}
  } // endfor {iPt}
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowPhi(const Short_t iEtaGapIndex, const Short_t iMassIndex)
{
  // Estimate the correlations for pT diff flow of phi candidates for all harmonics based on relevant flow vectors
  // *************************************************************

  FillPOIsVectors(iEtaGapIndex,kPhi,iMassIndex); // filling POIs (P,S) flow vectors

  const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;
  Float_t dEtaGap = fEtaGap[iEtaGapIndex];
  Short_t iHarmonics = 0;
  Double_t Dn2 = 0;
  TComplex vector = TComplex(0,0,kFALSE);
  Double_t dValue = 999;
  TProfile3D* prof = 0x0;
  Double_t dMassLow = fCutPhiInvMassMin + iMassIndex*(fCutPhiInvMassMax - fCutPhiInvMassMin)/fPhiNumBinsMass;
  Double_t dMassHigh = fCutPhiInvMassMin + (iMassIndex+1)*(fCutPhiInvMassMax - fCutPhiInvMassMin)/fPhiNumBinsMass;
  Double_t dMass = (dMassLow + dMassHigh)/2;

  for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
  {
    if(dEtaGap == -1) // no eta gap
    {
      // estimating <2'>
      Dn2 = TwoDiff(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          prof = fp3PhiCorrCor2Pos[iEtaGapIndex][iHarm];
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiff(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          {
            prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
            // printf("EtaGap %g (no | <2> | n %d) : Filling (mass %g | pt %g | value %g | weight %g) centrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
          }
        }
      }

      // estimating <4'>
      if(fCutFlowDoFourCorrelations)
      {
        Dn2 = FourDiff(0,0,0,0,iPt).Re();
        if(Dn2 != 0)
        {
          for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
          {
            prof = fp3PhiCorrCor4[iHarm];
            iHarmonics = fHarmonics[iHarm];
            vector = FourDiff(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics,iPt);
            dValue = vector.Re()/Dn2;
            // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
            if( TMath::Abs(dValue < 1) )
            {
              prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
              // printf("EtaGap %g (no | <4> | n %d) : Filling (mass %g | pt %g | value %g | weight %g) centrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
            }
          }
        }
      }
    }
    else // eta gap
    {
      // POIs in positive eta
      Dn2 = TwoDiffGapPos(0,0,iPt).Re();

      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          prof = fp3PhiCorrCor2Pos[iEtaGapIndex][iHarm];
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapPos(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Pos): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          {
            prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
            // printf("EtaGap %g (pos |  n %d): Filling (mass %g | pt %g | value %g | weight %g) centrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
          }
        }
      }

      // POIs in negative eta
      Dn2 = TwoDiffGapNeg(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          prof = fp3PhiCorrCor2Neg[iEtaGapIndex][iHarm];
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapNeg(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Neg): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          {
            prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
            // printf("EtaGap %g (neg <2>|  n %d): Filling (mass %g | pt %g | value %g | weight %g) cetrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
          }
        }
      }
    } // endif {dEtaGap}
  } // endfor {iPt}

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::DoFlowV0s(const Short_t iEtaGapIndex, const Short_t iMassIndex, const PartSpecies species)
{
  // Estimate <2> for pT diff flow of V0s for all harmonics based on relevant flow vectors
  // *************************************************************

  FillPOIsVectors(iEtaGapIndex,species,iMassIndex);  // filling POIs (P,S) flow vectors

  const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;
  Float_t dEtaGap = fEtaGap[iEtaGapIndex];
  Short_t iHarmonics = 0;
  Double_t Dn2 = 0;
  TComplex vector = TComplex(0,0,kFALSE);
  Double_t dValue = 999;
  TProfile3D* prof = 0x0;
  Double_t dMassLow = 0;
  Double_t dMassHigh = 0;
  // TProfile3D** profArr = 0x0;
  // TProfile3D** profArr4 = 0x0;

  // switch based on particle species
  switch (species)
  {
    case kK0s:
      // profArr = fp3V0sCorrK0sCor2[iEtaGapIndex];
      // profArr4 = fp3V0sCorrK0sCor4;
      dMassLow = fCutV0sInvMassK0sMin + iMassIndex*(fCutV0sInvMassK0sMax - fCutV0sInvMassK0sMin)/fV0sNumBinsMass;
      dMassHigh = fCutV0sInvMassK0sMin + (iMassIndex+1)*(fCutV0sInvMassK0sMax - fCutV0sInvMassK0sMin)/fV0sNumBinsMass;
      break;

    case kLambda:
      // profArr = fp3V0sCorrLambdaCor2[iEtaGapIndex];
      // profArr4 = fp3V0sCorrLambdaCor4;
      dMassLow = fCutV0sInvMassLambdaMin + iMassIndex*(fCutV0sInvMassLambdaMax - fCutV0sInvMassLambdaMin)/fV0sNumBinsMass;
      dMassHigh = fCutV0sInvMassLambdaMin + (iMassIndex+1)*(fCutV0sInvMassLambdaMax - fCutV0sInvMassLambdaMin)/fV0sNumBinsMass;
      break;

    default:
      ::Error("DoFlowV0s","Selected particles are not K0s nor Lambdas!");
      return;
  }

  Double_t dMass = (dMassHigh+dMassLow) / 2;

  for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
  {
    if(dEtaGap == -1) // no eta gap
    {
      // estimating <2'>
      Dn2 = TwoDiff(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          switch (species)
          {
            case kK0s:
              prof = fp3V0sCorrK0sCor2Pos[iEtaGapIndex][iHarm];
              break;

            case kLambda:
              prof = fp3V0sCorrLambdaCor2Pos[iEtaGapIndex][iHarm];
              break;

            default: return;
          }

          // prof = profArr[iHarm];
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiff(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          {
            prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
            // printf("EtaGap %g (no | <2> | n %d) : Filling (mass %g | pt %g | value %g | weight %g) centrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
          }
        }
      }

      // estimating <4'>
      if(fCutFlowDoFourCorrelations)
      {
        Dn2 = FourDiff(0,0,0,0,iPt).Re();
        if(Dn2 != 0)
        {
          for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
          {
            switch (species)
            {
              case kK0s:
              prof = fp3V0sCorrK0sCor4[iHarm];
              break;

              case kLambda:
              prof = fp3V0sCorrLambdaCor4[iHarm];
              break;

              default: return;
            }

            iHarmonics = fHarmonics[iHarm];
            vector = FourDiff(iHarmonics,iHarmonics,-iHarmonics,-iHarmonics,iPt);
            dValue = vector.Re()/Dn2;
            // printf("Gap (no): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt, Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
            if( TMath::Abs(dValue < 1) )
            {
              prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
              // printf("EtaGap %g (no | <4> | n %d) : Filling (mass %g | pt %g | value %g | weight %g) centrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
            }
          }
        }
      }
    }
    else // eta gap
    {
      // POIs in positive eta
      Dn2 = TwoDiffGapPos(0,0,iPt).Re();

      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          switch (species)
          {
            case kK0s:
              prof = fp3V0sCorrK0sCor2Pos[iEtaGapIndex][iHarm];
              break;

            case kLambda:
              prof = fp3V0sCorrLambdaCor2Pos[iEtaGapIndex][iHarm];
              break;

            default: return;
          }
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapPos(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Pos): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          {
            prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
            // printf("EtaGap %g (pos |  n %d): Filling (mass %g | pt %g | value %g | weight %g) centrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
          }
        }
      }

      // POIs in negative eta
      Dn2 = TwoDiffGapNeg(0,0,iPt).Re();
      if(Dn2 != 0)
      {
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
          switch (species)
          {
            case kK0s:
              prof = fp3V0sCorrK0sCor2Neg[iEtaGapIndex][iHarm];
              break;

            case kLambda:
              prof = fp3V0sCorrLambdaCor2Neg[iEtaGapIndex][iHarm];
              break;

            default: return;
          }
          iHarmonics = fHarmonics[iHarm];
          vector = TwoDiffGapNeg(iHarmonics,-iHarmonics,iPt);
          dValue = vector.Re()/Dn2;
          // printf("Gap (Neg): %g Harm %d Pt %g | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fFlowVecPpos[0][0]: %g | fFlowVecPneg[0][0]: %g | fFlowVecS[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,dPt,Dn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fFlowVecPpos[0][0].Re(),fFlowVecPneg[0][0].Re(),fFlowVecS[0][0].Re(),fIndexCentrality);
          if( TMath::Abs(dValue < 1) )
          {
            prof->Fill(fIndexCentrality, iPt*dPtBinWidth, dMass, dValue, Dn2);
            // printf("EtaGap %g (neg <2>|  n %d): Filling (mass %g | pt %g | value %g | weight %g) cetrality %d\n",dEtaGap,iHarmonics,dMass,iPt*dPtBinWidth,dValue,Dn2,fIndexCentrality);
          }
        }
      }
    } // endif {dEtaGap}
  } // endfor {iPt}
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillRefsVectors(const Short_t iEtaGapIndex)
{
  // Filling Q flow vector with RFPs
  // return kTRUE if succesfull (i.e. no error occurs), kFALSE otherwise
  // *************************************************************
  const Float_t dEtaGap = fEtaGap[iEtaGapIndex];
  TH2D* h2Weights = 0x0;
  Double_t dWeight = 1.;
  if(fFlowUseWeights)
  {
    h2Weights = fh2WeightRefs;
    if(!h2Weights) { ::Error("FillRefsVectors","Histogtram with weights not found."); return; }
  }

  // clearing output (global) flow vectors
  ResetRFPsVector(fFlowVecQpos);
  ResetRFPsVector(fFlowVecQneg);

  Double_t dQcosPos, dQcosNeg, dQsinPos, dQsinNeg;

  // Double_t dQcosPos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
  // Double_t dQcosNeg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
  // Double_t dQsinPos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
  // Double_t dQsinNeg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};

  for (auto part = fVectorCharged->begin(); part != fVectorCharged->end(); part++)
  {
    // checking species of used particles (just for double checking purpose)
    if( part->species != kCharged)
    {
      ::Warning("FillRefsVectors","Unexpected part. species (%d) in selected sample (expected %d)",part->species,kCharged);
      continue;
    }

    // RFPs pT check
    if(fCutFlowRFPsPtMin > 0. && part->pt < fCutFlowRFPsPtMin)
      continue;

    if(fCutFlowRFPsPtMax > 0. && part->pt > fCutFlowRFPsPtMax)
      continue;

    // 0-ing variables
    dQcosPos = 0;
    dQcosNeg = 0;
    dQsinPos = 0;
    dQsinNeg = 0;

    // loading weights if needed
    if(fFlowUseWeights && h2Weights)
    {
      dWeight = h2Weights->GetBinContent(h2Weights->FindBin(part->eta,part->phi));
      if(dWeight <= 0) dWeight = 1.;
      // if(iEtaGapIndex == 0) fh3AfterWeightsRefs->Fill(part->phi,part->eta,part->pt, dWeight);
    }

    // RPF candidate passing all criteria: start filling flow vectors

    if(dEtaGap == -1) // no eta gap
    {
      for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
          dQcosPos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
          dQsinPos = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dQcosPos,dQsinPos,kFALSE);
        } // endfor {iPower}
    }
    else
    {
      if(part->eta > dEtaGap / 2 )
      {
        // RFP in positive eta acceptance
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            dQcosPos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
            dQsinPos = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dQcosPos,dQsinPos,kFALSE);
          }
      }
      if(part->eta < -dEtaGap / 2 )
      {
        // RFP in negative eta acceptance
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            dQcosNeg = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
            dQsinNeg = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dQcosNeg,dQsinNeg,kFALSE);
          }
      }
    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  // // filling local flow vectors to global flow vector arrays
  // for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  //   for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
  //   {
  //     fFlowVecQpos[iHarm][iPower] = TComplex(dQcosPos[iHarm][iPower],dQsinPos[iHarm][iPower],kFALSE);
  //     if(dEtaGap > -1)
  //       fFlowVecQneg[iHarm][iPower] = TComplex(dQcosNeg[iHarm][iPower],dQsinNeg[iHarm][iPower],kFALSE);
  //   }

  // printf("RFPs EtaGap %g : number %g (pos) %g (neg) \n", dEtaGap,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re());
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillPOIsVectors(const Short_t iEtaGapIndex, const PartSpecies species, const Short_t iMassIndex)
{
  // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
  // *************************************************************

  if(species == kUnknown) return;

  Double_t dWeight = 1.;  // for generic framework != 1

  // clearing output (global) flow vectors
  ResetPOIsVector(fFlowVecPpos);
  ResetPOIsVector(fFlowVecPneg);
  ResetPOIsVector(fFlowVecS);

  std::vector<FlowPart>* vector = 0x0;
  TH3D* hist = 0x0;
  Double_t dMassLow = 0, dMassHigh = 0;
  TH2D* h2Weights = 0x0;

  // swich based on species
  switch (species)
  {
    case kCharged:
      vector = fVectorCharged;
      if(fFlowUseWeights) { h2Weights = fh2WeightCharged; }
      break;

    case kPion:
      vector = fVectorPion;
      if(fFlowUseWeights) { h2Weights = fh2WeightPion; }
      break;

    case kKaon:
      vector = fVectorKaon;
      if(fFlowUseWeights) { h2Weights = fh2WeightKaon; }

      break;

    case kProton:
      vector = fVectorProton;
      if(fFlowUseWeights) { h2Weights = fh2WeightProton; }
      break;

    case kK0s:
      vector = fVectorK0s;
      hist = fh3V0sEntriesK0s[iEtaGapIndex];
      if(fFlowUseWeights) { h2Weights = fh2WeightK0s; }
      dMassLow = fCutV0sInvMassK0sMin + iMassIndex*(fCutV0sInvMassK0sMax - fCutV0sInvMassK0sMin)/fV0sNumBinsMass;
      dMassHigh = fCutV0sInvMassK0sMin + (iMassIndex+1)*(fCutV0sInvMassK0sMax - fCutV0sInvMassK0sMin)/fV0sNumBinsMass;
      break;

    case kLambda: // if a Lambda/ALambda candidates: first go through Lambda array and then goes through ALambda array
      vector = fVectorLambda;
      hist = fh3V0sEntriesLambda[iEtaGapIndex];
      if(fFlowUseWeights) { h2Weights = fh2WeightLambda; }
      dMassLow = fCutV0sInvMassLambdaMin + iMassIndex*(fCutV0sInvMassLambdaMax - fCutV0sInvMassLambdaMin)/fV0sNumBinsMass;
      dMassHigh = fCutV0sInvMassLambdaMin + (iMassIndex+1)*(fCutV0sInvMassLambdaMax - fCutV0sInvMassLambdaMin)/fV0sNumBinsMass;
      break;

    case kPhi:
      vector = fVectorPhi;
      hist = fh3PhiEntriesSignal[iEtaGapIndex];
      if(fFlowUseWeights) { h2Weights = fh2WeightPhi; }
      dMassLow = fCutPhiInvMassMin + iMassIndex*(fCutPhiInvMassMax - fCutPhiInvMassMin)/fPhiNumBinsMass;
      dMassHigh = fCutPhiInvMassMin + (iMassIndex+1)*(fCutPhiInvMassMax - fCutPhiInvMassMin)/fPhiNumBinsMass;
      break;

    default:
      ::Error("FillPOIsVectors","Selected species unknown.");
      return;
  }

  if(fFlowUseWeights && !h2Weights) { ::Error("FillPOIsVectors","Histogtram with weights not found."); return; }

  const Double_t dEtaGap = fEtaGap[iEtaGapIndex];
  const Double_t dMass = (dMassLow+dMassHigh)/2;

  Short_t iPtBin = 0;
  Double_t dCos = 0, dSin = 0;

  for (auto part = vector->begin(); part != vector->end(); part++)
  {
    // checking species of used particles (just for double checking purpose)
    if( part->species != species)
    {
      ::Warning("FillPOIsVectors","Unexpected part. species (%d) in selected sample (expected %d)",part->species,species);
      continue;
    }

    // POIs mass bin check for V0s candidates
    if(species == kK0s || species == kLambda || species == kPhi)
    {
      if(part->mass < dMassLow || part->mass >= dMassHigh)
        continue;
    }

    // assign iPtBin based on particle momenta
    iPtBin = GetPOIsPtBinIndex(part->pt);

    // 0-ing variables
    dCos = 0;
    dSin = 0;

    // POIs candidate passing all criteria: start filling flow vectors

    // loading weights if needed
    if(fFlowUseWeights && h2Weights)
    {
      dWeight = h2Weights->GetBinContent(h2Weights->FindBin(part->eta,part->phi));
      if(dWeight <= 0) dWeight = 1.;
    }

    if(dEtaGap == -1) // no eta gap
    {
      if(species == kK0s || species == kLambda || species == kPhi)
      {
        hist->Fill(fIndexCentrality,part->pt,dMass,1);
      }

      for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
          dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
          dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
          fFlowVecPpos[iHarm][iPower][iPtBin] += TComplex(dCos,dSin,kFALSE);

          // check if track (passing criteria) is overlapping with RFPs pT region; if so, fill S (q) vector
          if(species == kCharged || species == kPion || species == kKaon || species == kProton)
          {
            if(part->pt > fCutFlowRFPsPtMin && part->pt < fCutFlowRFPsPtMax)
            {
              dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
              dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
              fFlowVecS[iHarm][iPower][iPtBin] += TComplex(dCos,dSin,kFALSE);
            }
          }
        }

    }
    else // with eta gap
    {
      if(part->eta > dEtaGap / 2 )
      {
        // particle in positive eta acceptance
        if(species == kK0s || species == kLambda || species == kPhi)
        {
          hist->Fill(fIndexCentrality,part->pt,dMass,1);
        }

        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
            dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
            fFlowVecPpos[iHarm][iPower][iPtBin] += TComplex(dCos,dSin,kFALSE);
          }
       }
       if(part->eta < -dEtaGap / 2 )
       {
         // particle in negative eta acceptance
         if(hist && (species == kK0s || species == kLambda || species == kPhi))
         {
           hist->Fill(fIndexCentrality,part->pt,dMass,1);
         }

         for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
           for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
           {
             dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * part->phi);
             dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * part->phi);
             fFlowVecPneg[iHarm][iPower][iPtBin] += TComplex(dCos,dSin,kFALSE);
           }
       }
     } // endif {dEtaGap}
   } // endfor {tracks}
   return;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetPOIsPtBinIndex(const Double_t pt)
{
  // Return POIs pT bin index based on pT value
  // *************************************************************
  const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;
  // printf("Pt %g | index %d\n",pt,(Short_t) (pt / dPtBinWidth) );
  return (Short_t) (pt / dPtBinWidth);
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ResetRFPsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
      array[iHarm][iPower] = TComplex(0,0,kFALSE);
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ResetPOIsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins])
{
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
      for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
        array[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
  // List all values of given flow vector TComplex array
  // *************************************************************
  printf(" ### Listing (TComplex) flow vector array ###########################\n");
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    printf("Harm %d (power):",iHarm);
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
        printf("|(%d) %g+%g(i)",iPower, array[iHarm][iPower].Re(), array[iHarm][iPower].Im());
    }
    printf("\n");
  }
  return;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetSamplingIndex()
{
  // Assessing sampling index based on generated random number
  // returns centrality index
  // *************************************************************

  Short_t index = 0x0;

  if(fSampling && fNumSamples > 1)
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
  // Estimating centrality percentile based on selected estimator.
  // (Default) If no multiplicity estimator is specified (fMultEstimator == '' || Charged), percentile is estimated as number of selected / filtered charged tracks.
  // If a valid multiplicity estimator is specified, centrality percentile is estimated via AliMultSelection
  // otherwise -1 is returned (and event is skipped)
  // *************************************************************
  fMultEstimator.ToUpper();

  if(
      fMultEstimator.EqualTo("V0A") || fMultEstimator.EqualTo("V0C") || fMultEstimator.EqualTo("V0M") ||
      fMultEstimator.EqualTo("CL0") || fMultEstimator.EqualTo("CL1") ||
      fMultEstimator.EqualTo("ZNA") || fMultEstimator.EqualTo("ZNC")
    )
  {
    // some of supported AliMultSelection estimators (listed above)
    Float_t dPercentile = 300;

    // checking AliMultSelection
    AliMultSelection* multSelection = 0x0;
    multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1;}

    dPercentile = multSelection->GetMultiplicityPercentile(fMultEstimator.Data());
    if(dPercentile > 100 || dPercentile < 0) { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1;}
    else {return dPercentile;}
  }
  else if(fMultEstimator.EqualTo("") || fMultEstimator.EqualTo("CHARGED"))
  {
    // assigning centrality based on number of selected charged tracks
    return fVectorCharged->size();
  }
  else
  {
    AliWarning(Form("Multiplicity estimator '%s' not supported. Returning -1\n",fMultEstimator.Data()));
    return -1;
  }

  return -1;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::Terminate(Option_t* option)
{
  // called on end of task, after all events are processed
  // *************************************************************

  return;
}
//_____________________________________________________________________________
// Set of methods returning given complex flow vector based on flow harmonics (n) and weight power indexes (p)
// a la General Framework implementation.
// Q: flow vector of RFPs (with/out eta gap)
// P: flow vector of POIs (with/out eta gap) (in usual notation p)
// S: flow vector of overlaping RFPs and POIs (in usual notation q)

TComplex AliAnalysisTaskUniFlow::Q(const Short_t n, const Short_t p)
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapPos(const Short_t n, const Short_t p)
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapNeg(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
  else return fFlowVecQneg[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::P(const Short_t n, const Short_t p, const Short_t pt)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p][pt]);
  else return fFlowVecPpos[n][p][pt];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::PGapPos(const Short_t n, const Short_t p, const Short_t pt)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p][pt]);
  else return fFlowVecPpos[n][p][pt];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::PGapNeg(const Short_t n, const Short_t p, const Short_t pt)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPneg[-n][p][pt]);
  else return fFlowVecPneg[n][p][pt];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::S(const Short_t n, const Short_t p, const Short_t pt)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecS[-n][p][pt]);
  else return fFlowVecS[n][p][pt];
}
//____________________________________________________________________

// Set of flow calculation methods for cumulants of different orders with/out eta gap

TComplex AliAnalysisTaskUniFlow::Two(const Short_t n1, const Short_t n2)
{
  TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoGap(const Short_t n1, const Short_t n2)
{
  TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiff(const Short_t n1, const Short_t n2, const Short_t pt)
{
  TComplex formula = P(n1,1,pt)*Q(n2,1) - S(n1+n2,1,pt);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiffGapPos(const Short_t n1, const Short_t n2, const Short_t pt)
{
  TComplex formula = PGapPos(n1,1,pt)*QGapNeg(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t pt)
{
  TComplex formula = PGapNeg(n1,1,pt)*QGapPos(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
  return formula;
}
// //____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t pt)
{
  TComplex formula = P(n1,1,pt)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2,pt)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2,pt)*Q(n4,1)
                    - P(n1,1,pt)*Q(n2+n3,2)*Q(n4,1)+2.*S(n1+n2+n3,3,pt)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2,pt)
                    + Q(n2+n3,2)*S(n1+n4,2,pt)-P(n1,1,pt)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2,pt)*Q(n2+n4,2)
                    + 2.*Q(n3,1)*S(n1+n2+n4,3,pt)-P(n1,1,pt)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2,pt)*Q(n3+n4,2)
                    + 2.*Q(n2,1)*S(n1+n3+n4,3,pt)+2.*P(n1,1,pt)*Q(n2+n3+n4,3)-6.*S(n1+n2+n3+n4,4,pt);
    return formula;
}
//____________________________________________________________________
// TComplex* AliAnalysisTaskUniFlow::FourGap(int n1, int n2, int n3, int n4)
// {
//
//   TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
//                     - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
//                     + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
//                     + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
//                     + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
//   TComplex *out = (TComplex*) &formula;
//   return out;
//
// }
//____________________________________________________________________
AliAnalysisTaskUniFlow::FlowPart AliAnalysisTaskUniFlow::FlowPart::MakeMother(const FlowPart* part1, const FlowPart* part2, const PartSpecies species)
{
  // Reconstructing mother particle from two daughters (FlowPart type) and fill its properties.
  // return created mother particle
  // *************************************************************

  if(!part1 || !part2) return FlowPart();
  // printf(" part1: "); part1->PrintPart();
  // printf(" part2: "); part2->PrintPart();

  if(species == kPhi && (part1->species != kKaon || part2->species != kKaon)) return FlowPart();

  TVector3 vec1 = TVector3(part1->px, part1->py, part1->pz);
  TVector3 vec2 = TVector3(part2->px, part2->py, part2->pz);
  TVector3 mother = vec1 + vec2;

  Double_t dMass = InvMass(part1,part2);
  Short_t iCharge = part1->charge + part2->charge;

  return FlowPart(mother.Pt(), mother.Phi(), mother.PseudoRapidity(), iCharge, species, dMass, mother.Px(), mother.Py(), mother.Pz());
}
//____________________________________________________________________
Double_t AliAnalysisTaskUniFlow::FlowPart::InvMass(const FlowPart* part1, const FlowPart* part2)
{
  // Calculating invariant mass of potential mother from two daughter particles given their properties and species.
  // return calculated invariant mass or -999 if an error occurs.
  // *************************************************************

  // checking if both poiters are valid
  if(!part1 || !part2) return -999;

  // checking daughter properties
  if(part1->species == kUnknown || part2->species == kUnknown) return -999;

  Double_t dP1Sq = TMath::Power(part1->px,2) + TMath::Power(part1->py,2) + TMath::Power(part1->pz,2);
  Double_t dP2Sq = TMath::Power(part2->px,2) + TMath::Power(part2->py,2) + TMath::Power(part2->pz,2);
  Double_t dE1 = TMath::Sqrt( dP1Sq + TMath::Power(part1->mass,2) );
  Double_t dE2 = TMath::Sqrt( dP2Sq + TMath::Power(part2->mass,2) );
  Double_t dMassSq = TMath::Power((dE1+dE2),2) - ( TMath::Power((part1->px + part2->px),2) + TMath::Power((part1->py + part2->py),2) + TMath::Power((part1->pz + part2->pz),2) );

  if(dMassSq < 0) return -999;
  return TMath::Sqrt(dMassSq);
}
//____________________________________________________________________
