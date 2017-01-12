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
#include "AliAnalysisTaskFlowPID.h"
#include "AliLog.h"

class AliAnalysisTaskFlowPID;

ClassImp(AliAnalysisTaskFlowPID) // classimp: necessary for root

Double_t AliAnalysisTaskFlowPID::fPtBinEdges[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0}; // You, Katarina binning
//Double_t AliAnalysisTaskFlowPID::fPtBinEdges[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.,5.5,6}; // PID flow v2 JHEP paper

//Double_t AliAnalysisTaskFlowPID::fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.}; # PbPb
Double_t AliAnalysisTaskFlowPID::fCentBinEdges[] = {0.,30.,50.,150.};
//Double_t AliAnalysisTaskFlowPID::fMinvFlowBinEdgesK0s[] = {0.4,0.42,0.44,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.54,0.56,0.58,0.6};
Double_t AliAnalysisTaskFlowPID::fMinvFlowBinEdgesK0s[] = {0.4,0.425,0.45,0.47,0.49,0.495,0.5,0.505,0.51,0.53,0.55,0.575,0.6};
Double_t AliAnalysisTaskFlowPID::fMinvFlowBinEdgesLambda[] = {1.08,1.09,1.10,1.105,1.11,1.115,1.12,1.125,1.13,1.14,1.15,1.16};
Int_t AliAnalysisTaskFlowPID::fHarmonics[AliAnalysisTaskFlowPID::fNumHarmonics] = {2,3,4};
Double_t AliAnalysisTaskFlowPID::fEtaGap[AliAnalysisTaskFlowPID::fNumEtaGap] = {-1.,0.,0.5,1.0};
Short_t AliAnalysisTaskFlowPID::fTracksScanFB[AliAnalysisTaskFlowPID::fNumScanFB] = {1,2,16,32,64,96,128,256,512,768,1024};

AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID() : AliAnalysisTaskSE(),
  fAOD(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fTPCPIDResponse(),
  fEtaCutFlag(0),
  fHarmFlag(0),
  fV0candK0s(0),
  fV0candLambda(0),
  fV0candALambda(0),
  fV0MinMassK0s(fMinvFlowBinEdgesK0s[0]),
  fV0MaxMassK0s(fMinvFlowBinEdgesK0s[fNumMinvFlowBinsK0s]),
  fV0MinMassLambda(fMinvFlowBinEdgesLambda[0]),
  fV0MaxMassLambda(fMinvFlowBinEdgesLambda[fNumMinvFlowBinsLambda]),
  fCentBinIndex(0),
  fCentPercentile(0),
  fPtBinIndex(0),
  fMinvFlowBinIndex(0),
  fSampleBinIndex(0),
  fVecRefPos(0,0,kFALSE),
  fVecRefNeg(0,0,kFALSE),
  fCountRefPos(0),
  fCountRefNeg(0),
  fArrTracksFiltered("AliAODTrack",10000),
  fArrPionFiltered("AliAODTrack",5000),
  fArrKaonFiltered("AliAODTrack",5000),
  fArrProtonFiltered("AliAODTrack",5000),
  fArrV0sK0sFiltered("AliAODv0",5000),
  fArrV0sLambdaFiltered("AliAODv0",5000),
  fArrV0sALambdaFiltered("AliAODv0",5000),

  fOutListCumulants(0),
  fOutListEvents(0),
  fOutListTracks(0),
  fOutListPID(0),
  fOutListV0s(0),
  //fOutListQA(0),

  fAODAnalysis(kTRUE),
  fPbPb(kFALSE),
  fPP(kFALSE),
  fPPb(kFALSE),
  fLHC10h(kTRUE),
  fTrigger(kFALSE),
  fRejectPileFromSPD(kFALSE),
  fUseIsPileUpFromSPD(kFALSE),
  fUsePlpMV(kFALSE),
  fRejectOutOfBunchPU(kFALSE),
  fCentFlag(0),
  fSampling(0),
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

  fEventCounter(0),
  fSampleCounter(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentDistUnitBin(0),
  fCentSPDvsV0M(0),

  fTracksCounter(0x0),
  fTracksScanCounter(0x0),
  fTracksScanPt(0x0),
  fTracksScanPhi(0x0),
  fTracksScanEta(0x0),

  fPionsCounter(0x0),
  fKaonsCounter(0x0),
  fProtonsCounter(0x0),
  fPionsMult(0x0),
  fPionsPt(0x0),
  fPionsEta(0x0),
  fPionsPhi(0x0),
  fPionsTPCdEdx(0x0),
  fPionsTOFbeta(0x0),
  fPionsNsigmasTPCTOF(0x0),
  fPionsBayesAsPion(0x0),
  fPionsBayesAsKaon(0x0),
  fPionsBayesAsProton(0x0),
  fPionsBayes(0x0),
  fKaonsMult(0x0),
  fKaonsPt(0x0),
  fKaonsEta(0x0),
  fKaonsPhi(0x0),
  fKaonsTPCdEdx(0x0),
  fKaonsTOFbeta(0x0),
  fKaonsNsigmasTPCTOF(0x0),
  fKaonsBayesAsPion(0x0),
  fKaonsBayesAsKaon(0x0),
  fKaonsBayesAsProton(0x0),
  fKaonsBayes(0x0),
  fProtonsMult(0x0),
  fProtonsPt(0x0),
  fProtonsEta(0x0),
  fProtonsPhi(0x0),
  fProtonsTPCdEdx(0x0),
  fProtonsTOFbeta(0x0),
  fProtonsNsigmasTPCTOF(0x0),
  fProtonsBayesAsPion(0x0),
  fProtonsBayesAsKaon(0x0),
  fProtonsBayesAsProton(0x0),
  fProtonsBayes(0x0),
  fQAPIDTOFbetaNoTOF(0x0),

  fQAV0sCounter(0),
  fQAV0sCounterK0s(0),
  fQAV0sCounterLambda(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fTPCPIDResponse(),
  fEtaCutFlag(0),
  fHarmFlag(0),

  fV0candK0s(0),
  fV0candLambda(0),
  fV0candALambda(0),
  fV0MinMassK0s(fMinvFlowBinEdgesK0s[0]),
  fV0MaxMassK0s(fMinvFlowBinEdgesK0s[fNumMinvFlowBinsK0s]),
  fV0MinMassLambda(fMinvFlowBinEdgesLambda[0]),
  fV0MaxMassLambda(fMinvFlowBinEdgesLambda[fNumMinvFlowBinsLambda]),

  fCentBinIndex(0),
  fCentPercentile(0),
  fPtBinIndex(0),
  fMinvFlowBinIndex(0),
  fSampleBinIndex(0),

  fVecRefPos(0),
  fVecRefNeg(0),
  fCountRefPos(0),
  fCountRefNeg(0),
  fArrTracksFiltered("AliAODTrack",10000),
  fArrPionFiltered("AliAODTrack",5000),
  fArrKaonFiltered("AliAODTrack",5000),
  fArrProtonFiltered("AliAODTrack",5000),
  fArrV0sK0sFiltered("AliAODv0",5000),
  fArrV0sLambdaFiltered("AliAODv0",5000),
  fArrV0sALambdaFiltered("AliAODv0",5000),

  fOutListCumulants(0),
  fOutListEvents(0),
  fOutListTracks(0),
  fOutListPID(0),
  fOutListV0s(0),
  //fOutListQA(0),

  fAODAnalysis(kTRUE),
  fPbPb(kFALSE),
  fPP(kFALSE),
  fPPb(kFALSE),
  fLHC10h(kTRUE),
  fTrigger(kFALSE),
  fRejectPileFromSPD(kFALSE),
  fUseIsPileUpFromSPD(kFALSE),
  fRejectOutOfBunchPU(kFALSE),
  fUsePlpMV(kFALSE),
  fCentFlag(0),
  fDoFlow(0),
  fSampling(0),
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

  fEventCounter(0),
  fSampleCounter(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentDistUnitBin(0),
  fCentSPDvsV0M(0),

  fTracksCounter(0x0),
  fTracksScanCounter(0x0),
  fTracksScanPt(0x0),
  fTracksScanPhi(0x0),
  fTracksScanEta(0x0),

  fPionsCounter(0x0),
  fKaonsCounter(0x0),
  fProtonsCounter(0x0),
  fPionsMult(0x0),
  fPionsPt(0x0),
  fPionsEta(0x0),
  fPionsPhi(0x0),
  fPionsTPCdEdx(0x0),
  fPionsTOFbeta(0x0),
  fPionsNsigmasTPCTOF(0x0),
  fPionsBayesAsPion(0x0),
  fPionsBayesAsKaon(0x0),
  fPionsBayesAsProton(0x0),
  fPionsBayes(0x0),
  fKaonsMult(0x0),
  fKaonsPt(0x0),
  fKaonsEta(0x0),
  fKaonsPhi(0x0),
  fKaonsTPCdEdx(0x0),
  fKaonsTOFbeta(0x0),
  fKaonsNsigmasTPCTOF(0x0),
  fKaonsBayesAsPion(0x0),
  fKaonsBayesAsKaon(0x0),
  fKaonsBayesAsProton(0x0),
  fKaonsBayes(0x0),
  fProtonsMult(0x0),
  fProtonsPt(0x0),
  fProtonsEta(0x0),
  fProtonsPhi(0x0),
  fProtonsTPCdEdx(0x0),
  fProtonsTOFbeta(0x0),
  fProtonsNsigmasTPCTOF(0x0),
  fProtonsBayesAsPion(0x0),
  fProtonsBayesAsKaon(0x0),
  fProtonsBayesAsProton(0x0),
  fProtonsBayes(0x0),
  fQAPIDTOFbetaNoTOF(0x0),

  fQAV0sCounter(0),
  fQAV0sCounterK0s(0),
  fQAV0sCounterLambda(0)
{
  // constructor

  for(Int_t i(0); i < fNumEtaGap; i++)
  {
    fListCumRef[i] = 0x0;
    fListCumTracks[i] = 0x0;
    fListCumPions[i] = 0x0;
    fListCumKaons[i] = 0x0;
    fListCumProtons[i] = 0x0;
    fListCumK0s[i] = 0x0;
    fListCumLambda[i] = 0x0;
  }

  // KAtarina gen fram
  for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
  {
    for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
    {
      Qvector[iharm][ipow] = TComplex(0,0,kFALSE);

      for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
      {
        QvectorGap[iVec][iharm][ipow] = TComplex(0,0,kFALSE);
      }

      pvector[iharm][ipow] = TComplex(0,0,kFALSE);
      qvector[iharm][ipow] = TComplex(0,0,kFALSE);

      for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
      {
        pvectorGap[iVec][iharm][ipow] = TComplex(0,0,kFALSE);
        //qvectorGap[iVec][iharm][ipow][iPt] = TComplex(0,0,kFALSE);

      }

    }
  }

  for(Int_t k(0); k < fNumEtaGap; k++)
  {
    for(Int_t i(0); i < fNumHarmonics; i++)
    {
      for(Int_t j(0); j < fNumSampleBins; j++)
      {
        fcn2Tracks[k][i][j] = 0x0;
        fcn4Tracks[k][i][j] = 0x0;
      }
    }
  }
  for(Int_t m(0); m < fNumEtaGap; m++)
  {
    for(Int_t i(0); i < fNumCentBins; i++)
    {
      for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
      {
        fInvMassPtK0s[iVec][m][i] = 0x0;
        fInvMassPtLambda[iVec][m][i] = 0x0;
      }

      for(Int_t k(0); k < fNumHarmonics; k++)
      {
        for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
        {
          fdn2K0s[iVec][m][k][i] = 0x0;
          fdn2Lambda[iVec][m][k][i] = 0x0;
        }

        fdn4K0s[m][k][i] = 0x0;
        fdn4Lambda[m][k][i] = 0x0;

        for(Int_t j(0); j < fNumSampleBins; j++)
        {
          for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
          {
            fdn2Tracks[iVec][m][k][i][j] = 0x0;
            fdn2Pion[iVec][m][k][i][j] = 0x0;
            fdn2Kaon[iVec][m][k][i][j] = 0x0;
            fdn2Proton[iVec][m][k][i][j] = 0x0;
          }

          fdn4Tracks[m][k][i][j] = 0x0;
          fdn4Pion[m][k][i][j] = 0x0;
          fdn4Kaon[m][k][i][j] = 0x0;
          fdn4Proton[m][k][i][j] = 0x0;

        }
      }
    }
  }


  //end of katarina gen fram

  for(Int_t i(0); i < fNumCentBins; i++)
  {
    for(Int_t j(0); j < fNumHarmonics; j++)
    {
      for(Int_t k(0); k < fNumEtaGap; k++)
      {
        fV0sPtInvMassK0s[i][j][k] = 0x0;
        fV0sPtInvMassLambda[i][j][k] = 0x0;

        for(Int_t m(0); m < fNumSampleBins; m++)
        {
          fV0sDiffTwoPos_K0s[i][j][k][m] = 0x0;
          fV0sDiffTwoNeg_K0s[i][j][k][m] = 0x0;
          fV0sDiffTwoPos_Lambda[i][j][k][m] = 0x0;
          fV0sDiffTwoNeg_Lambda[i][j][k][m] = 0x0;

          fTracksDiffTwoPos[i][j][k][m] = 0x0;
          fTracksDiffTwoNeg[i][j][k][m] = 0x0;

          fPionsDiffTwoPos[i][j][k][m] = 0x0;
          fPionsDiffTwoNeg[i][j][k][m] = 0x0;
          fKaonsDiffTwoPos[i][j][k][m] = 0x0;
          fKaonsDiffTwoNeg[i][j][k][m] = 0x0;
          fProtonsDiffTwoPos[i][j][k][m] = 0x0;
          fProtonsDiffTwoNeg[i][j][k][m] = 0x0;
        }
      }
    }
  }

  for(Int_t k(0); k < fNumEtaGap; k++)
  {
    fV0sInvMassK0s[k] = 0x0;
    fV0sInvMassLambda[k] = 0x0;
  }

  for(Int_t i(0); i < fNumHarmonics; i++)
  {
    for(Int_t j(0); j < fNumEtaGap; j++)
    {
      for(Int_t m(0); m < fNumSampleBins; m++)
      {
        fTracksRefTwo[i][j][m] = 0x0;
      }
    }
  }

  // QA plots
  for (Int_t i(0); i < fQANumSteps; i++)
  {
    // Events
    fQAEventsPVz[i] = 0x0;
    fQAEventsNumContrPV[i] = 0x0;
    fQAEventsNumSPDContrPV[i] = 0x0;
    fQAEventsDistPVSPD[i] = 0x0;
    fQAEventsSPDresol[i] = 0x0;
    fQAEventsTriggers[i] = 0x0;
    fQAEventsTriggerSelection[i] = 0x0;
    // Tracks
    fQATracksMult[i] = 0x0;
    fQATracksPt[i] = 0x0;
    fQATracksEta[i] = 0x0;
    fQATracksPhi[i] = 0x0;
    fQATracksFilterMap[i] = 0x0;
    fQATracksFilterMapBit[i] = 0x0;
    fQATracksNumTPCcls[i] = 0x0;
    fQATracksDCAxy[i] = 0x0;
    fQATracksDCAz[i] = 0x0;
    fQATracksTPCstatus[i] = 0x0;
    fQATracksTOFstatus[i] = 0x0;
    fQATracksTPCdEdx[i] = 0x0;
    fQATracksTOFbeta[i] = 0x0;
    fQATracksTOF[i] = 0x0;
    // PID
    fQAPIDTPCdEdx[i] = 0x0;
    fQAPIDTOFbeta[i] = 0x0;
    fQAPIDTOFbetaWithTOF[i] = 0x0;
    fQAPIDNsigmasTPCasPion[i] = 0x0;
    fQAPIDNsigmasTOFasPion[i] = 0x0;
    fQAPIDNsigmasTPCasKaon[i] = 0x0;
    fQAPIDNsigmasTOFasKaon[i] = 0x0;
    fQAPIDNsigmasTPCasProton[i] = 0x0;
    fQAPIDNsigmasTOFasProton[i] = 0x0;
    fQAPIDBayesProbPion[i] = 0x0;
    fQAPIDBayesProbKaon[i] = 0x0;
    fQAPIDBayesProbProton[i] = 0x0;
    fQAPIDBayesProb[i] = 0x0;

    // V0s
  	fQAV0sRecoMethod[i] = 0;
		fQAV0sTPCRefit[i] = 0;
		fQAV0sKinks[i] = 0;
		fQAV0sDCAtoPV[i] = 0;
		fQAV0sDCADaughters[i] = 0;
		fQAV0sDecayRadius[i] = 0;
    fQAV0sDaughterPt[i] = 0;
		fQAV0sDaughterPhi[i] = 0;
		fQAV0sDaughterEta[i] = 0;
    fQAV0sDaughterTPCdEdxPt[i] = 0;
    fQAV0sMotherPt[i] = 0;
		fQAV0sMotherPhi[i] = 0;
		fQAV0sMotherEta[i] = 0;
		fQAV0sMotherRapK0s[i] = 0;
		fQAV0sMotherRapLambda[i] = 0;
    fQAV0sInvMassK0s[i] = 0;
    fQAV0sInvMassLambda[i] = 0;
		fQAV0sCPAK0s[i] = 0;
		fQAV0sCPALambda[i] = 0;
		fQAV0sNumTauK0s[i] = 0;
		fQAV0sNumTauLambda[i] = 0;
		fQAV0sArmenterosK0s[i] = 0;
		fQAV0sArmenterosLambda[i] = 0;
		fQAV0sProtonNumSigmaPtLambda[i] = 0;
  }

  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                      // this chain is created by the analysis manager, so no need to worry about it,
                                      // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                      // you can add more output objects by calling DefineOutput(2, classname::Class())
                                      // if you add more output objects, make sure to call PostData for all of them, and to
                                      // make changes to your AddTask macro!
  DefineOutput(2, TList::Class());

  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  //DefineOutput(6, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::~AliAnalysisTaskFlowPID()
{
  // destructor
  if(fPIDCombined)
  {
    delete fPIDCombined;
  }

  if(fOutListCumulants)
  {
    delete fOutListCumulants;
  }

  if(fOutListEvents)
  {
    delete fOutListEvents;
  }

  if(fOutListTracks)
  {
    delete fOutListTracks;
  }

  if(fOutListPID)
  {
    delete fOutListPID;
  }

  if(fOutListV0s)
  {
    delete fOutListV0s;
  }

  /*
  if(fOutListQA)
  {
    delete fOutListQA;
  }
  */
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you create the histograms that you want to use
  // the histograms are in this case added to a TList, this list is in the end saved to an output file

  ListParameters();

  if((fPP && fPPb) || (fPP && fPbPb) || (fPPb && fPbPb) )
  {
    ::Error("UserCreateOutputObjects","Analysis switch on for more than one colliding system. Terminating!");
    return;
  }

  if(fUseBayesPID)
  {
    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();
    fPIDCombined->SetSelectedSpecies(5);

    if(fTrackFilterBit == 2)
    {
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetITS);
    }
    else
    {
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask
    }
  }


  fOutListCumulants = new TList();
  fOutListCumulants->SetOwner(kTRUE);

  fOutListEvents = new TList();
  fOutListEvents->SetOwner(kTRUE);

  fOutListTracks = new TList();
  fOutListTracks->SetOwner(kTRUE);

  fOutListPID = new TList();
  fOutListPID->SetOwner(kTRUE);

  fOutListV0s = new TList();
  fOutListV0s->SetOwner(kTRUE);

  //fOutListQA = new TList();
  //fOutListQA->SetOwner(kTRUE);

  if(fDoFlow)
  {

    for(Int_t i(0); i < fNumEtaGap; i++)
    {
      fListCumRef[i] = new TList();
      fListCumRef[i]->SetOwner(kTRUE);
      fListCumRef[i]->SetName(Form("fListRef_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumRef[i]);

      fListCumTracks[i] = new TList();
      fListCumTracks[i]->SetOwner(kTRUE);
      fListCumTracks[i]->SetName(Form("fListTracks_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumTracks[i]);

      fListCumPions[i] = new TList();
      fListCumPions[i]->SetOwner(kTRUE);
      fListCumPions[i]->SetName(Form("fListPions_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumPions[i]);

      fListCumKaons[i] = new TList();
      fListCumKaons[i]->SetOwner(kTRUE);
      fListCumKaons[i]->SetName(Form("fListKaons_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumKaons[i]);

      fListCumProtons[i] = new TList();
      fListCumProtons[i]->SetOwner(kTRUE);
      fListCumProtons[i]->SetName(Form("fListProtons_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumProtons[i]);

      fListCumK0s[i] = new TList();
      fListCumK0s[i]->SetOwner(kTRUE);
      fListCumK0s[i]->SetName(Form("fListK0s_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumK0s[i]);

      fListCumLambda[i] = new TList();
      fListCumLambda[i]->SetOwner(kTRUE);
      fListCumLambda[i]->SetName(Form("fListLambda_Gap%02.2g",10*fEtaGap[i]));
      fOutListCumulants->Add(fListCumLambda[i]);
    }

    // Katarina GF

    for(Int_t m(0); m < fNumEtaGap; m++)
    {
      for(Int_t i(0); i < fNumHarmonics; i++)
      {
        for(Int_t j(0); j < fNumSampleBins; j++)
        {
          if(!fSampling && j > 0)
            continue;

          fcn2Tracks[m][i][j] = new TProfile(Form("fTracksRef_n%d2_gap%02.2g_number%d",fHarmonics[i], fEtaGap[m]*10, j), Form("Tracks: <<2>> ref harmonics %d Gap %g sample %d; cent",fHarmonics[i],fEtaGap[m],j), fNumCentBins,fCentBinEdges);
          fcn2Tracks[m][i][j]->Sumw2();
          fListCumRef[m]->Add(fcn2Tracks[m][i][j]);
          // four particle cumulants does not have gap
          if(fEtaGap[m] >= -0.) // case with gap
            continue;

          fcn4Tracks[m][i][j] = new TProfile(Form("fTracksRef_n%d4_gap%02.2g_number%d",fHarmonics[i], fEtaGap[m]*10,j), Form("Tracks: <<4>> ref harmonics %d Gap %g sample %d; cent",fHarmonics[i],fEtaGap[m],j), fNumCentBins,fCentBinEdges);
          fcn4Tracks[m][i][j]->Sumw2();
          fListCumRef[m]->Add(fcn4Tracks[m][i][j]);
        }
      }
    }

    TString sVec[fGFKNumVectors] = {"Pos","Neg"};
    for(Int_t m(0); m < fNumEtaGap; m++)
    {
      for(Int_t i(0); i < fNumCentBins; i++)
      {
        for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
        {
          if(fEtaGap[m] < 0. && iVec == 1) //only positive vectors for No Gap situation
            continue;

          fInvMassPtK0s[iVec][m][i] = new TH2D(Form("fInvMassPtK0s_%s_Gap%02.2g_Cent%d",sVec[iVec].Data(),10*fEtaGap[m],i),Form("K^{0}_{S} candidates |#Delta#it{#eta}| > %g POI %s Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{m}^{V0}_{inv} (GeV/#it{c}^{2});",fEtaGap[m],sVec[iVec].Data(),fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins, fPtBinEdges, 200, fV0MinMassK0s, fV0MaxMassK0s);
          fInvMassPtK0s[iVec][m][i]->Sumw2();
          fListCumK0s[m]->Add(fInvMassPtK0s[iVec][m][i]);

          fInvMassPtLambda[iVec][m][i] = new TH2D(Form("fInvMassPtLambda_%s_Gap%02.2g_Cent%d",sVec[iVec].Data(),10*fEtaGap[m],i),Form("#Lambda candidates |#Delta#it{#eta}| > %g POI %s Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{m}^{V0}_{inv} (GeV/#it{c}^{2});",fEtaGap[m],sVec[iVec].Data(),fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins, fPtBinEdges, 200, fV0MinMassLambda, fV0MaxMassLambda);
          fInvMassPtLambda[iVec][m][i]->Sumw2();
          fListCumLambda[m]->Add(fInvMassPtLambda[iVec][m][i]);
        }

        for(Int_t k(0); k < fNumHarmonics; k++)
        {
          for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
          {
            if(fEtaGap[m] < 0. && iVec == 1) //only positive vectors for No Gap situation
              continue;

            fdn2K0s[iVec][m][k][i] = new TProfile2D(Form("fK0s_n%d2_%s_gap%02.2g_cent%d", fHarmonics[k],sVec[iVec].Data(),fEtaGap[m]*10, i), Form("K_{S}^{0}: <<2'>> harmonics %d diff Gap %g POI %s cent %d; pT; Minv;",fHarmonics[k],fEtaGap[m],sVec[iVec].Data(),i), fNumPtBins,fPtBinEdges,fNumMinvFlowBinsK0s,fMinvFlowBinEdgesK0s);
            fdn2K0s[iVec][m][k][i]->Sumw2();
            fListCumK0s[m]->Add(fdn2K0s[iVec][m][k][i]);

            fdn2Lambda[iVec][m][k][i] = new TProfile2D(Form("fLambda_n%d2_%s_gap%02.2g_cent%d", fHarmonics[k],sVec[iVec].Data(),fEtaGap[m]*10, i), Form("#Lambda: <<2'>> harmonics %d diff Gap %g POI %s cent %d; pT; Minv;",fHarmonics[k],fEtaGap[m],sVec[iVec].Data(),i), fNumPtBins,fPtBinEdges,fNumMinvFlowBinsLambda,fMinvFlowBinEdgesLambda);
            fdn2Lambda[iVec][m][k][i]->Sumw2();
            fListCumLambda[m]->Add(fdn2Lambda[iVec][m][k][i]);
          }

          if(fEtaGap[m] < 0.) // case with no gap (otherw. no need for 4 cummulants)
          {
            fdn4K0s[m][k][i] = new TProfile2D(Form("fK0s_n%d4_gap%02.2g_cent%d", fHarmonics[k],fEtaGap[m]*10, i), Form("K_{S}^{0}: <<4'>> harmonics %d diff Gap %g cent %d; pT; Minv;",fHarmonics[k],fEtaGap[m],i), fNumPtBins,fPtBinEdges,fNumMinvFlowBinsK0s,fMinvFlowBinEdgesK0s);
            fdn4K0s[m][k][i]->Sumw2();
            fListCumK0s[m]->Add(fdn4K0s[m][k][i]);

            fdn4Lambda[m][k][i] = new TProfile2D(Form("fLambda_n%d4_gap%02.2g_cent%d", fHarmonics[k],fEtaGap[m]*10, i), Form("#Lambda: <<4'>> harmonics %d diff Gap %g cent %d; pT; Minv;",fHarmonics[k],fEtaGap[m],i), fNumPtBins,fPtBinEdges,fNumMinvFlowBinsLambda,fMinvFlowBinEdgesLambda);
            fdn4Lambda[m][k][i]->Sumw2();
            fListCumLambda[m]->Add(fdn4Lambda[m][k][i]);
          }

          for(Int_t j(0); j < fNumSampleBins; j++)
          {
            if(!fSampling && j > 0)
              continue;

            for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
            {
              if(fEtaGap[m] < 0. && iVec == 1) //only positive vectors for No Gap situation
                continue;

              fdn2Tracks[iVec][m][k][i][j] = new TProfile(Form("fTracks_n%d2_%s_gap%02.2g_cent%d_number%d", fHarmonics[k],sVec[iVec].Data(),fEtaGap[m]*10, i, j), Form("#Tracks: <<2'>> harmonics %d diff Gap %g POI %s cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],sVec[iVec].Data(),i,j), fNumPtBins,fPtBinEdges);
              fdn2Tracks[iVec][m][k][i][j]->Sumw2();
              fListCumTracks[m]->Add(fdn2Tracks[iVec][m][k][i][j]);

              fdn2Pion[iVec][m][k][i][j] = new TProfile(Form("fPion_n%d2_%s_gap%02.2g_cent%d_number%d", fHarmonics[k],sVec[iVec].Data(),fEtaGap[m]*10, i, j), Form("#pi: <<2'>> harmonics %d diff Gap %g POI %s cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],sVec[iVec].Data(),i,j), fNumPtBins,fPtBinEdges);
              fdn2Pion[iVec][m][k][i][j]->Sumw2();
              fListCumPions[m]->Add(fdn2Pion[iVec][m][k][i][j]);

              fdn2Kaon[iVec][m][k][i][j] = new TProfile(Form("fKaon_n%d2_%s_gap%02.2g_cent%d_number%d", fHarmonics[k],sVec[iVec].Data(),fEtaGap[m]*10,i, j), Form("K: <<2'>> harmonics %d diff Gap %g POI %s cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],sVec[iVec].Data(),i,j), fNumPtBins,fPtBinEdges);
              fdn2Kaon[iVec][m][k][i][j]->Sumw2();
              fListCumKaons[m]->Add(fdn2Kaon[iVec][m][k][i][j]);

              fdn2Proton[iVec][m][k][i][j] = new TProfile(Form("fProton_n%d2_%s_gap%02.2g_cent%d_number%d", fHarmonics[k],sVec[iVec].Data(),fEtaGap[m]*10,i, j), Form("p: <<2'>> harmonics %d diff Gap %g POI %s cent %d sample %d ; pT",fHarmonics[k],fEtaGap[m],sVec[iVec].Data(),i,j), fNumPtBins,fPtBinEdges);
              fdn2Proton[iVec][m][k][i][j]->Sumw2();
              fListCumProtons[m]->Add(fdn2Proton[iVec][m][k][i][j]);

            }

            // four particle cumulants does not have gap
            if(fEtaGap[m] >= -0.) // case with gap
              continue;

            fdn4Tracks[m][k][i][j] = new TProfile(Form("fTracks_n%d4_gap%02.2g_cent%d_number%d", fHarmonics[k],fEtaGap[m]*10,i, j), Form("Tracks: <<4'>> harmonics %d Gap %g cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],i,j), fNumPtBins,fPtBinEdges);
            fdn4Tracks[m][k][i][j]->Sumw2();
            fListCumTracks[m]->Add(fdn4Tracks[m][k][i][j]);

            fdn4Pion[m][k][i][j] = new TProfile(Form("fPion_n%d4_gap%02.2g_cent%d_number%d", fHarmonics[k],fEtaGap[m]*10,i, j), Form("#pi: <<4'>> harmonics %d Gap %g cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],i,j), fNumPtBins,fPtBinEdges);
            fdn4Pion[m][k][i][j]->Sumw2();
            fListCumPions[m]->Add(fdn4Pion[m][k][i][j]);

            fdn4Kaon[m][k][i][j] = new TProfile(Form("fKaon_n%d4_gap%02.2g_cent%d_number%d", fHarmonics[k],fEtaGap[m]*10,i, j), Form("K: <<4'>> harmonics %d diff Gap %g cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],i,j), fNumPtBins,fPtBinEdges);
            fdn4Kaon[m][k][i][j]->Sumw2();
            fListCumKaons[m]->Add(fdn4Kaon[m][k][i][j]);

            fdn4Proton[m][k][i][j] = new TProfile(Form("fProton_n%d4_gap%02.2g_cent%d_number%d", fHarmonics[k],fEtaGap[m]*10,i, j), Form("p: <<4'>> harmonics %d diff Gap %0g cent %d sample %d; pT",fHarmonics[k],fEtaGap[m],i,j), fNumPtBins,fPtBinEdges);
            fdn4Proton[m][k][i][j]->Sumw2();
            fListCumProtons[m]->Add(fdn4Proton[m][k][i][j]);


          }
        }
      }
    }

    // end of Katarina GF

    if(fOldFlow)
    {
      // reference
      for(Int_t i(0); i < fNumHarmonics; i++)
      {
        for(Int_t j(0); j < fNumEtaGap; j++)
        {
          for(Int_t k(0); k < fNumSampleBins; k++)
          {
            if(!fSampling && k > 0)
              continue;

            fTracksRefTwo[i][j][k] = new TProfile(Form("fTracksRefTwo_n%d_Gap%02.2g_sample%d",fHarmonics[i],fEtaGap[j]*10,k),Form("#LT#LT2#GT#GT_{%d,|#Delta#it{#eta}| > %g} sample %d (ref. flow); centrality;",fHarmonics[i],fEtaGap[j],k),fNumCentBins,fCentBinEdges);
            fTracksRefTwo[i][j][k]->Sumw2();
            fOutListTracks->Add(fTracksRefTwo[i][j][k]);
          }
        }
      }

      // tracks diff
      if(fDiffFlow) // do differential flow switch
      {
        for(Int_t j(0); j < fNumHarmonics; j++)
        {
          for(Int_t k(0); k < fNumEtaGap; k++)
          {
            for(Int_t i = 0; i < fNumCentBins; i++)
            {
              for(Int_t m(0); m < fNumSampleBins; m++)
              {
                if(!fSampling && m > 0)
                  continue;

                fTracksDiffTwoPos[i][j][k][m] = new TProfile(Form("fTracksDiffTwoPos_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("#LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                fTracksDiffTwoPos[i][j][k][m]->Sumw2();
                fOutListTracks->Add(fTracksDiffTwoPos[i][j][k][m]);

                if(fPID)
                {
                  fPionsDiffTwoPos[i][j][k][m] = new TProfile(Form("fPionsDiffTwoPos_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("#pi: #LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                  fPionsDiffTwoPos[i][j][k][m]->Sumw2();
                  fOutListPID->Add(fPionsDiffTwoPos[i][j][k][m]);

                  fKaonsDiffTwoPos[i][j][k][m] = new TProfile(Form("fKaonsDiffTwoPos_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("K: #LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                  fKaonsDiffTwoPos[i][j][k][m]->Sumw2();
                  fOutListPID->Add(fKaonsDiffTwoPos[i][j][k][m]);

                  fProtonsDiffTwoPos[i][j][k][m] = new TProfile(Form("fProtonsDiffTwoPos_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("p: #LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                  fProtonsDiffTwoPos[i][j][k][m]->Sumw2();
                  fOutListPID->Add(fProtonsDiffTwoPos[i][j][k][m]);
                }

                if(!(fEtaGap[k] < 0.)) // make only with eta gap
                {
                  fTracksDiffTwoNeg[i][j][k][m] = new TProfile(Form("fTracksDiffTwoNeg_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("#LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} < -%g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                  fTracksDiffTwoNeg[i][j][k][m]->Sumw2();
                  fOutListTracks->Add(fTracksDiffTwoNeg[i][j][k][m]);

                  if(fPID)
                  {
                    fPionsDiffTwoNeg[i][j][k][m] = new TProfile(Form("fPionsDiffTwoNeg_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("#pi: #LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} < -%g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                    fPionsDiffTwoNeg[i][j][k][m]->Sumw2();
                    fOutListPID->Add(fPionsDiffTwoNeg[i][j][k][m]);

                    fKaonsDiffTwoNeg[i][j][k][m] = new TProfile(Form("fKaonsDiffTwoNeg_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("K: #LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} < -%g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                    fKaonsDiffTwoNeg[i][j][k][m]->Sumw2();
                    fOutListPID->Add(fKaonsDiffTwoNeg[i][j][k][m]);

                    fProtonsDiffTwoNeg[i][j][k][m] = new TProfile(Form("fProtonsDiffTwoNeg_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m),Form("p: #LT#LT2'#GT#GT_{%d, #Delta#it{#eta}| > %g, #it{#eta}^{POI} < -%g} Cent %g-%g%% sample %d (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges);
                    fProtonsDiffTwoNeg[i][j][k][m]->Sumw2();
                    fOutListPID->Add(fProtonsDiffTwoNeg[i][j][k][m]);
                  }

                }

              }
            }
          }
        }
      }
    }

  } // end of if(fDoFlow)

  Short_t iMaxMult = 1500;
  Short_t iNumBinsMult = 1500;
  if(fPbPb)
  {
    iMaxMult = 5000;
    iNumBinsMult = 500;
  }

  // Events
  fEventMult = new TH1D("fEventMult","Track multiplicity (all tracks in selected events); tracks;",100,0,10000);
  fOutListEvents->Add(fEventMult);
  fCentralityDis = new TH1D("fCentralityDis", "Event centrality distribution (final binning); centrality percentile; Events;", fNumCentBins,fCentBinEdges);
  fOutListEvents->Add(fCentralityDis);
  fCentDistUnitBin = new TH1D("fCentDistUnitBin", "Event centrality distibution (unit binning); centrality percentile; Events", fCentBinEdges[fNumCentBins], 0, fCentBinEdges[fNumCentBins]);
  fOutListEvents->Add(fCentDistUnitBin);


  // PID tracks histos
  if(fPID)
  {
    // counters
    Int_t iNPIDCounterBins = 6;
    TString sPIDCounterLabel[] = {"Input","PID TPC OK","PID TOF OK","TPC && TOF PID OK","#it{n#sigma}","Selected"};

    fPionsCounter = new TH1D("fPionsCounter", "#pi: Counter", iNPIDCounterBins,0,iNPIDCounterBins);
    fKaonsCounter = new TH1D("fKaonsCounter", "K: Counter", iNPIDCounterBins,0,iNPIDCounterBins);
    fProtonsCounter = new TH1D("fProtonsCounter", "p: Counter", iNPIDCounterBins,0,iNPIDCounterBins);
    for(Int_t i = 0; i < iNPIDCounterBins; i++)
    {
      fPionsCounter->GetXaxis()->SetBinLabel(i+1, sPIDCounterLabel[i].Data() );
      fKaonsCounter->GetXaxis()->SetBinLabel(i+1, sPIDCounterLabel[i].Data() );
      fProtonsCounter->GetXaxis()->SetBinLabel(i+1, sPIDCounterLabel[i].Data() );
    }
    fOutListPID->Add(fPionsCounter);
    fOutListPID->Add(fKaonsCounter);
    fOutListPID->Add(fProtonsCounter);

    fPionsMult = new TH1D("fPionsMult","#pi: Event multiplicity (selected); multiplicity",iNumBinsMult,0,iMaxMult);
    fOutListPID->Add(fPionsMult);
    fPionsPt = new TH1D("fPionsPt", "#pi: #it{p}_{T} (selected); #it{p}_{T} (GeV/#it{c});", 100, 0, 10);
    fOutListPID->Add(fPionsPt);
    fPionsEta = new TH1D("fPionsEta", "#pi: #it{#eta} (selected); #it{#eta};", 400, -2, 2);
    fOutListPID->Add(fPionsEta);
    fPionsPhi = new TH1D("fPionsPhi", "#pi: #it{#varphi} (selected); #it{#varphi};", 360, 0., TMath::TwoPi());
    fOutListPID->Add(fPionsPhi);
    fPionsTOFbeta = new TH2D("fPionsTOFbeta","#pi: TOF #beta (selected); #it{p} (GeV/#it{c}); TOF #beta", 1000,0,10, 100,0,1.5);
    fOutListPID->Add(fPionsTOFbeta);
    fPionsTPCdEdx = new TH2D("fPionsTPCdEdx","#pi: TPC #it{dEdx} (selected); #it{p} (GeV/#it{c}); TPC #it{dEdx}", 1000,0,10, 1000,0,1000);
    fOutListPID->Add(fPionsTPCdEdx);
    fPionsNsigmasTPCTOF = new TH2D("fPionsNsigmasTPCTOF","#pi: #it{n#sigma} TPC & TOF (selected); #it{n#sigma}^{TPC}; #it{n#sigma}^{TOF}", 20,-10,10, 20,-10,10);
    fOutListPID->Add(fPionsNsigmasTPCTOF);
    fPionsBayesAsPion = new TH1D("fPionsBayesAsPion","#pi: Bayes PID (as #pi); probability",100,0,1);
    fOutListPID->Add(fPionsBayesAsPion);
    fPionsBayesAsKaon = new TH1D("fPionsBayesAsKaon","#pi: Bayes PID (as K); probability",100,0,1);
    fOutListPID->Add(fPionsBayesAsKaon);
    fPionsBayesAsProton = new TH1D("fPionsBayesAsProton","#pi: Bayes PID (as p); probability",100,0,1);
    fOutListPID->Add(fPionsBayesAsProton);
    fPionsBayes = new TH3D("fPionsBayes","#pi: Bayes PID; probability #pi; probability K; probability p",100,0,1,100,0,1,100,0,1);
    fOutListPID->Add(fPionsBayes);

    fKaonsMult = new TH1D("fKaonsMult","K: Event multiplicity (selected); multiplicity",iNumBinsMult,0,iMaxMult);
    fOutListPID->Add(fKaonsMult);
    fKaonsPt = new TH1D("fKaonsPt", "K: #it{p}_{T} (selected); #it{p}_{T} (GeV/#it{c});", 100, 0, 10);
    fOutListPID->Add(fKaonsPt);
    fKaonsEta = new TH1D("fKaonsEta", "K: #it{#eta} (selected); #it{#eta};", 400, -2, 2);
    fOutListPID->Add(fKaonsEta);
    fKaonsPhi = new TH1D("fKaonsPhi", "K: #it{#varphi} (selected); #it{#varphi};", 360, 0., TMath::TwoPi());
    fOutListPID->Add(fKaonsPhi);
    fKaonsTOFbeta = new TH2D("fKaonsTOFbeta","K: TOF #beta (selected); #it{p} (GeV/#it{c}); TOF #beta", 1000,0,10, 100,0,1.5);
    fOutListPID->Add(fKaonsTOFbeta);
    fKaonsTPCdEdx = new TH2D("fKaonsTPCdEdx","K: TPC #it{dEdx} (selected); #it{p} (GeV/#it{c}); TPC #it{dEdx}", 1000,0,10, 1000,0,1000);
    fOutListPID->Add(fKaonsTPCdEdx);
    fKaonsNsigmasTPCTOF = new TH2D("fKaonsNsigmasTPCTOF","K: #it{n#sigma} TPC & TOF (selected); #it{n#sigma}^{TPC}; #it{n#sigma}^{TOF}", 20,-10,10, 20,-10,10);
    fOutListPID->Add(fKaonsNsigmasTPCTOF);
    fKaonsBayesAsPion = new TH1D("fKaonsBayesAsPion","K: Bayes PID (as #pi); probability",100,0,1);
    fOutListPID->Add(fKaonsBayesAsPion);
    fKaonsBayesAsKaon = new TH1D("fKaonsBayesAsKaon","K: Bayes PID (as K); probability",100,0,1);
    fOutListPID->Add(fKaonsBayesAsKaon);
    fKaonsBayesAsProton = new TH1D("fKaonsBayesAsProton","K: Bayes PID (as p); probability",100,0,1);
    fOutListPID->Add(fKaonsBayesAsProton);
    fKaonsBayes = new TH3D("fKaonsBayes","K: Bayes PID; probability #pi; probability K; probability p",100,0,1,100,0,1,100,0,1);
    fOutListPID->Add(fKaonsBayes);

    fProtonsMult = new TH1D("fProtonsMult","p: Event multiplicity (selected); multiplicity",iNumBinsMult,0,iMaxMult);
    fOutListPID->Add(fProtonsMult);
    fProtonsPt = new TH1D("fProtonsPt", "p: #it{p}_{T} (selected); #it{p}_{T} (GeV/#it{c});", 100, 0, 10);
    fOutListPID->Add(fProtonsPt);
    fProtonsEta = new TH1D("fProtonsEta", "p: #it{#eta} (selected); #it{#eta};", 400, -2, 2);
    fOutListPID->Add(fProtonsEta);
    fProtonsPhi = new TH1D("fProtonsPhi", "p: #it{#varphi} (selected); #it{#varphi};", 360, 0., TMath::TwoPi());
    fOutListPID->Add(fProtonsPhi);
    fProtonsTOFbeta = new TH2D("fProtonsTOFbeta","p: TOF #beta (selected); #it{p} (GeV/#it{c}); TOF #beta", 1000,0,10, 100,0,1.5);
    fOutListPID->Add(fProtonsTOFbeta);
    fProtonsTPCdEdx = new TH2D("fProtonsTPCdEdx","p: TPC #it{dEdx} (selected); #it{p} (GeV/#it{c}); TPC #it{dEdx}", 1000,0,10, 1000,0,1000);
    fOutListPID->Add(fProtonsTPCdEdx);
    fProtonsNsigmasTPCTOF = new TH2D("fProtonsNsigmasTPCTOF","p: #it{n#sigma} TPC & TOF (selected); #it{n#sigma}^{TPC}; #it{n#sigma}^{TOF}", 20,-10,10, 20,-10,10);
    fOutListPID->Add(fProtonsNsigmasTPCTOF);
    fProtonsBayesAsPion = new TH1D("fProtonsBayesAsPion","p: Bayes PID (as #pi); probability",100,0,1);
    fOutListPID->Add(fProtonsBayesAsPion);
    fProtonsBayesAsKaon = new TH1D("fProtonsBayesAsKaon","p: Bayes PID (as K); probability",100,0,1);
    fOutListPID->Add(fProtonsBayesAsKaon);
    fProtonsBayesAsProton = new TH1D("fProtonsBayesAsProton","p: Bayes PID (as p); probability",100,0,1);
    fOutListPID->Add(fProtonsBayesAsProton);
    fProtonsBayes = new TH3D("fProtonsBayes","p: Bayes PID; probability #pi; probability K; probability p",100,0,1,100,0,1,100,0,1);
    fOutListPID->Add(fProtonsBayes);


  }

  if(fOldFlow)
  {
    // V0s histos
    if(fPID)
    {
      for(Int_t k(0); k < fNumEtaGap; k++)
      {
        fV0sInvMassK0s[k] = new TH1D(Form("fV0sInvMass_K0s_Gap%02.2g",10*fEtaGap[k]),Form("K^{0}_{S} InvMass |#Delta#it{#eta}| > %g (selected); #it{m}_{inv} (GeV/#it{c}^{2});",fEtaGap[k]), 200,fV0MinMassK0s,fV0MaxMassK0s);
        fOutListV0s->Add(fV0sInvMassK0s[k]);

        fV0sInvMassLambda[k] = new TH1D(Form("fV0sInvMass_Lambda_Gap%02.2g",10*fEtaGap[k]),Form("#Lambda+#bar{#Lambda} InvMass |#Delta#it{#eta}| > %g (selected); #it{m}_{inv} (GeV/#it{c}^{2});",fEtaGap[k]), 200,fV0MinMassLambda,fV0MaxMassLambda);
        fOutListV0s->Add(fV0sInvMassLambda[k]);

      }

      for(Int_t j(0); j < fNumHarmonics; j++)
      {
        for(Int_t k(0); k < fNumEtaGap; k++)
        {
          for(Int_t i(0); i < fNumCentBins; i++)
          {
            fV0sPtInvMassK0s[i][j][k] = new TH2D(Form("fV0sPtInvMass_K0s_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i),Form("K^{0}_{S} candidates n=%d |#Delta#it{#eta}| > %g Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{m}^{V0}_{inv} (GeV/#it{c}^{2});",fHarmonics[j],fEtaGap[k],fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins, fPtBinEdges, 200, fV0MinMassK0s, fV0MaxMassK0s);
            fV0sPtInvMassK0s[i][j][k]->Sumw2();
            fOutListV0s->Add(fV0sPtInvMassK0s[i][j][k]);

            fV0sPtInvMassLambda[i][j][k] = new TH2D(Form("fV0sPtInvMass_Lambda_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i),Form("#Lambda+#bar{#Lambda} candidates n=%d |#Delta#it{#eta}| > %g Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{m}^{V0}_{inv} (GeV/#it{c}^{2});",fHarmonics[j],fEtaGap[k],fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins, fPtBinEdges, 200, fV0MinMassLambda, fV0MaxMassLambda);
            fV0sPtInvMassLambda[i][j][k]->Sumw2();
            fOutListV0s->Add(fV0sPtInvMassLambda[i][j][k]);

            for(Int_t m(0); m < fNumSampleBins; m++)
            {
              if(!fSampling && m > 0)
                continue;

              fV0sDiffTwoPos_K0s[i][j][k][m] = new TProfile2D(Form("fV0sDiffTwoPos_K0s_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m), Form("K^{0}_{S} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%% sample %d; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsK0s,fMinvFlowBinEdgesK0s);
              fV0sDiffTwoPos_K0s[i][j][k][m]->Sumw2();
              fOutListV0s->Add(fV0sDiffTwoPos_K0s[i][j][k][m]);

              if(!(fEtaGap[k] < 0.)) // make only with eta gap
              {
                fV0sDiffTwoNeg_K0s[i][j][k][m] = new TProfile2D(Form("fV0sDiffTwoNeg_K0s_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m), Form("K^{0}_{S} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, #it{#eta}^{POI} < -%g} Cent %g-%g%% sample %d; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsK0s,fMinvFlowBinEdgesK0s);
                fV0sDiffTwoNeg_K0s[i][j][k][m]->Sumw2();
                fOutListV0s->Add(fV0sDiffTwoNeg_K0s[i][j][k][m]);
              }


              fV0sDiffTwoPos_Lambda[i][j][k][m] = new TProfile2D(Form("fV0sDiffTwoPos_Lambda_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m), Form("#Lambda+#bar{#Lambda} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%% sample %d; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsLambda,fMinvFlowBinEdgesLambda);
              fV0sDiffTwoPos_Lambda[i][j][k][m]->Sumw2();
              fOutListV0s->Add(fV0sDiffTwoPos_Lambda[i][j][k][m]);

              if(!(fEtaGap[k] < 0.)) // make only with eta gap
              {
                fV0sDiffTwoNeg_Lambda[i][j][k][m] = new TProfile2D(Form("fV0sDiffTwoNeg_Lambda_n%d_Gap%02.2g_Cent%d_sample%d",fHarmonics[j],10*fEtaGap[k],i,m), Form("#Lambda+#bar{#Lambda} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, it{#eta}^{POI} < -%g} Cent %g-%g%% sample %d; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1],m),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsLambda,fMinvFlowBinEdgesLambda);
                fV0sDiffTwoNeg_Lambda[i][j][k][m]->Sumw2();
                fOutListV0s->Add(fV0sDiffTwoNeg_Lambda[i][j][k][m]);
              }
            }
          }
        }
      }
    }
  }

  Int_t iNEventCounterBins = 1;

	TString sQAlabel[fQANumSteps] = {"Before","After"/*,"Test"*/};
  // QA events output
  if(fPP)
  {
    if(fRejectOutOfBunchPU)
    {
      iNEventCounterBins = 14;
      TString sEventCounterLabel[] = {"Input","AOD OK","Physics selection OK","PV OK","SPD Vtx OK","Pileup SPD OK","Pileup MV OK","OOBPU OK","SPDclTrBg OK","SPDPU OK","Utils OK","PV #it{z} OK","At least 2 selected tracks","At least 1 V0 candidate"};
      fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
      for(Int_t i = 0; i < iNEventCounterBins; i++)
        fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    }
    else
    {
      iNEventCounterBins = 10;
      TString sEventCounterLabel[] = {"Input","AOD OK","Physics selection OK","PV OK","SPD Vtx OK","Pileup SPD OK","Pileup MV OK","PV #it{z} OK","At least 2 selected tracks","At least 1 V0 candidate"};
      fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
      for(Int_t i = 0; i < iNEventCounterBins; i++)
        fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    }
    fOutListEvents->Add(fEventCounter);
  }
  else
  {
    iNEventCounterBins = 9;
    TString sEventCounterLabel[] = {"Input","AOD OK","Pile-up OK","PV OK","SPD Vtx OK","PV #it{z} OK","Centrality OK","At least 2 selected tracks","At least 1 V0 candidate"};
    fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
    for(Int_t i = 0; i < iNEventCounterBins; i++)
      fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    fOutListEvents->Add(fEventCounter);
  }

  fSampleCounter = new TH2D("fSampleCounter","Event distribution if sampling bins; sampling bin index; centrality bin index; events",fNumSampleBins,0,fNumSampleBins, fNumCentBins, 0, fNumCentBins);
  for(Int_t i = 0; i < fNumSampleBins; i++)
    fSampleCounter->GetXaxis()->SetBinLabel(i+1,Form("%d",i));
  fOutListEvents->Add(fSampleCounter);

  fCentSPDvsV0M = new TH2D("fCentSPDvsV0M", "V0M-cent vs SPD-cent; V0M; SPD-cent", 100, 0, 100, 100, 0, 100);
  fOutListEvents->Add(fCentSPDvsV0M);

  Short_t iNEventTriggerBins = 4;
  TString sEventTriggerLabel[] = {"kINT7","kHighMultV0","kHighMultSPD","Other"};

  Short_t iNEventTriggersPresentBins = 4;
  TString sEventTriggersPresentLabel[] = {"kINT7","kHighMultV0","CINT7-B-NOPF-CENT","CVHMV0M-B-SPD2-CENT"};

  for(Int_t i(0); i < fQANumSteps; i++)
  {
    fQAEventsPVz[i] = new TH1D(Form("fQAEventsPVz_%s",sQAlabel[i].Data()),Form("QA Events: PV #it{z} (%s cuts); #it{z} (cm);",sQAlabel[i].Data()),101,-50,50);
    fOutListEvents->Add(fQAEventsPVz[i]);
    fQAEventsNumContrPV[i] = new TH1D(Form("fQAEventsNumContrPV_%s",sQAlabel[i].Data()),Form("QA Events: Number of contributors to AOD PV (%s cuts); Number of contributors;",sQAlabel[i].Data()),20,0,20);
    fOutListEvents->Add(fQAEventsNumContrPV[i]);
    fQAEventsNumSPDContrPV[i] = new TH1D(Form("fQAEventsNumSPDContrPV_%s",sQAlabel[i].Data()),Form("QA Events: Number of contributors to SPD PV (%s cuts); Number of contributors;",sQAlabel[i].Data()),20,0,20);
    fOutListEvents->Add(fQAEventsNumSPDContrPV[i]);
    fQAEventsDistPVSPD[i] = new TH1D(Form("fQAEventsDistPVSPD_%s",sQAlabel[i].Data()),Form("QA Events: #it{z}-istance between SPD vertex & PV (%s cuts); #it{z} (cm);",sQAlabel[i].Data()),50,0,5);
    fOutListEvents->Add(fQAEventsDistPVSPD[i]);

    if(fPP || fPPb)
    {
      fQAEventsSPDresol[i] = new TH1D(Form("fQAEventsSPDresol_%s",sQAlabel[i].Data()),Form("QA Events: SPD vertexer z resolution (%s cuts); #it{z} resolution (?cm);",sQAlabel[i].Data()),150,0,15);
      fOutListEvents->Add(fQAEventsSPDresol[i]);
      fQAEventsTriggerSelection[i] = new TH1D(Form("fQAEventsTriggerSelection_%s",sQAlabel[i].Data()),Form("QA Events: 2016 pp trigger selection (%s cuts); #it{z} (cm);",sQAlabel[i].Data()),iNEventTriggerBins,0,iNEventTriggerBins);
      for(Short_t j(0); j < iNEventTriggerBins; j++)
      {
        fQAEventsTriggerSelection[i]->GetXaxis()->SetBinLabel(j+1,sEventTriggerLabel[j].Data() );
      }
      fOutListEvents->Add(fQAEventsTriggerSelection[i]);


      fQAEventsTriggers[i] = new TH1D(Form("fQAEventsTriggers_%s",sQAlabel[i].Data()),Form("QA Events: Triggers present (%s Physics selection)",sQAlabel[i].Data()),iNEventTriggersPresentBins,0,iNEventTriggersPresentBins);
      for(Short_t j(0); j < iNEventTriggersPresentBins; j++)
      {
        fQAEventsTriggers[i]->GetXaxis()->SetBinLabel(j+1,sEventTriggersPresentLabel[j].Data() );
      }
      fOutListEvents->Add(fQAEventsTriggers[i]);
    }
  }

  // Tracks scan output
  if(fTracksScan)
  {
    const Short_t iNumBinsScanCounter = 5;
    TString sTracksScanCounterLabel[] = {"All","Selected","ITS PID OK","TPC PID OK","TOF PID OK"};

    fTracksScanCounter = new TH2D("fTracksScanCounter","Tracks: Scan counter; FB; ", fNumScanFB,0,fNumScanFB, iNumBinsScanCounter,0,iNumBinsScanCounter);
    fTracksScanPt = new TH2D("fTracksScanPt","Tracks: #it{p}_T scan; FB; #it{p}_T (GeV/#it{c})", fNumScanFB,0,fNumScanFB, 100,0,10);
    fTracksScanPhi = new TH2D("fTracksScanPhi","Tracks: #it{#varphi} scan; FB; #it{#varphi}", fNumScanFB,0,fNumScanFB, 100,0,TMath::TwoPi());
    fTracksScanEta = new TH2D("fTracksScanEta","Tracks: #it{#eta} scan; FB; #it{#eta}", fNumScanFB,0,fNumScanFB, 401,-2,2);

    for(Short_t i(0); i < fNumScanFB; i++)
    {
      fTracksScanCounter->GetXaxis()->SetBinLabel(i+1, Form("%d",fTracksScanFB[i]));
      fTracksScanPt->GetXaxis()->SetBinLabel(i+1, Form("%d",fTracksScanFB[i]));
      fTracksScanPhi->GetXaxis()->SetBinLabel(i+1, Form("%d",fTracksScanFB[i]));
      fTracksScanEta->GetXaxis()->SetBinLabel(i+1, Form("%d",fTracksScanFB[i]));
    }

    for(Short_t i(0); i < iNumBinsScanCounter; i++)
    {
      fTracksScanCounter->GetYaxis()->SetBinLabel(i+1, Form("%s",sTracksScanCounterLabel[i].Data()) );
    }

    fOutListTracks->Add(fTracksScanCounter);
    fOutListTracks->Add(fTracksScanPt);
    fOutListTracks->Add(fTracksScanPhi);
    fOutListTracks->Add(fTracksScanEta);
  }

  // QA tracks output

  Short_t iNTracksCounterBins = 8;
  TString sTracksCounterLabel[] = {"Input","Track OK","FB OK","TPC Ncls OK","Eta OK","Pt OK","DCAz OK","Selected"};
  fTracksCounter = new TH1D("fTrackCounter","Track Counter",iNTracksCounterBins,0,iNTracksCounterBins);
  for(Int_t i = 0; i < iNTracksCounterBins; i++)
    fTracksCounter->GetXaxis()->SetBinLabel(i+1, sTracksCounterLabel[i].Data() );
  fOutListTracks->Add(fTracksCounter);

  Short_t iNFilterMapBinBins = 12; // filter map bins (passing bit)

  for(Int_t i(0); i < fQANumSteps; i++)
  {
    fQATracksMult[i] = new TH1D(Form("fQATracksMult_%s",sQAlabel[i].Data()),Form("QA Tracks: Number of tracks in selected events (%s cuts); #it{N}^{tracks}",sQAlabel[i].Data()), iNumBinsMult,0,iMaxMult);
    fOutListTracks->Add(fQATracksMult[i]);
    fQATracksFilterMap[i] = new TH1D(Form("fQATracksFilterMap_%s",sQAlabel[i].Data()),Form("QA Tracks: Track filter bit map (%s cuts); filter bit",sQAlabel[i].Data()), 5000,0,5000);
    fOutListTracks->Add(fQATracksFilterMap[i]);

    fQATracksFilterMapBit[i] = new TH1D(Form("fQATracksFilterMapBit_%s",sQAlabel[i].Data()),Form("QA Tracks: Track filter bits (%s cuts); filter bit",sQAlabel[i].Data()), iNFilterMapBinBins,0,iNFilterMapBinBins);
    for(Int_t j = 0; j < iNFilterMapBinBins; j++)
      fQATracksFilterMapBit[i]->GetXaxis()->SetBinLabel(j+1, Form("%d",j));
    fOutListTracks->Add(fQATracksFilterMapBit[i]);

    fQATracksNumTPCcls[i] = new TH1D(Form("fQATracksNumTPCcls_%s",sQAlabel[i].Data()),Form("QA Tracks: Track number of TPC clusters (%s cuts); #it{N}^{TPC clusters}",sQAlabel[i].Data()), 200,0,200);
    fOutListTracks->Add(fQATracksNumTPCcls[i]);
    fQATracksPt[i] = new TH1D(Form("fQATracksPt_%s",sQAlabel[i].Data()),Form("QA Tracks: Track #it{p}_{T} (%s cuts); #it{p}_{T} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,10.);
    fOutListTracks->Add(fQATracksPt[i]);
    fQATracksEta[i] = new TH1D(Form("fQATracksEta_%s",sQAlabel[i].Data()),Form("QA Tracks: Track #it{#eta} (%s cuts); #it{#eta}",sQAlabel[i].Data()), 401,-2,2);
    fOutListTracks->Add(fQATracksEta[i]);
    fQATracksPhi[i] = new TH1D(Form("fQATracksPhi_%s",sQAlabel[i].Data()),Form("QA Tracks: Track #it{#varphi} (%s cuts); #it{#varphi}",sQAlabel[i].Data()), 100,0.,TMath::TwoPi());
    fOutListTracks->Add(fQATracksPhi[i]);
    fQATracksDCAxy[i] = new TH1D(Form("fQATracksDCAxy_%s",sQAlabel[i].Data()),Form("QA Tracks: Track DCA_#it{xy} (%s cuts); DCA_#it{xy} (cm)",sQAlabel[i].Data()), 100,-10.,10);
    fOutListTracks->Add(fQATracksDCAxy[i]);
    fQATracksDCAz[i] = new TH1D(Form("fQATracksDCAz_%s",sQAlabel[i].Data()),Form("QA Tracks: Track DCA_#it{z} (%s cuts); DCA_#it{z} (cm)",sQAlabel[i].Data()), 100,-5.,5);
    fOutListTracks->Add(fQATracksDCAz[i]);
    fQATracksTPCdEdx[i] = new TH2D(Form("fQATracksTPCdEdx_%s",sQAlabel[i].Data()),Form("QA Tracks: TPC PID information (%s cuts); #it{p} (GeV/#it{c}); TPC dEdx (au)",sQAlabel[i].Data()), 100,0,10, 2000,-10,1000);
    fOutListTracks->Add(fQATracksTPCdEdx[i]);
    fQATracksTOF[i] = new TH2D(Form("fQATracksTOF_%s",sQAlabel[i].Data()),Form("QA Tracks: TOF #it{t-t_{0}} information (%s cuts); #it{p} (GeV/#it{c}); TOF #it{t-t_{0}} (ps?)  ",sQAlabel[i].Data()), 100,0,10, 20000,-10,1000);
    fOutListTracks->Add(fQATracksTOF[i]);
    fQATracksTOFbeta[i] = new TH2D(Form("fQATracksTOFbeta_%s",sQAlabel[i].Data()),Form("QA Tracks: TOF #beta information (%s cuts); #it{p} (GeV/#it{c}); TOF #beta",sQAlabel[i].Data()), 1000,0,10, 500,0, 1.5);
    fOutListTracks->Add(fQATracksTOFbeta[i]);

    Int_t iNPIDstatusBins = 4;
    TString sPIDstatus[] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};

    fQATracksTPCstatus[i] = new TH1D(Form("fQATracksTPCstatus_%s",sQAlabel[i].Data()),Form("QA Tracks: PID status: TPC (%s cuts); isTOFok?",sQAlabel[i].Data()), iNPIDstatusBins,0,iNPIDstatusBins);
    fOutListTracks->Add(fQATracksTPCstatus[i]);
    fQATracksTOFstatus[i] = new TH1D(Form("fQATracksTOFstatus_%s",sQAlabel[i].Data()),Form("QA Tracks: PID status: TOF (%s cuts); isTOFok?",sQAlabel[i].Data()), iNPIDstatusBins,0,iNPIDstatusBins);
    fOutListTracks->Add(fQATracksTOFstatus[i]);

    for(Int_t j = 0; j < iNPIDstatusBins; j++)
    {
      fQATracksTOFstatus[i]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
      fQATracksTPCstatus[i]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
    }
  }

  // QA tracks PID output
  if(fPID)
  {
    fQAPIDTOFbetaNoTOF = new TH2D(Form("fQAPIDTOFbetaNoTOF"),Form("#pi,K,p: TOF #beta with NO TOF; #it{p} (GeV/#it{c}); TOF #beta"), 1000,0,10, 100,0,1.5);
    fOutListPID->Add(fQAPIDTOFbetaNoTOF);

    for(Int_t i(0); i < fQANumSteps; i++)
    {
      fQAPIDTOFbeta[i] = new TH2D(Form("fQAPIDTOFbeta_%s",sQAlabel[i].Data()),Form("#pi,K,p: TOF #beta (%s cuts); #it{p} (GeV/#it{c}); TOF #beta",sQAlabel[i].Data()), 1000,0,10, 100,0,1.5);
      fOutListPID->Add(fQAPIDTOFbeta[i]);
      fQAPIDTOFbetaWithTOF[i] = new TH2D(Form("fQAPIDTOFbetaWithTOF_%s",sQAlabel[i].Data()),Form("#pi,K,p: TOF #beta rejecting tracks with no TOF info (%s cuts); #it{p} (GeV/#it{c}); TOF #beta",sQAlabel[i].Data()), 1000,0,10, 100,0,1.5);
      fOutListPID->Add(fQAPIDTOFbetaWithTOF[i]);
      fQAPIDTPCdEdx[i] = new TH2D(Form("fQAPIDTPCdEdx_%s",sQAlabel[i].Data()),Form("#pi,K,p: TPC #it{dEdx} (%s cuts); #it{p} (GeV/#it{c}); TPC #it{dEdx}",sQAlabel[i].Data()), 1000,0,10, 1000,0,1000);
      fOutListPID->Add(fQAPIDTPCdEdx[i]);
      fQAPIDNsigmasTPCasPion[i] = new TH2D(Form("fQAPIDNsigmasTPCasPion_%s",sQAlabel[i].Data()),Form("#pi,K,p: #it{n#sigma} TPC as #pi (%s cuts); #it{p} (GeV/#it{c}); #it{n#sigma}^{TPC}",sQAlabel[i].Data()), 100,0,10, 100,-10,10);
      fOutListPID->Add(fQAPIDNsigmasTPCasPion[i]);
      fQAPIDNsigmasTOFasPion[i] = new TH2D(Form("fQAPIDNsigmasTOFasPion_%s",sQAlabel[i].Data()),Form("#pi,K,p: #it{n#sigma} TOF as #pi (%s cuts); #it{p} (GeV/#it{c}); #it{n#sigma}^{TOF}",sQAlabel[i].Data()), 100,0,10, 100,-10,10);
      fOutListPID->Add(fQAPIDNsigmasTOFasPion[i]);
      fQAPIDNsigmasTPCasKaon[i] = new TH2D(Form("fQAPIDNsigmasTPCasKaon_%s",sQAlabel[i].Data()),Form("#pi,K,p: #it{n#sigma} TPC as K (%s cuts); #it{p} (GeV/#it{c}); #it{n#sigma}^{TPC}",sQAlabel[i].Data()), 100,0,10, 100,-10,10);
      fOutListPID->Add(fQAPIDNsigmasTPCasKaon[i]);
      fQAPIDNsigmasTOFasKaon[i] = new TH2D(Form("fQAPIDNsigmasTOFasKaon_%s",sQAlabel[i].Data()),Form("#pi,K,p: #it{n#sigma} TOF as K (%s cuts); #it{p} (GeV/#it{c}); #it{n#sigma}^{TOF}",sQAlabel[i].Data()), 100,0,10, 100,-10,10);
      fOutListPID->Add(fQAPIDNsigmasTOFasKaon[i]);
      fQAPIDNsigmasTPCasProton[i] = new TH2D(Form("fQAPIDNsigmasTPCasProton_%s",sQAlabel[i].Data()),Form("#pi,K,p: #it{n#sigma} TPC as p (%s cuts); #it{p} (GeV/#it{c}); #it{n#sigma}^{TPC}",sQAlabel[i].Data()), 100,0,10, 100,-10,10);
      fOutListPID->Add(fQAPIDNsigmasTPCasProton[i]);
      fQAPIDNsigmasTOFasProton[i] = new TH2D(Form("fQAPIDNsigmasTOFasProton_%s",sQAlabel[i].Data()),Form("#pi,K,p: #it{n#sigma} TOF as p (%s cuts); #it{p} (GeV/#it{c}); #it{n#sigma}^{TOF}",sQAlabel[i].Data()), 100,0,10, 100,-10,10);
      fOutListPID->Add(fQAPIDNsigmasTOFasProton[i]);
      fQAPIDBayesProbPion[i] = new TH1D(Form("fQAPIDBayesProbPion_%s",sQAlabel[i].Data()),Form("#pi,K,p: Bayesian probablibity as #pi (%s cuts); probability",sQAlabel[i].Data()), 200, 0., 1.);
      fOutListPID->Add(fQAPIDBayesProbPion[i]);
      fQAPIDBayesProbKaon[i] = new TH1D(Form("fQAPIDBayesProbKaon_%s",sQAlabel[i].Data()), Form("#pi,K,p: Bayesian probablibity as K (%s cuts); probability",sQAlabel[i].Data()), 200, 0., 1.);
      fOutListPID->Add(fQAPIDBayesProbKaon[i]);
      fQAPIDBayesProbProton[i] = new TH1D(Form("fQAPIDBayesProbProton_%s",sQAlabel[i].Data()),Form("#pi,K,p: Bayesian probablibity as proton (%s cuts); probability",sQAlabel[i].Data()), 200, 0., 1.);
      fOutListPID->Add(fQAPIDBayesProbProton[i]);
      fQAPIDBayesProb[i] = new TH3D(Form("fQAPIDBayesProb_%s",sQAlabel[i].Data()), Form("#pi,K,p: Bayesian probability (%s cuts); probability #pi; probability K; probability p",sQAlabel[i].Data()), 100,0,1, 100,0,1, 100,0,1 );
      fOutListPID->Add(fQAPIDBayesProb[i]);
    }
  }

  // QA V0s output
  if(fDoV0s)
  {
    Int_t iNV0sCounterBins = 18;
    TString sV0sCounterLabel[] = {"Input","Daughters exists","Reco. method","Mother acceptance (#it{#eta},#it{p}_{T})","Decay radius","TPC refit","not Kink","Daughters charge","Daughter acceptance (#it{#eta},#it{p}_{T})","DCA to PV","Daughters DCA","Selected","K^{0}_{S}","#Lambda","#bar{#Lambda}","#Lambda && #bar{#Lambda}","K0s && (#Lambda || #bar{#Lambda})","K0s && #Lambda && #bar{#Lambda}"};
    fQAV0sCounter = new TH1D("fQAV0sCounter","QA V^{0}_{S} counter",iNV0sCounterBins,0,iNV0sCounterBins);
    for(Int_t i = 0; i < iNV0sCounterBins; i++)
      fQAV0sCounter->GetXaxis()->SetBinLabel(i+1, sV0sCounterLabel[i].Data() );
    fOutListV0s->Add(fQAV0sCounter);

    Int_t iNV0sCounterK0sBins = 7;
    TString sV0sCounterK0sLabel[] = {"Input","Inv. Mass","Mother #it{y} acceptance","CPA","Life-time","Armenteros","Selected"};
    fQAV0sCounterK0s = new TH1D("fQAV0sCounterK0s","QA K^{0}_{S} counter",iNV0sCounterK0sBins,0,iNV0sCounterK0sBins);
    for(Int_t i = 0; i < iNV0sCounterK0sBins; i++)
      fQAV0sCounterK0s->GetXaxis()->SetBinLabel(i+1, sV0sCounterK0sLabel[i].Data() );
    fOutListV0s->Add(fQAV0sCounterK0s);

    Int_t iNV0sCounterLambdaBins = 7;
    TString sV0sCounterLambdaLabel[] = {"Input","Inv. Mass","Mother #it{y} acceptance","CPA","Life-time","proton PID","Selected"};
    fQAV0sCounterLambda = new TH1D("fQAV0sCounterLambda","QA #Lambda/#bar{#Lambda} counter",iNV0sCounterLambdaBins,0,iNV0sCounterLambdaBins);
    for(Int_t i = 0; i < iNV0sCounterLambdaBins; i++)
      fQAV0sCounterLambda->GetXaxis()->SetBinLabel(i+1, sV0sCounterLambdaLabel[i].Data() );
    fOutListV0s->Add(fQAV0sCounterLambda);

  	for(Int_t i(0); i < fQANumSteps; i++)
    {
    	fQAV0sRecoMethod[i] = new TH1D(Form("fQAV0sRecoMethod_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Reconstruction method (%s cuts)",sQAlabel[i].Data()), 2,-0.5,1.5);
    	fQAV0sRecoMethod[i]->GetXaxis()->SetBinLabel(1, "offline");
    	fQAV0sRecoMethod[i]->GetXaxis()->SetBinLabel(2, "online (on-the-fly)");
    	fOutListV0s->Add(fQAV0sRecoMethod[i]);
  		fQAV0sTPCRefit[i] = new TH1D(Form("fQAV0sTPCRefit_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: TPC refit (%s cuts)",sQAlabel[i].Data()), 2,-0.5,1.5);
  		fQAV0sTPCRefit[i]->GetXaxis()->SetBinLabel(1, "not refit");
    	fQAV0sTPCRefit[i]->GetXaxis()->SetBinLabel(2, "refit");
    	fOutListV0s->Add(fQAV0sTPCRefit[i]);
  		fQAV0sKinks[i] = new TH1D(Form("fQAV0sKinks_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Kinks (%s cuts)",sQAlabel[i].Data()), 2,-0.5,1.5);
  		fQAV0sKinks[i]->GetXaxis()->SetBinLabel(1, "NOT AliAODVertex::kKink");
    	fQAV0sKinks[i]->GetXaxis()->SetBinLabel(2, "AliAODVertex:kKink");
    	fOutListV0s->Add(fQAV0sKinks[i]);
  		fQAV0sDCAtoPV[i] = new TH1D(Form("fQAV0sDCAtoPV_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter DCA to PV (%s cuts); daughter DCA^{PV} (cm)",sQAlabel[i].Data()), 100,0.,10.);
    	fOutListV0s->Add(fQAV0sDCAtoPV[i]);
  		fQAV0sDCADaughters[i] = new TH1D(Form("fQAV0sDCADaughters_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: DCA among daughters (%s cuts); DCA^{daughters} (cm)",sQAlabel[i].Data()), 100,0.,10.);
    	fOutListV0s->Add(fQAV0sDCADaughters[i]);
  		fQAV0sDecayRadius[i] = new TH1D(Form("fQAV0sDecayRadius_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Decay radius (%s cuts); #it{r_{xy}}^{decay} (cm)",sQAlabel[i].Data()), 200,0.,500.);
    	fOutListV0s->Add(fQAV0sDecayRadius[i]);
      fQAV0sDaughterPt[i] = new TH1D(Form("fQAV0sDaughterPt_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter #it{p}_{T} (%s cuts); #it{p}_{T}^{daughter} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,10.);
      fOutListV0s->Add(fQAV0sDaughterPt[i]);
      fQAV0sDaughterPhi[i] = new TH1D(Form("fQAV0sDaughterPhi_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter #it{#varphi} (%s cuts); #it{#varphi}^{daughter} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,TMath::TwoPi());
      fOutListV0s->Add(fQAV0sDaughterPhi[i]);
      fQAV0sDaughterTPCdEdxPt[i] = new TH2D(Form("fQAV0sDaughterTPCdEdxPt_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: TPC dEdx daughters (%s cuts); #it{p}_{T}^{daughter} (GeV/#it{c}); TPC dEdx (au);",sQAlabel[i].Data()), 100,0.,10,200,0.,200.);
    	fOutListV0s->Add(fQAV0sDaughterTPCdEdxPt[i]);
      fQAV0sDaughterEta[i] = new TH1D(Form("fQAV0sDaughterEta_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter #it{#eta} (%s cuts); #it{#eta}^{daugter}",sQAlabel[i].Data()), 401,-2,2);
      fOutListV0s->Add(fQAV0sDaughterEta[i]);
      fQAV0sMotherPt[i] = new TH1D(Form("fQAV0sMotherPt_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{p}_{T} (%s cuts); #it{p}_{T}^{V0} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,10.);
      fOutListV0s->Add(fQAV0sMotherPt[i]);
      fQAV0sMotherPhi[i] = new TH1D(Form("fQAV0sMotherPhi_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{#varphi} (%s cuts); #it{#varphi}^{V0} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,TMath::TwoPi());
      fOutListV0s->Add(fQAV0sMotherPhi[i]);
      fQAV0sMotherEta[i] = new TH1D(Form("fQAV0sMotherEta_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{#eta} (%s cuts); #it{#eta}^{V0}",sQAlabel[i].Data()), 401,-2,2);
      fOutListV0s->Add(fQAV0sMotherEta[i]);
      fQAV0sMotherRapK0s[i] = new TH1D(Form("fQAV0sMotherRapK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{y} (K^{0}_{S} hypo) (%s cuts); #it{y}^{V0,K0s}",sQAlabel[i].Data()), 400,-2,2);
      fOutListV0s->Add(fQAV0sMotherRapK0s[i]);
      fQAV0sMotherRapLambda[i] = new TH1D(Form("fQAV0sMotherRapLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{y} (Lambda/#bar{#Lambda} hypo) (%s cuts); #it{y}^{V0,#Lambda}",sQAlabel[i].Data()),400,-2,2);
      fOutListV0s->Add(fQAV0sMotherRapLambda[i]);
      fQAV0sInvMassK0s[i] = new TH1D(Form("fQAV0sInvMassK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: K^{0}_{S}: InvMass (%s cuts); #it{m}_{inv} (GeV/#it{c}^{2});",sQAlabel[i].Data()), 200,fV0MinMassK0s,fV0MaxMassK0s);
      fOutListV0s->Add(fQAV0sInvMassK0s[i]);
      fQAV0sInvMassLambda[i] = new TH1D(Form("fQAV0sInvMassLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: InvMass (%s cuts); #it{m}_{inv} (GeV/#it{c}^{2});",sQAlabel[i].Data()), 80,fV0MinMassLambda,fV0MaxMassLambda);
      fOutListV0s->Add(fQAV0sInvMassLambda[i]);
      fQAV0sCPAK0s[i] = new TH1D(Form("fQAV0sCPAK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: K^{0}_{S}: CPA (%s cuts); CPA^{K0s}",sQAlabel[i].Data()), 100,0.9,1.);
      fOutListV0s->Add(fQAV0sCPAK0s[i]);
      fQAV0sCPALambda[i] = new TH1D(Form("fQAV0sCPALambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: CPA (%s cuts); CPA^{#Lambda}",sQAlabel[i].Data()), 100, 0.9,1.);
      fOutListV0s->Add(fQAV0sCPALambda[i]);
      fQAV0sNumTauK0s[i] = new TH1D(Form("fQAV0sNumTauK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}:  K^{0}_{S}: Number of #it{c#tau} (%s cuts); #it{c#tau}^{K0s} (cm)",sQAlabel[i].Data()), 100, 0.,10.);
      fOutListV0s->Add(fQAV0sNumTauK0s[i]);
      fQAV0sNumTauLambda[i] = new TH1D(Form("fQAV0sNumTauLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Number of #it{c#tau} (%s cuts); #it{c#tau}^{#Lambda} (cm)",sQAlabel[i].Data()), 300, 0.,30);
      fOutListV0s->Add(fQAV0sNumTauLambda[i]);
      fQAV0sArmenterosK0s[i] = new TH2D(Form("fQAV0sArmenterosK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}:  K^{0}_{S}: Armenteros-Podolaski plot (%s cuts); #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});",sQAlabel[i].Data()), 100,-1.,1., 100,0.,0.3);
      fOutListV0s->Add(fQAV0sArmenterosK0s[i]);
      fQAV0sArmenterosLambda[i] = new TH2D(Form("fQAV0sArmenterosLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: Armenteros-Podolaski plot (%s cuts); #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});",sQAlabel[i].Data()), 100,-1.,1., 100,0.,0.3);
      fOutListV0s->Add(fQAV0sArmenterosLambda[i]);
      fQAV0sProtonNumSigmaPtLambda[i] = new TH2D(Form("fQAV0sProtonNumSigmaPtLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: (anti-)proton PID (%s cuts); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});",sQAlabel[i].Data()), 100,0.,10,100,0.,10.);
      fOutListV0s->Add(fQAV0sProtonNumSigmaPtLambda[i]);
  	}
  }

  PostData(1, fOutListCumulants);
  PostData(2, fOutListEvents);
  PostData(3, fOutListTracks);
  PostData(4, fOutListPID);
	PostData(5, fOutListV0s);
	//PostData(6, fOutListQA);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserExec(Option_t *)
{
  // this function is called once for each event
  if(!fPP && !fPPb && !fPbPb)
  {
    ::Error("UserExec","Neighter pp, pPb nor PbPb analysis switch is on. Terminating!");
    return;
  }

	fEventCounter->Fill(0); // input event

  if(fAODAnalysis) // AOD analysis
  {
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
  }

  if( !fAODAnalysis || dynamic_cast<AliESDEvent*>(InputEvent()) ) // ESD analysis
  {
    ::Warning("UserExec","ESD event: not implemented. Terminating!");
    return;
  }

  fEventCounter->Fill(1); // event AOD ok

  FillEventQA(fAOD,0); // before cuts event QA

  // event selection

  if(fPbPb)
  {
    if(fUseOldCent) // as by You Run1
    {
      if(!OldIsEventSelected(fAOD))
      {
        //printf("Rejected\n");
        return;
      }
    }
    else // Run2 cent selection based on AliMult
    {
      if(!IsEventSelected(fAOD))
      {
        return;
      }
    }
  }

  if(fPP || fPPb)
  {
    if(!IsEventSelectedPP(fAOD))
      return;
  }

  //printf("Event selection passed\n");

  // from here only selected events survive

  FillEventQA(fAOD,1); // after cuts events QA

  // loading PID response for protons
  if(fTracksScan || fPID || ( (fCutV0ProtonNumSigmaMax > 0.) && (fCutV0ProtonPIDPtMax > 0.) ) )
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse)
    {
      ::Error("UserExec","No PID response object found");
      return;
    }

    fTPCPIDResponse = fPIDResponse->GetTPCResponse();
  }

  if(fSampling) // randomly assign sampling bin index for stat. uncertanity estimation
  {
    TRandom3 rr(0);
    Double_t ranNum = rr.Rndm();

    if(fGFKNumSamples == 5)
    {
      // sampling for 5 samples
      if (ranNum <= 0.2) fSampleBinIndex = 0;
      else if (ranNum <= 0.4 && ranNum > 0.2) fSampleBinIndex = 1;
      else if (ranNum <= 0.6 && ranNum > 0.4) fSampleBinIndex = 2;
      else if (ranNum <= 0.8 && ranNum > 0.6) fSampleBinIndex = 3;
      else fSampleBinIndex = 4;
    }
    else
    {
      // sampling for 10 samples
      if (ranNum <= 0.1) fSampleBinIndex = 0;
      else if (ranNum <= 0.2 && ranNum > 0.1) fSampleBinIndex = 1;
      else if (ranNum <= 0.3 && ranNum > 0.2) fSampleBinIndex = 2;
      else if (ranNum <= 0.4 && ranNum > 0.3) fSampleBinIndex = 3;
      else if (ranNum <= 0.5 && ranNum > 0.4) fSampleBinIndex = 4;
      else if (ranNum <= 0.6 && ranNum > 0.5) fSampleBinIndex = 5;
      else if (ranNum <= 0.7 && ranNum > 0.6) fSampleBinIndex = 6;
      else if (ranNum <= 0.8 && ranNum > 0.7) fSampleBinIndex = 7;
      else if (ranNum <= 0.9 && ranNum > 0.8) fSampleBinIndex = 8;
      else fSampleBinIndex = 9;
    }
  }
  else
  {
    fSampleBinIndex = 0;
  }


  const Int_t iTracks(fAOD->GetNumberOfTracks());
  fEventMult->Fill(iTracks);

  // === tracks & PID species filtering ====

  // tracks scanning (pt,eta,phi) of various FBs
  if(fTracksScan)
    ScanTracks();


  // cleaning all TClonesArray containers
  fArrTracksFiltered.Clear("C");
  if(fPID)
  {
    fArrPionFiltered.Clear("C");
    fArrKaonFiltered.Clear("C");
    fArrProtonFiltered.Clear("C");
    if(fDoV0s)
    {
      fArrV0sK0sFiltered.Clear("C");
      fArrV0sLambdaFiltered.Clear("C");
      fArrV0sALambdaFiltered.Clear("C");
    }
  }

  // filtering objects and filling relevant containers
  FilterTracks(); // filter reference tracks

  if(fPID)
  {
    if(fUseBayesPID) // filter PID tracks Using Bayes PID
    {
      FilterPIDTracksBayesPID();
    }
    else
    {
      FilterPIDTracks(); // filter PID tracks (pi,K,p)
    }

    if(fDoV0s)
      FilterV0s(); // filter V0s
  }

  // all the possible tracks & candidates (passing selection criteria) are filtered and filled and relevant TClonesArray

  // put it here since in pp analysus, fCentBinIndex is returned durinf track filtering
  fSampleCounter->Fill(fSampleBinIndex,fCentBinIndex);

  // flow analysis
  if(fDoFlow)
  {

    if(fDoGenFramKat) // do flow according to Gen Framework implemented by Kat.
    {
      DoGenFramKatarina();
    }

    if(fOldFlow)
    {
      // ==== FLOW ANALYSIS ====
      // begin of loop over EtaGap & Harmonics
      Int_t iHarm = 0;
      Double_t dEtaGap = 0;
      TProfile2D* profV0sPos = 0x0;
      TProfile2D* profV0sNeg = 0x0;


      for(Int_t indexHarm(0); indexHarm < fNumHarmonics; indexHarm++)
      {
        for(Int_t indexEtaGap(0); indexEtaGap < fNumEtaGap; indexEtaGap++)
        {
          iHarm = fHarmonics[indexHarm];
          dEtaGap = fEtaGap[indexEtaGap];

          FillRefFlowVectors(dEtaGap, iHarm);
          EstimateRefCumulant(dEtaGap, iHarm, fTracksRefTwo[indexHarm][indexEtaGap][fSampleBinIndex]);

          if(fDiffFlow)
          {
            // tracks
            EstimatePtDiffCumulant(fArrTracksFiltered,dEtaGap,iHarm,fTracksDiffTwoPos[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex],fTracksDiffTwoNeg[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex]);
            // PID tracks
            EstimatePtDiffCumulant(fArrPionFiltered,dEtaGap,iHarm,fPionsDiffTwoPos[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex],fPionsDiffTwoNeg[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex]);
            EstimatePtDiffCumulant(fArrKaonFiltered,dEtaGap,iHarm,fKaonsDiffTwoPos[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex],fKaonsDiffTwoNeg[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex]);
            EstimatePtDiffCumulant(fArrProtonFiltered,dEtaGap,iHarm,fProtonsDiffTwoPos[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex],fProtonsDiffTwoNeg[fCentBinIndex][indexHarm][indexEtaGap][fSampleBinIndex]);
          }

          if(fPID)
          {
            if(fDoV0s)
              EstimateV0Cumulant(indexEtaGap,indexHarm,fSampleBinIndex);
          }
        } // end of loop over eta gap
      } // end of loop over harmonics
    } // end of doOldFlow

  } // end of if(fDoFlow)


  PostData(1, fOutListCumulants);
  PostData(2, fOutListEvents);
  PostData(3, fOutListTracks);
  PostData(4, fOutListPID);
  PostData(5, fOutListV0s);
  //PostData(6, fOutListQA);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::Terminate(Option_t *)
{


  // terminate
  // called at the END of the analysis (when all events are processed)

}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::ListParameters()
{
  printf("\n======= List of input parameters ======\n");

  printf("--------- Analysis & events setters --------\n");
  printf("SetAODAnalysis(Bool_t %d)\n",fAODAnalysis);
  printf("SetPbPbAnalysis(Bool_t %d)\n",fPbPb);
  printf("SetPPAnalysis(Bool_t %d)\n",fPP);
  printf("SetPPbAnalysis(Bool_t %d)\n",fPPb);
  printf("SetTrigger(Short_t %d)\n",fTrigger);
  printf("SetPeriod10h(Bool_t %d)\n",fLHC10h);
  printf("SetRejectPileUpSPD(Bool_t %d)\n",fRejectPileFromSPD);
  printf("SetUsePlpMV(Bool_t %d)\n",fUsePlpMV);
  printf("SetRejectOutOfBunchPU(Bool_t %d)\n",fRejectOutOfBunchPU);
  printf("SetUseIsPileUpFromSPD(Bool_t %d)\n",fUseIsPileUpFromSPD);
  printf("SetCentFlag(Short_t %d)\n",fCentFlag);
  printf("SetPVtxZMax(Double_t %f)\n",fPVtxCutZ);
  printf("SetSampling(Bool_t %d)\n",fSampling);
  printf("SetDoFlow(Bool_t %d)\n",fDoFlow);
  printf("SetDiffFlow(Bool_t %d)\n",fDiffFlow);
  printf("SetPID(Bool_t %d)\n",fPID);
  printf("SetUseBayesPID(Bool_t %d)\n",fUseBayesPID);
  printf("SetTracksScan(Bool_t %d)",fTracksScan);
  printf("SetDoV0s(Bool_t %d)\n",fDoV0s);
  printf("SetDoFlowGenFramKatarina(Bool_t %d)\n",fDoGenFramKat);
  printf("SetDoOldFlow(Bool_t %d)\n",fOldFlow);
  printf("SetUseOldCent(Bool_t %d)\n",fUseOldCent);

  printf("--------- Tracks & PID setters --------\n");
  printf("SetTrackEtaMax(Double_t %f)\n",fTrackEtaMax);
  printf("SetTrackPtMax(Double_t %f)\n", fTrackPtMax);
  printf("SetTrackPtMin(Double_t %f)\n", fTrackPtMin);
  printf("SetNumTPCclsMin(UShort_t %d)\n", fNumTPCclsMin);
  printf("SetTrackFilterBit(UInt_t %d)\n", fTrackFilterBit);
  printf("SetPionNumSigmasMax(Double_t %f)\n", fCutPionNumSigmaMax);
  printf("SetKaonNumSigmasMax(Double_t %f)\n", fCutKaonNumSigmaMax);
  printf("SetProtonNumSigmasMax(Double_t %f)\n", fCutProtonNumSigmaMax);
  printf("SetPIDBayesProbPionMin(Double_t %f)\n", fCutBayesPIDPionMin);
  printf("SetPIDBayesProbKaonMin(Double_t %f)\n", fCutBayesPIDKaonMin);
  printf("SetPIDBayesProbProtonMin(Double_t %f)\n", fCutBayesPIDProtonMin);

  printf("--------- V0s setters --------\n");
  printf("SetV0sOnFly(Bool_t %d)\n", fCutV0onFly);
  printf("SetV0sTPCRefit(Bool_t %d)\n", fCutV0refitTPC);
  printf("SetV0sRejectKinks(Bool_t %d)\n", fCutV0rejectKinks);
  printf("SetV0sDCAPVMin(Double_t %f)\n", fCutV0MinDCAtoPV);
  printf("SetV0sDCAPVMax(Double_t %f)\n", fCutV0MaxDCAtoPV);
  printf("SetV0sDCADaughtersMax(Double_t %f)\n", fCutV0MaxDCADaughters);
  printf("SetV0sDecayRadiusMin(Double_t %f)\n", fCutV0MinDecayRadius);
  printf("SetV0sDecayRadiusMax(Double_t %f)\n", fCutV0MaxDecayRadius);
  printf("SetV0sDaughterPtMin(Double_t %f)\n", fCutV0DaughterPtMin);
  printf("SetV0sDaughterEtaMax(Double_t %f)\n", fCutV0DaughterEtaMax);
  printf("SetV0sMotherEtaMax(Double_t %f)\n",  fCutV0MotherEtaMax);
  printf("SetV0sMotherRapMax(Double_t %f)\n", fCutV0MotherRapMax);
  printf("SetV0sMotherPtMin(Double_t %f)\n", fCutV0MotherPtMin);
  printf("SetV0sMotherPtMax(Double_t %f)\n", fCutV0MotherPtMax);
  printf("SetV0sK0sCPAMin(Double_t %f)\n", fCutV0MinCPAK0s);
  printf("SetV0sLambdaCPAMin(Double_t %f)\n", fCutV0MinCPALambda);
  printf("SetV0sK0sNumTauMax(Double_t %f)\n", fCutV0NumTauK0sMax);
  printf("SetV0sLambdaNumTauMax(Double_t %f)\n", fCutV0NumTauLambdaMax);
  printf("SetV0sK0sArmenterosAlphaMin(Double_t %f)\n", fCutV0K0sArmenterosAlphaMin);
  printf("SetV0sProtonNumSigmaMax(Double_t %f)\n", fCutV0ProtonNumSigmaMax);
  printf("SetV0sProtonPIDPtMax(Double_t %f)\n", fCutV0ProtonPIDPtMax);
  printf("====================================\n\n");
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::EstimatePtDiffCumulant(TClonesArray &array, const Float_t dEtaGap, const Short_t iHarm,TProfile* profilePos, TProfile* profileNeg)
{
  // checking if pointers are valid
  if(!profilePos)
    return;

  if(!profileNeg && !(dEtaGap < 0.)) // if with eta gap, negative profile has to be valid as well
    return;

  // check if RFPs flow vectors are empty
  if(!AreRefFlowVectorsFilled(dEtaGap,iHarm)) // if yes fill ref flow vectors
    FillRefFlowVectors(dEtaGap, iHarm);

  // POIs
  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for POI
  const Int_t iNumPOIs = array.GetEntries();

  TComplex vecPpos[fNumPtBins]; // flow vector for POIs in positive eta (including eta gap cut) [or all POIs with no eta gap]
  TComplex vecPneg[fNumPtBins]; // flow vector for POIs in negative eta (including eta gap cut)
  Int_t iCountPpos[fNumPtBins]; // counter of POIs in positive eta [or all POIs with no eta gap]
  Int_t iCountPneg[fNumPtBins]; // counter of POIs in negative eta

  AliAODTrack* track = 0x0;
  Double_t dPhi = 0, dEta = 0;
  Double_t dWeight = 0, dAmp = 0, dVal = 0;
  Short_t iPtBinIndex = -1;

  for(Int_t i(0); i < fNumPtBins; i++) // 0-ing staff
  {
    vecPpos[i] = TComplex(0,0,kFALSE);
    vecPneg[i] = TComplex(0,0,kFALSE);
    iCountPpos[i] = 0;
    iCountPneg[i] = 0;
  }

  for(Int_t i(0); i < iNumPOIs; i++) // loop over POIs
  {
    track = static_cast<AliAODTrack*>(array.At(i));

    if(!track)
      continue;

    iPtBinIndex = GetPtBinIndex(track->Pt());

    if(iPtBinIndex == -1)
      continue;

    dPhi = track->Phi();
    dEta = track->Eta();

    // filling the flow vectors for POIs
    if(dEtaGap < 0.) // no eta gap
    {
      vecPpos[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
      iCountPpos[iPtBinIndex]++;
    }
    else // with eta gap
    {
      if(dEta > dEtaGapCut)
      {
        vecPpos[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
        iCountPpos[iPtBinIndex]++;
      }

      if(dEta < -dEtaGapCut)
      {
        vecPneg[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
        iCountPneg[iPtBinIndex]++;
      }
    }
  } // end of loop over POIs entries

  // Filling the profiles
  for(Int_t i(0); i < fNumPtBins; i++)
  {
    if(dEtaGap < 0.) // with no eta gap
    {
      dWeight = iCountPpos[i] * (fCountRefPos - 1);
      dAmp = (vecPpos[i] * (TComplex::Conjugate(fVecRefPos))).Re();
      dVal = (dAmp - iCountPpos[i]) / dWeight;
      if( TMath::Abs(dVal < 1) && (dWeight > 1))
        profilePos->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal, dWeight);
    }
    else // with eta gap
    {
      dWeight = iCountPpos[i] * (fCountRefNeg);
      dAmp = (vecPpos[i] * (TComplex::Conjugate(fVecRefNeg))).Re();
      dVal = dAmp / dWeight;
      if( TMath::Abs(dVal < 1) && (dWeight > 0))
        profilePos->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal, dWeight);

      dWeight = iCountPneg[i] * (fCountRefPos);
      dAmp = (vecPneg[i] * (TComplex::Conjugate(fVecRefPos))).Re();
      dVal = dAmp / dWeight;
      if( TMath::Abs(dVal < 1) && (dWeight > 0))
        profileNeg->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, dVal, dWeight);
    }
  } // end of loop over pt bins
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::EstimateV0Cumulant(const Short_t iEtaGapIndex, const Short_t iHarmonicsIndex, const Short_t iSampleIndex)
{
  const Int_t iHarm = fHarmonics[iHarmonicsIndex];
  const Double_t dEtaGap = fEtaGap[iEtaGapIndex];
  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for POI

  // check if RFPs flow vectors are empty
  if(!AreRefFlowVectorsFilled(dEtaGap,iHarm)) // if yes fill ref flow vectors
    FillRefFlowVectors(dEtaGap, iHarm);

  TComplex vecVpos[fNumPtBins]; // flow vector for POIs in positive eta (including eta gap cut) [or all POIs with no eta gap]
  TComplex vecVneg[fNumPtBins]; // flow vector for POIs in negative eta (including eta gap cut)
  Int_t iCountVpos[fNumPtBins]; // counter of POIs in positive eta [or all POIs with no eta gap]
  Int_t iCountVneg[fNumPtBins]; // counter of POIs in negative eta

  AliAODv0* v0 = 0x0;
  Double_t dPhi = 0, dEta = 0, dMass = 0;
  Double_t dWeight = 0, dAmp = 0, dVal = 0;
  Short_t iMinvFlowBinIndex = -1, iPtBinIndex = -1;

  TClonesArray arrayV0s = TClonesArray("AliAODv0",5000);
  Int_t iNumV0s = 0;
  TProfile2D* profV0sPos;
  TProfile2D* profV0sNeg;
  TH2D* hPtInvMass;
  TH1D* hInvMass;
  Int_t iNumMinvFlowBins = 0;
  Double_t* dInvMassBinEdges;

  for(Int_t k(0); k < 2; k++) // loop over particle species (0: K0s / 1: Lambda / 2: ALambda)
  {
    switch(k)
    {
      case 0: // K0s
        arrayV0s = fArrV0sK0sFiltered;
        profV0sPos = fV0sDiffTwoPos_K0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex][iSampleIndex];
        profV0sNeg = fV0sDiffTwoNeg_K0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex][iSampleIndex];
        hPtInvMass = fV0sPtInvMassK0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hInvMass = fV0sInvMassK0s[iEtaGapIndex];
        iNumMinvFlowBins = fNumMinvFlowBinsK0s;
        dInvMassBinEdges = fMinvFlowBinEdgesK0s;
      break;

      case 1: // Lambda
        arrayV0s = fArrV0sLambdaFiltered;
        profV0sPos = fV0sDiffTwoPos_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex][iSampleIndex];
        profV0sNeg = fV0sDiffTwoNeg_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex][iSampleIndex];
        hPtInvMass = fV0sPtInvMassLambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hInvMass = fV0sInvMassLambda[iEtaGapIndex];
        iNumMinvFlowBins = fNumMinvFlowBinsLambda;
        dInvMassBinEdges = fMinvFlowBinEdgesLambda;
      break;

      case 2: // Anti-Lambda
        arrayV0s = fArrV0sALambdaFiltered;
        profV0sPos = fV0sDiffTwoPos_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex][iSampleIndex];
        profV0sNeg = fV0sDiffTwoNeg_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex][iSampleIndex];
        hPtInvMass = fV0sPtInvMassLambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hInvMass = fV0sInvMassLambda[iEtaGapIndex];
        iNumMinvFlowBins = fNumMinvFlowBinsLambda;
        dInvMassBinEdges = fMinvFlowBinEdgesLambda;
      break;

      default:
      break;
    }

    iNumV0s = arrayV0s.GetEntries();

    // checking if pointers are valid
    if(!profV0sPos)
      return;

    if(!profV0sNeg && !(dEtaGap < 0.0)) // if with eta gap, negative profile has to be valid as well
      return;

    for(Int_t j(0); j < iNumMinvFlowBins; j++) // loop over Minv bins
    {
      for(Int_t i(0); i < fNumPtBins; i++) // 0-ing staff
      {
        vecVpos[i] = TComplex(0,0,kFALSE);
        vecVneg[i] = TComplex(0,0,kFALSE);
        iCountVpos[i] = 0;
        iCountVneg[i] = 0;
      }

      for(Int_t i(0); i < iNumV0s; i++) // loop over POIs
      {
        v0 = static_cast<AliAODv0*>(arrayV0s.At(i));

        if(!v0)
          continue;

        switch(k)
        {
          case 0: // K0s
            dMass = v0->MassK0Short();
            iMinvFlowBinIndex = GetMinvFlowBinIndexK0s(dMass);
          break;

          case 1: // Lambda
            dMass = v0->MassLambda();
            iMinvFlowBinIndex = GetMinvFlowBinIndexLambda(dMass);
          break;

          case 2: // Anti-Lambda
            dMass = v0->MassAntiLambda();
            iMinvFlowBinIndex = GetMinvFlowBinIndexLambda(dMass);
          break;

          default:
          break;
        }

        iPtBinIndex = GetPtBinIndex(v0->Pt());

        if(iPtBinIndex == -1 || iMinvFlowBinIndex != j)
          continue;

        dPhi = v0->Phi();
        dEta = v0->Eta();

        // filling the flow vectors for POIs
        if(dEtaGap < 0.) // no eta gap
        {
          vecVpos[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
          iCountVpos[iPtBinIndex]++;
          hPtInvMass->Fill(v0->Pt(), dMass);
          hInvMass->Fill(dMass);
        }
        else // with eta gap
        {
          if(dEta > dEtaGapCut)
          {
            vecVpos[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
            iCountVpos[iPtBinIndex]++;
            hPtInvMass->Fill(v0->Pt(), dMass);
            hInvMass->Fill(dMass);
          }

          if(dEta < -dEtaGapCut)
          {
            vecVneg[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
            iCountVneg[iPtBinIndex]++;
            hPtInvMass->Fill(v0->Pt(), dMass);
            hInvMass->Fill(dMass);
          }
        }
      } // end of loop over POIs entries

      // filling the TProfiles
      for(Int_t i(0); i < fNumPtBins; i++)
      {
        if(dEtaGap < 0.) // with no eta gap
        {
          dWeight = iCountVpos[i] * fCountRefPos;
          dAmp = (vecVpos[i] * (TComplex::Conjugate(fVecRefPos))).Re();
          dVal = dAmp / dWeight;
          if( TMath::Abs(dVal < 1) && (dWeight > 0))
            profV0sPos->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, (dInvMassBinEdges[j+1] + dInvMassBinEdges[j])/2 ,dVal, dWeight);
        }
        else // with eta gap
        {
          dWeight = iCountVpos[i] * (fCountRefNeg);
          dAmp = (vecVpos[i] * (TComplex::Conjugate(fVecRefNeg))).Re();
          dVal = dAmp / dWeight;

          if( TMath::Abs(dVal < 1) && (dWeight > 0))
            profV0sPos->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, (dInvMassBinEdges[j+1] + dInvMassBinEdges[j])/2 ,dVal, dWeight);


          dWeight = iCountVneg[i] * (fCountRefPos);
          dAmp = (vecVneg[i] * (TComplex::Conjugate(fVecRefPos))).Re();
          dVal = dAmp / dWeight;

          if( TMath::Abs(dVal < 1) && (dWeight > 0))
            profV0sNeg->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, (dInvMassBinEdges[j+1] + dInvMassBinEdges[j])/2 ,dVal, dWeight);

        }
      }
    } // end of loop over Minv bins
  }// end of loop over particle species
  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::AreRefFlowVectorsFilled(const Float_t dEtaGap, const Short_t iHarm)
{
  if(dEtaGap != fEtaCutFlag || iHarm != fHarmFlag)
    return kFALSE;

  if(fVecRefPos.Re() == 0 && fVecRefPos.Im() == 0)
    return kFALSE;

  if(dEtaGap != -1 && fVecRefNeg.Re() == 0 && fVecRefNeg.Im() == 0) // if eta gap == -1 -> no eta gap only fVecRefPos & fCountRefPos are filled
    return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FillRefFlowVectors(const Float_t dEtaGap, const Short_t iHarm)
{
  //::Info("FillRefFlowVectors",Form("Filling RFPs flow vectors with eta gap %g & harmonics %d",dEtaGap,iHarm));

  fVecRefPos = TComplex(0,0,kFALSE); // flow vector for RFPs in positive eta (including eta gap cut) [or all RFPs with no eta gap]
  fVecRefNeg = TComplex(0,0,kFALSE); // flow vector for RFPs in negative eta (including eta gap cut)
  fCountRefPos = 0; // counter of POIs in positive eta [or all RFPs with no eta gap]
  fCountRefNeg = 0; // counter of POIs in negative eta

  // setting value of flags
  fEtaCutFlag = dEtaGap;
  fHarmFlag = iHarm;

  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for RFPs
  const Int_t iNumEntries = fArrTracksFiltered.GetEntries(); // number of entries in RFPs array

  AliAODTrack* track = 0x0;
  Double_t dPhi = 0, dEta = 0;

  for(Int_t i(0); i < iNumEntries; i++) // loop over RFPs
  {
    track = static_cast<AliAODTrack*>(fArrTracksFiltered.At(i));

    if(!track)
      continue;

    dPhi = track->Phi();
    dEta = track->Eta();

    if(dEtaGap < 0.) // with no gap -> only positive quantities are filled
    {
      //::Info("FillRefFlowVectors","Filling with NO eta gap");
      fVecRefPos += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
      fCountRefPos++;
    }
    else if(!(dEtaGap < 0.)) // with eta gap
    {
      //::Info("FillRefFlowVectors",Form("Filling with eta gap %f",dEtaGap));
      if(dEta > dEtaGapCut)
      {
        fVecRefPos += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
        fCountRefPos++;
      }

      if(dEta < -dEtaGapCut)
      {
        fVecRefNeg += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
        fCountRefNeg++;
      }
    }
  } // end of loop over RFPs entries

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::EstimateRefCumulant(const Float_t dEtaGap, const Short_t iHarm, TProfile* profile)
{
  if(!profile) // invalid TProfile pointer
    return;

  if(!AreRefFlowVectorsFilled(dEtaGap,iHarm))
    FillRefFlowVectors(dEtaGap,iHarm);

  Double_t dWeight = 0, dAmp = 0, dVal = 0;

  if(dEtaGap < 0.) // no eta gap
  {
    dWeight = fCountRefPos * (fCountRefPos - 1);
    dAmp = (fVecRefPos * (TComplex::Conjugate(fVecRefPos))).Re();
    dVal = (dAmp - fCountRefPos) / dWeight;
    if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCentPercentile > 0) )
      profile->Fill(fCentPercentile, dVal, dWeight);
  }
  else // with eta gap
  {
    dWeight = fCountRefPos * fCountRefNeg;
    dAmp = (fVecRefPos * (TComplex::Conjugate(fVecRefNeg))).Re();
    dVal = (dAmp) / dWeight;
    if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCentPercentile > 0) )
      profile->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::OldIsEventSelected(const AliAODEvent* event)
{
  //printf("Inside OldIsEventSelected\n");
  //GetVertex()
  Float_t vtxz = -999.;

  if (!fAODAnalysis)
  {
    return kFALSE;
  }

  AliAODEvent* aod = (AliAODEvent*)event;

  const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0)
      return kFALSE;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks"))
      return kFALSE;
  const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
  if (!spdVtx || spdVtx->GetNContributors()<=0)
      return kFALSE;
  if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5)
      return kFALSE;

  vtxz = trkVtx->GetZ();


  if(vtxz < -990)
  {
    //printf("Vtx 1\n");
    return kFALSE;
  }
  else
  {
    //fVtx->Fill(1);
    //fVtxBeforeCuts->Fill(zvtx);
    if (TMath::Abs(vtxz) < fPVtxCutZ)
    {
      //fMultCorBeforeCuts->Fill(GetGlobalMult(fAOD), GetTPCMult(fAOD));
      if(fRejectPileFromSPD)
      {

        if(fUseIsPileUpFromSPD)
        {
          Short_t isPileup = fAOD->IsPileupFromSPD(3);
          if (isPileup != 0)
            return kFALSE;
        }

        //if (fAOD->GetHeader()->GetRefMultiplicityComb08() < 0)
        //    return;

        //New cut from Ruben for pileup hybrid
        if (fUsePlpMV && plpMV(fAOD))
            return kFALSE;


        if ( (fLHC10h) && ((Float_t(GetTPCMult(fAOD)) < (-40.3+1.22*GetGlobalMult(fAOD))) || (Float_t(GetTPCMult(fAOD)) > (32.1+1.59*GetGlobalMult(fAOD))) ) )
            return kFALSE;

        //if ( (!fLHC10h) && ( (Float_t(GetTPCMult(fAOD)) < (-36.73+1.48*GetGlobalMult(fAOD))) || (Float_t(GetTPCMult(fAOD)) > (62.87+1.78*GetGlobalMult(fAOD))) ) )
            //return kFALSE;

      }

      //printf("Vtx & pileup passed\n");

      // GetCentrCode()
      AliAODEvent* aod = (AliAODEvent*)event;
      AliCentrality* centrality = ((AliAODHeader*)aod->GetHeader())->GetCentralityP();

      Short_t centrCode = -1;
      Float_t cent = -1.;

      if (fCentFlag == 0)
        cent = centrality->GetCentralityPercentile("V0M");
      else if (fCentFlag == 1)
        cent = centrality->GetCentralityPercentile("TRK");
      else if (fCentFlag == 2)
        cent = centrality->GetCentralityPercentile("CL1");

      Float_t centTrk = centrality->GetCentralityPercentile("TRK");
      Float_t centV0 = centrality->GetCentralityPercentile("V0M");

      if (fLHC10h)
      {
        if (TMath::Abs(centV0 - centTrk) < 5.0 && cent <= 80 && cent > 0)
        {
          if ((cent > 0) && (cent <= 5.0))
            centrCode = 0;
          else if ((cent > 5.0) && (cent <= 10.0))
            centrCode = 1;
          else if ((cent > 10.0) && (cent <= 20.0))
            centrCode = 2;
          else if ((cent > 20.0) && (cent <= 30.0))
            centrCode = 3;
          else if ((cent > 30.0) && (cent <= 40.0))
            centrCode = 4;
          else if ((cent > 40.0) && (cent <= 50.0))
            centrCode = 5;
          else if ((cent > 50.0) && (cent <= 60.0))
            centrCode = 6;
          else if ((cent > 60.0) && (cent <= 70.0))
            centrCode = 7;
          else if ((cent > 70.0) && (cent <= 80.0))
            centrCode = 8;
        }
      }

      if (!fLHC10h)
      {
        if (cent <= 80 && cent > 0)
        {
          if ((cent > 0) && (cent <= 5.0))
             centrCode = 0;
          else if ((cent > 5.0) && (cent <= 10.0))
             centrCode = 1;
          else if ((cent > 10.0) && (cent <= 20.0))
             centrCode = 2;
          else if ((cent > 20.0) && (cent <= 30.0))
             centrCode = 3;
          else if ((cent > 30.0) && (cent <= 40.0))
             centrCode = 4;
          else if ((cent > 40.0) && (cent <= 50.0))
             centrCode = 5;
          else if ((cent > 50.0) && (cent <= 60.0))
             centrCode = 6;
          else if ((cent > 60.0) && (cent <= 70.0))
             centrCode = 7;
          else if ((cent > 70.0) && (cent <= 80.0))
            centrCode = 8;
        }
      }

      Short_t cenAOD = centrCode;

      if (cenAOD >= 0)
      {
        fCentralityDis->Fill(fCentPercentile);
        fCentDistUnitBin->Fill(fCentPercentile);
        fEventCounter->Fill(6);

        fCentPercentile = cent;
        fCentBinIndex = cenAOD;

        return kTRUE;
      }
    }
  }

  return kFALSE;
}
//____________________________________________________________________
Int_t AliAnalysisTaskFlowPID::GetTPCMult(AliVEvent* ev) const
{

  // function to return the tpc only multiplicity as in the flow package; it is used to cut on the event multiplcity to remove outliers in the centrality

  Int_t multTPC = 0;

  if (fAODAnalysis) {

    AliAODEvent* aod = (AliAODEvent*)ev;
    const Int_t nGoodTracks = aod->GetNumberOfTracks();
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
      AliAODTrack* trackAOD =(AliAODTrack*) aod->GetTrack(iTracks);

      if (!trackAOD)
  continue;

      if (!(trackAOD->TestFilterBit(1)))
  continue;

      if ((trackAOD->Pt() < fTrackPtMin) || (trackAOD->Pt() > fTrackPtMax) || (TMath::Abs(trackAOD->Eta()) > fTrackEtaMax) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2))
  continue;

      multTPC++;
    }//track loop

  }//analysis type

  return multTPC;
}

//____________________________________________________________________
Int_t AliAnalysisTaskFlowPID::GetGlobalMult(AliVEvent* ev) const
{

  // function to return the global multiplicity as in the flow package; it is used to cut on the event multiplcity to remove outliers in the centrality

  Int_t multGlobal = 0;
  if (fAODAnalysis) {

    AliAODEvent* aod = (AliAODEvent*)ev;
    const Int_t nGoodTracks = aod->GetNumberOfTracks();
    for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
      AliAODTrack* trackAOD = (AliAODTrack*)aod->GetTrack(iTracks);

      if (!trackAOD)
  continue;

      if (!(trackAOD->TestFilterBit(16)))
  continue;

      if ((trackAOD->Pt() < fTrackPtMin) || (trackAOD->Pt() > fTrackPtMax) || (TMath::Abs(trackAOD->Eta()) > fTrackEtaMax) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1))
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
Bool_t AliAnalysisTaskFlowPID::IsEventSelected(const AliAODEvent* event)
{
  //printf("Inside IsEventSelected\n");
  // event selection criteria

  // Pileup rejection
  if( fRejectPileFromSPD && event->IsPileupFromSPD(3) ) //min contributors ???
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
    EstimateCentrality(fAOD);
    if (fCentBinIndex < 0 || fCentBinIndex > fNumCentBins || fCentPercentile < 0 || fCentPercentile > 80)
      return kFALSE;
  }

  fCentralityDis->Fill(fCentPercentile);
  fCentDistUnitBin->Fill(fCentPercentile);
  fEventCounter->Fill(6);

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsEventSelectedPP(AliVEvent* event)
{
//  TString sEventCounterLabel[] = {"Input","AOD OK","Pile-up OK","PV OK","SPD Vtx OK","PV #it{z} OK","Centrality OK","At least 2 selected tracks","At least 1 V0 candidate"};
  if(!fPP && !fPPb) // not pp or pPb analysis
    return kFALSE;

  if(!event) // not valid event
    return kFALSE;


  if(!fAODAnalysis) // not impelemted for ESD
    return kFALSE;

  //fEventCounter->Fill("AOD OK",1);

  // trigggers QA check present : suggested by Evgeny
  TString firedClasses = fInputEvent->GetFiredTriggerClasses();
  Bool_t isINT7 = firedClasses.Contains("CINT7-B-NOPF-CENT");
  Bool_t isV0HM = firedClasses.Contains("CVHMV0M-B-SPD2-CENT");

  if(isINT7)
  {
    fQAEventsTriggers[0]->Fill("CINT7-B-NOPF-CENT",1);
  }

  if(isV0HM)
  {
    fQAEventsTriggers[0]->Fill("CVHMV0M-B-SPD2-CENT",1);
  }

  // event / physics selection criteria
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();

  if(fSelectMask& AliVEvent::kINT7)
  {
    fQAEventsTriggers[0]->Fill("kINT7",1);
  }

  if(fSelectMask& AliVEvent::kHighMultV0)
  {
    fQAEventsTriggers[0]->Fill("kHighMultV0",1);
  }

  Bool_t isTriggerSelected = kFALSE;

  switch(fTrigger) // check for high multiplicity trigger
  {
    case 0:
      isTriggerSelected = fSelectMask& AliVEvent::kINT7;
      break;

    case 1:
      isTriggerSelected = fSelectMask& AliVEvent::kHighMultV0;
      break;

    case 2: isTriggerSelected = fSelectMask& AliVEvent::kHighMultSPD; break;
    default: isTriggerSelected = kFALSE;
  }

  if(!isTriggerSelected)
    return kFALSE;

  if(isINT7)
  {
    fQAEventsTriggers[1]->Fill("CINT7-B-NOPF-CENT",1);
  }

  if(isV0HM)
  {
    fQAEventsTriggers[1]->Fill("CVHMV0M-B-SPD2-CENT",1);
  }

  if(fSelectMask& AliVEvent::kINT7)
  {
    fQAEventsTriggers[1]->Fill("kINT7",1);
  }

  if(fSelectMask& AliVEvent::kHighMultV0)
  {
    fQAEventsTriggers[1]->Fill("kHighMultV0",1);
  }

  // events passing physics selection

  fEventCounter->Fill("Physics selection OK",1);

  // primary vertex selection
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(event->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;

  fEventCounter->Fill("PV OK",1);

  // SPD vertex selection
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(event->GetPrimaryVertexSPD());

  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol))
  {
    return kFALSE;
  }

  fEventCounter->Fill("SPD Vtx OK",1);

  // PileUp rejection included in Physics selection
  // but with values for high mult pp (> 5 contrib) => for low ones: do manually (> 3 contrib)

  /*
  if(fTrigger == 0 && fAOD->IsPileupFromSPD(3,0.8) )
  {
    return kFALSE;
  }
  */

  fEventCounter->Fill("Pileup SPD OK",1);

  // pileup rejection from multivertexer
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5);
  utils.SetMaxPlpChi2MV(5);
  utils.SetMinWDistMV(15);
  utils.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t isPileupFromMV = utils.IsPileUpMV(event);

  if(isPileupFromMV)
  {
    return kFALSE;
  }

  fEventCounter->Fill("Pileup MV OK",1);

  if(fRejectOutOfBunchPU) // out-of-bunch rejection (provided by Christian)
  {
    //out-of-bunch 11 BC
    if (utils.IsOutOfBunchPileUp(event))
    {
      return kFALSE;
    }
    fEventCounter->Fill("OOBPU OK",1);

    if (utils.IsSPDClusterVsTrackletBG(event))
    {
      return kFALSE;
    }

    fEventCounter->Fill("SPDClTrBG OK",1);

    // SPD pileup
    if (utils.IsPileUpSPD(event))
    {
      return kFALSE;
    }

    fEventCounter->Fill("SPDPU OK",1);
  }

  fEventCounter->Fill("Utils OK",1);

  // cutting on PV z-distance
  const Double_t aodVtxZ = vtx->GetZ();
  if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
  {
    return kFALSE;
  }
  fEventCounter->Fill("PV #it{z} OK",1);

  // centrality estimation (pPb analysis)

  if(fPPb)
  {
    Float_t lPercentile = 300;
    AliMultSelection* MultSelection = 0x0;
    MultSelection = (AliMultSelection*) event->FindListObject("MultSelection");

    if( !MultSelection )
    {
      //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
      ::Warning("IsEventSelectedPP","AliMultSelection object not found!");
    }
    else
    {
      lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    }

    fCentPercentile = lPercentile;

  }



 	return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsTrackSelected(const AliAODTrack* track, const Bool_t bTestFB)
{
	// track selection procedure excluding FilterBit
  //{"Input","Track OK","FB OK","TPC Ncls OK","Eta OK","Pt OK","DCAz OK","Selected"};

  fTracksCounter->Fill("Input",1);

  if(!track)
  {
    return kFALSE;
  }

  fTracksCounter->Fill("Track OK",1);

  if(!bTestFB)
  {
    if( !track->TestFilterBit(fTrackFilterBit) )
    {
      return kFALSE;
    }

    fTracksCounter->Fill("FB OK",1);
  }

  if(fTrackFilterBit != 2) // not stand-alone track  => no TPC info
  {
    if(track->GetTPCNcls() < fNumTPCclsMin)
    {
      return kFALSE;
    }

    fTracksCounter->Fill("TPC Ncls OK",1);
  }

  if(TMath::Abs(track->Eta()) > fTrackEtaMax)
  {
    return kFALSE;
  }

  fTracksCounter->Fill("Eta OK",1);

  if( (track->Pt() > fTrackPtMax) || (track->Pt() < fTrackPtMin) )
  {
    return kFALSE;
  }

  fTracksCounter->Fill("Pt OK",1);

  if( fTracksDCAzMax > 0. && TMath::Abs(track->ZAtDCA()) > fTracksDCAzMax)
  {
    return kFALSE;
  }

  fTracksCounter->Fill("DCAz OK",1);

  fTracksCounter->Fill("Selected",1);
	return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsTrackPion(const AliAODTrack* track)
{
  if(!fPIDResponse)
  {
    return kFALSE;
  }

  fPionsCounter->Fill(0); // input

  Float_t dNumSigmaTPC = 0;
  Float_t dNumSigmaTOF = 0;

  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);

  if( pidStatusTPC == AliPIDResponse::kDetPidOk)
    fPionsCounter->Fill(1); // TPC PID OK

  if( pidStatusTOF == AliPIDResponse::kDetPidOk )
    fPionsCounter->Fill(2); // TOF PID OK

  if( pidStatusTOF != AliPIDResponse::kDetPidOk || pidStatusTPC != AliPIDResponse::kDetPidOk)
    return kFALSE;

  fPionsCounter->Fill(3); // TOF && TPC PID OK

  dNumSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  dNumSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

  Float_t dNumSigmaCombined = TMath::Sqrt(TMath::Power(dNumSigmaTOF,2) + TMath::Power(dNumSigmaTPC,2));

  if(dNumSigmaCombined >= fCutPionNumSigmaMax)
  {
    return kFALSE;
  }
  fPionsCounter->Fill(4); // number of sigmas

  fPionsCounter->Fill(5); // selected
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsTrackKaon(const AliAODTrack* track)
{
 if(!fPIDResponse)
    return kFALSE;

  fKaonsCounter->Fill(0); // input

  Float_t dNumSigmaTPC = 0;
  Float_t dNumSigmaTOF = 0;

  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);

  if(pidStatusTPC == AliPIDResponse::kDetPidOk)
    fKaonsCounter->Fill(1); // TOF PID OK

  if(pidStatusTOF == AliPIDResponse::kDetPidOk)
    fKaonsCounter->Fill(2); // TOF PID OK

  if( pidStatusTOF != AliPIDResponse::kDetPidOk || pidStatusTPC != AliPIDResponse::kDetPidOk)
    return kFALSE;

  fKaonsCounter->Fill(3); // TOF && TPC PID OK

  dNumSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  dNumSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Float_t dNumSigmaCombined = TMath::Sqrt(TMath::Power(dNumSigmaTOF,2) + TMath::Power(dNumSigmaTPC,2));

  if(dNumSigmaCombined >= fCutKaonNumSigmaMax)
  {
    return kFALSE;
  }

  fKaonsCounter->Fill(4); // number of sigmas

  fKaonsCounter->Fill(5); // selected
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsTrackProton(const AliAODTrack* track)
{
  if(!fPIDResponse)
    return kFALSE;

  fProtonsCounter->Fill(0); // input

  Float_t dNumSigmaTPC = 0;
  Float_t dNumSigmaTOF = 0;

  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);

  if(pidStatusTPC == AliPIDResponse::kDetPidOk)
    fProtonsCounter->Fill(1); // TOF PID OK

  if(pidStatusTOF == AliPIDResponse::kDetPidOk)
    fProtonsCounter->Fill(2); // TOF PID OK

  if( pidStatusTOF != AliPIDResponse::kDetPidOk || pidStatusTPC != AliPIDResponse::kDetPidOk)
    return kFALSE;

  fProtonsCounter->Fill(3); // TOF && TPC PID OK

  dNumSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  dNumSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

  Float_t dNumSigmaCombined = TMath::Sqrt(TMath::Power(dNumSigmaTOF,2) + TMath::Power(dNumSigmaTPC,2));

  if(dNumSigmaCombined >= fCutProtonNumSigmaMax)
  {
    return kFALSE;
  }

  fProtonsCounter->Fill(4); // number of sigmas

  fProtonsCounter->Fill(5); // selected
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::ScanTracks()
{
  if(!fTracksScan)
    return;

  const Int_t iNumTracks = fAOD->GetNumberOfTracks();

  Short_t iFB = 0;
  AliAODTrack* track = 0x0;
  for(Int_t i(0); i < iNumTracks; i++)
  {
    track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track) continue;

    // test filter bit
    for(Short_t j(0); j < fNumScanFB; j++)
    {
      iFB = fTracksScanFB[j];

      if(track->TestFilterBit(iFB))
      {
        fTracksScanPt->Fill(Form("%d",iFB),track->Pt(), 1.);
        fTracksScanPhi->Fill(Form("%d",iFB),track->Phi(),1.);
        fTracksScanEta->Fill(Form("%d",iFB),track->Eta(),1.);

        fTracksScanCounter->Fill(Form("%d",iFB),"All",1.);

        if(IsTrackSelected(track,kTRUE)) // without testing FB for scanning purpose
          fTracksScanCounter->Fill(Form("%d",iFB),"Selected",1.);

        if(!fPIDResponse)
        {
          ::Warning("ScanTracks","No AliPIDResponse object found! Skipping PID scanning");
          return;
        }

        if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,track) == AliPIDResponse::kDetPidOk) // checking TOF
          fTracksScanCounter->Fill(Form("%d",iFB),"ITS PID OK",1.);

        if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,track) == AliPIDResponse::kDetPidOk) // checking TOF
          fTracksScanCounter->Fill(Form("%d",iFB),"TPC PID OK",1.);

        if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track) == AliPIDResponse::kDetPidOk) // checking TOF
          fTracksScanCounter->Fill(Form("%d",iFB),"TOF PID OK",1.);

      }
    }
  } // end of loop over tracks

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FilterTracks()
{
  // Filter all input tracks and clones the ones which passes selection to relevant TClonesArrays
  // Filter both reference particles and PID particles (pi,K,p) / with different filter bits

  const Int_t iNumTracks = fAOD->GetNumberOfTracks();
  Int_t iNumSelected = 0;

  if(iNumTracks < 1)
    return;

  AliAODTrack* track = 0x0;
  for(Int_t i(0); i < iNumTracks; i++)
  {
    if(iNumSelected >= 10000) // checking TClonesArray overflow
      return;

    track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));

    if(!track)
      continue;


    FillTrackQA(track,0); // tracks QA before cuts

    if(IsTrackSelected(track))
    {
      // (ref.) tracks passing all criteria
      new(fArrTracksFiltered[iNumSelected]) AliAODTrack(*track);
      iNumSelected++;
      FillTrackQA(track,1); // tracks QA after cuts
    }
  } // end of loop over all tracks

  fQATracksMult[0]->Fill(iNumTracks);
  fQATracksMult[1]->Fill(iNumSelected);

  // estimating 'centrality' based on multiplicity of charged tracks passing selection criteria
  if(fPP)
  {
    Short_t iMultIndex = GetCentBinIndexPP(iNumSelected);
    if(iMultIndex != -1)
    {
      fCentBinIndex = iMultIndex;
      fCentPercentile = iNumSelected;

      fCentralityDis->Fill(fCentPercentile);
      fCentDistUnitBin->Fill(fCentPercentile);
    }
  }


  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FilterPIDTracks()
{
  const Int_t iNumTracks = fArrTracksFiltered.GetEntriesFast();

  if(iNumTracks == 0) // no filter tracks
    return;

  Int_t iNumPions = 0;
  Int_t iNumKaons = 0;
  Int_t iNumProtons = 0;

  Bool_t bIsEither = kFALSE;

  Double_t dTOF[5] = {0};
  Double_t dBetaTOF = 0;

  AliAODTrack* track = 0x0;

  for(Int_t i(0); i < iNumTracks; i++)
  {
    track = static_cast<AliAODTrack*>(fArrTracksFiltered.At(i));

    bIsEither = kFALSE;

    FillPIDQA(track,0);

    if(IsTrackPion(track))
    {
      bIsEither = kTRUE;
      new(fArrPionFiltered[iNumPions]) AliAODTrack(*track);
      iNumPions++;

      fPionsPt->Fill(track->Pt());
      fPionsEta->Fill(track->Eta());
      fPionsPhi->Fill(track->Phi());
      fPionsNsigmasTPCTOF->Fill( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion) , fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion) );
      fPionsTPCdEdx->Fill(track->P(),track->GetTPCsignal());

      track->GetIntegratedTimes(dTOF);
      dBetaTOF = dTOF[0] / track->GetTOFsignal();
      fPionsTOFbeta->Fill(track->P(), dBetaTOF);
    }

    if(IsTrackKaon(track))
    {
      bIsEither = kTRUE;
      new(fArrKaonFiltered[iNumKaons]) AliAODTrack(*track);
      iNumKaons++;

      fKaonsPt->Fill(track->Pt());
      fKaonsEta->Fill(track->Eta());
      fKaonsPhi->Fill(track->Phi());
      fKaonsNsigmasTPCTOF->Fill( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon) , fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon) );
      fKaonsTPCdEdx->Fill(track->P(),track->GetTPCsignal());

      track->GetIntegratedTimes(dTOF);
      dBetaTOF = dTOF[0] / track->GetTOFsignal();
      fKaonsTOFbeta->Fill(track->P(), dBetaTOF);

    }

    if(IsTrackProton(track))
    {
      bIsEither = kTRUE;
      new(fArrProtonFiltered[iNumProtons]) AliAODTrack(*track);
      iNumProtons++;

      fProtonsPt->Fill(track->Pt());
      fProtonsEta->Fill(track->Eta());
      fProtonsPhi->Fill(track->Phi());
      fProtonsNsigmasTPCTOF->Fill( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton) , fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton) );
      fProtonsTPCdEdx->Fill(track->P(),track->GetTPCsignal());

      track->GetIntegratedTimes(dTOF);
      dBetaTOF = dTOF[0] / track->GetTOFsignal();
      fProtonsTOFbeta->Fill(track->P(), dBetaTOF);
    }

    if(bIsEither) // either PID as pi,K or p
    {
      FillPIDQA(track,1);
    }
  } // end of loop over filtered tracks

  fPionsMult->Fill(iNumPions);
  fKaonsMult->Fill(iNumKaons);
  fProtonsMult->Fill(iNumProtons);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FilterPIDTracksBayesPID()
{
  const Int_t iNumTracks = fArrTracksFiltered.GetEntriesFast();

  if(iNumTracks == 0) // no filter tracks
    return;

  Int_t iNumPions = 0;
  Int_t iNumKaons = 0;
  Int_t iNumProtons = 0;

  Bool_t bIsEither = kFALSE;

  Double_t dTOF[5] = {-900};
  Double_t dBetaTOF = 0.0001;

  Bool_t bIsTOF = kFALSE;

  Double_t dPropPID[5] = {0}; // array 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Bool_t bIsPion = kFALSE;
  Bool_t bIsKaon = kFALSE;
  Bool_t bIsProton = kFALSE;


  UInt_t iDetUsed = 0;
  AliAODTrack* track = 0x0;
  for(Int_t i(0); i < iNumTracks; i++)
  {
    track = static_cast<AliAODTrack*>(fArrTracksFiltered.At(i));

    FillPIDQA(track,0);

    bIsPion = kFALSE;
    bIsKaon = kFALSE;
    bIsProton = kFALSE;

    bIsEither = kFALSE;

    bIsTOF = (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME); // checking TOF
    iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dPropPID);

    // checking the PID probability
    if(dPropPID[2] > fCutBayesPIDPionMin) // is Pion?
      bIsPion = kTRUE;

    if(dPropPID[3] > fCutBayesPIDKaonMin) // is Kaon?
      bIsKaon= kTRUE;

    if(dPropPID[4] > fCutBayesPIDProtonMin) // is Proton?
      bIsProton = kTRUE;

    //printf("PID: pion %g (%d) / kaon %g (%d) / proton %g (%d) \n", dPropPID[2], bIsPion, dPropPID[3], bIsKaon, dPropPID[4], bIsProton);

    if(bIsTOF)
    {
      track->GetIntegratedTimes(dTOF);
      dBetaTOF = dTOF[0] / track->GetTOFsignal();
    }

    if(bIsPion)
    {
      bIsEither = kTRUE;
      new(fArrPionFiltered[iNumPions]) AliAODTrack(*track);
      iNumPions++;

      fPionsPt->Fill(track->Pt());
      fPionsEta->Fill(track->Eta());
      fPionsPhi->Fill(track->Phi());
      fPionsNsigmasTPCTOF->Fill( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion) , fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion) );
      fPionsTPCdEdx->Fill(track->P(),track->GetTPCsignal());
      fPionsBayesAsPion->Fill(dPropPID[2]);
      fPionsBayesAsKaon->Fill(dPropPID[3]);
      fPionsBayesAsProton->Fill(dPropPID[4]);
      fPionsBayes->Fill(dPropPID[2],dPropPID[3],dPropPID[4]);

      if(bIsTOF)
        fPionsTOFbeta->Fill(track->P(), dBetaTOF);
    }

    if(bIsKaon)
    {
      bIsEither = kTRUE;
      new(fArrKaonFiltered[iNumKaons]) AliAODTrack(*track);
      iNumKaons++;

      fKaonsPt->Fill(track->Pt());
      fKaonsEta->Fill(track->Eta());
      fKaonsPhi->Fill(track->Phi());
      fKaonsNsigmasTPCTOF->Fill( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon) , fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon) );
      fKaonsTPCdEdx->Fill(track->P(),track->GetTPCsignal());
      fKaonsBayesAsPion->Fill(dPropPID[2]);
      fKaonsBayesAsKaon->Fill(dPropPID[3]);
      fKaonsBayesAsProton->Fill(dPropPID[4]);
      fKaonsBayes->Fill(dPropPID[2],dPropPID[3],dPropPID[4]);

      if(bIsTOF)
        fKaonsTOFbeta->Fill(track->P(), dBetaTOF);

    }

    if(bIsProton)
    {
      bIsEither = kTRUE;
      new(fArrProtonFiltered[iNumProtons]) AliAODTrack(*track);
      iNumProtons++;

      fProtonsPt->Fill(track->Pt());
      fProtonsEta->Fill(track->Eta());
      fProtonsPhi->Fill(track->Phi());
      fProtonsNsigmasTPCTOF->Fill( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton) , fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton) );
      fProtonsTPCdEdx->Fill(track->P(),track->GetTPCsignal());
      fProtonsBayesAsPion->Fill(dPropPID[2]);
      fProtonsBayesAsKaon->Fill(dPropPID[3]);
      fProtonsBayesAsProton->Fill(dPropPID[4]);
      fProtonsBayes->Fill(dPropPID[2],dPropPID[3],dPropPID[4]);

      if(bIsTOF)
        fProtonsTOFbeta->Fill(track->P(), dBetaTOF);
    }

    if(bIsEither) // either PID as pi,K or p
    {
      FillPIDQA(track,1);

      if(!bIsTOF) // if no TOF information
        fQAPIDTOFbetaNoTOF->Fill(track->P(), dBetaTOF);
    }
  } // end of loop over filtered tracks

  fPionsMult->Fill(iNumPions);
  fKaonsMult->Fill(iNumKaons);
  fProtonsMult->Fill(iNumProtons);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FilterV0s()
{
  const Int_t iNumV0s = fAOD->GetNumberOfV0s();
  Int_t iNumK0sSelected = 0, iNumLambdaSelected = 0, iNumALambdaSelected = 0;

  if(iNumV0s < 1)
    return;

  AliAODv0* v0 = 0x0;
  for(Int_t i(0); i < iNumV0s; i++)
  {
    if(iNumK0sSelected >= 5000 || iNumLambdaSelected >= 5000 || iNumALambdaSelected >= 5000) // checking TClonesArray overflow
      return;

    v0 = static_cast<AliAODv0*>(fAOD->GetV0(i));

    if(!v0)
      continue;

    fV0candK0s = kTRUE;
    fV0candLambda = kTRUE;
    fV0candALambda = kTRUE;

    FillV0sQA(v0,0); // Filling QA histograms 'Before cuts' for V0s in selected events

    if(!IsV0Selected(v0))
      continue;

    // !!! after this point value of V0s flags (fV0cand...) should not be changed !!!

    FillV0sQA(v0,1); // Filling QA histos after cuts

    // Filling the Arrays
    if(fV0candK0s)
    {
      new(fArrV0sK0sFiltered[iNumK0sSelected]) AliAODv0(*v0);
      iNumK0sSelected++;
    }

    if(fV0candLambda)
    {
      new(fArrV0sLambdaFiltered[iNumLambdaSelected]) AliAODv0(*v0);
      iNumLambdaSelected++;
    }

    if(fV0candALambda)
    {
      new(fArrV0sALambdaFiltered[iNumALambdaSelected]) AliAODv0(*v0);
      iNumALambdaSelected++;
    }
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowPID::IsV0Selected(const AliAODv0* v0)
{
	Short_t iCounterIndex = 1;

  // invalid V0 pointer
  if(!v0)
  {
    //::Warning("IsV0Selected","Invalid pointer to V0!");
    return kFALSE;
  }

  // daughter track check
  const AliAODTrack* trackDaughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* trackDaughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  // invalid daughter track pointers
  if(!trackDaughterPos || !trackDaughterNeg)
  {
    //::Warning("IsV0Selected","Invalid pointer to V0 daughters!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // daughter tracks exists
  iCounterIndex++;

  // reconstruction method: online (on-the-fly) OR offline
  if(v0->GetOnFlyStatus() != fCutV0onFly)
  {
    //::Warning("IsV0Selected","Wrong reconstruction method!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // reconstruction method
  iCounterIndex++;

  if( fCutV0DaughterEtaMax > 0. && ( (TMath::Abs(trackDaughterNeg->Eta()) > fCutV0DaughterEtaMax) || (TMath::Abs(trackDaughterPos->Eta()) > fCutV0DaughterEtaMax) ) )
  {
    return kFALSE;
  }

  if(GetPtBinIndex(v0->Pt()) == -1) // check if the Pt is withing range of interest (histos binning)
    return kFALSE;

  if( fCutV0MotherPtMin > 0. && (v0->Pt() < fCutV0MotherPtMin) )
  {
    return kFALSE;
  }

  if( fCutV0MotherPtMax > 0. && (v0->Pt() > fCutV0MotherPtMax) )
  {
    return kFALSE;
  }

  fQAV0sCounter->Fill(iCounterIndex); // mother acceptance pT, eta
  iCounterIndex++;

  // radius of decay vertex in transverse plane
  Double_t dSecVtxCoor[3] = {0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  if( fCutV0MinDecayRadius > 0. && (dDecayRadius < fCutV0MinDecayRadius) )
  {
    //::Warning("IsV0Selected","Invalid vertex decay radius!");
    return kFALSE;
  }
  if( fCutV0MaxDecayRadius > 0. && (dDecayRadius > fCutV0MaxDecayRadius) )
  {
    //::Warning("IsV0Selected","Invalid vertex decay radius!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // decay radius
  iCounterIndex++;

	// TPC refit
  if( fCutV0refitTPC && ( !trackDaughterPos->IsOn(AliAODTrack::kTPCrefit) || !trackDaughterNeg->IsOn(AliAODTrack::kTPCrefit) ) )
  {
    //::Warning("IsV0Selected","TPC refit rejection!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // TPC refit
  iCounterIndex++;

  // Kinks
  const AliAODVertex* prodVtxDaughterPos = (AliAODVertex*) trackDaughterPos->GetProdVertex(); // production vertex of the positive daughter track
  const AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*) trackDaughterNeg->GetProdVertex(); // production vertex of the negative daughter track
  if( fCutV0rejectKinks && ( (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) || (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) ) )
  {
    //::Warning("IsV0Selected","Kink rejection!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // kinks OK
  iCounterIndex++;

  // charge of daughters
  if( trackDaughterPos->Charge() == trackDaughterNeg->Charge() ) // same charge
    return kFALSE;

  if( (trackDaughterPos->Charge() != 1) || (trackDaughterNeg->Charge() != -1) ) // expected charge
  {
    //::Warning("IsV0Selected","Bad charge!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // daughter charge
  iCounterIndex++;

  if( fCutV0DaughterPtMin > 0. && ( (trackDaughterNeg->Pt() < fCutV0DaughterPtMin) || (trackDaughterPos->Pt() < fCutV0DaughterPtMin) ) )
  {
    return kFALSE;
  }

  if( fCutV0DaughterEtaMax > 0. && ( (TMath::Abs(trackDaughterNeg->Eta()) > fCutV0DaughterEtaMax) || (TMath::Abs(trackDaughterPos->Eta()) > fCutV0DaughterEtaMax) ) )
  {
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // daughters acceptance pT, eta
  iCounterIndex++;


  // Daughters DCA to PV
  Double_t dDCAPosToPV = TMath::Abs(v0->DcaPosToPrimVertex());
  Double_t dDCANegToPV = TMath::Abs(v0->DcaNegToPrimVertex());
  if(fCutV0MinDCAtoPV > 0. && ( dDCAPosToPV < fCutV0MinDCAtoPV || dDCANegToPV < fCutV0MinDCAtoPV ) )
  {
    //::Warning("IsV0Selected","Wrong daughters DCA to PV!");
    return kFALSE;
  }
  if(fCutV0MaxDCAtoPV > 0. && ( dDCAPosToPV > fCutV0MaxDCAtoPV || dDCANegToPV > fCutV0MaxDCAtoPV ) )
  {
    //::Warning("IsV0Selected","Wrong daughters DCA to PV!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // daughters DCA to PV
  iCounterIndex++;

  // Daughter DCA among themselves
  Double_t dDCA = TMath::Abs(v0->DcaV0Daughters());
  if(fCutV0MaxDCADaughters > 0. && dDCA > fCutV0MaxDCADaughters)
  {
    //::Warning("IsV0Selected","Invalid daughters DCA!");
    return kFALSE;
  }
  fQAV0sCounter->Fill(iCounterIndex); // DCA among daughters
  iCounterIndex++;

  // re-setting the flags (to be sure)
  fV0candK0s = kTRUE;
  fV0candLambda = kTRUE;
  fV0candALambda = kTRUE;

  // is V0 either K0s or (A)Lambda candidate
  IsV0aK0s(v0);
  IsV0aLambda(v0);

  if( !fV0candK0s && !fV0candLambda && !fV0candALambda)
  {
    //::Warning("IsV0Selected","V0 is not K0s nor (A)Lambda!");
    return kFALSE; // V0 is neither K0s nor Lambda nor ALambda candidate
  }

  fQAV0sCounter->Fill(iCounterIndex); // Selected (K0s or Lambda candidates)
  iCounterIndex++;

  if(fV0candK0s)
		fQAV0sCounter->Fill(iCounterIndex); // K0s candidate

  iCounterIndex++;

  if(fV0candLambda)
    fQAV0sCounter->Fill(iCounterIndex); // Lambda candidate

  iCounterIndex++;

  if(fV0candALambda)
		fQAV0sCounter->Fill(iCounterIndex); // ALambda candidate

  iCounterIndex++;

  if(fV0candLambda && fV0candALambda)
    fQAV0sCounter->Fill(iCounterIndex); // Lambda AND ALambda candidate

  iCounterIndex++;

  if(fV0candK0s && (fV0candLambda || fV0candALambda))
    fQAV0sCounter->Fill(iCounterIndex); //  K0s AND (Lambda OR ALambda) candidates

  iCounterIndex++;

 	if(fV0candK0s && fV0candLambda && fV0candALambda)
		fQAV0sCounter->Fill(iCounterIndex); //  K0s AND Lambda AND ALambda candidates

	iCounterIndex++;

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::IsV0aK0s(const AliAODv0* v0)
{
	Short_t iCounterIndex = 0;

  if(!v0)
  {
    ::Warning("IsV0aK0s","Invalid V0 pointer!");
    fV0candK0s = kFALSE;
    return;
  }

  fQAV0sCounterK0s->Fill(iCounterIndex); // input
  iCounterIndex++;

  // inv. mass window
  Double_t dMass = v0->MassK0Short();
  if( dMass < fV0MinMassK0s || dMass > fV0MaxMassK0s )
  {
    fV0candK0s = kFALSE;
    return;
  }

	fQAV0sCounterK0s->Fill(iCounterIndex); // Inv mass window
  iCounterIndex++;


  if(fCutV0MotherRapMax > 0. && ( TMath::Abs(v0->RapK0Short()) > fCutV0MotherRapMax ) )
  {
    fV0candK0s = kFALSE;
    return;
  }

  fQAV0sCounterK0s->Fill(iCounterIndex); // V0 mother rapidity acceptance
  iCounterIndex++;

  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();

  // CPA
  Double_t dCPA = v0->CosPointingAngle(primVtx);
  if( fCutV0MinCPAK0s > 0. && dCPA < fCutV0MinCPAK0s )
  {
    fV0candK0s = kFALSE;
    return;
  }

  fQAV0sCounterK0s->Fill(iCounterIndex); // CPA
  iCounterIndex++;

  // proper life-time

  if( fCutV0NumTauK0sMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0}; // decay vector coor {xyz}
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 2; i++)
    {
      dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];
    }

    Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();

    Double_t dPropLife = ( (dMassPDGK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    if( dPropLife > (fCutV0NumTauK0sMax * 2.68) )
    {
      fV0candK0s = kFALSE;
      return;
    }
  }

  fQAV0sCounterK0s->Fill(iCounterIndex); // Proper lifetime
  iCounterIndex++;


  // Armenteros-Podolaski plot
  if(fCutV0K0sArmenterosAlphaMin > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if( dPtArm < (fCutV0K0sArmenterosAlphaMin * TMath::Abs(dAlpha)) )
    {
      fV0candK0s = kFALSE;
      return;
    }
  }

  fQAV0sCounterK0s->Fill(iCounterIndex); // Armernteros Podolanski
  iCounterIndex++;


  // cross-contamination

  fQAV0sCounterK0s->Fill(iCounterIndex); // Selected
  iCounterIndex++;

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::IsV0aLambda(const AliAODv0* v0)
{
	Short_t iCounterIndex = 0;

  fQAV0sCounterLambda->Fill(iCounterIndex); // input
  iCounterIndex++;

  if(!v0)
  {
    fV0candLambda = kFALSE;
    fV0candALambda = kFALSE;
    return;
  }

  // inv. mass window
  Double_t dMassLambda = v0->MassLambda();
  Double_t dMassALambda = v0->MassAntiLambda();

  if(dMassLambda < fV0MinMassLambda || dMassLambda > fV0MaxMassLambda)
    fV0candLambda = kFALSE;

  if(dMassALambda < fV0MinMassLambda || dMassALambda > fV0MaxMassLambda)
    fV0candALambda = kFALSE;

  if(!fV0candLambda && !fV0candALambda)
    return;

  fQAV0sCounterLambda->Fill(iCounterIndex); // Inv. mass window
  iCounterIndex++;

  if(fCutV0MotherRapMax > 0. && ( TMath::Abs(v0->RapLambda()) > fCutV0MotherRapMax ) )
  {
    fV0candLambda = kFALSE;
    fV0candALambda = kFALSE;
    return;
  }

	fQAV0sCounterLambda->Fill(iCounterIndex); // V0 mother rapidity window
  iCounterIndex++;

  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();

  // CPA
  Double_t dCPA = v0->CosPointingAngle(primVtx);
  if( fCutV0MinCPALambda > 0. && dCPA < fCutV0MinCPALambda )
  {
    fV0candLambda = kFALSE;
    fV0candALambda = kFALSE;
    return;
  }

	fQAV0sCounterLambda->Fill(iCounterIndex); // CPA
  iCounterIndex++;

  // proper life-time
  if( fCutV0NumTauLambdaMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0}; // decay vector coor {xyz}
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 2; i++)
    {
      dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];
    }

    Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

    Double_t dPropLife = ( (dMassPDGLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    if( dPropLife > (fCutV0NumTauLambdaMax * 7.89) )
    {
      fV0candLambda = kFALSE;
      fV0candALambda = kFALSE;
      return;
    }
  }

  fQAV0sCounterLambda->Fill(iCounterIndex); // proper life-time
  iCounterIndex++;

 	// proton PID of Lambda Candidates
  if( (fCutV0ProtonNumSigmaMax > 0.) && (fCutV0ProtonPIDPtMax > 0.) && fPIDResponse )
  {
  	const AliAODTrack* trackDaughterPos = (AliAODTrack*) v0->GetDaughter(0);
  	const AliAODTrack* trackDaughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  	Double_t dSigmaProtPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kProton));
  	Double_t dSigmaProtNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kProton));

    if(fV0candLambda && (trackDaughterPos->Pt() < fCutV0ProtonPIDPtMax) && (dSigmaProtPos > fCutV0ProtonNumSigmaMax)) // positive daughter is not a proton (within TPC sigmas)
    {
      fV0candLambda = kFALSE;
    }

    if(fV0candALambda && (trackDaughterNeg->Pt() < fCutV0ProtonPIDPtMax) && (dSigmaProtNeg > fCutV0ProtonNumSigmaMax)) // negative daughter is not a proton (within TPC sigmas)
    {
      fV0candALambda = kFALSE;
    }

    if(!fV0candLambda && !fV0candALambda)
    {
      return;
    }
  }

	fQAV0sCounterLambda->Fill(iCounterIndex); // proton pid
  iCounterIndex++;

  // cross-contamination

	fQAV0sCounterLambda->Fill(iCounterIndex); // selected
  iCounterIndex++;

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FillEventQA(const AliAODEvent* event, const Short_t iQAindex)
{
  // event QA procedure
	const AliAODVertex* aodVtx = event->GetPrimaryVertex();
  const Double_t dVtxZ = aodVtx->GetZ();
  const Int_t iNumContr = aodVtx->GetNContributors();
  const AliAODVertex* spdVtx = event->GetPrimaryVertexSPD();
  const Int_t iNumContrSPD = spdVtx->GetNContributors();
  const Double_t spdVtxZ = spdVtx->GetZ();

  fQAEventsPVz[iQAindex]->Fill(dVtxZ);
  fQAEventsNumContrPV[iQAindex]->Fill(iNumContr);
  fQAEventsNumSPDContrPV[iQAindex]->Fill(iNumContrSPD);
  fQAEventsDistPVSPD[iQAindex]->Fill(TMath::Abs(dVtxZ - spdVtxZ));

  if(fPP || fPPb) // pp OR pPb analysis
  {
    // event / physics selection criteria
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
    UInt_t fSelectMask = inputHandler->IsEventSelected();

    if( fSelectMask& AliVEvent::kINT7 ) { fQAEventsTriggerSelection[iQAindex]->Fill("kINT7",1); }
    else if (fSelectMask& AliVEvent::kHighMultV0) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultV0",1); }
    else if (fSelectMask& AliVEvent::kHighMultSPD) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultSPD",1); }
    else { fQAEventsTriggerSelection[iQAindex]->Fill("Other",1); }

    // SPD vertexer resolution
    Double_t cov[6] = {0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    fQAEventsSPDresol[iQAindex]->Fill(zRes);
  }

	return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FillTrackQA(const AliAODTrack* track, const Short_t iQAindex)
{
  fQATracksPt[iQAindex]->Fill(track->Pt());
  fQATracksPhi[iQAindex]->Fill(track->Phi());
  fQATracksEta[iQAindex]->Fill(track->Eta());
  fQATracksFilterMap[iQAindex]->Fill(track->GetFilterMap());
  fQATracksNumTPCcls[iQAindex]->Fill(track->GetTPCNcls());

  for(Short_t i(0); i < 12; i++)
  {
    if(track->TestFilterBit(TMath::Power(2.,i)))
      fQATracksFilterMapBit[iQAindex]->Fill(i);
  }

  const Float_t dDCAx = track->XAtDCA();
  const Float_t dDCAy = track->YAtDCA();
  const Float_t dDCAz = track->ZAtDCA();

  fQATracksDCAxy[iQAindex]->Fill(TMath::Sqrt(dDCAy*dDCAy + dDCAx*dDCAx));
  fQATracksDCAz[iQAindex]->Fill(dDCAz);

  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);

  fQATracksTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
  fQATracksTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );

  if(pidStatusTPC == AliPIDResponse::kDetPidOk)
  {
    fQATracksTPCdEdx[iQAindex]->Fill(track->P(), track->GetTPCsignal());
    //printf("signal TPC (in OK) %g\n",track->GetTPCsignal());
  }

  Bool_t bIsTOF = (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME); // checking TOF

  if(bIsTOF)
  {
    fQATracksTOF[iQAindex]->Fill(track->P(), track->GetTOFsignal()/1000);
    //printf("signal TOF (in OK) %g\n",track->GetTOFsignal());
    Double_t dTOF[5];
    track->GetIntegratedTimes(dTOF);
    Double_t dBetaTOF = dTOF[0] / track->GetTOFsignal();
    fQATracksTOFbeta[iQAindex]->Fill(track->P(),dBetaTOF);
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FillPIDQA(const AliAODTrack* track, const Short_t iQAindex)
{
  Double_t dTOF[5] = {-990};

  Double_t dBetaTOF = 0.001;

  Bool_t bIsTOF = (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME); // checking TOF

  if(bIsTOF)
  {
    track->GetIntegratedTimes(dTOF);
    dBetaTOF = dTOF[0] / track->GetTOFsignal();
    fQAPIDTOFbetaWithTOF[iQAindex]->Fill(track->P(), dBetaTOF);

    fQAPIDNsigmasTOFasPion[iQAindex]->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
    fQAPIDNsigmasTOFasKaon[iQAindex]->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
    fQAPIDNsigmasTOFasProton[iQAindex]->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
  }

  fQAPIDTOFbeta[iQAindex]->Fill(track->P(), dBetaTOF);
  fQAPIDTPCdEdx[iQAindex]->Fill(track->P(), track->GetTPCsignal());

  fQAPIDNsigmasTPCasPion[iQAindex]->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
  fQAPIDNsigmasTPCasKaon[iQAindex]->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
  fQAPIDNsigmasTPCasProton[iQAindex]->Fill(track->P(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));

  // Fill bayes probability
  Double_t dPropPID[5] = {0};
  UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dPropPID);
  fQAPIDBayesProbPion[iQAindex]->Fill(dPropPID[2]);
  fQAPIDBayesProbKaon[iQAindex]->Fill(dPropPID[3]);
  fQAPIDBayesProbProton[iQAindex]->Fill(dPropPID[4]);
  fQAPIDBayesProb[iQAindex]->Fill(dPropPID[2],dPropPID[3],dPropPID[4]);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::FillV0sQA(const AliAODv0* v0, const Short_t iQAindex)
{
	const Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
	const Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
	AliAODVertex* primVtx = fAOD->GetPrimaryVertex();

	Double_t dPrimVtxCoor[3] = {0};
	primVtx->GetXYZ(dPrimVtxCoor);

	AliAODTrack* trackDaughter[2] = {0x0, 0x0};
	AliAODVertex* prodVtxDaughter[2] = {0x0, 0x0};
	Double_t dSecVtxCoor[3] = {0};
	Double_t dDecayRadius = 0;
	Double_t dDecayCoor[3] = {0};
	Double_t dPropLife = 0;

	// universal cuts

	for(Int_t i(0); i < 2; i++) // daughters loop
	{
		trackDaughter[i] = (AliAODTrack*) v0->GetDaughter(i);
		if(!trackDaughter[i])
			continue;

		fQAV0sTPCRefit[iQAindex]->Fill(trackDaughter[i]->IsOn(AliAODTrack::kTPCrefit));
    fQAV0sDaughterPt[iQAindex]->Fill(trackDaughter[i]->Pt());
		fQAV0sDaughterPhi[iQAindex]->Fill(trackDaughter[i]->Phi());
		fQAV0sDaughterEta[iQAindex]->Fill(trackDaughter[i]->Eta());

		// kinks
		prodVtxDaughter[i] = (AliAODVertex*) trackDaughter[i]->GetProdVertex();
		fQAV0sKinks[iQAindex]->Fill( (prodVtxDaughter[i]->GetType() == AliAODVertex::kKink) );

    // TPC dEdx
    fQAV0sDaughterTPCdEdxPt[iQAindex]->Fill(trackDaughter[i]->Pt(),fTPCPIDResponse.GetTrackdEdx(trackDaughter[i]));
	}

	fQAV0sRecoMethod[iQAindex]->Fill(v0->GetOnFlyStatus());
  fQAV0sMotherPt[iQAindex]->Fill(v0->Pt());
	fQAV0sMotherPhi[iQAindex]->Fill(v0->Phi());
	fQAV0sMotherEta[iQAindex]->Fill(v0->Eta());
	fQAV0sDCADaughters[iQAindex]->Fill(TMath::Abs(v0->DcaV0Daughters()));
	fQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaPosToPrimVertex());
	fQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaNegToPrimVertex());

	//Decay radius
	v0->GetSecondaryVtx(dSecVtxCoor);
	dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
	fQAV0sDecayRadius[iQAindex]->Fill(dDecayRadius);


	for(Int_t i(0); i < 2; i++)
  {
    dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];
  }

	// Particle dependent cuts (K0s/(A)Lambda)
	if(fV0candK0s)
	{
		fQAV0sMotherRapK0s[iQAindex]->Fill(v0->RapK0Short());
		fQAV0sInvMassK0s[iQAindex]->Fill(v0->MassK0Short());
    fQAV0sCPAK0s[iQAindex]->Fill(v0->CosPointingAngle(primVtx));
		dPropLife = ( (dMassPDGK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
		fQAV0sNumTauK0s[iQAindex]->Fill(dPropLife);
		fQAV0sArmenterosK0s[iQAindex]->Fill(v0->AlphaV0(),v0->PtArmV0());
	}

	if(fV0candLambda || fV0candALambda)
	{
    if(fV0candLambda)
    {
      fQAV0sInvMassLambda[iQAindex]->Fill(v0->MassLambda());
      if(fPIDResponse)
      {
        fQAV0sProtonNumSigmaPtLambda[iQAindex]->Fill(trackDaughter[0]->Pt(),TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kProton)));
      }
    }

    if(fV0candALambda)
    {
      fQAV0sInvMassLambda[iQAindex]->Fill(v0->MassAntiLambda());
      // PID number of sigmas proton
      if(fPIDResponse)
      {
        fQAV0sProtonNumSigmaPtLambda[iQAindex]->Fill(trackDaughter[1]->Pt(),TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kProton)));
      }
    }

    fQAV0sMotherRapLambda[iQAindex]->Fill(v0->RapLambda());
		fQAV0sCPALambda[iQAindex]->Fill(v0->CosPointingAngle(primVtx));
		dPropLife = ( (dMassPDGLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
		fQAV0sNumTauLambda[iQAindex]->Fill(dPropLife);
		fQAV0sArmenterosLambda[iQAindex]->Fill(v0->AlphaV0(),v0->PtArmV0());
	}
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
      if(!MultSelection->IsEventSelected())
      {
        fCentPercentile = -100.;
        fCentBinIndex = -1;
        return;
      }

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

  if(centrCode > fNumCentBins)
  {
    fCentPercentile = -100.;
    fCentBinIndex = -1;
    return;
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
//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowPID::GetCentBinIndexPP(const Int_t iMult)
{
  if(!fPP)
  {
    ::Error("GetCentBinIndexPP","Not a pp analysis! Returning -1!");
    return -1;
  }
  for(Int_t i(0); i < fNumCentBins; i++)
  {
    if( (iMult >= fCentBinEdges[i]) && (iMult < fCentBinEdges[i+1]) )
      return i;
  }

  return -1;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowPID::GetMinvFlowBinIndexK0s(const Double_t dMass)
{
  for(Int_t i(0); i < fNumMinvFlowBinsK0s; i++)
  {
    if( (dMass >= fMinvFlowBinEdgesK0s[i]) && (dMass < fMinvFlowBinEdgesK0s[i+1]) )
      return i;
  }

  return -1;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowPID::GetMinvFlowBinIndexLambda(const Double_t dMass)
{
  for(Int_t i(0); i < fNumMinvFlowBinsLambda; i++)
  {
    if( (dMass >= fMinvFlowBinEdgesLambda[i]) && (dMass < fMinvFlowBinEdgesLambda[i+1]) )
      return i;
  }

  return -1;
}
//_____________________________________________________________________________
//_____________________________________________________________________________

//____________________________________________________________________
// Below is implementation of GF by Katarina wrapped in DoGenFramKatarina method
//_____________________________________________________________________
void AliAnalysisTaskFlowPID::DoGenFramKatarina()
{
  Double_t dEtaGap = 0;
  Double_t dMass = 0;

  for(Int_t iEtaGap(0); iEtaGap < fNumEtaGap; iEtaGap++)
  {
    //printf("Ref flow \n");
    dEtaGap = fEtaGap[iEtaGap];

    GFKFillRefVectors(fArrTracksFiltered,dEtaGap);

    for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
    {
      GFKDoRefFlow(fcn2Tracks[iEtaGap][iHarm][fSampleBinIndex],fcn4Tracks[iEtaGap][iHarm][fSampleBinIndex],fHarmonics[iHarm],dEtaGap);
    }

    //printf("Diff flow \n");
    // Diff flow
    if(fDiffFlow)
    {
      for(Int_t iPt(0); iPt < fNumPtBins; iPt++)
      {
        // fill p,q vectors

        GFKFillVectors(fArrTracksFiltered,iEtaGap,iPt,0);
        for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
        {
          GFKDoDiffFlow(fdn2Tracks[0][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn2Tracks[1][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn4Tracks[iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fHarmonics[iHarm],dEtaGap,iPt);
        }

        if(fPID)
        {
          // pions
          GFKFillVectors(fArrPionFiltered,iEtaGap,iPt,0);
          for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
          {
            GFKDoDiffFlow(fdn2Pion[0][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn2Pion[1][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn4Pion[iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fHarmonics[iHarm],dEtaGap,iPt);
          }

          // kaons
          GFKFillVectors(fArrKaonFiltered,iEtaGap,iPt,0);
          for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
          {
            GFKDoDiffFlow(fdn2Kaon[0][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn2Kaon[1][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn4Kaon[iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fHarmonics[iHarm],dEtaGap,iPt);
          }

          // protons
          GFKFillVectors(fArrProtonFiltered,iEtaGap,iPt,0);
          for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
          {
            GFKDoDiffFlow(fdn2Proton[0][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn2Proton[1][iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fdn4Proton[iEtaGap][iHarm][fCentBinIndex][fSampleBinIndex],fHarmonics[iHarm],dEtaGap,iPt);
          }


          // K0s
          for(Int_t iMass(0); iMass < fNumMinvFlowBinsK0s; iMass++)
          {
            dMass = (fMinvFlowBinEdgesK0s[iMass+1]+fMinvFlowBinEdgesK0s[iMass])/2;
            //printf("fill Kos pre Gap %g pT %d \n",dEtaGap,iPt);
            GFKFillVectors(fArrV0sK0sFiltered,iEtaGap,iPt,1,iMass);
            //printf("fill Kos post\n");

            for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
            {
              GFKDoDiffFlowV0s(fdn2K0s[0][iEtaGap][iHarm][fCentBinIndex],fdn2K0s[1][iEtaGap][iHarm][fCentBinIndex],fdn4K0s[iEtaGap][iHarm][fCentBinIndex],fHarmonics[iHarm],dEtaGap,iPt,dMass);
            }
          }

          //Lambda
          for(Int_t iMass(0); iMass < fNumMinvFlowBinsLambda; iMass++)
          {
            dMass = (fMinvFlowBinEdgesLambda[iMass+1]+fMinvFlowBinEdgesLambda[iMass])/2;

            // filling Lambda
            GFKFillVectors(fArrV0sLambdaFiltered,iEtaGap,iPt,2,iMass);
            for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
            {
              GFKDoDiffFlowV0s(fdn2Lambda[0][iEtaGap][iHarm][fCentBinIndex],fdn2Lambda[1][iEtaGap][iHarm][fCentBinIndex],fdn4Lambda[iEtaGap][iHarm][fCentBinIndex],fHarmonics[iHarm],dEtaGap,iPt,dMass);
            }

            //Filling antilambdas
            GFKFillVectors(fArrV0sALambdaFiltered,iEtaGap,iPt,3,iMass);
            for(Int_t iHarm = 0; iHarm < fNumHarmonics; iHarm++)
            {
              GFKDoDiffFlowV0s(fdn2Lambda[0][iEtaGap][iHarm][fCentBinIndex],fdn2Lambda[1][iEtaGap][iHarm][fCentBinIndex],fdn4Lambda[iEtaGap][iHarm][fCentBinIndex],fHarmonics[iHarm],dEtaGap,iPt,dMass);
            }
          }
        } // end of DoPID flag

      } // end of pT bins loop

    } // end of DoDiffFlow flag

  } // end of loop over eta gaps

  return;
}
//_____________________________________________________________________
void AliAnalysisTaskFlowPID::GFKFillRefVectors(TClonesArray &array,const Double_t dEtaGap)
{

  //printf("GFKFillRefVectors \n");

  double weight = 1.;
  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for RFPs

  Double_t Qcos[fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t Qsin[fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t QcosGap[fGFKNumVectors][fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t QsinGap[fGFKNumVectors][fMaxNumHarmonics][fMaxNumWeights] = {0};

  AliAODTrack* aodTrk = 0x0;
  Double_t dEta = 0;
  //Refs
  for(Int_t i(0); i < array.GetEntriesFast(); i++)
  {
    aodTrk = static_cast<AliAODTrack*>(array.At(i));
    dEta = aodTrk->Eta();

    for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
    {
      for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
      {
        if(dEtaGap < 0.) // with no gap -> only positive quantities are filled
        {
          Qcos[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          Qsin[iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
        else if(!(dEtaGap < 0.)) // with eta gap
        {
          if(dEta > dEtaGapCut)
          {
            QcosGap[0][iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap[0][iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }

          if(dEta < -dEtaGapCut)
          {
            QcosGap[1][iharm][ipow] += TMath::Power(weight, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap[1][iharm][ipow] += TMath::Power(weight, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
    }
  } // end of loop over particles

  for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
  {
    for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
    {
      Qvector[iharm][ipow] = TComplex(Qcos[iharm][ipow], Qsin[iharm][ipow]);

      for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
      {
        QvectorGap[iVec][iharm][ipow] = TComplex(QcosGap[iVec][iharm][ipow], QsinGap[iVec][iharm][ipow]);
        //printf("Qvector[iVec %d][iharm %d][ipow %d]: Re %g / Im %g\n",iVec, iharm, ipow, QvectorGap[iVec][iharm][ipow].Re(), QvectorGap[iVec][iharm][ipow].Im() );
      }
    }
  }


  return;
}
//_____________________________________________________________________
void AliAnalysisTaskFlowPID::GFKFillVectors(TClonesArray &array, const Int_t iEtaGapIndex, const Int_t iPtBin, const Int_t iV0s, const Int_t iMass)
{

  Double_t weightPt = 1.;
  const Double_t dEtaGap = fEtaGap[iEtaGapIndex];
  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for RFPs

  const Int_t iNumTracks = array.GetEntriesFast();

  Double_t pcosPt[fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t psinPt[fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t qcosPt[fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t qsinPt[fMaxNumHarmonics][fMaxNumWeights] = {0};

  Double_t pcosPtGap[fGFKNumVectors][fMaxNumHarmonics][fMaxNumWeights] = {0};
  Double_t psinPtGap[fGFKNumVectors][fMaxNumHarmonics][fMaxNumWeights] = {0};
  //Double_t qcosPtGap[fGFKNumVectors][fMaxNumHarmonics][fMaxNumWeights][fNumPtBins] = {0};
  //Double_t qsinPtGap[fGFKNumVectors][fMaxNumHarmonics][fMaxNumWeights][fNumPtBins] = {0};

  //AliAODTrack* aodTrk = 0x0;
  AliVTrack* aodTrk = 0x0;
  AliAODv0* v0 = 0x0; // only for extracting the mass
  Double_t dEta = 0;
  Double_t dMass = 0;

  Int_t iPtBinIndex = -1;
  Int_t iMassBinIndex = -1;

  //POIs
  for(Int_t i(0); i < iNumTracks; i++)
  {
    if(iV0s == 0) // tracks
    {
      aodTrk = static_cast<AliAODTrack*>(array.At(i));
    }
    else
    {
      aodTrk = static_cast<AliAODv0*>(array.At(i));
      v0 = static_cast<AliAODv0*>(array.At(i));

      if(iV0s == 1) // Kos
      {
        dMass = v0->MassK0Short();
        iMassBinIndex = GetMinvFlowBinIndexK0s(dMass);
      }
      else if(iV0s == 2) // Lambda
      {
        dMass = v0->MassLambda();
        iMassBinIndex = GetMinvFlowBinIndexLambda(dMass);
      }
      else // anti - lambda
      {
        dMass = v0->MassAntiLambda();
        iMassBinIndex = GetMinvFlowBinIndexLambda(dMass);
      }
    }

    dEta = aodTrk->Eta();

    iPtBinIndex = GetPtBinIndex(aodTrk->Pt());
    if(iPtBinIndex != iPtBin)
      continue;

    if(iV0s && iMass != iMassBinIndex)
      continue;


    if(dEtaGap < 0.) // with no gap -> only positive quantities are filled
    {
      // filling V0s Inv Mass distributions

      if(iV0s == 1) //K0s
        fInvMassPtK0s[0][iEtaGapIndex][fCentBinIndex]->Fill(aodTrk->Pt(),dMass);

      if(iV0s > 1) //Lambda & ALambda
        fInvMassPtLambda[0][iEtaGapIndex][fCentBinIndex]->Fill(aodTrk->Pt(),dMass);

      for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
      {
        for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
        {
          pcosPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          psinPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          qcosPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          qsinPt[iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }
    else if(!(dEtaGap < 0.)) // with eta gap
    {
      if(dEta > dEtaGapCut)
      {
        // filling V0s Inv Mass distributions
        if(iV0s == 1) //K0s // **** HERE IS THE BASTARD
        {
          fInvMassPtK0s[0][iEtaGapIndex][fCentBinIndex]->Fill(aodTrk->Pt(),dMass);
         // printf("KosPtMass:: Fill pt %g & mass %g in (eta: %d / cent: %d)\n",aodTrk->Pt(),dMass,iEtaGapIndex,fCentBinIndex);
        }
        if(iV0s > 1) //Lambda & ALambda
          fInvMassPtLambda[0][iEtaGapIndex][fCentBinIndex]->Fill(aodTrk->Pt(),dMass);

        for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
        {
          for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
          {
              pcosPtGap[0][iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
              psinPtGap[0][iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
              //qcosPtGapPos[iharm][ipow][iPtBinIndex] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
              //qsinPtGapPos[iharm][ipow][iPtBinIndex] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }

      if(dEta < -dEtaGapCut)
      {
        // filling V0s Inv Mass distributions
        if(iV0s == 1) //K0s
          fInvMassPtK0s[1][iEtaGapIndex][fCentBinIndex]->Fill(aodTrk->Pt(),dMass);
        if(iV0s > 1) //Lambda & ALambda
          fInvMassPtLambda[1][iEtaGapIndex][fCentBinIndex]->Fill(aodTrk->Pt(),dMass);
        for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
        {
          for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
          {
            pcosPtGap[1][iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            psinPtGap[1][iharm][ipow] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
            //qcosPtGapNeg[iharm][ipow][iPtBinIndex] += TMath::Power(weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            //qsinPtGapNeg[iharm][ipow][iPtBinIndex] += TMath::Power(weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
    }
  }

  for(Int_t iharm = 0; iharm < fMaxNumHarmonics; iharm++)
  {
    for(Int_t ipow = 0; ipow < fMaxNumWeights; ipow++)
    {
      pvector[iharm][ipow] = TComplex(pcosPt[iharm][ipow], psinPt[iharm][ipow]);
      qvector[iharm][ipow] = TComplex(qcosPt[iharm][ipow], qsinPt[iharm][ipow]);

      for(Int_t iVec(0); iVec < fGFKNumVectors; iVec++)
      {
        pvectorGap[iVec][iharm][ipow] = TComplex(pcosPtGap[iVec][iharm][ipow], psinPtGap[iVec][iharm][ipow]);
        //printf("pvector[iVec %d][iharm %d][ipow %d]: Re %g / Im %g\n",iVec, iharm, ipow, pvectorGap[iVec][iharm][ipow][iPt].Re(), pvectorGap[iVec][iharm][ipow][iPt].Im() );
        //qvectorGap[iVec][iharm][ipow][iPt] = TComplex(qcosPtGap[iVec][iharm][ipow][iPt], qsinPtGap[iVec][iharm][ipow][iPt]);
      }
    }
  }

  return;
}
//_____________________________________________________________________
void AliAnalysisTaskFlowPID::GFKDoRefFlow(TProfile* prof2,TProfile* prof4, const Short_t iHarm, const Double_t dEtaGap)
{
  Double_t dValue = 0;
  TComplex* vector = 0x0;

  if(dEtaGap < 0.) // no gap
  {
    Double_t Dn2 = Two(0,0)->Re();
    Double_t Dn4 = Four(0,0,0,0)->Re();

    //printf("DDn2: %g / DDn4: %g\n",DDn2, DDn4);

    if(Dn2 != 0)
    {
      //..v2{2} = <cos2(phi1 - phi2)>
      vector = Two(iHarm, -iHarm);
      dValue = vector->Re()/Dn2;
      if( TMath::Abs(dValue < 1) )
        prof2->Fill(fCentPercentile, dValue, Dn2);
    }

    if(Dn4 != 0)
    {
      //..v2{4} = <cos2(phi1 + phi2 - phi3 - phi4)>
      vector = Four(iHarm, iHarm, -iHarm, -iHarm);
      dValue = vector->Re()/Dn4;
      if( TMath::Abs(dValue < 1) )
        prof4->Fill(fCentPercentile, dValue, Dn4);
    }
  }
  else // with gap
  {
    Double_t Dn2 = TwoGap(0, 0)->Re();
    Double_t Dn4 = 0; // not implemented
    //double Dn4 = FourGap(0,0,0,0)->Re();

    //printf("DDn2: %g / DDn4: %g\n",Dn2, Dn4);

    if(Dn2 != 0)
    {
      //..v2{2} = <cos2(phi1 - phi2)>
      vector = TwoGap(iHarm, -iHarm);
      dValue = vector->Re()/Dn2;
      if( TMath::Abs(dValue < 1) )
        prof2->Fill(fCentPercentile, dValue, Dn2);
    }

    if(Dn4 != 0)
    {
      //..v2{4} = <cos2(phi1 + phi2 - phi3 - phi4)>
      vector = Four(iHarm, iHarm, -iHarm, -iHarm);
      dValue = vector->Re()/Dn4;
      if( TMath::Abs(dValue < 1) )
        prof4->Fill(fCentPercentile, dValue, Dn4);
    }
  }

  return;
}
//_____________________________________________________________________
void AliAnalysisTaskFlowPID::GFKDoDiffFlow(TProfile* prof2Pos, TProfile* prof2Neg, TProfile* prof4, const Short_t iHarm, const Double_t dEtaGap, const Int_t iPtBin)
{
  Double_t dValue = 0;
  TComplex* vector = 0x0;

  if(dEtaGap < 0.) // no gap
  {
    Double_t DDn2 = TwoDiff(0,0)->Re();
    Double_t DDn4 = FourDiff(0,0,0,0)->Re();

    //printf("DDn2: %g / DDn4: %g\n",DDn2, DDn4);

    //if(DDn2 > 1)
    if(DDn2 != 0)
    {
      //..v2{2}
      vector = TwoDiff(iHarm, -iHarm);
      dValue = vector->Re()/DDn2;
      if( TMath::Abs(dValue < 1) )
        prof2Pos->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dValue, DDn2);

    }

    if(DDn4 != 0)
    {
      //..v2{4}
      vector = FourDiff(iHarm, iHarm, -iHarm, -iHarm);
      dValue = vector->Re()/DDn4;
      if( TMath::Abs(dValue < 1) )
        prof4->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dValue, DDn4);
    }
  }
  else
  { // with gap
    Double_t DDn2Pos = TwoDiffGapPos(0,0)->Re();
    Double_t DDn2Neg = TwoDiffGapNeg(0,0)->Re();
    Double_t DDn4 = 0;

    //printf("DDn2Pos: %g / DDn2Neg: %g\n",DDn2Pos, DDn2Neg);

    //if(DDn2 > 1)
    if(DDn2Pos != 0)
    {
      //..v2{2}
      vector = TwoDiffGapPos(iHarm, -iHarm);
      dValue = vector->Re()/DDn2Pos;
      if( TMath::Abs(dValue < 1))
        prof2Pos->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dValue, DDn2Pos);

    }

    if(DDn2Neg != 0)
    {
      //..v2{2}
      vector = TwoDiffGapNeg(iHarm, -iHarm);
      dValue = vector->Re()/DDn2Neg;
      if( TMath::Abs(dValue < 1))
        prof2Neg->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dValue, DDn2Neg);

    }

    if(DDn4 != 0)
    {
      //..v2{4}
      vector = FourDiff(iHarm, iHarm, -iHarm, -iHarm);
      dValue = vector->Re()/DDn4;
      if( TMath::Abs(dValue < 1))
        prof4->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dValue, DDn4);
    }

  }


}
//______________________________________________
void AliAnalysisTaskFlowPID::GFKDoDiffFlowV0s(TProfile2D* prof2Pos, TProfile2D* prof2Neg, TProfile2D* prof4, const Short_t iHarm, const Double_t dEtaGap, const Int_t iPtBin, const Double_t dMass)
{
  Double_t dValue = 0;
  TComplex* vector = 0x0;

  if(dEtaGap < 0.) // no gap
  {
    Double_t DDn2 = TwoDiff(0,0)->Re();
    Double_t DDn4 = FourDiff(0,0,0,0)->Re();

    //printf("DDn2: %g / DDn4: %g\n",DDn2, DDn4);

    //if(DDn2 > 1)
    if(DDn2 != 0)
    {
      //..v2{2}
      vector = TwoDiff(iHarm, -iHarm);
      dValue = vector->Re()/DDn2;
      //printf("Filling:[%g][%g]: %g, %g\n", (fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2,dMass, dValue, DDn2);
      if( TMath::Abs(dValue < 1) )
        prof2Pos->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dMass, dValue, DDn2);


    }

    if(DDn4 != 0)
    {
      //..v2{4}
      vector = FourDiff(iHarm, iHarm, -iHarm, -iHarm);
      dValue = vector->Re()/DDn4;
      //printf("Filling:[%g][%g]: %g, %g\n", (fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2,dMass, dValue, DDn4);
      if( TMath::Abs(dValue < 1) )
        prof4->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dMass, dValue, DDn4);

    }
  }
  else
  { // with gap
    Double_t DDn2Pos = TwoDiffGapPos(0,0)->Re();
    Double_t DDn2Neg = TwoDiffGapNeg(0,0)->Re();
    Double_t DDn4 = 0;

    //printf("DDn2Pos: %g / DDn2Neg: %g\n",DDn2Pos, DDn2Neg);

    //if(DDn2 > 1)
    if(DDn2Pos != 0)
    {
      //..v2{2}
      vector = TwoDiffGapPos(iHarm, -iHarm);
      dValue = vector->Re()/DDn2Pos;
      //printf("Filling:[%g][%g]: %g, %g\n", (fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2,dMass, dValue, DDn2Pos);
      if( TMath::Abs(dValue < 1))
        prof2Pos->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dMass, dValue, DDn2Pos);

    }

    if(DDn2Neg != 0)
    {
      //..v2{2}
      vector = TwoDiffGapNeg(iHarm, -iHarm);
      dValue = vector->Re()/DDn2Neg;
      //printf("Filling:[%g][%g]: %g, %g\n", (fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2,dMass, dValue, DDn2Neg);
      if( TMath::Abs(dValue < 1))
        prof2Neg->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dMass, dValue, DDn2Neg);

    }

    if(DDn4 != 0)
    {
      //..v2{4}
      vector = FourDiff(iHarm, iHarm, -iHarm, -iHarm);
      dValue = vector->Re()/DDn4;
      //printf("Filling:[%g][%g]: %g, %g\n", (fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2,dMass, dValue, DDn4);
      if( TMath::Abs(dValue < 1))
        prof4->Fill((fPtBinEdges[iPtBin+1] + fPtBinEdges[iPtBin])/2, dMass, dValue, DDn4);
    }

  }


}
//_____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::Q(int n, int p)
{

  if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::QGapPos(int n, int p)
{

  if(n>=0) return QvectorGap[0][n][p];
  else return TComplex::Conjugate(QvectorGap[0][-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::QGapNeg(int n, int p)
{

  if(n>=0) return QvectorGap[1][n][p];
  else return TComplex::Conjugate(QvectorGap[1][-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::p(int n, int p)
{

  if(n>=0) return pvector[n][p];
  else return TComplex::Conjugate(pvector[n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::pGapPos(int n, int p)
{

  if(n>=0) return pvectorGap[0][n][p];
  else return TComplex::Conjugate(pvectorGap[0][n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::pGapNeg(int n, int p)
{

  if(n>=0) return pvectorGap[1][n][p];
  else return TComplex::Conjugate(pvectorGap[1][n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowPID::q(int n, int p)
{

  if(n>=0) return qvector[n][p];
  else return TComplex::Conjugate(qvector[n][p]);

}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::Two(int n1, int n2)
{
  TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  TComplex *out = (TComplex*) &formula;
  return out;
}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::TwoGap(int n1, int n2)
{
  TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
  TComplex *out = (TComplex*) &formula;
  return out;
}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::TwoDiff(int n1, int n2)
{

  TComplex formula = p(n1,1)*Q(n2,1) - q(n1+n2,1);
  TComplex *out = (TComplex*) &formula;
  return out;
}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::TwoDiffGapPos(int n1, int n2)
{

  TComplex formula = pGapPos(n1,1)*QGapNeg(n2,1);
  TComplex *out = (TComplex*) &formula;
  return out;
}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::TwoDiffGapNeg(int n1, int n2)
{

  TComplex formula = pGapNeg(n1,1)*QGapPos(n2,1);
  TComplex *out = (TComplex*) &formula;
  return out;
}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::Four(int n1, int n2, int n3, int n4)
{

  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
  TComplex *out = (TComplex*) &formula;
  return out;

}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::FourGap(int n1, int n2, int n3, int n4)
{

  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
  TComplex *out = (TComplex*) &formula;
  return out;

}
//____________________________________________________________________
TComplex* AliAnalysisTaskFlowPID::FourDiff(int n1, int n2, int n3, int n4)
{


  TComplex formula = p(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*q(n1+n3,2)*Q(n4,1)
                    - p(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*q(n1+n4,2)
                    + Q(n2+n3,2)*q(n1+n4,2)-p(n1,1)*Q(n3,1)*Q(n2+n4,2)+q(n1+n3,2)*Q(n2+n4,2)
                    + 2.*Q(n3,1)*q(n1+n2+n4,3)-p(n1,1)*Q(n2,1)*Q(n3+n4,2)+q(n1+n2,2)*Q(n3+n4,2)
                    + 2.*Q(n2,1)*q(n1+n3+n4,3)+2.*p(n1,1)*Q(n2+n3+n4,3)-6.*q(n1+n2+n3+n4,4);
  TComplex *out = (TComplex*) &formula;
  return out;

}
//____________________________________________________________________
