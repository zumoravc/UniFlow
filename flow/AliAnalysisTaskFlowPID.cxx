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

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskFlowPID.h"
#include "AliLog.h" 
class AliAnalysisTaskFlowPID;    

ClassImp(AliAnalysisTaskFlowPID) // classimp: necessary for root

//Double_t AliAnalysisTaskFlowPID::fPtBinEdges[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0}; // You, Katarina binning
Double_t AliAnalysisTaskFlowPID::fPtBinEdges[] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; // PID flow v2 JHEP paper
Double_t AliAnalysisTaskFlowPID::fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
Double_t AliAnalysisTaskFlowPID::fMinvFlowBinEdgesK0s[] = {0.4,0.42,0.44,0.46,0.48,0.49,0.5,0.51,0.52,0.54,0.56,0.58,0.6};
Double_t AliAnalysisTaskFlowPID::fMinvFlowBinEdgesLambda[] = {1.08,1.09,1.10,1.105,1.11,1.115,1.12,1.13,1.14,1.15,1.16};
Int_t AliAnalysisTaskFlowPID::fHarmonics[] = {2};
Double_t AliAnalysisTaskFlowPID::fEtaGap[] = {0.9};

AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID() : AliAnalysisTaskSE(), 
  fAOD(0),
  fPIDResponse(0),
  fTPCPIDResponse(),
  fEtaCutFlag(0),
  fHarmFlag(0), 
  fTrack(0),
  fTrackPt(0),
  fTrackEta(0),
  fTrackPhi(0),
  fNumV0s(0),
  fV0(0),
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
  fVecRefPos(0,0,kFALSE),
  fVecRefNeg(0,0,kFALSE),
  fCountRefPos(0),
  fCountRefNeg(0),
  fQvec2(0),
  fQvec3(0),
  fQvec4(0),
  fQvec5(0),
  fQvec2Gap00P(0),
  fQvec2Gap04P(0),
  fQvec2Gap08P(0),
  fQvec2Gap09P(0),
  fQvec2Gap10P(0),
  fQvec2Gap00N(0),
  fQvec2Gap04N(0),
  fQvec2Gap08N(0),
  fQvec2Gap09N(0),
  fQvec2Gap10N(0),
  fArrTracksFiltered("AliAODTrack",10000),
  fArrV0sK0sFiltered("AliAODv0",5000),
  fArrV0sLambdaFiltered("AliAODv0",5000),
  fArrV0sALambdaFiltered("AliAODv0",5000),

  fOutList(0),
  fOutListV0s(0),
  fOutListQA(0),

  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fCentFlag(0),
  fPVtxCutZ(0.),
  fTrackEtaMax(0),
  fTrackPtMax(0),
  fTrackPtMin(0),
  fNumTPCclsMin(0),
  fTrackFilterBit(0),
  fDiffFlow(0),
  fPID(0),  
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
  fV0sMult(0),
  fV0sInvMassK0sGap00(0),
  fV0sInvMassK0sGap09(0),
  fV0sInvMassLambdaGap00(0),
  fV0sInvMassLambdaGap09(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentSPDvsV0M(0),
  fMultTracksSelected(0),
  fTracksPtCent(0),
  fTracksPt(0),
  fTracksEta(0),
  fTracksPhi(0),
  fTracksCharge(0),
  fRefCorTwo2(0),
  fRefCorTwo3(0),
  fRefCorTwo4(0),
  fRefCorTwo5(0),
  fRefCorTwo2Gap00(0),
  fRefCorTwo2Gap04(0),
  fRefCorTwo2Gap08(0),
  fRefCorTwo2Gap09(0),
  fRefCorTwo2Gap09_test(0),
  fRefCorTwo2Gap10(0),
  fQAPVz(0),
  fQANumTracks(0),
  fQATrackPt(0),
  fQATrackEta(0),
  fQATrackPhi(0),
  fQATrackFilterMap(0),
  fQAV0sCounter(0),
  fQAV0sCounterK0s(0),
  fQAV0sCounterLambda(0),

  fTestTracksPt(0x0),
  fTestTracksMult(0x0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::AliAnalysisTaskFlowPID(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),
  fPIDResponse(0),
  fTPCPIDResponse(),
  fEtaCutFlag(0),
  fHarmFlag(0), 
  fTrack(0),
  fTrackPt(0),
  fTrackEta(0),
  fTrackPhi(0),

  fNumV0s(0),
  fV0(0),
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
  fVecRefPos(0),
  fVecRefNeg(0),
  fCountRefPos(0),
  fCountRefNeg(0),
  fQvec2(0),
  fQvec3(0),
  fQvec4(0),
  fQvec5(0),
  fQvec2Gap00P(0),
  fQvec2Gap04P(0),
  fQvec2Gap08P(0),
  fQvec2Gap09P(0),
  fQvec2Gap10P(0),
  fQvec2Gap00N(0),
  fQvec2Gap04N(0),
  fQvec2Gap08N(0),
  fQvec2Gap09N(0),
  fQvec2Gap10N(0),
  fArrTracksFiltered("AliAODTrack",10000),
  fArrV0sK0sFiltered("AliAODv0",5000),
  fArrV0sLambdaFiltered("AliAODv0",5000),
  fArrV0sALambdaFiltered("AliAODv0",5000),
  
  fAODAnalysis(kTRUE),
  fPbPb(kTRUE),
  fLHC10h(kTRUE),
  fCentFlag(0),
  fPVtxCutZ(0.),
  fTrackEtaMax(0),
  fTrackPtMax(0),
  fTrackPtMin(0),
  fNumTPCclsMin(0),
  fTrackFilterBit(0),
  fDiffFlow(0),
  fPID(0),  
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

  fOutList(0),
  fOutListV0s(0),
  fOutListQA(0),
  fEventCounter(0),

  fV0sMult(0),

  fV0sInvMassK0sGap00(0),
  fV0sInvMassK0sGap09(0),
  fV0sInvMassLambdaGap00(0),
  fV0sInvMassLambdaGap09(0),
  fEventMult(0),
  fCentralityDis(0),
  fCentSPDvsV0M(0),
  fMultTracksSelected(0),
  fTracksPtCent(0),
  fTracksPt(0),
  fTracksEta(0),
  fTracksPhi(0),
  fTracksCharge(0),
  fRefCorTwo2(0),
  fRefCorTwo3(0),
  fRefCorTwo4(0),
  fRefCorTwo5(0),
  fRefCorTwo2Gap00(0),
  fRefCorTwo2Gap04(0),
  fRefCorTwo2Gap08(0),
  fRefCorTwo2Gap09(0),
  fRefCorTwo2Gap09_test(0),
  fRefCorTwo2Gap10(0),
  fQAPVz(0),
  fQANumTracks(0),
  fQATrackPt(0),
  fQATrackEta(0),
  fQATrackPhi(0),
  fQATrackFilterMap(0),
  fQAV0sCounter(0),
  fQAV0sCounterK0s(0),
  fQAV0sCounterLambda(0),

  fTestTracksPt(0x0),
  fTestTracksMult(0x0)
{
  // constructor
  for(Int_t i(0); i < fNumCentBins; i++)
  {
    for(Int_t j(0); j < fNumHarmonics; j++)
    {
      for(Int_t k(0); k < fNumEtaGap; k++)
      {
        fV0sDiffTwoPos_K0s[i][j][k] = 0x0;
        fV0sDiffTwoNeg_K0s[i][j][k] = 0x0;
        fV0sDiffTwoPos_Lambda[i][j][k] = 0x0;
        fV0sDiffTwoNeg_Lambda[i][j][k] = 0x0;
        fV0sK0s[i][j][k] = 0x0;
        fV0sLambda[i][j][k] = 0x0;

        fV0sInvMassK0s[i][j][k] = 0x0;
        fV0sInvMassLambda[i][j][k] = 0x0;

      }
    }
    
    fDiffCorTwo2[i] = 0;
    fDiffCorTwo2Gap00[i] = 0;
    fDiffCorTwo2Gap04[i] = 0;
    fDiffCorTwo2Gap08[i] = 0;
    fDiffCorTwo2Gap10[i] = 0;
    fDiffCorTwo3[i] = 0;
  }

  // QA plots
  for (Int_t i(0); i < 2; i++)
  {
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
}
//_____________________________________________________________________________
AliAnalysisTaskFlowPID::~AliAnalysisTaskFlowPID()
{
  // destructor
  if(fOutList) 
  {
      delete fOutList;     // at the end of your task, it is deleted from memory by calling this function
  }

  if(fOutListV0s)
  {
    delete fOutListV0s;
  }

  if(fOutListQA)
  {
    delete fOutListQA;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you create the histograms that you want to use 
  // the histograms are in this case added to a TList, this list is in the end saved to an output file

  fOutList = new TList();          // this is a list which will contain all of your histograms at the end of the analysis, the contents of this list are written to the output file
  fOutList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested (dont worry about this now)

  fOutListV0s = new TList();
  fOutListV0s->SetOwner(kTRUE);

  fOutListQA = new TList();
  fOutListQA->SetOwner(kTRUE);

  // test histos
  fTestTracksPt = new TH1D("fTestTracksPt", "Tracks #it{p}_{T} (selected); #it{p}^{track}_{T} (GeV/#it{c});", 100, 0, 10);    
  fOutList->Add(fTestTracksPt);
  fTestTracksMult = new TH1D("fTestTracksMult","Track multiplicity (selected tracks in selected events); tracks;",100,0,5000);
  fOutList->Add(fTestTracksMult);
  
 
  // main output
  fEventMult = new TH1D("fEventMult","Track multiplicity (all tracks in selected events); tracks;",100,0,10000);
  fOutList->Add(fEventMult);
  fCentralityDis = new TH1D("fCentralityDis", "centrality distribution; centrality;", fNumCentBins,fCentBinEdges);
  fOutList->Add(fCentralityDis);
  fMultTracksSelected = new TH1D("fMultTracksSelected","Track multiplicity (selected tracks in selected events); tracks;",100,0,5000);
  fOutList->Add(fMultTracksSelected);
  fTracksPtCent = new TH2D("fTracksPtCent", "Tracks #it{p}_{T} vs. centrality (selected); #it{p}^{track}_{T} (GeV/#it{c}); centrality;", fNumPtBins,fPtBinEdges,fNumCentBins,fCentBinEdges);    
  fOutList->Add(fTracksPtCent);          
  fTracksPt = new TH1D("fTracksPt", "Tracks #it{p}_{T} (selected); #it{p}^{track}_{T} (GeV/#it{c});", 100, 0, 10);    
  fOutList->Add(fTracksPt);          
  fTracksEta = new TH1D("fTracksEta", "Tracks #it{#eta} (selected); #it{#eta}^{track};", 400, -2, 2);    
  fOutList->Add(fTracksEta);          
  fTracksPhi = new TH1D("fTracksPhi", "Tracks #it{#varphi} (selected); #it{#varphi}^{track};", 360, 0., TMath::TwoPi());    
  fOutList->Add(fTracksPhi);          
  fTracksCharge = new TH1D("fTracksCharge", "Track charge (selected); charge^{track};", 3,-1.5,1.5);    
  fOutList->Add(fTracksCharge);          
  
  // V0s histos

  fV0sMult = new TH1D("fV0sMult","V0s multiplicity (in selected events); V0s;",200,0,1000);
  fOutListV0s->Add(fV0sMult);
  fV0sInvMassK0sGap00 = new TH1D("fV0sInvMassK0sGap00","K^{0}_{S} InvMass |#Delta#it{#eta}| > 0 (selected); #it{m}_{inv} (GeV/#it{c}^{2});", 200,fV0MinMassK0s,fV0MaxMassK0s);          
  fOutListV0s->Add(fV0sInvMassK0sGap00);
  fV0sInvMassK0sGap09 = new TH1D("fV0sInvMassK0sGap09","K^{0}_{S} InvMass |#Delta#it{#eta}| > 0.9 (selected); #it{m}_{inv} (GeV/#it{c}^{2});", 200,fV0MinMassK0s,fV0MaxMassK0s);          
  fOutListV0s->Add(fV0sInvMassK0sGap09);
  fV0sInvMassK0sGap09_test = new TH1D("fV0sInvMassK0sGap09_test","K^{0}_{S} InvMass |#Delta#it{#eta}| > 0.9 (selected); #it{m}_{inv} (GeV/#it{c}^{2});", 200,fV0MinMassK0s,fV0MaxMassK0s);          
  fOutListV0s->Add(fV0sInvMassK0sGap09_test);
  fV0sInvMassLambdaGap00 = new TH1D("fV0sInvMassLambdaGap00","#Lambda+#bar{#Lambda} InvMass |#Delta#it{#eta}| > 0 (selected); #it{m}_{inv} (GeV/#it{c}^{2});", 80,fV0MinMassLambda,fV0MaxMassLambda);          
  fOutListV0s->Add(fV0sInvMassLambdaGap00);
  fV0sInvMassLambdaGap09 = new TH1D("fV0sInvMassLambdaGap09","#Lambda+#bar{#Lambda} InvMass |#Delta#it{#eta}| > 0.9 (selected); #it{m}_{inv} (GeV/#it{c}^{2});", 80,fV0MinMassLambda,fV0MaxMassLambda);          
  fOutListV0s->Add(fV0sInvMassLambdaGap09);
  
  for(Int_t i(0); i < fNumCentBins; i++)
  {
    for(Int_t j(0); j < fNumHarmonics; j++)
    {
      for(Int_t k(0); k < fNumEtaGap; k++)
      {
        fV0sDiffTwoPos_K0s[i][j][k] = new TProfile2D(Form("fV0sDiffTwoPos_K0s_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i), Form("K^{0}_{S} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsK0s,fMinvFlowBinEdgesK0s);
        fV0sDiffTwoPos_K0s[i][j][k]->Sumw2();
        fOutListV0s->Add(fV0sDiffTwoPos_K0s[i][j][k]);
  
        fV0sDiffTwoNeg_K0s[i][j][k] = new TProfile2D(Form("fV0sDiffTwoNeg_K0s_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i), Form("K^{0}_{S} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, #it{#eta}^{POI} < -%g} Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsK0s,fMinvFlowBinEdgesK0s);
        fV0sDiffTwoNeg_K0s[i][j][k]->Sumw2();
        fOutListV0s->Add(fV0sDiffTwoNeg_K0s[i][j][k]);

        fV0sK0s[i][j][k] = new TH2D(Form("fV0s_K0s_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i),Form("K^{0}_{S} candidates n=%d |#Delta#it{#eta}| > %g Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{m}^{V0}_{inv} (GeV/#it{c}^{2});",fHarmonics[j],fEtaGap[k],fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins, fPtBinEdges, 200, fV0MinMassK0s, fV0MaxMassK0s);
        fV0sK0s[i][j][k]->Sumw2();
        fOutListV0s->Add(fV0sK0s[i][j][k]);

        fV0sInvMassK0s[i][j][k] = new TH1D(Form("fV0sInvMass_K0s_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i),Form("K^{0}_{S} InvMass n=%d |#Delta#it{#eta}| > %g Cent %g-%g%% (selected); #it{m}_{inv} (GeV/#it{c}^{2});",fHarmonics[j],fEtaGap[k],fCentBinEdges[i],fCentBinEdges[i+1]), 200,fV0MinMassK0s,fV0MaxMassK0s);          
        fOutListV0s->Add(fV0sInvMassK0s[i][j][k]);

        fV0sDiffTwoPos_Lambda[i][j][k] = new TProfile2D(Form("fV0sDiffTwoPos_Lambda_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i), Form("#Lambda+#bar{#Lambda} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, #it{#eta}^{POI} > %g} Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsLambda,fMinvFlowBinEdgesLambda);
        fV0sDiffTwoPos_Lambda[i][j][k]->Sumw2();
        fOutListV0s->Add(fV0sDiffTwoPos_Lambda[i][j][k]);
  
        fV0sDiffTwoNeg_Lambda[i][j][k] = new TProfile2D(Form("fV0sDiffTwoNeg_Lambda_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i), Form("#Lambda+#bar{#Lambda} #LT#LT2'#GT#GT_{%d,|#Delta#it{#eta}| > %g, it{#eta}^{POI} < -%g} Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{M}_{inv}^{V0} (GeV/#it{c}^{2})",fHarmonics[j],fEtaGap[k],fEtaGap[k]/2,fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges,fNumMinvFlowBinsLambda,fMinvFlowBinEdgesLambda);
        fV0sDiffTwoNeg_Lambda[i][j][k]->Sumw2();
        fOutListV0s->Add(fV0sDiffTwoNeg_Lambda[i][j][k]);

        fV0sLambda[i][j][k] = new TH2D(Form("fV0s_Lambda_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i),Form("#Lambda+#bar{#Lambda} candidates n=%d |#Delta#it{#eta}| > %g Cent %g-%g%%; #it{p}^{V0}_{T} (GeV/#it{c}); #it{m}^{V0}_{inv} (GeV/#it{c}^{2});",fHarmonics[j],fEtaGap[k],fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins, fPtBinEdges, 200, fV0MinMassLambda, fV0MaxMassLambda);
        fV0sLambda[i][j][k]->Sumw2();
        fOutListV0s->Add(fV0sLambda[i][j][k]);

        fV0sInvMassLambda[i][j][k] = new TH1D(Form("fV0sInvMass_Lambda_n%d_Gap%02.2g_Cent%d",fHarmonics[j],10*fEtaGap[k],i),Form("#Lambda+#bar{#Lambda} InvMass n=%d |#Delta#it{#eta}| > %g Cent %g-%g%% (selected); #it{m}_{inv} (GeV/#it{c}^{2});",fHarmonics[j],fEtaGap[k],fCentBinEdges[i],fCentBinEdges[i+1]), 200,fV0MinMassLambda,fV0MaxMassLambda);          
        fOutListV0s->Add(fV0sInvMassLambda[i][j][k]);
      }
    }
  }

  
  
  fRefCorTwo2 = new TProfile("fRefCorTwo2","#LT#LT2#GT#GT_{2} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2->Sumw2();
  fOutList->Add(fRefCorTwo2);
  fRefCorTwo3 = new TProfile("fRefCorTwo3","#LT#LT2#GT#GT_{3} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo3->Sumw2();
  fOutList->Add(fRefCorTwo3);
  fRefCorTwo4 = new TProfile("fRefCorTwo4","#LT#LT2#GT#GT_{4} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo4->Sumw2();
  fOutList->Add(fRefCorTwo4);
  fRefCorTwo5 = new TProfile("fRefCorTwo5","#LT#LT2#GT#GT_{5} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo5->Sumw2();
  fOutList->Add(fRefCorTwo5);

  fRefCorTwo2Gap00 = new TProfile("fRefCorTwo2_Gap00","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap00->Sumw2();
  fOutList->Add(fRefCorTwo2Gap00);
  fRefCorTwo2Gap04 = new TProfile("fRefCorTwo2_Gap04","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.4} (ref. flow); centrality",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap04->Sumw2();
  fOutList->Add(fRefCorTwo2Gap04);
  fRefCorTwo2Gap08 = new TProfile("fRefCorTwo2_Gap08","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.8} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap08->Sumw2();
  fOutList->Add(fRefCorTwo2Gap08);
  fRefCorTwo2Gap09 = new TProfile("fRefCorTwo2_Gap09","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.9} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap09->Sumw2();
  fOutList->Add(fRefCorTwo2Gap09);
    fRefCorTwo2Gap09_test = new TProfile("fRefCorTwo2_Gap09_test","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 0.9} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap09_test->Sumw2();
  fOutList->Add(fRefCorTwo2Gap09_test);
  fRefCorTwo2Gap10 = new TProfile("fRefCorTwo2_Gap10","#LT#LT2#GT#GT_{2,|#Delta#it{#eta}| > 1} (ref. flow); centrality;",fNumCentBins,fCentBinEdges);
  fRefCorTwo2Gap10->Sumw2();
  fOutList->Add(fRefCorTwo2Gap10);
  
  if(fDiffFlow) // do differential flow switch
  {
    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2[i] = new TProfile(Form("fDiffCorTwo2_Cent%d",i),Form("#LT#LT2'#GT#GT_{2} Cent %g-%g%% (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2[i]->Sumw2();
      fOutList->Add(fDiffCorTwo2[i]);
    }
    
    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap00[i] = new TProfile(Form("fDiffCorTwo2_Gap00_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0} Cent %g-%g%% |#it{#eta}^{POI}|>0 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap00[i]->Sumw2();
      fOutList->Add(fDiffCorTwo2Gap00[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap04[i] = new TProfile(Form("fDiffCorTwo2_Gap04_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0.4} Cent %g-%g%% #it{#eta}^{POI}>0.2 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap04[i]->Sumw2();
      fOutList->Add(fDiffCorTwo2Gap04[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap08[i] = new TProfile(Form("fDiffCorTwo2_Gap08_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 0.8} Cent %g-%g%% #it{#eta}^{POI}>0.4 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap08[i]->Sumw2();
      fOutList->Add(fDiffCorTwo2Gap08[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo2Gap10[i] = new TProfile(Form("fDiffCorTwo2_Gap10_Cent%d",i),Form("#LT#LT2'#GT#GT_{2,|#Delta#it{#eta}| > 1} Cent %g-%g%% #it{#eta}^{POI}>0.5 (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo2Gap10[i]->Sumw2();
      fOutList->Add(fDiffCorTwo2Gap10[i]);
    }

    for(Int_t i = 0; i < fNumCentBins; i++)
    {
      fDiffCorTwo3[i] = new TProfile(Form("fDiffCorTwo3_Cent%d",i),Form("#LT#LT2'#GT#GT_{3} Cent %g-%g%% (diff. flow); #it{p}^{track}_{T} (GeV/#it{c})",fCentBinEdges[i],fCentBinEdges[i+1]),fNumPtBins,fPtBinEdges);
      fDiffCorTwo3[i]->Sumw2();
      fOutList->Add(fDiffCorTwo3[i]);
    }
    
  }
  

  // QA output
  Int_t iNEventCounterBins = 9;
  TString sEventCounterLabel[] = {"Input","AOD OK","Pile-up OK","PV OK","SPD Vtx OK","PV #it{z} OK","Centrality OK","At least 2 selected tracks","At least 1 V0 candidate"};
  fEventCounter = new TH1D("fEventCounter","Event Counter",iNEventCounterBins,0,iNEventCounterBins);
  for(Int_t i = 0; i < iNEventCounterBins; i++)
    fEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
  fOutListQA->Add(fEventCounter);
  

  fQAPVz = new TH1D("fQAPVz","QA: PV #it{z}; #it{z} (cm);",100,-50,50);
  fOutListQA->Add(fQAPVz);
  fCentSPDvsV0M = new TH2D("fCentSPDvsV0M", "V0M-cent vs SPD-cent; V0M; SPD-cent", 100, 0, 100, 100, 0, 100);
  fOutListQA->Add(fCentSPDvsV0M);
  fQANumTracks = new TH1D("fQANumTracks","QA: Number of AOD tracks; tracks;",100,0,10000);
  fOutListQA->Add(fQANumTracks);
  fQATrackPt = new TH1D("fQATrackPt","QA: Track #it{p}_{T} (all); #it{p}^{track}_{T} (GeV/#it{c});",100,0,10);
  fOutListQA->Add(fQATrackPt);
  fQATrackEta = new TH1D("fQATrackEta","QA: Track #it{#eta} (all); #it{#eta}^{track};",300,-1.5,1.5);
  fOutListQA->Add(fQATrackEta);
  fQATrackPhi = new TH1D("fQATrackPhi","QA: Track #it{#varphi} (all); #it{#varphi}^{track};",300,0,TMath::TwoPi());
  fOutListQA->Add(fQATrackPhi);
	fQATrackFilterMap = new TH1D("fQATrackFilterMap","QA: Tracks filter map (all); filter bit;",1000,0,1000);
	fOutListQA->Add(fQATrackFilterMap);

  // QA V0s output
  Int_t iNV0sCounterBins = 18;
  TString sV0sCounterLabel[] = {"Input","Daughters exists","Reco. method","Mother acceptance (#it{#eta},#it{p}_{T})","Decay radius","TPC refit","not Kink","Daughters charge","Daughter acceptance (#it{#eta},#it{p}_{T})","DCA to PV","Daughters DCA","Selected","K^{0}_{S}","#Lambda","#bar{#Lambda}","#Lambda && #bar{#Lambda}","K0s && (#Lambda || #bar{#Lambda})","K0s && #Lambda && #bar{#Lambda}"};
  fQAV0sCounter = new TH1D("fQAV0sCounter","QA V^{0}_{S} counter",iNV0sCounterBins,0,iNV0sCounterBins);
  for(Int_t i = 0; i < iNV0sCounterBins; i++)
    fQAV0sCounter->GetXaxis()->SetBinLabel(i+1, sV0sCounterLabel[i].Data() );
  fOutListQA->Add(fQAV0sCounter);

	Int_t iNV0sCounterK0sBins = 7;
  TString sV0sCounterK0sLabel[] = {"Input","Inv. Mass","Mother #it{y} acceptance","CPA","Life-time","Armenteros","Selected"};
  fQAV0sCounterK0s = new TH1D("fQAV0sCounterK0s","QA K^{0}_{S} counter",iNV0sCounterK0sBins,0,iNV0sCounterK0sBins);
  for(Int_t i = 0; i < iNV0sCounterK0sBins; i++)
    fQAV0sCounterK0s->GetXaxis()->SetBinLabel(i+1, sV0sCounterK0sLabel[i].Data() );
  fOutListQA->Add(fQAV0sCounterK0s);

	Int_t iNV0sCounterLambdaBins = 7;
  TString sV0sCounterLambdaLabel[] = {"Input","Inv. Mass","Mother #it{y} acceptance","CPA","Life-time","proton PID","Selected"};
  fQAV0sCounterLambda = new TH1D("fQAV0sCounterLambda","QA #Lambda/#bar{#Lambda} counter",iNV0sCounterLambdaBins,0,iNV0sCounterLambdaBins);
  for(Int_t i = 0; i < iNV0sCounterLambdaBins; i++)
    fQAV0sCounterLambda->GetXaxis()->SetBinLabel(i+1, sV0sCounterLambdaLabel[i].Data() );
  fOutListQA->Add(fQAV0sCounterLambda);


	TString sQAlabel[fQAV0sNumSteps] = {"Before","After"/*,"Test"*/};
	for(Int_t i(0); i < fQAV0sNumSteps; i++)
	{
  	fQAV0sRecoMethod[i] = new TH1D(Form("fQAV0sRecoMethod_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Reconstruction method (%s cuts)",sQAlabel[i].Data()), 2,-0.5,1.5);
  	fQAV0sRecoMethod[i]->GetXaxis()->SetBinLabel(1, "offline");
  	fQAV0sRecoMethod[i]->GetXaxis()->SetBinLabel(2, "online (on-the-fly)");
  	fOutListQA->Add(fQAV0sRecoMethod[i]);	
		fQAV0sTPCRefit[i] = new TH1D(Form("fQAV0sTPCRefit_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: TPC refit (%s cuts)",sQAlabel[i].Data()), 2,-0.5,1.5);
		fQAV0sTPCRefit[i]->GetXaxis()->SetBinLabel(1, "not refit");
  	fQAV0sTPCRefit[i]->GetXaxis()->SetBinLabel(2, "refit");
  	fOutListQA->Add(fQAV0sTPCRefit[i]);	
		fQAV0sKinks[i] = new TH1D(Form("fQAV0sKinks_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Kinks (%s cuts)",sQAlabel[i].Data()), 2,-0.5,1.5);
		fQAV0sKinks[i]->GetXaxis()->SetBinLabel(1, "NOT AliAODVertex::kKink");
  	fQAV0sKinks[i]->GetXaxis()->SetBinLabel(2, "AliAODVertex:kKink");
  	fOutListQA->Add(fQAV0sKinks[i]);	
		fQAV0sDCAtoPV[i] = new TH1D(Form("fQAV0sDCAtoPV_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter DCA to PV (%s cuts); daughter DCA^{PV} (cm)",sQAlabel[i].Data()), 100,0.,10.);	
  	fOutListQA->Add(fQAV0sDCAtoPV[i]);	
		fQAV0sDCADaughters[i] = new TH1D(Form("fQAV0sDCADaughters_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: DCA among daughters (%s cuts); DCA^{daughters} (cm)",sQAlabel[i].Data()), 100,0.,10.);	
  	fOutListQA->Add(fQAV0sDCADaughters[i]);	
		fQAV0sDecayRadius[i] = new TH1D(Form("fQAV0sDecayRadius_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Decay radius (%s cuts); #it{r_{xy}}^{decay} (cm)",sQAlabel[i].Data()), 200,0.,500.);	
  	fOutListQA->Add(fQAV0sDecayRadius[i]);	
    fQAV0sDaughterPt[i] = new TH1D(Form("fQAV0sDaughterPt_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter #it{p}_{T} (%s cuts); #it{p}_{T}^{daughter} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,10.);  
    fOutListQA->Add(fQAV0sDaughterPt[i]);
    fQAV0sDaughterPhi[i] = new TH1D(Form("fQAV0sDaughterPhi_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter #it{#varphi} (%s cuts); #it{#varphi}^{daughter} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,TMath::TwoPi()); 
    fOutListQA->Add(fQAV0sDaughterPhi[i]);   
    fQAV0sDaughterTPCdEdxPt[i] = new TH2D(Form("fQAV0sDaughterTPCdEdxPt_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: TPC dEdx daughters (%s cuts); #it{p}_{T}^{daughter} (GeV/#it{c}); TPC dEdx (au);",sQAlabel[i].Data()), 100,0.,10,200,0.,200.);
  	fOutListQA->Add(fQAV0sDaughterTPCdEdxPt[i]);	
    fQAV0sDaughterEta[i] = new TH1D(Form("fQAV0sDaughterEta_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Daughter #it{#eta} (%s cuts); #it{#eta}^{daugter}",sQAlabel[i].Data()), 400,-2,2);  
    fOutListQA->Add(fQAV0sDaughterEta[i]);  
    fQAV0sMotherPt[i] = new TH1D(Form("fQAV0sMotherPt_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{p}_{T} (%s cuts); #it{p}_{T}^{V0} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,10.);  
    fOutListQA->Add(fQAV0sMotherPt[i]); 
    fQAV0sMotherPhi[i] = new TH1D(Form("fQAV0sMotherPhi_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{#varphi} (%s cuts); #it{#varphi}^{V0} (GeV/#it{c})",sQAlabel[i].Data()), 100,0.,TMath::TwoPi()); 
    fOutListQA->Add(fQAV0sMotherPhi[i]);  
    fQAV0sMotherEta[i] = new TH1D(Form("fQAV0sMotherEta_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{#eta} (%s cuts); #it{#eta}^{V0}",sQAlabel[i].Data()), 400,-2,2); 
    fOutListQA->Add(fQAV0sMotherEta[i]);  
    fQAV0sMotherRapK0s[i] = new TH1D(Form("fQAV0sMotherRapK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{y} (K^{0}_{S} hypo) (%s cuts); #it{y}^{V0,K0s}",sQAlabel[i].Data()), 400,-2,2);
    fOutListQA->Add(fQAV0sMotherRapK0s[i]);
    fQAV0sMotherRapLambda[i] = new TH1D(Form("fQAV0sMotherRapLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Mother #it{y} (Lambda/#bar{#Lambda} hypo) (%s cuts); #it{y}^{V0,#Lambda}",sQAlabel[i].Data()),400,-2,2);
    fOutListQA->Add(fQAV0sMotherRapLambda[i]);     
    fQAV0sInvMassK0s[i] = new TH1D(Form("fQAV0sInvMassK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: K^{0}_{S}: InvMass (%s cuts); #it{m}_{inv} (GeV/#it{c}^{2});",sQAlabel[i].Data()), 200,fV0MinMassK0s,fV0MaxMassK0s);
    fOutListQA->Add(fQAV0sInvMassK0s[i]);  
    fQAV0sInvMassLambda[i] = new TH1D(Form("fQAV0sInvMassLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: InvMass (%s cuts); #it{m}_{inv} (GeV/#it{c}^{2});",sQAlabel[i].Data()), 80,fV0MinMassLambda,fV0MaxMassLambda);
    fOutListQA->Add(fQAV0sInvMassLambda[i]);  
    fQAV0sCPAK0s[i] = new TH1D(Form("fQAV0sCPAK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: K^{0}_{S}: CPA (%s cuts); CPA^{K0s}",sQAlabel[i].Data()), 100,0.9,1.);  
    fOutListQA->Add(fQAV0sCPAK0s[i]); 
    fQAV0sCPALambda[i] = new TH1D(Form("fQAV0sCPALambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: CPA (%s cuts); CPA^{#Lambda}",sQAlabel[i].Data()), 100, 0.9,1.); 
    fOutListQA->Add(fQAV0sCPALambda[i]);  
    fQAV0sNumTauK0s[i] = new TH1D(Form("fQAV0sNumTauK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}:  K^{0}_{S}: Number of #it{c#tau} (%s cuts); #it{c#tau}^{K0s} (cm)",sQAlabel[i].Data()), 100, 0.,10.);
    fOutListQA->Add(fQAV0sNumTauK0s[i]);  
    fQAV0sNumTauLambda[i] = new TH1D(Form("fQAV0sNumTauLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: Number of #it{c#tau} (%s cuts); #it{c#tau}^{#Lambda} (cm)",sQAlabel[i].Data()), 300, 0.,30);
    fOutListQA->Add(fQAV0sNumTauLambda[i]);
    fQAV0sArmenterosK0s[i] = new TH2D(Form("fQAV0sArmenterosK0s_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}:  K^{0}_{S}: Armenteros-Podolaski plot (%s cuts); #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});",sQAlabel[i].Data()), 100,-1.,1., 100,0.,0.3);
    fOutListQA->Add(fQAV0sArmenterosK0s[i]);
    fQAV0sArmenterosLambda[i] = new TH2D(Form("fQAV0sArmenterosLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: Armenteros-Podolaski plot (%s cuts); #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});",sQAlabel[i].Data()), 100,-1.,1., 100,0.,0.3);
    fOutListQA->Add(fQAV0sArmenterosLambda[i]);
    fQAV0sProtonNumSigmaPtLambda[i] = new TH2D(Form("fQAV0sProtonNumSigmaPtLambda_%s",sQAlabel[i].Data()),Form("QA V^{0}_{S}: #Lambda/#bar{#Lambda}: (anti-)proton PID (%s cuts); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});",sQAlabel[i].Data()), 100,0.,10,100,0.,10.);
    fOutListQA->Add(fQAV0sProtonNumSigmaPtLambda[i]);

	}

  PostData(1, fOutList);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
	PostData(2, fOutListV0s);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
	PostData(3, fOutListQA);           // postdata will notify the analysis manager of changes / updates to the fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::UserExec(Option_t *)
{
  // this function is called once for each event

	fEventCounter->Fill(0); // input event

	// loading PID response for protons
  if((fCutV0ProtonNumSigmaMax > 0.) && (fCutV0ProtonPIDPtMax > 0.))
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
  
  // cleaning all TClonesArray containers
  fArrTracksFiltered.Clear("C");
  fArrV0sK0sFiltered.Clear("C");
  fArrV0sLambdaFiltered.Clear("C");
  fArrV0sALambdaFiltered.Clear("C");

  // filtering objects and filling relevant containers
  FilterTracks();
  FilterV0s();

 
  //Double_t dEtaGap2 = 1.2;
  //TString sTest = Form("%02.2g %02.2g ",10*dEtaGap, 10*dEtaGap2);
  //printf("Test EtaGap (double) %f / (string) %s / (string *10) \n",dEtaGap,sTest.Data());





  // begin of loop over EtaGap & Harmonics
  Int_t iHarm = 0;
  Double_t dEtaGap = 0;
  TProfile2D* profV0sPos = 0x0;
  TProfile2D* profV0sNeg = 0x0;
  for(Int_t i(0); i < fNumHarmonics; i++)
  {
    for(Int_t j(0); j < fNumEtaGap; j++)
    {
      iHarm = fHarmonics[i];
      dEtaGap = fEtaGap[j];

      FillRefFlowVectors(dEtaGap, iHarm);
      EstimateRefCumulant(dEtaGap, iHarm, fRefCorTwo2Gap09_test);

      if(fDiffFlow)
      {
        //EstimateRefPtDiffCumulant(0.4,2,fDiffCorTwo2Gap04[fCentBinIndex],fDiffCorTwo2Gap00[fCentBinIndex]);  
      }
      
      if(fPID)
      {
        // output & profile finding testing
        //TString sNameProfPos = Form("fV0sDiffTwo%d_Gap%02.2gP_K0s_Cent%d",iHarm,10*0.,fCentBinIndex);
        //TString sNameProfNeg = Form("fV0sDiffTwo%d_Gap%02.2gN_K0s_Cent%d",iHarm,10*dEtaGap,fCentBinIndex);
        //printf("%s\n",sNameProfPos.Data());
        //profV0sPos = (TProfile2D*) fOutListV0s->FindObject(sNameProfPos.Data());
        //profV0sNeg = (TProfile2D*) fOutListV0s->FindObject(sNameProfNeg.Data());
        
        //profV0sPos = (TProfile2D*) fOutListV0s->FindObject(Form("fV0sDiffTwo%d_Gap%02.2gP_K0s_Cent%d",iHarm,10*dEtaGap,fCentBinIndex));
        //profV0sNeg = (TProfile2D*) fOutListV0s->FindObject(Form("fV0sDiffTwo%d_Gap%02.2gN_K0s_Cent%d",iHarm,10*dEtaGap,fCentBinIndex));
        //if(!profV0sPos || !profV0sNeg)
        //  continue;
        
        //EstimateV0Cumulant(dEtaGap, iHarm, fV0sDiffTwoPos_K0s, fV0);
        EstimateV0Cumulant(i,j);
      }
  
    } // end of loop over eta gap
  } // end of loop over harmonics
 
  // end of potential loop over EtaGap & Harmonics

  PostData(1, fOutList);  // stream the results the analysis of this event to the output manager which will take care of writing it to a file
  PostData(2, fOutListV0s);
  PostData(3, fOutListQA);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::EstimateRefPtDiffCumulant(const Float_t dEtaGap, const Short_t iHarm, TProfile* profilePos, TProfile* profileNeg)
{
  // checking if pointers are valid
  if(!profilePos)
    return;

  if(!profileNeg && dEtaGap >= 0.0) // if with eta gap, negative profile has to be valid as well
    return;

  // check if RFPs flow vectors are empty
  if(!AreRefFlowVectorsFilled(dEtaGap,iHarm)) // if yes fill ref flow vectors
    FillRefFlowVectors(dEtaGap, iHarm);

  // POIs
  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for POI
  const Int_t iNumPOIs = fArrTracksFiltered.GetEntries();
  
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
    track = static_cast<AliAODTrack*>(fArrTracksFiltered.At(i));

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

/*
void AliAnalysisTaskFlowPID::EstimateV0Cumulant(const Float_t dEtaGap, const Short_t iHarm, const Short_t iSpecies, TProfile2D* profilePos, TProfile2D* profileNeg)
{
  // checking if pointers are valid
  if(!profilePos)
    return;

  if(!profileNeg && dEtaGap >= 0.0) // if with eta gap, negative profile has to be valid as well
    return;

  // check if RFPs flow vectors are empty
  if(!AreRefFlowVectorsFilled(dEtaGap,iHarm)) // if yes fill ref flow vectors
    FillRefFlowVectors(dEtaGap, iHarm);

  // POIs

  TClonesArray arrayV0;
  if(iSpecies == 1) // K0s
  {
    arrayV0 = fArrV0sK0sFiltered;
  }
  else if (iSpecies == 2) // (anti)Lambda
  {
    arrayV0 = fArrV0sLambdaFiltered;
  }
  else 
    return;

  const Double_t dEtaGapCut = dEtaGap / 2; // eta cut limit for POI
  const Int_t iNumV0s = arrayV0.GetEntries();
  
  TComplex vecVpos[fNumPtBins]; // flow vector for POIs in positive eta (including eta gap cut) [or all POIs with no eta gap]
  TComplex vecVneg[fNumPtBins]; // flow vector for POIs in negative eta (including eta gap cut)
  Int_t iCountVpos[fNumPtBins]; // counter of POIs in positive eta [or all POIs with no eta gap]
  Int_t iCountVneg[fNumPtBins]; // counter of POIs in negative eta
  
  AliAODv0* v0 = 0x0;
  Double_t dPhi = 0, dEta = 0;
  Double_t dWeight = 0, dAmp = 0, dVal = 0;
  Short_t iMinvFlowBinIndexK0s = -1, iPtBinIndex = -1;

  for(Int_t j(0); j < fNumMinvFlowBinsK0s; j++) // loop over Minv bins
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
      v0 = static_cast<AliAODv0*>(arrayV0.At(i));

      if(!v0)
        continue;

      iMinvFlowBinIndexK0s = GetMinvFlowBinIndexK0s(v0->MassK0Short());
      iPtBinIndex = GetPtBinIndex(v0->Pt());

      if(iPtBinIndex == -1 || iMinvFlowBinIndexK0s != j)
        continue;

      dPhi = v0->Phi();
      dEta = v0->Eta();

      // filling the flow vectors for POIs
      if(dEtaGap < 0.) // no eta gap
      {
        vecVpos[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
        iCountVpos[iPtBinIndex]++;
      }
      else // with eta gap
      {
        if(dEta > dEtaGapCut)
        {
          vecVpos[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
          iCountVpos[iPtBinIndex]++;
          fV0sK0sGap09_test[fCentBinIndex]->Fill(v0->Pt(), v0->MassK0Short());     
          fV0sInvMassK0sGap09_test->Fill(v0->MassK0Short());
        }

        if(dEta < -dEtaGapCut)
        {
          vecVneg[iPtBinIndex] += TComplex(TMath::Cos(iHarm*dPhi),TMath::Sin(iHarm*dPhi),kFALSE);
          iCountVneg[iPtBinIndex]++;
          fV0sK0sGap09_test[fCentBinIndex]->Fill(v0->Pt(), v0->MassK0Short()); 
          fV0sInvMassK0sGap09_test->Fill(v0->MassK0Short());
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
          profilePos->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, (fMinvFlowBinEdgesK0s[j+1] + fMinvFlowBinEdgesK0s[j])/2 ,dVal, dWeight);
      }
      else // with eta gap
      {
        dWeight = iCountVpos[i] * (fCountRefNeg);
        dAmp = (vecVpos[i] * (TComplex::Conjugate(fVecRefNeg))).Re();
        dVal = dAmp / dWeight;
        if( TMath::Abs(dVal < 1) && (dWeight > 0))
          profilePos->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, (fMinvFlowBinEdgesK0s[j+1] + fMinvFlowBinEdgesK0s[j])/2 ,dVal, dWeight);
        
        dWeight = iCountVneg[i] * (fCountRefPos);
        dAmp = (vecVneg[i] * (TComplex::Conjugate(fVecRefPos))).Re();
        dVal = dAmp / dWeight;
        if( TMath::Abs(dVal < 1) && (dWeight > 0))
          profileNeg->Fill( (fPtBinEdges[i+1] + fPtBinEdges[i])/2, (fMinvFlowBinEdgesK0s[j+1] + fMinvFlowBinEdgesK0s[j])/2 ,dVal, dWeight);
      }
    }
  } // end of loop over Minv bins

  return; 
}
*/
//_____________________________________________________________________________
void AliAnalysisTaskFlowPID::EstimateV0Cumulant(Short_t iEtaGapIndex, Short_t iHarmonicsIndex)
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
        profV0sPos = fV0sDiffTwoPos_K0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        profV0sNeg = fV0sDiffTwoNeg_K0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hPtInvMass = fV0sK0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hInvMass = fV0sInvMassK0s[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        iNumMinvFlowBins = fNumMinvFlowBinsK0s;
        dInvMassBinEdges = fMinvFlowBinEdgesK0s;
      break;
        
      case 1: // Lambda
        arrayV0s = fArrV0sLambdaFiltered;
        profV0sPos = fV0sDiffTwoPos_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        profV0sNeg = fV0sDiffTwoNeg_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hPtInvMass = fV0sLambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hInvMass = fV0sInvMassLambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        iNumMinvFlowBins = fNumMinvFlowBinsLambda;
        dInvMassBinEdges = fMinvFlowBinEdgesLambda;
      break;

      case 2: // Anti-Lambda
        arrayV0s = fArrV0sALambdaFiltered;
        profV0sPos = fV0sDiffTwoPos_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        profV0sNeg = fV0sDiffTwoNeg_Lambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hPtInvMass = fV0sLambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
        hInvMass = fV0sInvMassLambda[fCentBinIndex][iHarmonicsIndex][iEtaGapIndex];
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

    if(!profV0sNeg && dEtaGap >= 0.0) // if with eta gap, negative profile has to be valid as well
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
    else if(dEtaGap >= 0.) // with eta gap
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
    //::Info("EstimateRefCumulant","Estimating with NO eta gap");
    dWeight = fCountRefPos * (fCountRefPos - 1);
    dAmp = (fVecRefPos * (TComplex::Conjugate(fVecRefPos))).Re();
    dVal = (dAmp - fCountRefPos) / dWeight;
    if( TMath::Abs(dVal < 1) && (dWeight > 1) && (fCentPercentile > 0) )
      profile->Fill(fCentPercentile, dVal, dWeight);
  }
  else // with eta gap 
  {
    //::Info("EstimateRefCumulant","Estimating with eta gap");
    dWeight = fCountRefPos * fCountRefNeg;
    dAmp = (fVecRefPos * (TComplex::Conjugate(fVecRefNeg))).Re();
    dVal = (dAmp) / dWeight;
    if( TMath::Abs(dVal < 1) && (dWeight > 0) && (fCentPercentile > 0) )
      profile->Fill(fCentPercentile, dVal, dWeight); // ! Fill always just VALUE and WEIGHT separately (not like value*weight) ->see testProfile
  }

  return;
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
void AliAnalysisTaskFlowPID::FilterTracks()
{
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

    if(IsTrackSelected(track))
    {  
      new(fArrTracksFiltered[iNumSelected]) AliAODTrack(*track);
      iNumSelected++;
    }
  }
  
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

  if(fPIDResponse)
  {
   //AliTPCPIDResponse responseTPC = (AliTPCPIDResponse) fPIDResponse->GetTPCResponse();
  }

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