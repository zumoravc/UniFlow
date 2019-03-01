#include "FlowTask.h"
#include "ProcessUniFlow.h"

#include <vector>
#include "TString.h"
#include "TList.h"
#include "TH1.h"

ClassImp(FlowTask);

//_____________________________________________________________________________
FlowTask::FlowTask(PartSpecies species, const char* name) :
  fHarmonics{0},
  fEtaGap{-1.0},
  fNumPtBins{-1},
  fPtBinsEdges{},
  fNumSamples{10},
  fConsCorr{kFALSE},
  fShowMult{kFALSE},
  fRebinning{kTRUE},
  fSampleMerging{kFALSE},
  fDesampleUseRMS{kFALSE},
  fCumOrderMax{kNon},
  fDoCorrMixed{kFALSE},
  fMergePosNeg{kFALSE},
  fFlowFitPhiSubtLS{kFALSE},
  fRebinFlowMass{0},
  fRebinInvMass{0},
  fFlowFitRangeLow{-1.0},
  fFlowFitRangeHigh{-1.0},
  fNumParMassSig{0},
  fNumParMassBG{0},
  fNumParFlowBG{0},
  fVecHistInvMass{},
  fVecHistInvMassBG{},
  fVecHistFlowMass{},
  fVecHistFlowMassFour{},
  fFlowFitMassSig{},
  fFlowFitMassBG{},
  fFlowFitFlowBG{},
  fListProfiles{nullptr},
  fListHistos{nullptr},
  fName{name},
  fTaskTag{},
  fSpecies{species},
  fInputTag{},
  fMixedDiff{},
  fMixedRefs{}
{
  fListProfiles = new TList();
  fListProfiles->SetOwner(kTRUE);
  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);

  fTaskTag = ProcessUniFlow::GetSpeciesName(species);
  if(!fName.EqualTo("")) fTaskTag.Append(Form("_%s",fName.Data()));
}
//_____________________________________________________________________________
FlowTask::~FlowTask()
{
  if(fListProfiles) { fListProfiles->Clear(); delete fListProfiles; }
  if(fListHistos) { fListHistos->Clear(); delete fListHistos; }
  if(fVecHistFlowMass) { delete fVecHistFlowMass; }
  if(fVecHistFlowMassFour) { delete fVecHistFlowMassFour; }
  if(fVecHistInvMass) { delete fVecHistInvMass; }
  if(fVecHistInvMassBG) { delete fVecHistInvMassBG; }
}
//_____________________________________________________________________________
void FlowTask::SetFitParDefaults(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { ProcessUniFlow::Error("Wrong size of parameters array.","SetFitParDefaults"); return; }
  if(!array) { ProcessUniFlow::Error("Wrong array.","SetFitParDefaults"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParDefaults[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::SetFitParLimitsLow(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { ProcessUniFlow::Error("Wrong size of parameters array.","SetFitParLimitsLow"); return; }
  if(!array) { ProcessUniFlow::Error("Wrong array.","SetFitParLimitsLow"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParLimLow[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::SetFitParLimitsHigh(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { ProcessUniFlow::Error("Wrong size of parameters array.","SetFitParLimitsHigh"); return; }
  if(!array) { ProcessUniFlow::Error("Wrong array.","SetFitParLimitsHigh"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParLimHigh[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::PrintTask()
{
  printf("----- Printing task info ------\n");
  printf("   fName: \"%s\"\n",fName.Data());
  printf("   fTaskTag: \"%s\"\n",fTaskTag.Data());
  printf("   fInputTag: \"%s\"\n",fInputTag.Data());
  printf("   fSpecies: %s (%d)\n",ProcessUniFlow::GetSpeciesName(fSpecies).Data(),fSpecies);
  printf("   fHarmonics: %d\n",fHarmonics);
  printf("   fEtaGap: %g\n",fEtaGap);
  printf("   fNumPtBins: %d\n",fNumPtBins);
  if(fNumPtBins > 1) { printf("   fPtBinsEdges: "); for(Int_t i(0); i < (Int_t) fPtBinsEdges.size() ; ++i) printf("%g ",fPtBinsEdges[i]); printf("\n"); }
  printf("   fRebinning: %s\n", fRebinning ? "true" : "false");
  printf("   fMergePosNeg: %s\n", fMergePosNeg ? "true" : "false");
  printf("   fNumSamples: %d\n", fNumSamples);
  printf("   fDesampleUseRMS: %s\n", fDesampleUseRMS ? "true" : "false");
  printf("   fSampleMerging: %s\n", fSampleMerging ? "true" : "false");
  printf("   fConsCorr: %s\n", fConsCorr ? "true" : "false");
  printf("   fShowMult: %s\n", fShowMult ? "true" : "false");
  printf("   fCumOrderMax: %d\n", fCumOrderMax);
  printf("   fDoCorrMixed: %s\n", fDoCorrMixed ? "true" : "false");
  printf("   fMixedRefs: %s\n",fMixedRefs.Data());
  printf("   fMixedDiff: %s\n",fMixedDiff.Data());
  printf("   fFlowFitRangeLow: %g\n",fFlowFitRangeLow);
  printf("   fFlowFitRangeHigh: %g\n",fFlowFitRangeHigh);
  printf("   fFlowFitPhiSubtLS: %s\n", fFlowFitPhiSubtLS ? "true" : "false");
  printf("   fRebinFlowMass: %d\n",fRebinFlowMass);
  printf("   fRebinInvMass: %d\n",fRebinInvMass);
  printf("------------------------------\n");

  return;
}
