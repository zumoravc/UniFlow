#include "FlowTask.h"

#include <vector>
#include "TROOT.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

ClassImp(FlowTask);

//_____________________________________________________________________________
FlowTask::FlowTask(PartSpecies species, const char* name) :
  fHarmonics(0),
  fEtaGap(0),
  fNumPtBins(-1),
  fPtBinsEdges()
{
  fName = name;
  fSpecies = species;
  fInputTag = "";
  fNumSamples = 10;
  fDoFour = kFALSE;
  fConsCorr = kFALSE;
  fShowMult = kFALSE;
  fRebinning = kTRUE;
  fSampleMerging = kFALSE;
  fDesampleUseRMS = kFALSE;
  fProcessMixed = kFALSE;
  fMixedDiff = "";
  fMixedRefs = "";
  fMergePosNeg = kFALSE;
  fFlowFitPhiSubtLS = kFALSE;
  fRebinFlowMass = 0;
  fRebinInvMass = 0;
  fFlowFitRangeLow = -1.0;
  fFlowFitRangeHigh = -1.0;
  fFlowFitMassSig = TString();
  fFlowFitMassBG = TString();
  fFlowFitFlowBG = TString();
  fNumParMassSig = 0;
  fNumParMassBG = 0;
  fNumParFlowBG = 0;
  fListProfiles = new TList();
  fListProfiles->SetOwner(kTRUE);
  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);
  fVecHistInvMass = new std::vector<TH1D*>;
  fVecHistInvMassBG = new std::vector<TH1D*>;
  fVecHistFlowMass = new std::vector<TH1D*>;
  fVecHistFlowMassFour = new std::vector<TH1D*>;

  fTaskTag = this->GetSpeciesName();
  if(!fName.EqualTo("")) fTaskTag.Append(Form("_%s",fName.Data()));
}
//_____________________________________________________________________________
FlowTask::~FlowTask()
{
  if(fListProfiles) { fListProfiles->Clear(); delete fListProfiles; }
  if(fListHistos) { fListHistos->Clear(); delete fListHistos; }
  if(fVecHistFlowMass) delete fVecHistFlowMass;
  if(fVecHistFlowMassFour) delete fVecHistFlowMassFour;
  if(fVecHistInvMass) delete fVecHistInvMass;
  if(fVecHistInvMassBG) delete fVecHistInvMassBG;
  if(fCanvas) delete fCanvas;
}
//_____________________________________________________________________________
void FlowTask::SetFitParDefaults(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { Error("Wrong size of parameters array.","SetFitParDefaults"); return; }
  if(!array) { Error("Wrong array.","SetFitParDefaults"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParDefaults[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::SetFitParLimitsLow(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { Error("Wrong size of parameters array.","SetFitParLimitsLow"); return; }
  if(!array) { Error("Wrong array.","SetFitParLimitsLow"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParLimLow[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::SetFitParLimitsHigh(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { Error("Wrong size of parameters array.","SetFitParLimitsHigh"); return; }
  if(!array) { Error("Wrong array.","SetFitParLimitsHigh"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParLimHigh[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
TString FlowTask::GetSpeciesName()
{
  TString name = TString();
  switch (fSpecies)
  {
    case kRefs : name.Append("Refs"); break;
    case kCharged : name.Append("Charged"); break;
    case kPion : name.Append("Pion"); break;
    case kKaon : name.Append("Kaon"); break;
    case kProton : name.Append("Proton"); break;
    case kPhi : name.Append("Phi"); break;
    case kK0s : name.Append("K0s"); break;
    case kLambda : name.Append("Lambda"); break;
    default: name.Append("Unknown");
  }
  return name;
}
//_____________________________________________________________________________
TString FlowTask::GetSpeciesLabel()
{
  TString label = TString();
  switch (fSpecies)
  {
    case kRefs : label.Append("RFP"); break;
    case kCharged : label.Append("h^{#pm}"); break;
    case kPion : label.Append("#pi^{#pm}"); break;
    case kKaon : label.Append("K^{#pm}"); break;
    case kProton : label.Append("p/{#bar{p}}"); break;
    case kPhi : label.Append("#phi"); break;
    case kK0s : label.Append("K_{S}^{0}"); break;
    case kLambda : label.Append("#Lambda/#bar{#Lambda}"); break;
    default: label.Append("Non");
  }
  return label;
}
//_____________________________________________________________________________
void FlowTask::PrintTask()
{
  printf("----- Printing task info ------\n");
  printf("   fTaskTag: \"%s\"\n",fTaskTag.Data());
  printf("   fName: \"%s\"\n",fName.Data());
  printf("   fInputTag: \"%s\"\n",fInputTag.Data());
  printf("   fSpecies: %s (%d)\n",GetSpeciesName().Data(),fSpecies);
  printf("   fHarmonics: %d\n",fHarmonics);
  printf("   fEtaGap: %g\n",fEtaGap);
  printf("   fProcessMixed: %s\n", fProcessMixed ? "true" : "false");
  printf("   fDoFour: %s\n", fDoFour ? "true" : "false");
  printf("   fShowMult: %s\n", fShowMult ? "true" : "false");
  printf("   fMergePosNeg: %s\n", fMergePosNeg ? "true" : "false");
  printf("   fFlowFitRangeLow: %g\n",fFlowFitRangeLow);
  printf("   fFlowFitRangeHigh: %g\n",fFlowFitRangeHigh);
  printf("   fNumPtBins: %d\n",fNumPtBins);
  if(fNumPtBins > 1) { printf("   fPtBinsEdges: "); for(Int_t i(0); i < (Int_t) fPtBinsEdges.size() ; ++i) printf("%g ",fPtBinsEdges[i]); printf("\n"); }
  printf("------------------------------\n");
  return;
}
