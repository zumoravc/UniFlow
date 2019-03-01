#ifndef FLOWTASK_H
#define FLOWTASK_H

#include <vector>
#include "TString.h"
#include "ProcessUniFlow.h"

class TList;
class TCanvas;
class TH1D;

class FlowTask
{
  public:
                FlowTask(PartSpecies species = kUnknown, const char* name = "");
                ~FlowTask(); // default destructor

    void        PrintTask(); // listing values of internal properties

    TString     GetEtaGapString() { return TString(Form("%02.2g",10*fEtaGap)); } // used for "character-safe" names

    void        SetHarmonics(Int_t harm) { fHarmonics = harm; }
    void        SetEtaGap(Float_t eta) { fEtaGap = eta; }
    void        SetNumSamples(Short_t num) { fNumSamples = num; }
    void        SetInputTag(const char* name) { fInputTag = name; }
    void        SetPtBins(std::vector<Double_t> array) { fPtBinsEdges = array; fNumPtBins = (Int_t) array.size() - 1; } // setup the pt binning for this task using std::vectors. NB: possible with {}
    void        SetShowMultDist(Bool_t show) { fShowMult = show; }
    void        SetConsiderCorrelations(Bool_t cor = kTRUE) { fConsCorr = cor; }
    void        SetRebinning(Bool_t rebin = kTRUE) { fRebinning = rebin; }
    void        SetMergePosNeg(Bool_t merge = kTRUE) { fMergePosNeg = merge; }
    void        SetDesamplingUseRMS(Bool_t use = kTRUE) { fDesampleUseRMS = use; }
    void        DoCorrMixed(TString nameDiff, TString nameRefs) { fDoCorrMixed = kTRUE; fMixedDiff = nameDiff; fMixedRefs = nameRefs; }
    void        DoCumTwo(Bool_t use = kTRUE) { fDoCumTwo = use; }
    void        DoCumFour(Bool_t use = kTRUE) { fDoCumFour = use; }

    // fitting
    void        SetInvMassRebin(Short_t rebin = 2) { fRebinInvMass = rebin; }
    void        SetFlowMassRebin(Short_t rebin = 2) { fRebinFlowMass = rebin; }
    // fitting parameters
    void        SetFitPhiSubtLS(Bool_t sub = kTRUE) { fFlowFitPhiSubtLS = sub; }
    void        SetFitMassRange(Double_t dLow, Double_t dHigh) { fFlowFitRangeLow = dLow; fFlowFitRangeHigh = dHigh; }
    void        SetFitMassSig(TString func, Int_t pars) { fFlowFitMassSig = func; fNumParMassSig = pars; }
    void        SetFitMassBG(TString func, Int_t pars) { fFlowFitMassBG = func; fNumParMassBG = pars; }
    void        SetFitFlowBG(TString func, Int_t pars) { fFlowFitFlowBG = func; fNumParFlowBG = pars; }
    void        SetFitParDefaults(Double_t* array, Int_t size);
    void        SetFitParLimitsLow(Double_t* array, Int_t size);
    void        SetFitParLimitsHigh(Double_t* array, Int_t size);

    TString     fTaskTag; // "unique" tag used primarily for storing output
    TString     fName; // task name
    PartSpecies fSpecies; // species involved
    TString     fInputTag; // alterinative tag appended to name of input histos & profiles
    Int_t       fHarmonics; // harmonics
    Double_t    fEtaGap; // eta gap
    Bool_t      fConsCorr; // consider correlations in cumulant / flow calculations
    Int_t       fNumSamples; // [10] number of samples
    Int_t       fNumPtBins; // actual number of pT bins (not size of array) for rebinning
    std::vector<Double_t>   fPtBinsEdges; // pt binning
    Bool_t      fShowMult; // show multiplicity distribution
    Bool_t      fSampleMerging; // [kFALSE] flag for merging TProfiles (good for refs)
    Bool_t      fRebinning; // [kTRUE] flag for rebinning prior to desampling
    Bool_t      fDesampleUseRMS; // [kFALSE] flag for using RMS as uncertainty during desampling
    Bool_t      fMergePosNeg; // [kFALSE] flag for merging results corresponding to positive and negative POIs
    Bool_t      fDoCumTwo; // [kTRUE] flag for processing 2-part. cumulants
    Bool_t      fDoCumFour; // [kFALSE] flag for processing 4-part. cumulants
    Bool_t      fDoCorrMixed; // [kFALSE] flag for processing mixed harmonics (non-linear flow modes)
    TString     fMixedDiff; // name (tag) for diff. profile for mixed harmonics
    TString     fMixedRefs; // name (tag) for reference profile for mixed harmonics
    // Reconstructed fitting
    Bool_t      fFlowFitPhiSubtLS; // [kFALSE] flag for subtraction of like-sign background from the unlike-sign one
    Int_t       fRebinInvMass; // flag for rebinning inv-mass (and BG) histo
    Int_t       fRebinFlowMass; // flag for rebinning flow-mass profile

    Double_t    fFlowFitRangeLow; // lower edge for fitting during flow extraction
    Double_t    fFlowFitRangeHigh; // high edge for fitting during flow extraction
    TString     fFlowFitMassSig; // formula for signal contribution of inv.mass fit
    TString     fFlowFitMassBG; // formula for bg contribution of inv.mass fit
    TString     fFlowFitFlowBG; // formula for bg contribution of flow-mass fit
    Int_t       fNumParMassSig; // number of parameters in 'fFlowFitMassSig'
    Int_t       fNumParMassBG; // number of parameters in 'fFlowFitMassBG'
    Int_t       fNumParFlowBG; // number of parameters in 'fFlowFitFlowBG'
    static const Int_t fNumParsMax = 20;
    Double_t    fFitParDefaults[fNumParsMax]; // default values for all parameters in all formulas (i.e. master flow-mass formula)
    Double_t    fFitParLimLow[fNumParsMax]; // low limit values for all parameters in all formulas (i.e. master flow-mass formula)
    Double_t    fFitParLimHigh[fNumParsMax]; // high limit values for all parameters in all formulas (i.e. master flow-mass formula)
    TList*      fListProfiles; // TList with profiles slices
    TList*      fListHistos; // TList with histos slices


    std::vector<TH1D*>* fVecHistInvMass; // container for sliced inv. mass projections
    std::vector<TH1D*>* fVecHistInvMassBG; // container for sliced inv. mass projections for BG candidates (phi)
    std::vector<TH1D*>* fVecHistFlowMass; // container for sliced flow-mass projections
    std::vector<TH1D*>* fVecHistFlowMassFour; // container for sliced flow-mass projections

  protected:
  private:

    ClassDef(FlowTask,1);
};

#endif
