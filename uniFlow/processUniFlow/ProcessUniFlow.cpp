/* ProcessUniFlow class
 *
 * Class implemented for processing results of AliAnalysisTaskUniFlow task.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */
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

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();


class FlowTask
{
  public:
    enum    PartSpecies {kUnknown=0, kRefs, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi}; // list of all particle species of interest

                FlowTask(PartSpecies species = kUnknown, const char* name = "");
                ~FlowTask(); // default destructor

    void        PrintTask(); // listing values of internal properties

    TString     GetSpeciesName();
    TString     GetEtaGapString() { return TString(Form("%02.2g",10*fEtaGap)); }

    void        SetHarmonics(Int_t harm) { fHarmonics = harm; }
    void        SetEtaGap(Float_t eta) { fEtaGap = eta; }
    void        SetNumSamples(Short_t num) { fNumSamples = num; }
    void        SetInputTag(const char* name) { fInputTag = name; }
    void        SetPtBins(Double_t* array, const Short_t size); // setup the pt binning for this task, where size is number of elements in array
    void        SetShowMultDist(Bool_t show) { fShowMult = show; }
    void        SetRebinning(Bool_t rebin = kTRUE) { fRebinning = rebin; }
    void        SetMergePosNeg(Bool_t merge = kTRUE) {fMergePosNeg = merge; }
    void        SetDesamplingUseRMS(Bool_t use = kTRUE) { fDesampleUseRMS = use; }
    void        SuggestPtBinning(Bool_t bin = kTRUE, Double_t entries = 20000) { fSuggestPtBins = bin; fSuggestPtBinEntries = entries; } // suggest pt binning based on number of candidates

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

  protected:
  private:

    TString     fTaskTag; // "unique" tag used primarily for storing output
    TString     fName; // task name
    PartSpecies fSpecies; // species involved
    TString     fInputTag; // alterinative tag appended to name of input histos & profiles
    Int_t       fHarmonics; // harmonics
    Double_t    fEtaGap; // eta gap
    Short_t     fNumSamples; // [10] number of samples
    static const Short_t fNumPtBinsMax = 100; // initialization (maximum) number of pt bins
    Double_t    fPtBinsEdges[fNumPtBinsMax]; // pt binning
    Short_t     fNumPtBins; // actual number of pT bins (not size of array) for rebinning
    Bool_t      fShowMult; // show multiplicity distribution
    Bool_t      fSampleMerging; // [kFALSE] flag for merging TProfiles (good for refs)
    Bool_t      fRebinning; // [kTRUE] flag for rebinning prior to desampling
    Bool_t      fDesampleUseRMS; // [kFALSE] flag for using RMS as uncertainty during desampling
    Bool_t      fMergePosNeg; // [kFALSE] flag for merging results corresponding to positive and negative POIs
    Bool_t      fSuggestPtBins; // suggest pt binning
    Double_t    fSuggestPtBinEntries; // suggest pt binning
    // Reconstructed fitting
    Bool_t      fFlowFitPhiSubtLS; // [kFALSE] flag for subtraction of like-sign background from the unlike-sign one
    Short_t     fRebinInvMass; // flag for rebinning inv-mass (and BG) histo
    Short_t     fRebinFlowMass; // flag for rebinning flow-mass profile

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

    std::vector<TH1D*>* fVecHistInvMass; // container for sliced inv. mass projections
    std::vector<TH1D*>* fVecHistInvMassBG; // container for sliced inv. mass projections for BG candidates (phi)
    std::vector<TH1D*>* fVecHistFlowMass; // container for sliced flow-mass projections
    TCanvas*     fCanvas; // temporary canvas for mass plotting
};

//_____________________________________________________________________________
FlowTask::FlowTask(PartSpecies species, const char* name)
{
  fName = name;
  fSpecies = species;
  fHarmonics = 0;
  fEtaGap = 0;
  fInputTag = "";
  fNumSamples = 10;
  fNumPtBins = -1;
  fShowMult = kFALSE;
  fRebinning = kTRUE;
  fSampleMerging = kFALSE;
  fDesampleUseRMS = kFALSE;
  fMergePosNeg = kFALSE;
  fSuggestPtBins = kFALSE;
  fSuggestPtBinEntries = 20000;
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
  fVecHistInvMass = new std::vector<TH1D*>;
  fVecHistInvMassBG = new std::vector<TH1D*>;
  fVecHistFlowMass = new std::vector<TH1D*>;

  fTaskTag = this->GetSpeciesName();
  if(!fName.EqualTo("")) fTaskTag.Append(Form("_%s",fName.Data()));
}
//_____________________________________________________________________________
FlowTask::~FlowTask()
{
  if(fVecHistFlowMass) delete fVecHistFlowMass;
  if(fVecHistInvMass) delete fVecHistInvMass;
  if(fVecHistInvMassBG) delete fVecHistInvMassBG;
  if(fCanvas) delete fCanvas;
}
//_____________________________________________________________________________
void FlowTask::SetPtBins(Double_t* array, const Short_t size)
{
  if(size < 0 || size > fNumPtBinsMax) { Error("Wrong size of pt binning array.","SetPtBins"); return; }
  if(!array) { Error("Wrong array.","SetPtBins"); return; }

  fNumPtBins = size - 1;

  for(Short_t i(0); i < size; i++)
  {
    fPtBinsEdges[i] = array[i];
  }

  return;
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
void FlowTask::PrintTask()
{
  printf("----- Printing task info ------\n");
  printf("   fTaskTag: \"%s\"\n",fTaskTag.Data());
  printf("   fName: \"%s\"\n",fName.Data());
  printf("   fInputTag: \"%s\"\n",fInputTag.Data());
  printf("   fSpecies: %s (%d)\n",GetSpeciesName().Data(),fSpecies);
  printf("   fHarmonics: %d\n",fHarmonics);
  printf("   fEtaGap: %g\n",fEtaGap);
  printf("   fShowMult: %s\n", fShowMult ? "true" : "false");
  printf("   fSuggestPtBins: %s\n", fSuggestPtBins ? "true" : "false");
  printf("   fMergePosNeg: %s\n", fMergePosNeg ? "true" : "false");
  printf("   fFlowFitRangeLow: %g\n",fFlowFitRangeLow);
  printf("   fFlowFitRangeHigh: %g\n",fFlowFitRangeHigh);
  printf("   fNumPtBins: %d (limit %d)\n",fNumPtBins,fNumPtBinsMax);
  if(fNumPtBins > -1) { printf("   fPtBinsEdges: "); for(Short_t i(0); i < fNumPtBins+1; i++) printf("%g ",fPtBinsEdges[i]); printf("\n"); }
  printf("------------------------------\n");
  return;
}
//#############################################################################
//----------------- End of FlowTask class -------------------------------------
//#############################################################################

class ProcessUniFlow
{
  public:

                ProcessUniFlow(); // default constructor
                ~ProcessUniFlow(); // default destructor

    void        SetInputFilePath(const char* path) { fsInputFilePath = path; }
    void        SetInputFileName(const char* name) { fsInputFileName = name; }
    void        SetOutputFilePath(const char* path) { fsOutputFilePathRoot = path; }
    void        SetOutputFileName(const char* name) { fsOutputFileName = name; }
    void        SetOutputFileMode(const char* mode = "RECREATE") { fsOutputFileMode = mode; }
    void        SetTaskName(const char* name) { fsTaskName = name; }
    void        SetGlobalProfNameLabel(const char* label = "") { fsGlobalProfNameLabel = label; } // add global profile label for all tasks NOTE: for the purpose of Flow sub
    void        SetSaveMult(Bool_t bSave = kTRUE) { fbSaveMult = bSave; } // save reference multiplicity
    void        SetMultiplicityBins(Double_t* array, const Short_t size); // setup the global multiplicity binning, where size is number of elements in array
    void        SetFitCumulants(Bool_t cum = kTRUE) { fFlowFitCumulants = cum; } // use cn{2} vs m_inv instead of vn{2} vs. m_inv
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }
    void        AddTask(FlowTask* task = 0x0); // add task to internal lists of all tasks
    void        Run(); // running the task (main body of the class)
    void        Clear(); // clearing (removing tasks, etc.) after running
  protected:

  private:
    Bool_t      Initialize(); // initialization task
    Bool_t      LoadLists(); // loading flow lists from input file

    Bool_t      ProcessTask(FlowTask* task = 0x0); // process FlowTask according to it setting
    Bool_t      ProcessRefs(FlowTask* task = 0x0); // process reference flow task
    Bool_t      ProcessDirect(FlowTask* task = 0x0, Short_t iMultBin = 0); // process PID (pion,kaon,proton) flow task
    Bool_t      ProcessReconstructed(FlowTask* task = 0x0, Short_t iMultBin = 0); // process  V0s flow
    Bool_t      PrepareSlices(const Short_t multBin, FlowTask* task = 0x0, TProfile3D* p3Cor = 0x0, TH3D* h3Entries = 0x0, TH3D* h3EntriesBG = 0x0); // prepare

    Bool_t 	    ExtractFlowOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowPhiOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowK0sOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowLambdaOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for Lambda candidates


    void        SuggestMultBinning(const Short_t numFractions);
    void        SuggestPtBinning(TH3D* histEntries = 0x0, TProfile3D* profFlowOrig = 0x0, FlowTask* task = 0x0, Short_t binMult = 0); //
    TH1D*       DesampleList(TList* list = 0x0, FlowTask* task = 0x0, Short_t iMultBin = 0); // Desample list of samples for estimating the uncertanity
    TH1D*       TestRebin(TH1D* hOrig = 0x0, FlowTask* task = 0x0); // testing desample - manual rebin

    void        TestProjections(); // testing projection of reconstructed particles
    TProfile2D* Project3DProfile(const TProfile3D* prof3dorig = 0x0); // making projection out of TProfile3D
    TProfile2D* DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,bool originalRange, bool useUF, bool useOF) const;
    TH2D*       DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY, bool computeErrors, bool originalRange, bool useUF, bool useOF) const;


    // printing output methods
    void        Fatal(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Error(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Warning(TString sMsg, TString sMethod = ""); // printf the msg as warning
    void        Info(TString sMsg, TString sMethod = ""); // printf the msg as info
    void        Debug(TString sMsg, TString sMethod = ""); // printf the msg as info

    static const Short_t fiNumMultBinsGlobal = 200; // global initialization number of bins
    Double_t    fdMultBins[fiNumMultBinsGlobal]; // global multiplicity/centrality binning
    Short_t     fiNumMultBins; // number of multiplicity bins (not size of array)

    TString     fsInputFilePath; // path to the input folder with input file
    TString     fsInputFileName; // name of input file
    TString     fsOutputFilePath; // current (working) path to the output folder (process/task dependent)
    TString     fsOutputFilePathRoot; // root (global/top-level) path to the output folder
    TString     fsOutputFileName; // name of output file
    TString     fsOutputFileMode; // [RECREATE] mode of output file
    TString     fsTaskName; // name of task (inchluded in data structure names)
    TString     fsOutputFileFormat; // [pdf] format of output files (pictures)
    TString     fsGlobalProfNameLabel; // global profile label for all task
    Bool_t      fbSaveMult; // [kFALSE]
    Bool_t      fFlowFitCumulants; // [kFALSE]

    Bool_t      fbInit; // flag for initialization status
    Bool_t      fbDebug; // flag for debugging : if kTRUE Debug() messages are displayed
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! output file container
    TFile*      ffDesampleFile; //! output file for results of desampling
    TFile*      ffFitsFile; //! output file for fitting procedure
    TList*      flFlowRefs; //! TList from input file with RFPs flow profiles
    TList*      flFlowCharged; //! TList from input file with Charged flow profiles
    TList*      flFlowPID; //! TList from input file with PID (pi,K,p) flow profiles
    TList*      flFlowPhi; //! TList from input file with Phi flow profiles
    TList*      flFlowK0s; //! TList from input file with K0s flow profiles
    TList*      flFlowLambda; //! TList from input file with Lambda flow profiles
    TList*      flQACharged; //! TList from input file with Charged QA plots / profiles
    TList*      flQAPID; //! TList from input file with PID (pi,K,p) QA plots / profiles
    TList*      flQAPhi; //! TList from input file with Phi QA plots / profiles
    TList*      flQAV0s; //! TList from input file with K0s QA plots / profiles
    std::vector<FlowTask*> fvTasks; // vector of task for individual species proccesing

};
//_____________________________________________________________________________
ProcessUniFlow::ProcessUniFlow() :
  fbDebug(kFALSE),
  fbInit(kFALSE),
  fbSaveMult(kFALSE),
  fFlowFitCumulants(kFALSE),
  ffInputFile(0x0),
  ffOutputFile(0x0),
  ffFitsFile(0x0),
  ffDesampleFile(0x0),
  flFlowRefs(0x0),
  flFlowCharged(0x0),
  flFlowPID(0x0),
  flFlowPhi(0x0),
  flFlowK0s(0x0),
  flFlowLambda(0x0),
  flQACharged(0x0),
  flQAPID(0x0),
  flQAPhi(0x0),
  flQAV0s(0x0)
{
  // default constructor
  fsInputFilePath = TString("");
  fsInputFileName = TString("AnalysisResults.root");
  fsOutputFilePathRoot = TString("");
  fsOutputFilePath = TString("");
  fsOutputFileName = TString("UniFlow.root");
  fsOutputFileMode = TString("RECREATE");
  fsTaskName = TString("UniFlow");
  fsOutputFileFormat = TString("pdf");
  fsGlobalProfNameLabel = TString("");
  fvTasks = std::vector<FlowTask*>();

  for(Short_t i(0); i < fiNumMultBinsGlobal; i++) fdMultBins[i] = 0;
}
//_____________________________________________________________________________
ProcessUniFlow::~ProcessUniFlow()
{
  // default destructor
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;
  if(ffFitsFile) delete ffFitsFile;

  if(flFlowRefs) delete flFlowRefs;
  if(flFlowCharged) delete flFlowCharged;
  if(flFlowPID) delete flFlowPID;
  if(flFlowPhi) delete flFlowPhi;
  if(flFlowK0s) delete flFlowK0s;
  if(flFlowLambda) delete flFlowLambda;

  // deleting the FlowTasks
  const Short_t iNumTasks = fvTasks.size();
  for(Short_t index(0); index < iNumTasks; index++)
  {
    delete fvTasks.at(index);
  }
}
//_____________________________________________________________________________
void ProcessUniFlow::Clear()
{
  Info("Cleaning ProcessUniFlow instance","Clear");
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;

  flFlowRefs = 0x0;
  flFlowCharged = 0x0;
  flFlowPID = 0x0;
  flFlowPhi = 0x0;
  flFlowK0s = 0x0;
  flFlowLambda = 0x0;

  const Short_t iNumTasks = fvTasks.size();
  for(Short_t index(0); index < iNumTasks; index++)
  {
    if(fvTasks.at(index)) delete fvTasks.at(index);
  }
  Info("Cleaning done!","Clear");
}
//_____________________________________________________________________________
void ProcessUniFlow::Run()
{
  gStyle->SetOptFit(1100);

  // main body of the class
  if(!Initialize()) { Fatal("Task not initialized","Run"); return; }

  const Short_t iNumTasks = fvTasks.size();

  Info("===== Running over tasks ======","Run");
  Info(Form("  Number of tasks: %d\n",iNumTasks),"Run");
  for(Short_t iTask(0); iTask < iNumTasks; iTask++)
  {
    FlowTask* currentTask = fvTasks.at(iTask);
    if(!currentTask) continue;
    if(!ProcessTask(currentTask)) {return;}
  }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Initialize()
{
  // initialization of all necessery prerequisits
  Info("Initializating task","Initialize");

  // opening input file
  ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(),fsInputFileName.Data()),"READ");
  if(!ffInputFile || !ffInputFile->IsOpen())
  {
    Fatal(Form("Input file %s/%s not open",fsInputFilePath.Data(),fsInputFileName.Data()),"Initialize");
    return kFALSE;
  }

  // Setting current output path
  fsOutputFilePath = fsOutputFilePathRoot;

  // checking specified output folder & required sub-folders
  gSystem->mkdir(fsOutputFilePath.Data(),kTRUE);

  // opening output file
  ffOutputFile = TFile::Open(Form("%s/%s",fsOutputFilePath.Data(),fsOutputFileName.Data()),fsOutputFileMode.Data());
  if(!ffOutputFile || !ffOutputFile->IsOpen())
  {
    Fatal(Form("Output file %s/%s not open",fsOutputFilePath.Data(),fsOutputFileName.Data()),"Initialize");
    return kFALSE;
  }

  // creating output file for Desampling
  ffDesampleFile = TFile::Open(Form("%s/desampling.root",fsOutputFilePath.Data()),"RECREATE");
  if(!ffDesampleFile) { Fatal(Form("Output desampling file '%s/desampling.root' not open!","Initialize")); return kFALSE; }

  // creating output file for fits
  ffFitsFile = TFile::Open(Form("%s/fits.root",fsOutputFilePath.Data()),"RECREATE");
  if(!ffFitsFile) { Fatal(Form("Output desampling file '%s/fits.root' not open!","Initialize")); return kFALSE; }

  Info("Files loaded","Initialize");

  if(!LoadLists()) return kFALSE;
  Info("Flow lists loaded","Initialize");

  // initialization succesfull
  Info("Initialization succesfull","Initialize");
  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::SetMultiplicityBins(Double_t* array, const Short_t size)
{
  if(size < 0 || size > fiNumMultBinsGlobal) { Fatal("Wrong number of bins selected","SetMultiplicityBins"); return; }
  if(!array) { Fatal("Multiplicity array not valid","SetMultiplicityBins"); return; }

  fiNumMultBins = size-1;

  for(Short_t i(0); i < size; i++)
  {
    fdMultBins[i] = array[i];
  }
  return;
}
//_____________________________________________________________________________
void ProcessUniFlow::SuggestMultBinning(const Short_t numFractions)
{
  if(numFractions < 1) { Error("Suggested number of fractions is too low","SuggestMultBinning"); return; }
  if(!fbInit) Initialize();

  // loading multiplicity dist.
  ffInputFile->cd(fsTaskName.Data());
  gDirectory->ls();

  TList* lQACharged = (TList*) gDirectory->Get(Form("QA_Charged_%s",fsTaskName.Data()));
  if(!lQACharged) return;

  lQACharged->ls();
  TH1D* hMult = (TH1D*) lQACharged->FindObject("fhQAChargedMult_After");

  TCanvas* canSuggestMult = new TCanvas("canSuggestMult","SuggestMultBinning",400,400);
  canSuggestMult->cd(1);
  gPad->SetLogy();
  hMult->Draw();

  Double_t dCumul = 0;
  Double_t dContent = 0;
  Double_t dWeight = 0;
  Double_t dSum = 0;
  Double_t dSumWeight = 0;

  for(Short_t bin(1); bin < hMult->GetNbinsX()+1; bin++)
  {
    dContent = hMult->GetBinContent(bin);
    dWeight = 1/(bin - 1);
    dSum += dContent*(dWeight); // weighed average
    dSumWeight += dWeight;
  }

  Double_t dEntries = hMult->GetEntries();
  // Double_t dSegment = dEntries / numFractions;
  Double_t dSegment = dSum / dSumWeight / numFractions;

  printf("entries %g | Sum %g | sumweight %g | segmend %g\n",dEntries,dSum,dSumWeight,dSegment);

  for(Short_t bin(1); bin < hMult->GetNbinsX()+1; bin++)
  {
    dContent = hMult->GetBinContent(bin);
    // dWeight = bin - 1;
    dWeight = 1/(bin-1);
    dCumul += dContent*(dWeight);
    // dSumWeight += dWeight;
    if(dCumul >= dSegment)
    {
      printf("bin %d : mult %g (cumul %g)\n",bin,hMult->GetBinLowEdge(bin),dCumul);
      dCumul = 0;
    }
  }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::LoadLists()
{
  // loading TLists into task
  if(!ffInputFile) { Fatal("Input file does not exists!","LoadLists"); return kFALSE; }
  ffInputFile->cd(fsTaskName.Data());

  flFlowRefs = (TList*) gDirectory->Get(Form("Flow_Refs_%s",fsTaskName.Data()));
  if(!flFlowRefs) { Fatal("flFlow_Refs list does not exists!","LoadLists"); ffInputFile->ls(); return kFALSE; }
  flFlowCharged = (TList*) gDirectory->Get(Form("Flow_Charged_%s",fsTaskName.Data()));
  if(!flFlowCharged) { Fatal("flFlow_Charged list does not exists!","LoadLists"); return kFALSE; }
  flFlowPID = (TList*) gDirectory->Get(Form("Flow_PID_%s",fsTaskName.Data()));
  if(!flFlowPID) { Fatal("flFlow_PID list does not exists!","LoadLists"); return kFALSE; }
  flFlowPhi = (TList*) gDirectory->Get(Form("Flow_Phi_%s",fsTaskName.Data()));
  if(!flFlowPhi) { Fatal("flFlow_Phi list does not exists!","LoadLists"); return kFALSE; }
  flFlowK0s = (TList*) gDirectory->Get(Form("Flow_K0s_%s",fsTaskName.Data()));
  if(!flFlowK0s) { Fatal("flFlow_K0s list does not exists!","LoadLists"); return kFALSE; }
  flFlowLambda = (TList*) gDirectory->Get(Form("Flow_Lambda_%s",fsTaskName.Data()));
  if(!flFlowLambda) { Fatal("flFlow_Lambda list does not exists!","LoadLists"); return kFALSE; }

  flQACharged = (TList*) gDirectory->Get(Form("QA_Charged_%s",fsTaskName.Data()));
  if(!flQACharged) { Fatal("flQACharged list does not exists!","LoadLists"); return kFALSE; }
  flQAPID = (TList*) gDirectory->Get(Form("QA_PID_%s",fsTaskName.Data()));
  if(!flQAPID) { Fatal("flQAPID list does not exists!","LoadLists"); return kFALSE; }
  flQAPhi = (TList*) gDirectory->Get(Form("QA_Phi_%s",fsTaskName.Data()));
  if(!flQAPhi) { Fatal("flQAPhi list does not exists!","LoadLists"); return kFALSE; }
  flQAV0s = (TList*) gDirectory->Get(Form("QA_V0s_%s",fsTaskName.Data()));
  if(!flQAV0s) { Fatal("flQAV0s list does not exists!","LoadLists"); return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessTask(FlowTask* task)
{
  fsOutputFilePath = Form("%s/%s",fsOutputFilePathRoot.Data(),task->fTaskTag.Data());
  gSystem->mkdir(fsOutputFilePath.Data(),kTRUE);

  Info(Form("Processing task: %s",task->fName.Data()),"ProcessTask");
  if(!task) { Error("Task not valid!","ProcessTask"); return kFALSE; }

  task->PrintTask();

  switch (task->fSpecies)
  {
    case FlowTask::kRefs:
      if(!ProcessRefs(task)) { Error(Form("Task '%s' not processed correctly!",task->fName.Data()),"ProcessTask"); return kFALSE; }
      break;

    case FlowTask::kCharged:
    case FlowTask::kPion:
    case FlowTask::kKaon:
    case FlowTask::kProton:
      for(Short_t binMult(0); binMult < fiNumMultBins; binMult++) { if(!ProcessDirect(task,binMult)) { Error(Form("Task '%s' not processed correctly!",task->fName.Data()),"ProcessTask"); return kFALSE; }}
      break;

    case FlowTask::kPhi:
    case FlowTask::kK0s:
    case FlowTask::kLambda:
      for(Short_t binMult(0); binMult < fiNumMultBins; binMult++) { if(!ProcessReconstructed(task,binMult)) { Error(Form("Task '%s' not processed correctly!",task->fName.Data()),"ProcessTask"); return kFALSE; } }
      break;
    default: break;
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessRefs(FlowTask* task)
{
  Info("Processing Refs task","ProcesRefs");
  if(!task) { Error("Task not valid!","ProcessRefs"); return kFALSE; }
  if(task->fSpecies != FlowTask::kRefs) { Error("Task species not kRefs!","ProcessRefs"); return kFALSE; }

  // preparing for desampling
  TProfile* prof = 0x0;
  TProfile* profRebin  = 0x0;
  TH1D* histProj = 0x0;
  TList* list = new TList();
  TList* listMerge = new TList();

  // rebinning <multiplicity>
  if(fbSaveMult)
  {
    TProfile* profMult = (TProfile*) flQACharged->FindObject(Form("fpRefsMult"));
    if(!profMult) { Error("MeanMult profile not found!"); flFlowRefs->ls(); return kFALSE; }
    TProfile* profMult_rebin = (TProfile*) profMult->Rebin(fiNumMultBins,Form("%s_rebin",profMult->GetName()),fdMultBins);

    ffOutputFile->cd();
    profMult_rebin->Write(profMult_rebin->GetName());
  }


  for(Short_t i(0); i < task->fNumSamples; i++)
  {
    prof = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_%s<2>_harm%d_gap%02.2g_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,i));
    if(!prof) { Warning(Form("Profile sample %d does not exits. Skipping",i),"ProcesRefs"); continue; }



    if(task->fRebinning)
    {
      // rebinning the profiles
      profRebin = (TProfile*) prof->Rebin(fiNumMultBins,Form("%s_rebin",prof->GetName()),fdMultBins);
      histProj = profRebin->ProjectionX(Form("%s_proj",profRebin->GetName()));
      listMerge->Add(profRebin);
    }
    else
    {
      // no rebinning
      histProj = prof->ProjectionX(Form("%s_proj",prof->GetName()));
      listMerge->Add(prof);
    }
    list->Add(histProj);
  }

  // desampling (similarly to PID & Charged tracks)
  TH1D* hDesampled = DesampleList(list,task);
  if(!hDesampled) { Error("Desampling unsuccesfull","ProcessRefs"); return kFALSE; }
  hDesampled->SetName(Form("hCum2_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));
  hDesampled->SetTitle(Form("%s c_{%d}{2} | Gap %s ",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TH1D* hDesampledFlow = (TH1D*) hDesampled->Clone(Form("%s_flow",hDesampled->GetName()));
  hDesampledFlow->SetName(Form("hFlow2_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));
  hDesampledFlow->SetTitle(Form("%s v_{%d}{2} | Gap %s ",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  // estimating vn out of cn
  Double_t dContent = 0., dError = 0.;
  for(Short_t iBin(1); iBin < hDesampled->GetNbinsX()+1; iBin++)
  {
    dContent = hDesampled->GetBinContent(iBin);
    dError = hDesampled->GetBinError(iBin);

    if(dContent > 0. && dError >= 0.)
    {
      hDesampledFlow->SetBinContent(iBin, TMath::Sqrt(dContent));
      hDesampledFlow->SetBinError(iBin, TMath::Sqrt( TMath::Power(dError,2) / (4*dContent) ));
    }
    else
    {
      hDesampledFlow->SetBinContent(iBin, 9.);
      hDesampledFlow->SetBinError(iBin, 9.);
    }
  }


  // saving to output file
  ffOutputFile->cd();
  hDesampled->Write(Form("%s",hDesampled->GetName()));
  hDesampledFlow->Write(Form("%s",hDesampledFlow->GetName()));

  // merging all samples together (NOTE: good for Refs)
  TH1D* hMerged = 0x0;
  if(task->fSampleMerging)
  {
    TProfile* merged = (TProfile*) listMerge->At(0)->Clone();
    merged->Reset();
    Double_t mergeStatus = merged->Merge(listMerge);
    merged->SetName(Form("hCum2_%s_harm%d_gap%s_merged",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));
    if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessRefs"); return kFALSE; }

    hMerged = (TH1D*) merged->ProjectionX();
    hMerged->SetName(Form("hFlow2_%s_harm%d_gap%s_merged",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

    // getting vn out of cn
    for(Short_t iBin(1); iBin < merged->GetNbinsX()+1; iBin++)
    {
      dContent = merged->GetBinContent(iBin);
      dError = merged->GetBinError(iBin);

      if(dContent > 0. && dError >= 0.)
      {
        hMerged->SetBinContent(iBin, TMath::Sqrt(dContent));
        hMerged->SetBinError(iBin, TMath::Sqrt( TMath::Power(dError,2) / (4*dContent) ));
      }
      else
      {
        hMerged->SetBinContent(iBin, 9.);
        hMerged->SetBinError(iBin, 9.);
      }
    }
    merged->Write(Form("%s",merged->GetName()));
    hMerged->Write(Form("%s",hMerged->GetName()));
  }



  if(!task->fRebinning)
  {
    // no rebinning
    TH1D* hNoRebin_rebinned = TestRebin(hDesampledFlow,task);
    hNoRebin_rebinned->Write(Form("%s",hNoRebin_rebinned->GetName()));

    if(task->fSampleMerging)
    {
      TH1D* hMerged_rebinned = TestRebin(hMerged,task);
      hMerged_rebinned->Write(Form("%s",hMerged_rebinned->GetName()));
    }
  }

  // TCanvas* canTest = new TCanvas("canTest","canTest");
  // canTest->Divide(2,1);
  // canTest->cd(1);
  // hDesampled->Draw();
  // canTest->cd(2);
  // hDesampledFlow->Draw();

  if(hDesampled) delete hDesampled;
  if(hDesampledFlow) delete hDesampledFlow;
  if(profRebin) delete profRebin;
  if(histProj) delete histProj;
  if(list) delete list;
  if(listMerge) delete listMerge;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessDirect(FlowTask* task, Short_t iMultBin)
{
  Info("Processing direct task","ProcesDirect");
  if(!task) { Error("Task not valid!","ProcessDirect"); return kFALSE; }

  TList* listInput = 0x0;
  switch (task->fSpecies)
  {
    case FlowTask::kCharged:
      listInput = flFlowCharged;
      break;

    case FlowTask::kPion:
    case FlowTask::kKaon:
    case FlowTask::kProton:
      listInput = flFlowPID;
      break;

    default:
      Error("Task species not direct!","ProcessDirect");
      return kFALSE;
  }

  // preparing vn' samples
  TList* listFlow = new TList();
  TList* listCum = new TList();
  TList* listMerge = new TList();
  TProfile2D* prof2 = 0x0;
  TProfile2D* prof2pos = 0x0;
  TProfile2D* prof2neg = 0x0;
  TProfile* prof = 0x0;
  TProfile* profRebin = 0x0;
  TProfile* profRef = 0x0;
  TProfile* profRefRebin = 0x0;
  TH1D* hFlow = 0x0;

  Short_t binMultLow = 0;
  Short_t binMultHigh = 0;
  Double_t dRef = 0;

  for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
  {
    if(task->fMergePosNeg)
    {
      // loading Pos & Neg if fMergePosNeg is ON
      prof2pos = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_Pos_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample));
      prof2neg = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_Neg_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample));

      if(!prof2pos || !prof2neg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessDirect"); return kFALSE; }
      listMerge->Clear();
      listMerge->Add(prof2pos);
      listMerge->Add(prof2neg);

      prof2 = (TProfile2D*) listMerge->At(0)->Clone();
      prof2->Reset();
      Double_t mergeStatus = prof2->Merge(listMerge);
      if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessDirect"); return kFALSE; }
    }
    else
    {
      // loading single profile
      if(task->fInputTag.EqualTo(""))
      { prof2 = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_Pos_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample)); }
      else
      { prof2 = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_%s_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data(),iSample)); }
    }
    if(!prof2) { Error(Form("Profile sample %d does not exists.",iSample),"ProcessDirect"); return kFALSE; }

    // preparing Refs
    profRef = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_%s<2>_harm%d_gap%02.2g_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample));
    if(!profRef) { Error(Form("Profile sample %d does not exists",iSample),"ProcesPID"); return kFALSE; }

    profRefRebin = (TProfile*) profRef->Rebin(fiNumMultBins,Form("%s_rebin",profRef->GetName()),fdMultBins);
    dRef = profRefRebin->GetBinContent(iMultBin+1);
    profRefRebin->Draw();
    // printf("dRef %g\n",dRef);
    // NOTE: complains about Sumw2

    // rebinning according in mult bin
    binMultLow = prof2->GetXaxis()->FindFixBin(fdMultBins[iMultBin]);
    binMultHigh = prof2->GetXaxis()->FindFixBin(fdMultBins[iMultBin+1]) - 1;
    prof = prof2->ProfileY(Form("%s_cent%d",prof2->GetName(),iMultBin),binMultLow,binMultHigh);

    if(task->fNumPtBins > 0)
    {
      // rebinning according to pt bins
      profRebin = (TProfile*) prof->Rebin(task->fNumPtBins,Form("%s_rebin",prof->GetName()),task->fPtBinsEdges);
    }
    else
    {
      profRebin = (TProfile*) prof->Clone(Form("%s_rebin",prof->GetName()));
    }



    // making TH1D projection (to avoid handling of bin entries)
    TH1D* hCum = (TH1D*) profRebin->ProjectionX();
    hCum->SetName(Form("%s_Cum",prof->GetName()));
    listCum->Add(hCum);

    hFlow = (TH1D*) profRebin->ProjectionX();
    // printf("Pre %g\n",hFlow->GetBinContent(20));
    // printf("Scale %g\n",1/TMath::Sqrt(dRef));
    hFlow->SetName(Form("%s_Flow",prof->GetName()));
    hFlow->Scale(1/TMath::Sqrt(dRef));
    // printf("Post %g\n",hFlow->GetBinContent(20));
    listFlow->Add(hFlow);
  }
  delete listMerge;

  Debug(Form("Number of samples in list pre merging %d",listFlow->GetEntries()),"ProcessDirect");

  TH1D* hDesampled = DesampleList(listFlow,task,iMultBin);
  if(!hDesampled) { Error("Desampling unsuccesfull","ProcessDirect"); return kFALSE; }

  hDesampled->SetName(Form("hFlow2_%s_harm%d_gap%s_cent%d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
  hDesampled->SetTitle(Form("%s v_{%d}{2} | Gap %s | Cent %d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));

  // desampling cumulants
  TH1D* hDesampled_Cum = DesampleList(listCum,task,iMultBin);
  if(!hDesampled_Cum) { Error("Desampling unsuccesfull","ProcessDirect"); return kFALSE; }

  hDesampled_Cum->SetName(Form("hCum2_%s_harm%d_gap%s_cent%d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
  hDesampled_Cum->SetTitle(Form("%s d_{%d}{2} | Gap %s | Cent %d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));


  // saving to output file
  ffOutputFile->cd();
  hDesampled->Write();
  hDesampled_Cum->Write();

  delete listFlow;
  delete prof;
  delete profRebin;
  delete profRefRebin;
  delete hFlow;
  delete hDesampled;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessReconstructed(FlowTask* task,Short_t iMultBin)
{
  Info("Processing task","ProcessReconstructed");
  if(!task) { Error("Task not valid!","ProcessReconstructed"); return kFALSE; }
  // if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessReconstructed"); return kFALSE; }


  TList* listMerge = 0x0;
  // preparing particle dependent variables for switch
  //  -- input histos / profiles with entries and correlations
  TH3D* histEntries = 0x0;
  TH3D* histEntriesPos = 0x0;
  TH3D* histEntriesNeg = 0x0;
  TH3D* histBG = 0x0; // entries for BG (phi)
  TH3D* histBGPos = 0x0; // entries for BG (phi)
  TH3D* histBGNeg = 0x0; // entries for BG (phi)
  TProfile3D* profFlow = 0x0;
  TProfile3D* profFlowPos = 0x0;
  TProfile3D* profFlowNeg = 0x0;
  //  -- naming variables
  TString sSpeciesName; // in objects name
  TString sSpeciesLabel; // LaTeX for titles

  // checking particles species and assigning particle dependent variables
  switch (task->fSpecies)
  {
    case FlowTask::kPhi :
      sSpeciesName = TString("Phi");
      sSpeciesLabel = TString("#phi");

      if(task->fMergePosNeg)
      {
        // loading Pos & Neg if fMergePosNeg is ON
        // merging profiles
        profFlowPos = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        profFlowNeg = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_%s<2>_harm%d_gap%02.2g_Neg",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));

        if(!profFlowPos || !profFlowNeg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessDirect"); return kFALSE; }
        listMerge = new TList();
        listMerge->Add(profFlowPos);
        listMerge->Add(profFlowNeg);

        profFlow = (TProfile3D*) listMerge->At(0)->Clone();
        profFlow->Reset();
        Double_t mergeStatus = profFlow->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        // merging histos
        histEntriesPos = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        histEntriesNeg = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_Neg",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        if(!histEntriesPos || !histEntriesNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histEntriesPos);
        listMerge->Add(histEntriesNeg);

        histEntries = (TH3D*) listMerge->At(0)->Clone();
        histEntries->Reset();
        mergeStatus = histEntries->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        histBGPos = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_Pos",10*task->fEtaGap));
        histBGNeg = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_Neg",10*task->fEtaGap));
        if(!histBGPos || !histBGNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histBGPos);
        listMerge->Add(histBGNeg);

        histBG = (TH3D*) listMerge->At(0)->Clone();
        histBG->Reset();
        mergeStatus = histBG->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;
      }
      else
      {
        // loading single profile
        if(task->fInputTag.EqualTo(""))
        {
          histEntries = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
          histBG = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_Pos",10*task->fEtaGap));
          profFlow = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_<2>_%sharm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        }
        else
        {
          histEntries = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_%s",fsGlobalProfNameLabel.Data(),10*task->fEtaGap,task->fInputTag.Data()));
          histBG = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_%s",10*task->fEtaGap,task->fInputTag.Data()));
          profFlow = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_%s<2>_harm%d_gap%02.2g_%s",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data()));
        }
      }
    break;

    case FlowTask::kK0s :
      sSpeciesName = TString("K0s");
      sSpeciesLabel = TString("K^{0}_{S}");

      if(task->fMergePosNeg)
      {
        // loading Pos & Neg if fMergePosNeg is ON
        // merging TProfiles
        profFlowPos = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        profFlowNeg = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_Neg",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        if(!profFlowPos || !profFlowNeg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(profFlowPos);
        listMerge->Add(profFlowNeg);

        profFlow = (TProfile3D*) listMerge->At(0)->Clone();
        profFlow->Reset();
        Double_t mergeStatus = profFlow->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging profiles unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        // merging histos
        histEntriesPos = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        histEntriesNeg = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_Neg",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        if(!histEntriesPos || !histEntriesNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histEntriesPos);
        listMerge->Add(histEntriesNeg);

        histEntries = (TH3D*) listMerge->At(0)->Clone();
        histEntries->Reset();
        mergeStatus = histEntries->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;
      }
      else
      {
        // loading single profile
        if(task->fInputTag.EqualTo(""))
        {
          histEntries = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
          profFlow = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        }
        else
        {
          histEntries = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_%s",fsGlobalProfNameLabel.Data(),10*task->fEtaGap,task->fInputTag.Data()));
          profFlow = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_%s",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data()));
        }
      }
    break;

    case FlowTask::kLambda :
      sSpeciesName = TString("Lambda");
      sSpeciesLabel = TString("#Lambda/#bar{#Lambda}");

      if(task->fMergePosNeg)
      {
        // loading Pos & Neg if fMergePosNeg is ON
        // merging profiles
        profFlowPos = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        profFlowNeg = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_Neg",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));

        if(!profFlowPos || !profFlowNeg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }
        listMerge = new TList();
        listMerge->Add(profFlowPos);
        listMerge->Add(profFlowNeg);

        profFlow = (TProfile3D*) listMerge->At(0)->Clone();
        profFlow->Reset();
        Double_t mergeStatus = profFlow->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        // merging histos
        histEntriesPos = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        histEntriesNeg = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_Neg",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        if(!histEntriesPos || !histEntriesNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histEntriesPos);
        listMerge->Add(histEntriesNeg);

        histEntries = (TH3D*) listMerge->At(0)->Clone();
        histEntries->Reset();
        mergeStatus = histEntries->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;
      }
      else
      {
        // loading single profile
        if(task->fInputTag.EqualTo(""))
        {
          histEntries = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
          profFlow = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        }
        else
        {
          histEntries = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_%s",fsGlobalProfNameLabel.Data(),10*task->fEtaGap,task->fInputTag.Data()));
          profFlow = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_%s",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data()));
        }
      }
    break;

    default:
      Error("Task species not V0s nor Phi!","ProcessReconstructed");
      return kFALSE;
  }

  if(!histEntries) { Error("Entries histos not found!","ProcessReconstructed"); return kFALSE; }
  if(!profFlow) { Error("Cumulant histos not found!","ProcessReconstructed"); return kFALSE; }

  // check if suggest pt binning flag is on if of Pt binning is not specified
  if(task->fSuggestPtBins || task->fNumPtBins < 1)
  {
    SuggestPtBinning(histEntries,profFlow,task,iMultBin);
  }

  if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessReconstructed"); return kFALSE; }

  // task->PrintTask();

  if(!PrepareSlices(iMultBin,task,profFlow,histEntries,histBG)) return kFALSE;

  TH1D* hFlow = 0x0;
  if(!fFlowFitCumulants)
  {
    hFlow = new TH1D(Form("hFlow2_%s_harm%d_gap%02.2g_cent%d",sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin),Form("%s: v_{%d}{2 | Gap %g} (%g - %g); #it{p}_{T} (GeV/#it{c})",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1]), task->fNumPtBins,task->fPtBinsEdges);
  }
  else
  {
    hFlow = new TH1D(Form("hCum2_%s_harm%d_gap%02.2g_cent%d",sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin),Form("%s: d_{%d}{2 | Gap %g} (%g - %g); #it{p}_{T} (GeV/#it{c})",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1]), task->fNumPtBins,task->fPtBinsEdges);
  }

  TH1D* hInvMass = 0x0;
  TH1D* hInvMassBG = 0x0;
  TH1D* hFlowMass = 0x0;
  Double_t dFlow = 0, dFlowError = 0; // containers for flow extraction results
  TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1600,1200); // canvas for fitting results

  TCanvas* canFlowAll = new TCanvas("canFlowAll","canFlowAll",1600,1200);
  TCanvas* canInvMassAll = new TCanvas("canInvMassAll","canInvMassAll",1600,1200);
  canFlowAll->Divide(3,ceil(task->fNumPtBins/3.));
  canInvMassAll->Divide(3,ceil(task->fNumPtBins/3.));

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.1);
  if(task->fSpecies == FlowTask::kPhi) {  latex->SetTextSize(0.08); }
  latex->SetNDC();

  TLatex* latex2 = new TLatex();
  // latex2->SetTextFont(43);
  // latex2->SetTextSize(40);
  latex2->SetNDC();

  for(Short_t binPt(0); binPt < task->fNumPtBins; binPt++)
  {
    hInvMass = task->fVecHistInvMass->at(binPt);
    hFlowMass = task->fVecHistFlowMass->at(binPt);

    hInvMass->SetTitle(Form("%s: InvMass dist (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
    hFlowMass->SetTitle(Form("%s: FlowMass (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));

    hInvMass->SetMarkerStyle(kFullCircle);
    hFlowMass->SetMarkerStyle(kFullCircle);

    TList* listFits = new TList();

    // extracting flow
    switch (task->fSpecies)
    {
      case FlowTask::kPhi :
        hInvMassBG = task->fVecHistInvMassBG->at(binPt);
        if( !ExtractFlowOneGo(task,hInvMass,hInvMassBG,hFlowMass,dFlow,dFlowError,canFitInvMass,listFits) ) { Warning("Flow extraction unsuccesfull","ProcessReconstructed"); return kFALSE; }
        // if( !ExtractFlowPhiOneGo(task,hInvMass,hInvMassBG,hFlowMass,dFlow,dFlowError,canFitInvMass,listFits) ) { Warning("Flow extraction unsuccesfull","ProcessReconstructed"); return kFALSE; }
      break;

      case FlowTask::kK0s :
        if( !ExtractFlowOneGo(task,hInvMass,0x0,hFlowMass,dFlow,dFlowError,canFitInvMass,listFits) ) { Warning("Flow extraction unsuccesfull (one go)","ProcessReconstructed"); return kFALSE; }
      break;

      case FlowTask::kLambda :
        if( !ExtractFlowOneGo(task,hInvMass,0x0,hFlowMass,dFlow,dFlowError,canFitInvMass,listFits) ) { Warning("Flow extraction unsuccesfull (one go)","ProcessReconstructed"); return kFALSE; }
      break;

      default :
        Error("Uknown species","ProcessReconstructed");
        return kFALSE;
    }

    ffFitsFile->cd();
    listFits->Write(Form("fits_%s_cent%d_pt%d",sSpeciesName.Data(),iMultBin,binPt),TObject::kSingleKey);


    gSystem->mkdir(Form("%s/fits/",fsOutputFilePath.Data()));


    TF1* fitInvMass2 = (TF1*) listFits->At(1);
    // printf("InvMass chi2 %g\n",fitInvMass2->GetChisquare());
    TF1* fitFlowMass2 = (TF1*) listFits->At(5);
    // printf("FlowMass chi2 %g\n",fitFlowMass2->GetChisquare());

    canFitInvMass->cd(1);
    // if(task->fSpecies == FlowTask::kPhi) canFitInvMass->cd(2);
    latex2->DrawLatex(0.17,0.85,Form("#color[9]{pt %g-%g GeV/c (%g-%g%%)}",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
    canFitInvMass->cd(2);
    latex2->DrawLatex(0.17,0.85,Form("#color[9]{pt %g-%g GeV/c (%g-%g%%)}",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
    canFitInvMass->SaveAs(Form("%s/fits/Fit_%s_n%d2_gap%02.2g_cent%d_pt%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,binPt,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

    canFlowAll->cd(binPt+1);
    hFlowMass->SetLabelFont(43,"XY");
    hFlowMass->SetLabelSize(18,"XY");
    hFlowMass->DrawCopy();
    TF1* fitVn = (TF1*) listFits->FindObject("fitVn");
    fitVn->DrawCopy("same");
    latex->DrawLatex(0.13,0.8,Form("#color[9]{%1.1f-%1.1f GeV/c (%g-%g%%)}",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
    latex->DrawLatex(0.13,0.2,Form("#color[9]{#chi2/ndf = %.1f/%d = %.1f (p=%.3f)}",fitFlowMass2->GetChisquare(), fitFlowMass2->GetNDF(),fitFlowMass2->GetChisquare()/fitFlowMass2->GetNDF(),fitFlowMass2->GetProb()));
    latex->DrawLatex(0.13,0.33,Form("#color[9]{d_{2} = %.2g +- %.2g }",dFlow,dFlowError));

    canInvMassAll->cd(binPt+1);
    // gPad->SetLogy();
    hInvMass->SetLabelFont(43,"XY");
    hInvMass->SetLabelSize(18,"XY");
    // hInvMass->SetMinimum(1);
    hInvMass->DrawCopy();
    TF1* fitInvMass = (TF1*) listFits->FindObject("fitMass");
    fitInvMass->DrawCopy("same");
    latex->DrawLatex(0.13,0.2,Form("#color[9]{#chi2/ndf = %.1f/%d = %.1f (p=%.3f)}",fitInvMass2->GetChisquare(), fitInvMass2->GetNDF(),fitInvMass2->GetChisquare()/fitInvMass2->GetNDF(),fitInvMass2->GetProb()));
    latex->DrawLatex(0.13,0.8,Form("#color[9]{%1.1f-%1.1f GeV/c (%g-%g%%)}",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));


    if(TMath::Abs(dFlow) > 1 )
    {
      hFlow->SetBinContent(binPt+1,0);
      hFlow->SetBinError(binPt+1,0);
    }
    else
    {
      hFlow->SetBinContent(binPt+1,dFlow);
      hFlow->SetBinError(binPt+1,dFlowError);
    }

  } // endfor {binPt}

  // task->fCanvas->Draw();
  canFlowAll->SaveAs(Form("%s/FlowMassFits_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
  canInvMassAll->SaveAs(Form("%s/InvMassFits_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());


  TCanvas* cFlow = new TCanvas("cFlow","cFlow");
  cFlow->cd();
  hFlow->SetStats(0);
  hFlow->Draw();
  cFlow->SaveAs(Form("%s/Flow_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

  ffOutputFile->cd();
  hFlow->Write();

  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::DesampleList(TList* list, FlowTask* task, Short_t iMultBin)
{
  if(!task) { Error("FlowTask does not exists","DesampleList"); return 0x0; }
  if(!list) { Error("List does not exists","DesampleList"); return 0x0; }
  if(list->GetEntries() < 1) { Error("List is empty","DesampleList"); return 0x0; }
  if(list->GetEntries() != task->fNumSamples) { Warning("Number of list entries is different from task number of samples","DesampleList"); }

  Debug(Form("Number of samples in list pre-desampling: %d",list->GetEntries()),"DesampleList");

  TH1D* hTempSample = (TH1D*) list->At(0);
  TH1D* hDesampled = (TH1D*) hTempSample->Clone(Form("%s_Desampled",hTempSample->GetName()));
  hDesampled->Reset();

  Double_t content = 0;
  Double_t error = 0;

  std::vector<Double_t> vecContents;
  std::vector<Double_t> vecErrors;
  Double_t dMean = 0;
  Double_t dRMS = 0;

  Double_t dSum = 0;
  Double_t dW = 0;
  Double_t dAverage = 0;
  Double_t dAve_err = 0;

  for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++)
  {
    dSum = 0;
    dW = 0;
    dAverage = 9.;
    dAve_err = 9.;

    vecContents.clear();
    vecErrors.clear();

    for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
    {
      hTempSample = (TH1D*) list->At(iSample);
      if(!hTempSample) { Warning(Form("Sample %d not found! Skipping!",iSample),"DesampleList"); continue; }

      content = hTempSample->GetBinContent(bin);
      error = hTempSample->GetBinError(bin);

      vecContents.push_back(content);
      vecErrors.push_back(error);

      if(error <= 0.) continue;

      dSum += content / TMath::Power(error,2);
      dW += 1 / TMath::Power(error,2);
      Debug(Form("Sample: %d | bin %d | %g +- %g",iSample,bin,content,error),"DesampleList");
    }

    Debug(Form(" --- bin %d | Sum %g +- %g",bin,dSum,dW),"DesampleList");

    if(dSum == 0. && dW == 0.) continue; // skipping empty bins

    if(dW > 0.)
    {
      dAverage = dSum / dW;
      dAve_err = TMath::Sqrt(1/dW);
    }

    Debug(Form("W average | bin %d | %g +- %g",bin,dAverage,dAve_err),"DesampleList");

    // working on weighted mean and RMS
    if(vecContents.size() != vecErrors.size()) { Warning("Something went wrong: Vector with contents and error have not the same size","DesampleList"); }

    dSum = 0;
    dW = 0;

    for(ULong_t iSample(0); iSample < vecContents.size(); iSample++)
    {
      content = vecContents.at(iSample);
      error = vecErrors.at(iSample);

      // weighted mean
      // dW += TMath::Power(error, -2);
      // dSum += content * TMath::Power(error, -2);

      // non-weighted mean
      dW++;
      dSum += content;
    }

    dMean = dSum/dW;

    dSum = 0;
    for(ULong_t iSample(0); iSample < vecContents.size(); iSample++)
    {
      content = vecContents.at(iSample);
      // error = vecErrors.at(iSample);
      error = 1.;

      dSum += TMath::Power( (dMean - content) / error ,2);
    }
    dRMS = TMath::Sqrt(dSum / (dW-1));
    Debug(Form("Mean | bin %d | %g +- %g",bin,dMean,dRMS),"DesampleList");

    //ratio->SetBinContent(bin, dAverage / merged->GetBinContent(bin));
    //ratioErr->SetBinContent(bin, dAve_err / merged->GetBinError(bin));

    if(task->fDesampleUseRMS)
    {
      hDesampled->SetBinContent(bin,dMean);
      // hDesampled->SetBinError(bin,dRMS);
      hDesampled->SetBinError(bin,dRMS/TMath::Sqrt(task->fNumSamples));
    }
    else
    {
      hDesampled->SetBinContent(bin,dAverage);
      hDesampled->SetBinError(bin,dAve_err);
    }
    Debug(Form("Desampled: bin %d | %g +- %g\n",bin,hDesampled->GetBinContent(bin),hDesampled->GetBinError(bin)),"DesampleList");
  }

  // at this point, hDesampled is ready

  // getting copy which does not affect histo which is returned
  TH1D* hDesampledClone = (TH1D*) hDesampled->Clone(Form("%sClone",hDesampled->GetName()));

  TList* listOutput = new TList(); // list for collecting all QA histos

  // doing QA plots with spread, etc.
  TCanvas* canDesample = new TCanvas(Form("canDesample_%s",task->fName.Data()),Form("canDesample_%s",task->fName.Data()),1200,400);
  canDesample->Divide(3,1);

  TH1D* hTempRatio = 0x0;
  TH1D* hTempError = 0x0;

  TLine* lineUnity = new TLine();
  lineUnity->SetLineColor(kRed);
  lineUnity->SetLineWidth(3);

  canDesample->cd(1);
  hDesampledClone->SetStats(kFALSE);
  hDesampledClone->SetFillColor(kBlue);
  hDesampledClone->SetStats(kFALSE);
  hDesampledClone->SetMarkerStyle(20);
  hDesampledClone->SetMarkerSize(0.5);
  hDesampledClone->SetMarkerColor(kRed);
  hDesampledClone->DrawCopy("E2");

  for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
  {
    hTempSample = (TH1D*) list->At(iSample);
    if(!hTempSample) { Warning(Form("Sample %d not found during plotting QA! Skipping!",iSample),"DesampleList"); continue; }

    canDesample->cd(1);
    hTempSample->SetStats(kFALSE);
    hTempSample->SetLineColor(30+2*iSample);
    hTempSample->SetMarkerColor(30+2*iSample);
    hTempSample->SetMarkerStyle(24);
    hTempSample->SetMarkerSize(0.5);
    hTempSample->DrawCopy("hist p same");

    hTempRatio = (TH1D*) hTempSample->Clone(Form("%s_ratio",hTempSample->GetName()));
    hTempRatio->Divide(hDesampled);
    hTempRatio->SetYTitle("Value: final / sample");
    hTempRatio->SetTitleOffset(1.2,"Y");

    canDesample->cd(2);
    hTempRatio->SetMinimum(0.6);
    hTempRatio->SetMaximum(1.4);
    hTempRatio->Draw("hist p same");

    hTempError = (TH1D*) hTempSample->Clone(Form("%s_error",hTempSample->GetName()));
    for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hTempSample->GetBinError(bin)); }

    canDesample->cd(3);
    hTempError->SetMinimum(0.);
    hTempError->SetMaximum(1.5*hTempError->GetMaximum());
    hTempError->SetYTitle("Uncertainty");
    hTempError->SetTitleOffset(1.2,"Y");

    hTempError->Draw("hist p same");

    listOutput->Add(hTempSample);
    listOutput->Add(hTempRatio);
    listOutput->Add(hTempError);
  }

  canDesample->cd(1);
  hDesampledClone->DrawCopy("hist p same");

  canDesample->cd(2);
  lineUnity->DrawLine(hTempRatio->GetXaxis()->GetXmin(),1,hTempRatio->GetXaxis()->GetXmax(),1);

  hTempError = (TH1D*) hDesampledClone->Clone(Form("%s_error",hDesampled->GetName()));
  for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hDesampledClone->GetBinError(bin)); }
  listOutput->Add(hTempError);

  canDesample->cd(3);
  hTempError->Draw("hist p same");

  // saving QA plots
  canDesample->SaveAs(Form("%s/Desampling_%s_harm%d_gap%g_cent%d_%s.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,task->fName.Data(),fsOutputFileFormat.Data()));

  Info("Saving desampling QA into output file","DesampleList");
  ffDesampleFile->cd();
  listOutput->Add(canDesample);
  listOutput->Write(Form("Desampling_%s_cent%d_%s",task->GetSpeciesName().Data(),iMultBin,task->fName.Data()),TObject::kSingleKey);

  // deleting created stuff
  delete listOutput;
  // delete canDesample;
  delete lineUnity;
  // if(hTempSample) delete hTempSample;
  // delete hTempRatio;
  // delete hTempError;
  // delete hDesampledClone;

  return hDesampled;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor, TH3D* h3Entries, TH3D* h3EntriesBG)
{

  if(!task) { Error("Input task not found!","PrepareSlices"); return kFALSE; }
  if(!h3Entries) { Error("Input hist with entries not found!","PrepareSlices"); return kFALSE; }
  if(!p3Cor) { Error("Input profile with correlations not found!","PrepareSlices"); return kFALSE; }
  if(multBin < 0 || multBin > fiNumMultBins) { Error("Wrong multiplicity bin index (not in range)!","PrepareSlices"); return kFALSE; }

  // cleaning the vectros with flow-mass and inv. mass plots
  if(task->fVecHistInvMass->size() > 0) task->fVecHistInvMass->clear();
  if(task->fVecHistInvMassBG->size() > 0) task->fVecHistInvMassBG->clear();
  if(task->fVecHistFlowMass->size() > 0) task->fVecHistFlowMass->clear();

  const Short_t binMultLow = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin]);
  const Short_t binMultHigh = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin+1]) - 1;
  // printf("Mult: %g(%d) -  %g(%d)\n",fdMultBins[multBin],binMultLow,fdMultBins[multBin+1],binMultHigh);

  Double_t dRefFlow = 1.;
  Double_t dRefFlowErr = 1.;

  if(!fFlowFitCumulants)
  {
    // loading reference flow, if not found, it will be prepared
    TH1D* hRefFlow = 0x0;
    hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
    if(!hRefFlow)
    {
      Warning("Relevant Reference flow not found within output file.","PrepareSlices");
      ffOutputFile->ls();
      // return kFALSE;

      Info("Creating relevant reference flow task.","PrepareSlices");
      FlowTask* taskRef = new FlowTask(FlowTask::kRefs,"Ref");
      taskRef->SetHarmonics(task->fHarmonics);
      taskRef->SetEtaGap(task->fEtaGap);
      taskRef->SetNumSamples(task->fNumSamples);
      taskRef->SetInputTag(task->fInputTag);
      if(ProcessRefs(taskRef))
      {
        hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
        if(!hRefFlow) {  Error("Automated Refs task completed, but RefFlow not found!","PrepareSlices"); return kFALSE; }
      }
      else { Error("Something went wrong when running automatic refs flow task:","PrepareSlices"); taskRef->PrintTask(); return kFALSE; }
    }
    dRefFlow = hRefFlow->GetBinContent(multBin+1);
    dRefFlowErr = hRefFlow->GetBinError(multBin+1);
    Debug(Form("Ref (bin %d): %g +- %g\n",multBin,dRefFlow,dRefFlowErr),"PrepareSlices");
  }

  // loop over pt
  Short_t binPtLow = 0;
  Short_t binPtHigh = 0;
  TH1D* hInvMass_temp = 0x0;
  TH1D* hInvMassBG_temp = 0x0;
  TProfile3D* prof3Flow_temp = 0x0;
  TProfile2D* prof2FlowMass_temp = 0x0;
  TProfile* profFlowMass_temp = 0x0;
  TH1D* hFlowMass_temp = 0x0;

  prof3Flow_temp = (TProfile3D*) p3Cor->Clone(Form("prof3Flow_temp_cent%d",multBin));
  prof3Flow_temp->GetXaxis()->SetRange(binMultLow,binMultHigh);
  prof2FlowMass_temp = Project3DProfile(prof3Flow_temp);

  Short_t iNumPtBins = task->fNumPtBins;

  TCanvas* canFlowMass = new TCanvas("canFlowMass","FlowMass",1400,600);
  TCanvas* canInvMass = new TCanvas("canInvMass","InvMass",1400,600);
  TCanvas* canInvMassBG = new TCanvas("canInvMassBG","InvMassBG",1400,600);
  canFlowMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canInvMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canInvMassBG->Divide(5,std::ceil(iNumPtBins/5)+1);

  Double_t dContent = 0, dError = 0;

  for(Short_t binPt(0); binPt < iNumPtBins; binPt++)
  {
    // estimating pt edges
    binPtLow = h3Entries->GetYaxis()->FindFixBin(task->fPtBinsEdges[binPt]);
    binPtHigh = h3Entries->GetYaxis()->FindFixBin(task->fPtBinsEdges[binPt+1]) - 1; // for rebin both bins are included (so that one needs to lower)
    printf("   Pt: %g(%d) -  %g(%d)\n",task->fPtBinsEdges[binPt],binPtLow,task->fPtBinsEdges[binPt+1],binPtHigh);

    // rebinning entries based on mult & pt binning
    hInvMass_temp = (TH1D*) h3Entries->ProjectionZ(Form("hInvMass_cent%d_pt%d",multBin,binPt),binMultLow,binMultHigh,binPtLow,binPtHigh,"e");
    if(h3EntriesBG) hInvMassBG_temp = (TH1D*) h3EntriesBG->ProjectionZ(Form("hInvMassBG_cent%d_pt%d",multBin,binPt),binMultLow,binMultHigh,binPtLow,binPtHigh,"e");

    // checking if rebinning inv mass hist
    if(task->fRebinInvMass > 1) { hInvMass_temp->Rebin(task->fRebinInvMass); if(h3EntriesBG){ hInvMassBG_temp->Rebin(task->fRebinInvMass); } }


    // projection of flow-mass profile
    profFlowMass_temp = (TProfile*) prof2FlowMass_temp->ProfileX(Form("profFlowMass_cent%d_pt%d",multBin,binPt),binPtLow,binPtHigh);

    // checking for rebinning the flow-mass profile
    if(task->fRebinFlowMass > 1) { profFlowMass_temp->Rebin(task->fRebinFlowMass); }

    hFlowMass_temp = (TH1D*) profFlowMass_temp->ProjectionX(Form("hFlowMass_cent%d_pt%d",multBin,binPt));
    // NOTE: this is the ONLY (for some freaking reason) way how to get proper TH1 wth <<2>> out of TProfile3D

    // scaling flow-mass with reference flow
    if(!fFlowFitCumulants)
    {
      for(Short_t bin(1); bin < hFlowMass_temp->GetNbinsX()+1; bin++)
      {
        dContent = hFlowMass_temp->GetBinContent(bin);
        dError = hFlowMass_temp->GetBinError(bin);

        if(dContent == 0. || dError == 0.) continue;

        hFlowMass_temp->SetBinContent(bin,dContent/dRefFlow);
        hFlowMass_temp->SetBinError(bin, TMath::Sqrt( TMath::Power(dError/dRefFlow,2) + TMath::Power(dRefFlowErr*dContent/(dRefFlow*dRefFlow),2) - 2*(dError*dRefFlowErr*dContent*TMath::Power(dRefFlow,-3))) );
        // printf("%g | %g \n", TMath::Sqrt( TMath::Power(dError/dRefFlow,2) + TMath::Power(dRefFlowErr*dContent/(dRefFlow*dRefFlow),2) - 2*(dError*dRefFlowErr*dContent*TMath::Power(dRefFlow,-3))), TMath::Sqrt( TMath::Power(dError/dRefFlow,2) + TMath::Power(dRefFlowErr*dContent/(dRefFlow*dRefFlow),2)) );
      }
    }

    // hInvMass_temp->SetTitle(Form("%s: Inv. Mass (|#Delta#eta| > %g)",sSpeciesLabel.Data(),0.05));
    // hFlowMass_temp->SetTitle("TEST");

    // ready to fitting
    task->fVecHistFlowMass->push_back(hFlowMass_temp);
    task->fVecHistInvMass->push_back(hInvMass_temp);
    if(h3EntriesBG) task->fVecHistInvMassBG->push_back(hInvMassBG_temp);

    canFlowMass->cd(binPt+1);
    hFlowMass_temp->Draw();

    canInvMass->cd(binPt+1);
    hInvMass_temp->Draw();

    if(h3EntriesBG)
    {
      canInvMassBG->cd(binPt+1);
      hInvMassBG_temp->Draw();
    }

  } // endfor {binPt}: over Pt bins

  printf(" # of slices: InvMass: %lu | InvMassBG %lu | FlowMass %lu\n",task->fVecHistInvMass->size(),task->fVecHistInvMassBG->size(),task->fVecHistFlowMass->size());

  gSystem->mkdir(Form("%s/slices/",fsOutputFilePath.Data()));
  canFlowMass->SaveAs(Form("%s/slices/Slices_FlowMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  canInvMass->SaveAs(Form("%s/slices/Slices_InvMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  if(h3EntriesBG) canInvMassBG->SaveAs(Form("%s/slices/Slices_InvMassBG_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));

  if(task->fVecHistInvMass->size() < 1 || task->fVecHistFlowMass->size() < 1 || task->fVecHistFlowMass->size() != task->fVecHistInvMass->size()) { Error("Output vector empty. Something went wrong","PrepareSlices"); return kFALSE; }
  if(h3EntriesBG && (task->fVecHistInvMassBG->size() < 1 || task->fVecHistInvMassBG->size() != task->fVecHistInvMass->size()) ) { Error("Output vector empty. Something went wrong with BG histograms","PrepareSlices"); return kFALSE; }
  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::TestRebin(TH1D* hOrig, FlowTask* task)
{
  if(!hOrig) { Error("Original histogram not found!","TestRebin"); return 0x0; }
  if(!task) { Error("Task not found!","TestRebin"); return 0x0; }
  TH1D* hRebin = (TH1D*) hOrig->Rebin(fiNumMultBins,Form("%s_testRebin",hOrig->GetName()),fdMultBins);

  // rebinning
  Short_t numBins = fiNumMultBins;
  Double_t* multBins = fdMultBins;
  Short_t binIndex = 0;


  const Short_t iNumBins = hOrig->GetNbinsX();

  Double_t dSum = 0;
  Double_t dSumWeights = 0;
  Double_t dContent = 0;
  Double_t dWeight = 0;
  for(Short_t mult = 1; mult < iNumBins+1; mult++)
  {
    if(hOrig->GetBinLowEdge(mult) < multBins[binIndex]) continue;

    dContent = hOrig->GetBinContent(mult);
    dWeight =  TMath::Power(hOrig->GetBinError(mult),-2);

    dSumWeights += dWeight;
    dSum += dContent*dWeight;

    if( hOrig->GetBinLowEdge(mult+1) == multBins[binIndex+1] )
    {
      hRebin->SetBinContent(binIndex+1, dSum / dSumWeights);
      hRebin->SetBinError(binIndex+1, TMath::Sqrt(1/dSumWeights));

      dSumWeights = 0;
      dSum = 0;
      binIndex++;
    }
  }

  return hRebin;
}
//_____________________________________________________________________________
void ProcessUniFlow::AddTask(FlowTask* task)
{
  if(!task) return;

  if(task->fSpecies == FlowTask::kRefs) fvTasks.insert(fvTasks.begin(), task);
  else fvTasks.push_back(task);

  return;
}
//_____________________________________________________________________________
void ProcessUniFlow::SuggestPtBinning(TH3D* histEntries, TProfile3D* profFlowOrig, FlowTask* task, Short_t binMult)
{
  if(!histEntries) return;
  if(!profFlowOrig) return;
  if(!task) return;

  Short_t binMultLow = histEntries->GetXaxis()->FindFixBin(fdMultBins[binMult]);
  Short_t binMultHigh = histEntries->GetXaxis()->FindFixBin(fdMultBins[binMult+1]) - 1; // for rebin both bins are included (so that one needs to lower)

  TH1D* hPtProj = (TH1D*) histEntries->ProjectionY("hPtProj",binMultLow,binMultHigh);

  const Double_t dMinEntries = task->fSuggestPtBinEntries;
  std::vector<Double_t> vecBins;
  std::vector<Double_t> vecContents;

  printf("vector pre: %lu\n", vecBins.size());
  Double_t dCount = 0;
  const Short_t iNBins = hPtProj->GetNbinsX();
  Short_t iBin = 1;

  // if some of the pt bins are set, the suggestion started with first bin after the last edge
  if(task->fNumPtBins > 0 && task->fNumPtBins < iNBins)
  {
    iBin = hPtProj->FindFixBin(task->fPtBinsEdges[task->fNumPtBins]) + 1;
    printf("LastBin %d (%g) | iBin %d (%g)\n",task->fNumPtBins,task->fPtBinsEdges[task->fNumPtBins], iBin, hPtProj->GetBinLowEdge(iBin));
  }
  else { vecBins.push_back(hPtProj->GetXaxis()->GetXmin()); } // pushing first edge

  for(Short_t bin(iBin); bin < iNBins+1; bin++)
  {
    dCount += hPtProj->GetBinContent(bin);
    if(dCount > dMinEntries)
    {
      vecBins.push_back(hPtProj->GetXaxis()->GetBinUpEdge(bin));
      vecContents.push_back(dCount);
      dCount = 0;
    }
  }
  vecContents.push_back(dCount); // pushing last (remaining) bin content
  vecBins.push_back(hPtProj->GetXaxis()->GetXmax()); // pushing last edge (low edge of underflowbin)
  const Short_t iNumBinEdges = vecBins.size();
  printf("vector post: %hd\n", iNumBinEdges);

  // filling internal pt binning with suggestion result
  for(Short_t index(0); index < iNumBinEdges; index++)
  {
    task->fPtBinsEdges[index + task->fNumPtBins] = vecBins.at(index); // NOTE +1 ???
  }
  task->fNumPtBins += iNumBinEdges-1;


  // TODO hBinsContent reimplmente
  // TH1D* hBinsContent = new TH1D("hBinsContent","BinContent", task->fNumPtBins, task->fPtBinsEdges);
  // for(Short_t iBin(0); iBin < task->-1; iBin++)
  // {
  //   hBinsContent->SetBinContent(iBin+1,vecContents.at(iBin));
  // }

  TLine* line = new TLine();
  line->SetLineColor(kRed);

  TCanvas* cPtBins = new TCanvas("cPtBins","cPtBins",600,600);
  cPtBins->Divide(1,2);
  cPtBins->cd(1);
  hPtProj->Draw();
  cPtBins->cd(2);
  // hBinsContent->SetMinimum(0);
  // hBinsContent->Draw();
  // line->DrawLine(task->fPtBinsEdges[0],dMinEntries, task->fPtBinsEdges[task->fNumPtBins], dMinEntries);
  cPtBins->cd(1);

  Double_t dPt = 0;
  printf("Suggested binning: ");
  for(Short_t index(0); index < task->fNumPtBins+1; index++)
  {
    dPt = task->fPtBinsEdges[index];
    printf("%hd (%g) | ",index,dPt);
    line->DrawLine(dPt,0,dPt,1.05*hPtProj->GetMaximum());
  }
  printf("\n");

  cPtBins->SaveAs(Form("%s/suggestBins/SuggestPtBins_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,binMult,fsOutputFileFormat.Data()));
}
//_____________________________________________________________________________
void ProcessUniFlow::TestProjections()
{
  Info("Testing profile projections");

  ffInputFile->cd("UniFlow");
  // ffInputFile->ls();

  // TList* lFlow = (TList*) gDirectory->Get("UniFlow");
  // lFlow->ls();

  //
  // TList* lRef = (TList*) lFlow->FindObject("fFlowRefs");
  // TProfile* Ref = (TProfile*) lRef->FindObject("fpRefs_<2>_harm2_gap-10_sample0");
  // if(!Ref) { Error("NotFound"); return; }
  // // Ref->Draw();
  // // NOTE: reference flow works
  //
  // // projections charged flow
  // TList* lCharged = (TList*) lFlow->FindObject("fFlowCharged");
  // // lCharged->ls();
  // TProfile2D* p2Charged = (TProfile2D*) lCharged->FindObject("fp2Charged_<2>_harm2_gap-10_sample0");
  // if(!p2Charged) { Error("NotFound"); return; }
  // TProfile* p2ChargedProjY = p2Charged->ProfileY("p2ChargedProjY",10,10);

  // TCanvas* cCharged = new TCanvas("cCharged","cCharged",1000,1000);
  // cCharged->Divide(1,2);
  // cCharged->cd(1);
  // p2Charged->Draw("colz");
  // cCharged->cd(2);
  // p2ChargedProjY->Draw();

  // projections V0s
  // TList* lV0s = (TList*) gDirectory->Get("Flow_V0s_UniFlow");
  TList* lV0s = (TList*) gDirectory->Get("Flow_Phi_UniFlow");
  if(!lV0s) { Error("NotFound"); return; }
  lV0s->ls();

  // entries
  TH3D* h3K0sEntries = (TH3D*) lV0s->FindObject("fh3PhiEntriesSignal_gap-10")->Clone("h3K0sEntries");
  // TH3D* h3K0sEntries = (TH3D*) lV0s->FindObject("fh3V0sEntriesK0s_gap-10")->Clone("h3K0sEntries");
  // if(!h3K0sEntries) { Error("NotFound"); return; }
  // TH1D* h3K0sEntriesProjX = h3K0sEntries->ProjectionX("h3K0sEntriesProjX"); // whole projeciton
  // TH1D* h3K0sEntriesProjY = h3K0sEntries->ProjectionY("h3K0sEntriesProjY"); // whole projeciton
  TH1D* h3K0sEntriesProjZ = h3K0sEntries->ProjectionZ("h3K0sEntriesProjZ"); // whole projeciton
  // TH1D* h3K0sEntriesProjZsub = h3K0sEntries->ProjectionZ("h3K0sEntriesProjZsub",10,12,10,12); // projection in certain cent & pt range
  //
  // TCanvas* cK0s = new TCanvas("cK0s","cK0s",1000,1000);
  // cK0s->Divide(2,3);
  // cK0s->cd(1);
  // h3K0sEntries->Draw();
  // cK0s->cd(2);
  // h3K0sEntriesProjZ->Draw();
  // cK0s->cd(3);
  // h3K0sEntriesProjX->Draw();
  // cK0s->cd(4);
  // h3K0sEntriesProjY->Draw();
  // cK0s->cd(5);
  // h3K0sEntriesProjZsub->Draw();
  // cK0s->cd(6);
  // NOTE seems to work properly

  // correlations
  // TProfile3D* p3K0sCor = (TProfile3D*) lV0s->FindObject("fp3V0sCorrK0s_<2>_harm2_gap-10")->Clone("p3K0sCor");
  TProfile3D* p3K0sCor = (TProfile3D*) lV0s->FindObject("fp3PhiCorr_<2>_harm2_gap00")->Clone("p3K0sCor");
  if(!p3K0sCor) { Error("NotFound"); return; }

  p3K0sCor->GetXaxis()->SetRange(20,50);
  // p3K0sCor->GetZaxis()->SetRange(1,10);

  TProfile2D* profROOT = p3K0sCor->Project3DProfile("yz");
  TProfile2D* prof = Project3DProfile(p3K0sCor);

  TProfile* profROOTy = profROOT->ProfileX("profROOTy");
  TProfile* profy = prof->ProfileX("profy",40,70);
  // profy->Rebin(30);

  // TCanvas* cK0sCor = new TCanvas("cK0sCor","cK0sCor",1000,1000);
  // cK0sCor->Divide(2,3);
  // cK0sCor->cd(1);
  // p3K0sCor->Draw();
  // cK0sCor->cd(2);
  // h3K0sEntriesProjZ->Draw();
  // cK0sCor->cd(3);
  // profROOT->Draw("colz");
  // cK0sCor->cd(4);
  // prof->Draw("colz");
  // cK0sCor->cd(5);
  // profROOTy->Draw();
  // cK0sCor->cd(6);
  // profy->Draw();




  return;
}
//_____________________________________________________________________________
TProfile2D* ProcessUniFlow::Project3DProfile(const TProfile3D* prof3dorig)
{
  if(!prof3dorig) return 0x0;
  TProfile3D* prof3d = (TProfile3D*) prof3dorig->Clone();

  Int_t iBinFirst = prof3d->GetXaxis()->GetFirst();
  Int_t iBinLast = prof3d->GetXaxis()->GetLast();
  Int_t iNumBins = prof3d->GetNbinsX();
  Int_t iNumBinsAxis = prof3d->GetXaxis()->GetNbins();
  // printf("Bins:  %d - %d (%d | %d) \n", iBinFirst,iBinLast,iNumBins,iNumBinsAxis);

  // // making 3d hist from 3d profile
  // TH3D* hist3d = prof3d->ProjectionXYZ();   //NOTE do not care about range !!!
  // TH3D* hist3d_entry = prof3d->ProjectionXYZ("hist3d_entry","B");   //NOTE do not care about range !!!
  // TH3D* hist3d_weight = prof3d->ProjectionXYZ("hist3d_weight","W");   //NOTE do not care about range !!!

  TProfile2D* prof2d_test = DoProjectProfile2D(prof3d,"prof2d_test","",prof3d->GetYaxis(),prof3d->GetZaxis(),1,0,0);
  // prof2d_test->Draw("colz");

  return prof2d_test;

  // resulting profile
  // TProfile* result = new TProfile("result","result",100,0,10);
  // for(Int_t i(0); i < 10; i++) result->Fill(i,1);

  //
  // TCanvas* canTest = new TCanvas("canTest");
  // canTest->Divide(2,2);
  // canTest->cd(1);
  // prof3d->Draw("box");
  // canTest->cd(2);
  // hist3d->Draw("box");
  // canTest->cd(3);
  // hist3d_entry->Draw("box");
  // canTest->cd(4);
  // hist3d_weight->Draw("box");


  // return result;
}
//_____________________________________________________________________________
TProfile2D * ProcessUniFlow::DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,
                                           bool originalRange, bool useUF, bool useOF) const
{
// internal method to project to a 2D Profile
 // called from TH3::Project3DProfile but re-implemented in case of the TPRofile3D since what is done is different

 // projX, projY: axes of the orifinal histogram to which the projection is done (e.g. xy)

 // Get the ranges where we will work.
 Int_t ixmin = projX->GetFirst();
 Int_t ixmax = projX->GetLast();
 Int_t iymin = projY->GetFirst();
 Int_t iymax = projY->GetLast();
 if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
 if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
 Int_t nx = ixmax-ixmin+1;
 Int_t ny = iymax-iymin+1;

 // Create the projected profiles
 TProfile2D *p2 = 0;
 // Create always a new TProfile2D (not as in the case of TH3 projection)

 const TArrayD *xbins = projX->GetXbins();
 const TArrayD *ybins = projY->GetXbins();
 // assume all axis have variable bins or have fixed bins
 if ( originalRange ) {
    if (xbins->fN == 0 && ybins->fN == 0) {
       p2 = new TProfile2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                           ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
    } else {
       p2 = new TProfile2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
    }
 } else {
    if (xbins->fN == 0 && ybins->fN == 0) {
       p2 = new TProfile2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                           ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
    } else {
       p2 = new TProfile2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
    }
 }

 // new profile p2 is set according to axis ranges (keeping originals or not)

 // weights
 bool useWeights = (h3->fBinSumw2.fN != 0); //array elements
 if (useWeights) p2->Sumw2();

 // make projection in a 3D first // from 3D profile -> TH3
 TH3D * h3dW = h3->ProjectionXYZ("h3temp-W","W"); // getbincontent*getBinentries
 TH3D * h3dN = h3->ProjectionXYZ("h3temp-N","B"); // bin content is original profile = GetEntriesBin

 // fix ???
 h3dW->GetXaxis()->SetRange(h3->GetXaxis()->GetFirst(),h3->GetXaxis()->GetLast());
 h3dW->GetYaxis()->SetRange(h3->GetYaxis()->GetFirst(),h3->GetYaxis()->GetLast());
 h3dW->GetZaxis()->SetRange(h3->GetZaxis()->GetFirst(),h3->GetZaxis()->GetLast());
 h3dN->GetXaxis()->SetRange(h3->GetXaxis()->GetFirst(),h3->GetXaxis()->GetLast());
 h3dN->GetYaxis()->SetRange(h3->GetYaxis()->GetFirst(),h3->GetYaxis()->GetLast());
 h3dN->GetZaxis()->SetRange(h3->GetZaxis()->GetFirst(),h3->GetZaxis()->GetLast());




 h3dW->SetDirectory(0); h3dN->SetDirectory(0); // istograms does not bellow to any directorz ???

 // note that h3dW is always a weighted histogram - so we need to compute error in the projection
 TAxis * projX_hW = h3dW->GetXaxis();
 TAxis * projX_hN = h3dN->GetXaxis();
 if (projX == h3->GetYaxis() ) {  projX_hW =  h3dW->GetYaxis();  projX_hN =  h3dN->GetYaxis(); }
 if (projX == h3->GetZaxis() ) {  projX_hW =  h3dW->GetZaxis();  projX_hN =  h3dN->GetZaxis(); }
 TAxis * projY_hW = h3dW->GetYaxis();
 TAxis * projY_hN = h3dN->GetYaxis();
 if (projY == h3->GetXaxis() ) {  projY_hW =  h3dW->GetXaxis();  projY_hN =  h3dN->GetXaxis(); }
 if (projY == h3->GetZaxis() ) {  projY_hW =  h3dW->GetZaxis();  projY_hN =  h3dN->GetZaxis(); }
 // checking the axes

 // TH3 -> TH2
 TH2D * h2W = DoProject2D(h3dW,"htemp-W","",projX_hW, projY_hW, true, originalRange, useUF, useOF);
 TH2D * h2N = DoProject2D(h3dN,"htemp-N","",projX_hN, projY_hN, useWeights, originalRange, useUF, useOF);
 h2W->SetDirectory(0); h2N->SetDirectory(0);


 // fill the bin content
 R__ASSERT( h2W->fN == p2->fN );
 R__ASSERT( h2N->fN == p2->fN );
 R__ASSERT( h2W->GetSumw2()->fN != 0); // h2W should always be a weighted histogram since h3dW is weighted


 // filling the new tprofile2D
 for (int i = 0; i < p2->fN ; ++i) {
    //std::cout << " proj bin " << i << "  " <<  h2W->fArray[i] << "  " << h2N->fArray[i] << std::endl;
    p2->fArray[i] = h2W->fArray[i];   // array of profile is sum of all values
    p2->GetSumw2()->fArray[i]  = h2W->GetSumw2()->fArray[i];   // array of content square of profile is weight square of the W projected histogram
    p2->SetBinEntries(i, h2N->fArray[i] );
    if (useWeights) p2->GetBinSumw2()->fArray[i] = h2N->GetSumw2()->fArray[i];    // sum of weight squares are stored to compute errors in h1N histogram
 }
 // delete the created histograms
 delete h3dW;
 delete h3dN;
 delete h2W;
 delete h2N;

 // Also we need to set the entries since they have not been correctly calculated during the projection
 // we can only set them to the effective entries
 p2->SetEntries( p2->GetEffectiveEntries() );

 return p2;
}
//_____________________________________________________________________________
TH2D* ProcessUniFlow::DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY,
                    bool computeErrors, bool originalRange,
                    bool useUF, bool useOF) const
{
  // internal method performing the projection to a 2D histogram
     // called from TH3::Project3D

     TH2D *h2 = 0;

     // Get range to use as well as bin limits
     Int_t ixmin = projX->GetFirst();
     Int_t ixmax = projX->GetLast();
     Int_t iymin = projY->GetFirst();
     Int_t iymax = projY->GetLast();
     if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
     if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
     Int_t nx = ixmax-ixmin+1;
     Int_t ny = iymax-iymin+1;

      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange )
      {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1]
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
         }
      } else {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1]
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
         }
      }

    //  // Copy the axis attributes and the axis labels if needed.
    //  THashList* labels1 = 0;
    //  THashList* labels2 = 0;
    //  // "xy"
    //  h2->GetXaxis()->ImportAttributes(projY);
    //  h2->GetYaxis()->ImportAttributes(projX);
    //  labels1 = projY->GetLabels();
    //  labels2 = projX->GetLabels();
    //  if (labels1) {
    //     TIter iL(labels1);
    //     TObjString* lb;
    //     Int_t i = 1;
    //     while ((lb=(TObjString*)iL())) {
    //        h2->GetXaxis()->SetBinLabel(i,lb->String().Data());
    //        i++;
    //     }
    //  }
    //  if (labels2) {
    //     TIter iL(labels2);
    //     TObjString* lb;
    //     Int_t i = 1;
    //     while ((lb=(TObjString*)iL())) {
    //        h2->GetYaxis()->SetBinLabel(i,lb->String().Data());
    //        i++;
    //     }
    //  }
    //  h2->SetLineColor(this->GetLineColor());
    //  h2->SetFillColor(this->GetFillColor());
    //  h2->SetMarkerColor(this->GetMarkerColor());
    //  h2->SetMarkerStyle(this->GetMarkerStyle());

     // Activate errors
     if ( computeErrors) h2->Sumw2();

     // Set references to the axis, so that the bucle has no branches.
     TAxis* out = 0;
     if ( projX != h3->GetXaxis() && projY != h3->GetXaxis() ) {
        out = h3->GetXaxis();
     } else if ( projX != h3->GetYaxis() && projY != h3->GetYaxis() ) {
        out = h3->GetYaxis();
     } else {
        out = h3->GetZaxis();
     }

     Int_t *refX = 0, *refY = 0, *refZ = 0;
     Int_t ixbin, iybin, outbin;
     if ( projX == h3->GetXaxis() && projY == h3->GetYaxis() ) { refX = &ixbin;  refY = &iybin;  refZ = &outbin; }
     if ( projX == h3->GetYaxis() && projY == h3->GetXaxis() ) { refX = &iybin;  refY = &ixbin;  refZ = &outbin; }
     if ( projX == h3->GetXaxis() && projY == h3->GetZaxis() ) { refX = &ixbin;  refY = &outbin; refZ = &iybin;  }
     if ( projX == h3->GetZaxis() && projY == h3->GetXaxis() ) { refX = &iybin;  refY = &outbin; refZ = &ixbin;  }
     if ( projX == h3->GetYaxis() && projY == h3->GetZaxis() ) { refX = &outbin; refY = &ixbin;  refZ = &iybin;  }
     if ( projX == h3->GetZaxis() && projY == h3->GetYaxis() ) { refX = &outbin; refY = &iybin;  refZ = &ixbin;  }
     R__ASSERT (refX != 0 && refY != 0 && refZ != 0);

     // Fill the projected histogram excluding underflow/overflows if considered in the option
     // if specified in the option (by default they considered)
     Double_t totcont  = 0;

     Int_t outmin = out->GetFirst();
     Int_t outmax = out->GetLast();
     // GetFirst(), GetLast() can return (0,0) when the range bit is set artifically (see TAxis::SetRange)
     if (outmin == 0 && outmax == 0) { outmin = 1; outmax = out->GetNbins(); }
     // correct for underflow/overflows
     if (useUF && !out->TestBit(TAxis::kAxisRange) )  outmin -= 1;
     if (useOF && !out->TestBit(TAxis::kAxisRange) )  outmax += 1;

     for (ixbin=0;ixbin<=1+projX->GetNbins();ixbin++){
        if ( projX->TestBit(TAxis::kAxisRange) && ( ixbin < ixmin || ixbin > ixmax )) continue;
        Int_t ix = h2->GetYaxis()->FindBin( projX->GetBinCenter(ixbin) );

        for (iybin=0;iybin<=1+projY->GetNbins();iybin++){
           if ( projY->TestBit(TAxis::kAxisRange) && ( iybin < iymin || iybin > iymax )) continue;
           Int_t iy = h2->GetXaxis()->FindBin( projY->GetBinCenter(iybin) );

           Double_t cont = 0;
           Double_t err2 = 0;

           // loop on the bins to be integrated (outbin should be called inbin)
           for (outbin = outmin; outbin <= outmax; outbin++){

              Int_t bin = h3->GetBin(*refX,*refY,*refZ);

              // sum the bin contents and errors if needed
              cont += h3->GetBinContent(bin);
              if (computeErrors) {
                 Double_t exyz = h3->GetBinError(bin);
                 err2 += exyz*exyz;
              }

           }

           // remember axis are inverted
           h2->SetBinContent(iy , ix, cont);
           if (computeErrors) h2->SetBinError(iy, ix, TMath::Sqrt(err2) );
           // sum all content
           totcont += cont;

        }
     }

     // since we use fill we need to reset and recalculate the statistics (see comment in DoProject1D )
     // or keep original statistics if consistent sumw2
     bool resetStats = true;
     double eps = 1.E-12;
     if (h3->IsA() == TH3F::Class() ) eps = 1.E-6;
     if (h3->fTsumw != 0 && TMath::Abs( h3->fTsumw - totcont) <  TMath::Abs(h3->fTsumw) * eps) resetStats = false;

     bool resetEntries = resetStats;
     // entries are calculated using underflow/overflow. If excluded entries must be reset
     resetEntries |= !useUF || !useOF;

     if (!resetStats) {
        Double_t stats[TH1::kNstat];
        Double_t oldst[TH1::kNstat]; // old statistics
        for (Int_t i = 0; i < TH1::kNstat; ++i) { oldst[i] = 0; }
        h3->GetStats(oldst);
        std::copy(oldst,oldst+TH1::kNstat,stats);
        // not that projX refer to Y axis and projX refer to the X axis of projected histogram
        // nothing to do for projection in Y vs X
        if ( projY == h3->GetXaxis() && projX == h3->GetZaxis() ) {  // case XZ
           stats[4] = oldst[7];
           stats[5] = oldst[8];
           stats[6] = oldst[9];
        }
        if ( projY == h3->GetYaxis() ) {
           stats[2] = oldst[4];
           stats[3] = oldst[5];
           if ( projX == h3->GetXaxis() )  { // case YX
              stats[4] = oldst[2];
              stats[5] = oldst[3];
           }
           if ( projX == h3->GetZaxis() )  { // case YZ
              stats[4] = oldst[7];
              stats[5] = oldst[8];
              stats[6] = oldst[10];
           }
        }
        else if  ( projY == h3->GetZaxis() ) {
           stats[2] = oldst[7];
           stats[3] = oldst[8];
           if ( projX == h3->GetXaxis() )  { // case ZX
              stats[4] = oldst[2];
              stats[5] = oldst[3];
              stats[6] = oldst[9];
           }
           if ( projX == h3->GetYaxis() )  { // case ZY
              stats[4] = oldst[4];
              stats[5] = oldst[5];
              stats[6] = oldst[10];
           }
        }
        // set the new statistics
        h2->PutStats(stats);
     }
     else {
        // recalculate the statistics
        h2->ResetStats();
     }

     if (resetEntries) {
        // use the effective entries for the entries
        // since this  is the only way to estimate them
        Double_t entries =  h2->GetEffectiveEntries();
        if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
        h2->SetEntries( entries );
     }
     else {
        h2->SetEntries( h3->fEntries );
     }


     return h2;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowOneGo"); return kFALSE; }
  if(listFits->GetEntries() > 0) { Error("TList not empty!","ExtractFlowOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowOneGo"); return kFALSE; }
  if(task->fSpecies == FlowTask::kPhi && task->fFlowFitPhiSubtLS && !hInvMassBG) { Error("Inv. Mass (BG) histogram does not exists!","ExtractFlowOneGo"); return kFALSE; }
  if(task->fFlowFitPhiSubtLS) { Error("Phi like-sign subtraction not implemented ATM. Please turn the switch off.","ExtractFlowOneGo"); return kFALSE; }

  // === Fitting parametrisation (species dependent default) ===

  Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  Double_t dMaximum = hInvMass->GetMaximum();

  Int_t iNpx = 10000;
  TString sFitOptMass = "RNL";
  TString sFitOptFlow = "RN";

  TString sMassBG = TString(); Int_t iNumParsMassBG = 0; // function for inv. mass dist. (BG component)
  TString sMassSig = TString();  Int_t iNumParsMassSig = 0; // function for inv. mass dist. (sig component)
  TString sFlowBG = TString();  Int_t iNumParsFlowBG = 0; // function for flow-mass (BG component)

  Int_t iParMass = 0;
  Int_t iParWidth = 0;
  Int_t iParWidth_2 = 0;

  std::vector<Double_t> dParDef;
  std::vector<Double_t> dParLimLow;
  std::vector<Double_t> dParLimHigh;

  if(task->fSpecies == FlowTask::kPhi)
  {
    dMassRangeLow = 0.994;
    // dMassRangeHigh = 1.134;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::BreitWigner(x,[5],[6])"; iNumParsMassSig = 3;
    sFlowBG = "[7]*x+[8]"; iNumParsFlowBG = 2;

    iParMass = 5;
    iParWidth = 6;

    Double_t dDef[] =      {1.0,1.0,1.0,1.0,   dMaximum,1.019445,0.0046, 1.0,1.0};
    Double_t dLimLow[] =   {-1,-1,-1,-1,    0.0,1.018,0.001, -1,-1};
    Double_t dLimHigh[] =  {-1,-1,-1,-1,  2.0*dMaximum,1.022,0.006,  -1,-1};

    // assignment to external arrays
    for(Int_t par(0); par < (Int_t) (sizeof(dDef)/sizeof(dDef[0])); ++par) { dParDef.push_back(dDef[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimLow)/sizeof(dLimLow[0])); ++par) { dParLimLow.push_back(dLimLow[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimHigh)/sizeof(dLimHigh[0])); ++par) { dParLimHigh.push_back(dLimHigh[par]); }
  }
  else if(task->fSpecies == FlowTask::kK0s)
  {
    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;
    sFlowBG = "[9]*x+[10]"; iNumParsFlowBG = 2;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;

    Double_t dDef[] =      {1.0,1.0,1.0,1.0,   dMaximum,0.4976,0.003,dMaximum,0.01, 1.0,1.0};
    Double_t dLimLow[] =   {-1,-1,-1,-1,    0.0,0.48,0.003,0.0,0.003, -1,-1};
    Double_t dLimHigh[] =  {-1,-1,-1,-1,  2.0*dMaximum,0.52,0.006,2.0*dMaximum,0.01,  -1,-1};

    // assignment to external arrays
    for(Int_t par(0); par < (Int_t) (sizeof(dDef)/sizeof(dDef[0])); ++par) { dParDef.push_back(dDef[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimLow)/sizeof(dLimLow[0])); ++par) { dParLimLow.push_back(dLimLow[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimHigh)/sizeof(dLimHigh[0])); ++par) { dParLimHigh.push_back(dLimHigh[par]); }
  }
  else if(task->fSpecies == FlowTask::kLambda)
  {
    dMassRangeLow = 1.096;
    dMassRangeHigh = 1.150;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;
    sFlowBG = "[9]*x+[10]"; iNumParsFlowBG = 2;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;

    Double_t dDef[] =      {1.0,1.0,1.0,1.0,   dMaximum,1.115, 0.001,dMaximum,0.01, 1.0,1.0};
    Double_t dLimLow[] =   {-1,-1,-1,-1,    0.0,1.10,0.001,0.0,0.001, -1,-1};
    Double_t dLimHigh[] =  {-1,-1,-1,-1,  2.0*dMaximum,1.13,0.008,2.0*dMaximum,0.01,  -1,-1};

    // assignment to external arrays
    for(Int_t par(0); par < (Int_t) (sizeof(dDef)/sizeof(dDef[0])); ++par) { dParDef.push_back(dDef[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimLow)/sizeof(dLimLow[0])); ++par) { dParLimLow.push_back(dLimLow[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimHigh)/sizeof(dLimHigh[0])); ++par) { dParLimHigh.push_back(dLimHigh[par]); }
  }
  else { Error("Invalid species","ExtractFlowOneGo"); return kFALSE; }

  // check if parametrisation is setup manually
  if(task->fFlowFitRangeLow > 0.0) { dMassRangeLow = task->fFlowFitRangeLow; }
  if(task->fFlowFitRangeHigh > 0.0) { dMassRangeLow = task->fFlowFitRangeHigh; }
  if(task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow >= task->fFlowFitRangeHigh) { Error("Wrong fitting ranges set!","ExtractFlowOneGo"); return kFALSE; }

  if(task->fNumParMassSig > 0) { sMassSig = task->fFlowFitMassSig; iNumParsMassSig = task->fNumParMassSig; Debug(" Task massSig set","ExtractFlowOneGo"); }
  if(task->fNumParMassBG > 0) { sMassBG = task->fFlowFitMassBG; iNumParsMassBG = task->fNumParMassBG; Debug(" Task massBG set","ExtractFlowOneGo"); }
  if(task->fNumParFlowBG > 0) { sFlowBG = task->fFlowFitFlowBG; iNumParsFlowBG = task->fNumParFlowBG; Debug(" Task flowBG set","ExtractFlowOneGo"); }

  // TODO paramatetrisaion

  Int_t iNumParDefs = dParDef.size();
  Int_t iNumParLimLow = dParLimLow.size();
  Int_t iNumParLimHigh = dParLimHigh.size();

  // check the output of the vector assigment
  if(fbDebug)
  {
    Debug("Post ifs","ExtractFlowOneGo");
    for(Int_t par(0); par < iNumParDefs; ++par) { printf("  par %d: %g (%g<%g)\n",par, dParDef.at(par), dParLimLow.at(par), dParLimHigh.at(par)); }
  }

  if(iNumParDefs != iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG) { Error(Form("Length of dParDef array does not match number of parameters (%d != %d)",iNumParDefs,iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG),"ExtractFlowOneGo"); return kFALSE; }
  if(iNumParDefs != iNumParLimLow) { Error(Form("Different length of arrays with parameter defauls and low limit values (%d != %d).",iNumParDefs,iNumParLimLow),"ExtractFlowOneGo"); return kFALSE; }
  if(iNumParDefs != iNumParLimHigh) { Error(Form("Different length of arrays with parameter defauls and high limit values (%d != %d).",iNumParDefs,iNumParLimHigh),"ExtractFlowOneGo"); return kFALSE; }

  // === Initialision ===
  Int_t iNumParMass = iNumParsMassSig+iNumParsMassBG;
  Int_t iParFlow = iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG; // index of Flow (vn/dn) parameter

  // master formula used in the fitting procedure
  if(!fbDebug) { sFitOptFlow += "Q"; sFitOptMass += "Q"; } // quite fitting option if NOT in debug

  TString sFuncMass = Form("%s + %s",sMassBG.Data(),sMassSig.Data());
  TString sFuncVn = Form("[%d]*(%s)/(%s + %s) + (%s)*(%s)/(%s + %s)", iParFlow, sMassSig.Data(), sMassSig.Data(), sMassBG.Data(), sFlowBG.Data(),sMassBG.Data(),sMassSig.Data(),sMassBG.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh), "ExtractFlowOneGo");
  Debug(Form("Fit Dist :\n    %s",sFuncMass.Data()), "ExtractFlowOneGo");
  Debug(Form("Fit Flow :\n    %s\n",sFuncVn.Data()), "ExtractFlowOneGo");

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // === Fitting procedure ===

  // fitting invariant mass distribution
  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetNpx(iNpx);

  for(Int_t par(0); par < iNumParMass; ++par)
  {
    fitMass->SetParameter(par, dParDef.at(par));

    Double_t dLimLow = dParLimLow.at(par);
    Double_t dLimHigh = dParLimHigh.at(par);

    if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
    else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"ExtractFlowOneGo"); return kFALSE; }
  }

  hInvMass->Fit(fitMass, sFitOptMass.Data());

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 15))
  {
    fitMass->SetParameter(0, fitMass->GetParameter(0)/nfitsA);
    for(Int_t par(0); par < iNumParMass; ++par)
    {
      fitMass->SetParameter(par, fitMass->GetParameter(par));

      Double_t dLimLow = dParLimLow.at(par);
      Double_t dLimHigh = dParLimHigh.at(par);

      if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
      else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"ExtractFlowOneGo"); return kFALSE; }
    }

    hInvMass->Fit(fitMass, sFitOptMass.Data());

    statusA = gMinuit->fCstatu.Data();
    nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv.mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("\n Inv.mass distribution fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g; %d iterations)\n",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF(),fitMass->GetProb(),nfitsA), "ExtractFlowOneGo");

  // fitting invariant mass distribution
  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  // fixing Nsig & Nbg terms extracted in previous step
  for(Int_t par(0); par < iNumParMass; ++par) { fitVn->FixParameter(par, fitMass->GetParameter(par)); }
  for(Int_t par(iNumParMass); par < iParFlow; ++par)
  {
    fitVn->SetParameter(par, dParDef.at(par));

    Double_t dLimLow = dParLimLow.at(par);
    Double_t dLimHigh = dParLimHigh.at(par);

    if(dLimLow > -1.0 && dLimHigh > -1.0) { fitVn->SetParLimits(par, dLimLow,dLimHigh); }
    else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Flow-mass: Only one of the parameter limits is set (par %d). Fix this!",par),"ExtractFlowOneGo"); return kFALSE; }
  }
  hFlowMass->Fit(fitVn, sFitOptFlow.Data());

  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Flow-mass fit does not converged!"), "ExtractFlowOneGo"); return kFALSE; }
  Info(Form("\nFlow-mass fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g)\n",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF(),fitVn->GetProb()), "ExtractFlowOneGo");

  // saving flow to output
  dFlow = fitVn->GetParameter(iParFlow);
  dFlowError = fitVn->GetParError(iParFlow);
  Info(Form("Final flow: %g +- %g\n=================================\n", dFlow,dFlowError), "ExtractFlowOneGo");

  // === Extracting fitting components to separated TF1's ===

  TF1* fitBg = new TF1("fitMassBG",sMassBG.Data(),dMassRangeLow,dMassRangeHigh);
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar)
  {
    fitBg->SetParameter(iPar, fitMass->GetParameter(iPar));
    fitBg->SetParError(iPar, fitMass->GetParError(iPar));
  }

  TF1* fitSig = new TF1("fitMassSig", sMassSig.Data(), dMassRangeLow,dMassRangeHigh);
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsMassBG); iPar < iNumParsMassBG+iNumParsMassSig; ++iPar)
  {
    fitSig->SetParameter(iPar, fitMass->GetParameter(iPar));
    fitSig->SetParError(iPar, fitMass->GetParError(iPar));
  }

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("(%s)*(%s)/(%s + %s)", sFlowBG.Data(), sMassBG.Data(), sMassSig.Data(), sMassBG.Data()), dMassRangeLow,dMassRangeHigh);
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  for(Int_t iPar(0); iPar < iParFlow; ++iPar)
  {
    fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar));
    fitFlowBg->SetParError(iPar, fitVn->GetParError(iPar));
  }

  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[%d]*(%s)/(%s + %s)", iParFlow, sMassSig.Data(), sMassSig.Data(), sMassBG.Data()), dMassRangeLow,dMassRangeHigh);
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);
  fitFlowSig->SetParameter(iParFlow, fitVn->GetParameter(iParFlow));
  for(Int_t iPar(0); iPar < iNumParsMassBG+iNumParsMassSig; ++iPar)
  {
    fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar));
    fitFlowSig->SetParError(iPar, fitVn->GetParError(iPar));
  }
  for(Int_t iPar(iNumParsMassBG+iNumParsMassSig); iPar < iParFlow; ++iPar) { fitFlowSig->SetParameter(iPar, 0.0); }

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);

  // === Drawing stuff to canvas ===

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  latex->SetNDC();

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hInvMass->SetMarkerStyle(20);
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitBg->DrawCopy("same");
  fitSig->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(iParMass),fitMass->GetParError(iParMass)));

  if(task->fSpecies == FlowTask::kPhi)
  {
    latex->DrawLatex(0.17,0.70,Form("#color[9]{#Gamma = %.6f #pm %.6f}",fitMass->GetParameter(iParWidth),fitMass->GetParError(iParWidth)));
  }
  else if(task->fSpecies == FlowTask::kK0s || task->fSpecies == FlowTask::kLambda)
  {
    latex->DrawLatex(0.17,0.70,Form("#color[9]{#sigma_{1} = %.6f #pm %.6f}",fitMass->GetParameter(iParWidth),fitMass->GetParError(iParWidth)));
    latex->DrawLatex(0.17,0.65,Form("#color[9]{#sigma_{2} = %.6f #pm %.6f}",fitMass->GetParameter(iParWidth_2),fitMass->GetParError(iParWidth_2)));
  }

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hFlowMass->SetMarkerStyle(20);
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  // fitFlowSig->DrawCopy("same");
  // fitFlowBg->DrawCopy("same");

  TString sResult = "v_{2}"; if(fFlowFitCumulants) { sResult = "d_{2}"; }
  latex->DrawLatex(0.17,0.80,Form("#color[9]{%s = %.4f #pm %.4f}",sResult.Data(),dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{p = %.3g}",fitVn->GetProb()));

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowPhiOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowPhiOneGo"); return kFALSE; }

  // subtraction of LS?

  // fitting parametrisation
  TString sFuncBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; Int_t iNumParsFuncBG = 4;
  TString sFuncSig = "[4]*TMath::BreitWigner(x,[5],[6])";  Int_t iNumParsFuncSig = 3;

  TString sFuncVn = Form("[9]*(%s)/(%s + %s) + ([7]*x+[8])*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data(), sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data());
  TString sFuncMass = Form("%s + %s",sFuncBG.Data(),sFuncSig.Data());

  // Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  Double_t dMassRangeLow = 0.994;
  // Double_t dMassRangeHigh =1.134;
  Double_t dMaximum = hInvMass->GetMaximum();

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh),"ExtractFlowPhiOneGo");
  Debug(Form("Fit func invmass:\n%s",sFuncVn.Data()));
  Debug(Form("Fit func vn:\n%s",sFuncVn.Data()));


  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetParameters(dMaximum/10.0, 1.0, 1.0, 1.0, dMaximum, 1.019445, 0.0046);
  fitMass->SetNpx(5000);
  fitMass->SetParLimits(5,1.018,1.022);
  fitMass->SetParLimits(6,0.001,0.006);
  hInvMass->Fit(fitMass, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 10))
  {
    fitMass->SetParameters(fitMass->GetParameter(0)/nfitsA, fitMass->GetParameter(1), fitMass->GetParameter(2), fitMass->GetParameter(3), fitMass->GetParameter(4), 1.019445, fitMass->GetParameter(5)*nfitsA);
    fitMass->SetParLimits(5,1.018,1.022);
    fitMass->SetParLimits(6,0.001,0.006);
    hInvMass->Fit(fitMass, "RNL");

    statusA = gMinuit->fCstatu.Data();
    nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Number of iterations: %d\n",nfitsA));

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  fitVn->SetParameter(7, 1.0);
  fitVn->SetParameter(8, 1.0);
  fitVn->SetParameter(9, 0.1);

  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitVn->FixParameter(iPar, fitMass->GetParameter(iPar)); }
  hFlowMass->Fit(fitVn, "RN");
  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged!")); return kFALSE; }

  // saving flow to output
  dFlow = fitVn->GetParameter(9);
  dFlowError = fitVn->GetParError(9);
  Info(Form("=================================\n Final flow: %g +- %g\n =================================\n", dFlow,dFlowError));

  TF1* fitBg = new TF1("fitMassBG",sFuncBG.Data(),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitBg->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);

  TF1* fitSig = new TF1("fitMassSig", sFuncSig.Data(), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsFuncBG); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitSig->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("([7]*x+[8])*(%s)/(%s + %s)",sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data()),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowBg->SetParameter(7, fitVn->GetParameter(7));
  fitFlowBg->SetParameter(8, fitVn->GetParameter(8));
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  //
  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[9]*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data()), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowSig->SetParameter(7, 0.0);
  fitFlowSig->SetParameter(8, 0.0);
  fitFlowSig->SetParameter(9, fitVn->GetParameter(9));
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);

  // Drawing stuff
  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hInvMass->SetMarkerStyle(20);
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitBg->DrawCopy("same");
  fitSig->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(5),fitMass->GetParError(5)));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{#Gamma = %.6f #pm %.6f}",fitMass->GetParameter(6),fitMass->GetParError(6)));

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hFlowMass->SetMarkerStyle(20);
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  // fitFlowSig->DrawCopy("same");
  // fitFlowBg->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{v_{2} = %.4f #pm %.4f}",dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowK0sOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowK0sOneGo"); return kFALSE; }


  Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  // Double_t dMassRangeLow = 0.41;
  // Double_t dMassRangeHigh = 0.59;
  Double_t dMaximum = hInvMass->GetMaximum();

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // fitting parametrisation
  TString sFuncBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; Int_t iNumParsFuncBG = 4;
  TString sFuncSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; Int_t iNumParsFuncSig = 5;

  TString sFuncVn = Form("[11]*(%s)/(%s + %s) + ([9]*x+[10])*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data(), sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data());
  TString sFuncMass = Form("%s + %s",sFuncBG.Data(),sFuncSig.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh),"ExtractFlowPhiOneGo");
  Debug(Form("Fit func invmass:\n%s",sFuncVn.Data()));
  Debug(Form("Fit func vn:\n%s",sFuncVn.Data()));

  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetNpx(10000);
  fitMass->SetParameters(dMaximum/10.0, 1.0, 1.0, 1.0, dMaximum, 0.4976, 0.001,dMaximum,0.001);
  fitMass->SetParLimits(4, 0, dMaximum*2.0);
  fitMass->SetParLimits(5, 0.48, 0.52);
  fitMass->SetParLimits(6, 0.003, 0.006);
  fitMass->SetParLimits(7, 0, dMaximum*2.0);
  fitMass->SetParLimits(8, 0.003, 0.01);
  hInvMass->Fit(fitMass, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 15)){

      fitMass->SetParameters(fitMass->GetParameter(0)/nfitsA, fitMass->GetParameter(1), fitMass->GetParameter(2), fitMass->GetParameter(3), fitMass->GetParameter(4), 0.4976, fitMass->GetParameter(6)*nfitsA,fitMass->GetParameter(7), fitMass->GetParameter(8)*nfitsA);
      // fitMass->SetParLimits(4, 0, dMaximum*2.0);
      fitMass->SetParLimits(5, 0.48, 0.52);
      fitMass->SetParLimits(6, 0.003, 0.006);
      // fitMass->SetParLimits(7, 0, dMaximum*2.0);
      fitMass->SetParLimits(8, 0.003, 0.01);
      hInvMass->Fit(fitMass, "RNL");

      statusA = gMinuit->fCstatu.Data();
      nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Number of iterations: %d\n",nfitsA));

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  fitVn->SetParameter(9, 1.0);
  fitVn->SetParameter(10, 1.0);
  fitVn->SetParameter(11, 0.1);

  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitVn->FixParameter(iPar, fitMass->GetParameter(iPar)); }
  hFlowMass->Fit(fitVn, "RN");
  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged!")); return kFALSE; }

  // saving flow to output
  dFlow = fitVn->GetParameter(11);
  dFlowError = fitVn->GetParError(11);
  Info(Form("=================================\n Final flow: %g +- %g\n =================================\n", dFlow,dFlowError));

  // Drawing stuff
  TF1* fitBg = new TF1("fitBG",sFuncBG.Data(),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitBg->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);

  TF1* fitSig = new TF1("fitSig", sFuncSig.Data(), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsFuncBG); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitSig->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("([7]*x+[8])*(%s)/(%s + %s)",sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data()),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowBg->SetParameter(9, fitVn->GetParameter(9));
  fitFlowBg->SetParameter(10, fitVn->GetParameter(10));
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  //
  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[9]*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data()), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowSig->SetParameter(9, 0.0);
  fitFlowSig->SetParameter(10, 0.0);
  fitFlowSig->SetParameter(11, fitVn->GetParameter(11));
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);


  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitSig->DrawCopy("same");
  fitBg->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(5),fitMass->GetParError(5)));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{#sigma_{1} = %.6f #pm %.6f}",fitMass->GetParameter(6),fitMass->GetParError(6)));
  latex->DrawLatex(0.17,0.65,Form("#color[9]{#sigma_{2} = %.6f #pm %.6f}",fitMass->GetParameter(8),fitMass->GetParError(8)));

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{K^{0}} (GeV/c^{2})");
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{v_{2} = %.4f #pm %.4f}",dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowLambdaOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowLambdaOneGo"); return kFALSE; }


  TString sFuncBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; Int_t iNumParsFuncBG = 4;
  TString sFuncSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; Int_t iNumParsFuncSig = 5;
  // TString sFuncSig = "gaus(4)+gaus(7)"; Int_t iNumParsFuncSig = 6;

  TString sFuncMass = Form("%s + %s",sFuncBG.Data(),sFuncSig.Data());
  TString sFuncVn = Form("[11]*(%s)/(%s + %s) + ([9]*x+[10])*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data(), sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data());

  // Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  // Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  Double_t dMassRangeLow = 1.096;
  Double_t dMassRangeHigh =1.150;
  Double_t dMaximum = hInvMass->GetMaximum();

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);


  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh),"ExtractFlowPhiOneGo");
  Debug(Form("Fit func invmass:\n%s",sFuncVn.Data()));
  Debug(Form("Fit func vn:\n%s",sFuncVn.Data()));

  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetParameters(dMaximum/10.0, 1.0, 1.0, 1.0, dMaximum, 1.115, 0.001,dMaximum, 0.001);
  fitMass->SetNpx(5000);
  fitMass->SetParLimits(4, 0, dMaximum*2.0);
  fitMass->SetParLimits(5, 1.10, 1.13);
  fitMass->SetParLimits(6, 0.001,0.006);
  fitMass->SetParLimits(7, 0, dMaximum*2.0);
  fitMass->SetParLimits(8, 0.001,0.01);
  hInvMass->Fit(fitMass, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 15))
  {
    fitMass->SetParameters(fitMass->GetParameter(0)/nfitsA, fitMass->GetParameter(1), fitMass->GetParameter(2), fitMass->GetParameter(3), fitMass->GetParameter(4), 1.115, fitMass->GetParameter(6)*nfitsA, fitMass->GetParameter(7), fitMass->GetParameter(9)*nfitsA);
    fitMass->SetParLimits(4, 0, dMaximum*2.0);
    fitMass->SetParLimits(5, 1.10, 1.13);
    fitMass->SetParLimits(6, 0.001,0.008);
    fitMass->SetParLimits(7, 0, dMaximum*2.0);
    fitMass->SetParLimits(8, 0.001,0.01);
    hInvMass->Fit(fitMass, "RNL");

    statusA = gMinuit->fCstatu.Data();
    nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Number of iterations: %d\n",nfitsA));

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  fitVn->SetParameter(9, 1.0);
  fitVn->SetParameter(10, 1.0);
  fitVn->SetParameter(11, 0.1);

  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitVn->FixParameter(iPar, fitMass->GetParameter(iPar)); }
  hFlowMass->Fit(fitVn, "RN");
  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged!")); return kFALSE; }

  // saving flow to output
  dFlow = fitVn->GetParameter(11);
  dFlowError = fitVn->GetParError(11);
  Info(Form("=================================\n Final flow: %g +- %g\n =================================\n", dFlow,dFlowError));

  // Drawing stuff
  TF1* fitBg = new TF1("fitBG",sFuncBG.Data(),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitBg->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);

  TF1* fitSig = new TF1("fitSig", sFuncSig.Data(), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsFuncBG); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitSig->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("([9]*x+[10])*(%s)/(%s + %s)",sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data()),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowBg->SetParameter(9, fitVn->GetParameter(9));
  fitFlowBg->SetParameter(10, fitVn->GetParameter(10));
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  //
  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[12]*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data()), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  // fitFlowSig->SetParameter(10, 0.0);
  // fitFlowSig->SetParameter(11, 0.0);
  fitFlowSig->SetParameter(11, fitVn->GetParameter(11));
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->GetXaxis()->SetTitle("M_{#Lambda} (GeV/c^{2})");
  hInvMass->SetMarkerStyle(20);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitSig->DrawCopy("same");
  fitBg->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(5),fitMass->GetParError(5)));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{#sigma_{1} = %.6f #pm %.6f}",fitMass->GetParameter(6),fitMass->GetParError(6)));
  latex->DrawLatex(0.17,0.65,Form("#color[9]{#sigma_{2} = %.6f #pm %.6f}",fitMass->GetParameter(8),fitMass->GetParError(8)));

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{K^{0}} (GeV/c^{2})");
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{v_{2} = %.4f #pm %.4f}",dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));

  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::Fatal(TString sMsg, TString sMethod)
{
	printf("Fatal::%s  %s. Terminating!\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Error(TString sMsg, TString sMethod)
{
	printf("Error::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Info(TString sMsg, TString sMethod)
{
	printf("Info::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Warning(TString sMsg, TString sMethod)
{
	printf("Warning::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Debug(TString sMsg, TString sMethod)
{
	if(fbDebug) printf("Debug::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}
//_____________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
//_____________________________________________________________________________
void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2();
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2();
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

}
//_____________________________________________________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}
//_____________________________________________________________________________
void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"
  TString sText;

  switch (logo)
  {
    case 0:
      sText = TString("ALICE");
      break;

    case 1:
      sText = TString("ALICE Preliminary");
      break;

    case 2:
      sText = TString("ALICE This thesis");
      break;
  }


  TLatex *   tex = new TLatex(xmin,ymin, sText.Data());
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();

  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}//_____________________________________________________________________________
