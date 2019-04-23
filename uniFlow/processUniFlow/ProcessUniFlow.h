#ifndef PROCESSUNIFLOW_H
#define PROCESSUNIFLOW_H

#include <vector>
#include "TString.h"

class FlowTask;

class TFile;
class TF1;
class TAxis;
class TCanvas;
class THnSparse;
class TH1;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
class TProfile3D;
class TFitResultPtr;

enum        PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest
enum        Cumulants {kNon = 0, kTwo = 2, kFour = 4}; // Cumulants order

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
    void        SetMultiplicityBins(std::vector<Double_t> array) { fdMultBins = array; fiNumMultBins = (Int_t) array.size() - 1; } // setup the global multiplicity binning, where size is number of elements in array
    void        SetFitCumulants(Bool_t cum = kTRUE) { fFlowFitCumulants = cum; } // use cn{2} vs m_inv instead of vn{2} vs. m_inv
    void        SetSaveInterSteps(Bool_t save = kTRUE) { fSaveInterSteps = save; }
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }

    static TString     GetSpeciesName(PartSpecies species); // system species name (Charged, K0s, ...)
    static TString     GetSpeciesLabel(PartSpecies species); // readable species name (h^{#pm}, K^{0}_{S}, ...)
    static Bool_t      IsSpeciesDirect(PartSpecies sp) { return (sp == kCharged || sp == kPion || sp == kKaon || sp == kProton); } //
    static Bool_t      IsSpeciesReconstructed(PartSpecies sp) { return (sp == kK0s || sp == kLambda || sp == kPhi); } //

    void        AddTask(FlowTask* task = 0x0); // add task to internal lists of all tasks
    void        Run(); // running the task (main body of the class)
    void        Clear(); // clearing (removing tasks, etc.) after running

    // printing output methods
    static void        Fatal(TString sMsg, TString sMethod = ""); // printf the msg as error
    static void        Error(TString sMsg, TString sMethod = ""); // printf the msg as error
    static void        Warning(TString sMsg, TString sMethod = ""); // printf the msg as warning
    static void        Info(TString sMsg, TString sMethod = ""); // printf the msg as info
    void               Debug(TString sMsg, TString sMethod = ""); // printf the msg as info
  protected:

  private:
    Bool_t      Initialize(); // initialization task
    Bool_t      LoadLists(); // loading flow lists from input file

    Bool_t      InitTask(FlowTask* task); // initialize FlowTask
    Bool_t      ProcessTask(FlowTask* task); // process FlowTask according to it setting
    Bool_t      ProcessMixed(FlowTask* task); // prepare FlowTask input for mixed harmonics
    Bool_t      ProcessRefs(FlowTask* task); // process reference flow task
    Bool_t      ProcessDirect(FlowTask* task, Short_t iMultBin = 0); // process PID (pion,kaon,proton) flow task
    Bool_t      ProcessReconstructed(FlowTask* task, Short_t iMultBin = 0); // process  V0s flow
    Bool_t      PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor = 0x0, TH3D* h3Entries = 0x0, TH3D* h3EntriesBG = 0x0, TProfile3D* p3CorFour = 0x0); // prepare
    Bool_t      PrepareSlicesNew(FlowTask* task, TString histName, Bool_t bDoCand = kTRUE); // wrapper for making/preparing per-task slices
    Bool_t      MakeProfileSlices(FlowTask* task, TH1* inputProf, TList* outList); // prepare slices out of inputHist
    Bool_t      MakeSparseSlices(FlowTask* task, THnSparse* inputSparse, TList* outList, const char* outName = "hInvMass"); // prepare slices out of 'inputSparse'

    TH1D*       CalcRefCumTwo(TProfile* hTwoRef, FlowTask* task); // calculate cn{2} out of correlation
    TH1D*       CalcRefCumFour(TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate cn{4} out of correlation
    TH1D*       CalcDifCumTwo(TH1D* hTwoDif, FlowTask* task); // calculate dn{2} out of correlation
    TH1D*       CalcDifCumTwo(TProfile* hTwoDif, FlowTask* task); // calculate dn{2} out of correlation
    TH1D*       CalcDifCumFour(TH1D* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate dn{4} out of correlation
    TH1D*       CalcDifCumFour(TProfile* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate dn{4} out of correlation

    TH1D*       CalcRefFlowTwo(TH1D* hTwoRef, FlowTask* task); // calculate vn{2} out of cn{2}
    TH1D*       CalcRefFlowFour(TH1D* hFourRef, FlowTask* task); // calculate vn{4} out of cn{4}
    TH1D*       CalcDifFlowTwo(TH1D* hTwoDif, TH1D* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate vn'{2} out of dn{2} & vn{2}
    TH1D*       CalcDifFlowFour(TH1D* hFourDif, TH1D* hFourRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate vn'{4} out of dn{4} and vn{4}

    void        PrintFitFunction(const TF1* func);
    Bool_t      SetFuncParameters(TF1* func, Double_t* dVec, const std::vector<Double_t>& vecLow, const std::vector<Double_t>& vecHigh, const std::vector<TString> vecNames = {}); // set func parameters & limits (including fixed paramters)
    Bool_t      SetFuncParameters(TF1* func, const std::vector<Double_t>& vecVal, const std::vector<Double_t>& vecLow, const std::vector<Double_t>& vecHigh, const std::vector<TString> vecNames = {}); // set func parameters & limits (including fixed paramters)
    Bool_t      CheckFitResult(TFitResultPtr result, Bool_t bIgnorePOSDEF = kFALSE);
    TH1*        SubtractInvMassBg(TH1* hInvMass, TH1* hInvMassBg, FlowTask* task);
    Bool_t      FitInvMass(TH1* hist, FlowTask* task, TF1& fitOut, TF1& fitOutSig, TF1& fitOutBg, TList* outList, TH1* histBg = nullptr);
    Bool_t      FitCorrelations(TH1* hist, FlowTask* task, TF1& fitOut, TF1& fitOutSig, TF1& fitOutBg, TF1& fitInSig, TF1& fitInBg, TList* outList,  Bool_t bIsFour = kFALSE);

    TList*      LoadSamples(TList* list, TString sHistName, Int_t iNumSamples); // find all samples histos in list
    TH1*        MergeListProfiles(TList* list); // merge list of TProfiles into single TProfile
    TH1*        Merge(TH1* a, TH1* b); // merge two histogram
    TH1D*       DesampleList(TList* list, TH1D* merged, FlowTask* task, TString name, Bool_t bSkipDesampling = kFALSE); // Desample list of samples for estimating the uncertanity
    Bool_t      PlotDesamplingQA(TList* list, TH1D* hDesampled, FlowTask* task); // produce QA plots for result of desampling procedure
    TH1D*       TestRebin(TH1D* hOrig = 0x0, FlowTask* task = 0x0); // testing desample - manual rebin

    void        TestProjections(); // testing projection of reconstructed particles
    TProfile2D* Project3DProfile(const TProfile3D* prof3dorig = 0x0); // making projection out of TProfile3D
    TProfile2D* DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,bool originalRange, bool useUF, bool useOF) const;
    TH2D*       DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY, bool computeErrors, bool originalRange, bool useUF, bool useOF) const;

    std::vector<Double_t>    fdMultBins; // global multiplicity/centrality binning
    Int_t     fiNumMultBins; // number of multiplicity bins (not size of array)

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
    Bool_t      fSaveInterSteps; // flag for saving intermediate steps (correlations, cumulants) into a main output file
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! output file container
    TFile*      ffDesampleFile; //! output file for results of desampling
    TFile*      ffFitsFile; //! output file for fitting procedure
    TList*      flFlow[kUnknown]; //! TList array for input flow profiles
    TList*      flQACharged; //! TList from input file with Charged QA plots / profiles
    TList*      flQAPID; //! TList from input file with PID (pi,K,p) QA plots / profiles
    TList*      flQAPhi; //! TList from input file with Phi QA plots / profiles
    TList*      flQAV0s; //! TList from input file with K0s QA plots / profiles
    std::vector<FlowTask*> fvTasks; // vector of task for individual species proccesing

    ClassDef(ProcessUniFlow,1);
};

#endif
