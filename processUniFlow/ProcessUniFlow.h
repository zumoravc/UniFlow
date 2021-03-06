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

enum        PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kMyProton, kK0s, kLambda, kPhi, kMyUnknown}; // list of all particle species of interest
enum        Cumulants {kNon = 0, kTwo = 2, kFour = 4, kSix = 6, kEight = 8, kTen = 10, kTwelve = 12, kFourteen = 14}; // Cumulants order

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

    static TString     GetSpeciesString(PartSpecies species); // system species string  (kCharged, kK0s, ...)
    static TString     GetSpeciesName(PartSpecies species); // system species name (Charged, K0s, ...)
    static TString     GetSpeciesLabel(PartSpecies species); // readable species name (h^{#pm}, K^{0}_{S}, ...)
    static Bool_t      IsSpeciesDirect(PartSpecies sp) { return (sp == kCharged || sp == kPion || sp == kKaon || sp == kMyProton); } //
    static Bool_t      IsSpeciesReconstructed(PartSpecies sp) { return (sp == kK0s || sp == kLambda || sp == kPhi); } //

    void        AddTask(FlowTask* task = 0x0); // add task to internal lists of all tasks
    Bool_t        Run(); // running the task (main body of the class)
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
    Bool_t      ProcessReconstructed(FlowTask* task); // process  V0s flow
    Bool_t      ProcessSubtraction(FlowTask* task);
    Bool_t      PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor = 0x0, TH3D* h3Entries = 0x0, TH3D* h3EntriesBG = 0x0, TProfile3D* p3CorFour = 0x0); // prepare
    Bool_t      ProcessFMC(FlowTask* task); // process FMC multiparticle correlations
    Bool_t      PrepareSlicesNew(FlowTask* task, TString histName, Bool_t bDoCand = kTRUE); // wrapper for making/preparing per-task slices
    Bool_t      MakeProfileSlices(FlowTask* task, TH1* inputProf, TList* outList); // prepare slices out of inputHist
    Bool_t      MakeSparseSlices(FlowTask* task, THnSparse* inputSparse, TList* outList, const char* outName = "hInvMass"); // prepare slices out of 'inputSparse'
    Int_t       ReturnThird(const Int_t first, const Int_t second); //needed for <<4'>> (permutations)
    Int_t       ReturnIndex3sub(const Int_t first); //needed for <<4'>> (permutations)

    TH1D*       CalcRefCumTwo(TProfile* hTwoRef, FlowTask* task, Int_t rf1Pos = 0, Int_t rf2Pos = 0); // calculate cn{2} out of correlation
    TH1D*       CalcRefCumFour(TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate cn{4} out of correlation
    TH1D*       CalcRefCumFour3sub(TProfile* hFourRef, TProfile* hTwoRef_sub1, TProfile* hTwoRef_sub2, FlowTask* task, Int_t side); // calculate cn{4} out of correlation
    TH1D*       CalcRefCumSix(TProfile* hSixRef, TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task); // calculate cn{6} out of correlation
    TH1D*       CalcRefCumEight(TProfile* hEightRef, TProfile* hSixRef, TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task); // calculate cn{8} out of correlation
    TH1D*       CalcRefCumTen(TProfile* hTenRef, TProfile* hEightRef, TProfile* hSixRef, TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task); // calculate cn{10} out of correlation
    TH1D*       CalcRefCumTwelve(TProfile* hTwelveRef, TProfile* hTenRef, TProfile* hEightRef, TProfile* hSixRef, TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task); // calculate cn{12} out of correlation
    TH1D*       CalcRefCumFourteen(TProfile* hFourteenRef, TProfile* hTwelveRef, TProfile* hTenRef, TProfile* hEightRef, TProfile* hSixRef, TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task); // calculate cn{14} out of correlation

    TH1D*       CalcDifCumTwo(TH1D* hTwoDif, FlowTask* task); // calculate dn{2} out of correlation
    TH1D*       CalcDifCumTwo(TProfile* hTwoDif, FlowTask* task); // calculate dn{2} out of correlation
    TH1D*       CalcDifCumFour(TH1D* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate dn{4} out of correlation
    TH1D*       CalcDifCumFour(TProfile* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate dn{4} out of correlation
    TH1D*       CalcDifCumFour3sub(TH1D* hFourDif, TH1* hTwoDif1, TH1* hTwoDif2, TH1* hTwoRef1, TH1* hTwoRef2, Int_t iRefBin, FlowTask* task); // calculate dn{4} out of correlation
    TH1D*       CalcDifCumFour3sub(TProfile* hFourDif, TH1* hTwoDif1, TH1* hTwoDif2, TH1* hTwoRef1, TH1* hTwoRef2, Int_t iRefBin, FlowTask* task); // calculate dn{4} out of correlation

    TH1D*       CalcRefFlowTwo(TH1D* hTwoRef, FlowTask* task, Int_t rf1Pos = 0, Int_t rf2Pos = 0); // calculate vn{2} out of cn{2}
    TH1D*       CalcRefFlowFour(TH1D* hFourRef, FlowTask* task); // calculate vn{4} out of cn{4}
    TH1D*       CalcRefFlowSix(TH1D* hSixRef, FlowTask* task); // calculate vn{6} out of cn{6}
    TH1D*       CalcRefFlowEight(TH1D* hEightRef, FlowTask* task); // calculate vn{8} out of cn{8}
    TH1D*       CalcRefFlowTen(TH1D* hTenRef, FlowTask* task); // calculate vn{10} out of cn{10}
    TH1D*       CalcRefFlowTwelve(TH1D* hTwelveRef, FlowTask* task); // calculate vn{12} out of cn{12}
    TH1D*       CalcRefFlowFourteen(TH1D* hFourteenRef, FlowTask* task); // calculate vn{14} out of cn{12}

    TH1D*       CalcDifFlowTwo(TH1D* hTwoDif, TH1D* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate vn'{2} out of dn{2} & vn{2}
    TH1D*       CalcDifFlowFour(TH1D* hFourDif, TH1D* hFourRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate vn'{4} out of dn{4} and vn{4}

    //FMCs
    TH1D*       CalcFourFMC(TProfile* hFour, TProfile* hTwo_1, TProfile* hTwo_2, TProfile* hTwo_1_gap, TProfile* hTwo_2_gap, FlowTask* task); // calculate FMC from 4 particle correlations with two different harmonics
    TH1D*       CalcSixTwoDif(TProfile* hSix, TProfile* hFour_12, TProfile* hFour_13, TProfile* hTwo_1, TProfile* hTwo_3, TProfile* hTwo_3_gap, FlowTask* task); // calculate FMC from 6 particle correlations with two different harmonics
    TH1D*       CalcSixThreeDif(TProfile* hSix, TProfile* hFour_12, TProfile* hFour_13, TProfile* hFour_23, TProfile* hTwo_1, TProfile* hTwo_2, TProfile* hTwo_3, TProfile* hTwo_1_gap, TProfile* hTwo_2_gap, TProfile* hTwo_3_gap, FlowTask* task, TProfile* hThree = nullptr); // calculate FMC from 6 particle correlations with three different harmonics
    TH1D*       CalcEight_ThreeOne(TProfile* hEight, TProfile* hSix_123, TProfile* hSix_124, TProfile* hFour_12, TProfile* hFour_14, TProfile* hTwo_1, TProfile* hTwo_4, TProfile* hTwo_4_gap, FlowTask* task); // calculate FMC from 8 particle correlations with two different harmonics -- 3 same, 1 different
    TH1D*       CalcEight_TwoTwo(TProfile* hEight, TProfile* hSix_123, TProfile* hSix_134, TProfile* hFour_12, TProfile* hFour_13, TProfile* hFour_34, TProfile* hTwo_1, TProfile* hTwo_4, FlowTask* task); // calculate FMC from 8 particle correlations with two different harmonics -- 2 & 2 same

    //fluctuations
    TH1D*       CalcMeanFromV(TH1D* hTwo, TH1D* hFour, FlowTask* task); // calculating mean value from vn{2} and vn{4}
    TH1D*       CalcSigmaFromV(TH1D* hTwo, TH1D* hFour, FlowTask* task); // calculating variance value from vn{2} and vn{4}
    TH1D*       CalcRatioV(TH1D* hTwo, TH1D* hFour, FlowTask* task); // calculating variance value from vn{2} and vn{4}
    TH1D*       CalcFlowFluctations(TH1D* hMean, TH1D* hSigma, FlowTask* task); // calculating flow fluctiatons out of mean and sigma

    //subtracted
    TH1D*       CalcSubtracted(FlowTask* task, Int_t iMultBin, TH2D* base, TH2D* raw, TH1D* refCum, TH1D* profDiffBase, TH1D* profDiffRaw); // calculate non-flow subtraction of pt diff with peripheral collisions

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
    TH1D*       DoJackknife(TList* list, TH1D* merged, FlowTask* task, TString name); // Jackknife procedure
    TH1D*       DoBootstrapping(TList* list, TH1D* merged, FlowTask* task, TString name); // Jackknife procedure
    Bool_t      PlotDesamplingQA(TList* list, TH1D* hDesampled, FlowTask* task); // produce QA plots for result of desampling procedure
    TH1D*       TestRebin(TH1D* hOrig = 0x0, FlowTask* task = 0x0); // testing desample - manual rebin
    TH1D*       GetMeanPOI(TH1D** flow, TString name); //get mean value of geometrical combinations
    TH1D*       GetMean(TH1D** flow, TString name, const Int_t max); //get mean value of geometrical combinations

    void        TestProjections(); // testing projection of reconstructed particles
    TProfile2D* Project3DProfile(const TProfile3D* prof3dorig = 0x0); // making projection out of TProfile3D
    TProfile2D* DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,bool originalRange, bool useUF, bool useOF) const;
    TH2D*       DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY, bool computeErrors, bool originalRange, bool useUF, bool useOF) const;

    char sides[4]; //name of sides for 3 sub

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
    TFile*      ffJackFile; //! output file for results of jackknife
    TFile*      ffBSFile; //! output file for results of bootstrap
    TFile*      ffFitsFile; //! output file for fitting procedure
    TList*      flFlow[kMyUnknown]; //! TList array for input flow profiles
    TList*      flQACharged; //! TList from input file with Charged QA plots / profiles
    TList*      flQAPID; //! TList from input file with PID (pi,K,p) QA plots / profiles
    TList*      flQAPhi; //! TList from input file with Phi QA plots / profiles
    TList*      flQAV0s; //! TList from input file with K0s QA plots / profiles
    std::vector<FlowTask*> fvTasks; // vector of task for individual species proccesing

    ClassDef(ProcessUniFlow,1);
};

#endif
