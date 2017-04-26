/* ProcessUniFlow class
 *
 * Class implemented for processing results of AliAnalysisTaskUniFlow task.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"



class FlowTask
{
  public:
    enum    PartSpecies {kUnknown=0, kRefs, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi}; // list of all particle species of interest

                FlowTask(); // default constructor
                FlowTask(const char* name, PartSpecies species = kUnknown); // named constructor
                ~FlowTask(); // default destructor

    void        SetHarmonics(Int_t harm) { fHarmonics = harm; }
    void        SetEtaGap(Float_t eta) { fEtaGap = eta; }
    void        SetNumSamples(Short_t num) { fNumSamples = num; }
    void        SetPtBins(Double_t* array, const Short_t size); // setup the pt binning for this task, where size is number of elements in array

  protected:
  private:
    void        PrintTask(); // listing values of internal properties

    TString     fName; // task name
    PartSpecies fSpecies; // species involved
    Int_t       fHarmonics; // harmonics
    Double_t    fEtaGap; // eta gap
    Short_t     fNumSamples; // [10] number of samples
    static const Short_t fNumPtBinsMax = 100; // initialization (maximum) number of pt bins
    Double_t    fPtBinsEdges[fNumPtBinsMax]; // pt binning
    Short_t     fNumPtBins; // actual number of pT bins (not size of array) for rebinning
};

FlowTask::FlowTask() { fHarmonics = 0; fEtaGap = 0; fNumSamples = 10; fNumPtBins = -1; }
FlowTask::FlowTask(const char* name, PartSpecies species) : FlowTask() { fName = name; fSpecies = species; }
FlowTask::~FlowTask() {}
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

void FlowTask::PrintTask()
{
  printf("----- Printing task info ------\n");
  printf("   fName: %s\n",fName.Data());
  printf("   fSpecies: %d\n",fSpecies);
  printf("   fHarmonics: %d\n",fHarmonics);
  printf("   fEtaGap: %g\n",fEtaGap);
  printf("   fNumPtBins: %d (limit %d)\n",fNumPtBins,fNumPtBinsMax);
  if(fNumPtBins > -1) { printf("   fPtBinsEdges: "); for(Short_t i(0); i < fNumPtBins+1; i++) printf("%g ",fPtBinsEdges[i]); printf("\n"); }
  printf("------------------------------\n");
  return;
}



class ProcessUniFlow
{
  public:

                ProcessUniFlow(); // default constructor
                ~ProcessUniFlow(); // default destructor

    void        SetInputFilePath(const char* path) { fsInputFilePath = path; }
    void        SetInputFileName(const char* name) { fsInputFileName = name; }
    void        SetOutputFilePath(const char* path) { fsOutputFilePath = path; }
    void        SetOutputFileName(const char* name) { fsOutputFileName = name; }
    void        SetTaskName(const char* name) { fsTaskName = name; }
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }
    void        AddTask(FlowTask* task = 0x0); // add task to internal lists of all tasks
    void        Run(); // running the task (main body of the class)
    void        SetMultiplicityBins(Double_t* array, const Short_t size); // setup the global multiplicity binning, where size is number of elements in array

  protected:

  private:
    Bool_t      Initialize(); // initialization task
    Bool_t      LoadLists(); // loading flow lists from input file
    void        TestProjections(); // testing projection of reconstructed particles

    TProfile2D* Project3DProfile(const TProfile3D* prof3dorig = 0x0); // making projection out of TProfile3D
    TProfile2D* DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,bool originalRange, bool useUF, bool useOF) const;
    TH2D*       DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY, bool computeErrors, bool originalRange, bool useUF, bool useOF) const;

    void        ProcessTask(FlowTask* task = 0x0); // process FlowTask according to it setting
    Bool_t      ProcessRefs(FlowTask* task = 0x0); // process reference flow task
    Bool_t      ProcessV0s(FlowTask* task = 0x0); // process  V0s flow
    Bool_t 	    ExtractFlowPhi(TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowK0s(TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for K0s candidates
		Bool_t 	    ExtractFlowLambda(TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for Lambda candidates

    // printing output methods
    void        Fatal(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Error(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Warning(TString sMsg, TString sMethod = ""); // printf the msg as warning
    void        Info(TString sMsg, TString sMethod = ""); // printf the msg as info
    void        Debug(TString sMsg, TString sMethod = ""); // printf the msg as info

    static const Short_t fiNumMultBinsGlobal = 20; // global initialization number of bins
    Double_t    fdMultBins[fiNumMultBinsGlobal]; // global multiplicity/centrality binning
    Short_t     fiNumMultBins; // number of multiplicity bins (not size of array)

    TString     fsInputFilePath; // path to the input folder with input file
    TString     fsInputFileName; // name of input file
    TString     fsOutputFilePath; // path to the ouput folder
    TString     fsOutputFileName; // name of output file
    TString     fsTaskName; // name of task (inchluded in data structure names)
    TString     fsOutputFileFormat; // [pdf] format of output files (pictures)

    Bool_t      fbInit; // flag for initialization status
    Bool_t      fbDebug; // flag for debugging : if kTRUE Debug() messages are displayed
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! input file container
    TList*      flFlowRefs; //! TList from input file with RFPs flow profiles
    TList*      flFlowCharged; //! TList from input file with Charged flow profiles
    TList*      flFlowPID; //! TList from input file with PID (pi,K,p) flow profiles
    TList*      flFlowPhi; //! TList from input file with Phi flow profiles
    TList*      flFlowK0s; //! TList from input file with K0s flow profiles
    TList*      flFlowLambda; //! TList from input file with Lambda flow profiles
    std::vector<FlowTask*> fvTasks; // vector of task for individual species proccesing

};
//_____________________________________________________________________________
ProcessUniFlow::ProcessUniFlow() :
  fbDebug(kFALSE),
  fbInit(kFALSE),
  ffInputFile(0x0),
  ffOutputFile(0x0),
  flFlowRefs(0x0),
  flFlowCharged(0x0),
  flFlowPID(0x0),
  flFlowPhi(0x0),
  flFlowK0s(0x0),
  flFlowLambda(0x0)
{
  // default constructor
  fsInputFilePath = TString("");
  fsInputFileName = TString("AnalysisResults.root");
  fsOutputFilePath = TString("");
  fsOutputFileName = TString("UniFlow.root");
  fsTaskName = TString("UniFlow");
  fsOutputFileFormat = TString("pdf");
  fvTasks = std::vector<FlowTask*>();

  for(Short_t i(0); i < fiNumMultBinsGlobal; i++) fdMultBins[i] = 0;
}
//_____________________________________________________________________________
ProcessUniFlow::~ProcessUniFlow()
{
  // default destructor
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;
}
//_____________________________________________________________________________
void ProcessUniFlow::Run()
{
  // main body of the class
  if(!Initialize()) { Fatal("Task not initialized","Run"); return; }

  Debug("Initialized");

  // flFlow_K0s->ls();
  // TestProjections();

  const Short_t iNumTasks = fvTasks.size();
  FlowTask* currentTask = 0x0;
  Info("===== Running over tasks ======","Run");
  Info(Form("  Number of tasks: %d\n",iNumTasks),"Run");
  for(Short_t iTask(0); iTask < iNumTasks; iTask++)
  {
    currentTask = fvTasks.at(iTask);
    if(!currentTask) continue;
    ProcessTask(currentTask);
  }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Initialize()
{
  // initialization of all necessery prerequisits
  Info("Initializating task","Initialize");
  fbInit = kFALSE;

  // opening input file
  ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(),fsInputFileName.Data()),"READ");
  if(!ffInputFile || !ffInputFile->IsOpen())
  {
    Fatal(Form("Input file %s/%s not open",fsInputFilePath.Data(),fsInputFileName.Data()),"Initialize");
    return fbInit;
  }

  // opening output file
  ffOutputFile = new TFile(Form("%s/%s",fsOutputFilePath.Data(),fsOutputFileName.Data()),"RECREATE");
  if(!ffOutputFile || !ffOutputFile->IsOpen())
  {
    Fatal(Form("Output file %s/%s not open",fsOutputFilePath.Data(),fsOutputFileName.Data()),"Initialize");
    return fbInit;
  }
  Info("Files loaded","Initialize");

  if(!LoadLists()) return fbInit;
  Info("Flow lists loaded","Initialize");

  // initialization succesfull
  fbInit = kTRUE;
  Info("Initialization succesfull","Initialize");
  return fbInit;
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
Bool_t ProcessUniFlow::LoadLists()
{
  // loading TLists into task
  if(!ffInputFile) { Fatal("Input file does not exists!","LoadLists"); return kFALSE; }
  ffInputFile->cd(fsTaskName.Data());

  flFlowRefs = (TList*) gDirectory->Get(Form("Flow_Refs_%s",fsTaskName.Data()));
  if(!flFlowRefs) { Fatal("flFlow_Refs list does not exists!","LoadLists"); return kFALSE; }
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

  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::ProcessTask(FlowTask* task)
{
  Info(Form("Processing task: %s",task->fName.Data()),"ProcessTask");
  if(!task) { Error("Task not valid!","ProcessTask"); return; }

  task->PrintTask();

  Bool_t bProcessed = kFALSE;
  switch (task->fSpecies)
  {
    case FlowTask::kRefs:
      bProcessed = ProcessRefs(task);
      break;

    case FlowTask::kPhi:
    case FlowTask::kK0s:
    case FlowTask::kLambda:
      bProcessed = ProcessV0s(task);
      break;
    default: break;
  }

  if(!bProcessed) { Error("Task not processed correctly!","ProcessTask"); return; }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessRefs(FlowTask* task)
{
  Info("Processing Refs task","ProcesRefs");
  if(!task) { Error("Task not valid!","ProcessRefs"); return kFALSE; }
  if(task->fSpecies != FlowTask::kRefs) { Error("Task species not kRefs!","ProcessRefs"); return kFALSE; }

  // Info("Listing input refs profiles","ProcesRefs");
  // flFlowRefs->ls();

  // merging samples together
  TProfile* prof = 0x0;
  TList* list = new TList();
  // for(Short_t i(0); i < 2; i++)
  for(Short_t i(0); i < task->fNumSamples; i++)
  {
    prof = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_<2>_harm%d_gap%02.2g_sample%d",task->fHarmonics,10*task->fEtaGap,i));
    if(!prof) { Warning(Form("Profile sample %d does not exits. Skipping",i),"ProcesRefs"); continue; }
    list->Add(prof);
  }
  Debug(Form("Number of samples in list pre merging %d",list->GetEntries()));


  TProfile* merged = (TProfile*) prof->Clone();
  Double_t mergeStatus = merged->Merge(list);
  if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessRefs"); return kFALSE; }
  // merged->Draw();

  // merged->SetName(Form("fpRefs_<2>_harm%d_gap%g",task->fHarmonics,task->fEtaGap));
  // merged->SetTitle(Form("Ref: <<2>> | n=%d | Gap %02.2g",task->fHarmonics,task->fEtaGap));

  // rebinning: multiplicity
  // Double_t xbins[] = {1,14,16,60};
  // Int_t iNumBins = sizeof(xbins)/sizeof(xbins[0]) - 1;
  // printf("bins: %d\n",iNumBins);

  TProfile* rebin = (TProfile*) merged->Rebin(fiNumMultBins,Form("%s_rebin",merged->GetName()),fdMultBins);
  TH1D* histRebin = rebin->ProjectionX();

  // at this point correlations are processed
  // start doing vns out of them

  TH1D* hFlow = (TH1D*) histRebin->Clone(Form("hFlow_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
  hFlow->SetTitle(Form("Ref: v_{%d}{2 | Gap %g}",task->fHarmonics,task->fEtaGap));
  // TH1D* hFlow = (TH1D*) histRebin->Clone("hFlow");
  // hFlow->Reset();

  const Short_t iBinsX = hFlow->GetNbinsX();

  Double_t dContent = 0, dValue = 0;
  for(Short_t iBin(1); iBin < iBinsX+1; iBin++)
  {
    dContent = histRebin->GetBinContent(iBin);
    if(dContent < 0) hFlow->SetBinContent(iBin, -9.);
    else hFlow->SetBinContent(iBin,TMath::Sqrt(dContent));
    // printf("%g | %g\n",dContent, TMath::Sqrt(dContent));
  }

  // printf("%g±%g\n", merged->GetBinContent(14), merged->GetBinError(14));
  // printf("%g±%g\n", merged->GetBinContent(15), merged->GetBinError(15));
  // printf("%g±%g\n", histRebin->GetBinContent(2), histRebin->GetBinError(2));



  // TCanvas* canTest = new TCanvas("canTestRefs","canTestReffs",1000,1000);
  // canTest->Divide(2,2);
  // canTest->cd(1);
  // merged->Draw();
  // canTest->cd(2);
  // rebin->Draw();
  // canTest->cd(3);
  // histRebin->Draw();
  // canTest->cd(4);
  // hFlow->Draw();


  ffOutputFile->cd();
  hFlow->Write();


  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessV0s(FlowTask* task)
{
  Info("Processing V0s task","ProcessV0s");
  if(!task) { Error("Task not valid!","ProcessV0s"); return kFALSE; }
  if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessV0s"); return kFALSE; }

  // preparing particle dependent variables for switch
  //  -- input histos / profiles with entries and correlations
  TH3D* histEntries = 0x0;
  TH3D* histBG = 0x0; // etries for BG (phi)
  TProfile3D* profFlow = 0x0;
  //  -- naming variables
  TString sSpeciesName; // in objects name
  TString sSpeciesLabel; // LaTeX for titles

  // checking particles species and assigning particle dependent variables
  switch (task->fSpecies)
  {
    case FlowTask::kPhi :
      histEntries = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_gap%02.2g",10*task->fEtaGap));
      histBG = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g",10*task->fEtaGap));
      if(!histBG) { Error("Histo with BG entries not found","ProcessV0s"); return kFALSE; }
      profFlow = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_<2>_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
      sSpeciesName = TString("Phi");
      sSpeciesLabel = TString("#phi");
    break;

    case FlowTask::kK0s :
      histEntries = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_gap%02.2g",10*task->fEtaGap));
      profFlow = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_<2>_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
      sSpeciesName = TString("K0s");
      sSpeciesLabel = TString("K^{0}_{S}");
    break;

    case FlowTask::kLambda :
      histEntries = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_gap%02.2g",10*task->fEtaGap));
      profFlow = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_<2>_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
      sSpeciesName = TString("Lambda");
      sSpeciesLabel = TString("#Lambda/#bar{#Lambda}");
    break;

    default:
      Error("Task species not V0s nor Phi!","ProcessV0s");
      return kFALSE;
  }

  if(!histEntries) { Error("Entries histos not found!","ProcessV0s"); return kFALSE; }
  if(!profFlow) { Error("Cumulant histos not found!","ProcessV0s"); return kFALSE; }

  TH1D* hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
  if(!hRefFlow)
  {
    Warning("Relevant Reference flow not found within output file.","ProcessV0s");
    Info("Creating relevant reference flow task.","ProcessV0s");

    FlowTask* taskRef = new FlowTask("Ref",FlowTask::kRefs);
    taskRef->SetHarmonics(task->fHarmonics);
    taskRef->SetEtaGap(task->fEtaGap);
    if(ProcessRefs(taskRef)) return ProcessV0s(task);
    else { Error("Something went wrong when running automatic refs flow task:","ProcessV0s"); taskRef->PrintTask(); return kFALSE; }
  }

  // preparing for loops
  TH1D* hFlow = 0x0;
  Double_t dRefs = 0, dRefsError = 0;
  Short_t binMultLow = 0, binMultHigh = 0;
  Short_t binPtLow = 0, binPtHigh = 0;
  TH1D* hInvMass = 0x0;
  TH1D* hInvMassBG = 0x0;
  TProfile2D* prof2DFlowMass = 0x0;
  TProfile* profFlowMass = 0x0;
  TH1D* hFlowMass = 0x0;
  Double_t dContent = 0, dError = 0;
  Double_t dFlow = 0, dFlowError = 0; // containers for flow extraction results
  TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200); // canvas for fitting results
  TCanvas* cFlow = new TCanvas("cFlow","Flow",400,400); // canvas for resulting flow vs pt

  for(Short_t binMult(0); binMult < fiNumMultBins; binMult++)
  {
    // preparing resulting flow vs pt (mult) histo
    hFlow = new TH1D(Form("hFlow_%s_harm%d_gap%02.2g_mult%d",sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,binMult),Form("%s: v_{%d}{2 | Gap %g} (%g - %g); #it{p}_{T} (GeV/#it{c})",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[binMult],fdMultBins[binMult+1]), task->fNumPtBins,task->fPtBinsEdges);

    // estimating multiplicity edges and selection axis range
    binMultLow = histEntries->GetXaxis()->FindFixBin(fdMultBins[binMult]);
    binMultHigh = histEntries->GetXaxis()->FindFixBin(fdMultBins[binMult+1]) - 1; // for rebin both bins are included (so that one needs to lower)
    profFlow->GetXaxis()->SetRange(binMultLow,binMultHigh); // setting multiplicity range

    // loading relevant reference flow
    dRefs = hRefFlow->GetBinContent(binMult+1);
    dRefsError = hRefFlow->GetBinError(binMult+1);


    for(Short_t binPt(0); binPt < task->fNumPtBins; binPt++)
    {
      // estimating pt edges
      binPtLow = histEntries->GetYaxis()->FindFixBin(task->fPtBinsEdges[binPt]);
      binPtHigh = histEntries->GetYaxis()->FindFixBin(task->fPtBinsEdges[binPt+1]) - 1; // for rebin both bins are included (so that one needs to lower)

      // rebinning entries based on mult & pt binning
      hInvMass = (TH1D*) histEntries->ProjectionZ("hInvMass",binMultLow,binMultHigh,binPtLow,binPtHigh,"e");
      if(histBG) hInvMassBG = (TH1D*) histBG->ProjectionZ("hInvMassBG",binMultLow,binMultHigh,binPtLow,binPtHigh,"e");
      // here can be rebinning of inv mass if needed
      // hInvMass ready to fitting

      // projection of flow-mass profile
      prof2DFlowMass = Project3DProfile(profFlow); // doing projection with modified function
      profFlowMass = (TProfile*) prof2DFlowMass->ProfileX("profFlowMass",binPtLow,binPtHigh);
      hFlowMass = (TH1D*) profFlowMass->ProjectionX("hFlowMass");
      // NOTE: this is the ONLY (for some freaking reason) way how to get proper TH1 wth <<2>> out of TProfile3D

      // make flow out of corralation
      for(Short_t iMass(1); iMass < profFlowMass->GetNbinsX()+1; iMass++)
      {
        dContent = profFlowMass->GetBinContent(iMass);
        dError = profFlowMass->GetBinError(iMass);
        hFlowMass->SetBinContent(iMass, dContent / dRefs);
        hFlowMass->SetBinError(iMass, TMath::Sqrt(TMath::Power(dError/dRefs,2) + TMath::Power(dContent*dRefsError/(dRefs*dRefs),2)) );
        // printf("bin %d | %g±%g\n",iMass, dContent,dError);
      }
      // hFlowMass is ready for fitting

      // extracting flow
      switch (task->fSpecies)
      {
        case FlowTask::kPhi :
          if( !ExtractFlowPhi(hInvMass,hInvMassBG,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull","ProcessV0s"); return kFALSE; }
        break;

        case FlowTask::kK0s :
          if( !ExtractFlowK0s(hInvMass,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull","ProcessV0s"); return kFALSE; }
        break;

        case FlowTask::kLambda :
          if( !ExtractFlowLambda(hInvMass,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull","ProcessV0s"); return kFALSE; }
        break;


        default : break;
      }
      Info(Form("Extracted flow: %g ± %g\n",dFlow,dFlowError),"ProcessV0s");

      // saving pt bin resulting flow to flow-vs-pt hist
      hFlow->SetBinContent(binPt+1, dFlow);
      hFlow->SetBinError(binPt+1, dFlowError);

      // printing fitting results
      canFitInvMass->SaveAs(Form("%s/FlowMass_%s_n%d2_gap%02.2g_cent%d_pt%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,binMult,binPt,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

    } // endfor {binPt}

    cFlow->cd();
    hFlow->Draw();
    cFlow->SaveAs(Form("%s/Flow_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,binMult,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

    ffOutputFile->cd();
    hFlow->Write();

  } // endfor {binMult}

  // TCanvas* canTest = new TCanvas("canTest","Test",600,600);
  // canTest->Divide(2,3);
  // canTest->cd(1);
  // histEntries->Draw();
  // canTest->cd(2);
  // profFlow->Draw();
  // canTest->cd(3);
  // hInvMass->Draw();
  // canTest->cd(4);
  // profFlowMass->Draw();
  // canTest->cd(5);
  // hFlowMass->Draw();
  // canTest->cd(6);
  // prof2DFlowMass->Draw("colz");


  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::AddTask(FlowTask* task)
{
  if(!task) return;

  fvTasks.push_back(task);
  return;
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
Bool_t ProcessUniFlow::ExtractFlowPhi(TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
  if(!hInvMass) { Error("InvMass (signal) histogram does not exists!","ExtractFlowPhi"); return kFALSE; }
  if(!hInvMassBG) { Error("InvMass (BG) histogram does not exists!","ExtractFlowPhi"); return kFALSE; }
	if(!hFlowMass) { Error("Flow-mass histogram does not exists!","ExtractFlowPhi"); return kFALSE;	}
	if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowPhi"); return kFALSE; }

	// Reseting the canvas (removing drawn things)
	canFitInvMass->Clear();

	TString sOutputFormat = fsOutputFileFormat;
	const Short_t iNumSigmas = 7;
	Double_t dMeanShot = 0;
	Double_t dSigmaShot = 0;
	Double_t dMassLow = 0;
	Double_t dMassHigh = 0;





  // subtraction of BG
  // normalisation BG to Signal in region of interest
  const Double_t dNormMassLow = 1.03;
  const Double_t dNormMassHigh = 1.06;

  Short_t iBinNormLow = hInvMassBG->FindFixBin(dNormMassLow);
  Short_t iBinNormHigh = hInvMassBG->FindFixBin(dNormMassHigh);
  printf("Bin low %g | high %g\n",hInvMassBG->GetBinCenter(iBinNormLow),hInvMassBG->GetBinCenter(iBinNormHigh));

  Double_t dIntSignal = 0;
  Double_t dIntBG = 0;
  for(Short_t iBin(iBinNormLow); iBin < iBinNormHigh+1; iBin++)
  {
    dIntSignal += hInvMass->GetBinContent(iBin);
    dIntBG += hInvMassBG->GetBinContent(iBin);
  }
  Double_t dScaleFact = dIntSignal/dIntBG;
  printf("integral: sig %g | bg %g | scale %g\n",dIntSignal,dIntBG,dScaleFact);

  // TODO which one?
  TH1D* hInvMassBG_scaled = (TH1D*) hInvMassBG->Clone("hInvMassBG_scaled");
  // hInvMass->Scale(1/dIntSignal);
  // hInvMassBG_scaled->Scale(dScaleFact);
  // hInvMassBG_scaled->Scale(1/dScaleFact);
  hInvMassBG_scaled->SetLineColor(kRed);

  TH1D* hInvMass_subs = (TH1D*) hInvMass->Clone("hInvMass_subs");
  hInvMass_subs->Add(hInvMassBG_scaled,-1);


  // TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200);
	canFitInvMass->Divide(3,3);
  canFitInvMass->cd(7);
  hInvMass->Draw();
  hInvMassBG->Draw("same");
  hInvMassBG_scaled->Draw("same");
  canFitInvMass->cd(8);
  hInvMass_subs->Draw();

  return kFALSE;

	TH1D* hInvMass_side = 0x0;
	TH1D* hInvMass_residual = 0x0;
	TH1D* hInvMass_ratio = 0x0;

	TH1D* hFlowMass_side = 0x0;

	const Short_t iNumBinsMassFlow = hFlowMass->GetNbinsX();

	// inv mass fitts
	TF1* fitShot = 0x0;
	TF1* fitSide = 0x0;
	TF1* fitRatio = 0x0;

	// flow mass fits
	TF1* fitFlowSide = 0x0;
	TF1* fitFlowTot = 0x0;


	hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side");
	hInvMass_residual = (TH1D*) hInvMass->Clone("hInvMass_residual");


	printf("\n====== K0s: Fitting InvMass first shot ========\n");
	canFitInvMass->cd(1);
	fitShot = new TF1("fitShot","gaus(0)+pol2(3)",0.4,0.6);
	fitShot->SetNpx(10000);
	fitShot->SetParameter(1,0.5);
	fitShot->SetParLimits(1,0.485,0.515);
	fitShot->SetParameter(2,0.01);
	fitShot->SetParLimits(2,0.,0.05);
	hInvMass->Fit("fitShot","R");

	// TODO checking the fitting results

	// extract mean & sigma for sidebands fitting reagion
	dMeanShot = fitShot->GetParameter(1);
	dSigmaShot = fitShot->GetParameter(2);
	dMassLow = dMeanShot - iNumSigmas*dSigmaShot;
	dMassHigh = dMeanShot + iNumSigmas*dSigmaShot;
	printf("=========================\nFitting region: %f - %f \n==========================\n", dMassLow,dMassHigh);

	const Short_t iNumBinsMass = hInvMass->GetNbinsX();
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hInvMass_side->GetBinCenter(iMass) > dMassLow && hInvMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hInvMass_side->SetBinError(iMass,9999999999999);
		}
	}
	canFitInvMass->cd(2);
	hInvMass_side->SetMaximum(hInvMass->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_side->SetMinimum(0);
	hInvMass_side->Draw();

	// fitting background in sidebands
	printf("\n====== K0s: Fitting InvMass side bands ========\n");
	fitSide = new TF1("fitSide","pol2(0)",0.4,0.6);
	fitSide->SetNpx(10000);
	hInvMass_side->Fit("fitSide","R");

	Double_t dContent = 0;
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		dContent = hInvMass_side->GetBinContent(iMass) - fitSide->Eval(hInvMass_side->GetBinCenter(iMass));
		hInvMass_residual->SetBinContent(iMass,dContent);
	}

	canFitInvMass->cd(3);
	hInvMass_residual->Draw();

	hInvMass_ratio = (TH1D*) hInvMass_residual->Clone("hInvMass_ratio");
	hInvMass_ratio->Sumw2();
	hInvMass_ratio->Divide(hInvMass);

	printf("\n====== K0s: Fitting InvMass sig/tot ratio ========\n");
	canFitInvMass->cd(4);

	//hInvMass_ratio->Draw();
	fitRatio = new TF1("fitRatio","gaus(0)+pol3(3)",0.4,0.6);
	fitRatio->SetNpx(1000);
	fitRatio->SetParameter(0,0.98);
	fitRatio->SetParameter(1,0.5);
	fitRatio->SetParLimits(1,0.48,0.51);
	fitRatio->SetParameter(2,0.01);
	fitRatio->SetParLimits(2,0.,0.05);
	hInvMass_ratio->Fit("fitRatio","R");

	// flow mass Fitting
	hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");
	hFlowMass_side->SetMaximum(1.5*hFlowMass->GetMaximum());
	hFlowMass_side->SetMinimum(0.5*hFlowMass->GetMinimum());

	// fitting side bands
	for(Short_t iMass(1); iMass < iNumBinsMassFlow+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hFlowMass_side->GetBinCenter(iMass) > dMassLow && hFlowMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hFlowMass_side->SetBinError(iMass,9999999999999);
		}
	}

	printf("\n====== K0s: Fitting FlowMass sidebands ========\n");
	canFitInvMass->cd(6);
	fitFlowSide = new TF1("fitFlowSide","pol2(0)",0.4,0.6);
	hFlowMass_side->Fit("fitFlowSide","R");

	canFitInvMass->cd(5);
	printf("\n====== K0s: Fitting FlowMass total flow ========\n");
	fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol2(8)",0.4,0.6);
	// Inv mass ratio signal/total
	fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
	fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
	fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
	fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
	fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
	fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
	fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
	// FlowMass backround / sidebands
	fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(0));
	fitFlowTot->FixParameter(9,fitFlowSide->GetParameter(1));
	fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
	hFlowMass->Fit("fitFlowTot","R");

	dFlow = fitFlowTot->GetParameter(0);
	dFlowError = fitFlowTot->GetParError(0);

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowK0s(TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
	if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowK0s"); return kFALSE; }
	if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowK0s"); return kFALSE; }
	if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowK0s"); return kFALSE; }

	// Reseting the canvas (removing drawn things)
	canFitInvMass->Clear();

	// Fitting K0s
	const TString sOutputFormat = fsOutputFileFormat;
	const Short_t iNumSigmas = 7;
	Double_t dMeanShot = 0;
	Double_t dSigmaShot = 0;
	Double_t dMassLow = 0;
	Double_t dMassHigh = 0;

	//TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200);
	canFitInvMass->Divide(3,2);

	TH1D* hInvMass_side = 0x0;
	TH1D* hInvMass_residual = 0x0;
	TH1D* hInvMass_ratio = 0x0;

	TH1D* hFlowMass_side = 0x0;

	const Short_t iNumBinsMassFlow = hFlowMass->GetNbinsX();

	// inv mass fitts
	TF1* fitShot = 0x0;
	TF1* fitSide = 0x0;
	TF1* fitRatio = 0x0;

	// flow mass fits
	TF1* fitFlowSide = 0x0;
	TF1* fitFlowTot = 0x0;


	hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side");
	hInvMass_residual = (TH1D*) hInvMass->Clone("hInvMass_residual");

	printf("\n====== K0s: Fitting InvMass first shot ========\n");
	canFitInvMass->cd(1);
	fitShot = new TF1("fitShot","gaus(0)+pol2(3)",0.4,0.6);
	fitShot->SetNpx(10000);
	fitShot->SetParameter(1,0.5);
	fitShot->SetParLimits(1,0.485,0.515);
	fitShot->SetParameter(2,0.01);
	fitShot->SetParLimits(2,0.,0.05);
	hInvMass->Fit("fitShot","R");

	// TODO checking the fitting results

	// extract mean & sigma for sidebands fitting reagion
	dMeanShot = fitShot->GetParameter(1);
	dSigmaShot = fitShot->GetParameter(2);
	dMassLow = dMeanShot - iNumSigmas*dSigmaShot;
	dMassHigh = dMeanShot + iNumSigmas*dSigmaShot;
	printf("=========================\nFitting region: %f - %f \n==========================\n", dMassLow,dMassHigh);

	const Short_t iNumBinsMass = hInvMass->GetNbinsX();
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hInvMass_side->GetBinCenter(iMass) > dMassLow && hInvMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hInvMass_side->SetBinError(iMass,9999999999999);
		}
	}
	canFitInvMass->cd(2);
	hInvMass_side->SetMaximum(hInvMass->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_side->SetMinimum(0);
	hInvMass_side->Draw();

	// fitting background in sidebands
	printf("\n====== K0s: Fitting InvMass side bands ========\n");
	fitSide = new TF1("fitSide","pol2(0)",0.4,0.6);
	fitSide->SetNpx(10000);
	hInvMass_side->Fit("fitSide","R");

	Double_t dContent = 0;
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		dContent = hInvMass_side->GetBinContent(iMass) - fitSide->Eval(hInvMass_side->GetBinCenter(iMass));
		hInvMass_residual->SetBinContent(iMass,dContent);
	}

	canFitInvMass->cd(3);
	hInvMass_residual->Draw();

	hInvMass_ratio = (TH1D*) hInvMass_residual->Clone("hInvMass_ratio");
	hInvMass_ratio->Sumw2();
	hInvMass_ratio->Divide(hInvMass);

	printf("\n====== K0s: Fitting InvMass sig/tot ratio ========\n");
	canFitInvMass->cd(4);

	//hInvMass_ratio->Draw();
	fitRatio = new TF1("fitRatio","gaus(0)+pol3(3)",0.4,0.6);
	fitRatio->SetNpx(1000);
	fitRatio->SetParameter(0,0.98);
	fitRatio->SetParameter(1,0.5);
	fitRatio->SetParLimits(1,0.48,0.51);
	fitRatio->SetParameter(2,0.01);
	fitRatio->SetParLimits(2,0.,0.05);
	hInvMass_ratio->Fit("fitRatio","R");

	// flow mass Fitting
	hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");
	hFlowMass_side->SetMaximum(1.5*hFlowMass->GetMaximum());
	hFlowMass_side->SetMinimum(0.5*hFlowMass->GetMinimum());

	// fitting side bands
	for(Short_t iMass(1); iMass < iNumBinsMassFlow+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hFlowMass_side->GetBinCenter(iMass) > dMassLow && hFlowMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hFlowMass_side->SetBinError(iMass,9999999999999);
		}
	}

	printf("\n====== K0s: Fitting FlowMass sidebands ========\n");
	canFitInvMass->cd(6);
	fitFlowSide = new TF1("fitFlowSide","pol2(0)",0.4,0.6);
	hFlowMass_side->Fit("fitFlowSide","R");

	canFitInvMass->cd(5);
	printf("\n====== K0s: Fitting FlowMass total flow ========\n");
	fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol2(8)",0.4,0.6);
	// Inv mass ratio signal/total
	fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
	fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
	fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
	fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
	fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
	fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
	fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
	// FlowMass backround / sidebands
	fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(0));
	fitFlowTot->FixParameter(9,fitFlowSide->GetParameter(1));
	fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
	hFlowMass->Fit("fitFlowTot","R");

	dFlow = fitFlowTot->GetParameter(0);
	dFlowError = fitFlowTot->GetParError(0);

	return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowLambda(TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
	if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowLambda"); return kFALSE; }
	if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowLambda");	return kFALSE; }
	if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowLambda"); return kFALSE; }

	// Reseting the canvas (removing drawn things)
	canFitInvMass->Clear();

	// Fitting K0s
	const TString sOutputFormat = fsOutputFileFormat;
	const Short_t iNumSigmas = 6;
	const Double_t fitLimitLow = 1.095;
	const Double_t fitLimitHigh = 1.15;

	Double_t dMeanShot = 0;
	Double_t dSigmaShot = 0;
	Double_t dMassLow = 0;
	Double_t dMassHigh = 0;

	//TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200);
	canFitInvMass->Divide(3,2);

	TH1D* hInvMass_side = 0x0;
	TH1D* hInvMass_residual = 0x0;
	TH1D* hInvMass_ratio = 0x0;

	TH1D* hFlowMass_side = 0x0;

	const Short_t iNumBinsMassFlow = hFlowMass->GetNbinsX();

	// inv mass fitts
	TF1* fitShot = 0x0;
	TF1* fitSide = 0x0;
	TF1* fitRatio = 0x0;

	// flow mass fits
	TF1* fitFlowSide = 0x0;
	TF1* fitFlowTot = 0x0;


	hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side");
	hInvMass_residual = (TH1D*) hInvMass->Clone("hInvMass_residual");

	canFitInvMass->cd(1);
	printf("\n====== Lambda: Fitting InvMass first shot ========\n");
	fitShot = new TF1("fitShot","gaus(0)+pol2(3)",1.1,1.13);
	fitShot->SetNpx(10000);
	fitShot->SetParameter(1,1.115);
	fitShot->SetParLimits(1,1.113,1.12);
	fitShot->SetParameter(2,0.01);
  fitShot->SetParLimits(2,0.,0.002);
	hInvMass->Fit("fitShot","RI");

	// TODO checking the fitting results

	// extract mean & sigma for sidebands fitting reagion
	dMeanShot = fitShot->GetParameter(1);
	dSigmaShot = fitShot->GetParameter(2);
	dMassLow = dMeanShot - iNumSigmas*dSigmaShot;
	dMassHigh = dMeanShot + iNumSigmas*dSigmaShot;
	printf("=========================\nFitting region: %f - %f \n==========================\n", dMassLow,dMassHigh);

	// return kTRUE; // testing

	const Short_t iNumBinsMass = hInvMass->GetNbinsX();
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hInvMass_side->GetBinCenter(iMass) > dMassLow && hInvMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hInvMass_side->SetBinError(iMass,9999999999999);
		}
	}
	canFitInvMass->cd(2);
	hInvMass_side->SetMaximum(hInvMass->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_side->SetMinimum(0);
	hInvMass_side->Draw();

	// fitting background in sidebands
	printf("\n====== Lambda: Fitting InvMass sidebands ========\n");
	fitSide = new TF1("fitSide","pol3(0)",fitLimitLow,fitLimitHigh);
	fitSide->SetNpx(10000);
	hInvMass_side->Fit("fitSide","RI");

	Double_t dContent = 0;
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		dContent = hInvMass_side->GetBinContent(iMass) - fitSide->Eval(hInvMass_side->GetBinCenter(iMass));
		hInvMass_residual->SetBinContent(iMass,dContent);
	}

	canFitInvMass->cd(3);
	hInvMass_residual->Draw();

	hInvMass_ratio = (TH1D*) hInvMass_residual->Clone("hInvMass_ratio");
	hInvMass_ratio->Sumw2();
	hInvMass_ratio->Divide(hInvMass);
	hInvMass_ratio->SetMinimum(-0.5);
	hInvMass_ratio->SetMaximum(1.);

	canFitInvMass->cd(4);
	printf("\n====== Lambda: Fitting InvMass sig/total ratio ========\n");
	//hInvMass_ratio->Draw();
	fitRatio = new TF1("fitRatio","gaus(0)+pol3(3)",fitLimitLow,fitLimitHigh);
	fitRatio->SetNpx(1000);
	fitRatio->SetParameter(0,0.98);
	fitRatio->SetParameter(1,1.115);
	fitRatio->SetParameter(2,0.001);
	hInvMass_ratio->Fit("fitRatio","RI");

	// flow mass Fitting
	hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");
	hFlowMass_side->SetMaximum(1.5*hFlowMass->GetMaximum());
	hFlowMass_side->SetMinimum(0.5*hFlowMass->GetMinimum());

	// fitting side bands
	for(Short_t iMass(1); iMass < iNumBinsMassFlow+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hFlowMass_side->GetBinCenter(iMass) > dMassLow && hFlowMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hFlowMass_side->SetBinError(iMass,9999999999999);
		}
	}

	canFitInvMass->cd(6);
	printf("\n====== Lambda: Fitting FlowMass sidebands ========\n");
	fitFlowSide = new TF1("fitFlowSide","pol3(0)",fitLimitLow,fitLimitHigh);
	hFlowMass_side->Fit("fitFlowSide","RI");

	canFitInvMass->cd(5);
	printf("\n====== Lambda: Fitting FlowMass total flow ========\n");
	fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol3(8)",fitLimitLow,fitLimitHigh);

	// Inv mass ratio signal/total
	fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
	fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
	fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
	fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
	fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
	fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
	fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
	// FlowMass backround / sidebands
	fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(0));
	fitFlowTot->FixParameter(9,fitFlowSide->GetParameter(1));
	fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
	fitFlowTot->FixParameter(11,fitFlowSide->GetParameter(3));
	hFlowMass->Fit("fitFlowTot","RIB");

	dFlow = fitFlowTot->GetParameter(0);
	dFlowError = fitFlowTot->GetParError(0);

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
