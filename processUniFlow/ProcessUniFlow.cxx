#include "ProcessUniFlow.h"
#include "FlowTask.h"

#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
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
#include "THnSparse.h"

ClassImp(ProcessUniFlow);

//_____________________________________________________________________________
ProcessUniFlow::ProcessUniFlow() :
  fbDebug{kFALSE},
  fbInit{kFALSE},
  fbSaveMult{kFALSE},
  fFlowFitCumulants{kFALSE},
  fSaveInterSteps{kFALSE},
  sides{"LMR"},
  fdMultBins{},
  fiNumMultBins{0},
  ffInputFile{nullptr},
  ffOutputFile{nullptr},
  ffFitsFile{nullptr},
  ffDesampleFile{nullptr},
  fsInputFilePath{},
  fsInputFileName{"AnalysisResults.root"},
  fsOutputFilePathRoot{},
  fsOutputFilePath{},
  fsOutputFileName{"UniFlow.root"},
  fsOutputFileMode{"RECREATE"},
  fsTaskName{"UniFlow"},
  fsOutputFileFormat{"pdf"},
  fsGlobalProfNameLabel{},
  flFlow{},
  flQACharged{nullptr},
  flQAPID{nullptr},
  flQAPhi{nullptr},
  flQAV0s{nullptr},
  fvTasks{}
{
  // default constructor

}
//_____________________________________________________________________________
ProcessUniFlow::~ProcessUniFlow()
{
  // default destructor
  if(ffInputFile) { delete ffInputFile; }
  if(ffOutputFile) { delete ffOutputFile; }
  if(ffFitsFile) { delete ffFitsFile; }

  for(Int_t i(0); i < kUnknown; i++) { if(flFlow[i]) delete flFlow[i]; }

  // deleting the FlowTasks
  const Int_t iNumTasks = fvTasks.size();
  for(Int_t index(0); index < iNumTasks; ++index) { delete fvTasks.at(index); }
}
//_____________________________________________________________________________
void ProcessUniFlow::Clear()
{
  Info("Cleaning ProcessUniFlow instance","Clear");
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;

  for(Int_t i(0); i < kUnknown; ++i) { flFlow[i] = nullptr; }

  const Short_t iNumTasks = fvTasks.size();
  for(Short_t index(0); index < iNumTasks; index++)
  {
    if(fvTasks.at(index)) delete fvTasks.at(index);
  }
  Info("Cleaning done!","Clear");
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Run()
{
  gStyle->SetOptFit(1100);

  // main body of the class
  if(!Initialize()) { Fatal("Task not initialized","Run"); return kFALSE; }

  const Int_t iNumTasks = fvTasks.size();

  Info("===== Running over tasks ======","Run");
  Info(Form("  Number of tasks: %d\n",iNumTasks),"Run");
  for(Int_t iTask(0); iTask < iNumTasks; iTask++)
  {
    FlowTask* currentTask = fvTasks.at(iTask);
    if(!currentTask) { continue; }
    if(!InitTask(currentTask)) { return kFALSE; }
    if(!ProcessTask(currentTask)) { return kFALSE; }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Initialize()
{
  // initialization of all necessery prerequisits
  // if(fiNumMultBins < 1) { Warning(Form("Not enough mult. bins: %d (at least 1 needed)! Switching to fRebinning = 0",fiNumMultBins),"Initialize"); return kFALSE; }

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
  ffDesampleFile = TFile::Open(Form("%s/desampling.root",fsOutputFilePath.Data()),fsOutputFileMode.Data());
  if(!ffDesampleFile) { Fatal(Form("Output desampling file '%s/desampling.root' not open!","Initialize")); return kFALSE; }

  // creating output file for fits
  ffFitsFile = TFile::Open(Form("%s/fits.root",fsOutputFilePath.Data()),fsOutputFileMode.Data());
  if(!ffFitsFile) { Fatal(Form("Output desampling file '%s/fits.root' not open!","Initialize")); return kFALSE; }

  Info("Files loaded","Initialize");

  if(!LoadLists()) return kFALSE;
  Info("Flow lists loaded","Initialize");

  // initialization succesfull
  Info("Initialization succesfull","Initialize");
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::LoadLists()
{
  // loading TLists into task
  if(!ffInputFile) { Fatal("Input file does not exists!","LoadLists"); return kFALSE; }
  ffInputFile->cd(fsTaskName.Data());

  for(Int_t spec(0); spec < kUnknown; ++spec) {
    flFlow[spec] = (TList*) gDirectory->Get(Form("Flow_%s_%s",ProcessUniFlow::GetSpeciesName(PartSpecies(spec)).Data(),fsTaskName.Data()));
    if(!flFlow[spec]) { Fatal(Form("flFlow_%s list does not exists!",GetSpeciesName(PartSpecies(spec)).Data()),"LoadLists"); ffInputFile->ls(); return kFALSE; }
  }

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
TString ProcessUniFlow::GetSpeciesString(PartSpecies species)
{
    TString name = TString("k");
    name.Append(GetSpeciesName(species));
    return name;
}
//_____________________________________________________________________________
TString ProcessUniFlow::GetSpeciesName(PartSpecies species)
{
  TString name = TString();
  switch (species)
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
TString ProcessUniFlow::GetSpeciesLabel(PartSpecies species)
{
  TString label = TString();
  switch (species)
  {
    case kRefs : label.Append("RFP"); break;
    case kCharged : label.Append("h^{#pm}"); break;
    case kPion : label.Append("#pi^{#pm}"); break;
    case kKaon : label.Append("K^{#pm}"); break;
    case kProton : label.Append("p/#bar{p}"); break;
    case kPhi : label.Append("#phi"); break;
    case kK0s : label.Append("K_{S}^{0}"); break;
    case kLambda : label.Append("#Lambda/#bar{#Lambda}"); break;
    default: label.Append("N/A");
  }

  return label;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::InitTask(FlowTask* task)
{
    if(!task) { Error("Task not valid!","InitTask"); return kFALSE; }

    // checking fNumSamplesRefs
    if(task->fNumSamplesRefs != task->fNumSamples && !IsSpeciesReconstructed(task->fSpecies)) {
        Warning(Form("fNumSamplesRefs set to fNumSamples (was %d | now %d)",task->fNumSamplesRefs,task->fNumSamples),"InitTask");
        task->fNumSamplesRefs = task->fNumSamples;
    }

    if(fiNumMultBins < 1) { task->fRebinning = kFALSE; }

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

  PartSpecies spec = task->fSpecies;

  // task checks & initialization
  if(spec != kRefs && task->fNumPtBins < 1) { Error(Form("Too small number of bins: %d (at least 1 needed)!",task->fNumPtBins),"ProcesTask"); return kFALSE; }
  if(!task->HasGap() && task->fMergePosNeg) { task->fMergePosNeg = kFALSE; Warning("Merging Pos&Neg 'fMergePosNeg' switch off (no gap)","ProcessTask"); }

  // processing mixed
  if(task->fDoCorrMixed) {
    TList* listSlicesProfiles = task->fListProfiles;
    TList* listSlicesHistos = task->fListHistos;

    if(!PrepareSlicesNew(task,task->fMixedDiff)) { Error("Preparing slices failed!","ProcessTask"); return kFALSE; }

    ffOutputFile->cd();
    if(fSaveInterSteps) {
      listSlicesProfiles->Write("MakeProfileSlices",TObject::kSingleKey);
      listSlicesHistos->Write("MakeHistosSlices",TObject::kSingleKey);
    }

    if(!ProcessMixed(task)) { Error(Form("ProcessMixed '%s' failed!",task->fMixedDiff.Data()),"ProcessTask"); return kFALSE; }
  }

  if(task->fDoFour || task->fDoSix || task->fDoEight)
  {
    if(spec == kRefs && !ProcessFMC(task)) {
      Error(Form("Task '%s' (%s) not processed correctly!",task->fName.Data(), GetSpeciesName(spec).Data()),"ProcessTask");
      return kFALSE;
    }
  }
  // processing standard cumulants
  if(task->fCumOrderMax > 0) {
    if(spec == kRefs && !ProcessRefs(task)) {
      Error(Form("Task '%s' (%s) not processed correctly!",task->fName.Data(), GetSpeciesName(spec).Data()),"ProcessTask");
      return kFALSE;
    }

    if(IsSpeciesDirect(spec)) {
      for(Int_t binMult(0); binMult < fiNumMultBins; ++binMult) {
        if(!ProcessDirect(task,binMult)) {
          Error(Form("Task '%s' (%s; mult. bin %d) not processed correctly!",task->fName.Data(),GetSpeciesName(spec).Data(),binMult),"ProcessTask");
          return kFALSE;
        }
      }
    }

    if(IsSpeciesReconstructed(spec) && !ProcessReconstructed(task,0)) {
      Error(Form("Task '%s' (%s) not processed correctly!",task->fName.Data(),GetSpeciesName(spec).Data()),"ProcessTask");
      return kFALSE;
    }
  }

  if(task->fBaseCentBin > -1.0 && !ProcessSubtraction(task)){
    Error(Form("Subtraction in task '%s' (%s) not processed correctly!",task->fName.Data(),GetSpeciesName(spec).Data()),"ProcessTask");
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessFMC(FlowTask* task)
{
  if(!task) { Error("Task not valid!","ProcessFMC"); return kFALSE; }
  if(task->fNumSamplesRefs < 2) { Error("Implemented only for more samples. Terminating!","ProcessFMC"); return kFALSE; }

  Bool_t bDoSix = task->fDoSix;
  Bool_t bDoEight = task->fDoEight;
  Bool_t bDoFour = task->fDoFour;

  Int_t nSameIdx = 0;

  TList* listCorFour_12 = new TList();
  TString nameCorFour_12 = Form("Refs_pCor4_harm%d%d",task->fHarm[0],task->fHarm[1]);
  TString sProfFourName_12 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[0],task->fHarm[1],task->fHarm[0],task->fHarm[1]);
  if(task->HasGap()) {
    sProfFourName_12 += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TList* listCorFour_13 = new TList();
  TString nameCorFour_13 = Form("Refs_pCor4_harm%d%d",task->fHarm[0],task->fHarm[2]);
  TString sProfFourName_13 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[0],task->fHarm[2],task->fHarm[0],task->fHarm[2]);
  if(task->HasGap()) {
    sProfFourName_13 += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TList* listCorFour_23 = new TList();
  TString nameCorFour_23 = Form("Refs_pCor4_harm%d%d",task->fHarm[1],task->fHarm[2]);
  TString sProfFourName_23 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[1],task->fHarm[2],task->fHarm[1],task->fHarm[2]);
  if(task->HasGap()) {
    sProfFourName_23 += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TList* listCorTwo_1 = new TList();
  TString nameCorTwo_1 = Form("Refs_pCor2_harm%d",task->fHarm[0]);
  TString sProfTwoName_1 = Form("<<2>>(%d,-%d)",task->fHarm[0],task->fHarm[0]);
  if(task->HasGap()) {
    sProfTwoName_1 += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TList* listCorTwo_2 = new TList();
  TString nameCorTwo_2 = Form("Refs_pCor2_harm%d",task->fHarm[1]);
  TString sProfTwoName_2 = Form("<<2>>(%d,-%d)",task->fHarm[1],task->fHarm[1]);
  if(task->HasGap()) {
    sProfTwoName_2 += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TList* listCorTwo_3 = new TList();
  TString nameCorTwo_3 = Form("Refs_pCor2_harm%d",task->fHarm[2]);
  TString sProfTwoName_3 = Form("<<2>>(%d,-%d)",task->fHarm[2],task->fHarm[2]);
  if(task->HasGap()) {
    sProfTwoName_3 += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TList* listCorTwo_1_gap = new TList();
  TString nameCorTwo_1_gap = Form("Refs_pCor2_harm%d_gap",task->fHarm[0]);
  TString sProfTwoName_1_gap = Form("<<2>>(%d,-%d)_2sub(%.2g)",task->fHarm[0],task->fHarm[0],1.0);

  TList* listCorTwo_2_gap = new TList();
  TString nameCorTwo_2_gap = Form("Refs_pCor2_harm%d_gap",task->fHarm[1]);
  TString sProfTwoName_2_gap = Form("<<2>>(%d,-%d)_2sub(%.2g)",task->fHarm[1],task->fHarm[1],1.0);

  TList* listCorTwo_3_gap = new TList();
  TString nameCorTwo_3_gap = Form("Refs_pCor2_harm%d_gap",task->fHarm[2]);
  TString sProfTwoName_3_gap = Form("<<2>>(%d,-%d)_2sub(%.2g)",task->fHarm[2],task->fHarm[2],1.0);

  TH1D* fourFmc = nullptr;
  TList* listFmcFour = new TList();
  TString nameFmcFour = Form("Refs_FMC4_harm%d%d",task->fHarm[0],task->fHarm[2]);

  if(bDoFour)
  {
      if( task->fHarm[0] != -(task->fHarm[2]) || task->fHarm[1] != -(task->fHarm[3]) ) return kFALSE;
      if( task->fHarm[0] >= task->fHarm[1] ) { Error("Implemented only for different moments. Terminating!","ProcessFMC - 4pc"); return kFALSE; }

      nameFmcFour = Form("Refs_FMC4_harm%d%d",task->fHarm[0],task->fHarm[1]);
      if(task->HasGap()) {
        nameFmcFour += Form("_2sub(%.2g)",task->fEtaGap);
      }
      for(Short_t iSample(0); iSample < task->fNumSamplesRefs; ++iSample)
      {
        TProfile* pCorFour_12 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample));
        if(!pCorFour_12) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorFour_12 = (TProfile*) pCorFour_12->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_12.Data(), iSample),fdMultBins.data());
        TProfile* pCorTwo_1 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_1.Data(), iSample));
        TProfile* pCorTwo_2 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_2.Data(), iSample));
        if(!pCorTwo_1 || !pCorTwo_2) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfTwoName_1.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        TProfile* pCorTwo_1_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_1_gap.Data(), iSample));
        TProfile* pCorTwo_2_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_2_gap.Data(), iSample));
        if(!pCorTwo_1_gap || !pCorTwo_2_gap) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",nameCorTwo_1_gap.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorTwo_1 = (TProfile*) pCorTwo_1->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_1.Data(), iSample),fdMultBins.data());
        if(task->fRebinning) pCorTwo_2 = (TProfile*) pCorTwo_2->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_2.Data(), iSample),fdMultBins.data());
        if(task->fRebinning) pCorTwo_1_gap = (TProfile*) pCorTwo_1_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_1_gap.Data(), iSample),fdMultBins.data());
        if(task->fRebinning) pCorTwo_2_gap = (TProfile*) pCorTwo_2_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_2_gap.Data(), iSample),fdMultBins.data());
        listCorFour_12->Add(pCorFour_12);
        listCorTwo_1->Add(pCorTwo_1);
        listCorTwo_2->Add(pCorTwo_2);
        listCorTwo_1_gap->Add(pCorTwo_1_gap);
        listCorTwo_2_gap->Add(pCorTwo_2_gap);
        fourFmc = CalcFourFMC(pCorFour_12, pCorTwo_1, pCorTwo_2, pCorTwo_1_gap, pCorTwo_2_gap, task);
        fourFmc->SetName(Form("%s_sample%d", nameFmcFour.Data(), iSample));
        listFmcFour->Add(fourFmc);
      }
      TProfile* pCorFour_12Merged = (TProfile*) MergeListProfiles(listCorFour_12);
      TProfile* pCorTwo_1Merged = (TProfile*) MergeListProfiles(listCorTwo_1);
      TProfile* pCorTwo_2Merged = (TProfile*) MergeListProfiles(listCorTwo_2);
      TProfile* pCorTwo_1_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_1_gap);
      TProfile* pCorTwo_2_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_2_gap);
      if(!pCorFour_12Merged || !pCorTwo_1Merged || !pCorTwo_2Merged || !pCorTwo_1_gapMerged || !pCorTwo_2_gapMerged) { Error("Merging of 'pCorFourMerged' failed!","ProcessFMC"); return kFALSE; }
      pCorFour_12Merged->SetName(Form("%s_merged", nameCorFour_12.Data()));
      pCorTwo_1Merged->SetName(Form("%s_merged", nameCorTwo_1.Data()));
      pCorTwo_2Merged->SetName(Form("%s_merged", nameCorTwo_2.Data()));
      pCorTwo_1_gapMerged->SetName(Form("%s_merged", nameCorTwo_1.Data()));
      pCorTwo_2_gapMerged->SetName(Form("%s_merged", nameCorTwo_2.Data()));
      TH1D* hFourFmcMerged = CalcFourFMC(pCorFour_12Merged, pCorTwo_1Merged, pCorTwo_2Merged, pCorTwo_1_gapMerged, pCorTwo_2_gapMerged, task);
      TH1D* hCorFourDesampled = DesampleList(listFmcFour, hFourFmcMerged, task, nameFmcFour, kFALSE);
      if(!hCorFourDesampled) { Error("Desampling 'hCorFourDesampled' failed","ProcessRefs"); return kFALSE; }
      hCorFourDesampled->SetName(nameFmcFour.Data());
      ffOutputFile->cd();
      hCorFourDesampled->Write();
      return kTRUE;
  }

  if(bDoSix)
  {
    if( task->fHarm[0] != -(task->fHarm[3]) || task->fHarm[1] != -(task->fHarm[4]) || task->fHarm[2] != -(task->fHarm[5]) ) return kFALSE;
    if( task->fHarm[0] > task->fHarm[1] || task->fHarm[1] > task->fHarm[2] ) return kFALSE;

    if( task->fHarm[0] == task->fHarm[1] ) { nSameIdx += 1; }
    if( task->fHarm[1] == task->fHarm[2] ) { nSameIdx += 2; }
    if( nSameIdx == 3) { Error("Implemented only for different moments. Terminating!","ProcessFMC - 6pc"); return kFALSE; }

    TList* listCorSix = new TList();
    TString nameCorSix = Form("Refs_pCor6_harm%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2]);
    TList* listFmcSix = new TList();
    TString nameFmcSix = Form("Refs_fMC6_harm%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2]);
    if(task->HasGap()) {
      nameFmcSix += Form("_2sub(%.2g)",task->fEtaGap);
    }
    TString sProfSixName = Form("<<6>>(%d,%d,%d,-%d,-%d,-%d)",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[0],task->fHarm[1],task->fHarm[2]);
    if(task->HasGap()) {
      sProfSixName += Form("_2sub(%.2g)",task->fEtaGap);
    }

    for(Short_t iSample(0); iSample < task->fNumSamplesRefs; ++iSample)
    {
      TProfile* pCorSix = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfSixName.Data(), iSample));
      if(!pCorSix) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfSixName.Data(), iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(task->fRebinning) pCorSix = (TProfile*) pCorSix->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorSix.Data(), iSample),fdMultBins.data());

      TProfile* pCorFour_13 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_13.Data(),iSample));
      TProfile* pCorTwo_1 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_1.Data(), iSample));
      TProfile* pCorTwo_3 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_3.Data(), iSample));
      TProfile* pCorTwo_1_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_1_gap.Data(), iSample));
      TProfile* pCorTwo_2_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_2_gap.Data(), iSample));
      TProfile* pCorTwo_3_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_3_gap.Data(), iSample));
      if(!pCorFour_13 || !pCorTwo_1 || !pCorTwo_3) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_13.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(!pCorTwo_1_gap || !pCorTwo_2_gap || !pCorTwo_3_gap) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",nameCorTwo_1_gap.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(task->fRebinning) pCorFour_13 = (TProfile*) pCorFour_13->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_13.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_1 = (TProfile*) pCorTwo_1->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_1.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_3 = (TProfile*) pCorTwo_3->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_3.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_1_gap = (TProfile*) pCorTwo_1_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_1_gap.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_2_gap = (TProfile*) pCorTwo_2_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_2_gap.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_3_gap = (TProfile*) pCorTwo_3_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_3_gap.Data(), iSample),fdMultBins.data());
      listCorSix->Add(pCorSix);
      listCorFour_13->Add(pCorFour_13);
      listCorTwo_1->Add(pCorTwo_1);
      listCorTwo_3->Add(pCorTwo_3);
      listCorTwo_1_gap->Add(pCorTwo_1_gap);
      listCorTwo_2_gap->Add(pCorTwo_2_gap);
      listCorTwo_3_gap->Add(pCorTwo_3_gap);

      if(bDoFour && nSameIdx != 0)
      {
        fourFmc = CalcFourFMC(pCorFour_13, pCorTwo_1, pCorTwo_3, pCorTwo_1_gap, pCorTwo_3_gap, task);
        fourFmc->SetName(Form("%s_sample%d", nameFmcFour.Data(), iSample));
        listFmcFour->Add(fourFmc);
      }

      TH1D* sixFmc = nullptr;
      if(nSameIdx == 1)
      {
        TProfile* pCorFour_12 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample));
        if(!pCorFour_12) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorFour_12 = (TProfile*) pCorFour_12->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_12.Data(), iSample),fdMultBins.data());
        listCorFour_12->Add(pCorFour_12);
        sixFmc = CalcSixTwoDif(pCorSix, pCorFour_12, pCorFour_13, pCorTwo_1, pCorTwo_3, pCorTwo_3_gap, task);
        sixFmc->SetName(Form("%s_sample%d", nameFmcSix.Data(), iSample));
        listFmcSix->Add(sixFmc);
      }
      if(nSameIdx == 2)
      {
        TProfile* pCorFour_23 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_23.Data(),iSample));
        if(!pCorFour_23) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_23.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorFour_23 = (TProfile*) pCorFour_23->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_23.Data(), iSample),fdMultBins.data());
        listCorFour_23->Add(pCorFour_23);
        sixFmc = CalcSixTwoDif(pCorSix, pCorFour_23, pCorFour_13, pCorTwo_3, pCorTwo_1, pCorTwo_1_gap, task);
        sixFmc->SetName(Form("%s_sample%d", nameFmcSix.Data(), iSample));
        listFmcSix->Add(sixFmc);
      }
      if(nSameIdx == 0)
      {
        TProfile* pCorFour_23 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_23.Data(),iSample));
        TProfile* pCorFour_12 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample));
        TProfile* pCorTwo_2 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_2.Data(), iSample));
        if(!pCorFour_23 || !pCorFour_23 || !pCorTwo_2) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_23.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorTwo_2 = (TProfile*) pCorTwo_2->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_2.Data(), iSample),fdMultBins.data());
        if(task->fRebinning) pCorFour_23 = (TProfile*) pCorFour_23->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_23.Data(), iSample),fdMultBins.data());
        if(task->fRebinning) pCorFour_12 = (TProfile*) pCorFour_12->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_12.Data(), iSample),fdMultBins.data());
        listCorTwo_2->Add(pCorTwo_2);
        listCorFour_12->Add(pCorFour_12);
        listCorFour_23->Add(pCorFour_23);
        sixFmc = CalcSixThreeDif(pCorSix, pCorFour_12, pCorFour_13, pCorFour_23, pCorTwo_1, pCorTwo_2, pCorTwo_3, pCorTwo_1_gap, pCorTwo_2_gap, pCorTwo_3_gap, task);
        sixFmc->SetName(Form("%s_sample%d", nameFmcSix.Data(), iSample));
        listFmcSix->Add(sixFmc);
      }

    }// end-for {iSample}: samples
    Debug("Samples processing done!","ProcessFMC");

    TProfile* pCorSixMerged = (TProfile*) MergeListProfiles(listCorSix);
    TProfile* pCorFour_13Merged = (TProfile*) MergeListProfiles(listCorFour_13);
    TProfile* pCorTwo_1Merged = (TProfile*) MergeListProfiles(listCorTwo_1);
    TProfile* pCorTwo_3Merged = (TProfile*) MergeListProfiles(listCorTwo_3);
    TProfile* pCorTwo_1_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_1_gap);
    TProfile* pCorTwo_3_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_3_gap);
    if(!pCorSixMerged || !pCorFour_13Merged || !pCorTwo_1Merged || !pCorTwo_3Merged || !pCorTwo_1_gapMerged || !pCorTwo_3_gapMerged) { Error("Merging of 'pCorSixMerged' failed!","ProcessFMC"); return kFALSE; }
    pCorSixMerged->SetName(Form("%s_merged", nameCorSix.Data()));
    pCorFour_13Merged->SetName(Form("%s_merged", nameCorFour_13.Data()));
    pCorTwo_1Merged->SetName(Form("%s_merged", nameCorTwo_1.Data()));
    pCorTwo_3Merged->SetName(Form("%s_merged", nameCorTwo_3.Data()));
    pCorTwo_1_gapMerged->SetName(Form("%s_merged", nameCorTwo_1.Data()));
    pCorTwo_3_gapMerged->SetName(Form("%s_merged", nameCorTwo_3.Data()));
    TH1D* hSixFMC_Merged = nullptr;
    if(nSameIdx == 1)
    {
      TProfile* pCorFour_12Merged = (TProfile*) MergeListProfiles(listCorFour_12);
      if(!pCorFour_12Merged) { Error("Merging of 'pCorFour_12Merged' failed!","ProcessFMC"); return kFALSE; }
      pCorFour_12Merged->SetName(Form("%s_merged", nameCorFour_12.Data()));
      hSixFMC_Merged = CalcSixTwoDif(pCorSixMerged, pCorFour_12Merged, pCorFour_13Merged, pCorTwo_1Merged, pCorTwo_3Merged, pCorTwo_3_gapMerged, task);
    }
    if(nSameIdx == 2)
    {
      TProfile* pCorFour_23Merged = (TProfile*) MergeListProfiles(listCorFour_23);
      if(!pCorFour_23Merged) { Error("Merging of 'pCorFour_23Merged' failed!","ProcessFMC"); return kFALSE; }
      pCorFour_23Merged->SetName(Form("%s_merged", nameCorFour_23.Data()));
      hSixFMC_Merged = CalcSixTwoDif(pCorSixMerged, pCorFour_23Merged, pCorFour_13Merged, pCorTwo_3Merged, pCorTwo_1Merged, pCorTwo_1_gapMerged, task);
    }
    if(nSameIdx == 0)
    {
      TProfile* pCorFour_23Merged = (TProfile*) MergeListProfiles(listCorFour_23);
      TProfile* pCorFour_12Merged = (TProfile*) MergeListProfiles(listCorFour_12);
      TProfile* pCorTwo_2Merged = (TProfile*) MergeListProfiles(listCorTwo_2);
      TProfile* pCorTwo_2_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_2_gap);
      if(!pCorFour_23Merged || !pCorFour_12Merged || !pCorTwo_2Merged || !pCorTwo_2_gapMerged) { Error("Merging of 'pCorTwo_2Merged' failed!","ProcessFMC"); return kFALSE; }
      pCorFour_12Merged->SetName(Form("%s_merged", nameCorFour_12.Data()));
      pCorFour_23Merged->SetName(Form("%s_merged", nameCorFour_23.Data()));
      pCorTwo_2Merged->SetName(Form("%s_merged", nameCorTwo_2.Data()));
      pCorTwo_2_gapMerged->SetName(Form("%s_merged", nameCorTwo_2_gap.Data()));
      hSixFMC_Merged = CalcSixThreeDif(pCorSixMerged, pCorFour_12Merged, pCorFour_13Merged, pCorFour_23Merged, pCorTwo_1Merged, pCorTwo_2Merged, pCorTwo_3Merged, pCorTwo_1_gapMerged, pCorTwo_2_gapMerged, pCorTwo_3_gapMerged, task);
    }
    // desampling
    Debug("Desampling","ProcessFMC");
    TH1D* hCorSixDesampled = DesampleList(listFmcSix, hSixFMC_Merged, task, nameFmcSix, kFALSE);
    if(!hCorSixDesampled) { Error("Desampling 'hCorSixDesampled' failed","ProcessRefs"); return kFALSE; }
    hCorSixDesampled->SetName(nameFmcSix.Data());

    ffOutputFile->cd();
    hCorSixDesampled->Write();

    if(bDoFour && nSameIdx != 0)
    {
      TH1D* hFourFmcMerged = CalcFourFMC(pCorFour_13Merged, pCorTwo_1Merged, pCorTwo_3Merged, pCorTwo_1_gapMerged, pCorTwo_3_gapMerged, task);
      TH1D* hCorFourDesampled = DesampleList(listFmcFour, hFourFmcMerged, task, nameFmcFour, kFALSE);
      if(!hCorFourDesampled) { Error("Desampling 'hCorFourDesampled' failed","ProcessRefs"); return kFALSE; }
      hCorFourDesampled->SetName(nameFmcFour.Data());
      ffOutputFile->cd();
      hCorFourDesampled->Write();
    }

    return kTRUE;

  }//end doSix

  if(bDoEight)
  {
    if( task->fHarm[0] != -(task->fHarm[4]) || task->fHarm[1] != -(task->fHarm[5]) || task->fHarm[2] != -(task->fHarm[6]) || task->fHarm[3] != -(task->fHarm[7])) return kFALSE;
    if( task->fHarm[0] > task->fHarm[1] || task->fHarm[1] > task->fHarm[2] || task->fHarm[2] > task->fHarm[3]) return kFALSE;

    if( task->fHarm[0] == task->fHarm[1] ) { nSameIdx += 1; }
    if( task->fHarm[1] == task->fHarm[2] ) { nSameIdx += 2; }
    if( task->fHarm[2] == task->fHarm[3] ) { nSameIdx += 4; }
    if( nSameIdx == 7 ) { Error("Implemented only for different moments. Terminating!","ProcessFMC - 8pc"); return kFALSE; }
    if( nSameIdx < 3 || nSameIdx == 4 ) { Error("Not implemented for this combinations so far. Terminating!","ProcessFMC - 8pc"); return kFALSE; }

    TList* listCorEight = new TList();
    TString nameCorEight = Form("Refs_pCor8_harm%d%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]);
    TString sProfEightName = Form("<<8>>(%d,%d,%d,%d,-%d,-%d,-%d,-%d)",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3],task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]);
    if(task->HasGap()) {
      sProfEightName += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listFmcEight = new TList();
    TString nameFmcEight = Form("Refs_fMC8_harm%d%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]);
    if(task->HasGap()) {
      nameFmcEight += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorSix_123 = new TList();
    TString nameCorSix_123 = Form("Refs_pCor6_harm%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2]);
    TString sProfSixName_123 = Form("<<6>>(%d,%d,%d,-%d,-%d,-%d)",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[0],task->fHarm[1],task->fHarm[2]);
    if(task->HasGap()) {
      sProfSixName_123 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorSix_124 = new TList();
    TString nameCorSix_124 = Form("Refs_pCor6_harm%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[3]);
    TString sProfSixName_124 = Form("<<6>>(%d,%d,%d,-%d,-%d,-%d)",task->fHarm[0],task->fHarm[1],task->fHarm[3],task->fHarm[0],task->fHarm[1],task->fHarm[3]);
    if(task->HasGap()) {
      sProfSixName_124 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorSix_234 = new TList();
    TString nameCorSix_234 = Form("Refs_pCor6_harm%d%d%d",task->fHarm[1],task->fHarm[2],task->fHarm[3]);
    TString sProfSixName_234 = Form("<<6>>(%d,%d,%d,-%d,-%d,-%d)",task->fHarm[1],task->fHarm[2],task->fHarm[3],task->fHarm[1],task->fHarm[2],task->fHarm[3]);
    if(task->HasGap()) {
      sProfSixName_234 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorSix_134 = new TList();
    TString nameCorSix_134 = Form("Refs_pCor6_harm%d%d%d",task->fHarm[0],task->fHarm[2],task->fHarm[3]);
    TString sProfSixName_134 = Form("<<6>>(%d,%d,%d,-%d,-%d,-%d)",task->fHarm[0],task->fHarm[2],task->fHarm[3],task->fHarm[0],task->fHarm[2],task->fHarm[3]);
    if(task->HasGap()) {
      sProfSixName_134 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorFour_14 = new TList();
    TString nameCorFour_14 = Form("Refs_pCor4_harm%d%d",task->fHarm[0],task->fHarm[3]);
    TString sProfFourName_14 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[0],task->fHarm[3],task->fHarm[0],task->fHarm[3]);
    if(task->HasGap()) {
      sProfFourName_14 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorFour_13 = new TList();
    TString nameCorFour_13 = Form("Refs_pCor4_harm%d%d",task->fHarm[0],task->fHarm[2]);
    TString sProfFourName_13 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[0],task->fHarm[2],task->fHarm[0],task->fHarm[2]);
    if(task->HasGap()) {
      sProfFourName_13 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorFour_34 = new TList();
    TString nameCorFour_34 = Form("Refs_pCor4_harm%d%d",task->fHarm[2],task->fHarm[3]);
    TString sProfFourName_34 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[2],task->fHarm[3],task->fHarm[2],task->fHarm[3]);
    if(task->HasGap()) {
      sProfFourName_34 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorFour_23 = new TList();
    TString nameCorFour_23 = Form("Refs_pCor4_harm%d%d",task->fHarm[1],task->fHarm[2]);
    TString sProfFourName_23 = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarm[1],task->fHarm[1],task->fHarm[1],task->fHarm[2]);
    if(task->HasGap()) {
      sProfFourName_23 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorTwo_4 = new TList();
    TString nameCorTwo_4 = Form("Refs_pCor2_harm%d",task->fHarm[3]);
    TString sProfTwoName_4 = Form("<<2>>(%d,-%d)",task->fHarm[3],task->fHarm[3]);
    if(task->HasGap()) {
      sProfTwoName_4 += Form("_2sub(%.2g)",task->fEtaGap);
    }

    TList* listCorTwo_4_gap = new TList();
    TString nameCorTwo_4_gap = Form("Refs_pCor2_harm%d_gap",task->fHarm[3]);
    TString sProfTwoName_4_gap = Form("<<2>>(%d,-%d)_2sub(%.2g)",task->fHarm[3],task->fHarm[3],1.0);


    for(Short_t iSample(0); iSample < task->fNumSamplesRefs; ++iSample)
    {
      TProfile* pCorEight = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfEightName.Data(), iSample));
      if(!pCorEight) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfEightName.Data(), iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(task->fRebinning) pCorEight = (TProfile*) pCorEight->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorEight.Data(), iSample),fdMultBins.data());

      TProfile* pCorSix_123 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfSixName_123.Data(), iSample));
      if(!pCorSix_123) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfSixName_123.Data(), iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(task->fRebinning) pCorSix_123 = (TProfile*) pCorSix_123->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorSix_123.Data(), iSample),fdMultBins.data());

      TProfile* pCorFour_12 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample));
      TProfile* pCorTwo_1 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_1.Data(), iSample));
      TProfile* pCorTwo_1_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_1_gap.Data(), iSample));
      TProfile* pCorTwo_4 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_4.Data(), iSample));
      TProfile* pCorTwo_4_gap = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName_4_gap.Data(), iSample));
      if(!pCorFour_12) { Warning(Form("Profile '%s' not valid! 8FMC",Form("%s_Pos_sample%d",sProfFourName_12.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(!pCorTwo_1) { Warning(Form("Profile '%s' not valid! 8FMC",Form("%s_Pos_sample%d",sProfTwoName_1.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(!pCorTwo_1_gap) { Warning(Form("Profile '%s' not valid! 8FMC",Form("%s_Pos_sample%d",sProfTwoName_1_gap.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(!pCorTwo_4) { Warning(Form("Profile '%s' not valid! 8FMC",Form("%s_Pos_sample%d",sProfTwoName_4.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(!pCorTwo_4_gap) { Warning(Form("Profile '%s' not valid! 8FMC",Form("%s_Pos_sample%d",sProfTwoName_4_gap.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
      if(task->fRebinning) pCorFour_12 = (TProfile*) pCorFour_12->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_12.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_1 = (TProfile*) pCorTwo_1->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_1.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_1_gap = (TProfile*) pCorTwo_1_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_1_gap.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_4 = (TProfile*) pCorTwo_4->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_4.Data(), iSample),fdMultBins.data());
      if(task->fRebinning) pCorTwo_4_gap = (TProfile*) pCorTwo_4_gap->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo_4_gap.Data(), iSample),fdMultBins.data());
      listCorEight->Add(pCorEight);
      listCorSix_123->Add(pCorSix_123);
      listCorFour_12->Add(pCorFour_12);
      listCorTwo_1->Add(pCorTwo_1);
      listCorTwo_1_gap->Add(pCorTwo_1_gap);
      listCorTwo_4->Add(pCorTwo_4);
      listCorTwo_4_gap->Add(pCorTwo_4_gap);

      TH1D* eightFmc = nullptr;
      if(nSameIdx == 3)
      {
        TProfile* pCorSix_124 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfSixName_124.Data(), iSample));
        if(!pCorSix_124) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfSixName_124.Data(), iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorSix_124 = (TProfile*) pCorSix_124->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorSix_124.Data(), iSample),fdMultBins.data());

        TProfile* pCorFour_14 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_14.Data(),iSample));
        if(!pCorFour_14) { Warning(Form("Profile '%s' not valid! 8FMC",Form("%s_Pos_sample%d",sProfFourName_14.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorFour_14 = (TProfile*) pCorFour_14->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_14.Data(), iSample),fdMultBins.data());

        listCorSix_124->Add(pCorSix_124);
        listCorFour_14->Add(pCorFour_14);
        eightFmc = CalcEight_ThreeOne(pCorEight, pCorSix_123, pCorSix_124, pCorFour_12, pCorFour_14, pCorTwo_1, pCorTwo_4, pCorTwo_4_gap, task);
      }
      if(nSameIdx == 6)
      {
        TProfile* pCorSix_234 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfSixName_234.Data(), iSample));
        if(!pCorSix_234) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfSixName_234.Data(), iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorSix_234 = (TProfile*) pCorSix_234->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorSix_234.Data(), iSample),fdMultBins.data());

        TProfile* pCorFour_23 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_23.Data(),iSample));
        if(!pCorFour_23) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_23.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorFour_23 = (TProfile*) pCorFour_23->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_23.Data(), iSample),fdMultBins.data());

        listCorSix_234->Add(pCorSix_234);
        listCorFour_23->Add(pCorFour_23);
        eightFmc = CalcEight_ThreeOne(pCorEight, pCorSix_234, pCorSix_123, pCorFour_23, pCorFour_12, pCorTwo_4, pCorTwo_1, pCorTwo_1_gap, task);
      }
      if(nSameIdx == 5)
      {
        TProfile* pCorSix_134 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfSixName_134.Data(), iSample));
        if(!pCorSix_134) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfSixName_134.Data(), iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorSix_134 = (TProfile*) pCorSix_134->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorSix_134.Data(), iSample),fdMultBins.data());

        TProfile* pCorFour_13 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_13.Data(),iSample));
        TProfile* pCorFour_34 = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName_34.Data(),iSample));
        if(!pCorFour_13 || !pCorFour_34) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName_13.Data(),iSample)),"ProcessFMC"); flFlow[kRefs]->ls(); return kFALSE; }
        if(task->fRebinning) pCorFour_13 = (TProfile*) pCorFour_13->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_13.Data(), iSample),fdMultBins.data());
        if(task->fRebinning) pCorFour_34 = (TProfile*) pCorFour_34->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour_34.Data(), iSample),fdMultBins.data());

        listCorSix_134->Add(pCorSix_134);
        listCorFour_13->Add(pCorFour_13);
        listCorFour_34->Add(pCorFour_34);
        eightFmc = CalcEight_TwoTwo(pCorEight, pCorSix_123, pCorSix_134, pCorFour_12, pCorFour_13, pCorFour_34, pCorTwo_1, pCorTwo_4, task);
      }
      eightFmc->SetName(Form("%s_sample%d", nameFmcEight.Data(), iSample));
      listFmcEight->Add(eightFmc);
    }// end-for {iSample}: samples
    Debug("Samples processing done!","ProcessFMC");

    TProfile* pCorEightMerged = (TProfile*) MergeListProfiles(listCorEight);
    TProfile* pCorSix_123Merged = (TProfile*) MergeListProfiles(listCorSix_123);
    TProfile* pCorFour_12Merged = (TProfile*) MergeListProfiles(listCorFour_12);
    TProfile* pCorTwo_1Merged = (TProfile*) MergeListProfiles(listCorTwo_1);
    TProfile* pCorTwo_4Merged = (TProfile*) MergeListProfiles(listCorTwo_4);
    TProfile* pCorTwo_1_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_1_gap);
    TProfile* pCorTwo_4_gapMerged = (TProfile*) MergeListProfiles(listCorTwo_4_gap);
    if(!pCorEightMerged || !pCorSix_123Merged || !pCorFour_12Merged || !pCorTwo_1Merged || !pCorTwo_4Merged || !pCorTwo_1_gapMerged || !pCorTwo_4_gapMerged) { Error("Merging of 'pCorEightMerged' failed!","ProcessFMC"); return kFALSE; }
    pCorEightMerged->SetName(Form("%s_merged", nameCorEight.Data()));
    pCorSix_123Merged->SetName(Form("%s_merged", nameCorSix_123.Data()));
    pCorFour_12Merged->SetName(Form("%s_merged", nameCorFour_12.Data()));
    pCorTwo_1Merged->SetName(Form("%s_merged", nameCorTwo_1.Data()));
    pCorTwo_4Merged->SetName(Form("%s_merged", nameCorTwo_4.Data()));
    pCorTwo_1_gapMerged->SetName(Form("%s_merged", nameCorTwo_1_gap.Data()));
    pCorTwo_4_gapMerged->SetName(Form("%s_merged", nameCorTwo_4_gap.Data()));

    TH1D* hEightFMC_Merged = nullptr;
    if(nSameIdx == 3)
    {
      TProfile* pCorSix_124Merged = (TProfile*) MergeListProfiles(listCorSix_124);
      TProfile* pCorFour_14Merged = (TProfile*) MergeListProfiles(listCorFour_14);
      if(!pCorSix_124Merged || !pCorFour_14Merged) { Error("Merging of 'pCorSix_124Merged' failed!","ProcessFMC"); return kFALSE; }
      pCorSix_124Merged->SetName(Form("%s_merged", nameCorSix_124.Data()));
      pCorFour_14Merged->SetName(Form("%s_merged", nameCorFour_14.Data()));
      hEightFMC_Merged = CalcEight_ThreeOne(pCorEightMerged, pCorSix_123Merged, pCorSix_124Merged, pCorFour_12Merged, pCorFour_14Merged, pCorTwo_1Merged, pCorTwo_4Merged, pCorTwo_4_gapMerged, task);
    }
    if(nSameIdx == 6)
    {
      TProfile* pCorSix_234Merged = (TProfile*) MergeListProfiles(listCorSix_234);
      TProfile* pCorFour_23Merged = (TProfile*) MergeListProfiles(listCorFour_23);
      if(!pCorSix_234Merged || !pCorFour_23Merged) { Error("Merging of 'pCorSix_234Merged' failed!","ProcessFMC"); return kFALSE; }
      pCorSix_234Merged->SetName(Form("%s_merged", nameCorSix_234.Data()));
      pCorFour_23Merged->SetName(Form("%s_merged", nameCorFour_23.Data()));
      hEightFMC_Merged = CalcEight_ThreeOne(pCorEightMerged, pCorSix_234Merged, pCorSix_123Merged, pCorFour_23Merged, pCorFour_12Merged, pCorTwo_4Merged, pCorTwo_1Merged, pCorTwo_1_gapMerged, task);
    }
    if(nSameIdx == 5)
    {
      TProfile* pCorSix_134Merged = (TProfile*) MergeListProfiles(listCorSix_134);
      TProfile* pCorFour_13Merged = (TProfile*) MergeListProfiles(listCorFour_13);
      TProfile* pCorFour_34Merged = (TProfile*) MergeListProfiles(listCorFour_34);
      if(!pCorSix_134Merged || !pCorFour_13Merged || !pCorFour_34Merged) { Error("Merging of 'pCorSix_134Merged' failed!","ProcessFMC"); return kFALSE; }
      pCorSix_134Merged->SetName(Form("%s_merged", nameCorSix_134.Data()));
      pCorFour_13Merged->SetName(Form("%s_merged", nameCorFour_13.Data()));
      pCorFour_34Merged->SetName(Form("%s_merged", nameCorFour_34.Data()));
      hEightFMC_Merged = CalcEight_TwoTwo(pCorEightMerged, pCorSix_123Merged, pCorSix_134Merged, pCorFour_12Merged, pCorFour_13Merged, pCorFour_34Merged, pCorTwo_1Merged, pCorTwo_4Merged, task);
    }
    // desampling
    Debug("Desampling","ProcessFMC");
    TH1D* hCorEightDesampled = DesampleList(listFmcEight, hEightFMC_Merged, task, nameFmcEight, kFALSE);
    if(!hCorEightDesampled) { Error("Desampling 'hCorEightDesampled' failed","ProcessRefs"); return kFALSE; }
    hCorEightDesampled->SetName(nameFmcEight.Data());

    ffOutputFile->cd();
    hCorEightDesampled->Write();
    return kTRUE;

  }//end doEight

  if(listCorFour_12) delete listCorFour_12;
  if(listCorFour_13) delete listCorFour_13;
  if(listCorFour_23) delete listCorFour_23;
  if(listCorTwo_1) delete listCorTwo_1;
  if(listCorTwo_2) delete listCorTwo_2;
  if(listCorTwo_3) delete listCorTwo_3;
  if(listCorTwo_1_gap) delete listCorTwo_1_gap;
  if(listCorTwo_2_gap) delete listCorTwo_2_gap;
  if(listCorTwo_3_gap) delete listCorTwo_3_gap;

  return kFALSE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessMixed(FlowTask* task)
{
  if(!task) { Error("Task not valid!","ProcessMixed"); return kFALSE; }
  Debug("Processing mixed","ProcessMixed");

  PartSpecies species = task->fSpecies;
  Bool_t bReco = IsSpeciesReconstructed(species);

    if(bReco && task->fNumSamples > 1) {
        Error("Implemented only for 1 sample. Terminating!","ProcessMixed");
        Warning("TO DO","ProcessMixed");
        return kFALSE;
    }

  Int_t iSample = 0;
  TString sNameRefs = Form("%s_Pos_sample%d", task->fMixedRefs.Data(),iSample);
  TString sNamePOIs = Form("%s_Pos_sample%d", task->fMixedDiff.Data(), iSample);
  TString sNamePOIsNeg = Form("%s_Neg_sample%d", task->fMixedDiff.Data(), iSample);

  TList trashCol;
  trashCol.SetOwner(kTRUE);

  // ### Preparing Refs ###
  TList* listRefs = LoadSamples(flFlow[kRefs], Form("%s_Pos", task->fMixedRefs.Data()), task->fNumSamplesRefs);
  if(!listRefs) { Error("Loading list with Refs samples failed","ProcessMixed"); return kFALSE; }

  // TProfile* profRef_preRebin = (TProfile*) flFlow[kRefs]->FindObject(sNameRefs.Data());
  TProfile* profRef_preRebin = (TProfile*) MergeListProfiles(listRefs);
  if(!profRef_preRebin) { Error(Form("Refs profile '%s' pre-rebin not found!",sNameRefs.Data()),"ProcessMixed"); flFlow[kRefs]->ls(); return kFALSE; }
  TProfile* profRef = (TProfile*) profRef_preRebin->Rebin(fiNumMultBins,Form("%s_rebin",sNameRefs.Data()),fdMultBins.data());
  if(!profRef) { Error("Refs profile rebinning failed!","ProcessMixed"); return kFALSE; }
  trashCol.Add(profRef);

   // ### Preparing POIs ###
  for(Int_t iMultBin(0); iMultBin < fiNumMultBins; ++iMultBin) {

    TH1D* histFlow = nullptr; // histo with final results
    // TString sHistFlowName = Form("%s_%s_mult%d",GetSpeciesName(task->fSpecies).Data(),sNamePOIs.Data(),iMultBin);
    TString sHistFlowName = Form("%s_%s_mult%d",GetSpeciesName(task->fSpecies).Data(),task->fMixedDiff.Data(),iMultBin);

    // direct species
    if(!bReco) {
      TString sName = Form("%s_Pos_sample0_mult%d",task->fMixedDiff.Data(),iMultBin);
      TProfile* profVn = (TProfile*) task->fListProfiles->FindObject(sName.Data());
      if(!profVn) { Error("Loading slice failed!","ProcessMixed"); task->fListProfiles->ls(); return kFALSE; }

      // Making vn out of cn,dn
      TH1D* histVn = (TH1D*) profVn->ProjectionX();
      trashCol.Add(histVn);

      // dividing POIS / sqrt(refs)
      Double_t dRefCont = profRef->GetBinContent(iMultBin+1);
      Double_t dRefErr = profRef->GetBinError(iMultBin+1);

      for(Int_t bin(0); bin < profVn->GetNbinsX()+1; ++bin) {
        if(dRefCont < 0.0) {
          histVn->SetBinContent(bin, 9999.9);
          histVn->SetBinError(bin, 9999.9);
          continue;
        }

        Double_t dOldCont = histVn->GetBinContent(bin);
        Double_t dOldErr = histVn->GetBinError(bin);

        Double_t dNewCont = dOldCont / TMath::Sqrt(dRefCont);
        Double_t dNewErrSq = dOldErr*dOldErr/dRefCont + 0.25*TMath::Power(dRefCont,-3.0)*dOldCont*dOldCont*dRefErr*dRefErr;

        histVn->SetBinContent(bin, dNewCont);
        histVn->SetBinError(bin, TMath::Sqrt(dNewErrSq));
      }

      histFlow = histVn;
    } else { // end-if {!bReco}
      // recnstructed
      // gSystem->mkdir(Form("%s/fits/",fsOutputFilePath.Data()));

      histFlow = new TH1D(sHistFlowName.Data(),Form("%s: %s; #it{p}_{T} (GeV/#it{c});",GetSpeciesLabel(task->fSpecies).Data(),sHistFlowName.Data()), task->fNumPtBins,task->fPtBinsEdges.data());
      if(!histFlow) { Error("Creation of 'histFlow' failed!","ProcessMixed"); return kFALSE; }
      trashCol.Add(histFlow);

      // dividing POIS / sqrt(refs)
      Double_t dRefCont = profRef->GetBinContent(iMultBin+1);
      Double_t dRefErr = profRef->GetBinError(iMultBin+1);

      for(Int_t iPtBin(0); iPtBin < task->fNumPtBins; ++iPtBin) {

        TList* listFits = new TList();
        // listFits->SetOwner(kTRUE); // NB: when on, seg fault happen

        // Adding description strings
        TNamed* corr = new TNamed("corr",Form("%s",task->fMixedDiff.Data()));
        listFits->Add(corr);
        TNamed* spec = new TNamed("spec",Form("%s",GetSpeciesLabel(task->fSpecies).Data()));
        listFits->Add(spec);
        TNamed* nCent = new TNamed("cent",Form("%d-%d",(Int_t)fdMultBins[iMultBin],(Int_t)fdMultBins[iMultBin+1]));
        listFits->Add(nCent);
        TNamed* nPt = new TNamed("pt",Form("%g-%g",task->fPtBinsEdges[iPtBin],task->fPtBinsEdges[iPtBin+1]));
        listFits->Add(nPt);

        TH1D* hInvMass = (TH1D*) task->fListHistos->FindObject(Form("hInvMass_mult%d_pt%d",iMultBin,iPtBin));
        if(!hInvMass) { Error("Loading inv. mass slice failed!","ProcessMixed"); task->fListHistos->ls(); return kFALSE; }

        TH1D* hInvMassBg = nullptr;
        if(species == kPhi) {
          hInvMassBg = (TH1D*) task->fListHistos->FindObject(Form("hInvMassBg_mult%d_pt%d",iMultBin,iPtBin));
          if(!hInvMassBg) { Error("Loading inv. mass (Bg) slice failed!","ProcessMixed"); task->fListHistos->ls(); return kFALSE; }
        }

        TString sName = Form("%s_Pos_sample0_mult%d_pt%d",task->fMixedDiff.Data(),iMultBin,iPtBin);
        TProfile* profVn = (TProfile*) task->fListProfiles->FindObject(sName.Data());
        if(!profVn) { Error("Loading correlation slice failed!","ProcessMixed"); task->fListProfiles->ls(); return kFALSE; }

        // Making vn out of cn,dn
        TH1D* histVn = (TH1D*) profVn->ProjectionX();

        // for(Int_t bin(0); bin < profVn->GetNbinsX()+1; ++bin) {
        //   if(dRefCont < 0.0) {
        //     histVn->SetBinContent(bin, 9999.9);
        //     histVn->SetBinError(bin, 9999.9);
        //     continue;
        //   }
        //
        //   Double_t dOldCont = histVn->GetBinContent(bin);
        //   Double_t dOldErr = histVn->GetBinError(bin);
        //
        //   Double_t dNewCont = dOldCont / TMath::Sqrt(dRefCont);
        //   Double_t dNewErrSq = dOldErr*dOldErr/dRefCont + 0.25*TMath::Power(dRefCont,-3.0)*dOldCont*dOldCont*dRefErr*dRefErr;
        //
        //   histVn->SetBinContent(bin, dNewCont);
        //   histVn->SetBinError(bin, TMath::Sqrt(dNewErrSq));
        // }

        // Here ready for fitting

        TCanvas* canFitInvMass = new TCanvas("canFitInvMass","canFitInvMass",1600,1200); // canvas for fitting results

        Bool_t bFitMass = kFALSE;
        Bool_t bFitFlow = kFALSE;

        Double_t dFlow = 0.0;
        Double_t dFlowError = 0.0;

        TF1 fitMass, fitMassSig, fitMassBg, fitFlow, fitFlowSig, fitFlowBg;

        bFitMass = FitInvMass(hInvMass, task, fitMass, fitMassSig, fitMassBg, listFits, hInvMassBg);

        if(bFitMass) {
            bFitFlow = FitCorrelations(histVn, task, fitFlow, fitFlowSig, fitFlowBg, fitMassSig, fitMassBg, listFits);
        }

        // if either FitInvMass or FitCorrelations fails, terminate here!
        // NB: It is important to save output listFist first for debugging

        ffFitsFile->cd();
        listFits->Write(Form("%s_%s_cent%d_pt%d",GetSpeciesName(task->fSpecies).Data(),task->fMixedDiff.Data(),iMultBin,iPtBin),TObject::kSingleKey);

        if(!bFitMass) {
          Error(Form("Fitting inv.mass unsuccesfull (mult %d | pt %d)",iMultBin,iPtBin),"ProcessMixed");
          delete canFitInvMass;
          // delete listFits;
          return kFALSE;
        }

        if(!bFitFlow) {
          Error(Form("Fitting vn-mass unsuccesfull (mult %d | pt %d)",iMultBin,iPtBin),"ProcessMixed");
          delete canFitInvMass;
          // delete listFits;
          return kFALSE;
        }

        Int_t iParFlow = fitFlowSig.GetNpar() - 1;
        dFlow = fitFlowSig.GetParameter(iParFlow);
        dFlowError = fitFlowSig.GetParError(iParFlow);
        // histFlow->SetBinContent(iPtBin+1,dFlow);
        // histFlow->SetBinError(iPtBin+1,dFlowError);

        Double_t dFlowRel = -999.9; if(TMath::Abs(dFlow) > 0.0) { dFlowRel = dFlowError / dFlow; }
        Info(Form("Final corr(n,m,k): (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,iPtBin,dFlow,dFlowError,dFlowRel), "ProcessMixed");

        Double_t dNewCont = dFlow / TMath::Sqrt(dRefCont);
        Double_t dNewErrSq = dFlowError*dFlowError/dRefCont + 0.25*TMath::Power(dRefCont,-3.0)*dFlow*dFlow*dRefErr*dRefErr;
        Double_t dNewErr = TMath::Sqrt(dNewErrSq);

        histFlow->SetBinContent(iPtBin+1,dNewCont);
        histFlow->SetBinError(iPtBin+1,dNewErr);
        Info(Form("Final v(n,m,k): (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,iPtBin,dNewCont,dNewErr,dNewErr/dNewCont), "ProcessMixed");

        // // === Plotting fits ===
        // TLatex latex2;
        // // latex2.SetTextFont(43);
        // // latex2.SetTextSize(40);
        // latex2.SetNDC();
        //
        // canFitInvMass->cd(1);
        // // if(task->fSpecies == kPhi) canFitInvMass->cd(2);
        // latex2.DrawLatex(0.17,0.85,Form("#color[9]{pt %g-%g GeV/c (%g-%g%%)}",task->fPtBinsEdges[iPtBin],task->fPtBinsEdges[iPtBin+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
        // canFitInvMass->cd(2);
        // latex2.DrawLatex(0.17,0.85,Form("#color[9]{pt %g-%g GeV/c (%g-%g%%)}",task->fPtBinsEdges[iPtBin],task->fPtBinsEdges[iPtBin+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
        // canFitInvMass->SaveAs(Form("%s/fits/%s_%s_mult%d_pt%d.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),sNamePOIs.Data(),iMultBin,iPtBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
        //
        // delete canFitInvMass;
        // delete listFits;

      } // end-for {iPtBin}
    } // end-else {!bReco}

    histFlow->SetName(sHistFlowName.Data());
    histFlow->SetTitle(sHistFlowName.Data());

    ffOutputFile->cd();
    histFlow->Write();

    TCanvas* cFlow = new TCanvas("cFlow","cFlow");
    cFlow->cd();
    histFlow->SetStats(0);
    histFlow->DrawCopy();
    cFlow->SaveAs(Form("%s/Flow_%s_%s_mult%d.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),sNamePOIs.Data(),iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
    delete cFlow;
  } // end-for {iMultBin}

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessRefs(FlowTask* task)
{
  Info("Processing Refs task","ProcesRefs");
  if(!task) { Error("Task not valid!","ProcessRefs"); return kFALSE; }
  if(task->fSpecies != kRefs) { Error("Task species not kRefs!","ProcessRefs"); return kFALSE; }

  Bool_t bDoFour = (task->fCumOrderMax >= 4); // check if cn{4} should be processed
  Bool_t bCorrelated = task->fConsCorr; // check if correlated uncrt. are considered
  Bool_t bDo3sub = task->Has3sub();
  Int_t nOfSamples = task->fNumSamplesRefs;


  // for saving profiles into TList for merging (into single one) -> estimation for central values
  TList* listCorTwo = new TList(); TString nameCorTwo = Form("Refs_pCor2_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());
  TList* listCorFour = new TList(); TString nameCorFour = Form("Refs_pCor4_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());

  // 3 subevents
  TList* listCorTwo3sub[3][3] = {nullptr};
  TString nameCorTwo3sub[3][3] = {""};
  TList* listCorFour3sub[3] = {nullptr};
  TString nameCorFour3sub[3] = {""};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        nameCorTwo3sub[rf1Pos][rf2Pos] = Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]);
        listCorTwo3sub[rf1Pos][rf2Pos] = new TList();
      }
    }
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      nameCorFour3sub[twoPos] = Form("Refs_pCor4_harm%d_gap(%s,%s)_3sub_two_%c",task->fHarmonics, task->GetEtaGapString().Data(), task->GetEtaGapString().Data(), sides[twoPos]);
      listCorFour3sub[twoPos] = new TList();
    }
  }

  // for cumulants (applicable for diff. flow)
  TList* listCumTwo = new TList(); TString nameCumTwo = Form("Refs_hCum2_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());
  TList* listCumFour = new TList(); TString nameCumFour = Form("Refs_hCum4_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());

  // 3 subevents
  TList* listCumTwo3sub[3][3] = {nullptr};
  TString nameCumTwo3sub[3][3] = {""};
  TList* listCumFour3sub[3] = {nullptr};
  TString nameCumFour3sub[3] = {""};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        nameCumTwo3sub[rf1Pos][rf2Pos] = Form("Refs_hCum2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]);
        listCumTwo3sub[rf1Pos][rf2Pos] = new TList();
      }
    }
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      nameCumFour3sub[twoPos] = Form("Refs_hCum4_harm%d_gap(%s,%s)_3sub_two_%c",task->fHarmonics, task->GetEtaGapString().Data(), task->GetEtaGapString().Data(), sides[twoPos]);
      listCumFour3sub[twoPos] = new TList();
    }
  }

  TH1D* desAllCombi[10] = {nullptr};
  if(nOfSamples > 10 && bDoFour && bDo3sub) { Error(Form("Number of samples: %d is more than 10! Mixing of combinations when working with v24 3 sub implemented just for 10 samples! \n !!Change here!!",nOfSamples),"ProcessRefs"); return kFALSE; }

  // for vns desampling
  TList* listFlowTwo = new TList(); TString nameFlowTwo = Form("Refs_hFlow2_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());
  TList* listFlowFour = new TList(); TString nameFlowFour = Form("Refs_hFlow4_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());

  // 3 subevents
  TList* listFlowTwo3sub[3][3] = {nullptr};
  TString nameFlowTwo3sub[3][3] = {""};
  TList* listFlowFour3sub[3] = {nullptr};
  TString nameFlowFour3sub[3] = {""};
  TList* listMergedAllCombinations = new TList();
  TString nameFlowAllCombi = Form("Refs_hFlow4_harm%d_gap(%s,%s)_3sub",task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data());
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        nameFlowTwo3sub[rf1Pos][rf2Pos] = Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]);
        listFlowTwo3sub[rf1Pos][rf2Pos] = new TList();
      }
    }
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      nameFlowFour3sub[twoPos] =  Form("Refs_hFlow4_harm%d_gap(%s,%s)_3sub_two_%c",task->fHarmonics, task->GetEtaGapString().Data(), task->GetEtaGapString().Data(), sides[twoPos]);
      listFlowFour3sub[twoPos] = new TList();
    }
  }

  // estimating <multiplicity>
  if(fbSaveMult)
  {
    TProfile* profMult = (TProfile*) flQACharged->FindObject(Form("fpRefsMult"));
    if(!profMult) { Error("MeanMult profile not found!"); flFlow[kRefs]->ls(); return kFALSE; }
    TProfile* profMult_rebin = (TProfile*) profMult->Rebin(fiNumMultBins,Form("%s_rebin",profMult->GetName()),fdMultBins.data());

    ffOutputFile->cd();
    profMult_rebin->Write(profMult_rebin->GetName());
  }

  // new naming convention for input histos (from FlowTask)
  TString sProfTwoName = Form("<<2>>(%d,-%d)",task->fHarmonics, task->fHarmonics);
  TString sProfFourName = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarmonics, task->fHarmonics, task->fHarmonics, task->fHarmonics);
  if(task->HasGap()) {
    if(!task->Has3sub()){
      sProfTwoName += Form("_2sub(%.2g)",task->fEtaGap);
      sProfFourName += Form("_2sub(%.2g)",task->fEtaGap);
    }
    else {
      sProfTwoName += Form("_3sub(%.2g,%.2g)",task->fEtaGap,task->fEtaGapSecond);
      sProfFourName += Form("_3sub(%.2g,%.2g)",task->fEtaGap,task->fEtaGapSecond);
    }
  }

  Debug("Processing samples","ProcessRefs");
  //for(Short_t iSample(0); iSample < task->fNumSamples; ++iSample)
  for(Short_t iSample(0); iSample < nOfSamples; ++iSample)
  {
    TProfile* pCorTwo = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfTwoName.Data(), iSample));
    if(!pCorTwo) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfTwoName.Data(), iSample)),"ProcesRefs"); flFlow[kRefs]->ls(); return kFALSE; }
    // TProfile* pCorTwo = (TProfile*) flFlow[kRefs]->FindObject(Form("fpRefs_%s<2>_harm%d_gap%s_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iSample));
    // if(!pCorTwo) { Warning(Form("Profile 'pCorTwo' (sample %d) not valid",iSample),"ProcesRefs"); flFlow[kRefs]->ls(); return kFALSE; }

    // 3 subevents
    TProfile* pCorTwo3sub[3][3] = {nullptr};
    if(bDo3sub){
      for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
        if(rf1Pos > 1) break;
        for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
          if(rf1Pos >= rf2Pos) continue;
          pCorTwo3sub[rf1Pos][rf2Pos] = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d_rf1_%c_rf2_%c",sProfTwoName.Data(), iSample,sides[rf1Pos],sides[rf2Pos]));
          if(!pCorTwo3sub[rf1Pos][rf2Pos]) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d_rf1_%c_rf2_%c",sProfTwoName.Data(), iSample,sides[rf1Pos],sides[rf2Pos])),"ProcesRefs"); flFlow[kRefs]->ls(); return kFALSE; }
        }
      }
    }

    // Process 4-particle correlations
    TProfile* pCorFour = nullptr;
    TProfile* pCorFour3sub[3] = {nullptr};
    if(bDoFour)
    {
      pCorFour = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d",sProfFourName.Data(),iSample));
      if(!pCorFour) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName.Data(),iSample)),"ProcesRefs"); flFlow[kRefs]->ls(); return kFALSE; }
      // pCorFour = (TProfile*) flFlow[kRefs]->FindObject(Form("fpRefs_%s<4>_harm%d_gap%s_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iSample));
      // if(!pCorFour) { Warning(Form("Profile 'pCorFour' (sample %d) not valid!",iSample),"ProcesRefs"); flFlow[kRefs]->ls(); return kFALSE; }

      //3 subevents
      if(bDo3sub){
        for(Int_t twoPos(0); twoPos < 3; twoPos++){
          pCorFour3sub[twoPos] = (TProfile*) flFlow[kRefs]->FindObject(Form("%s_Pos_sample%d_two_%c",sProfFourName.Data(),iSample,sides[twoPos]));
          if(!pCorFour3sub[twoPos]) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d_two_%c",sProfFourName.Data(),iSample,sides[twoPos])),"ProcesRefs"); flFlow[kRefs]->ls(); return kFALSE; }
        }
      }
    }

    // rebinning the profiles
    if(task->fRebinning)
    {
      pCorTwo = (TProfile*) pCorTwo->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo.Data(), iSample),fdMultBins.data());
      if(bDoFour) { pCorFour = (TProfile*) pCorFour->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour.Data(), iSample),fdMultBins.data()); }
      if(task->Has3sub()){
        for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
          if(rf1Pos > 1) break;
          for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
            if(rf1Pos >= rf2Pos) continue;
            pCorTwo3sub[rf1Pos][rf2Pos] = (TProfile*) pCorTwo3sub[rf1Pos][rf2Pos]->Rebin(fiNumMultBins,Form("%s_Pos_sample%d_rebin",nameCorTwo3sub[rf1Pos][rf2Pos].Data(), iSample),fdMultBins.data());
          }
        }
        if(bDoFour){
          for(Int_t twoPos(0); twoPos < 3; twoPos++){
            pCorFour3sub[twoPos] = (TProfile*) pCorFour3sub[twoPos]->Rebin(fiNumMultBins,Form("%s_Pos_sample%d_two_%c",nameCorFour3sub[twoPos].Data(),iSample,sides[twoPos]),fdMultBins.data());
          }
        }
      }
    }

    // naming <<X>>
    TString sGap = TString();
    if(task->HasGap() && !bDo3sub) { sGap.Append(Form("{|#Delta#eta| > %g}",task->fEtaGap)); }
    pCorTwo->SetTitle(Form("%s: <<2>>_{%d} %s",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
    listCorTwo->Add(pCorTwo);


    if(bDoFour)
    {
      pCorFour->SetTitle(Form("%s: <<4>>_{%d} %s",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
      listCorFour->Add(pCorFour);
    }

    //3 subevents
    if(bDo3sub){
      sGap.Append(Form(" 3 subevents (#eta mid: -%g,%g)",task->fEtaGap/2,task->fEtaGapSecond/2));
      for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
        if(rf1Pos > 1) break;
        for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
          if(rf1Pos >= rf2Pos) continue;
          pCorTwo3sub[rf1Pos][rf2Pos]->SetTitle(Form("%s: <<2>>_{%d} %s (%c,%c)",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data(),sides[rf1Pos],sides[rf2Pos]));
          listCorTwo3sub[rf1Pos][rf2Pos]->Add(pCorTwo3sub[rf1Pos][rf2Pos]);
        }
      }
      if(bDoFour){
        for(Int_t twoPos(0); twoPos < 3; twoPos++){
          pCorFour3sub[twoPos]->SetTitle(Form("%s: <<4>>_{%d} %s (two %c)",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data(),sides[twoPos]));
          listCorFour3sub[twoPos]->Add(pCorFour3sub[twoPos]);
        }
      }
    }

    // Making cumulants out of correlations : <<N>>_n -> c_n{N} -> v_n{N}
    // cn{2}
    TH1D* hCumTwo = CalcRefCumTwo(pCorTwo,task);
    if(!hCumTwo) { Error(Form("cn{2} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
    hCumTwo->SetName(Form("%s_sample%d", nameCumTwo.Data(), iSample));
    listCumTwo->Add(hCumTwo);

    //3 subevents
    TH1D* hCumTwo3sub[3][3] = {nullptr};
    if(bDo3sub){
      for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
        if(rf1Pos > 1) break;
        for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
          if(rf1Pos >= rf2Pos) continue;
          hCumTwo3sub[rf1Pos][rf2Pos] = CalcRefCumTwo(pCorTwo3sub[rf1Pos][rf2Pos],task,rf1Pos,rf2Pos);
          if(!hCumTwo3sub[rf1Pos][rf2Pos]) { Error(Form("cn{2} rf1_%c_rf2_%c (sample %d) not processed correctly!",sides[rf1Pos],sides[rf2Pos],iSample),"ProcessRefs"); return kFALSE; }
          hCumTwo3sub[rf1Pos][rf2Pos]->SetName(Form("%s_sample%d", nameCumTwo3sub[rf1Pos][rf2Pos].Data(), iSample));
          listCumTwo3sub[rf1Pos][rf2Pos]->Add(hCumTwo3sub[rf1Pos][rf2Pos]);
        }
      }
    }

    // vn{2}
    TH1D* hFlowTwo = CalcRefFlowTwo(hCumTwo,task);
    if(!hFlowTwo) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
    hFlowTwo->SetName(Form("%s_sample%d", nameFlowTwo.Data(), iSample));
    listFlowTwo->Add(hFlowTwo);

    //3 subevents
    TH1D* hFlowTwo3sub[3][3] = {nullptr};
    if(bDo3sub){
      for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
        if(rf1Pos > 1) break;
        for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
          if(rf1Pos >= rf2Pos) continue;
          hFlowTwo3sub[rf1Pos][rf2Pos] = CalcRefFlowTwo(hCumTwo3sub[rf1Pos][rf2Pos],task,rf1Pos,rf2Pos);
          if(!hFlowTwo3sub[rf1Pos][rf2Pos]) { Error(Form("vn{2} rf1_%c_rf2_%c (sample %d) not processed correctly!",sides[rf1Pos],sides[rf2Pos],iSample),"ProcessRefs"); return kFALSE; }
          hFlowTwo3sub[rf1Pos][rf2Pos]->SetName(Form("%s_sample%d", nameFlowTwo3sub[rf1Pos][rf2Pos].Data(), iSample));
          listFlowTwo3sub[rf1Pos][rf2Pos]->Add(hFlowTwo3sub[rf1Pos][rf2Pos]);
        }
      }
    }

    if(task->fCumOrderMax >= 4)
    {
      // cn{4}
      TH1D* hCumFour = CalcRefCumFour(pCorFour, pCorTwo, task, bCorrelated);
      if(!hCumFour) { Error(Form("cn{4} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
      hCumFour->SetName(Form("%s_sample%d", nameCumFour.Data(), iSample));
      listCumFour->Add(hCumFour);

      // vn{4}
      TH1D* hFlowFour = CalcRefFlowFour(hCumFour, task);
      if(!hFlowFour) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
      hFlowFour->SetName(Form("%s_sample%d", nameFlowFour.Data(), iSample));
      listFlowFour->Add(hFlowFour);

      TH1D* hCumFour3sub[3] = {nullptr};
      TH1D* hFlowFour3sub[3] = {nullptr};
      TH1D* mergedAllCombi = nullptr;
      TList* listAllCombi = new TList();

      if(bDo3sub){
        for(Int_t twoPos(0); twoPos < 3; twoPos++){
          //find correct histograms for subtraction
          TProfile* profTwoSubtracting[2] = {nullptr};
          for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
            if(rf1Pos > 1) break;
            for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
              if(rf1Pos >= rf2Pos) continue;
              if(rf1Pos == twoPos || rf2Pos == twoPos){
                if(!profTwoSubtracting[0]) profTwoSubtracting[0] = pCorTwo3sub[rf1Pos][rf2Pos];
                else profTwoSubtracting[1] = pCorTwo3sub[rf1Pos][rf2Pos];

              }
            }
          }
          hCumFour3sub[twoPos] = CalcRefCumFour3sub(pCorFour3sub[twoPos], profTwoSubtracting[0], profTwoSubtracting[1], task, twoPos);
          if(!hCumFour3sub[twoPos]) { Error(Form("cn{4} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
          hCumFour3sub[twoPos]->SetName(Form("%s_sample%d", nameCumFour3sub[twoPos].Data(), iSample));
          listCumFour3sub[twoPos]->Add(hCumFour3sub[twoPos]);

          hFlowFour3sub[twoPos] = CalcRefFlowFour(hCumFour3sub[twoPos], task);
          if(!hFlowFour3sub[twoPos]) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
          hFlowFour3sub[twoPos] ->SetName(Form("%s_sample%d", nameFlowFour3sub[twoPos].Data(), iSample));
          listFlowFour3sub[twoPos]->Add(hFlowFour3sub[twoPos]);

          //making "merge" of all combinations
          listAllCombi->Add(hFlowFour3sub[twoPos]);
          if(!mergedAllCombi) {mergedAllCombi = (TH1D*) hFlowFour3sub[twoPos]->Clone(Form("%s_allC_sample%d",nameFlowAllCombi.Data(),iSample)); mergedAllCombi->Reset(); }
          mergedAllCombi->Add(hFlowFour3sub[twoPos]);
        } // end twoPos
        if(!mergedAllCombi) { Error("Merging of 'hFlowFourDif3sub' failed!","ProcessDirect"); return kFALSE; }

        mergedAllCombi->Scale((Double_t) 1./3);
        mergedAllCombi->SetName(Form("%s_merged_sample%d",nameFlowAllCombi.Data(),iSample));

        TString desName = Form("%s_sample%d", nameFlowAllCombi.Data(),iSample);
        desAllCombi[iSample] = DesampleList(listAllCombi, mergedAllCombi, task, desName);
        if(!desAllCombi[iSample]) { Error("Desampling vn{2} failed","ProcessDirect"); return kFALSE; }
        desAllCombi[iSample]->SetName(Form("%s_des_sample%d",nameFlowAllCombi.Data(),iSample));

        listMergedAllCombinations->Add(desAllCombi[iSample]);
        delete listAllCombi;
      }
    }
  } // end-for {iSample}: samples
  Debug("Samples processing done!","ProcessRefs");

  // merging correlation profiles to get central values
  Debug("Merging correlations for central values", "ProcessRefs");
  TProfile* pCorTwoMerged = (TProfile*) MergeListProfiles(listCorTwo);
  if(!pCorTwoMerged) { Error("Merging of 'pCorTwoMerged' failed!","ProcessRefs"); return kFALSE; }
  pCorTwoMerged->SetName(Form("%s_merged", nameCorTwo.Data()));

  //3 subevents
  TProfile* pCorTwoMerged3sub[3][3] = {nullptr};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        pCorTwoMerged3sub[rf1Pos][rf2Pos] = (TProfile*) MergeListProfiles(listCorTwo3sub[rf1Pos][rf2Pos]);
        if(!pCorTwoMerged3sub[rf1Pos][rf2Pos]) { Error("Merging of 'pCorTwoMerged' failed!","ProcessRefs"); return kFALSE; }
        pCorTwoMerged3sub[rf1Pos][rf2Pos]->SetName( Form("%s_merged", nameCorTwo3sub[rf1Pos][rf2Pos].Data()) );
      }
    }
  }

  TH1D* hCumTwoMerged = CalcRefCumTwo(pCorTwoMerged, task);
  if(!hCumTwoMerged) { Error(Form("cn{2} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
  hCumTwoMerged->SetName(Form("%s_merged", nameCumTwo.Data()));

  //3 subevents
  TH1D* hCumTwoMerged3sub[3][3] = {nullptr};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        hCumTwoMerged3sub[rf1Pos][rf2Pos] = CalcRefCumTwo(pCorTwoMerged3sub[rf1Pos][rf2Pos],task,rf1Pos,rf2Pos);
        if(!hCumTwoMerged3sub[rf1Pos][rf2Pos]) { Error(Form("cn{2} rf1_%c_rf2_%c (merged) not processed correctly!", sides[rf1Pos],sides[rf2Pos]),"ProcessRefs"); return kFALSE; }
        hCumTwoMerged3sub[rf1Pos][rf2Pos]->SetName(Form("%s_merged", nameCumTwo3sub[rf1Pos][rf2Pos].Data()));
      }
    }
  }

  TH1D* hFlowTwoMerged = CalcRefFlowTwo(hCumTwoMerged, task);
  if(!hFlowTwoMerged) { Error(Form("vn{2} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
  hFlowTwoMerged->SetName(Form("%s_merged", nameFlowTwo.Data()));

  //3 subevents
  TH1D* hFlowTwoMerged3sub[3][3] = {nullptr};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        hFlowTwoMerged3sub[rf1Pos][rf2Pos] = CalcRefFlowTwo(hCumTwoMerged3sub[rf1Pos][rf2Pos],task,rf1Pos,rf2Pos);
        if(!hFlowTwoMerged3sub[rf1Pos][rf2Pos]) { Error(Form("vn{2} rf1_%c_rf2_%c (merged) not processed correctly!", sides[rf1Pos],sides[rf2Pos]),"ProcessRefs"); return kFALSE; }
        hFlowTwoMerged3sub[rf1Pos][rf2Pos]->SetName(Form("%s_merged_rf1_%c_rf2_%c", nameFlowTwo3sub[rf1Pos][rf2Pos].Data(),sides[rf1Pos],sides[rf2Pos]));
      }
    }
  }

  TProfile* pCorFourMerged = nullptr;
  TH1D* hCumFourMerged = nullptr;
  TH1D* hFlowFourMerged = nullptr;
  if(bDoFour)
  {
    pCorFourMerged = (TProfile*) MergeListProfiles(listCorFour);
    if(!pCorFourMerged) { Error("Merging of 'pCorFourMerged' failed!","ProcessRefs"); return kFALSE; }
    pCorFourMerged->SetName(Form("%s_merged", nameCorFour.Data()));

    hCumFourMerged = CalcRefCumFour(pCorFourMerged, pCorTwoMerged, task, bCorrelated);
    if(!hCumFourMerged) { Error(Form("cn{4} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
    hCumFourMerged->SetName(Form("%s_merged", nameCumFour.Data()));

    hFlowFourMerged = CalcRefFlowFour(hCumFourMerged, task);
    if(!hFlowFourMerged) { Error(Form("vn{4} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
    hFlowFourMerged->SetName(Form("%s_merged", nameFlowFour.Data()));
  }

  TProfile* pCorFour3subMerged[3] = {nullptr};
  TH1D* hCumFour3subMerged[3] = {nullptr};
  TH1D* hFlowFour3subMerged[3] = {nullptr};
  if(bDoFour && bDo3sub){
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      pCorFour3subMerged[twoPos] = (TProfile*) MergeListProfiles(listCorFour3sub[twoPos]);
      if(!pCorFour3subMerged[twoPos]) { Error("Merging of 'pCorFourMerged' failed!","ProcessRefs"); return kFALSE; }
      pCorFour3subMerged[twoPos]->SetName(Form("%s_merged", nameCorFour3sub[twoPos].Data()));

      //find correct histograms for subtraction
      TProfile* profTwoSubtracting[2] = {nullptr};
      for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
        if(rf1Pos > 1) break;
        for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
          if(rf1Pos >= rf2Pos) continue;
          if(rf1Pos == twoPos || rf2Pos == twoPos){
            if(!profTwoSubtracting[0]) profTwoSubtracting[0] = pCorTwoMerged3sub[rf1Pos][rf2Pos];
            else profTwoSubtracting[1] = pCorTwoMerged3sub[rf1Pos][rf2Pos];
          }
        }
      }
      hCumFour3subMerged[twoPos] = CalcRefCumFour3sub(pCorFour3subMerged[twoPos], profTwoSubtracting[0], profTwoSubtracting[1], task, twoPos);
      if(!hCumFour3subMerged[twoPos]) { Error(Form("cn{4} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
      hCumFour3subMerged[twoPos]->SetName(Form("%s_merged", nameCumFour3sub[twoPos].Data()));

      hFlowFour3subMerged[twoPos] = CalcRefFlowFour(hCumFour3subMerged[twoPos], task);
      if(!hFlowFour3subMerged[twoPos]) { Error(Form("vn{4} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
      hFlowFour3subMerged[twoPos]->SetName(Form("%s_merged", nameFlowFour3sub[twoPos].Data()));
    }
  }

  // desampling
  Debug("Desampling","ProcessRefs");

  TH1D* hCorTwoDesampled = DesampleList(listCorTwo, pCorTwoMerged->ProjectionX(), task, nameCorTwo, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
  if(!hCorTwoDesampled) { Error("Desampling 'hCorTwoDesampled' failed","ProcessRefs"); return kFALSE; }
  hCorTwoDesampled->SetName(nameCorTwo.Data());

  TH1D* hCumTwoDesampled = DesampleList(listCumTwo, hCumTwoMerged, task, nameCumTwo, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
  if(!hCumTwoDesampled) { Error("Desampling 'hCumTwoDesampled' failed","ProcessRefs"); return kFALSE; }
  hCumTwoDesampled->SetName(nameCumTwo.Data());

  TH1D* hFlowTwoDesampled = DesampleList(listFlowTwo, hFlowTwoMerged, task, nameFlowTwo, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
  if(!hFlowTwoDesampled) { Error("Desampling 'hFlowTwoDesampled' failed","ProcessRefs"); return kFALSE; }
  hFlowTwoDesampled->SetName(nameFlowTwo.Data());

  //3 subevents
  TH1D* hCorTwoDesampled3sub[3][3] = {nullptr};
  TH1D* hCumTwoDesampled3sub[3][3] = {nullptr};
  TH1D* hFlowTwoDesampled3sub[3][3] = {nullptr};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        TString name3sub = Form("%s",nameCorTwo3sub[rf1Pos][rf2Pos].Data());
        hCorTwoDesampled3sub[rf1Pos][rf2Pos] = DesampleList(listCorTwo3sub[rf1Pos][rf2Pos], pCorTwoMerged3sub[rf1Pos][rf2Pos]->ProjectionX(), task, name3sub);
        if(!hCorTwoDesampled3sub[rf1Pos][rf2Pos]) { Error("Desampling 'hCorTwoDesampled' failed","ProcessRefs"); return kFALSE; }
        hCorTwoDesampled3sub[rf1Pos][rf2Pos]->SetName(name3sub.Data());

        hCumTwoDesampled3sub[rf1Pos][rf2Pos] = DesampleList(listCumTwo3sub[rf1Pos][rf2Pos], hCumTwoMerged3sub[rf1Pos][rf2Pos], task, nameCumTwo3sub[rf1Pos][rf2Pos]);
        if(!hCumTwoDesampled3sub[rf1Pos][rf2Pos]) { Error("Desampling 'hCumTwoDesampled' failed","ProcessRefs"); return kFALSE; }
        hCumTwoDesampled3sub[rf1Pos][rf2Pos]->SetName(nameCumTwo.Data());

        hFlowTwoDesampled3sub[rf1Pos][rf2Pos] = DesampleList(listFlowTwo3sub[rf1Pos][rf2Pos], hFlowTwoMerged3sub[rf1Pos][rf2Pos], task, nameFlowTwo3sub[rf1Pos][rf2Pos]);
        if(!hFlowTwoDesampled3sub[rf1Pos][rf2Pos]) { Error("Desampling 'hFlowTwoDesampled' failed","ProcessRefs"); return kFALSE; }
        hFlowTwoDesampled3sub[rf1Pos][rf2Pos]->SetName(nameFlowTwo3sub[rf1Pos][rf2Pos].Data());
      }
    }
  }

  ffOutputFile->cd();
  hCorTwoDesampled->Write();
  hCumTwoDesampled->Write();
  hFlowTwoDesampled->Write();
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      if(rf1Pos > 1) break;
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos >= rf2Pos) continue;
        hCorTwoDesampled3sub[rf1Pos][rf2Pos]->Write();
        hCumTwoDesampled3sub[rf1Pos][rf2Pos]->Write();
        hFlowTwoDesampled3sub[rf1Pos][rf2Pos]->Write();
      }
    }
  }

  if(bDoFour)
  {
    TH1D* hCorFourDesampled = DesampleList(listCorFour, pCorFourMerged->ProjectionX(), task, nameCorFour); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
    if(!hCorFourDesampled) { Error("Desampling 'hCorFourDesampled' failed","ProcessRefs"); return kFALSE; }
    hCorFourDesampled->SetName(nameCorFour.Data());

    TH1D* hCumFourDesampled = DesampleList(listCumFour, hCumFourMerged, task, nameCumFour);
    if(!hCumFourDesampled) { Error("Desampling 'hCumFourDesampled' failed","ProcessRefs"); return kFALSE; }
    hCumFourDesampled->SetName(nameCumFour.Data());

    TH1D* hFlowFourDesampled = DesampleList(listFlowFour, hFlowFourMerged, task, nameFlowFour);
    if(!hFlowFourDesampled) { Error("Desampling 'hFlowFourDesampled' failed","ProcessRefs"); return kFALSE; }
    hFlowFourDesampled->SetName(nameFlowFour.Data());

    ffOutputFile->cd();
    hCorFourDesampled->Write();
    hCumFourDesampled->Write();
    hFlowFourDesampled->Write();
  }

  if(bDoFour && bDo3sub){
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      TH1D* hCorFourDesampled3sub = DesampleList(listCorFour3sub[twoPos], pCorFour3subMerged[twoPos]->ProjectionX(), task, nameCorFour3sub[twoPos]);
      if(!hCorFourDesampled3sub) { Error("Desampling 'hCorFourDesampled3sub' failed","ProcessRefs"); return kFALSE; }
      hCorFourDesampled3sub->SetName(nameCorFour3sub[twoPos].Data());

      TH1D* hCumFourDesampled3sub = DesampleList(listCumFour3sub[twoPos], hCumFour3subMerged[twoPos], task, nameCumFour3sub[twoPos]);
      if(!hCumFourDesampled3sub) { Error("Desampling 'hCumFourDesampled3sub' failed","ProcessRefs"); return kFALSE; }
      hCumFourDesampled3sub->SetName(nameCumFour3sub[twoPos].Data());

      TH1D* hFlowFourDesampled3sub = DesampleList(listFlowFour3sub[twoPos], hFlowFour3subMerged[twoPos], task, nameFlowFour3sub[twoPos]);
      if(!hFlowFourDesampled3sub) { Error("Desampling 'hFlowFourDesampled3sub' failed","ProcessRefs"); return kFALSE; }
      hFlowFourDesampled3sub->SetName(nameFlowFour3sub[twoPos].Data());

      TH1D* mergedAllCombiAllSamples = nullptr;
      for(Int_t iSample(0); iSample < nOfSamples; iSample++){
        if(!mergedAllCombiAllSamples) {mergedAllCombiAllSamples = (TH1D*) desAllCombi[iSample]->Clone(Form("%s_allC",nameFlowAllCombi.Data())); mergedAllCombiAllSamples->Reset(); }
        mergedAllCombiAllSamples->Add(desAllCombi[iSample]);
      }
      mergedAllCombiAllSamples->Scale((Double_t) 1./nOfSamples);
      mergedAllCombiAllSamples->SetName(Form("%s_merged",nameFlowAllCombi.Data()));

      TH1D* hDesampledAllC = DesampleList(listMergedAllCombinations, mergedAllCombiAllSamples, task, nameFlowAllCombi);
      if(!hDesampledAllC) { Error("Desampling vn{4} failed","ProcessDirect"); return kFALSE; }
      hDesampledAllC->SetName(nameFlowAllCombi.Data());

      ffOutputFile->cd();
      hCorFourDesampled3sub->Write();
      hCumFourDesampled3sub->Write();
      hFlowFourDesampled3sub->Write();
      hDesampledAllC->Write();
    }
  }

  // Comment :: Not sure what this is about =====>
  // if(!task->fRebinning)
  // {
  //   // no rebinning
  //   TH1D* hNoRebin_rebinned = TestRebin(hDesampledFlow,task);
  //   hNoRebin_rebinned->Write(Form("%s",hNoRebin_rebinned->GetName()));
  //
  //   if(task->fSampleMerging)
  //   {
  //     TH1D* hMerged_rebinned = TestRebin(hMerged,task);
  //     hMerged_rebinned->Write(Form("%s",hMerged_rebinned->GetName()));
  //   }
  // }
  // <=====

  Debug("Processing done","ProcessRefs");

  if(listCorTwo) delete listCorTwo;
  if(listCumTwo) delete listCumTwo;
  if(listFlowTwo) delete listFlowTwo;
  if(listCorFour) delete listCorFour;
  if(listCumFour) delete listCumFour;
  if(listFlowFour) delete listFlowFour;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessSubtraction(FlowTask* task)
{
  // for small systems
  // subtraction of non-flow
  Info("Processing subtraction task","ProcessSubtraction");
  if(!task) { Error("FlowTask not found!","ProcessSubtraction"); return kFALSE; }
  if(fiNumMultBins > 9 ) { Error("Not implemented for more mult bin than 10!","ProcessSubtraction"); return kFALSE; }
  if(task->fCumOrderMax != 2) { Error("Not implemented for differemt cumulant order!","ProcessSubtraction"); return kFALSE; }
  if(!flQACharged) { Error("List 'flQACharged' not found!","ProcessSubtraction"); return kFALSE; }

  Debug("Checks done!","ProcessSubtraction");

  TList* listRefTwo = (TList*) ffDesampleFile->Get(Form("Refs_hCum2_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefTwo) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  Debug("Refs and multiplicity histograms loaded!","ProcessSubtraction");

  TH2D* mult[10] = {nullptr};
  TList* listCumTwo[10] = {nullptr};

  for(Int_t binMult(0); binMult < fiNumMultBins; ++binMult){
    listCumTwo[binMult] = (TList*) ffDesampleFile->Get(Form("%s_hCum2_harm%d_gap%s_cent%d_list", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics,task->GetEtaGapString().Data(),binMult));
    if(!listCumTwo[binMult]) { Error(Form("List 'listCumTwo' bin %d not found!",binMult),"ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
  }
  Debug(Form("Loaded lists : %s _ gap%s for all centralities.",GetSpeciesName(task->fSpecies).Data(),task->GetEtaGapString().Data()),"ProcessSubtraction");
  //different binning for mult histogram (fixed per 10)
  for(Int_t binMult(0); binMult < 10; ++binMult)
  {
    mult[binMult] = (TH2D*) flQACharged->FindObject(Form("fh2MeanMultCharged_Cent%d",binMult));
    if(!mult[binMult]) { Error(Form("Histogram 'MeanMultCharged_Cent%d' not found!",binMult),"ProcessSubtraction"); ffDesampleFile->ls(); return kFALSE; }
  }

  Debug("All set!","ProcessSubtraction");
  Debug("Base: 60-80% based on V0A. Some parts hard coded. Sorry.","ProcessSubtraction");
  Int_t lastBin = 5; //peripheral collisions

  if(fdMultBins[5] != 60 || fdMultBins[6] != 80) { Error("Problem with multiplicity bins!","ProcessSubtraction"); return kFALSE; }

  TH2D* base = (TH2D*) mult[6]->Clone("base");
  base->Add(mult[7]);

  for(Int_t iMultBin(0); iMultBin < lastBin; iMultBin++){
    Debug(Form("Working on mult. bin %d.",iMultBin),"ProcessSubtraction");

    //List for desampling
    TList* listFlowTwo = new TList();
    TString nameFlowTwo = Form("%s_hFlow2_harm%d_gap%s_cent%d_subtracted", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

    TH2D* raw = (TH2D*) mult[iMultBin]->Clone("raw");

    for(Short_t iSample(0); iSample < task->fNumSamples; iSample++){

      Debug(Form("Working on sample %d.",iSample),"ProcessSubtraction");

      TH1D* refCum = (TH1D*) listRefTwo->FindObject(Form("Refs_hCum2_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(), iSample));
      if(!refCum) { Error("Reference cn{2} not loaded!","ProcessSubtraction"); listRefTwo->ls(); return kFALSE; }

      TH1D* hDiffBase = (TH1D*) listCumTwo[lastBin]->FindObject(Form("%s_hCum2_harm%d_gap%s_cent%d_sample%d",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data(), lastBin, iSample));
      if(!hDiffBase) { Error("Base dn{2} not loaded!","ProcessSubtraction"); listRefTwo->ls(); return kFALSE; }

      TH1D* hDiffRaw = (TH1D*) listCumTwo[iMultBin]->FindObject(Form("%s_hCum2_harm%d_gap%s_cent%d_sample%d",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data(), iMultBin, iSample));
      if(!hDiffRaw) { Error(Form("Raw dn{2} (mult %d) not loaded!",iMultBin),"ProcessSubtraction"); listRefTwo->ls(); return kFALSE; }

      TH1D* hFlowTwoDif = CalcSubtracted(task,iMultBin,base,raw,refCum,hDiffBase,hDiffRaw);
      listFlowTwo->Add(hFlowTwoDif);
    }// end for samples

    Debug("Merging correlations for central values", "ProcessSubtraction");

    TH1D* refCumMerged = (TH1D*) ffDesampleFile->Get(Form("Refs_hCum2_harm%d_gap%s_merged",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!refCumMerged) { Error(Form("Reference cn{2} (merged) not loaded!"),"ProcessSubtraction"); return kFALSE; }

    TH1D* hDiffBaseMerged = (TH1D*) ffDesampleFile->Get(Form("%s_hCum2_harm%d_gap%s_cent%d_merged",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data(),lastBin));
    if(!hDiffBaseMerged) { Error(Form("Base dn{2} (merged) not loaded!"),"ProcessSubtraction"); return kFALSE; }

    TH1D* hDiffRawMerged = (TH1D*) ffDesampleFile->Get(Form("%s_hCum2_harm%d_gap%s_cent%d_merged",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
    if(!hDiffRawMerged) { Error(Form("Raw dn{2} (merged) not loaded!"),"ProcessSubtraction"); return kFALSE; }

    TH1D* hFlowTwoDifMerged = CalcSubtracted(task,iMultBin,base,raw,refCumMerged,hDiffBaseMerged,hDiffRawMerged);

    Debug("Desampling","ProcessSubtraction");

    TH1D* hDesampledTwo = DesampleList(listFlowTwo, hFlowTwoDifMerged, task, nameFlowTwo);
    if(!hDesampledTwo) { Error("Desampling vn{2} failed","ProcessSubtraction"); return kFALSE; }
    hDesampledTwo->SetName(nameFlowTwo.Data());

    ffOutputFile->cd();
    hDesampledTwo->Write();

    delete hDesampledTwo;

    delete listFlowTwo;
  } // end multiplicity bins

  delete listRefTwo;
  for(Int_t binMult(0); binMult < 10; ++binMult){
    delete mult[binMult];
    delete listCumTwo[binMult];
  }
  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcSubtracted(FlowTask* task, Int_t iMultBin, TH2D* base, TH2D* raw, TH1D* refCum, TH1D* hDiffBase, TH1D* hDiffRaw)
{
  //calculate subtraction using pt dependent mean values of charged particles

  if(!task) { Error("FlowTask not found!","CalcSubtracted"); return nullptr; }
  if(!base) { Error("Base multiplcity histogram not found!","CalcSubtracted"); return nullptr; }
  if(!raw) { Error("Raw multiplcity histogram not found!","CalcSubtracted"); return nullptr; }
  if(!refCum) { Error("RFP cumulant not found!","CalcSubtracted"); return nullptr; }
  if(!hDiffBase) { Error("POI base cumulant not found!","CalcSubtracted"); return nullptr; }
  if(!hDiffRaw) { Error("POI raw cumulant not found!","CalcSubtracted"); return nullptr; }
  if(hDiffBase->GetNbinsX() != hDiffRaw->GetNbinsX()) { Error("Different number of bins!","CalcSubtracted"); return nullptr; }

  TH1D* hFlowTwoDif = (TH1D*) hDiffBase->Clone("hFlowTwoDif");
  hFlowTwoDif->Reset();

  TH1D* multBaseIntegrated = (TH1D*) base->ProjectionY("multBase",0,-1);
  Double_t meanBaseIntegrated = multBaseIntegrated->GetMean();

  TH1D* multThisBinIntegrated = (TH1D*) raw->ProjectionY("multThisBinIntegrated",0,-1);
  Double_t meanRawIntegrated = multThisBinIntegrated->GetMean();

  Short_t iNumPtBins = task->fNumPtBins;
  Short_t binPtLow = 0;
  Short_t binPtHigh = 0;

  Double_t cumBase = refCum->GetBinContent(6);
  Double_t cumRaw = refCum->GetBinContent(iMultBin + 1);

  for(Short_t binPt(0); binPt < iNumPtBins; binPt++){
    binPtLow = base->GetXaxis()->FindFixBin(task->fPtBinsEdges[binPt]);
    binPtHigh = base->GetXaxis()->FindFixBin(task->fPtBinsEdges[binPt+1]) - 1;

    TH1D* multBase = (TH1D*) base->ProjectionY("multBase",binPtLow,binPtHigh);
    Double_t meanBase = multBase->GetMean();

    TH1D* multThisBin = (TH1D*) raw->ProjectionY("multThisBin",binPtLow,binPtHigh);
    Double_t meanRaw = multThisBin->GetMean();

    Double_t ratio = -9.0;
    if(meanRawIntegrated > 0) ratio = meanBaseIntegrated / meanRawIntegrated;
    else Warning("Problem with the ratio (mean raw integrated < 0)!","CalcSubtracted");
    Double_t denom = cumRaw - ratio * cumBase;
    if(denom > 0) denom = TMath::Sqrt(denom);
    else denom = -9.0;

    Double_t dContBase = hDiffBase->GetBinContent(binPt);
    Double_t dContRaw = hDiffRaw->GetBinContent(binPt);

    if(meanRaw > 0) ratio = (ratio * meanBase) / meanRaw;
    else Warning("Problem with the ratio (mean raw < 0)!","CalcSubtracted");
    if(ratio > 0) ratio = TMath::Sqrt(ratio);
    else Warning("Problem with the ratio (< 0)!","CalcSubtracted");
    Double_t dContOut = dContRaw - ratio * dContRaw;
    if(denom > 0) dContOut = dContOut/denom;
    // printf("ptbin %d: %f - %f , dContBase : %f dContRaw: %f \n", binPt, task->fPtBinsEdges[binPt], task->fPtBinsEdges[binPt+1], dContBase, dContRaw);
    // printf("Multiplicty: base : %f raw: %f \n", meanBase, meanRaw);
    // printf("ratio : %f denom: %f \n\n\n", ratio, denom);
    hFlowTwoDif->SetBinContent(binPt, dContOut);
  }

  return hFlowTwoDif;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefCumTwo(TProfile* hTwoRef, FlowTask* task, Int_t rf1Pos, Int_t rf2Pos)
{
  // Calculate reference c_n{2} out of correlations
  // NOTE: it is just a fancier Clone(): for consistency
  // cn{2} = <<2>>

  if(!hTwoRef) { Error("Profile 'hTwoRef' not valid!","CalcRefCumTwo"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcRefCumTwo"); return nullptr; }

  TH1D* histCum = nullptr;
  if(rf1Pos == rf2Pos)
    histCum = (TH1D*) hTwoRef->ProjectionX(Form("hCum2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
  else
    histCum = (TH1D*) hTwoRef->ProjectionX(Form("hCum2_Refs_harm%d_gap%s_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]));

  TString sGap = TString();
  if(task->HasGap() && !task->Has3sub()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  if(task->Has3sub()) { sGap.Append(Form(" 3 subevents (#eta mid: -%g,%g)",task->fEtaGap/2,task->fEtaGapSecond/2)); }
  histCum->SetTitle(Form("%s: c_{%d}{2%s}",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcSixTwoDif(TProfile* hSix, TProfile* hFour_12, TProfile* hFour_13, TProfile* hTwo_1, TProfile* hTwo_3, TProfile* hTwo_3_gap, FlowTask* task)
{
  if(!task) { Error("FlowTask not found!","CalcSixTwoDif"); return nullptr; }
  if(!hSix) { Error("Profile 'hSix' not valid!","CalcSixTwoDif"); return nullptr; }
  if(!hFour_12) { Error("Profile 'hFour_12' not valid!","CalcSixTwoDif"); return nullptr; }
  if(!hFour_13) { Error("Profile 'hFour_13' not valid!","CalcSixTwoDif"); return nullptr; }
  if(!hTwo_1) { Error("Profile 'hTwo_1' not valid!","CalcSixTwoDif"); return nullptr; }
  if(!hTwo_3 || !hTwo_3_gap) { Error("Profile 'hTwo_3' not valid!","CalcSixTwoDif"); return nullptr; }
  if(hSix->GetNbinsX() != hFour_12->GetNbinsX() || hSix->GetNbinsX() != hFour_13->GetNbinsX()) { Error("Different number of bins! 4-par corr.","CalcSixTwoDif"); return nullptr; }
  if(hSix->GetNbinsX() != hTwo_1->GetNbinsX() || hSix->GetNbinsX() != hTwo_3->GetNbinsX()) { Error("Different number of bins! 2-par corr.","CalcSixTwoDif"); return nullptr; }
  if(hSix->GetNbinsX() != hTwo_3_gap->GetNbinsX()) { Error("Different number of bins! 2-par corr. with gap.","CalcSixTwoDif"); return nullptr; }

  TH1D* hSixCor = (TH1D*) hSix->ProjectionX(Form("hCor6_Refs_harm%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2]));
  hSixCor->SetTitle(Form("%s: FMC6_{%d,%d,%d}",GetSpeciesName(task->fSpecies).Data(), task->fHarm[0],task->fHarm[1],task->fHarm[2]));
  hSixCor->Reset();

  for(Int_t iBin(0); iBin < hSix->GetNbinsX()+2; ++iBin)
  {
    Double_t dContSix = hSix->GetBinContent(iBin);
    Double_t dErrSix = hSix->GetBinError(iBin);

    Double_t dContFour_12 = hFour_12->GetBinContent(iBin);
    Double_t dErrFour_12 = hFour_12->GetBinError(iBin);

    Double_t dContFour_13 = hFour_13->GetBinContent(iBin);
    Double_t dErrFour_13 = hFour_13->GetBinError(iBin);

    Double_t dContTwo_1 = hTwo_1->GetBinContent(iBin);
    Double_t dErrTwo_1 = hTwo_1->GetBinError(iBin);

    Double_t dContTwo_3 = hTwo_3->GetBinContent(iBin);
    Double_t dErrTwo_3 = hTwo_3->GetBinError(iBin);

    Double_t dContTwo_3_gap = hTwo_3_gap->GetBinContent(iBin);
    Double_t dErrTwo_3_gap = hTwo_3_gap->GetBinError(iBin);

    Double_t dContOut = dContSix - 4.0 * dContFour_13 * dContTwo_1 - dContFour_12 * dContTwo_3 + 4.0 * dContTwo_1 * dContTwo_1 * dContTwo_3;
    Double_t norm = dContFour_12 * dContTwo_3_gap;
    if(task->fIsHijing) norm = 1.0;
    hSixCor->SetBinContent(iBin, dContOut/norm);

    Double_t dErrOut = TMath::Power(dErrSix, 2.0);
    Double_t dErrHelp = dErrFour_12 * dContTwo_3;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = 4.0 * dErrFour_13 * dContTwo_1;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrTwo_3 * ( dContFour_12 + 4.0 *  dContTwo_1);
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrTwo_1 * ( 4.0 *  dContFour_13 + 8.0 * dContTwo_1 * dContTwo_3);
    dErrOut += TMath::Power(dErrHelp, 2.0);

    Double_t dErNorm = TMath::Sqrt( TMath::Power( dErrFour_12 * dContTwo_3_gap, 2.0)
                      + TMath::Power(dErrTwo_3_gap * dContFour_12, 2.0) );
    if(task->fIsHijing) dErNorm = 0.0;
    dErrOut = TMath::Power((dErrOut / norm), 2.0) + TMath::Power(dErNorm * dContOut / (norm * norm), 2.0);

    hSixCor->SetBinError(iBin, TMath::Sqrt(dErrOut));
  }
  return hSixCor;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcSixThreeDif(TProfile* hSix, TProfile* hFour_12, TProfile* hFour_13, TProfile* hFour_23, TProfile* hTwo_1, TProfile* hTwo_2, TProfile* hTwo_3, TProfile* hTwo_1_gap, TProfile* hTwo_2_gap, TProfile* hTwo_3_gap, FlowTask* task)
{
  if(!task) { Error("FlowTask not found!","CalcSixThreeDif"); return nullptr; }
  if(!hSix) { Error("Profile 'hSix' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hFour_12) { Error("Profile 'hFour_12' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hFour_13) { Error("Profile 'hFour_13' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hFour_23) { Error("Profile 'hFour_23' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hTwo_1 || !hTwo_1_gap) { Error("Profile 'hTwo_1' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hTwo_2 || !hTwo_2_gap) { Error("Profile 'hTwo_2' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hTwo_3 || !hTwo_3_gap) { Error("Profile 'hTwo_3' not valid!","CalcSixThreeDif"); return nullptr; }
  if(hSix->GetNbinsX() != hFour_12->GetNbinsX() || hSix->GetNbinsX() != hFour_13->GetNbinsX() || hSix->GetNbinsX() != hFour_23->GetNbinsX()) { Error("Different number of bins! 4-par corr.","CalcSixThreeDif"); return nullptr; }
  if(hSix->GetNbinsX() != hTwo_1->GetNbinsX() || hSix->GetNbinsX() != hTwo_2->GetNbinsX() || hSix->GetNbinsX() != hTwo_3->GetNbinsX() ) { Error("Different number of bins! 2-par corr.","CalcSixThreeDif"); return nullptr; }
  if(hSix->GetNbinsX() != hTwo_1_gap->GetNbinsX() || hSix->GetNbinsX() != hTwo_2_gap->GetNbinsX() || hSix->GetNbinsX() != hTwo_3_gap->GetNbinsX() ) { Error("Different number of bins! 2-par corr. with gap.","CalcSixThreeDif"); return nullptr; }

  TH1D* hSixCor = (TH1D*) hSix->ProjectionX(Form("hCor6_Refs_harm%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2]));
  hSixCor->SetTitle(Form("%s: FMC6_{%d,%d,%d}",GetSpeciesName(task->fSpecies).Data(), task->fHarm[0],task->fHarm[1],task->fHarm[2]));
  hSixCor->Reset();

  for(Int_t iBin(0); iBin < hSix->GetNbinsX()+2; ++iBin)
  {
    Double_t dContSix = hSix->GetBinContent(iBin);
    Double_t dErrSix = hSix->GetBinError(iBin);

    Double_t dContFour_12 = hFour_12->GetBinContent(iBin);
    Double_t dErrFour_12 = hFour_12->GetBinError(iBin);

    Double_t dContFour_13 = hFour_13->GetBinContent(iBin);
    Double_t dErrFour_13 = hFour_13->GetBinError(iBin);

    Double_t dContFour_23 = hFour_23->GetBinContent(iBin);
    Double_t dErrFour_23 = hFour_23->GetBinError(iBin);

    Double_t dContTwo_1 = hTwo_1->GetBinContent(iBin);
    Double_t dErrTwo_1 = hTwo_1->GetBinError(iBin);

    Double_t dContTwo_2 = hTwo_2->GetBinContent(iBin);
    Double_t dErrTwo_2 = hTwo_2->GetBinError(iBin);

    Double_t dContTwo_3 = hTwo_3->GetBinContent(iBin);
    Double_t dErrTwo_3 = hTwo_3->GetBinError(iBin);

    Double_t dContTwo_1_gap = hTwo_1_gap->GetBinContent(iBin);
    Double_t dErrTwo_1_gap = hTwo_1_gap->GetBinError(iBin);

    Double_t dContTwo_2_gap = hTwo_2_gap->GetBinContent(iBin);
    Double_t dErrTwo_2_gap = hTwo_2_gap->GetBinError(iBin);

    Double_t dContTwo_3_gap = hTwo_3_gap->GetBinContent(iBin);
    Double_t dErrTwo_3_gap = hTwo_3_gap->GetBinError(iBin);

    Double_t dContOut = dContSix - dContFour_13 * dContTwo_2 - dContFour_12 * dContTwo_3 - dContFour_23 * dContTwo_1 + 2.0 * dContTwo_1 * dContTwo_2 * dContTwo_3;
    Double_t norm = dContTwo_1_gap * dContTwo_2_gap * dContTwo_3_gap;
    if(task->fIsHijing) norm = 1.0;
    hSixCor->SetBinContent(iBin, dContOut/norm);

    Double_t dErrOut = TMath::Power(dErrSix, 2.0);
    Double_t dErrHelp = dErrFour_12 * dContTwo_3;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrFour_13 * dContTwo_2;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrFour_23 * dContTwo_1;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrTwo_3 * ( dContFour_12 + 2.0 *  dContTwo_1 *  dContTwo_2);
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrTwo_2 * ( dContFour_13 + 2.0 *  dContTwo_1 *  dContTwo_3);
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrTwo_1 * ( dContFour_23 + 2.0 *  dContTwo_2 *  dContTwo_3);
    dErrOut += TMath::Power(dErrHelp, 2.0);
    //normalisation
    Double_t dErNorm = TMath::Sqrt( TMath::Power(dErrTwo_1_gap * dContTwo_2_gap * dContTwo_3_gap, 2.0)
                      + TMath::Power(dErrTwo_2_gap * dContTwo_1_gap * dContTwo_3_gap, 2.0)
                      + TMath::Power(dErrTwo_3_gap * dContTwo_1_gap * dContTwo_2_gap, 2.0) );
    if(task->fIsHijing) dErNorm = 0.0;
    dErrOut = TMath::Power((dErrOut / norm), 2.0) + TMath::Power(dErNorm * dContOut / (norm * norm), 2.0);
    hSixCor->SetBinError(iBin, TMath::Sqrt(dErrOut));
  }
  return hSixCor;

}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcFourFMC(TProfile* hFour, TProfile* hTwo_1, TProfile* hTwo_2, TProfile* hTwo_1_gap, TProfile* hTwo_2_gap, FlowTask* task)
{
  if(!task) { Error("FlowTask not found!","CalcSixThreeDif"); return nullptr; }
  if(!hFour) { Error("Profile 'hFour' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hTwo_1 || !hTwo_1_gap) { Error("Profile 'hTwo_1' not valid!","CalcSixThreeDif"); return nullptr; }
  if(!hTwo_2 || !hTwo_2_gap) { Error("Profile 'hTwo_2' not valid!","CalcSixThreeDif"); return nullptr; }
  if(hFour->GetNbinsX() != hTwo_1->GetNbinsX() || hFour->GetNbinsX() != hFour->GetNbinsX()) { Error("Different number of bins! 2-par corr.","CalcSixThreeDif"); return nullptr; }
  if(hFour->GetNbinsX() != hTwo_1_gap->GetNbinsX() || hFour->GetNbinsX() != hTwo_2_gap->GetNbinsX() ) { Error("Different number of bins! 2-par corr. with gap.","CalcSixThreeDif"); return nullptr; }

  TH1D* hFourCor = (TH1D*) hFour->ProjectionX(Form("hCor4_Refs_harm%d%d",task->fHarm[0],task->fHarm[1]));
  hFourCor->SetTitle(Form("%s: FMC4_{%d,%d}",GetSpeciesName(task->fSpecies).Data(), task->fHarm[0],task->fHarm[1]));
  hFourCor->Reset();

  for(Int_t iBin(0); iBin < hFourCor->GetNbinsX()+2; ++iBin)
  {
    Double_t dContFour = hFour->GetBinContent(iBin);
    Double_t dErrFour = hFour->GetBinError(iBin);

    Double_t dContTwo_1 = hTwo_1->GetBinContent(iBin);
    Double_t dErrTwo_1 = hTwo_1->GetBinError(iBin);

    Double_t dContTwo_2 = hTwo_2->GetBinContent(iBin);
    Double_t dErrTwo_2 = hTwo_2->GetBinError(iBin);

    Double_t dContTwo_1_gap = hTwo_1_gap->GetBinContent(iBin);
    Double_t dErrTwo_1_gap = hTwo_1_gap->GetBinError(iBin);

    Double_t dContTwo_2_gap = hTwo_2_gap->GetBinContent(iBin);
    Double_t dErrTwo_2_gap = hTwo_2_gap->GetBinError(iBin);

    Double_t dContOut = dContFour - dContTwo_1 * dContTwo_2;
    Double_t norm = dContTwo_1_gap * dContTwo_2_gap;
    if(task->fIsHijing) norm = 1.0;
    hFourCor->SetBinContent(iBin, dContOut/norm);

    Double_t dErrOut = TMath::Power(dErrFour, 2.0);
    Double_t dErrHelp = dErrTwo_2 * dContTwo_1;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    dErrHelp = dErrTwo_1 * dContTwo_2;
    dErrOut += TMath::Power(dErrHelp, 2.0);
    //normalisation
    Double_t dErNorm = TMath::Sqrt( TMath::Power(dErrTwo_1_gap * dContTwo_2_gap, 2.0)
                      + TMath::Power(dErrTwo_2_gap * dContTwo_1_gap, 2.0) );
    if(task->fIsHijing) dErNorm = 0.0;
    dErrOut = TMath::Power((dErrOut / norm), 2.0) + TMath::Power(dErNorm * dContOut / (norm * norm), 2.0);
    hFourCor->SetBinError(iBin, TMath::Sqrt(dErrOut));
  }
  return hFourCor;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcEight_ThreeOne(TProfile* hEight, TProfile* hSix_123, TProfile* hSix_124, TProfile* hFour_12, TProfile* hFour_14, TProfile* hTwo_1, TProfile* hTwo_4, TProfile* hTwo_4_gap, FlowTask* task)
{
  if(!task) { Error("FlowTask not found!","CalcEight_ThreeOne"); return nullptr; }
  if(!hEight) { Error("Profile 'hEight' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hSix_123) { Error("Profile 'hSix_123' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hSix_124) { Error("Profile 'hSix_124' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hFour_12) { Error("Profile 'hFour_12' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hFour_14) { Error("Profile 'hFour_14' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hTwo_1) { Error("Profile 'hTwo_1' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hTwo_4 || !hTwo_4_gap) { Error("Profile 'hTwo_4' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hSix_123->GetNbinsX() || hEight->GetNbinsX() != hSix_123->GetNbinsX()) { Error("Different number of bins! 6-par corr.","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hFour_12->GetNbinsX() || hEight->GetNbinsX() != hFour_14->GetNbinsX()) { Error("Different number of bins! 4-par corr.","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hTwo_1->GetNbinsX() || hEight->GetNbinsX() != hTwo_4->GetNbinsX()) { Error("Different number of bins! 2-par corr.","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hTwo_4_gap->GetNbinsX()) { Error("Different number of bins! 2-par corr. with gap.","CalcEight_ThreeOne"); return nullptr; }

  TH1D* hEightCor = (TH1D*) hEight->ProjectionX(Form("hCor8_Refs_harm%d%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]));
  hEightCor->SetTitle(Form("%s: FMC8_{%d,%d,%d,%d}",GetSpeciesName(task->fSpecies).Data(), task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]));
  hEightCor->Reset();

  for(Int_t iBin(0); iBin < hEightCor->GetNbinsX()+2; ++iBin)
  {
    Double_t dContEight = hEight->GetBinContent(iBin);
    Double_t dErrEight = hEight->GetBinError(iBin);

    Double_t dContSix_123 = hSix_123->GetBinContent(iBin);
    Double_t dErrSix_123 = hSix_123->GetBinError(iBin);

    Double_t dContSix_124 = hSix_124->GetBinContent(iBin);
    Double_t dErrSix_124 = hSix_124->GetBinError(iBin);

    Double_t dContFour_12 = hFour_12->GetBinContent(iBin);
    Double_t dErrFour_12 = hFour_12->GetBinError(iBin);

    Double_t dContFour_14 = hFour_14->GetBinContent(iBin);
    Double_t dErrFour_14 = hFour_14->GetBinError(iBin);

    Double_t dContTwo_1 = hTwo_1->GetBinContent(iBin);
    Double_t dErrTwo_1 = hTwo_1->GetBinError(iBin);

    Double_t dContTwo_4 = hTwo_4->GetBinContent(iBin);
    Double_t dErrTwo_4 = hTwo_4->GetBinError(iBin);

    Double_t dContTwo_4_gap = hTwo_4_gap->GetBinContent(iBin);
    Double_t dErrTwo_4_gap = hTwo_4_gap->GetBinError(iBin);

    Double_t dContOut = dContEight - 9.0 * dContSix_124 * dContTwo_1 - dContSix_123 * dContTwo_4
                        - 9.0 * dContFour_12 * dContFour_14
                        - 36.0 * TMath::Power(dContTwo_1, 3.0) * dContTwo_4
                        + 18.0 * dContTwo_1 * dContTwo_4 * dContFour_12
                        + 36.0 * dContTwo_1 * dContTwo_1 * dContFour_14;
    Double_t norm = dContSix_123 * dContTwo_4_gap;
    if(task->fIsHijing) norm = 1.0;
    // Double_t dContOut = dContEight;
    if(!(TMath::Abs(norm) > 0.0)) norm = 1.0;
    // Double_t norm = 1.0;
    hEightCor->SetBinContent(iBin, dContOut/norm);

    Double_t dErrOut = TMath::Power(dErrEight, 2.0);
    // need to be implemented
    // hEightCor->SetBinError(iBin, TMath::Sqrt(dErrOut/norm));
    hEightCor->SetBinError(iBin, 1.0);
  }
  return hEightCor;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcEight_TwoTwo(TProfile* hEight, TProfile* hSix_123, TProfile* hSix_134, TProfile* hFour_12, TProfile* hFour_13, TProfile* hFour_34, TProfile* hTwo_1, TProfile* hTwo_4, FlowTask* task)
{
  if(!task) { Error("FlowTask not found!","CalcEight_ThreeOne"); return nullptr; }
  if(!hEight) { Error("Profile 'hEight' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hSix_123) { Error("Profile 'hSix_123' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hSix_134) { Error("Profile 'hSix_134' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hFour_12) { Error("Profile 'hFour_12' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hFour_13) { Error("Profile 'hFour_13' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hFour_34) { Error("Profile 'hFour_34' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hTwo_1) { Error("Profile 'hTwo_1' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(!hTwo_4) { Error("Profile 'hTwo_4' not valid!","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hSix_123->GetNbinsX() || hEight->GetNbinsX() != hSix_134->GetNbinsX()) { Error("Different number of bins! 6-par corr.","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hFour_12->GetNbinsX() || hEight->GetNbinsX() != hFour_13->GetNbinsX() || hEight->GetNbinsX() != hFour_34->GetNbinsX() ) { Error("Different number of bins! 4-par corr.","CalcEight_ThreeOne"); return nullptr; }
  if(hEight->GetNbinsX() != hTwo_1->GetNbinsX() || hEight->GetNbinsX() != hTwo_4->GetNbinsX()) { Error("Different number of bins! 2-par corr.","CalcEight_ThreeOne"); return nullptr; }

  TH1D* hEightCor = (TH1D*) hEight->ProjectionX(Form("hCor8_Refs_harm%d%d%d%d",task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]));
  hEightCor->SetTitle(Form("%s: FMC8_{%d,%d,%d,%d}",GetSpeciesName(task->fSpecies).Data(), task->fHarm[0],task->fHarm[1],task->fHarm[2],task->fHarm[3]));
  hEightCor->Reset();

  for(Int_t iBin(0); iBin < hEightCor->GetNbinsX()+2; ++iBin)
  {
    Double_t dContEight = hEight->GetBinContent(iBin);
    Double_t dErrEight = hEight->GetBinError(iBin);

    Double_t dContSix_123 = hSix_123->GetBinContent(iBin);
    Double_t dErrSix_123 = hSix_123->GetBinError(iBin);

    Double_t dContSix_134 = hSix_134->GetBinContent(iBin);
    Double_t dErrSix_134 = hSix_134->GetBinError(iBin);

    Double_t dContFour_12 = hFour_12->GetBinContent(iBin);
    Double_t dErrFour_12 = hFour_12->GetBinError(iBin);

    Double_t dContFour_13 = hFour_13->GetBinContent(iBin);
    Double_t dErrFour_13 = hFour_13->GetBinError(iBin);

    Double_t dContFour_34 = hFour_34->GetBinContent(iBin);
    Double_t dErrFour_34 = hFour_34->GetBinError(iBin);

    Double_t dContTwo_1 = hTwo_1->GetBinContent(iBin);
    Double_t dErrTwo_1 = hTwo_1->GetBinError(iBin);

    Double_t dContTwo_4 = hTwo_4->GetBinContent(iBin);
    Double_t dErrTwo_4 = hTwo_4->GetBinError(iBin);

    Double_t dContOut = dContEight - 4.0 * dContSix_123 * dContTwo_4
                        - 4.0 * dContSix_134 * dContTwo_1
                        - dContFour_12 * dContFour_34
                        - 8.0 * dContFour_13 * dContFour_13
                        - 24.0 * dContTwo_1 * dContTwo_1 * dContTwo_4 * dContTwo_4
                        + 4.0 * dContTwo_1 * dContTwo_1 * dContFour_34
                        + 4.0 * dContTwo_4 * dContTwo_4 * dContFour_12
                        + 32.0 * dContTwo_1 * dContTwo_4 * dContFour_13;
    Double_t norm = dContFour_12 * dContFour_34;
    if(task->fIsHijing) norm = 1.0;
    if(!(TMath::Abs(norm) > 0.0)) norm = 1.0;
    // Double_t norm = 1;
    hEightCor->SetBinContent(iBin, dContOut/norm);

    // Double_t dErrOut = TMath::Power(dErrEight, 2.0);
    // need to be implemented
    // hEightCor->SetBinError(iBin, TMath::Sqrt(dErrOut/norm));
    hEightCor->SetBinError(iBin, 1.0);
  }
  return hEightCor;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefCumFour(TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task, Bool_t bCorrel)
{
  // Calculate reference c_n{4} out of correlations
  // cn{4} = <<4>> - 2*<<2>>^2

  if(!task) { Error("FlowTask not found!","CalcRefCumFour"); return nullptr; }
  if(!hFourRef) { Error("Profile 'hFourRef' not valid!","CalcRefCumFour"); return nullptr; }
  if(!hTwoRef) { Error("Profile 'hTwoRef' not valid!","CalcRefCumFour"); return nullptr; }
  if(hFourRef->GetNbinsX() != hTwoRef->GetNbinsX()) { Error("Different number of bins!","CalcRefCumFour"); return nullptr; }

  TH1D* histCum = (TH1D*) hFourRef->ProjectionX(Form("hCum4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: c_{%d}{4%s}",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histCum->Reset();

  for(Int_t iBin(0); iBin < hFourRef->GetNbinsX()+2; ++iBin)
  {
    Double_t dContInFour = hFourRef->GetBinContent(iBin);
    Double_t dErrInFour = hFourRef->GetBinError(iBin);

    Double_t dContInTwo = hTwoRef->GetBinContent(iBin);
    Double_t dErrInTwo = hTwoRef->GetBinError(iBin);

    Double_t dContOut = dContInFour - 2.0 * dContInTwo * dContInTwo;
    histCum->SetBinContent(iBin, dContOut);

    Double_t dErrOutFour = dErrInFour; // wrt. <4>
    Double_t dErrOutTwo = -4.0 * dContInTwo * dErrInTwo; // wrt. <2>

    Double_t dErrOut = TMath::Power(dErrOutFour, 2.0) + TMath::Power(dErrOutTwo, 2.0);
    if(bCorrel) { dErrOut += 2.0 * dErrOutFour * dErrOutTwo; }
    histCum->SetBinError(iBin, TMath::Sqrt(dErrOut));
  }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefCumFour3sub(TProfile* hFourRef, TProfile* hTwoRef_sub1, TProfile* hTwoRef_sub2, FlowTask* task, Int_t side)
{
  // Calculate reference c_n{4} out of correlations in 3 subevents case
  // cn{4} = <<4>> - 2*<<2>>*<<2>> (<<2>>s have to be from correct sub-event!)

  if(!task) { Error("FlowTask not found!","CalcRefCumFour3sub"); return nullptr; }
  if(!hFourRef) { Error("Profile 'hFourRef' not valid!","CalcRefCumFour3sub"); return nullptr; }
  if(!hTwoRef_sub1) { Error("Profile 'hTwoRef_sub1' not valid!","CalcRefCumFour3sub"); return nullptr; }
  if(!hTwoRef_sub2) { Error("Profile 'hTwoRef_sub2' not valid!","CalcRefCumFour3sub"); return nullptr; }
  if(hFourRef->GetNbinsX() != hTwoRef_sub1->GetNbinsX() || hFourRef->GetNbinsX() != hTwoRef_sub2->GetNbinsX()) { Error("Different number of bins!","CalcRefCumFour3sub"); return nullptr; }

  TH1D* histCum = (TH1D*) hFourRef->ProjectionX(Form("hCum4_Refs_harm%d_gap%s_two_%c",task->fHarmonics,task->GetEtaGapString().Data(),sides[side]));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",M: |#Delta#eta| < %g",task->fEtaGap/2)); }

  histCum->SetTitle(Form("%s: c_{%d}{4%s} (%c)",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data(), sides[side]));
  histCum->Reset();

  for(Int_t iBin(0); iBin < hFourRef->GetNbinsX()+2; ++iBin)
  {
    Double_t dContInFour = hFourRef->GetBinContent(iBin);
    Double_t dErrInFour = hFourRef->GetBinError(iBin);

    Double_t dContInTwo1 = hTwoRef_sub1->GetBinContent(iBin);
    Double_t dErrInTwo1 = hTwoRef_sub1->GetBinError(iBin);

    Double_t dContInTwo2 = hTwoRef_sub2->GetBinContent(iBin);
    Double_t dErrInTwo3 = hTwoRef_sub2->GetBinError(iBin);

    Double_t dContOut = dContInFour - 2.0 * dContInTwo1 * dContInTwo2;
    histCum->SetBinContent(iBin, dContOut);

    //not changed, errors from desampling!
    Double_t dErrOutFour = dErrInFour; // wrt. <4>
    // Double_t dErrOutTwo = -4.0 * dContInTwo * dErrInTwo; // wrt. <2>

    Double_t dErrOut = TMath::Power(dErrOutFour, 2.0);
    histCum->SetBinError(iBin, TMath::Sqrt(dErrOut));
  }

  return histCum;

}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumTwo(TProfile* hTwoDif, FlowTask* task)
{
  // Same as with TH1D argument.
  // First cast as TH1D and then call CalcDifCumTwo

  TH1D* histCum = CalcDifCumTwo((TH1D*)hTwoDif->ProjectionX("_temp"),task);
  if(!histCum) { Error("Failed!","CalcDifCumTwo"); return nullptr; }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumTwo(TH1D* hTwoDif, FlowTask* task)
{
  // Calculate reference d_n{2} out of correlations
  // NOTE: it is just a fancier Clone(): for consistency
  // dn{2} = <<2'>>

  if(!hTwoDif) { Error("Input 'hTwoDif' not valid!","CalcDifCumTwo"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcDifCumTwo"); return nullptr; }

  TH1D* histCum = (TH1D*) hTwoDif->Clone(Form("hCum2_%s_harm%d_gap%s",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: d_{%d}{2%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumFour(TProfile* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Same as with TH1D argument.
  // First cast as TH1D and then call CalcDifCumFour

  TH1D* histCum = CalcDifCumFour((TH1D*)hFourDif->ProjectionX("_temp"),hTwoDif, hTwoRef, iRefBin, task, bCorrel);
  if(!histCum) { Error("Failed!","CalcDifCumFour"); return nullptr; }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumFour(TH1D* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Calculate reference d_n{4} out of correlations
  // dn{4} = <<4'>> - 2<<2>><<2'>>

  if(!hFourDif) { Error("Input 'hFourDif' not valid!","CalcDifCumFour"); return nullptr; }
  if(!hTwoDif) { Error("Input 'hTwoDif' not valid!","CalcDifCumFour"); return nullptr; }
  if(!hTwoRef) { Error("Input 'hTwoRef' not valid!","CalcDifCumFour"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcDifCumFour"); return nullptr; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!'","CalcDifCumFour"); return nullptr; }
  if(hFourDif->GetNbinsX() != hTwoDif->GetNbinsX()) { Error("Different number of bins!","CalcDifCumFlow"); return nullptr; }

  TH1D* histCum = (TH1D*) hFourDif->Clone(Form("hCum4_%s_harm%d_gap%s",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: d_{%d}{4%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histCum->Reset();

  Double_t dContInTwoRef = hTwoRef->GetBinContent(iRefBin);
  Double_t dErrInTwoRef = hTwoRef->GetBinError(iRefBin);

  for(Int_t iBin(0); iBin < histCum->GetNbinsX()+2; ++iBin)
  {
    Double_t dContInFourDif = hFourDif->GetBinContent(iBin);
    Double_t dErrInFourDif = hFourDif->GetBinError(iBin);

    Double_t dContInTwoDif = hTwoDif->GetBinContent(iBin);
    Double_t dErrInTwoDif = hTwoDif->GetBinError(iBin);

    Double_t dContOut = dContInFourDif - 2.0 * dContInTwoDif * dContInTwoRef;
    histCum->SetBinContent(iBin, dContOut);

    Double_t dErrOutFour = dErrInFourDif; // wrt. <4>
    Double_t dErrOutTwoDif =  -2.0 * dContInTwoRef * dErrInTwoDif; // wrt. <2'>
    Double_t dErrOutTwoRef =  -2.0 * dContInTwoDif * dErrInTwoRef; // wrt. <2>

    Double_t dErrOutSq = TMath::Power(dErrOutFour, 2.0) + TMath::Power(dErrOutTwoDif, 2.0) + TMath::Power(dErrOutTwoRef, 2.0);
    if(bCorrel) { dErrOutSq += 2.0 * dErrOutFour * dErrOutTwoDif + 2.0 * dErrOutFour * dErrOutTwoRef + 2.0 * dErrOutTwoDif * dErrOutTwoRef; }
    histCum->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
  }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumFour3sub(TProfile* hFourDif, TH1* hTwoDif1, TH1* hTwoDif2, TH1* hTwoRef1, TH1* hTwoRef2, Int_t iRefBin, FlowTask* task)
{
  // Same as with TH1D argument.
  // First cast as TH1D and then call CalcDifCumFour

  TH1D* histCum = CalcDifCumFour3sub((TH1D*)hFourDif->ProjectionX("_temp"),hTwoDif1, hTwoDif2, hTwoRef1, hTwoRef2, iRefBin, task);
  if(!histCum) { Error("Failed!","CalcDifCumFour"); return nullptr; }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumFour3sub(TH1D* hFourDif, TH1* hTwoDif1, TH1* hTwoDif2, TH1* hTwoRef1, TH1* hTwoRef2, Int_t iRefBin, FlowTask* task)
{
  // Calculate reference d_n{4} out of correlations
  // dn{4} = <<4'>> - <<2>><<2'>> - <<2>><<2'>>
  // only for case 2-p.c. is within the same sub-event as POI

  if(!hFourDif) { Error("Input 'hFourDif' not valid!","CalcDifCumFour3sub"); return nullptr; }
  if(!hTwoDif1) { Error("Input 'hTwoDif1' not valid!","CalcDifCumFour3sub"); return nullptr; }
  if(!hTwoDif2) { Error("Input 'hTwoDif2' not valid!","CalcDifCumFour3sub"); return nullptr; }
  if(!hTwoRef1) { Error("Input 'hTwoRef1' not valid!","CalcDifCumFour3sub"); return nullptr; }
  if(!hTwoRef2) { Error("Input 'hTwoRef2' not valid!","CalcDifCumFour3sub"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcDifCumFour3sub"); return nullptr; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!'","CalcDifCumFour3sub"); return nullptr; }
  if(hFourDif->GetNbinsX() != hTwoDif1->GetNbinsX() || hFourDif->GetNbinsX() != hTwoDif2->GetNbinsX() ) { Error("Different number of bins!","CalcDifCumFour3sub"); return nullptr; }

  TH1D* histCum = (TH1D*) hFourDif->Clone(Form("hCum4_%s_harm%d_gap%s",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  Debug(Form("\n Names: \n  %s -   %s  *  %s   - %s  *  %s ",hFourDif->GetName(), hTwoDif1->GetName(), hTwoRef1->GetName(), hTwoDif2->GetName(), hTwoRef2->GetName()), "CalcDifCumFour3sub");

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: d_{%d}{4%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histCum->Reset();

  Double_t dContInTwoRef1 = hTwoRef1->GetBinContent(iRefBin);
  Double_t dErrInTwoRef1 = hTwoRef1->GetBinError(iRefBin);

  Double_t dContInTwoRef2 = hTwoRef2->GetBinContent(iRefBin);
  Double_t dErrInTwoRef2 = hTwoRef2->GetBinError(iRefBin);

  for(Int_t iBin(0); iBin < histCum->GetNbinsX()+2; ++iBin)
  {
    Double_t dContInFourDif = hFourDif->GetBinContent(iBin);
    Double_t dErrInFourDif = hFourDif->GetBinError(iBin);

    Double_t dContInTwoDif1 = hTwoDif1->GetBinContent(iBin);
    Double_t dErrInTwoDif1 = hTwoDif1->GetBinError(iBin);

    Double_t dContInTwoDif2 = hTwoDif2->GetBinContent(iBin);
    Double_t dErrInTwoDif2 = hTwoDif2->GetBinError(iBin);

    Double_t dContOut = dContInFourDif - dContInTwoDif1 * dContInTwoRef1 - dContInTwoDif2 * dContInTwoRef2;
    histCum->SetBinContent(iBin, dContOut);

    Double_t dErrOutFour = dErrInFourDif; // wrt. <4>
    Double_t dErrOutTwoDif1 =  dContInTwoRef1 * dErrInTwoDif1; // wrt. <2'>1
    Double_t dErrOutTwoRef1 =  dContInTwoDif1 * dErrInTwoRef1; // wrt. <2>1
    Double_t dErrOutTwoDif2 =  dContInTwoRef2 * dErrInTwoDif2; // wrt. <2'>2
    Double_t dErrOutTwoRef2 =  dContInTwoDif2 * dErrInTwoRef2; // wrt. <2>2

    Double_t dErrOutSq = TMath::Power(dErrOutFour, 2.0) + TMath::Power(dErrOutTwoDif1, 2.0) + TMath::Power(dErrOutTwoRef1, 2.0) + TMath::Power(dErrOutTwoDif2, 2.0) + TMath::Power(dErrOutTwoRef2, 2.0);
    histCum->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
  }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefFlowTwo(TH1D* hTwoRef, FlowTask* task, Int_t rf1Pos, Int_t rf2Pos)
{
  // Calculate reference v_n{2} out of c_n{2}
  // vn{2} = cn{2}^(1/2)

  if(!hTwoRef) { Error("Histo 'hTwoRef' not valid!","CalcRefFlowTwo"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcRefFlowTwo"); return nullptr; }

  TH1D* histFlow = nullptr;
  if(rf1Pos == rf2Pos)
    histFlow = (TH1D*) hTwoRef->Clone(Form("hFlow2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
  else
    histFlow = (TH1D*) hTwoRef->Clone(Form("hFlow2_Refs_harm%d_gap%s_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]));

  TString sGap = TString();
  if(task->HasGap() && !task->Has3sub()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  if(task->Has3sub()) { sGap.Append(Form(" 3 subevents (#eta mid: -%g,%g)",task->fEtaGap/2,task->fEtaGapSecond/2)); }
  histFlow->SetTitle(Form("%s: v_{%d}{2%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  for(Short_t iBin(0); iBin < hTwoRef->GetNbinsX()+2; ++iBin)
  {
    Double_t dContIn = hTwoRef->GetBinContent(iBin);
    Double_t dErrIn = hTwoRef->GetBinError(iBin);

    if(dContIn > 0.0 && dErrIn >= 0.0)
    {
      Double_t dContOut = TMath::Sqrt(dContIn);
      histFlow->SetBinContent(iBin, dContOut);
      Double_t dErrOutSq = 0.25 * dErrIn * dErrIn / dContIn;
      histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
    }
    else
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }
  }

  return histFlow;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefFlowFour(TH1D* hFourRef, FlowTask* task)
{
  // Calculate reference v_n{4} out of c_n{4}
  // vn{4} = (-cn{4})^(1/4)

  if(!hFourRef) { Error("Histo 'hFourRef' not valid!","CalcRefFlowFour"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcRefFlowFour"); return nullptr; }

  TH1D* histFlow = (TH1D*) hFourRef->Clone(Form("hFlow4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{4%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  for(Short_t iBin(0); iBin < hFourRef->GetNbinsX()+2; ++iBin)
  {
    Double_t dContIn = hFourRef->GetBinContent(iBin);
    Double_t dErrIn = hFourRef->GetBinError(iBin);

    if(dContIn < 0.0 && dErrIn >= 0.0)
    {
      Double_t dContOut = TMath::Power(-1.0 * dContIn, 0.25);
      histFlow->SetBinContent(iBin, dContOut);
      Double_t dErrOutSq = TMath::Power(0.25 * dErrIn * TMath::Power(-dContIn, -0.75), 2.0);
      histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
    }
    else
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }
  }

  return histFlow;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifFlowTwo(TH1D* hTwoDif, TH1D* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Calculate differential v_n{2} out of d_n{2} & v_n{2} (!)
  // vn'{2} = dn{2} / vn{2}

  if(!hTwoDif) { Error("Histo 'hTwoDif' not valid!","CalcDifFlowTwo"); return nullptr; }
  if(!hTwoRef) { Error("Histo 'hTwoRef' not valid!","CalcDifFlowTwo"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcDifFlowTwo"); return nullptr; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!'","CalcDifFlowTwo"); return nullptr; }

  TH1D* histFlow = (TH1D*) hTwoDif->Clone(Form("hFlow2_%s_harm%d_gap%s",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{2%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  Double_t dContInRef = hTwoRef->GetBinContent(iRefBin);
  Double_t dErrInRef = hTwoRef->GetBinError(iRefBin);

  // flow not real -> putting 'wrong' numbers in
  if(dContInRef < -9.0 || (dContInRef <= 0.0 && dErrInRef > 1000))
  {
    for(Short_t iBin(0); iBin < histFlow->GetNbinsX()+2; ++iBin)
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }

    return histFlow;
  }

  // flow real -> correct analytical calculation
  for(Short_t iBin(0); iBin < histFlow->GetNbinsX()+2; ++iBin)
  {
    Double_t dContInDif = hTwoDif->GetBinContent(iBin);
    Double_t dErrInDif = hTwoDif->GetBinError(iBin);

    Double_t dContOut = dContInDif / dContInRef;
    histFlow->SetBinContent(iBin, dContOut);

    Double_t dErrOutDif = dErrInDif / dContInRef;
    Double_t dErrOutRef = -1.0 * dContInDif * TMath::Power(dContInRef, -2.0) * dErrInRef;

    Double_t dErrOutSq = TMath::Power(dErrOutDif, 2.0) + TMath::Power(dErrOutRef, 2.0);
    if(bCorrel) { dErrOutSq += 2.0 * dErrOutDif * dErrOutRef; }
    histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
  }

  return histFlow;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifFlowFour(TH1D* hFourDif, TH1D* hFourRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Calculate differential v_n{4} out of d_n{4} & v_n{4} (!)
  // vn'{2} = - dn{4} / vn{4}^3

  if(!hFourDif) { Error("Histo 'hFourDif' not valid!","CalcDifFlowFour"); return nullptr; }
  if(!hFourRef) { Error("Histo 'hFourRef' not valid!","CalcDifFlowFour"); return nullptr; }
  if(!task) { Error("FlowTask not found!","CalcDifFlowFour"); return nullptr; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!","CalcDifFlowFour"); return nullptr; }

  TH1D* histFlow = (TH1D*) hFourDif->Clone(Form("hFlow4_%s_harm%d_gap%s",GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->HasGap()) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{4%s}",GetSpeciesLabel(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  Double_t dContInRef = hFourRef->GetBinContent(iRefBin);
  Double_t dErrInRef = hFourRef->GetBinError(iRefBin);

  // flow not real -> putting 'wrong' numbers in
  if(dContInRef <= 0.0 || (dContInRef < -9.0 && dErrInRef > 1000))
  {
    for(Short_t iBin(0); iBin < histFlow->GetNbinsX()+2; ++iBin)
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }

    return histFlow;
  }

  // flow real -> correct analytical calculation
  for(Short_t iBin(0); iBin < histFlow->GetNbinsX()+2; ++iBin)
  {
    Double_t dContInDif = hFourDif->GetBinContent(iBin);
    Double_t dErrInDif = hFourDif->GetBinError(iBin);

    Double_t dContOut = -1.0 * dContInDif * TMath::Power(dContInRef, -3.0);
    histFlow->SetBinContent(iBin, dContOut);

    Double_t dErrOutDif = -1.0 * dErrInDif * TMath::Power(dContInRef, -3.0); // wrt. dn
    Double_t dErrOutRef = 3.0 * dContInDif * dErrInRef * TMath::Power(dContInRef, -4.0); // wrt. vn

    Double_t dErrOutSq = TMath::Power(dErrOutDif, 2.0) + TMath::Power(dErrOutRef, 2.0);
    if(bCorrel) { dErrOutSq += 2.0 * dErrOutDif * dErrOutRef; }
    histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
  }

  return histFlow;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessDirect(FlowTask* task, Short_t iMultBin)
{
  Info(Form("Processing direct task (mult. bin %d)", iMultBin),"ProcesDirect");
  if(!task) { Error("Task not valid!","ProcessDirect"); return kFALSE; }

  PartSpecies spec = task->fSpecies;
  if(!IsSpeciesDirect(spec)) {
    Error("Task species not direct!","ProcessDirect");
    return kFALSE;
  }

  TList* listInput = flFlow[spec];
  if(!listInput) { Error("Input list not loaded!","ProcessDirect"); return kFALSE; }

  Bool_t bDoFour = task->fCumOrderMax >= 4; // check if cn{4} should be processed
  Bool_t bCorrelated = task->fConsCorr; // check if correlated uncrt. are considered
  Bool_t bDo3sub = task->Has3sub();

  // Loading list where reference flow samples are stored
  TList* listRefCorTwo = (TList*) ffDesampleFile->Get(Form("Refs_pCor2_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefCorTwo) { Error("List 'listRefCorTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefTwo = (TList*) ffDesampleFile->Get(Form("Refs_hFlow2_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefTwo) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  //3 sub-events
  TList* listRefCorTwo3sub[3][3] = {nullptr};
  TList* listRefTwo3sub[3][3] = {nullptr};
  if(bDo3sub){
    for(Int_t rf1Pos(0); rf1Pos < 3; rf1Pos++){
      for(Int_t rf2Pos(0); rf2Pos < 3; rf2Pos++ ){
        if(rf1Pos == rf2Pos) continue;
        if(rf1Pos < rf2Pos) {
          listRefCorTwo3sub[rf1Pos][rf2Pos] = (TList*) ffDesampleFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_list",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]));
          if(!listRefCorTwo3sub[rf1Pos][rf2Pos]) { Error("List 'listRefCorTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

          listRefTwo3sub[rf1Pos][rf2Pos] = (TList*) ffDesampleFile->Get(Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_list",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf1Pos],sides[rf2Pos]));
          if(!listRefTwo3sub[rf1Pos][rf2Pos]) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
        }
        else {
          listRefCorTwo3sub[rf1Pos][rf2Pos] = (TList*) ffDesampleFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_list",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf2Pos],sides[rf1Pos]));
          if(!listRefCorTwo3sub[rf1Pos][rf2Pos]) { Error("List 'listRefCorTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

          listRefTwo3sub[rf1Pos][rf2Pos] = (TList*) ffDesampleFile->Get(Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_list",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rf2Pos],sides[rf1Pos]));
          if(!listRefTwo3sub[rf1Pos][rf2Pos]) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
        }
      }
    }
  }

  TList* listRefFour = nullptr;
  TList* listRefFour3sub[3] = {nullptr};
  if(bDoFour)
  {
    listRefFour = (TList*) ffDesampleFile->Get(Form("Refs_hFlow4_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!listRefFour) { Error("List 'listRefFour' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

    if(bDo3sub){
      for(Int_t twoPos(0); twoPos < 3; twoPos++){
        listRefFour3sub[twoPos] = (TList*) ffDesampleFile->Get(Form("Refs_hFlow4_harm%d_gap(%s,%s)_3sub_two_%c_list",task->fHarmonics,task->GetEtaGapString().Data(), task->GetEtaGapString().Data(), sides[twoPos]));
        if(!listRefFour3sub[twoPos]) { Error(Form("List 'listRefFour3sub[%c]' not found!",sides[twoPos]),"ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
      }
    }
  }

  // List for desampling
  TList* listCorTwo = new TList(); TString nameCorTwo = Form("%s_pCor2_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listCumTwo = new TList(); TString nameCumTwo = Form("%s_hCum2_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listFlowTwo = new TList(); TString nameFlowTwo = Form("%s_hFlow2_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

  //3 sub-events
  TList* listCorTwo3sub[3][3] = {nullptr};
  TList* listCumTwo3sub[3][3] = {nullptr};
  TList* listFlowTwo3sub[3][3] = {nullptr};
  TString nameCorTwo3sub[3][3] = {""};
  TString nameCumTwo3sub[3][3] = {""};
  TString nameFlowTwo3sub[3][3] = {""};
  if(bDo3sub){
    for(Int_t poiPos(0); poiPos < 3; poiPos++){
      for(Int_t rfPos(0); rfPos < 3; rfPos++){
        if(poiPos == rfPos) continue;
        nameCorTwo3sub[poiPos][rfPos] =  Form("%s_pCor2_harm%d_gap(%s,%s)_cent%d_3sub_poi_%c_rf_%c", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),iMultBin,sides[poiPos],sides[rfPos]);
        nameCumTwo3sub[poiPos][rfPos] =  Form("%s_hCum2_harm%d_gap(%s,%s)_cent%d_3sub_poi_%c_rf_%c", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),iMultBin,sides[poiPos],sides[rfPos]);
        nameFlowTwo3sub[poiPos][rfPos] =  Form("%s_hFlow2_harm%d_gap(%s,%s)_cent%d_3sub_poi_%c_rf_%c", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),iMultBin,sides[poiPos],sides[rfPos]);
        listCorTwo3sub[poiPos][rfPos] = new TList();
        listCumTwo3sub[poiPos][rfPos] = new TList();
        listFlowTwo3sub[poiPos][rfPos] = new TList();
      }
    }
  }

  TList* listCorFour = new TList(); TString nameCorFour = Form("%s_pCor4_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listCumFour = new TList(); TString nameCumFour = Form("%s_hCum4_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listFlowFour = new TList(); TString nameFlowFour = Form("%s_hFlow4_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

  //3 sub-events
  TList* listCorFour3sub[3][3] = {nullptr};
  TList* listCumFour3sub[3][3] = {nullptr};
  TList* listFlowFour3sub[3][3] = {nullptr};
  TList* listMergedAllCombinations = new TList();
  TString nameCorFour3sub[3][3] = {""};
  TString nameCumFour3sub[3][3] = {""};
  TString nameFlowFour3sub[3][3] = {""};
  TString nameFlowAllCombi = Form("%s_hFlow4_harm%d_gap(%s,%s)_cent%d_3sub", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(), task->GetEtaGapString().Data(), iMultBin);
  if(bDoFour && bDo3sub){
    for(Int_t poiPos(0); poiPos < 3; poiPos++){
      for(Int_t twoPos(0); twoPos < 3; twoPos++){
        nameCorFour3sub[poiPos][twoPos] =  Form("%s_pCor4_harm%d_gap(%s,%s)_cent%d_3sub_poi_%c_two_%c", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),iMultBin,sides[poiPos],sides[twoPos]);
        nameCumFour3sub[poiPos][twoPos] =  Form("%s_hCum4_harm%d_gap(%s,%s)_cent%d_3sub_poi_%c_two_%c", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),iMultBin,sides[poiPos],sides[twoPos]);
        nameFlowFour3sub[poiPos][twoPos] =  Form("%s_hFlow4_harm%d_gap(%s,%s)_cent%d_3sub_poi_%c_two_%c", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),iMultBin,sides[poiPos],sides[twoPos]);
        listCorFour3sub[poiPos][twoPos] = new TList();
        listCumFour3sub[poiPos][twoPos] = new TList();
        listFlowFour3sub[poiPos][twoPos] = new TList();
      }
    }
  }

  Debug("Processing samples","ProcessDirect");

  // new naming convention for input histos (from FlowTask)
  TString sProfTwoName = Form("<<2>>(%d,-%d)",task->fHarmonics, task->fHarmonics);
  TString sProfFourName = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarmonics, task->fHarmonics, task->fHarmonics, task->fHarmonics);
  if(task->HasGap()){
    if(!bDo3sub){
      sProfTwoName += Form("_2sub(%.2g)",task->fEtaGap);
      sProfFourName += Form("_2sub(%.2g)",task->fEtaGap);
    }
    else {
      sProfTwoName += Form("_3sub(%.2g,%.2g)",task->fEtaGap,task->fEtaGapSecond);
      sProfFourName += Form("_3sub(%.2g,%.2g)",task->fEtaGap,task->fEtaGapSecond);
    }
  }

  Int_t nOfSamples = task->fNumSamples;

  TProfile2D* p2CorTwoDif = nullptr;
  TProfile2D* p2CorFourDif = nullptr;
  TProfile2D* p2CorTwoDif3sub[3][3] = {nullptr};
  TProfile* pCorTwoDif3sub[3][3] = {nullptr};
  TProfile2D* p2CorFourDif3sub[3][3] = {nullptr};
  TProfile* pCorFourDif3sub[3][3] = {nullptr};

  TH1D* desAllCombi[10] = {nullptr};
  if(nOfSamples > 10 && bDoFour && bDo3sub) { Error(Form("Number of samples: %d is more than 10! Mixing of combinations when working with v24 3 sub implemented just for 10 samples! \n !!Change here!!",nOfSamples),"ProcessDirect"); return kFALSE; }

  for(Short_t iSample(0); iSample < nOfSamples; iSample++)
  {
    Debug(Form("Processing sample %d",iSample), "ProcessDirect");
    // <<2'>>

    if(task->fMergePosNeg)
    {
      // loading pos & neg if fMergePosNeg is ON
      TProfile2D* prof2pos = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfTwoName.Data(),iSample));
      TProfile2D* prof2neg = (TProfile2D*) listInput->FindObject(Form("%s_Neg_sample%d",sProfTwoName.Data(),iSample));
      if(!prof2pos || !prof2neg) { Error("<2>: Pos & Neg profile merging: 'prof2pos' OR 'prof2neg' not found!.","ProcessDirect"); return kFALSE; }

      // merging pos & neg
      TList* listMerge = new TList();
      listMerge->Add(prof2pos);
      listMerge->Add(prof2neg);
      p2CorTwoDif = (TProfile2D*) MergeListProfiles(listMerge);
      delete listMerge; // first delete, then check (return)
      if(!p2CorTwoDif) { Error("<2>: Pos & Neg profile merging failed!","ProcessDirect"); return kFALSE; }
    }
    else
    {
      // loading single (Pos) profile
      if(task->fInputTag.EqualTo("")) // loading default-ly named profile
      { p2CorTwoDif = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfTwoName.Data(),iSample)); }
      else // loading "non-standardly" named profile
      { p2CorTwoDif = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%s_%s_sample%d",GetSpeciesName(task->fSpecies).Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),task->fInputTag.Data(),iSample)); }
    }
    if(!p2CorTwoDif) { Error(Form("Profile '%s' (sample %d) does not exists.",sProfTwoName.Data(),iSample),"ProcessDirect"); return kFALSE; }

    if(bDo3sub){
      for(Int_t poiPos(0); poiPos < 3; poiPos++){
        for(Int_t rfPos(0); rfPos < 3; rfPos++){
          if(poiPos == rfPos) continue;
          p2CorTwoDif3sub[poiPos][rfPos] = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d_poi_%c_rfp_%c",sProfTwoName.Data(),iSample,sides[poiPos],sides[rfPos]));
          if(!p2CorTwoDif3sub[poiPos][rfPos]) { Error(Form("Profile '%s' (sample %d) does not exists.",sProfTwoName.Data(),iSample),"ProcessDirect"); return kFALSE; }
        }
      }
    }
    // <<4'>>
    if(bDoFour)
    {
      if(task->fMergePosNeg)
      {
        // loading pos & neg if fMergePosNeg is ON
        TProfile2D* prof2pos = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfFourName.Data(),iSample));
        TProfile2D* prof2neg = (TProfile2D*) listInput->FindObject(Form("%s_Neg_sample%d",sProfFourName.Data(),iSample));
        if(!prof2pos || !prof2neg) { Error("<4>: Pos & Neg profile merging: 'prof2pos' OR 'prof2neg' not found!.","ProcessDirect"); return kFALSE; }

        // merging pos & neg
        TList* listMerge = new TList();
        listMerge->Add(prof2pos);
        listMerge->Add(prof2neg);
        p2CorFourDif = (TProfile2D*) MergeListProfiles(listMerge);
        delete listMerge; // first delete, then check (return)
        if(!p2CorFourDif) { Error("<4>: Pos & Neg profile merging failed!","ProcessDirect"); return kFALSE; }
      }
      else
      {
        // loading single (Pos) profile
        if(task->fInputTag.EqualTo("")) // loading default-ly named profile
        { p2CorFourDif = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfFourName.Data(),iSample)); }
        else // loading "non-standardly" named profile
        { p2CorFourDif = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<4>_harm%d_gap%s_%s_sample%d",GetSpeciesName(task->fSpecies).Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),task->fInputTag.Data(),iSample)); }
      }
      if(!p2CorFourDif) { Error(Form("Profile '%s' (sample %d) does not exists.",sProfFourName.Data(),iSample),"ProcessDirect"); listInput->ls(); return kFALSE; }
    }

    if(bDoFour && bDo3sub){
      for(Int_t poiPos(0); poiPos < 3; poiPos++){
        for(Int_t twoPos(0); twoPos < 3; twoPos++){
          p2CorFourDif3sub[poiPos][twoPos] = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d_poi_%c_two_%c",sProfFourName.Data(),iSample,sides[poiPos],sides[twoPos]));
          if(!p2CorFourDif3sub[poiPos][twoPos]) { Error(Form("Profile '%s' (sample %d, poi %c, two %c) does not exists.",sProfFourName.Data(),iSample,sides[poiPos],sides[twoPos]),"ProcessDirect"); return kFALSE; }
        }
      }
    }

    Debug("Rebinning profiles","ProcessDirect");
    // rebinning according in mult bin
    Short_t binMultLow = p2CorTwoDif->GetXaxis()->FindFixBin(fdMultBins[iMultBin]);
    Short_t binMultHigh = p2CorTwoDif->GetXaxis()->FindFixBin(fdMultBins[iMultBin+1]) - 1;

    TProfile* pCorTwoDif = p2CorTwoDif->ProfileY(nameCorTwo.Data(),binMultLow,binMultHigh);
    TProfile* pCorFourDif = nullptr; if(bDoFour) { pCorFourDif = p2CorFourDif->ProfileY(nameCorFour.Data(),binMultLow,binMultHigh); }

    if(bDo3sub){
      for(Int_t poiPos(0); poiPos < 3; poiPos++){
        for(Int_t rfPos(0); rfPos < 3; rfPos++){
          if(poiPos == rfPos) continue;
          pCorTwoDif3sub[poiPos][rfPos] = p2CorTwoDif3sub[poiPos][rfPos]->ProfileY(nameCorTwo3sub[poiPos][rfPos].Data(),binMultLow,binMultHigh);
        }
      }
      if(bDoFour){
        for(Int_t poiPos(0); poiPos < 3; poiPos++){
          for(Int_t twoPos(0); twoPos < 3; twoPos++){
            pCorFourDif3sub[poiPos][twoPos] = p2CorFourDif3sub[poiPos][twoPos]->ProfileY(nameCorFour3sub[poiPos][twoPos].Data(),binMultLow,binMultHigh);
          }
        }
      }
    }

    // rebinning according to pt bins
    if(task->fNumPtBins > 0)
    {
      pCorTwoDif = (TProfile*) pCorTwoDif->Rebin(task->fNumPtBins,Form("%s_sample%d_rebin", nameCorTwo.Data(), iSample), task->fPtBinsEdges.data());
      if(bDoFour) { pCorFourDif = (TProfile*) pCorFourDif->Rebin(task->fNumPtBins,Form("%s_sample%d_rebin", nameCorFour.Data(), iSample), task->fPtBinsEdges.data()); }
      if(task->Has3sub()){
        for(Int_t poiPos(0); poiPos < 3; poiPos++){
          for(Int_t rfPos(0); rfPos < 3; rfPos++){
            if(poiPos == rfPos) continue;
            pCorTwoDif3sub[poiPos][rfPos] = (TProfile*) pCorTwoDif3sub[poiPos][rfPos]->Rebin(task->fNumPtBins,Form("%s_sample%d_rebin", nameCorTwo3sub[poiPos][rfPos].Data(), iSample), task->fPtBinsEdges.data());
          }
        }
        if(bDoFour){
          for(Int_t poiPos(0); poiPos < 3; poiPos++){
            for(Int_t twoPos(0); twoPos < 3; twoPos++){
              pCorFourDif3sub[poiPos][twoPos] = (TProfile*) pCorFourDif3sub[poiPos][twoPos]->Rebin(task->fNumPtBins,Form("%s_sample%d_rebin", nameCorFour3sub[poiPos][twoPos].Data(), iSample), task->fPtBinsEdges.data());
            }
          }
        }
      } //task 3 sub
    }
    else
    {
      pCorTwoDif = (TProfile*) pCorTwoDif->Clone(Form("%s_sample%d_rebin", nameCorTwo.Data(), iSample));
      if(bDoFour) { pCorFourDif = (TProfile*) pCorFourDif->Clone(Form("%s_sample%d_rebin", nameCorFour.Data(), iSample)); }
      if(bDo3sub){
        for(Int_t poiPos(0); poiPos < 3; poiPos++){
          for(Int_t rfPos(0); rfPos < 3; rfPos++){
            if(poiPos == rfPos) continue;
            pCorTwoDif3sub[poiPos][rfPos] = (TProfile*) pCorTwoDif3sub[poiPos][rfPos]->Clone(Form("%s_sample%d_rebin", nameCorTwo3sub[poiPos][rfPos].Data(), iSample));
          }
        }
        if(bDoFour){
          for(Int_t poiPos(0); poiPos < 3; poiPos++){
            for(Int_t twoPos(0); twoPos < 3; twoPos++){
              pCorFourDif3sub[poiPos][twoPos] = (TProfile*) pCorFourDif3sub[poiPos][twoPos]->Clone(Form("%s_sample%d_rebin", nameCorFour3sub[poiPos][twoPos].Data(), iSample));
            }
          }
        }
      }
    }

    // renaming
    TString sGap = TString(); if(task->HasGap()  && !bDo3sub) { sGap.Append(Form("{|#Delta#eta| > %g}",task->fEtaGap)); }
    pCorTwoDif->SetTitle(Form("%s: <<2'>>_{%d} %s; #it{p}_{T} (GeV/#it{c});",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data()));
    if(bDoFour) { pCorFourDif->SetTitle(Form("%s: <<4'>>_{%d} %s; #it{p}_{T} (GeV/#it{c});",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data())); }

    // NOTE: Here the <X'> is ready & rebinned
    listCorTwo->Add(pCorTwoDif);
    if(bDoFour) { listCorFour->Add(pCorFourDif); }

    if(bDo3sub){
      sGap.Append(Form(" 3 subevents (#eta mid: -%g,%g)",task->fEtaGap/2,task->fEtaGapSecond/2));
      for(Int_t poiPos(0); poiPos < 3; poiPos++){
        for(Int_t rfPos(0); rfPos < 3; rfPos++){
          if(poiPos == rfPos) continue;
          pCorTwoDif3sub[poiPos][rfPos]->SetTitle(Form("%s: <<2'>>_{%d}  %s (poi %c, rfp %c); #it{p}_{T} (GeV/#it{c});",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data(), sides[poiPos], sides[rfPos]));
          listCorTwo3sub[poiPos][rfPos]->Add(pCorTwoDif3sub[poiPos][rfPos]);
        }
      }
      if(bDoFour){
        for(Int_t poiPos(0); poiPos < 3; poiPos++){
          for(Int_t twoPos(0); twoPos < 3; twoPos++){
            pCorFourDif3sub[poiPos][twoPos]->SetTitle(Form("%s: <<4'>>_{%d}  %s (poi %c, two %c); #it{p}_{T} (GeV/#it{c});",GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, sGap.Data(), sides[poiPos], sides[twoPos]));
            listCorFour3sub[poiPos][twoPos]->Add(pCorFourDif3sub[poiPos][twoPos]);
          }
        }
      }
    }

    Debug("Calculating flow","ProcessDirect");
    // loading reference vn{2}
    TH1D* hFlowRefTwo = (TH1D*) listRefTwo->FindObject(Form("Refs_hFlow2_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
    if(!hFlowRefTwo) { Error(Form("Histo 'hFlowRefTwo' (sample %d) does not exists",iSample),"ProcessDirect"); listRefTwo->ls(); ffDesampleFile->ls(); return kFALSE; }

    // dn{2}
    TH1D* hCumTwoDif = CalcDifCumTwo(pCorTwoDif, task);
    if(!hCumTwoDif) { Error(Form("dn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
    hCumTwoDif->SetName(Form("%s_sample%d", nameCumTwo.Data(), iSample));
    listCumTwo->Add(hCumTwoDif);

    // v'n{2}
    TH1D* hFlowTwoDif = CalcDifFlowTwo(hCumTwoDif, hFlowRefTwo, iMultBin+1, task, bCorrelated);
    if(!hFlowTwoDif) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
    hFlowTwoDif->SetName(Form("%s_sample%d", nameFlowTwo.Data(), iSample));
    listFlowTwo->Add(hFlowTwoDif);

    TH1D* hFlowRefTwo3sub = nullptr;
    TH1D* hCumTwoDif3sub = nullptr;
    TH1D* hFlowTwoDif3sub = nullptr;
    if(bDo3sub){
      for(Int_t poiPos(0); poiPos < 3; poiPos++){
        for(Int_t rfPos(0); rfPos < 3; rfPos++){
          if(poiPos == rfPos) continue;
          if(poiPos < rfPos) hFlowRefTwo3sub = (TH1D*) listRefTwo3sub[poiPos][rfPos]->FindObject(Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[poiPos],sides[rfPos],iSample));
          else hFlowRefTwo3sub = (TH1D*) listRefTwo3sub[poiPos][rfPos]->FindObject(Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rfPos],sides[poiPos],iSample));
          if(!hFlowRefTwo3sub) { Error(Form("Histo 'hFlowRefTwo3sub' (sample %d) does not exists",iSample),"ProcessDirect"); listRefTwo3sub[poiPos][rfPos]->ls(); ffDesampleFile->ls(); return kFALSE; }

          // dn{2}
          hCumTwoDif3sub = CalcDifCumTwo(pCorTwoDif3sub[poiPos][rfPos], task);
          if(!hCumTwoDif3sub) { Error(Form("dn{2} 3sub (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
          hCumTwoDif3sub->SetName(Form("%s_sample%d", nameCumTwo3sub[poiPos][rfPos].Data(), iSample));
          listCumTwo3sub[poiPos][rfPos]->Add(hCumTwoDif3sub);

          // v'n{2}
          hFlowTwoDif3sub = CalcDifFlowTwo(hCumTwoDif3sub, hFlowRefTwo3sub, iMultBin+1, task, bCorrelated);
          if(!hFlowTwoDif) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
          hFlowTwoDif3sub->SetName(Form("%s_sample%d", nameFlowTwo3sub[poiPos][rfPos].Data(), iSample));
          listFlowTwo3sub[poiPos][rfPos]->Add(hFlowTwoDif3sub);
        }
      }
    }

    if(bDoFour)
    {
      // loading reference <<2>>
      TProfile* pCorTwoRef = (TProfile*) listRefCorTwo->FindObject(Form("Refs_pCor2_harm%d_gap%s_sample%d_rebin",task->fHarmonics, task->GetEtaGapString().Data(), iSample));
      if(!pCorTwoRef) { Error(Form("Profile 'pCorTwoRef' (sample %d) does not exists",iSample),"ProcessDirect"); listRefCorTwo->ls(); return kFALSE; }

      // loading reference vn{4}
      TH1D* hFlowRefFour = (TH1D*) listRefFour->FindObject(Form("Refs_hFlow4_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
      if(!hFlowRefFour) { Error(Form("Histo 'hFlowRefFour' (sample %d) does not exists",iSample),"ProcessDirect"); listRefFour->ls(); return kFALSE; }

      // dn{4}
      TH1D* hCumFourDif = CalcDifCumFour(pCorFourDif, pCorTwoDif, pCorTwoRef, iMultBin+1, task, bCorrelated);
      if(!hCumFourDif) { Error(Form("dn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hCumFourDif->SetName(Form("%s_sample%d", nameCumFour.Data(), iSample));
      listCumFour->Add(hCumFourDif);

      // v'n{4}
      TH1D* hFlowFourDif = CalcDifFlowFour(hCumFourDif, hFlowRefFour, iMultBin+1, task, bCorrelated);
      if(!hFlowFourDif) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hFlowFourDif->SetName(Form("%s_sample%d", nameFlowFour.Data(), iSample));
      listFlowFour->Add(hFlowFourDif);
    }

    if(bDoFour && bDo3sub){
      TH1D* mergedAllCombi = nullptr;
      TList* listAllCombi = new TList();

      Debug("Processing bDoFour && bDo3sub", "ProcessDirect");
      for(Int_t poiPos(0); poiPos < 3; poiPos++){
        for(Int_t twoPos(0); twoPos < 3; twoPos++){
          // loading reference vn{4}
          TH1D* hFlowRefFour3sub = (TH1D*) listRefFour3sub[twoPos]->FindObject(Form("Refs_hFlow4_harm%d_gap(%s,%s)_3sub_two_%c_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[twoPos],iSample));
          if(!hFlowRefFour3sub) { Error(Form("Histo 'hFlowRefFour' (sample %d) does not exists for two %c",iSample,sides[twoPos]),"ProcessDirect"); listRefFour3sub[twoPos]->ls(); return kFALSE; }

          TH1D* hCumFourDif3sub = nullptr;

          if(poiPos != twoPos){
            Debug("Processing poiPos != twoPos", "ProcessDirect");

            Int_t third = ReturnThird(poiPos,twoPos);
            // loading reference <<2>>
            TProfile* pCorTwoRef3sub = nullptr;
            if(twoPos < third) pCorTwoRef3sub = (TProfile*) listRefCorTwo3sub[twoPos][third]->FindObject(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_Pos_sample%d_rebin",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[twoPos],sides[third],iSample));
            else pCorTwoRef3sub = (TProfile*) listRefCorTwo3sub[twoPos][third]->FindObject(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_Pos_sample%d_rebin",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[third],sides[twoPos],iSample));
            if(!pCorTwoRef3sub) { Error(Form("Profile 'pCorTwoRef' (sample %d, poi %c, two %c) does not exists",iSample, sides[twoPos],sides[third]),"ProcessDirect");  listRefCorTwo3sub[poiPos][twoPos]->ls(); return kFALSE; }

            // dn{4}
            hCumFourDif3sub = CalcDifCumFour(pCorFourDif3sub[poiPos][twoPos], pCorTwoDif3sub[poiPos][twoPos], pCorTwoRef3sub, iMultBin+1, task, bCorrelated);
            if(!hCumFourDif3sub) { Error(Form("dn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
            hCumFourDif3sub->SetName(Form("%s_sample%d", nameCumFour3sub[poiPos][twoPos].Data(), iSample));
            listCumFour3sub[poiPos][twoPos]->Add(hCumFourDif3sub);
          } // end poiPos != twoPos
          else {
            Debug("Processing poiPos == twoPos", "ProcessDirect");
            Int_t indexUp = ReturnIndex3sub(poiPos);
            Int_t third = ReturnThird(poiPos,indexUp);

            TProfile* pCorTwoRef1 = nullptr;
            TProfile* pCorTwoRef2 = nullptr;

            if(poiPos < indexUp) pCorTwoRef1 = (TProfile*) listRefCorTwo3sub[poiPos][indexUp]->FindObject(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_Pos_sample%d_rebin",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[poiPos],sides[indexUp],iSample));
            else pCorTwoRef1 = (TProfile*) listRefCorTwo3sub[poiPos][indexUp]->FindObject(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_Pos_sample%d_rebin",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[indexUp],sides[poiPos],iSample));
            if(!pCorTwoRef1) { Error(Form("Profile 'pCorTwoRef' (sample %d, poi %c, two %c) does not exists",iSample, sides[poiPos],sides[indexUp]),"ProcessDirect");  listRefCorTwo3sub[poiPos][indexUp]->ls(); return kFALSE; }

            if(poiPos < third) pCorTwoRef2 = (TProfile*) listRefCorTwo3sub[poiPos][third]->FindObject(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_Pos_sample%d_rebin",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[poiPos],sides[third],iSample));
            else pCorTwoRef2 = (TProfile*) listRefCorTwo3sub[poiPos][third]->FindObject(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c_Pos_sample%d_rebin",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[third],sides[poiPos],iSample));
            if(!pCorTwoRef2) { Error(Form("Profile 'pCorTwoRef' (sample %d, poi %c, two %c) does not exists",iSample, sides[poiPos],sides[third]),"ProcessDirect");  listRefCorTwo3sub[poiPos][third]->ls(); return kFALSE; }

            Debug("Loaded ref", "ProcessDirect");

            // dn{4}
            hCumFourDif3sub = CalcDifCumFour3sub(pCorFourDif3sub[poiPos][twoPos], pCorTwoDif3sub[poiPos][indexUp], pCorTwoDif3sub[poiPos][third], pCorTwoRef2, pCorTwoRef1, iMultBin+1, task);
            if(!hCumFourDif3sub) { Error(Form("dn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
            hCumFourDif3sub->SetName(Form("%s_sample%d", nameCumFour3sub[poiPos][twoPos].Data(), iSample));
            listCumFour3sub[poiPos][twoPos]->Add(hCumFourDif3sub);

            Debug("4p 3 sub: cumulant calculated", "ProcessDirect");
          } // end poiPos == twoPos

          // v'n{4}
          TH1D* hFlowFourDif3sub = CalcDifFlowFour(hCumFourDif3sub, hFlowRefFour3sub, iMultBin+1, task, bCorrelated);
          if(!hFlowRefFour3sub) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
          hFlowFourDif3sub->SetName(Form("%s_sample%d", nameFlowFour3sub[poiPos][twoPos].Data(), iSample));
          listFlowFour3sub[poiPos][twoPos]->Add(hFlowFourDif3sub);

          //combination of all positions of poi & 2pc
          listAllCombi->Add(hFlowFourDif3sub);
          if(!mergedAllCombi) {mergedAllCombi = (TH1D*) hFlowFourDif3sub->Clone(Form("%s_allC_sample%d",nameFlowAllCombi.Data(),iSample)); mergedAllCombi->Reset(); }
          mergedAllCombi->Add(hFlowFourDif3sub);
        }
      }
      if(!mergedAllCombi) { Error("Merging of 'hFlowFourDif3sub' failed!","ProcessDirect"); return kFALSE; }

      mergedAllCombi->Scale((Double_t) 1./9);
      mergedAllCombi->SetName(Form("%s_merged_sample%d",nameFlowAllCombi.Data(),iSample));
      TString desName = Form("%s_sample%d", nameFlowAllCombi.Data(),iSample);
      desAllCombi[iSample] = DesampleList(listAllCombi, mergedAllCombi, task, desName);
      if(!desAllCombi[iSample]) { Error("Desampling vn{2} failed","ProcessDirect"); return kFALSE; }
      desAllCombi[iSample]->SetName(Form("%s_des_sample%d",nameFlowAllCombi.Data(),iSample));

      listMergedAllCombinations->Add(desAllCombi[iSample]);
      delete listAllCombi;
    } //end v'n{4} with 3 sub-events
  } // end-for {iSample} : loop over samples

  Debug("Merging correlations for central values", "ProcessDirect");

  // loading reference vn{2} merged for vn{2} dif. merged
  TH1D* hFlowTwoRefMerged = (TH1D*) ffOutputFile->Get(Form("Refs_hFlow2_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!hFlowTwoRefMerged) { Error(Form("Reference vn{2} (merged) not loaded!"),"ProcessDirect"); return kFALSE; }

  // <<2>>
  TProfile* pCorTwoMerged = (TProfile*) MergeListProfiles(listCorTwo);
  if(!pCorTwoMerged) { Error("Merging of 'pCorTwoMerged' failed!","ProcessDirect"); return kFALSE; }
  pCorTwoMerged->SetName(Form("%s_merged", nameCorTwo.Data()));

  // dn{2}
  TH1D* hCumTwoMerged = CalcDifCumTwo(pCorTwoMerged, task);
  if(!hCumTwoMerged) { Error(Form("dn{2} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
  hCumTwoMerged->SetName(Form("%s_merged", nameCumTwo.Data()));

  // vn{2}
  TH1D* hFlowTwoMerged = CalcDifFlowTwo(hCumTwoMerged, hFlowTwoRefMerged, iMultBin+1, task);
  if(!hFlowTwoMerged) { Error(Form("vn{2} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
  hFlowTwoMerged->SetName(Form("%s_merged", nameFlowTwo.Data()));

  TH1D* hFlowTwoRefMerged3sub[3][3] = {nullptr};
  TProfile* pCorTwoMerged3sub[3][3] = {nullptr};
  TH1D* hCumTwoMerged3sub[3][3] = {nullptr};
  TH1D* hFlowTwoMerged3sub[3][3] = {nullptr};
  if(bDo3sub){
    for(Int_t poiPos(0); poiPos < 3; poiPos++){
      for(Int_t rfPos(0); rfPos < 3; rfPos++){
        if(poiPos == rfPos) continue;
        if(poiPos < rfPos) hFlowTwoRefMerged3sub[poiPos][rfPos] = (TH1D*) ffOutputFile->Get(Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[poiPos],sides[rfPos]));
        else hFlowTwoRefMerged3sub[poiPos][rfPos] = (TH1D*) ffOutputFile->Get(Form("Refs_hFlow2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[rfPos],sides[poiPos]));
        if(!hFlowTwoRefMerged3sub[poiPos][rfPos]) { Error(Form("Reference vn{2} (%c,%c) (merged) not loaded!",sides[poiPos],sides[rfPos]),"ProcessDirect"); return kFALSE; }

        //<<2>>
        pCorTwoMerged3sub[poiPos][rfPos] = (TProfile*) MergeListProfiles(listCorTwo3sub[poiPos][rfPos]);
        if(!pCorTwoMerged3sub[poiPos][rfPos]) { Error("Merging of 'pCorTwoMerged' failed!","ProcessDirect"); return kFALSE; }
        pCorTwoMerged3sub[poiPos][rfPos]->SetName(Form("%s_merged", nameCorTwo3sub[poiPos][rfPos].Data()));

        // dn{2}
        hCumTwoMerged3sub[poiPos][rfPos] = CalcDifCumTwo(pCorTwoMerged3sub[poiPos][rfPos], task);
        if(!hCumTwoMerged3sub[poiPos][rfPos]) { Error(Form("dn{2} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
        hCumTwoMerged3sub[poiPos][rfPos]->SetName(Form("%s_merged", nameCumTwo3sub[poiPos][rfPos].Data()));

        // vn{2}
        hFlowTwoMerged3sub[poiPos][rfPos] = CalcDifFlowTwo(hCumTwoMerged3sub[poiPos][rfPos], hFlowTwoRefMerged3sub[poiPos][rfPos], iMultBin+1, task);
        if(!hFlowTwoMerged3sub[poiPos][rfPos]) { Error(Form("vn{2} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
        hFlowTwoMerged3sub[poiPos][rfPos]->SetName(Form("%s_merged", nameFlowTwo3sub[poiPos][rfPos].Data()));
      }
    }
  }

  TProfile* pCorFourMerged = nullptr;
  TH1D* hCumFourMerged = nullptr;
  TH1D* hFlowFourMerged = nullptr;
  if(bDoFour)
  {
    // loading reference <<2>> merged
    TProfile* pCorTwoRefMerged = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap%s", task->fHarmonics, task->GetEtaGapString().Data() ));
    if(!pCorTwoRefMerged) { Error(Form("Reference <<2>> (merged) not loaded!"),"ProcessDirect"); ffOutputFile->ls(); return kFALSE; }

    // loading reference vn{4} merged
    TH1D* hFlowFourRefMerged = (TH1D*) ffOutputFile->Get(Form("Refs_hFlow4_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!hFlowFourRefMerged) { Error(Form("Reference vn{4} (merged) not loaded!"),"ProcessDirect"); ffOutputFile->ls(); return kFALSE; }

    // <<4>>
    pCorFourMerged = (TProfile*) MergeListProfiles(listCorFour);
    if(!pCorFourMerged) { Error("Merging of 'pCorFourMerged' failed!","ProcessDirect"); return kFALSE; }
    pCorFourMerged->SetName(Form("%s_merged", nameCorFour.Data()));

    // dn{4}
    hCumFourMerged = CalcDifCumFour(pCorFourMerged, pCorTwoMerged, pCorTwoRefMerged, iMultBin+1, task, bCorrelated);
    if(!hCumFourMerged) { Error(Form("cn{4} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
    hCumFourMerged->SetName(Form("%s_merged", nameCumFour.Data()));

    // vn{4}
    hFlowFourMerged = CalcDifFlowFour(hCumFourMerged, hFlowFourRefMerged, iMultBin+1, task, bCorrelated);
    if(!hFlowFourMerged) { Error(Form("vn{4} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
    hFlowFourMerged->SetName(Form("%s_merged", nameFlowFour.Data()));
  }

  TProfile* pCorFourMerged3sub[3][3] = {nullptr};
  TH1D* hCumFourMerged3sub[3][3] = {nullptr};
  TH1D* hFlowFourMerged3sub[3][3] = {nullptr};
  if(bDoFour && bDo3sub){
    for(Int_t poiPos(0); poiPos < 3; poiPos++){
      for(Int_t twoPos(0); twoPos < 3; twoPos++){
        // loading reference vn{4} merged
        TH1D* hFlowRefFourMerged3sub = (TH1D*) ffOutputFile->Get(Form("Refs_hFlow4_harm%d_gap(%s,%s)_3sub_two_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[twoPos]));
        if(!hFlowRefFourMerged3sub) { Error(Form("Histo 'hFlowRefFour' (merged) does not exists for two %c",sides[twoPos]),"ProcessDirect"); ffOutputFile->ls(); return kFALSE; }

        // merging <<4>>
        pCorFourMerged3sub[poiPos][twoPos] = (TProfile*) MergeListProfiles(listCorFour3sub[poiPos][twoPos]);
        if(!pCorFourMerged3sub[poiPos][twoPos]) { Error("Merging of 'pCorFourMerged' failed!","ProcessDirect"); return kFALSE; }
        pCorFourMerged3sub[poiPos][twoPos]->SetName(Form("%s_merged", nameCorFour3sub[poiPos][twoPos].Data()));


        if(poiPos != twoPos){
          Debug("Merging poiPos != twoPos", "ProcessDirect");

          Int_t third = ReturnThird(poiPos,twoPos);

          // loading reference <<2>> merged
          TProfile* pCorTwoRefMerged3sub = nullptr;
          if(twoPos < third) pCorTwoRefMerged3sub = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c", task->fHarmonics, task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[twoPos],sides[third]));
          else pCorTwoRefMerged3sub = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c", task->fHarmonics, task->GetEtaGapString().Data(), task->GetEtaGapString().Data(),sides[third],sides[twoPos]));
          if(!pCorTwoRefMerged3sub) { Error(Form("Reference <<2>> (%c,%c) (merged) not loaded!",sides[twoPos],sides[third]),"ProcessDirect"); ffOutputFile->ls(); return kFALSE; }

          // dn{4}
          hCumFourMerged3sub[poiPos][twoPos] = CalcDifCumFour(pCorFourMerged3sub[poiPos][twoPos], pCorTwoMerged3sub[poiPos][twoPos], pCorTwoRefMerged3sub, iMultBin+1, task, bCorrelated);
          if(!hCumFourMerged3sub[poiPos][twoPos]) { Error("dn{4} (merged) not processed correctly!","ProcessDirect"); return kFALSE; }
          hCumFourMerged3sub[poiPos][twoPos]->SetName(Form("%s_merged", nameCumFour3sub[poiPos][twoPos].Data()));
        } // end poiPos != twoPos
        else {
          Debug("Merging poiPos == twoPos", "ProcessDirect");
          Int_t indexUp = ReturnIndex3sub(poiPos);
          Int_t third = ReturnThird(poiPos,indexUp);

          printf("poiPos %d twoPos %d indexUp %d third %d \n", poiPos, twoPos, indexUp, third);

          TProfile* pCorTwoRef1Merged = nullptr;
          TProfile* pCorTwoRef2Merged = nullptr;

          if(poiPos < indexUp) pCorTwoRef1Merged = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[poiPos],sides[indexUp]));
          else pCorTwoRef1Merged = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[indexUp],sides[poiPos]));
          if(!pCorTwoRef1Merged) { Error(Form("Profile 'pCorTwoRef' (poi %c, two %c) does not exists", sides[poiPos],sides[indexUp]),"ProcessDirect");  ffOutputFile->ls(); return kFALSE; }

          Debug("Ref1 loaded", "ProcessDirect");


          if(poiPos < third) pCorTwoRef2Merged = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[poiPos],sides[third]));
          else pCorTwoRef2Merged = (TProfile*) ffOutputFile->Get(Form("Refs_pCor2_harm%d_gap(%s,%s)_3sub_rf1_%c_rf2_%c",task->fHarmonics,task->GetEtaGapString().Data(),task->GetEtaGapString().Data(),sides[third],sides[poiPos]));
          if(!pCorTwoRef2Merged) { Error(Form("Profile 'pCorTwoRef' (poi %c, two %c) does not exists", sides[poiPos],sides[third]),"ProcessDirect");  ffOutputFile->ls(); return kFALSE; }

          Debug("Ref2 loaded", "ProcessDirect");

          // dn{4}
          hCumFourMerged3sub[poiPos][twoPos] = CalcDifCumFour3sub(pCorFourMerged3sub[poiPos][twoPos], pCorTwoMerged3sub[poiPos][indexUp], pCorTwoMerged3sub[poiPos][third], pCorTwoRef2Merged, pCorTwoRef1Merged, iMultBin+1, task);
          if(!hCumFourMerged3sub[poiPos][twoPos]) { Error("dn{4} (merged) not processed correctly!","ProcessDirect"); return kFALSE; }
          hCumFourMerged3sub[poiPos][twoPos]->SetName(Form("%s_merged", nameCumFour3sub[poiPos][twoPos].Data()));

          Debug("Cum calculated", "ProcessDirect");
        } // end poiPos == twoPos
        // v'n{4}
        hFlowFourMerged3sub[poiPos][twoPos] = CalcDifFlowFour(hCumFourMerged3sub[poiPos][twoPos], hFlowRefFourMerged3sub, iMultBin+1, task, bCorrelated);            if(!hFlowFourMerged3sub[poiPos][twoPos]) { Error("vn{4} (merged) not processed correctly!","ProcessDirect"); return kFALSE; }
        hFlowFourMerged3sub[poiPos][twoPos]->SetName(Form("%s_merged", nameFlowFour3sub[poiPos][twoPos].Data()));
      }
    }
  }

  Debug("Desampling","ProcessDirect");

  Debug(Form("<<2'>>: Number of samples in list pre-merging: %d",listCorTwo->GetEntries()),"ProcessDirect");
  TH1D* hDesampledTwo_Cor = DesampleList(listCorTwo, pCorTwoMerged->ProjectionX(), task, nameCorTwo, kTRUE); // skipping desampling for output structure
  if(!hDesampledTwo_Cor) { Error("Desampling <<2'>> failed","ProcessDirect"); return kFALSE; }
  hDesampledTwo_Cor->SetName(nameCorTwo.Data());

  Debug(Form("dn{2}: Number of samples in list pre-merging: %d",listCumTwo->GetEntries()),"ProcessDirect");
  TH1D* hDesampledTwo_Cum = DesampleList(listCumTwo, hCumTwoMerged, task, nameCumTwo);
  if(!hDesampledTwo_Cum) { Error("Desampling dn{2} failed","ProcessDirect"); return kFALSE; }
  hDesampledTwo_Cum->SetName(nameCumTwo.Data());

  Debug(Form("vn{2}: Number of samples in list pre-merging: %d",listFlowTwo->GetEntries()),"ProcessDirect");
  TH1D* hDesampledTwo = DesampleList(listFlowTwo, hFlowTwoMerged, task, nameFlowTwo);
  if(!hDesampledTwo) { Error("Desampling vn{2} failed","ProcessDirect"); return kFALSE; }
  hDesampledTwo->SetName(nameFlowTwo.Data());

  // saving to output file & cleaning
  ffOutputFile->cd();
  if(fSaveInterSteps) {
    hDesampledTwo_Cor->Write();
    hDesampledTwo_Cum->Write();
  }
  hDesampledTwo->Write();

  delete hDesampledTwo_Cum;
  // delete hDesampledTwo;

  if(bDoFour)
  {
    Debug(Form("<<4'>>: Number of samples in list pre-merging: %d",listCorFour->GetEntries()),"ProcessDirect");
    TH1D* hDesampledFour_Cor = DesampleList(listCorFour, pCorFourMerged->ProjectionX(), task, nameCorFour, kTRUE); // skipping desampling for structure
    if(!hDesampledFour_Cor) { Error("Desampling dn{4} failed","ProcessDirect"); return kFALSE; }
    hDesampledFour_Cor->SetName(nameCorFour.Data());

    Debug(Form("dn{4}: Number of samples in list pre-merging: %d",listCumFour->GetEntries()),"ProcessDirect");
    TH1D* hDesampledFour_Cum = DesampleList(listCumFour, hCumFourMerged, task, nameCumFour);
    if(!hDesampledFour_Cum) { Error("Desampling dn{4} failed","ProcessDirect"); return kFALSE; }
    hDesampledFour_Cum->SetName(nameCumFour.Data());

    Debug(Form("vn{4}: Number of samples in list pre-merging: %d",listFlowFour->GetEntries()),"ProcessDirect");
    TH1D* hDesampledFour = DesampleList(listFlowFour, hFlowFourMerged, task, nameFlowFour);
    if(!hDesampledFour) { Error("Desampling vn{4} failed","ProcessDirect"); return kFALSE; }
    hDesampledFour->SetName(nameFlowFour.Data());

    // saving to output file & cleaning
    ffOutputFile->cd();
    if(fSaveInterSteps) {
      hDesampledFour_Cor->Write();
      hDesampledFour_Cum->Write();
    }
    hDesampledFour->Write();

    delete hDesampledFour_Cor;
    delete hDesampledFour_Cum;
    delete hDesampledFour;
  };

  if(bDo3sub){
    for(Int_t poiPos(0); poiPos < 3; poiPos++){
      for(Int_t rfPos(0); rfPos < 3; rfPos++){
        if(poiPos == rfPos) continue;
        Debug(Form("<<2'>> 3 sub (poi: %c, rf: %c): Number of samples in list pre-merging: %d",sides[poiPos],sides[rfPos],listCorTwo3sub[poiPos][rfPos]->GetEntries()),"ProcessDirect");
        TH1D* hDesampledTwo_Cor3sub = DesampleList(listCorTwo3sub[poiPos][rfPos], pCorTwoMerged3sub[poiPos][rfPos]->ProjectionX(), task, nameCorTwo3sub[poiPos][rfPos]);
        if(!hDesampledTwo_Cor3sub) { Error("Desampling <<2'>> failed","ProcessDirect"); return kFALSE; }
        hDesampledTwo_Cor3sub->SetName(nameCorTwo3sub[poiPos][rfPos].Data());

        TH1D* hDesampledTwo_Cum3sub = DesampleList(listCumTwo3sub[poiPos][rfPos], hCumTwoMerged3sub[poiPos][rfPos], task, nameCumTwo3sub[poiPos][rfPos]);
        if(!hDesampledTwo_Cum3sub) { Error("Desampling dn{2} failed","ProcessDirect"); return kFALSE; }
        hDesampledTwo_Cum3sub->SetName(nameCumTwo3sub[poiPos][rfPos].Data());

        TH1D* hDesampledTwo3sub = DesampleList(listFlowTwo3sub[poiPos][rfPos], hFlowTwoMerged3sub[poiPos][rfPos], task, nameFlowTwo3sub[poiPos][rfPos]);
        if(!hDesampledTwo3sub) { Error("Desampling vn{2} failed","ProcessDirect"); return kFALSE; }
        hDesampledTwo3sub->SetName(nameFlowTwo3sub[poiPos][rfPos].Data());

        ffOutputFile->cd();
        if(fSaveInterSteps) {
          hDesampledTwo_Cor3sub->Write();
          hDesampledTwo_Cum3sub->Write();
        }
        hDesampledTwo3sub->Write();

        delete hDesampledTwo_Cor3sub;
        delete hDesampledTwo_Cum3sub;
        delete hDesampledTwo3sub;
      }
    }
  }

  if(bDoFour && bDo3sub){
    for(Int_t poiPos(0); poiPos < 3 ; poiPos++){
      for(Int_t twoPos(0); twoPos < 3; twoPos++){
        Debug(Form("<<4'>> 3 sub (poi: %c, two: %c): Number of samples in list pre-merging: %d",sides[poiPos],sides[twoPos],listCorFour3sub[poiPos][twoPos]->GetEntries()),"ProcessDirect");

        TH1D* hDesampledFour_Cor3sub = DesampleList(listCorFour3sub[poiPos][twoPos], pCorFourMerged3sub[poiPos][twoPos]->ProjectionX(), task, nameCorFour3sub[poiPos][twoPos]);
        if(!hDesampledFour_Cor3sub) { Error("Desampling cn{4} failed","ProcessDirect"); return kFALSE; }
        hDesampledFour_Cor3sub->SetName(nameCorFour3sub[poiPos][twoPos].Data());
        Debug("Desampling pCor done","ProcessDirect");

        TH1D* hDesampledFour_Cum3sub = DesampleList(listCumFour3sub[poiPos][twoPos], hCumFourMerged3sub[poiPos][twoPos], task, nameCumFour3sub[poiPos][twoPos]);
        if(!hDesampledFour_Cum3sub) { Error("Desampling dn{4} failed","ProcessDirect"); return kFALSE; }
        hDesampledFour_Cum3sub->SetName(nameCumFour3sub[poiPos][twoPos].Data());

        Debug("Desampling cum 4 done","ProcessDirect");

        TH1D* hDesampledFour3sub = DesampleList(listFlowFour3sub[poiPos][twoPos], hFlowFourMerged3sub[poiPos][twoPos], task, nameFlowFour3sub[poiPos][twoPos]);
        if(!hDesampledFour3sub) { Error("Desampling vn{4} failed","ProcessDirect"); return kFALSE; }
        hDesampledFour3sub->SetName(nameFlowFour3sub[poiPos][twoPos].Data());

        Debug("Desampling flow 4 done","ProcessDirect");
        Debug("Merging & desampling flow for all combinations","ProcessDirect");

        TH1D* mergedAllCombiAllSamples = nullptr;
        for(Int_t iSample(0); iSample < nOfSamples; iSample++){
          if(!mergedAllCombiAllSamples) {mergedAllCombiAllSamples = (TH1D*) desAllCombi[iSample]->Clone(Form("%s_allC",nameFlowAllCombi.Data())); mergedAllCombiAllSamples->Reset(); }
          mergedAllCombiAllSamples->Add(desAllCombi[iSample]);
        }
        mergedAllCombiAllSamples->Scale((Double_t) 1./nOfSamples);
        mergedAllCombiAllSamples->SetName(Form("%s_merged",nameFlowAllCombi.Data()));

        TH1D* hDesampledAllC = DesampleList(listMergedAllCombinations, mergedAllCombiAllSamples, task, nameFlowAllCombi);
        if(!hDesampledAllC) { Error("Desampling vn{4} failed","ProcessDirect"); return kFALSE; }
        hDesampledAllC->SetName(nameFlowAllCombi.Data());

        ffOutputFile->cd();
        if(fSaveInterSteps) {
          hDesampledFour_Cor3sub->Write();
          hDesampledFour_Cum3sub->Write();
        }
        hDesampledFour3sub->Write();
        hDesampledAllC->Write();

        // ffOutputFile->cd();
        // if(fSaveInterSteps) {
        //   hDesampledFour_Cor3sub->Write();
        //   hDesampledFour_Cum3sub->Write();
        // }
        // hDesampledFour3sub->Write();

        delete hDesampledFour_Cor3sub;
        delete hDesampledFour_Cum3sub;
        delete hDesampledFour3sub;
      }
    }
  }

  delete listCorTwo;
  delete listCumTwo;
  delete listFlowTwo;

  delete listCorFour;
  delete listCumFour;
  delete listFlowFour;

  for(Int_t poiPos(0); poiPos < 3; poiPos++){
    for(Int_t rfPos(0); rfPos < 3; rfPos++){
      delete listRefTwo3sub[poiPos][rfPos];
      delete listCorTwo3sub[poiPos][rfPos];
      delete listCumTwo3sub[poiPos][rfPos];
      delete listFlowTwo3sub[poiPos][rfPos];
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessReconstructed(FlowTask* task,Short_t iMultBin)
{
  Debug("Processing task","ProcessReconstructed");
  if(!task) { Error("Task not valid!","ProcessReconstructed"); return kFALSE; }
  if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessReconstructed"); return kFALSE; }

  if(fFlowFitCumulants) { fFlowFitCumulants = kFALSE; Warning("Fitting cumulants currently not available! WIP! switching flag off"); }

  // new naming convention for input histos (from FlowTask)
  TString sProfTwoName = Form("<<2>>(%d,-%d)",task->fHarmonics, task->fHarmonics);
  TString sProfFourName = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarmonics, task->fHarmonics, task->fHarmonics, task->fHarmonics);
  if(task->HasGap()) {
    sProfTwoName += Form("_2sub(%.2g)",task->fEtaGap);
    sProfFourName += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TString sSpeciesName = GetSpeciesName(task->fSpecies);
  TString sSpeciesLabel = GetSpeciesLabel(task->fSpecies);

  // ### Preparing slices of pt
  // if(!PrepareSlices(iMultBin,task,profFlow,histEntries,histEntriesBg,profFlowFour)) { return kFALSE; }

  if(!PrepareSlicesNew(task,sProfTwoName,kTRUE)) { Error(Form("PrepareSlicesNew '%s' failed!",sProfTwoName.Data()),"ProcessReconstructed"); return kFALSE; }
  if(task->fCumOrderMax >= 4 && !PrepareSlicesNew(task,sProfFourName,kFALSE)) { Error(Form("PrepareSlicesNew for '%s' failed!",sProfFourName.Data()),"ProcessReconstructed"); return kFALSE; }


  // things for later

  Bool_t bCorrelated = 0;

  // Loading list where reference flow samples are stored
  TList* listRefCorTwo = (TList*) ffDesampleFile->Get(Form("Refs_pCor2_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefCorTwo) { Error("List 'listRefCorTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefTwo = (TList*) ffDesampleFile->Get(Form("Refs_hFlow2_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefTwo) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefFour = nullptr;
  if(task->fCumOrderMax >= 4)
  {
    listRefFour = (TList*) ffDesampleFile->Get(Form("Refs_hFlow4_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!listRefFour) { Error("List 'listRefFour' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
  }


  Int_t iSample = 0;

  for(Int_t binMult(0); binMult < fiNumMultBins; ++binMult)
  {
    iMultBin = binMult;

    // List for desampling : later
    TString nameCorTwo = Form("%s_pCor2_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameCumTwo = Form("%s_hCum2_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameFlowTwo = Form("%s_hFlow2_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

    TString nameCorFour = Form("%s_pCor4_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameCumFour = Form("%s_hCum4_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameFlowFour = Form("%s_hFlow4_harm%d_gap%s_cent%d", GetSpeciesName(task->fSpecies).Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    //

    TH1D* hCumTwoDif = nullptr;
    TH1D* hFlowTwoDif = nullptr;
    TH1D* hCumFourDif = nullptr;
    TH1D* hFlowFourDif = nullptr;

    Double_t dFlow = 0.0, dFlowError = 0.0; // containers for flow extraction results
    TH1D* hFlowMass = nullptr;
    TH1D* pCorTwoDif = nullptr;
    if(!fFlowFitCumulants) { pCorTwoDif = new TH1D(nameCorTwo.Data(),Form("%s: <<2>>_{%d}{|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); <<2>>_{%d}{|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }
    else { pCorTwoDif = new TH1D(nameCumTwo.Data(),Form("%s: d_{%d}{2,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); d_{%d}{2,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }


    Double_t dFlowFour = 0.0, dFlowFourError = 0.0; // containers for flow extraction results
    TH1D* hFlowMassFour = nullptr;
    TH1D* pCorFourDif = nullptr;
    if(!fFlowFitCumulants) { pCorFourDif = new TH1D(nameCorFour.Data(),Form("%s: <<4>>_{%d}{|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); <<4>>_{%d}{|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }
    else { pCorFourDif = new TH1D(nameCumFour.Data(),Form("%s: d_{%d}{4,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); d_{%d}{4,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }


    for(Short_t binPt(0); binPt < task->fNumPtBins; binPt++)
    {
      TH1D* hInvMass = (TH1D*) task->fListHistos->FindObject(Form("hInvMass_mult%d_pt%d",binMult,binPt));
      if(!hInvMass) { Error("hInvMass histo not found among slices!","ProcessReconstructed"); task->fListHistos->ls(); return kFALSE; }
      hInvMass->SetTitle(Form("%s: InvMass dist (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
      hInvMass->SetMarkerStyle(kFullCircle);

      TH1D* hInvMassBg = nullptr;
      if(task->fSpecies == kPhi) {
         hInvMassBg = (TH1D*) task->fListHistos->FindObject(Form("hInvMassBg_mult%d_pt%d",binMult,binPt));
         if(!hInvMassBg) { Error("Loading inv. mass (Bg) slice failed!","ProcessReconstructed"); task->fListHistos->ls(); return kFALSE; }
      }

      TH1D* hFlowMass = (TH1D*) task->fListProfiles->FindObject(Form("%s_Pos_sample0_mult%d_pt%d",sProfTwoName.Data(),binMult,binPt));
      if(!hFlowMass) { Error(Form("hFlowMass histo '%s' not found among slices!",sProfTwoName.Data()),"ProcessReconstructed"); task->fListProfiles->ls(); return kFALSE; }
      hFlowMass->SetTitle(Form("%s: FlowMass (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
      hFlowMass->SetMarkerStyle(kFullCircle);

      if(task->fCumOrderMax >= 4)
      {
        hFlowMassFour = (TH1D*) task->fListProfiles->FindObject(Form("%s_Pos_sample0_mult%d_pt%d",sProfFourName.Data(),binMult,binPt));
        if(!hFlowMassFour) { Error(Form("hFlowMassFour histo '%s' not found among slices!",sProfFourName.Data()),"ProcessReconstructed"); task->fListProfiles->ls(); return kFALSE; }
        hFlowMassFour->SetTitle(Form("%s: FlowMassFour (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
        hFlowMassFour->SetMarkerStyle(kFullCircle);
      }

      // ### Fitting the correlations

      Bool_t bFitMass = kFALSE;
      Bool_t bFitCor = kFALSE;
      Bool_t bFitCorFour = kFALSE;

      TF1 fitOut, fitOutSig, fitOutBg;
      TF1 fitCor, fitCorSig, fitCorBg;
      TF1 fitCorFour, fitCorSigFour, fitCorBgFour;

      TList* listMass = new TList();
      listMass->SetOwner(1);

      bFitMass = FitInvMass(hInvMass, task, fitOut, fitOutSig, fitOutBg, listMass, hInvMassBg);

      TList* listFits = (TList*) listMass->Clone("listFits");
      listFits->SetOwner(1);

      TList* listFitsFour = (TList*) listMass->Clone("listFitsFour");
      listFitsFour->SetOwner(1);

      if(bFitMass) {
          bFitCor = FitCorrelations(hFlowMass, task, fitCor, fitCorSig, fitCorBg, fitOutSig, fitOutBg, listFits, kFALSE);
          if(task->fCumOrderMax >= 4) {
              bFitCorFour = FitCorrelations(hFlowMassFour, task, fitCorFour, fitCorSigFour, fitCorBgFour, fitOutSig, fitOutBg, listFitsFour, kTRUE);
        }
      }


    // Adding description strings
    TNamed* corr = new TNamed("corr",Form("%s",sProfTwoName.Data()));
    listFits->Add(corr);
    TNamed* spec = new TNamed("spec",Form("%s",sSpeciesLabel.Data()));
    listFits->Add(spec);
    TNamed* nCent = new TNamed("cent",Form("%d-%d",(Int_t)fdMultBins[binMult],(Int_t)fdMultBins[binMult+1]));
    listFits->Add(nCent);
    TNamed* nPt = new TNamed("pt",Form("%g-%g",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1]));
    listFits->Add(nPt);

    if(task->fCumOrderMax >= 4)
    {
      TNamed* corrFour = new TNamed("corr",Form("%s",sProfFourName.Data()));
      listFitsFour->Add(corrFour);
      listFitsFour->Add(spec);
      listFitsFour->Add(nCent);
      listFitsFour->Add(nPt);
    }

      // Saving fits to output file
      ffFitsFile->cd();
      listFits->Write(Form("%s_%s_cent%d_pt%d",sSpeciesName.Data(),sProfTwoName.Data(),iMultBin,binPt),TObject::kSingleKey);
      if(task->fCumOrderMax >= 4) { listFitsFour->Write(Form("%s_%s_cent%d_pt%d",sSpeciesName.Data(),sProfFourName.Data(),iMultBin,binPt),TObject::kSingleKey); }

      if(!bFitMass) {
          Error(Form("Fitting inv.mass unsuccesfull (mult %d | pt %d)",binMult,binPt),"ProcessReconstructed");
          return kFALSE;
      }

      if(!bFitCor) {
        Error(Form("Fitting vn-mass unsuccesfull (mult %d | pt %d)",binMult,binPt),"ProcessReconstructed");
        return kFALSE;
      }

      if(task->fCumOrderMax >= 4 && !bFitCorFour) {
        Error(Form("Fitting vn-mass for Four unsuccesfull (mult %d | pt %d)",binMult,binPt),"ProcessReconstructed");
        return kFALSE;
      }


        Int_t iParFlow = fitCorSig.GetNpar() - 1;
        dFlow = fitCorSig.GetParameter(iParFlow);
        dFlowError = fitCorSig.GetParError(iParFlow);

        if(task->fCumOrderMax >= 4) {
          Int_t iParFlow = fitCorSigFour.GetNpar() - 1;
          dFlowFour = fitCorSigFour.GetParameter(iParFlow);
          dFlowFourError = fitCorSigFour.GetParError(iParFlow);
        }

      // setting the flow
      pCorTwoDif->SetBinContent(binPt+1,dFlow);
      pCorTwoDif->SetBinError(binPt+1,dFlowError);

      Double_t dFlowRel = -999.9; if(TMath::Abs(dFlow) > 0.0) { dFlowRel = dFlowError / dFlow; }
      Info(Form("Final vn{2}: (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,binPt,dFlow,dFlowError,dFlowRel), "ProcessReconstructed");

      if(task->fCumOrderMax >= 4)
      {
        pCorFourDif->SetBinContent(binPt+1,dFlowFour);
        pCorFourDif->SetBinError(binPt+1,dFlowFourError);

        Double_t dFlowRel = -999.9; if(TMath::Abs(dFlowFour) > 0.0) { dFlowRel = dFlowFourError / dFlowFour; }
        Info(Form("Final vn{4}: (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,binPt,dFlowFour,dFlowFourError,dFlowRel), "ProcessReconstructed");
      }

      // HERE THE CORRELATIONS (NOT FLOW) are ready!!!
      // TODO: Rename that stuff!!!

    } // endfor {binPt}

    pCorTwoDif->SetName(nameCorTwo.Data());

    Debug("Calculating flow","ProcessDirect");
    // loading reference vn{2}
    TH1D* hFlowRefTwo = (TH1D*) listRefTwo->FindObject(Form("Refs_hFlow2_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
    if(!hFlowRefTwo) { Error(Form("Histo 'hFlowRefTwo' (sample %d) does not exists",iSample),"ProcessReconstructed"); listRefTwo->ls(); ffDesampleFile->ls(); return kFALSE; }

    // dn{2}
    hCumTwoDif = CalcDifCumTwo(pCorTwoDif, task);
    if(!hCumTwoDif) { Error(Form("dn{2} (sample %d) not processed correctly!",iSample),"ProcessReconstructed"); return kFALSE; }
    hCumTwoDif->SetName(Form("%s", nameCumTwo.Data()));

    // v'n{2}
    hFlowTwoDif = CalcDifFlowTwo(hCumTwoDif, hFlowRefTwo, iMultBin+1, task, bCorrelated);
    if(!hFlowTwoDif) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
    hFlowTwoDif->SetName(Form("%s", nameFlowTwo.Data()));

    if(task->fCumOrderMax >= 4)
    {
      pCorFourDif->SetName(nameCorFour.Data());
      // loading reference <<2>>
      TProfile* pCorTwoRef = (TProfile*) listRefCorTwo->FindObject(Form("Refs_pCor2_harm%d_gap%s_sample%d_rebin",task->fHarmonics, task->GetEtaGapString().Data(), iSample));
      if(!pCorTwoRef) { Error(Form("Profile 'pCorTwoRef' (sample %d) does not exists",iSample),"ProcessDirect"); listRefCorTwo->ls(); return kFALSE; }

      // loading reference vn{4}
      TH1D* hFlowRefFour = (TH1D*) listRefFour->FindObject(Form("Refs_hFlow4_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
      if(!hFlowRefFour) { Error(Form("Histo 'hFlowRefFour' (sample %d) does not exists",iSample),"ProcessDirect"); listRefFour->ls(); return kFALSE; }

      // dn{4}
      hCumFourDif = CalcDifCumFour(pCorFourDif, pCorTwoDif, pCorTwoRef, iMultBin+1, task, bCorrelated);
      if(!hCumFourDif) { Error(Form("dn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hCumFourDif->SetName(Form("%s", nameCumFour.Data()));

      // v'n{4}
      hFlowFourDif = CalcDifFlowFour(hCumFourDif, hFlowRefFour, iMultBin+1, task, bCorrelated);
      if(!hFlowFourDif) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hFlowFourDif->SetName(Form("%s", nameFlowFour.Data()));
    }


    ffOutputFile->cd();

    if(fSaveInterSteps) {
      pCorTwoDif->Write();
      hCumTwoDif->Write();
    }

    hFlowTwoDif->Write();

    if(task->fCumOrderMax >= 4) {
      if(fSaveInterSteps) {
        pCorFourDif->Write();
        hCumFourDif->Write();
      }
      hFlowFourDif->Write();
    }

    // if(fFlowFitCumulants)
    // {
    //   TH1D* hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
    //   if(!hRefFlow) { Error("Something went wrong when running automatic refs flow task:","ProcessReconstructed"); return kFALSE; }
    //
    //   TH1D* hFlow_vn = CalcDifFlowTwo(hFlow, hRefFlow, iMultBin+1 ,task, task->fConsCorr);
    //   hFlow_vn->SetName(Form("hFlow2_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
    //   hFlow_vn->SetTitle(Form("%s: v_{%d}{2,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); v_{%d}{2,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap));
    //   hFlow_vn->Write();
    //
    //   if(task->fCumOrderMax >= 4)
    //   {
    //     TH1D* hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
    //     if(!hRefFlow) { Error("Something went wrong when running automatic refs flow task:","ProcessReconstructed"); return kFALSE; }
    //     TH1D* hFlow_vn = CalcDifFlowFour(hFlowFour, hRefFlow, iMultBin+1, task, task->fConsCorr);
    //     hFlow_vn->SetName(Form("hFlow4_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
    //     hFlow_vn->SetTitle(Form("%s: v_{%d}{4,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); v_{%d}{4,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap));
    //     hFlow_vn->Write();
    //   }
    // }

    TCanvas* cFlow = new TCanvas("cFlow","cFlow");
    cFlow->cd();
    hFlowTwoDif->SetStats(0);
    hFlowTwoDif->Draw();
    cFlow->SaveAs(Form("%s/%s_Flow_n%d2_gap%s_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

    if(task->fCumOrderMax >= 4)
    {
      TCanvas* cFlow = new TCanvas("cFlow","cFlow");
      cFlow->cd();
      hFlowFourDif->SetStats(0);
      hFlowFourDif->Draw();
      cFlow->SaveAs(Form("%s/%s_FlowFour_n%d2_gap%s_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
    }
  }


  return kTRUE;
}
//_____________________________________________________________________________
TList* ProcessUniFlow::LoadSamples(TList* list, TString sHistName, Int_t iNumSamples)
{
    if(!list) {
        Error("Input list not found!", "LoadSamples");
        return nullptr;
    }

    TList* outList = new TList();

    for(Int_t iSample(0); iSample < iNumSamples; ++iSample) {
        TString sSampleName = sHistName + Form("_sample%d",iSample);
        TH1* hist = (TH1*) list->FindObject(sSampleName.Data());
        if(!hist) {
            Error(Form("Sample '%s' not found in input list!", sSampleName.Data()),"LoadSamples");
            list->ls();
            delete outList;
            return nullptr;
        }
        outList->Add(hist);
    }

    if(outList->GetEntries() != iNumSamples) {
        Warning(Form("Unexpected number of samples (%d instead of %d) !!!",outList->GetEntries(),iNumSamples),"LoadSamples");
    }

    return outList;
}
//_____________________________________________________________________________
TH1* ProcessUniFlow::MergeListProfiles(TList* list)
{
  // merge list of TProfiles into single TProfile and return it
  if(!list || list->IsEmpty()) { Error("List not valid or empty","MergeListProfiles"); return nullptr; }

  TH1* merged = (TH1*) list->At(0)->Clone();
  // merged->SetName(Form("%s_merged",merged->GetName()));

  if(list->GetEntries() < 2) // only 1 entry
  {
    Warning("Only one entry for merging; returning it directly instead!","MergeListProfiles");
    return merged;
  }

  merged->Reset();
  Double_t mergeStatus = merged->Merge(list);
  if(mergeStatus == -1) { Error("Merging failed!","MergeListProfiles"); return nullptr; }

  return merged;
}
//_____________________________________________________________________________
TH1* ProcessUniFlow::Merge(TH1* a, TH1* b)
{
    if(!a) { Error("Histogram (first) not found","Merge"); return nullptr; }
    if(!b) { Error("Histogram (second) not found","Merge"); return nullptr; }

    TList* list = new TList();
    list->Add(a);
    list->Add(b);

    TH1* merged = MergeListProfiles(list);
    delete list;
    return merged;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::DesampleList(TList* list, TH1D* merged, FlowTask* task, TString name, Bool_t bSkipDesampling)
{
  // bSkipDesampling [kFALSE] is used to keep the output structure cosistent (for check) even if Desampling can be skipped
  // e.g. vn{2} ~ cn{2} = <<2>>

  if(!merged) { Error("Merged histogram not valid","DesampleList"); return nullptr; }
  if(!list) { Error("List does not valid","DesampleList"); return nullptr; }
  if(list->GetEntries() < 1) { Error("List is empty","DesampleList"); return nullptr; }
  if(list->GetEntries() != task->fNumSamples) { Warning("Number of list entries is different from task number of samples","DesampleList"); }
  if(!task) { Error("FlowTask does not exists","DesampleList"); return nullptr; }

  Debug(Form("Number of samples in list pre-desampling: %d",list->GetEntries()),"DesampleList");

  TH1D* hDesampled = (TH1D*) merged->Clone(Form("%s_desampled",name.Data()));
  if(!hDesampled) { Error("Histo 'hDesampled' cloning failed","DesampleList"); return nullptr; }

  // saving objects to separate desampling file
  ffDesampleFile->cd();
  list->SetName(Form("%s_list",name.Data()));
  list->Write(0,TObject::kSingleKey);
  merged->SetName(Form("%s_merged",name.Data()));
  merged->Write();

  if(bSkipDesampling || task->fNumSamples < 2 || list->GetEntries() < 2)  // only one sample -> no sampling needed
  {
    // Warning("Only 1 sample for desampling; returning merged instead!","DesampleList");
    ffDesampleFile->cd();
    hDesampled->Write();
    return hDesampled;
  }

  for(Int_t iBin(0); iBin < hDesampled->GetNbinsX()+2; ++iBin)
  {
    const Double_t dDesMean = hDesampled->GetBinContent(iBin);

    Double_t dSum = 0.0;
    Int_t iCount = 0;

    for(Short_t iSample(0); iSample < list->GetEntries(); ++iSample)
    {
      TH1D* hTemp = (TH1D*) list->At(iSample);
      if(!hTemp) { Error(Form("Histo 'hTemp' (bin %d, sample %d) not found in list",iBin,iSample),"DesampleList"); return nullptr; }

      Double_t dContent = hTemp->GetBinContent(iBin);
      Double_t dError = hTemp->GetBinError(iBin);

      // TODO Check content & error

      dSum += TMath::Power(dDesMean - dContent, 2.0);
      iCount++;
    } // end-for {iSample} : loop over samples in list

    Double_t dError = TMath::Sqrt(dSum / (iCount*iCount));
    hDesampled->SetBinError(iBin,dError);

  } // end-for {iBin} : loop over bins in histo

  PlotDesamplingQA(list, hDesampled, task);

  hDesampled->Write();
  return hDesampled;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PlotDesamplingQA(TList* list, TH1D* hDesampled, FlowTask* task)
{
  if(!list) { Error("Input samples list not found!","PlotDesamplingQA"); return kFALSE; }
  if(!hDesampled) { Error("Desampled result not found!","PlotDesamplingQA"); return kFALSE; }
  if(!task) { Error("FlowTask not found!","PlotDesamplingQA"); return kFALSE; }

  Warning("Not implemented yet!","PlotDesamplingQA");
    //
    //
    // // getting copy which does not affect histo which is returned
    // TH1D* hDesampledClone = (TH1D*) hDesampled->Clone(Form("%sClone",hDesampled->GetName()));
    //
    // TList* listOutput = new TList(); // list for collecting all QA histos
    //
    // // doing QA plots with spread, etc.
    // TCanvas* canDesample = new TCanvas(Form("canDesample_%s",task->fName.Data()),Form("canDesample_%s",task->fName.Data()),1200,400);
    // canDesample->Divide(3,1);
    //
    // TH1D* hTempRatio = nullptr;
    // TH1D* hTempError = nullptr;
    //
    // TLine* lineUnity = new TLine();
    // lineUnity->SetLineColor(kRed);
    // lineUnity->SetLineWidth(3);
    //
    // canDesample->cd(1);
    // hDesampledClone->SetStats(kFALSE);
    // hDesampledClone->SetFillColor(kBlue);
    // hDesampledClone->SetStats(kFALSE);
    // hDesampledClone->SetMarkerStyle(20);
    // hDesampledClone->SetMarkerSize(0.5);
    // hDesampledClone->SetMarkerColor(kRed);
    // hDesampledClone->DrawCopy("E2");
    //
    // for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
    // {
    //   hTempSample = (TH1D*) list->At(iSample);
    //   if(!hTempSample) { Warning(Form("Sample %d not found during plotting QA! Skipping!",iSample),"DesampleList"); continue; }
    //
    //   canDesample->cd(1);
    //   hTempSample->SetStats(kFALSE);
    //   hTempSample->SetLineColor(30+2*iSample);
    //   hTempSample->SetMarkerColor(30+2*iSample);
    //   hTempSample->SetMarkerStyle(24);
    //   hTempSample->SetMarkerSize(0.5);
    //   hTempSample->DrawCopy("hist p same");
    //
    //   hTempRatio = (TH1D*) hTempSample->Clone(Form("%s_ratio",hTempSample->GetName()));
    //   hTempRatio->Divide(hDesampled);
    //   hTempRatio->SetYTitle("Value: final / sample");
    //   hTempRatio->SetTitleOffset(1.2,"Y");
    //
    //   canDesample->cd(2);
    //   hTempRatio->SetMinimum(0.6);
    //   hTempRatio->SetMaximum(1.4);
    //   hTempRatio->Draw("hist p same");
    //
    //   hTempError = (TH1D*) hTempSample->Clone(Form("%s_error",hTempSample->GetName()));
    //   for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hTempSample->GetBinError(bin)); }
    //
    //   canDesample->cd(3);
    //   hTempError->SetMinimum(0.);
    //   hTempError->SetMaximum(1.5*hTempError->GetMaximum());
    //   hTempError->SetYTitle("Uncertainty");
    //   hTempError->SetTitleOffset(1.2,"Y");
    //
    //   hTempError->Draw("hist p same");
    //
    //   listOutput->Add(hTempSample);
    //   listOutput->Add(hTempRatio);
    //   listOutput->Add(hTempError);
    // }
    //
    // canDesample->cd(1);
    // hDesampledClone->DrawCopy("hist p same");
    //
    // canDesample->cd(2);
    // lineUnity->DrawLine(hTempRatio->GetXaxis()->GetXmin(),1,hTempRatio->GetXaxis()->GetXmax(),1);
    //
    // hTempError = (TH1D*) hDesampledClone->Clone(Form("%s_error",hDesampled->GetName()));
    // for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hDesampledClone->GetBinError(bin)); }
    // listOutput->Add(hTempError);
    //
    // canDesample->cd(3);
    // hTempError->Draw("hist p same");
    //
    // // saving QA plots
    // canDesample->SaveAs(Form("%s/Desampling_%s_harm%d_gap%g_cent%d_%s.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,task->fName.Data(),fsOutputFileFormat.Data()));
    //
    // Info("Saving desampling QA into output file","DesampleList");
    // ffDesampleFile->cd();
    // listOutput->Add(canDesample);
    // listOutput->Write(Form("Desampling_%s_cent%d_%s",GetSpeciesName(task->fSpecies).Data(),iMultBin,task->fName.Data()),TObject::kSingleKey);
    //
    // // deleting created stuff
    // delete listOutput;
    // // delete canDesample;
    // delete lineUnity;
    // // if(hTempSample) delete hTempSample;
    // // delete hTempRatio;
    // // delete hTempError;
    // // delete hDesampledClone;
    //
    //

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor, TH3D* h3Entries, TH3D* h3EntriesBG, TProfile3D* p3CorFour)
{

  if(!task) { Error("Input task not found!","PrepareSlices"); return kFALSE; }
  if(!h3Entries) { Error("Input hist with entries not found!","PrepareSlices"); return kFALSE; }
  if(!p3Cor) { Error("Input profile with correlations not found!","PrepareSlices"); return kFALSE; }
  if(task->fCumOrderMax >= 4 && !p3CorFour) { Error("Input profile with <<4'>> not found!","PrepareSlices"); return kFALSE; }
  if(multBin < 0 || multBin > fiNumMultBins) { Error("Wrong multiplicity bin index (not in range)!","PrepareSlices"); return kFALSE; }

  // cleaning the vectros with flow-mass and inv. mass plots
  if(task->fVecHistInvMass->size() > 0) { task->fVecHistInvMass->clear(); }
  if(task->fVecHistInvMassBG->size() > 0) { task->fVecHistInvMassBG->clear(); }
  if(task->fVecHistFlowMass->size() > 0) { task->fVecHistFlowMass->clear(); }
  if(task->fCumOrderMax >= 4 && task->fVecHistFlowMassFour->size() > 0) { task->fVecHistFlowMassFour->clear(); }

  const Short_t binMultLow = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin]);
  const Short_t binMultHigh = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin+1]) - 1;
  // printf("Mult: %g(%d) -  %g(%d)\n",fdMultBins[multBin],binMultLow,fdMultBins[multBin+1],binMultHigh);

  TH1D* hRefFlow = nullptr;
  if(!fFlowFitCumulants)
  {
    // loading reference flow, if not found, it will be prepared
    hRefFlow = (TH1D*) ffDesampleFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g_desampled",task->fHarmonics,10*task->fEtaGap));
    if(!hRefFlow)
    {
      Warning("Relevant Reference flow not found within output file.","PrepareSlices");
      ffOutputFile->ls();

      Info("Creating relevant reference flow task.","PrepareSlices");
      FlowTask* taskRef = new FlowTask(kRefs,"Ref");
      taskRef->SetHarmonics(task->fHarmonics);
      taskRef->SetEtaGap(task->fEtaGap);
      taskRef->SetNumSamples(task->fNumSamples);
      taskRef->SetInputTag(task->fInputTag);
      taskRef->DoCumOrderMax(task->fCumOrderMax);

      if(ProcessRefs(taskRef))
      {
        hRefFlow = (TH1D*) ffDesampleFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
        if(!hRefFlow) {  Error("Automated Refs task completed, but RefFlow not found!","PrepareSlices"); ffDesampleFile->ls(); return kFALSE; }
      }
      else { Error("Something went wrong when running automatic refs flow task:","PrepareSlices"); taskRef->PrintTask(); return kFALSE; }
    }
  }

  TH1D* hRefCor2 = nullptr;
  TH1D* hRefFlow4 = nullptr;
  if(task->fCumOrderMax >= 4)
  {
    // loading <<2>> refs for dn{4}
    const char* name = Form("pCor2_Refs_harm%d_gap%02.2g_desampled",task->fHarmonics,10*task->fEtaGap);
    hRefCor2 = (TH1D*) ffDesampleFile->Get(name);
    if(!hRefCor2) { Error(Form("Input '%s' (hRefCor2) not found!",name),"PrepareSlices"); ffDesampleFile->ls(); return kFALSE; }

    const char* nameFour = Form("hFlow4_Refs_harm%d_gap%02.2g_desampled",task->fHarmonics,10*task->fEtaGap);
    hRefFlow4 = (TH1D*) ffDesampleFile->Get(nameFour);
    if(!hRefFlow4) { Error(Form("Input '%s' (hRefFlow4) not found!",nameFour),"PrepareSlices"); ffDesampleFile->ls(); return kFALSE; }
  }

  // loop over pt
  Short_t binPtLow = 0;
  Short_t binPtHigh = 0;
  TH1D* hInvMass_temp = nullptr;
  TH1D* hInvMassBG_temp = nullptr;
  TProfile3D* prof3Flow_temp = nullptr;
  TProfile2D* prof2FlowMass_temp = nullptr;
  TProfile* profFlowMass_temp = nullptr;
  TH1D* hFlowMass_temp = nullptr;

  prof3Flow_temp = (TProfile3D*) p3Cor->Clone(Form("prof3Flow_temp_cent%d",multBin));
  prof3Flow_temp->GetXaxis()->SetRange(binMultLow,binMultHigh);
  // prof2FlowMass_temp = (TProfile2D*) prof3Flow_temp->Project3DProfile("yz"); // NOTE: standard ROOT way - working properly in ROOTv6-12 onwards
  prof2FlowMass_temp = Project3DProfile(prof3Flow_temp);

  TProfile3D* prof3FlowFour_temp = nullptr;
  TProfile2D* prof2FlowMassFour_temp = nullptr;
  TProfile* profFlowMassFour_temp = nullptr;
  TH1D* hFlowMassFour_temp = nullptr;
  if(task->fCumOrderMax >= 4)
  {
    prof3FlowFour_temp = (TProfile3D*) p3CorFour->Clone(Form("prof3FlowFour_temp_cent%d",multBin));
    prof3FlowFour_temp->GetXaxis()->SetRange(binMultLow,binMultHigh);
    // prof2FlowMassFour_temp = (TProfile2D*) prof3FlowFour_temp->Project3DProfile("yz"); // NOTE: standard ROOT way - working properly in ROOTv6-12 onwards
    prof2FlowMassFour_temp = Project3DProfile(prof3FlowFour_temp);
  }

  Short_t iNumPtBins = task->fNumPtBins;

  TCanvas* canInvMass = new TCanvas("canInvMass","InvMass",1400,600);
  TCanvas* canInvMassBG = new TCanvas("canInvMassBG","InvMassBG",1400,600);
  TCanvas* canFlowMass = new TCanvas("canFlowMass","FlowMass",1400,600);
  TCanvas* canFlowMassFour = new TCanvas("canFlowMassFour","FlowMassFour",1400,600);
  canInvMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canInvMassBG->Divide(5,std::ceil(iNumPtBins/5)+1);
  canFlowMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canFlowMassFour->Divide(5,std::ceil(iNumPtBins/5)+1);

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

    task->fVecHistInvMass->push_back(hInvMass_temp);
    if(h3EntriesBG) { task->fVecHistInvMassBG->push_back(hInvMassBG_temp); }
    // done with InvMass

    // ### projection of flow-mass profile
    profFlowMass_temp = (TProfile*) prof2FlowMass_temp->ProfileX(Form("profFlowMass_cent%d_pt%d",multBin,binPt),binPtLow,binPtHigh);
    if(task->fCumOrderMax >= 4) { profFlowMassFour_temp = (TProfile*) prof2FlowMassFour_temp->ProfileX(Form("profFlowMassFour_cent%d_pt%d",multBin,binPt),binPtLow,binPtHigh); }

    // checking for rebinning the flow-mass profile
    if(task->fRebinFlowMass > 1) {
      profFlowMass_temp->Rebin(task->fRebinFlowMass);
      if(task->fCumOrderMax >= 4) { profFlowMassFour_temp->Rebin(task->fRebinFlowMass); }
    }

    // prepare slices for <<2>>
    hFlowMass_temp = CalcDifCumTwo(profFlowMass_temp,task);
    if(!hFlowMass_temp) { Error("<<2'>> not ready!","PrepareSlices"); return kFALSE; }
    hFlowMass_temp->SetName(Form("hFlowMass_cent%d_pt%d",multBin,binPt));
    if(!fFlowFitCumulants) { hFlowMass_temp = CalcDifFlowTwo(hFlowMass_temp, hRefFlow, multBin+1, task, task->fConsCorr); }
    if(!hFlowMass_temp) { Error("hFlowMass_temp not ready! Something went wrong!","PrepareSlices"); return kFALSE; }
    task->fVecHistFlowMass->push_back(hFlowMass_temp);

    // prepare slices for <<4>>
    if(task->fCumOrderMax >= 4) {
      hFlowMassFour_temp = CalcDifCumFour(profFlowMassFour_temp, profFlowMass_temp, hRefCor2, multBin+1, task, task->fConsCorr);
      if(!hFlowMassFour_temp) { Error("<<4'>> not ready!","PrepareSlices"); return kFALSE; }
      hFlowMassFour_temp->SetName(Form("hFlowMassFour_cent%d_pt%d",multBin,binPt));
      if(!fFlowFitCumulants) { hFlowMassFour_temp = CalcDifFlowFour(hFlowMassFour_temp, hRefFlow4, multBin+1, task, task->fConsCorr); }
      if(!hFlowMassFour_temp) { Error("FlowMassFour_temp not ready! Something went wrong!","PrepareSlices"); return kFALSE; }
      task->fVecHistFlowMassFour->push_back(hFlowMassFour_temp);
    }

    canInvMass->cd(binPt+1);
    hInvMass_temp->Draw();

    if(h3EntriesBG)
    {
      canInvMassBG->cd(binPt+1);
      hInvMassBG_temp->Draw();
    }

    canFlowMass->cd(binPt+1);
    hFlowMass_temp->Draw();

    if(task->fCumOrderMax >= 4)
    {
      canFlowMassFour->cd(binPt+1);
      hFlowMassFour_temp->Draw();
    }

  } // endfor {binPt}: over Pt bins

  printf(" # of slices: InvMass: %lu | InvMassBG %lu | FlowMass %lu\n",task->fVecHistInvMass->size(),task->fVecHistInvMassBG->size(),task->fVecHistFlowMass->size());

  gSystem->mkdir(Form("%s/slices/",fsOutputFilePath.Data()));
  canInvMass->SaveAs(Form("%s/slices/Slices_InvMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  if(h3EntriesBG) canInvMassBG->SaveAs(Form("%s/slices/Slices_InvMassBG_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  canFlowMass->SaveAs(Form("%s/slices/Slices_FlowMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  if(task->fCumOrderMax >= 4) { canFlowMassFour->SaveAs(Form("%s/slices/Slices_FlowMassFour_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),GetSpeciesName(task->fSpecies).Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data())); }

  if(task->fVecHistInvMass->size() < 1 || task->fVecHistFlowMass->size() < 1 || task->fVecHistFlowMass->size() != task->fVecHistInvMass->size()) { Error("Output vector empty. Something went wrong","PrepareSlices"); return kFALSE; }
  if(h3EntriesBG && (task->fVecHistInvMassBG->size() < 1 || task->fVecHistInvMassBG->size() != task->fVecHistInvMass->size()) ) { Error("Output vector empty. Something went wrong with BG histograms","PrepareSlices"); return kFALSE; }
  if(task->fCumOrderMax >= 4 && (task->fVecHistFlowMassFour->size() < 1 || task->fVecHistFlowMass->size() != task->fVecHistInvMass->size()) ) { Error("Output vector for <<4>> empty. Something went wrong","PrepareSlices"); return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PrepareSlicesNew(FlowTask* task, TString histName, Bool_t bDoCand)
{
  // wrapper for making/preparing per-task slices
  if(!task) { Error("FlowTask does not exists!","PrepareSlicesNew"); return kFALSE; }

  PartSpecies species = task->fSpecies;
  Bool_t bReco = IsSpeciesReconstructed(species);

  TList* inputList = flFlow[species];
  if(!inputList) { Error("Input TList not found!","PrepareSlicesNew"); return kFALSE; }

  // preparing flow slices
  TH1* prof = nullptr;

  TList* listPos = LoadSamples(inputList, Form("%s_Pos",histName.Data()), task->fNumSamples);
  TH1* profPos = (TH1*) MergeListProfiles(listPos);
  if(!profPos) { Error(Form("Positive profile '%s_Neg' not found!",histName.Data()),"PrepareSlicesNew"); return kFALSE; }
  prof = profPos;

  if(task->fMergePosNeg) {
    TList* listNeg = LoadSamples(inputList, Form("%s_Neg",histName.Data()), task->fNumSamples);
    TH1* profNeg = (TH1*) MergeListProfiles(listNeg);
    if(!profNeg) { Error(Form("Negative profile '%s_Neg' not found!",histName.Data()),"PrepareSlicesNew"); return kFALSE; }

    prof = (TH1*) Merge(profPos, profNeg);
  }

  if(!prof) { Error(Form("Profile '%s' not found!",histName.Data()),"PrepareSlicesNew"); inputList->ls(); return kFALSE; }
  if(!MakeProfileSlices(task,prof,task->fListProfiles)) { Error("Profile Slices failed!","PrepareSlicesNew"); return kFALSE; };

  // preparing inv. mass slices (NB: merging pos/neg done in MakeSparseSlices() )
  if(bReco && bDoCand) {

    TString sNameCand = Form("fhsCand%s",GetSpeciesName(task->fSpecies).Data());

    THnSparseD* sparse = (THnSparseD*) inputList->FindObject(sNameCand.Data());
    if(!MakeSparseSlices(task,sparse,task->fListHistos)) { Error("Histo Slices failed!","PrepareSlicesNew"); return kFALSE; };

    if(species == kPhi) {
      THnSparseD* sparseBg = (THnSparseD*) inputList->FindObject(Form("%sBg",sNameCand.Data()));
      if(!MakeSparseSlices(task,sparseBg,task->fListHistos,"hInvMassBg")) { Error("Histo Slices for Phi BG failed!","PrepareSlicesNew"); return kFALSE; };
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t ProcessUniFlow::MakeProfileSlices(FlowTask* task, TH1* inputProf, TList* outList)
{
  // prepare slices out of inputHist
  if(!task) { Error("FlowTask does not exists!","MakeProfileSlices"); return kFALSE; }
  if(!inputProf) { Error("Input profile does not exists!","MakeProfileSlices"); return kFALSE; }
  if(!outList) { Error("Output TList does not exists!","MakeProfileSlices"); return kFALSE; }

  PartSpecies spec = task->fSpecies;
  if(spec == kRefs) { Error("Species is 'kRefs': no slicing required!","MakeProfileSlices"); return kFALSE; }

  Bool_t bReco = IsSpeciesReconstructed(spec);

  Int_t iNumBinsMult = fiNumMultBins;
  Int_t iNumBinsPt = task->fNumPtBins;

  TAxis* axisMult = inputProf->GetXaxis();
  TAxis* axisPt = inputProf->GetYaxis();

  TList trashCol;
  trashCol.SetOwner(kTRUE);

  for(Int_t iBinMult(0); iBinMult < iNumBinsMult; ++iBinMult) {
    const Int_t iBinMultLow = axisMult->FindFixBin(fdMultBins[iBinMult]);
    const Int_t iBinMultHigh = axisMult->FindFixBin(fdMultBins[iBinMult+1]) - 1;

    if(!bReco) {
      // direct species
      TProfile* prof1D_preRebin = ((TProfile2D*)inputProf)->ProfileY("",iBinMultLow,iBinMultHigh);
      if(!prof1D_preRebin) { Error("Profile 'prof1D_preRebin' failed!","MakeProfileSlices"); return kFALSE; }
      trashCol.Add(prof1D_preRebin); // to ensure to be deleted

      TProfile* prof1D = nullptr;
      if(iNumBinsPt > 0) { prof1D = (TProfile*) prof1D_preRebin->Rebin(iNumBinsPt,Form("%s_rebin", inputProf->GetName()), task->fPtBinsEdges.data()); }
      else { prof1D = (TProfile*) prof1D_preRebin->Clone(); }
      if(!prof1D) { Error("Profile 'prof1D' does not exists!","MakeProfileSlices"); return kFALSE; }

      prof1D->SetName(Form("%s_mult%d",inputProf->GetName(),iBinMult));
      prof1D->GetXaxis()->SetTitle(inputProf->GetYaxis()->GetTitle());
      outList->Add(prof1D);
    } else {
      // reconstructed species
      axisMult->SetRange(iBinMultLow,iBinMultHigh);
      TProfile2D* prof2D = Project3DProfile((TProfile3D*) inputProf);
      if(!prof2D) { Error("Mult projection failed!","MakeProfileSlices"); return kFALSE; }
      trashCol.Add(prof2D); // NB: to ensure that it will be deleted

      for(Int_t iBinPt(0); iBinPt < iNumBinsPt; ++iBinPt) {
        const Double_t dEdgePtLow = task->fPtBinsEdges[iBinPt];
        const Double_t dEdgePtHigh = task->fPtBinsEdges[iBinPt+1];

        const Int_t iBinPtLow = axisPt->FindFixBin(dEdgePtLow);
        const Int_t iBinPtHigh = axisPt->FindFixBin(dEdgePtHigh) - 1;

        TProfile* prof1D = prof2D->ProfileX("",iBinPtLow,iBinPtHigh);
        if(!prof1D) { Error("Profile 'prof1D' failed!","MakeProfileSlices"); return kFALSE; }

        if(task->fRebinFlowMass > 1) {
          prof1D->Rebin(task->fRebinFlowMass);
        }

        prof1D->SetName(Form("%s_mult%d_pt%d",inputProf->GetName(),iBinMult,iBinPt));
        prof1D->GetXaxis()->SetTitle(inputProf->GetZaxis()->GetTitle());
        outList->Add(prof1D);
      } // end-for {binPt}
    } // end-else {!bReco}
  } // end-for {binMult}

  Info("Successfull!","MakeProfileSlices");
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::MakeSparseSlices(FlowTask* task, THnSparse* inputSparse, TList* outList, const char* outName)
{
  // prepare slices out of inputHist
  if(!task) { Error("FlowTask does not exists!","MakeSparseSlices"); return kFALSE; }
  if(!inputSparse) { Error("Input THnSparse does not exists!","MakeSparseSlices"); return kFALSE; }
  if(!outList) { Error("Output TList does not exists!","MakeSparseSlices"); return kFALSE; }
  // if(outList->GetEntries() > 0) { Error("Output TList is not empty!","MakeSparseSlices"); return kFALSE; }

  PartSpecies species = task->fSpecies;
  if(!IsSpeciesReconstructed(species)) {
    Error("Not a reconstructed species!","MakeSparseSlices");
    return kFALSE;
  }

  TList trashCol;
  trashCol.SetOwner(kTRUE);

  Double_t dEtaGap = task->fEtaGap;
  Bool_t bHasGap = task->HasGap();

  TH3D* histEntries = nullptr;
  TH3D* histEntriesPos = nullptr;
  TH3D* histEntriesNeg = nullptr;

  if(!bHasGap) {
    histEntries = (TH3D*) inputSparse->Projection(1,2,0);
    trashCol.Add(histEntries);
  } // end-if {!bHasGap}
  else {
    TAxis* axisEta = inputSparse->GetAxis(3);
    if(!axisEta) { Error("TAxis 'axisEta' does not exists!","MakeSparseSlices"); return kFALSE; }

    // positive POIs
    axisEta->SetRangeUser(dEtaGap/2.0,axisEta->GetXmax());
    TH3D* histEntriesPos = (TH3D*) inputSparse->Projection(1,2,0);
    if(!histEntriesPos) { Error("Projection 'histEntriesPos' failed!","MakeSparseSlices"); return kFALSE; }
    trashCol.Add(histEntriesPos);

    // negative POIs
    axisEta->SetRangeUser(axisEta->GetXmin(),-dEtaGap/2.0);
    TH3D* histEntriesNeg = (TH3D*) inputSparse->Projection(1,2,0);
    if(!histEntriesNeg) { Error("Projection 'histEntriesNeg' failed!","MakeSparseSlices"); return kFALSE; }
    trashCol.Add(histEntriesNeg);

    if(task->fMergePosNeg) {
      TList* listMerge = new TList();
      // No need to add ownership since it is collected by trashCol TList
      listMerge->Add(histEntriesPos);
      listMerge->Add(histEntriesNeg);
      histEntries = (TH3D*) MergeListProfiles(listMerge);
      delete listMerge;
      trashCol.Add(histEntries);
    } // end-if {fMergePosNeg}
    else {
      // loading single histo (positive by default)
      if(task->fInputTag.EqualTo("")) {
        histEntries = histEntriesPos;
      }
      else if (task->fInputTag.EqualTo("Neg")) {
        histEntries = histEntriesNeg;
      }
      else {
        Error(Form("Invalid InputTag '%s'!",task->fInputTag.Data()),"ProcessReconstructed");
        return kFALSE;
      }
    } // end-else {fMergePosNeg}
  } // end-else {!bHasGap}

  if(!histEntries) { Error("Histo 'histEntries' failed!","MakeSparseSlices"); return kFALSE; }
  // histEntries (TH3D*) ready for slicing in mult & pt bins (NB should be put in trashCol)

  Int_t iNumBinsMult = fiNumMultBins;
  Int_t iNumBinsPt = task->fNumPtBins;

  TAxis* axisMult = histEntries->GetXaxis();
  TAxis* axisPt = histEntries->GetYaxis();

  for(Int_t iBinMult(0); iBinMult < iNumBinsMult; ++iBinMult) {
    const Int_t iBinMultLow = axisMult->FindFixBin(fdMultBins[iBinMult]);
    const Int_t iBinMultHigh = axisMult->FindFixBin(fdMultBins[iBinMult+1]) - 1;

    for(Int_t iBinPt(0); iBinPt < iNumBinsPt; ++iBinPt) {
      const Double_t dEdgePtLow = task->fPtBinsEdges[iBinPt];
      const Double_t dEdgePtHigh = task->fPtBinsEdges[iBinPt+1];
      const Int_t iBinPtLow = axisPt->FindFixBin(dEdgePtLow);
      const Int_t iBinPtHigh = axisPt->FindFixBin(dEdgePtHigh) - 1;

      TH1D* histInvMass = (TH1D*) histEntries->ProjectionZ(Form("%s_mult%d_pt%d",outName,iBinMult,iBinPt),iBinMultLow,iBinMultHigh,iBinPtLow,iBinPtHigh,"e");
      if(!histInvMass) { Error("Projection 'histInvMass' failed!","MakeSparseSlices"); return kFALSE; }

      if(task->fRebinInvMass > 1) {
        histInvMass->Rebin(task->fRebinInvMass);
      }

      outList->Add(histInvMass);
    } // end-for {binPt}
  } // end-for {binMult}

  Debug("Successfull!","MakeSparseSlices");
  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::TestRebin(TH1D* hOrig, FlowTask* task)
{
  if(!hOrig) { Error("Original histogram not found!","TestRebin"); return nullptr; }
  if(!task) { Error("Task not found!","TestRebin"); return nullptr; }
  TH1D* hRebin = (TH1D*) hOrig->Rebin(fiNumMultBins,Form("%s_testRebin",hOrig->GetName()),fdMultBins.data());

  // rebinning
  Short_t numBins = fiNumMultBins;
  Double_t* multBins = fdMultBins.data();
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

  if(task->fSpecies == kRefs) fvTasks.insert(fvTasks.begin(), task);
  else fvTasks.push_back(task);

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
  if(!prof3dorig) { Error("Input profile does not exists!","Project3DProfile"); return nullptr; }

  TProfile3D* prof3d = (TProfile3D*) prof3dorig->Clone();
  if(!prof3d) { Error("Cloning failed!","Project3DProfile"); return nullptr; }

  Int_t iBinFirst = prof3d->GetXaxis()->GetFirst();
  Int_t iBinLast = prof3d->GetXaxis()->GetLast();
  Int_t iNumBins = prof3d->GetNbinsX();
  Int_t iNumBinsAxis = prof3d->GetXaxis()->GetNbins();
  // printf("Bins:  %d - %d (%d | %d) \n", iBinFirst,iBinLast,iNumBins,iNumBinsAxis);

  // // making 3d hist from 3d profile
  // TH3D* hist3d = prof3d->ProjectionXYZ();   //NOTE do not care about range !!!
  // TH3D* hist3d_entry = prof3d->ProjectionXYZ("hist3d_entry","B");   //NOTE do not care about range !!!
  // TH3D* hist3d_weight = prof3d->ProjectionXYZ("hist3d_weight","W");   //NOTE do not care about range !!!

  TProfile2D* prof2d_test = DoProjectProfile2D(prof3d,Form("%s_px",prof3d->GetName()),prof3d->GetTitle(),prof3d->GetYaxis(),prof3d->GetZaxis(),1,0,0);
  if(!prof2d_test) { Error("DoProjectProfile2D failed!","Project3DProfile"); delete prof3d; return nullptr; }
  prof2d_test->GetXaxis()->SetTitle(prof3d->GetZaxis()->GetTitle());
  prof2d_test->GetYaxis()->SetTitle(prof3d->GetYaxis()->GetTitle());

  // prof2d_test->Draw("colz");
  delete prof3d;
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
 bool useWeights = (h3->GetBinSumw2()->fN != 0); //array elements
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

     Double_t stats[TH1::kNstat] = {0};
     h3->GetStats(stats);
     double dfTsumw = stats[0];

     if (h3->IsA() == TH3F::Class() ) eps = 1.E-6;
     if (dfTsumw != 0 && TMath::Abs( dfTsumw - totcont) <  TMath::Abs(dfTsumw) * eps) resetStats = false;

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
        // since this  is the only way to CalcCum them
        Double_t entries =  h2->GetEffectiveEntries();
        if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
        h2->SetEntries( entries );
     }
     else {
        h2->SetEntries( h3->GetEntries() );
     }


     return h2;
}
//_____________________________________________________________________________
void ProcessUniFlow::PrintFitFunction(const TF1* func)
{
    if(!func) { Error("Input function not found!","PrintFitFunction"); return; }
    printf("Name '%s'\n",func->GetName());
    printf("Title '%s'\n", func->GetTitle());
    // printf("Status '%s'\n",func->fCstatus.Data());
    for(Int_t iPar(0); iPar < func->GetNpar(); ++iPar) {
        Double_t dLimLow, dLimHigh;
        func->GetParLimits(iPar,dLimLow,dLimHigh);
        printf("  p%d  \t'%s' :    \t%g \t+- \t%g   \t(%g <=> %g)\n", iPar, func->GetParName(iPar), func->GetParameter(iPar), func->GetParError(iPar), dLimLow, dLimHigh);
    }
    printf("Chi2/ndf = %.3g/%d = %.3g\n", func->GetChisquare(), func->GetNDF(),func->GetChisquare()/func->GetNDF());
    printf("p-value = %g\n",func->GetProb());
    printf("##############################################################\n");
}
//_____________________________________________________________________________
TH1* ProcessUniFlow::SubtractInvMassBg(TH1* hInvMass, TH1* hInvMassBg, FlowTask* task)
{
    if(!hInvMass) { Error("Input inv. mass histo not found!","SubtractInvMassBg"); return nullptr; }
    if(!hInvMassBg) { Error("Input inv. mass Bg histo not found!","SubtractInvMassBg"); return nullptr; }
    if(!task) { Error("FlowTask not found!","SubtractInvMassBg"); return nullptr; }

    if(hInvMass->GetNbinsX() != hInvMassBg->GetNbinsX()) { Error("Different number of bins for signal & bg histo!","SubtractInvMassBg"); return nullptr; }

    Double_t dNorm = 1.0;

    // normalise
    if(task->fbNormLS) {
        Double_t dLow = task->fdNormLSLow;
        Double_t dHigh = task->fdNormLSHigh;

        Int_t binLow = hInvMass->FindBin(task->fdNormLSLow);
        Int_t binHigh = hInvMass->FindBin( task->fdNormLSHigh);

        if(binLow >= binHigh) { Error("Low norm. range is higher than high norm. range!","SubtractInvMassBg"); return nullptr; }
        if(binLow < 1) { Error("Low norm. range is lower than histo range!","SubtractInvMassBg"); return nullptr; }
        if(binHigh > hInvMass->GetNbinsX()) { Error("High norm. range is higher than histo range!","SubtractInvMassBg"); return nullptr; }
        if(binLow != (hInvMassBg->FindBin(task->fdNormLSLow)) ) { Error("Low norm. range bin is different in hInvMassBg histo!","SubtractInvMassBg"); return nullptr; }
        if(binHigh != (hInvMassBg->FindBin(task->fdNormLSHigh)) ) { Error("High norm. range bin is different in hInvMassBg histo!","SubtractInvMassBg"); return nullptr; }

        Double_t dIntegral = hInvMass->Integral(binLow,binHigh);
        Double_t dIntegralBg = hInvMassBg->Integral(binLow,binHigh);


        if(dIntegral > 0 && dIntegralBg > 0) {
            dNorm = dIntegral / dIntegralBg;
        } else {
            Error("Either of the integrals is 0!","SubtractInvMassBg");
            return nullptr;
        }

        Debug(Form("dNorm = (dIntegral = %g / dIntegralBg = %g) = %g",dIntegral,dIntegralBg,dNorm),"SubtractInvMassBg");
    }

    TH1* hInvMassSubt = (TH1*) hInvMass->Clone();
    hInvMassSubt->Sumw2();
    hInvMassSubt->Add(hInvMassBg, -1.0*dNorm);

    // for(Int_t iBin(0); iBin < iNumBins+1; ++iBin) {
    //     Double_t dContSig = hInvMass->GetBinContent(iBin);
    //     Double_t dContBg = hInvMassBg->GetBinContent(iBin);
    //     Double_t dContNew = dContSig - dContBg;
    //
    //     hInvMassSubt->SetBinContent(iBin, dContNew);
    //     hInvMassSubt->SetBinError(iBin, TMath::Sqrt(dContSig+dContBg));
    // }

    return hInvMassSubt;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::SetFuncParameters(TF1* func, Double_t* dVal, const std::vector<Double_t>& vecLow, const std::vector<Double_t>& vecHigh, const std::vector<TString> vecNames)
{
    if(!func) { Error("Input function not found!","SetFuncParameters"); return kFALSE; }
    Int_t iNumPar = func->GetNpar();
    std::vector<Double_t> dVec;

    for(Int_t iPar(0); iPar < iNumPar; ++iPar) {
        Double_t v = dVal[iPar];
        dVec.push_back(v);
        // printf("%d  :  %f    %f\n",iPar,v,dVec.at(iPar));
    }

    return SetFuncParameters(func, dVec, vecLow, vecHigh, vecNames);
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::SetFuncParameters(TF1* func, const std::vector<Double_t>& vecVal, const std::vector<Double_t>& vecLow, const std::vector<Double_t>& vecHigh,  const std::vector<TString> vecNames)
{
    if(!func) { Error("Input function not found!","SetFuncParameters"); return kFALSE; }

    Int_t iNumPar = func->GetNpar();
    Int_t iSizeVal = vecVal.size();
    Int_t iSizeLow = vecLow.size();
    Int_t iSizeHigh = vecHigh.size();
    Int_t iSizeNames = vecNames.size();

    if(iSizeVal != iSizeLow) { Error("Size of vector vecVal different from dLow!","SetFuncParameters"); return kFALSE; }
    if(iSizeVal != iSizeHigh) { Error("Size of vector vecVal different from dHigh!","SetFuncParameters"); return kFALSE; }
    if(iSizeVal < iNumPar) { Error("Size of vector vecVal smaller than number of func parameters!","SetFuncParameters"); return kFALSE; }
    if(iSizeNames && iSizeNames != iNumPar) { Error("Size of vector vecNames smaller than number of func parameters!","SetFuncParameters"); return kFALSE; }

    for(Int_t par(0); par < iNumPar; ++par) {
        Double_t dVal = vecVal.at(par);
        Double_t dLow = vecLow.at(par);
        Double_t dHigh = vecHigh.at(par);

        if(iSizeNames) { func->SetParName(par, vecNames.at(par).Data()); }

        if(dLow > -1.0 || dHigh > -1.0) {
            if(dLow > -1.0 && dHigh > -1.0) {
                // both limits == par val -> fix parameter
                // if(dLow == dHigh && dLow == dVal) {
                if(dLow == dHigh) {
                    printf("Fixing parameters %d : %f\n", par, dVal);
                    func->FixParameter(par,dVal);
                    printf("Done : %f\n", func->GetParameter(par));
                    continue;
                // both limits set and proper value
                } else if (dLow < dHigh && dLow <= dVal && dVal <= dHigh) {
                    func->SetParameter(par, dVal);
                    func->SetParLimits(par, dLow, dHigh);
                    continue;
                } else {
                    Error(Form("Badly setup parameter limits: parameter %d | dLow %g | dPar %g | dHigh %g", par, dLow, dVal, dHigh),"SetFuncParameters");
                    return kFALSE;
                }
            }
            // only one limit set
            Error(Form("Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!", par, dVal, dLow, dHigh),"SetFuncParameters");
            return kFALSE;
        } else {
            // no limits set, just set parameter
            func->SetParameter(par, dVal);
        }
    }

    return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::CheckFitResult(TFitResultPtr result, Bool_t bIgnorePOSDEF)
{
    Bool_t isOK = kTRUE;

    Int_t status = (Int_t) result;

    if(status != 0) {
        Warning(Form("Fit result status not zero (%d)!",status),"CheckFitResult");
        isOK = kFALSE;
    }

    if(!gMinuit) {
        Warning(Form("gMinuit not available! "),"CheckFitResult");
        isOK = kFALSE;
    }

    TString sMinuit = gMinuit->fCstatu;

    if(sMinuit.Contains("NOT POSDEF")) {
        if(!bIgnorePOSDEF) {
            Warning(Form("gMinuit status is '%s' while ignorePOSDEF is OFF! ",sMinuit.Data()),"CheckFitResult");
            isOK = kFALSE;
        } else {
            Warning(Form("gMinuit status is '%s'! Ignored.",sMinuit.Data()),"CheckFitResult");
        }
    } else if (!sMinuit.Contains("CONVERGED") && !sMinuit.Contains("OK") ) {
        Warning(Form("gMinuit status ('%s') does not converged! ",sMinuit.Data()),"CheckFitResult");
        isOK = kFALSE;
    }

    // NB: if Debug flag is ON, the fit statistics will be shown anyway
    if(!isOK && !fbDebug) {
        TVirtualFitter* fitter = TVirtualFitter::GetFitter();
        fitter->PrintResults(3,0.0);
    }

    return isOK;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::FitInvMass(TH1* hist, FlowTask* task, TF1& fitOut, TF1& fitOutSig, TF1& fitOutBg, TList* outList, TH1* histBg)
{
  if(!hist) { Error("Input histo not found!","FitInvMass"); return kFALSE; }
  if(!task) { Error("FlowTask not found!","FitInvMass"); return kFALSE; }
  if(!outList) { Error("Output TList outList not found!"); return kFALSE; }

  PartSpecies species = task->fSpecies;
  if(task->fFlowFitPhiSubtLS && species == kPhi && !histBg) { Error("Input histo bg not found!","FitInvMass"); return kFALSE; }
  if(!IsSpeciesReconstructed(species)) { Error("Invalid species!","FitInvMass"); return kFALSE; }

  Double_t dMassRangeLow = hist->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hist->GetXaxis()->GetXmax();
  Double_t dMaximum = hist->GetBinContent(hist->GetMaximumBin());

  hist->SetName("histMass");
  outList->Add(hist);

  Int_t iNpx = 10000;
  TString sFitOptMass = "RNBS";
  // TString sFitOptMass = "RNLB";

  TString sMassBG = TString(); Int_t iNumParsMassBG = 0; // function for inv. mass dist. (BG component)
  TString sMassSig = TString();  Int_t iNumParsMassSig = 0; // function for inv. mass dist. (sig component)

  Int_t iParMass = 0;
  Int_t iParMass_2 = 0;
  Int_t iParWidth = 0;
  Int_t iParWidth_2 = 0;

  Double_t dPeakLow = 0.0;
  Double_t dPeakHigh = 0.0;

  Double_t dFracLimLow = 0.0; // shift of fix fraction
  Double_t dFracLimHigh = 0.0;

  std::vector<Double_t> dParDef;
  std::vector<Double_t> dParLimLow;
  std::vector<Double_t> dParLimHigh;
  std::vector<TString> sParNames;

  if(species == kPhi) {
    Debug("Setting parameters for Phi","FitInvMass");

    dMassRangeLow = 1.00;
    dMassRangeHigh = 1.06;
    // dMassRangeLow = 0.994;
    // dMassRangeHigh = 1.134;

    dPeakLow = 1.007;
    dPeakHigh = 1.03;

    dFracLimLow = 0.05;
    dFracLimHigh = 0.05;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*(TMath::BreitWigner(x,[5],[6]))"; iNumParsMassSig = 3;

    iParMass = 5;
    iParWidth = 6;

    sParNames.push_back("bg0");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg1");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg2");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg3");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);

    sParNames.push_back("ampTot");      dParDef.push_back(1.2*dMaximum);    dParLimLow.push_back(0.0);          dParLimHigh.push_back(1.2*dMaximum);
    // sParNames.push_back("ampBW");       dParDef.push_back(1.0);         dParLimLow.push_back(1.0);          dParLimHigh.push_back(1.0);
    sParNames.push_back("mean");        dParDef.push_back(1.019445);    dParLimLow.push_back(1.0185);       dParLimHigh.push_back(1.021);
    sParNames.push_back("width");       dParDef.push_back(0.004);       dParLimLow.push_back(0.004);        dParLimHigh.push_back(0.009);
  }

  if(species == kK0s) {
    Debug("Setting parameters for K0s","FitInvMass");

    // dMassRangeLow = 0.994;
    // dMassRangeHigh = 1.134;

    dPeakLow = 0.44;
    dPeakHigh = 0.55;

    dFracLimLow = 0.05;
    dFracLimHigh = 0.07;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*([5]*TMath::Gaus(x,[6],[7])+(1.0-[5])*TMath::Gaus(x,[8],[9]))"; iNumParsMassSig = 6;

    iParMass = 6;
    iParMass_2 = 8;
    iParWidth = 7;
    iParWidth_2 = 9;

    sParNames.push_back("bg0");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg1");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg2");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg3");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);

    sParNames.push_back("ampTot");      dParDef.push_back(1.2*dMaximum);    dParLimLow.push_back(0.0);          dParLimHigh.push_back(1.2*dMaximum);
    sParNames.push_back("ampG1");       dParDef.push_back(0.55);         dParLimLow.push_back(0.55);          dParLimHigh.push_back(0.99);
    sParNames.push_back("meanG1");      dParDef.push_back(0.4976);      dParLimLow.push_back(0.48);         dParLimHigh.push_back(0.51);
    sParNames.push_back("sigmaG1");     dParDef.push_back(0.001);       dParLimLow.push_back(0.001);        dParLimHigh.push_back(0.05);
    // sParNames.push_back("ampG2");       dParDef.push_back(0.0);         dParLimLow.push_back(0.0);          dParLimHigh.push_back(0.4);
    sParNames.push_back("meanG2");      dParDef.push_back(0.48);      dParLimLow.push_back(0.48);         dParLimHigh.push_back(0.52);
    sParNames.push_back("sigmaG2");     dParDef.push_back(0.2);        dParLimLow.push_back(0.001);        dParLimHigh.push_back(0.2);
  }

  if(species == kLambda)
  {
    Debug("Setting parameters for Lambda","FitInvMass");

    dMassRangeLow = 1.099;
    // dMassRangeHigh = 1.15;

    dPeakLow = 1.106;
    dPeakHigh = 1.125;

    dFracLimLow = 0.05;
    dFracLimHigh = 0.07;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*([5]*TMath::Gaus(x,[6],[7])+(1.0-[5])*TMath::Gaus(x,[8],[9]))"; iNumParsMassSig = 6;

    iParMass = 6;
    iParMass_2 = 8;
    iParWidth = 7;
    iParWidth_2 = 9;

    sParNames.push_back("bg0");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg1");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg2");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);
    sParNames.push_back("bg3");         dParDef.push_back(0.0);         dParLimLow.push_back(-1);           dParLimHigh.push_back(-1);

    sParNames.push_back("ampTot");      dParDef.push_back(dMaximum);    dParLimLow.push_back(0.0);          dParLimHigh.push_back(1.2*dMaximum);
    sParNames.push_back("ampG1");       dParDef.push_back(0.6);         dParLimLow.push_back(0.6);          dParLimHigh.push_back(1.0);
    sParNames.push_back("meanG1");      dParDef.push_back(1.115);       dParLimLow.push_back(1.11);         dParLimHigh.push_back(1.12);
    sParNames.push_back("sigmaG1");     dParDef.push_back(0.001);       dParLimLow.push_back(0.001);        dParLimHigh.push_back(0.007);
    // sParNames.push_back("ampG2");       dParDef.push_back(0.2);         dParLimLow.push_back(0.0);          dParLimHigh.push_back(1.0);
    sParNames.push_back("meanG2");      dParDef.push_back(1.109);       dParLimLow.push_back(1.109);         dParLimHigh.push_back(1.125);
    sParNames.push_back("sigmaG2");     dParDef.push_back(0.001);        dParLimLow.push_back(0.001);        dParLimHigh.push_back(0.01);
  }

  // check if parametrisation is setup manually
  if(task->fFlowFitRangeLow > 0.0) { dMassRangeLow = task->fFlowFitRangeLow; }
  if(task->fFlowFitRangeHigh > 0.0) { dMassRangeHigh = task->fFlowFitRangeHigh; }
  if(task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow >= task->fFlowFitRangeHigh) { Error("Wrong fitting ranges set!","FitInvMass"); return kFALSE; }

  Bool_t bUserPars = kFALSE;
  if(task->fNumParMassSig > 0) { bUserPars = kTRUE; sMassSig = task->fFlowFitMassSig; iNumParsMassSig = task->fNumParMassSig; Debug(" Task massSig set","FitInvMass"); }
  if(task->fNumParMassBG > 0) { bUserPars = kTRUE; sMassBG = task->fFlowFitMassBG; iNumParsMassBG = task->fNumParMassBG; Debug(" Task massBG set","FitInvMass"); }
  if(bUserPars && (task->fNumParMassSig == 0 || task->fNumParMassBG == 0)) { Error("Only a subset of functions has been changed. Provide all, or non.","FitInvMass"); return kFALSE; }

  if(bUserPars) {
    Info("Setting UserParameters","FitInvMass");
    dParDef.clear();
    dParLimLow.clear();
    dParLimHigh.clear();

    Int_t iNumParTot = task->fNumParMassSig + task->fNumParMassBG;
    for(Int_t par(0); par < iNumParTot; ++par)
    {
      dParDef.push_back(task->fFitParDefaults[par]);
      dParLimLow.push_back(task->fFitParLimLow[par]);
      dParLimHigh.push_back(task->fFitParLimHigh[par]);
    }
  }

  Int_t iNumParTot = iNumParsMassSig+iNumParsMassBG;
  Int_t iNumParDefs = dParDef.size();
  Int_t iNumParLimLow = dParLimLow.size();
  Int_t iNumParLimHigh = dParLimHigh.size();

  if(iNumParDefs != iNumParTot) { Error(Form("Length of dParDef array does not match number of parameters (%d != %d)",iNumParDefs,iNumParTot),"FitInvMass"); return kFALSE; }
  if(iNumParDefs != iNumParLimLow) { Error(Form("Different length of arrays with parameter defauls and low limit values (%d != %d).",iNumParDefs,iNumParLimLow),"FitInvMass"); return kFALSE; }
  if(iNumParDefs != iNumParLimHigh) { Error(Form("Different length of arrays with parameter defauls and high limit values (%d != %d).",iNumParDefs,iNumParLimHigh),"FitInvMass"); return kFALSE; }

  // check the output of the vector assigment
  Debug("Fitting parameters setting done","FitInvMass");

  // === Initialision ===

  // master formula used in the fitting procedure
  if(!fbDebug) { sFitOptMass += "Q"; } // quite fitting option if NOT in debug

  TString sFuncMass = Form("%s + %s",sMassBG.Data(),sMassSig.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh), "FitInvMass");
  Debug(Form("Fit :\n    %s",sFuncMass.Data()), "FitInvMass");

  // changes the axis
  hist->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // === Fitting procedure ===

  //  fitting invariant mass distribution
  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetNpx(iNpx);

  if(!SetFuncParameters(fitMass, dParDef,dParLimLow,dParLimHigh,sParNames)) { Error("Setting default fitMass parameters failed!","FitInvMass"); return kFALSE; }

  // fitting
  TVirtualFitter::SetMaxIterations(10000);

  Debug("Fitting : 'fitMass'","FitInvMass");

  Bool_t bFitOK = kFALSE;
  Int_t nfitsA = 1;

  while(!bFitOK && (nfitsA < 30))
  {
    if(nfitsA > 1)
    {
        Double_t* dOldPars = fitMass->GetParameters();
        dOldPars[0] = fitMass->GetParameter(0)/nfitsA;

        if(!SetFuncParameters(fitMass, dOldPars,dParLimLow,dParLimHigh,sParNames)) { Error(Form("Setting fitMass parameters failed! (iteration %d)",nfitsA),"FitInvMass"); return kFALSE; }
    }

    bFitOK = CheckFitResult(hist->Fit(fitMass, sFitOptMass.Data()), nfitsA > 10);

    nfitsA++;
  }

  // === Extracting fitting components to separated TF1's ===
  TF1* fitMass_partBg = new TF1("fitMass_partBg",sMassBG.Data(),dMassRangeLow,dMassRangeHigh);
  fitMass_partBg->SetLineColor(kBlue);
  fitMass_partBg->SetLineStyle(2);

  TF1* fitMass_partSig = new TF1("fitMass_partSig", sMassSig.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass_partSig->SetLineColor(kGreen+2);
  fitMass_partSig->SetLineStyle(2);

  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) {
    fitMass_partBg->SetParameter(iPar, fitMass->GetParameter(iPar));
    fitMass_partBg->SetParError(iPar, fitMass->GetParError(iPar));

    fitMass_partSig->SetParameter(iPar, 0.0);
    fitMass_partSig->SetParError(iPar, 0.0);
  }

  for(Int_t iPar(iNumParsMassBG); iPar < iNumParTot; ++iPar) {
    fitMass_partSig->SetParameter(iPar, fitMass->GetParameter(iPar));
    fitMass_partSig->SetParError(iPar, fitMass->GetParError(iPar));
  }

  outList->Add(fitMass);
  outList->Add(fitMass_partSig);
  outList->Add(fitMass_partBg);

  // making "hard" fractions insted of individual components
  TF1* fitMass_fracSigHard = new TF1("fitMass_fracSigHard",Form("(%s)/(%s+%s)",sMassSig.Data(), sMassBG.Data(), sMassSig.Data()), dMassRangeLow, dMassRangeHigh);
  fitMass_fracSigHard->SetParameters(fitMass->GetParameters());
  fitMass_fracSigHard->SetParErrors(fitMass->GetParErrors());
  TF1* fitMass_fracBgHard = new TF1("fitMass_fracBgHard",Form("(%s)/(%s+%s)",sMassBG.Data(), sMassBG.Data(), sMassSig.Data()), dMassRangeLow, dMassRangeHigh);
  fitMass_fracBgHard->SetParameters(fitMass->GetParameters());
  fitMass_fracBgHard->SetParErrors(fitMass->GetParErrors());
  for(Int_t par(0); par < fitMass->GetNpar(); ++par) {
      Double_t dLow, dHigh;
      fitMass->GetParLimits(par,dLow,dHigh);
      fitMass_fracSigHard->SetParLimits(par,dLow,dHigh);
      fitMass_fracBgHard->SetParLimits(par,dLow,dHigh);
  }
  // PrintFitFunction(fitMass);
  // PrintFitFunction(fitMass_fracSigHard);
  // PrintFitFunction(fitMass_fracBgHard);

  fitOut = *fitMass;
  fitOutSig = *fitMass_fracSigHard;
  fitOutBg = *fitMass_fracBgHard;

  outList->Add(fitMass);
  outList->Add(fitMass_fracSigHard);
  outList->Add(fitMass_fracBgHard);

  if(!bFitOK) {
      Error("fitMass failed!","FitInvMass");
      PrintFitFunction(fitMass);
      // return kFALSE;
  }

  // === Preparing fractions ==
  Debug("Working on proper fractions!","FitInvMass");

  // Cloning original inv.mass distrubution : serves as baseline for fractions
  TH1D* hist_copy = nullptr;

  // Subtracting LS background out of US histo
  if(species == kPhi && task->fFlowFitPhiSubtLS) {
      histBg->SetName("histMass_bgLS");
      outList->Add(histBg);

      TH1D* histMass_subtLS = (TH1D*) SubtractInvMassBg(hist, histBg, task);
      if(!histMass_subtLS) { Error("Inv. mass BG subtraction failed!","FitInvMass"); return kFALSE; }
      histMass_subtLS->SetName("histMass_subtLS");
      outList->Add(histMass_subtLS);

      // setting subtLs as a baseline for all following operations
      hist_copy = histMass_subtLS;
  } else {
      hist_copy = (TH1D*) hist->Clone("hist_copy");
  }

  // excluding peak region
  TH1D* histMass_exclPeak = (TH1D*) hist_copy->Clone("histMass_exclPeak");
  for(Int_t iBin(0); iBin < histMass_exclPeak->GetNbinsX()+2; ++iBin) {
      Double_t x = histMass_exclPeak->GetBinCenter(iBin);
      if(x > dPeakLow && x < dPeakHigh) { histMass_exclPeak->SetBinError(iBin,2.0*dMaximum); }
  }
  outList->Add(histMass_exclPeak);

  // Fitting BG only in histMass_exclPeak
  Debug("Fitting : 'fitMass_exclPeak'","FitInvMass");
  TF1* fitMass_exclPeak = new TF1("fitMass_exclPeak",sMassBG.Data(),dMassRangeLow,dMassRangeHigh);
  if(!SetFuncParameters(fitMass_exclPeak, dParDef, dParLimLow, dParLimHigh)) { Error("Setting fitMass_exclPeak parameters failed!","FitInvMass"); return kFALSE; }
  if(bFitOK) { for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) { fitMass_exclPeak->SetParameter(iPar, fitMass->GetParameter(iPar)); } }
  for(Int_t iPar(iNumParsMassBG); iPar < iNumParTot; ++iPar) { fitMass_exclPeak->FixParameter(iPar,0.0); }
  outList->Add(fitMass_exclPeak);
  if(!CheckFitResult(histMass_exclPeak->Fit(fitMass_exclPeak, sFitOptMass.Data()), kTRUE) ) { Error("fitMass_exclPeak failed!","FitInvMass"); PrintFitFunction(fitMass_exclPeak); return kFALSE; }

  // Subtracting BG from total hist_copy (signal only)
  TH1D* histMass_sig = (TH1D*) hist_copy->Clone("histMass_sig");
  histMass_sig->Add(fitMass_exclPeak,-1.0);
  outList->Add(histMass_sig);

  // Subtracting Sig from total hist (bg only)  (NOT 'hist_copy')
  TH1D* histMass_bg = (TH1D*) hist->Clone("histMass_bg");
  histMass_bg->Add(histMass_sig,-1.0);
  outList->Add(histMass_bg);

  // Fitting Signal in BG free
  Debug("Fitting : 'fitMass_sig'","FitInvMass");
  TF1* fitMass_sig = new TF1("fitMass_sig",sMassSig.Data(),dMassRangeLow,dMassRangeHigh);
  if(!SetFuncParameters(fitMass_sig, dParDef,dParLimLow, dParLimHigh, sParNames)) { Error("Setting fitMass_sig parameters failed!","FitInvMass"); return kFALSE; }
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) { fitMass_sig->FixParameter(iPar,0.0); }
  if(bFitOK) { for(Int_t iPar(iNumParsMassBG); iPar < iNumParTot; ++iPar) { fitMass_sig->SetParameter(iPar, fitMass->GetParameter(iPar)); } }
  outList->Add(fitMass_sig);
  if(!CheckFitResult(histMass_sig->Fit(fitMass_sig, sFitOptMass.Data()), kTRUE)) { Error("fitMass_sig failed!","FitInvMass"); PrintFitFunction(fitMass_sig); return kFALSE; }

  // Cloning fitMass_bg for consistency reason (QA)
  TF1* fitMass_bg = (TF1*) fitMass_exclPeak->Clone("fitMass_bg");
  outList->Add(fitMass_bg);

  // Making sum of sig + bg (QA reasons)
  TF1* fitMass_tot = new TF1("fitMass_tot",sFuncMass.Data(), dMassRangeLow, dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) { fitMass_tot->SetParameter(iPar,fitMass_bg->GetParameter(iPar)); }
  for(Int_t iPar(iNumParsMassBG); iPar < iNumParTot; ++iPar) { fitMass_tot->SetParameter(iPar,fitMass_sig->GetParameter(iPar)); }
  outList->Add(fitMass_tot);

  // Making fraction // NB: Here 'hist' has to be in denominator (not 'hist_copy' !!!)
  TH1D* histMass_fracSig = (TH1D*) histMass_sig->Clone("histMass_fracSig");
  histMass_fracSig->Divide(hist);
  outList->Add(histMass_fracSig);

  TH1D* histMass_fracBg = (TH1D*) histMass_bg->Clone("histMass_fracBg");
  histMass_fracBg->Divide(hist);
  outList->Add(histMass_fracBg);

  // Fitting fractions
  Debug("Fitting : 'fitMass_fracSig'","FitInvMass");
  TF1* fitMass_fracSig = (TF1*) fitMass_sig->Clone("fitMass_fracSig");
  if(fbDebug) { Debug("Fitting: pre", "FitInvMass"); PrintFitFunction(fitMass_fracSig); }
  if(!SetFuncParameters(fitMass_fracSig, fitMass_sig->GetParameters(),dParLimLow, dParLimHigh, sParNames)) { Error("Setting fitMass_fracSig parameters failed!","FitInvMass"); return kFALSE; }
  if(bFitOK) { for(Int_t iPar(iNumParsMassBG); iPar < iNumParTot; ++iPar) { fitMass_fracSig->SetParameters(iPar, fitMass_fracSigHard->GetParameter(iPar)); } }
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) { fitMass_fracSig->FixParameter(iPar,0.0); }
  Double_t dFracMax = histMass_fracSig->GetBinContent(histMass_fracSig->GetMaximumBin());
  // fitMass_fracSig->SetParameter(iNumParsMassBG, dFracMax - dFracLimLow); // contribution from sMassSig == 0
  // fitMass_fracSig->SetParLimits(iNumParsMassBG, ((dFracMax - dFracLimLow) > 0.0 ? (dFracMax - dFracLimLow) : 0.0) , ((dFracMax + dFracLimHigh) < 1.0 ? (dFracMax + dFracLimHigh) : 1.0)); // contribution from sMassSig == 0
  fitMass_fracSig->SetParameter(iNumParsMassBG, dFracMax); // contribution from sMassSig == 0
  fitMass_fracSig->SetParLimits(iNumParsMassBG, 0.0,1.0); // contribution from sMassSig == 0
  if(fbDebug) { Debug("Fitting: post", "FitInvMass"); PrintFitFunction(fitMass_fracSig); }

  outList->Add(fitMass_fracSig);
  // if(!CheckFitResult(histMass_fracSig->Fit(fitMass_fracSig, sFitOptMass.Data()), kTRUE)) { Error("fitMass_fracSig failed!","FitInvMass"); PrintFitFunction(fitMass_fracSig); return kFALSE; }
  Bool_t bFitFracOK = CheckFitResult(histMass_fracSig->Fit(fitMass_fracSig, Form("IL%s",sFitOptMass.Data())), kTRUE);

  TF1* fitMass_fracBg = new TF1("fitMass_fracBg",Form("1.0-(%s)",sMassSig.Data()),dMassRangeLow,dMassRangeHigh);
  fitMass_fracBg->SetParameters(fitMass_fracSig->GetParameters());
  outList->Add(fitMass_fracBg);

  if(species != kPhi) {
      // splitiign two gauss
      TF1* fitMass_fracSig_partGauss1 = (TF1*) fitMass_fracSig->Clone("fitMass_fracSig_partGauss1");
      fitMass_fracSig_partGauss1->SetParameter(iParMass_2,0.0);

      TF1* fitMass_fracSig_partGauss2 = (TF1*) fitMass_fracSig->Clone("fitMass_fracSig_partGauss2");
      fitMass_fracSig_partGauss2->SetParameter(iParMass,0.0);

      outList->Add(fitMass_fracSig_partGauss1);
      outList->Add(fitMass_fracSig_partGauss2);
  }

   if(!bFitFracOK) { Error("fitMass_fracSig failed!","FitInvMass"); PrintFitFunction(fitMass_fracSig); return kFALSE; }

  fitOut = *fitMass;
  fitOutSig = *fitMass_fracSig;
  fitOutBg = *fitMass_fracBg;

  Info(Form("Inv.mass distribution fit (fitMass): SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g; %d iterations)",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF(),fitMass->GetProb(),nfitsA), "FitInvMass");
  Info(Form("Inv.mass distribution fit (fitMass_fracSig): SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g)",fitMass_fracSig->GetChisquare(), fitMass_fracSig->GetNDF(),fitMass_fracSig->GetChisquare()/fitMass_fracSig->GetNDF(),fitMass_fracSig->GetProb()), "FitInvMass");

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::FitCorrelations(TH1* hist, FlowTask* task, TF1& fitOut, TF1& fitOutSig, TF1& fitOutBg, TF1& fitInSig, TF1& fitInBg, TList* outList, Bool_t bIsFour)
{
  if(!hist) { Error("Input histo not found!","FitCorrelations"); return kFALSE; }
  if(!task) { Error("FlowTask not found!","FitCorrelations"); return kFALSE; }
  if(!outList) { Error("Output TList outList not found!","FitCorrelations"); return kFALSE; }

  PartSpecies species = task->fSpecies;
  if(!IsSpeciesReconstructed(species)) { Error("Invalid species!","FitCorrelations"); return kFALSE; }

  // === Fitting parametrisation (species dependent default) ===

  Double_t dMassRangeLow = fitInSig.GetXmin();
  Double_t dMassRangeHigh = fitInSig.GetXmax();
  Double_t dMaximum = hist->GetMaximum();
  hist->SetName("histCorr");
  outList->Add(hist);

  Int_t iNumParMass = fitInSig.GetNpar();

  Int_t iNpx = 10000;
  TString sFitOptFlow = "RNIS";

  std::vector<Double_t> dParDef;
  std::vector<Double_t> dParLimLow;
  std::vector<Double_t> dParLimHigh;

  TString sFlowBG = TString();  Int_t iNumParsFlowBG = 0; // function for flow-mass (BG component)

  // ========== SET FITTING FUNCTIONS ===========

  // Species (independent) flow shape
  sFlowBG = Form("[%d]*x+[%d]", iNumParMass,iNumParMass+1); iNumParsFlowBG = 2;

  if(bIsFour) {
      dParDef = {1e-2,1e-3,1e-3};
      dParLimLow = {-1,-1,-0.5};
      dParLimHigh = {-1,-1,1.0};

  } else {
      dParDef = {0.0,1.0,0.5};
      dParLimLow = {-1,-1,-0.5};
      dParLimHigh = {-1,-1,1.0};
  }

  // ========== END :: SET FITTING FUNCTIONS ===========

  if(task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow >= task->fFlowFitRangeHigh) { Error("Wrong fitting ranges set!","FitCorrelations"); return kFALSE; }
  if(task->fFlowFitRangeLow > 0.0) { dMassRangeLow = task->fFlowFitRangeLow; }
  if(task->fFlowFitRangeHigh > 0.0) { dMassRangeHigh = task->fFlowFitRangeHigh; }

  Bool_t bUserPars = kFALSE;
  if(task->fNumParFlowBG > 0) {
      Debug(" Task flowBG set","FitCorrelations");
      Warning("Setting User parameters for correlations not implemented yet","FitCorrelations");
      return kFALSE;
      bUserPars = kTRUE;
      sFlowBG = task->fFlowFitFlowBG;
      iNumParsFlowBG = task->fNumParFlowBG;
  }

  // TODO: TO be implemented (???)
  // if(bUserPars)
  // {
  //   Warning("Setting User parameters for correlations not implemented yet","FitCorrelations");
  //   // TODO: Due to way how default parameters are passed -> need reimplementation
  //   // using std::vector<> a = {} syntax: Setparames({2,1,2}) ...
  //
  //   // dParDef.clear();
  //   // dParLimLow.clear();
  //   // dParLimHigh.clear();
  //   //
  //   // Int_t iNumParTot = task->fNumParFlowBG;
  //   // for(Int_t par(0); par < iNumParTot; ++par)
  //   // {
  //   //   dParDef.push_back(task->fFitParDefaults[par]);
  //   //   dParLimLow.push_back(task->fFitParLimLow[par]);
  //   //   dParLimHigh.push_back(task->fFitParLimHigh[par]);
  //   // }
  // }

  Int_t iNumParDefs = dParDef.size();
  Int_t iNumParLimLow = dParLimLow.size();
  Int_t iNumParLimHigh = dParLimHigh.size();

  // check the output of the vector assigment
  if(fbDebug) {
    Debug("Fittin setting done","FitCorrelations");
    printf("Form: %s\n",sFlowBG.Data());
    for(Int_t par(0); par < iNumParDefs; ++par) { printf("  par %d: %g (%g<%g)\n",par, dParDef.at(par), dParLimLow.at(par), dParLimHigh.at(par)); }
  }

  // Fitting only BG first (only for Phi)
  TF1* fitBgOnly = nullptr;
  // if(species == kPhi) {
  if(1) {
    fitBgOnly = new TF1("fitFlowBg",sFlowBG.Data(),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    for(Int_t p(0); p < iNumParMass; ++p) {
        fitBgOnly->FixParameter(p, 0.0);
    }

    for(Int_t p(iNumParMass); p < iNumParMass+2; ++p) {
        fitBgOnly->SetParameter(p, 0.0);
    }

    Double_t dPeakLow = 0.0;
    Double_t dPeakHigh = 0.0;

    if(species == kK0s) {
        dPeakLow = 0.44;
        dPeakHigh = 0.55;
    }

    if(species == kLambda) {
        dPeakLow = 1.106;
        dPeakHigh = 1.134;
    }

    if(species == kPhi) {
        dPeakLow = 1.007;
        dPeakHigh = 1.03;
    }

    TProfile* prof_temp = (TProfile*) hist->Clone("prof_temp");
    TH1D* hist_Bg = (TH1D*) prof_temp->ProjectionX("hist_Bg");
    Double_t dMaximum = hist_Bg->GetMaximum();

    for(Int_t iBin(0); iBin < hist_Bg->GetNbinsX()+2; ++iBin) {

        Double_t x = hist_Bg->GetBinCenter(iBin);
        if(x > dPeakLow && x < dPeakHigh) { hist_Bg->SetBinError(iBin,10*dMaximum); }
    }

    outList->Add(hist_Bg);
    outList->Add(fitBgOnly);

    if(!CheckFitResult( hist_Bg->Fit(fitBgOnly,sFitOptFlow.Data()) )) { Error(Form("Flow-mass BG only fit does not converged within iterations limit (%d)!",1), "FitCorrelations"); PrintFitFunction(fitBgOnly); return kFALSE; }
  }

  // === Initialision ===
  Int_t iParFlow = iNumParMass +iNumParsFlowBG; // index of Flow (vn/dn) parameter

  // master formula used in the fitting procedure
  if(!fbDebug) { sFitOptFlow += "Q"; } // quite fitting option if NOT in debug

  TString sFuncVn = Form("[%d]*(%s) + (%s)*(%s)", iParFlow, fitInSig.GetTitle(), sFlowBG.Data(), fitInBg.GetTitle());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh), "FitCorrelations");
  Debug(Form("Fit Flow :\n    %s\n",sFuncVn.Data()), "FitCorrelations");

  // changes the axis
  // hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  // hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // === Fitting procedure ===

  TF1* fitCorr = new TF1(Form("fitCorr"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);

  // fixing fraction from input fits
  for(Int_t par(0); par < iNumParMass; ++par) {
      fitCorr->FixParameter(par, fitInSig.GetParameter(par));
  }

  for(Int_t par(iNumParMass); par < iParFlow+1; ++par) {
    // if(species == kPhi && fitBgOnly && par < iParFlow) {
    if(fitBgOnly && par < iParFlow) {
        fitCorr->FixParameter(par, fitBgOnly->GetParameter(par));
    } else {
        // Here par-iNumParMass is to account for a fact that dParDef takes only flow part (vector index != parameter index)
        fitCorr->SetParameter(par, dParDef.at(par-iNumParMass));
        Debug(Form("Parameter %d : %f",par,dParDef.at(par-iNumParMass) ),"FitCorrelations");
        Double_t dLimLow = dParLimLow.at(par-iNumParMass);
        Double_t dLimHigh = dParLimHigh.at(par-iNumParMass);

        if(dLimLow > -1.0 && dLimHigh > -1.0) { fitCorr->SetParLimits(par, dLimLow,dLimHigh); }
        else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Flow-mass: Only one of the parameter limits is set (par %d). Fix this!",par),"FitCorrelations"); return kFALSE; }
    }

  }

  // fitCorr->SetParameter(iParFlow, 0.5);
  // fitCorr->SetParLimits(iParFlow, -0.5,1.0);

  // fitting

  TVirtualFitter::SetMaxIterations(10000);

  // NB: Currently only one iteration
  // // fitting
  // Int_t nfitsA = 1;
  // Bool_t bFitOK = kFALSE;
  //
  // while(!bFitOK && (nfitsA < 15))
  // {
  //   if(nfitsA > 1)
  //   {
  //     fitMass->SetParameter(0, fitMass->GetParameter()/nfitsA);
  //
  //     for(Int_t par(1); par < iNumParMass; ++par)
  //     {
  //       fitMass->SetParameter(par, fitMass->GetParameter(par));
  //
  //       Double_t dLimLow = dParLimLow.at(par);
  //       Double_t dLimHigh = dParLimHigh.at(par);
  //
  //       if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
  //       else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"FitInvMass"); return kFALSE; }
  //     }
  //   }

  //   hist->Fit(fitCorr, sFitOptFlow.Data());
  //
  //   TString statusA = gMinuit->fCstatu.Data();
  //   if(statusA.Contains("CONVERGED")) { bFitOK = kTRUE; }
  //   nfitsA++;
  // }

  outList->Add(fitCorr);
  if(!CheckFitResult( hist->Fit(fitCorr, sFitOptFlow.Data()) )) { Error(Form("Flow-mass fit does not converged within iterations limit (%d)!",1), "FitCorrelations"); PrintFitFunction(fitCorr); return kFALSE; }

  // === Extracting fitting components to separated TF1's ===

  TF1* fitFlowSig = new TF1("fitCorSig", Form("[%d]*(%s)", iParFlow, fitInSig.GetTitle()), dMassRangeLow,dMassRangeHigh);
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitCorBg", Form("(%s)*(%s)", sFlowBG.Data(), fitInBg.GetTitle()), dMassRangeLow,dMassRangeHigh);
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);

  for(Int_t iPar(0); iPar < iNumParMass; ++iPar) {
      fitFlowSig->SetParameter(iPar, fitInSig.GetParameter(iPar));
      fitFlowSig->SetParError(iPar, fitInSig.GetParError(iPar));

      fitFlowBg->SetParameter(iPar, fitInSig.GetParameter(iPar));
      fitFlowBg->SetParError(iPar, fitInSig.GetParError(iPar));
  }

  for(Int_t iPar(iNumParMass); iPar < iParFlow; ++iPar) {
      fitFlowSig->SetParameter(iPar, 0.0);
      fitFlowSig->SetParError(iPar, 0.0);

      fitFlowBg->SetParameter(iPar, fitCorr->GetParameter(iPar));
      fitFlowBg->SetParError(iPar, fitCorr->GetParError(iPar));
  }

  fitFlowSig->SetParameter(iParFlow, fitCorr->GetParameter(iParFlow));
  fitFlowSig->SetParError(iParFlow, fitCorr->GetParError(iParFlow));

  outList->Add(fitFlowSig);
  outList->Add(fitFlowBg);

  fitOut = *fitCorr;
  fitOutSig = *fitFlowSig;
  fitOutBg = *fitFlowBg;

  Info(Form("Flow-mass fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g)",fitCorr->GetChisquare(), fitCorr->GetNDF(),fitCorr->GetChisquare()/fitCorr->GetNDF(),fitCorr->GetProb()), "FitCorrelations");

  return kTRUE;
}
//_____________________________________________________________________________
Int_t ProcessUniFlow::ReturnThird(Int_t first, Int_t second)
{
  //looking for the third index.. maybe there is a more elegant solution?
  Int_t third = -1;
  for(Int_t i(0); i < 3; i++){
    if(i != first && i != second) third = i;
  }
  return third;
}
//_____________________________________________________________________________
Int_t ProcessUniFlow::ReturnIndex3sub(Int_t index)
{
  //looking for the corect index.. maybe there is a more elegant solution?
  // reference??
  // check cpp book... or learn cpp...
  Int_t returnValue;
  if(index < 2)
    return returnValue = index + 1;
  else
    return 0;
}
//_____________________________________________________________________________
void ProcessUniFlow::Fatal(TString sMsg, TString sMethod)
{
	printf("\033[91mFatal::%s  %s. Terminating!\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Error(TString sMsg, TString sMethod)
{
	printf("\033[91mError::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Info(TString sMsg, TString sMethod)
{
	printf("\033[96mInfo::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Warning(TString sMsg, TString sMethod)
{
	printf("\033[93mWarning::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Debug(TString sMsg, TString sMethod)
{
	if(fbDebug) printf("\033[95mDebug::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
