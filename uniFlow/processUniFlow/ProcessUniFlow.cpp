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

class ProcessUniFlow
{
  public:
                ProcessUniFlow();
                ~ProcessUniFlow();

    void        SetInputFilePath(const char* path) { fsInputFilePath = path; }
    void        SetInputFileName(const char* name) { fsInputFileName = name; }
    void        SetOutputFilePath(const char* path) { fsOutputFilePath = path; }
    void        SetOutputFileName(const char* name) { fsOutputFileName = name; }
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }
    void        Run(); // running the task (main body of the class)

  protected:

  private:
    Bool_t      Initialize(); // initialization task
    void        TestProjections(); // testing projection of reconstructed particles
    // printing output methods
    void        Fatal(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Error(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Warning(TString sMsg, TString sMethod = ""); // printf the msg as warning
    void        Info(TString sMsg, TString sMethod = ""); // printf the msg as info
    void        Debug(TString sMsg, TString sMethod = ""); // printf the msg as info

    Bool_t      fbInit; // flag for initialization status
    Bool_t      fbDebug; // flag for debugging : if kTRUE Debug() messages are displayed
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! input file container
    TString     fsInputFilePath; // path to the input folder with input file
    TString     fsInputFileName; // name of input file
    TString     fsOutputFilePath; // path to the ouput folder
    TString     fsOutputFileName; // name of output file
};
//_____________________________________________________________________________
ProcessUniFlow::ProcessUniFlow() :
  fbDebug(kFALSE),
  fbInit(kFALSE),
  ffInputFile(0x0),
  ffOutputFile(0x0)
{
  // default constructor
  fsInputFilePath = TString("");
  fsInputFileName = TString("AnalysisResults.root");
  fsOutputFilePath = TString("");
  fsOutputFileName = TString("UniFlow.root");
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

  TestProjections();

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


  // initialization succesfull
  fbInit = kTRUE;
  Info("Initialization succesfull","Initialize");
  return fbInit;
}
//_____________________________________________________________________________
void ProcessUniFlow::TestProjections()
{
  Info("Testing profile projections");

  ffInputFile->cd("UniFlow_CENT_wSDD");
  // ffInputFile->ls();

  TList* lFlow = (TList*) gDirectory->Get("Flow_UniFlow_CENT_wSDD");
  // lFlow->ls();


  TList* lRef = (TList*) lFlow->FindObject("fFlowRefs");
  TProfile* Ref = (TProfile*) lRef->FindObject("fpRefs_<2>_harm2_gap-10_sample0");
  if(!Ref) { Error("NotFound"); return; }
  // Ref->Draw();
  // NOTE: reference flow works

  // projections charged flow
  TList* lCharged = (TList*) lFlow->FindObject("fFlowCharged");
  // lCharged->ls();
  TProfile2D* p2Charged = (TProfile2D*) lCharged->FindObject("fp2Charged_<2>_harm2_gap-10_sample0");
  if(!p2Charged) { Error("NotFound"); return; }
  TProfile* p2ChargedProjY = p2Charged->ProfileY("p2ChargedProjY",10,10);

  // TCanvas* cCharged = new TCanvas("cCharged","cCharged",1000,1000);
  // cCharged->Divide(1,2);
  // cCharged->cd(1);
  // p2Charged->Draw("colz");
  // cCharged->cd(2);
  // p2ChargedProjY->Draw();

  // projections V0s
  TList* lV0s = (TList*) lFlow->FindObject("fFlowV0s");
  if(!lV0s) { Error("NotFound"); return; }
  lV0s->ls();

  // entries
  TH3D* h3K0sEntries = (TH3D*) lV0s->FindObject("fh3V0sEntriesK0s_gap-10")->Clone("h3K0sEntries");
  if(!h3K0sEntries) { Error("NotFound"); return; }
  TH1D* h3K0sEntriesProjX = h3K0sEntries->ProjectionX("h3K0sEntriesProjX"); // whole projeciton
  TH1D* h3K0sEntriesProjY = h3K0sEntries->ProjectionY("h3K0sEntriesProjY"); // whole projeciton
  TH1D* h3K0sEntriesProjZ = h3K0sEntries->ProjectionZ("h3K0sEntriesProjZ"); // whole projeciton
  TH1D* h3K0sEntriesProjZsub = h3K0sEntries->ProjectionZ("h3K0sEntriesProjZsub",10,12,10,12); // projection in certain cent & pt range

  TCanvas* cK0s = new TCanvas("cK0s","cK0s",1000,1000);
  cK0s->Divide(2,3);
  cK0s->cd(1);
  h3K0sEntries->Draw();
  cK0s->cd(2);
  h3K0sEntriesProjZ->Draw();
  cK0s->cd(3);
  h3K0sEntriesProjX->Draw();
  cK0s->cd(4);
  h3K0sEntriesProjY->Draw();
  cK0s->cd(5);
  h3K0sEntriesProjZsub->Draw();
  cK0s->cd(6);
  // NOTE seems to work properly

  // correlations
  TProfile3D* p3K0sCor = (TProfile3D*) lV0s->FindObject("fp3V0sCorrK0s_<2>_harm2_gap-10")->Clone("p3K0sCor");
  if(!p3K0sCor) { Error("NotFound"); return; }

  TProfile2D* p3K0sCorProjZY = p3K0sCor->Project3DProfile("yz"); // projection over centrality
  TProfile* p3K0sCorProjZ = p3K0sCorProjZY->ProfileX("p3K0sCorProjZ"); // projection over pT

  TH3D* p3K0sCorHist = p3K0sCor->ProjectionXYZ("p3K0sCorHist"); // making TH3D from TProfile3D
  p3K0sCorHist->Sumw2();
  TProfile3D* p3K0sCorSub = (TProfile3D*) lV0s->FindObject("fp3V0sCorrK0s_<2>_harm2_gap-10")->Clone("p3K0sCorSub");
  if(!p3K0sCorSub) { Error("NotFound"); return; }

  // p3K0sCorSub->GetYaxis()->SetRangeUser(10,20);

  // TProfile2D* p3K0sCorProjZYsub = p3K0sCorSub->Project3DProfile("xz"); // projection over centrality
  // TProfile* p3K0sCorProjZsub = p3K0sCorProjZYsub->ProfileX("p3K0sCorProjZsub"); // projection over pT

  // p3K0sCorHist->GetXaxis()->SetRange(10,20);

  TProfile2D* p3K0sCorProjZYsub = p3K0sCorHist->Project3DProfile("yz"); // projection over centrality
  TProfile* p3K0sCorProjZsub = p3K0sCorProjZYsub->ProfileX("p3K0sCorProjZsub"); // projection over pT


  Int_t binx = 10, biny = 50, binz = 1;
  printf("Profile: %f ± %f | ", p3K0sCor->GetBinContent(binx,biny,binz), p3K0sCor->GetBinError(binx,biny,binz));
  printf("TH: %f ± %f\n", p3K0sCorHist->GetBinContent(binx,biny,binz), p3K0sCorHist->GetBinError(binx,biny,binz));

  printf("Projected 2D from Profile: %f ± %f | ", p3K0sCorProjZY->GetBinContent(1,5), p3K0sCorProjZY->GetBinError(1,5));
  printf(" from TH: %f ± %f\n", p3K0sCorProjZYsub->GetBinContent(1,5), p3K0sCorProjZYsub->GetBinError(1,5));

  TCanvas* cK0sCor = new TCanvas("cK0sCor","cK0sCor",1000,1000);
  cK0sCor->Divide(2,3);
  cK0sCor->cd(1);
  p3K0sCor->Draw();
  cK0sCor->cd(2);
  // p3K0sCorHist->Draw();
  p3K0sCorSub->Draw();
  cK0sCor->cd(3);
  p3K0sCorProjZY->Draw("colz");
  cK0sCor->cd(4);
  p3K0sCorProjZYsub->Draw("colz");
  cK0sCor->cd(5);
  p3K0sCorProjZ->Draw();
  cK0sCor->cd(6);
  p3K0sCorProjZsub->Draw();




  return;
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
