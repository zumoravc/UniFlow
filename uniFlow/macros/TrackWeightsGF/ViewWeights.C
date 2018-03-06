#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"


TH2D* LoadWeights(const TH3D* weights);


void ViewWeights()
{
  TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pp-16l";
  TString sTaskTag = "UniFlow";
  TString species[] = {"Phi"};
  Int_t iNumPart = sizeof(species)/sizeof(species[0]);



  TFile* fInput = TFile::Open(Form("%s/AnalysisResults_intermerge.root",sPath.Data()),"READ");
  fInput->cd(sTaskTag.Data());
  TList* list =  (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));

  for(Short_t part = 0; part < iNumPart; part++)
  {
    printf(" -part %d (out of %d)\n",part+1,iNumPart);
    // h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
    TH3D* h3Weights = (TH3D*) list->FindObject(Form("fh3Weights%s",species[part].Data()));
    if(!h3Weights) { printf("Hist 'fh3Weights%s' not found\n",species[part].Data()); continue; }
    TH2D* h2Weights = LoadWeights(h3Weights);
    // CheckWeights(h3Weights, iRunList[iRun]);
    h2Weights->SetName(Form("%s",species[part].Data()));
    h2Weights->SetTitle(Form("%s",species[part].Data()));
    // listRun->Add(h2Weights);

    h2Weights->Draw("colz");
  }

}

TH2D* LoadWeights(const TH3D* weights)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH2D* weights_2d = (TH2D*) weights->Project3D("xy");
  return weights_2d;
}
//_____________________________________________________________________________
