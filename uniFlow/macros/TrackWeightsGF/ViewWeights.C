// =================================================================================================
// ViewWeights.C
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016-2018
//
// Macro created for checking NUA (phi,eta) 2D-weights for UniFlow output
// =================================================================================================

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TSystem.h"

TH2D* LoadWeights(const TH3D* weights);

void ViewWeights()
{
  Bool_t bSavePlots = kTRUE;
  TString sWeightName = "fh3AfterWeights";
  // TString sWeightName = "fh3Weights";

  // TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/tracking/merged-16q-nua/";
  // TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/tracking/pPb-16q-FAST-nua/";
  // TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pPb-16q-woSDD-nua/";
  // TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pPb-16q-FAST-nua/";
  TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/";
  TString sFileName = "AnalysisResults.root";
  TString sTaskTag = "UniFlow";
  // TString sTaskTag = "UniFlow_cls";
  TString species[] = {"Refs","Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
  // TString species[] = {"K0s"};

 // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if(bSavePlots) { gSystem->mkdir(Form("%s/weights/",sPath.Data()),kTRUE); }

  Int_t iNumPart = sizeof(species)/sizeof(species[0]);

  TFile* fInput = TFile::Open(Form("%s%s",sPath.Data(),sFileName.Data()),"READ");
  fInput->cd(sTaskTag.Data());
  TList* list =  (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));



  TCanvas* can = new TCanvas("can","can",iNumPart*400,400);
  can->Divide(iNumPart,2);

  for(Short_t part = 0; part < iNumPart; part++)
  {
    printf(" -part %d (out of %d)\n",part+1,iNumPart);
    TH3D* h3Weights = (TH3D*) list->FindObject(Form("%s%s",sWeightName.Data(),species[part].Data()));
    if(!h3Weights) { printf("Hist '%s%s' not found\n",sWeightName.Data(),species[part].Data()); return; }

    TH2D* h2Weights = LoadWeights(h3Weights);
    if(!h2Weights) { printf("Hist h2Weights '%s%s' not found\n",sWeightName.Data(),species[part].Data()); return; }
    h2Weights->SetName(Form("%s",species[part].Data()));
    h2Weights->SetTitle(Form("%s",species[part].Data()));

    TH1D* hWeights = (TH1D*) h2Weights->ProjectionX();
    can->cd(part+1);
    // gPad->SetLogz();
    // h2Weights->SetStats(0);
    h2Weights->DrawCopy("colz");
    can->cd(part+1+iNumPart);
    // hWeights->SetStats(0);
    hWeights->DrawCopy("colz");

    TCanvas* canPart = new TCanvas("canPart","canPart",400,600);
    canPart->Divide(1,2);
    canPart->cd(1);
    h2Weights->SetStats(0);
    h2Weights->DrawCopy("colz");
    canPart->cd(2);
    hWeights->SetMinimum(0.0);
    hWeights->SetMarkerStyle(kFullCircle);
    hWeights->SetMarkerSize(0.5);
    hWeights->SetStats(0);
    hWeights->DrawCopy("colz");
    if(bSavePlots) { canPart->SaveAs(Form("%s/weights/%s_%s_%s.pdf",sPath.Data(),sWeightName.Data(),sTaskTag.Data(),species[part].Data())); }
  }

  if(bSavePlots) { can->SaveAs(Form("%s/weights/%s_%s.pdf",sPath.Data(),sWeightName.Data(), sTaskTag.Data())); }

}

TH2D* LoadWeights(const TH3D* weights)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH2D* weights_2d = (TH2D*) weights->Project3D("yx");
  return weights_2d;
}
//_____________________________________________________________________________
