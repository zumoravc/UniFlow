#include "TFile.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"


TH2D* LoadWeights(const TH3D* weights);

TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh);

void SliceWeights()
{
  TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pp-16l";
  TString sOutputPath = sPath+"/weights/";
  TString sTaskTag = "UniFlow";
  TString species[] = {"Refs", "Charged"};
  Int_t iNumPart = sizeof(species)/sizeof(species[0]);

  // setting output
  gSystem->mkdir(sOutputPath.Data(), kTRUE);
  TFile* fileOut = TFile::Open(Form("%s/weights.root",sOutputPath.Data()),"RECREATE");


  TFile* fInput = TFile::Open(Form("%s/AnalysisResults_intermerge.root",sPath.Data()),"READ");
  fInput->cd(sTaskTag.Data());
  TList* list =  (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));

  for(Short_t part = 0; part < iNumPart; part++)
  {
    printf(" -part %d (out of %d)\n",part+1,iNumPart);
    // h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
    TH3D* h3Weights = (TH3D*) list->FindObject(Form("fh3Weights%s",species[part].Data()));
    if(!h3Weights) { printf("Hist 'fh3Weights%s' not found\n",species[part].Data()); return; }


    TCanvas* canBins = new TCanvas("canBins","canBins", 1600,2000);
    canBins->Divide(5,4);

    for(Int_t iBinZ(1); iBinZ < h3Weights->GetNbinsZ()+1; ++iBinZ)
    {
      // TH2D* h2Weights = LoadWeights(h3Weights);

      Int_t iBinLow = iBinZ;
      Int_t iBinHigh = iBinZ;

      TH2D* h2Weights = SliceWeights(h3Weights,iBinZ,iBinHigh);

      // CheckWeights(h3Weights, iRunList[iRun]);
      h2Weights->SetName(Form("%s",species[part].Data()));
      h2Weights->SetTitle(Form("%s",species[part].Data()));
      // listRun->Add(h2Weights);

      TCanvas* canTemp = new TCanvas("canTemp","canTemp",1000,1000);
      canTemp->cd();
      h2Weights->Draw("colz");
      h2Weights->SaveAs(Form("%s/weight_%s_%d.pdf",sOutputPath.Data(),species[part].Data(),iBinLow),"pdf");

      canBins->cd(iBinZ);
      h2Weights->SetStats(0);
      h2Weights->Draw("colz");
    }
    canBins->SaveAs(Form("%s/weight_%s_all.pdf",sOutputPath.Data(),species[part].Data()),"pdf");
  }

  return;
}
//_____________________________________________________________________________
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH3D* weightsTemp = (TH3D*) weights->Clone();

  weightsTemp->GetZaxis()->SetRange(iBinLow,iBinHigh);

  TH2D* weights_2d = (TH2D*) weightsTemp->Project3D("xy");
  return weights_2d;
}
//_____________________________________________________________________________
TH2D* LoadWeights(const TH3D* weights)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH2D* weights_2d = (TH2D*) weights->Project3D("xy");
  return weights_2d;
}
//_____________________________________________________________________________
