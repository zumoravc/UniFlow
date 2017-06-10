/* PrepareWeights
 *
 * Macro for extracting the phi,eta,pt weights for Generic Framework calculations
 * and preparing source ROOT file used for running the analysis (on run-by-run basis).
 * Also contains QA tools for weights.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

 #include "TFile.h"
 #include "TObject.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "TH3.h"
 #include "TCanvas.h"

void CheckWeights(const TH3D* weight = 0x0, const Int_t runNumber = -1);
TH2D* ProcessWeights(const TH3D* weight = 0x0, const Int_t runNumber = -1);

void PrepareWeights()
{
  // const Short_t iNumRuns = 1;
  // const Int_t iRunList[iNumRuns] = {265387};
  const Short_t iNumRuns = 31;
  const Int_t iRunList[iNumRuns] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387,265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
  // const TString sOutputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_weights_V0A/FAST_16q"
  const TString sPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst/baseline/CENT_woSDD_16q/";
  const TString sTaskTag = "UniFlow";

  TFile* fOutput = new TFile(Form("%s/weights_Cor_CENTwoSDD_16q.root",sPath.Data()),"RECREATE");

  const Short_t iNumPart = 8;
  const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};

  TFile* fileTemp = 0x0;
  TList* listTemp = 0x0;
  TH3D* h3Weights = 0x0;
  TH2D* h2Weights= 0x0;

  TList* listRun = new TList();

  for(Short_t iRun = 0; iRun < iNumRuns; iRun++)
  {
    printf(" ### Run %d (out of %d)\n",iRun+1,iNumRuns);
    fileTemp = new TFile(Form("%s/merge/merge_%d/AnalysisResults.root",sPath.Data(),iRunList[iRun]),"READ");
    if(!fileTemp) { printf("Run %d | Input file with weights not found\n",iRunList[iRun]); continue; }
    fileTemp->cd(sTaskTag.Data());

    listTemp = (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));
    if(!listTemp) { printf("Run %d | TList with weights not found\n",iRunList[iRun]); continue; }

    // listTemp->ls();

    for(Short_t part = 0; part < iNumPart; part++)
    {
      printf(" -part %d (out of %d)\n",part+1,iNumPart);
      h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
      if(!h3Weights) { printf("Run %d | Hist 'fh3Weights%s' not found\n",iRunList[iRun],species[part].Data()); continue; }
      h2Weights = ProcessWeights(h3Weights, iRunList[iRun]);
      // if(part == 7)
      // {
      //   for(Short_t binX(1); binX < h2Weights->GetNbinsX()+1; binX++)
      //     for(Short_t binY(1); binY < h2Weights->GetNbinsX()+1; binY++)
      //     { h2Weights->SetBinContent(binX,binY,1.); }
      // }
      h2Weights->SetName(Form("%s",species[part].Data()));
      h2Weights->SetTitle(Form("%s | Run %d",species[part].Data(),iRunList[iRun]));
      listRun->Add(h2Weights);

    }

    fOutput->cd();
    listRun->Write(Form("%d",iRunList[iRun]),TObject::kSingleKey);
    listRun->Clear();
    delete fileTemp;
  }

  return;
}
//_____________________________________________________________________________
void CheckWeights(const TH3D* weights, const Int_t runNumber)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return;}
  if(runNumber == -1) { printf("PrepareWeights: run number not specified\n"); return;}

  const Short_t bins = 5;

  TH3D* weights_temp = 0x0;
  TH2D* weights_2d = 0x0;
  Short_t binLow = 0;
  Short_t binHigh = 0;

  TCanvas* canCheck = new TCanvas("canCheck","canCheck",1000,1000);
  canCheck->Divide(3,3);
  canCheck->cd(1);

  weights_temp = (TH3D*) weights->Clone("weights_temp_int");
  weights_2d = (TH2D*) weights_temp->Project3D("xy");
  weights_2d->Draw("colz");

  for(Short_t pt(0); pt < 8; pt++)
  {
    binLow = pt*bins+1;
    binHigh = (pt+1)*bins;
    printf("bin low %d | high %d\n",binLow,binHigh);
    weights_temp = (TH3D*) weights->Clone(Form("weights_temp_%d",pt));
    weights_temp->GetZaxis()->SetRange(binLow,binHigh);
    weights_2d = (TH2D*) weights_temp->Project3D("xy");
    weights_2d->Scale(1/weights_2d->GetEntries());

    canCheck->cd(pt+2);
    weights_2d->Draw("colz");
  }

  return;
}
//_____________________________________________________________________________
TH2D* ProcessWeights(const TH3D* weights, const Int_t runNumber)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}
  if(runNumber == -1) { printf("PrepareWeights: run number not specified\n"); return 0x0; }

  TH2D* weights_2d = (TH2D*) weights->Project3D("xy");

  // preparing weigths:  bin content = max / bincontent
  const Short_t iNumBinX = weights_2d->GetNbinsX();
  const Short_t iNumBinY = weights_2d->GetNbinsY();
  const Double_t dMax = weights_2d->GetMaximum();

  Double_t dContent = 0;

  for(Short_t binX(0); binX < iNumBinX+2; binX++)
    for(Short_t binY(0); binY < iNumBinY+2; binY++)
    {
      dContent = weights_2d->GetBinContent(binX,binY);
      if(dContent > 0.) { weights_2d->SetBinContent(binX,binY, dMax / dContent); }
      else {weights_2d->SetBinContent(binX,binY, 1); }

      weights_2d->SetBinError(binX,binY,0);
    }

  // weights_2d->Scale(1/weights_2d->GetMaximum());
  // weights_2d->Scale(1/weights_2d->GetEntries());
  // weights_2d->SetName(Form("%s_%d",species.Data(),runNumber));
  // weights_2d->SetTitle(Form("%s | Run %d",species.Data(),runNumber));
  // weights_2d->Draw("colz");

  return weights_2d;
}
//_____________________________________________________________________________
