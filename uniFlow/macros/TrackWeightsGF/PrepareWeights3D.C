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
 #include "TSystem.h"
 #include "TPad.h"

void CheckWeights(const TH3D* weight = 0x0, const Int_t runNumber = -1);
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh);
TH2D* ProcessWeights(const TH3D* weight = 0x0);
TH3D* ProcessWeights3D(const TH3D* dist);
TH1D* CalculateWeight(TH1D* proj);

const TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pp-16l";
const TString sOutputPath = sPath+"/weights/";
const Bool_t bRunByRun = kFALSE;

void PrepareWeights3D()
{
  const TString sTaskTag = "UniFlow";

  TFile* fOutput = new TFile(Form("%s/weights_Cor_CENT_woSDD_16q_FB768.root",sPath.Data()),"RECREATE");

  const Short_t iNumPart = 8;
  const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
  // const Short_t iNumPart = 1;
  // const TString species[iNumPart] = {"Refs"};

  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  TList* listRun = new TList();
  if(!bRunByRun)
  {
    TFile* fInput = TFile::Open(Form("%s/AnalysisResults_intermerge.root",sPath.Data()),"READ");
    if(!fInput) return;

    fInput->cd(sTaskTag.Data());
    TList* list =  (TList*) gDirectory->Get("Flow_Weights_UniFlow");
    if(!list) return;

    for(Short_t part = 0; part < iNumPart; part++)
    {
      printf(" -part %d (out of %d)\n",part+1,iNumPart);
      // h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
      TH3D* h3Weights = (TH3D*) list->FindObject(Form("fh3Weights%s",species[part].Data()));
      if(!h3Weights) { printf("Hist 'fh3Weights%s' not found\n",species[part].Data()); continue; }

      TCanvas* canDist = new TCanvas("canDist","canDist",1600,2000);
      canDist->Divide(5,4);
      for(Int_t iBinZ(1); iBinZ < h3Weights->GetNbinsZ()+1; ++iBinZ)
      {
        TH2D* slice = SliceWeights(h3Weights,iBinZ,iBinZ);
        slice->SetStats(0);
        // if(!slice) continue;

        canDist->cd(iBinZ);
        // gPad->SetLogz();
        slice->DrawCopy("colz");
      }

      canDist->SaveAs(Form("%s/dist_%s_all.pdf",sOutputPath.Data(), species[part].Data()), "pdf");

      TH3D* weights = ProcessWeights3D(h3Weights);
      if(!weights) { printf("Hist with 3D weights not processed!\n"); continue; }

      // CheckWeights(h3Weights, iRunList[iRun]);
      weights->SetName(Form("%s",species[part].Data()));
      weights->SetTitle(Form("%s",species[part].Data()));
      listRun->Add(weights);

      // Plot slices in z of weights
      TCanvas* canSlices = new TCanvas("canSlices","canSlices",1600,2000);
      canSlices->Divide(5,4);
      for(Int_t iBinZ(1); iBinZ < weights->GetNbinsZ()+1; ++iBinZ)
      {
        TH2D* slice = SliceWeights(weights,iBinZ,iBinZ);
        slice->SetStats(0);
        // if(!slice) continue;

        canSlices->cd(iBinZ);
        // gPad->SetLogz();
        slice->DrawCopy("colz");
      }

      canSlices->SaveAs(Form("%s/weights_%s_all.pdf",sOutputPath.Data(), species[part].Data()), "pdf");
    }

    fOutput->cd();
    listRun->Write("weights",TObject::kSingleKey);
  }
  else
  {
    printf("WARNING Run-By-Run need reimplementation !!!\n");
    return;

    // const Short_t iNumRuns = 1; const Int_t iRunList[iNumRuns] = {265387};
    const Short_t iNumRuns = 31; const Int_t iRunList[iNumRuns] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
    // const Short_t iNumRuns = 16; const Int_t iRunList[iNumRuns] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387};//,265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};

    // const Short_t iNumRuns = 4;  const Int_t iRunList[iNumRuns] = {267166, 267165, 267164, 267163};
    // const TString sOutputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_weights_V0A/FAST_16q"

    TFile* fileTemp = 0x0;
    TList* listTemp = 0x0;
    TH3D* h3Weights = 0x0;
    TH2D* h2Weights= 0x0;

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
        // h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
        h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
        if(!h3Weights) { printf("Run %d | Hist 'fh3Weights%s' not found\n",iRunList[iRun],species[part].Data()); continue; }
        h2Weights = ProcessWeights(h3Weights);
        // CheckWeights(h3Weights, iRunList[iRun]);

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
  weights_2d->SaveAs(Form("%s/weights_dist.pdf",sPath.Data()));

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
TH2D* ProcessWeights(const TH3D* weights)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH2D* weights_2d = (TH2D*) weights->Project3D("xy");

  // preparing weigths:  bin content = max / bincontent
  const Short_t iNumBinX = weights_2d->GetNbinsX();
  const Short_t iNumBinY = weights_2d->GetNbinsY();
  const Double_t dMax = weights_2d->GetMaximum();

  Double_t dContent = 0.0;

  for(Short_t binX(0); binX < iNumBinX+2; binX++)
    for(Short_t binY(0); binY < iNumBinY+2; binY++)
    {
      dContent = weights_2d->GetBinContent(binX,binY);
      if(dContent > 0.0) { weights_2d->SetBinContent(binX,binY, dMax / dContent); }
      else {weights_2d->SetBinContent(binX,binY, 1.0); }
      weights_2d->SetBinError(binX,binY,0.0);
    }

  return weights_2d;
}
//_____________________________________________________________________________
TH3D* ProcessWeights3D(const TH3D* dist)
{
  if(!dist) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH3D* distTemp = (TH3D*) dist->Clone("distTemp");
  // distTemp->Draw();

  // phi dist is on x-axis
  TH3D* weights = (TH3D*) dist->Clone("weights");
  weights->Reset();

  // preparing weigths:  bin content = max / bincontent
  for(Short_t binY(1); binY < dist->GetNbinsY()+1; binY++)
  {
    for(Short_t binZ(1); binZ < dist->GetNbinsZ()+1; binZ++)
    {
      distTemp->GetYaxis()->SetRange(binY,binY);
      distTemp->GetZaxis()->SetRange(binZ,binZ);

      TH1D* proj = (TH1D*) distTemp->Project3D(Form("x_%d_%d",binY,binZ));
      if(!proj) { printf("No projection! Something went wrong!\n"); return 0x0; }


      // calculating weights in this 1D projection
      TH1D* hWeight = CalculateWeight(proj);
      for(Short_t binPhi(1); binPhi < hWeight->GetNbinsX()+1; binPhi++)
      {
        weights->SetBinContent(binPhi,binY,binZ, hWeight->GetBinContent(binPhi));
        weights->SetBinError(binPhi,binY,binZ, 0.0);
      }

    }
  }

  return weights;
}
//_____________________________________________________________________________
TH1D* CalculateWeight(TH1D* proj)
{
  if(!proj) { printf("CalculateWeight-ERROR: Input proj not found!\n"); return 0x0; }

  TH1D* hWeight = (TH1D*) proj->Clone("weight");
  if(!hWeight) { printf("CalculateWeight-ERROR: hWeight not cloned!\n"); return 0x0; }

  Double_t dMax = proj->GetMaximum();
  // printf("dMax %f\n",dMax);

  for(Int_t iBin(0); iBin < proj->GetNbinsX()+2; ++iBin)
  {
    Double_t dContent = proj->GetBinContent(iBin);
    if(dContent > 0.0) { hWeight->SetBinContent(iBin, dMax / dContent); }
    else { hWeight->SetBinContent(iBin, 1.0); }
    hWeight->SetBinError(iBin, 0.0);
  }

  return hWeight;
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
