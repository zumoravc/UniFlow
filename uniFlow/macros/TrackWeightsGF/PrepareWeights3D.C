/* PrepareWeights
 *
 * Macro for extracting the (eta,phi,z) distribution and make them into 3D (eta,phi,z) and 2D (eta,phi) weights
 * for Generic Framework calculations and preparing source ROOT file used for running the analysis.
 *
 * This can be done on run-integrated results and by run-by-run basis (if bRunByRun = kTRUE)
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


const TString sTaskTag = "UniFlow";
const TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pp-16l-run2";
const TString sOutputPath = sPath+"/weights_16l/";
const TString sOutFileName = "weights_16l.root";

const Short_t iNumPart = 8;
const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
// const Short_t iNumPart = 1;
// const TString species[iNumPart] = {"Refs"};

const Bool_t bRunByRun = kTRUE;
// Run Lists

// pPb 5.02 TeV
// RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20171129_v2.txt [31 runs]
// Int_t iNumRuns = 31; Int_t iRunList[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};


// pp 13 TeV
// RunList_LHC16l_pass1_CentralBarrelTracking_hadronPID_20170509_v2.txt [70]
Int_t iNumRuns = 70; Int_t iRunList[] = {259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259713, 259711, 259705, 259704, 259703, 259700, 259697, 259668, 259650, 259649, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962, 258923, 258919};


TH2D* ProcessWeights2D(const TH3D* dist);
TH3D* ProcessWeights3D(const TH3D* dist);
TH1D* CalculateWeight(TH1D* proj);
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh);

void PrepareWeights3D()
{
  gSystem->mkdir(sOutputPath.Data(),kTRUE);
  TFile* fOutput = new TFile(Form("%s/%s",sPath.Data(),sOutFileName.Data()),"RECREATE");

  TList* listRun = new TList();

  TFile* fInput = TFile::Open(Form("%s/AnalysisResults.root",sPath.Data()),"READ");
  if(!fInput) return;

  fInput->cd(sTaskTag.Data());
  TList* list =  (TList*) gDirectory->Get("Flow_Weights_UniFlow");
  if(!list) return;

  for(Short_t part = 0; part < iNumPart; part++)
  {
    printf(" -part %d (out of %d)\n",part+1,iNumPart);
    TH3D* h3Weights = (TH3D*) list->FindObject(Form("fh3Weights%s",species[part].Data()));
    if(!h3Weights) { printf("Hist 'fh3Weights%s' not found\n",species[part].Data()); continue; }

    // Plotting overview plots with dist slices in z
    // TCanvas* canDist = new TCanvas("canDist","canDist",1600,2000);
    // canDist->Divide(5,4);
    // for(Int_t iBinZ(1); iBinZ < h3Weights->GetNbinsZ()+1; ++iBinZ)
    // {
    //   TH2D* slice = SliceWeights(h3Weights,iBinZ,iBinZ);
    //   slice->SetStats(0);
    //   // if(!slice) continue;
    //
    //   canDist->cd(iBinZ);
    //   // gPad->SetLogz();
    //   slice->DrawCopy("colz");
    // }
    // canDist->SaveAs(Form("%s/dist_%s_all.pdf",sOutputPath.Data(), species[part].Data()), "pdf");

    TH2D* weights = ProcessWeights2D(h3Weights);
    if(!weights) { printf("Hist with 2D weights not processed!\n"); continue; }

    // CheckWeights(h3Weights, iRunList[iRun]);
    weights->SetName(Form("%s",species[part].Data()));
    weights->SetTitle(Form("%s",species[part].Data()));
    listRun->Add(weights);

    TH3D* weights3D = ProcessWeights3D(h3Weights);
    if(!weights3D) { printf("Hist with 3D weights not processed!\n"); continue; }

    // CheckWeights(h3Weights, iRunList[iRun]);
    weights3D->SetName(Form("%s3D",species[part].Data()));
    weights3D->SetTitle(Form("%s3D",species[part].Data()));
    listRun->Add(weights3D);

    // Plotting overview plots with dist slices in z
    // TCanvas* canSlices = new TCanvas("canSlices","canSlices",1600,2000);
    // canSlices->Divide(5,4);
    // for(Int_t iBinZ(1); iBinZ < weights->GetNbinsZ()+1; ++iBinZ)
    // {
    //   TH2D* slice = SliceWeights(weights3D,iBinZ,iBinZ);
    //   slice->SetStats(0);
    //   // if(!slice) continue;
    //
    //   canSlices->cd(iBinZ);
    //   // gPad->SetLogz();
    //   slice->DrawCopy("colz");
    // }
    // canSlices->SaveAs(Form("%s/weights_%s_all.pdf",sOutputPath.Data(), species[part].Data()), "pdf");

  } // end-for{species}

  fOutput->cd();
  listRun->Write("weights",TObject::kSingleKey);

  if(bRunByRun)
  {
    for(Short_t iRun = 0; iRun < iNumRuns; iRun++)
    {
      printf(" ### Run %d (out of %d)\n",iRun+1,iNumRuns);
      TFile* fileTemp = new TFile(Form("%s/merge/merge_%d/AnalysisResults.root",sPath.Data(),iRunList[iRun]),"READ");
      if(!fileTemp) { printf("Run %d | Input file with weights not found\n",iRunList[iRun]); return; }

      fileTemp->cd(sTaskTag.Data());

      TList* listTemp = (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));
      if(!listTemp) { printf("Run %d | TList with weights not found\n",iRunList[iRun]); gDirectory->ls() ;return; }

      for(Short_t part = 0; part < iNumPart; part++)
      {
        printf(" -part %d (out of %d)\n",part+1,iNumPart);
        TH3D* h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
        if(!h3Weights) { printf("Run %d | Hist 'fh3Weights%s' not found\n",iRunList[iRun],species[part].Data()); listTemp->ls(); return; }

        TH2D* h2Weights = ProcessWeights2D(h3Weights);

        h2Weights->SetName(Form("%s",species[part].Data()));
        h2Weights->SetTitle(Form("%s | Run %d",species[part].Data(),iRunList[iRun]));
        listRun->Add(h2Weights);
      }

      fOutput->cd();
      listRun->Write(Form("%d",iRunList[iRun]),TObject::kSingleKey);
      listRun->Clear();
    }
  }

  return;
}
//_____________________________________________________________________________
TH2D* ProcessWeights2D(const TH3D* dist)
{
  if(!dist) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH2D* dist_2d = (TH2D*) dist->Project3D("xy");
  TH2D* weights_2d = (TH2D*) dist_2d->Clone();
  weights_2d->Reset();

  for(Short_t binEta(1); binEta < dist_2d->GetNbinsX()+1; ++binEta)
  {
    // projection onto phi (in eta bin)
    TH1D* proj = (TH1D*) dist_2d->ProjectionY(Form("y_%d",binEta), binEta,binEta);
    if(!proj) { printf("No projection! Something went wrong!\n"); return 0x0; }

    TH1D* hWeight = CalculateWeight(proj);

    for(Short_t binPhi(1); binPhi < weights_2d->GetNbinsY()+1; ++binPhi)
    {
      weights_2d->SetBinContent(binEta,binPhi, hWeight->GetBinContent(binPhi));
      weights_2d->SetBinError(binEta,binPhi, 0.0);
    }

  }

  return weights_2d;
}
//_____________________________________________________________________________
TH3D* ProcessWeights3D(const TH3D* dist)
{
  if(!dist) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH3D* distTemp = (TH3D*) dist->Clone("distTemp");
  // phi dist is on y-axis
  TH3D* weights = (TH3D*) dist->Clone("weights");
  weights->Reset();

  // preparing weigths:  bin content = max / bincontent

  // z-dimension
  for(Short_t binZ(1); binZ < dist->GetNbinsZ()+1; ++binZ)
  {
    // y-dimension : eta
    for(Short_t binEta(1); binEta < dist->GetNbinsX()+1; ++binEta)
    {
      distTemp->GetXaxis()->SetRange(binEta,binEta);
      distTemp->GetZaxis()->SetRange(binZ,binZ);

      // projection into phi
      TH1D* proj = (TH1D*) distTemp->Project3D(Form("y_%d_%d",binEta,binZ));
      if(!proj) { printf("No projection! Something went wrong!\n"); return 0x0; }

      // calculating weights in this 1D projection
      TH1D* hWeight = CalculateWeight(proj);
      for(Short_t binPhi(1); binPhi < hWeight->GetNbinsX()+1; ++binPhi)
      {
        weights->SetBinContent(binEta,binPhi,binZ, hWeight->GetBinContent(binPhi));
        weights->SetBinError(binEta,binPhi,binZ, 0.0);
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
