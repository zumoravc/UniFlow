/*
* Macro for plotting the comparison between FAST, CENTwSDD and CENTwoSDD trigger clusters
* QA and flow results.
* It contains comparison between flow from POIs in positive and negative eta as requested by Flow PAG
* Vojtech Pacik (vojtech.pacik@cern.ch) - 30 May 2017
*/
#include <vector>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

TH1D* DoRatio(TH1D* up, TH1D* low);

void PlotComparisonAlex()
{
  TString sInputFile = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/comparisonAlex/UniFlow_wSDD.root");
  // TString sInputFileRecon = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/merged/UniFlow.root");

  TString sInputFileAlexSDD = TString("/Users/vpacik/NBI/Flow/results/Alex/hist_v2_def_V0ACal_16qSDD.root");
  TString sInputFileAlexFAST = TString("/Users/vpacik/NBI/Flow/results/Alex/hist_v2_def_V0ACal_16qFAST.root");

  TString sOutputFilePath = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/comparisonAlex2");

  TString sOutputFormat = TString("png");

  // TString sGap = "08";
  const Short_t iNumCentrality = 5;
  TString sCent[iNumCentrality] = {"0-20","20-40","40-60","60-80","80-100"};

  // ALICE Preferred colors and markers (from figure template)
  const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

  // colors setting
  Color_t colCharged = kBlack;
  Color_t colPion = kRed;
  Color_t colKaon = kBlue;
  Color_t colProton = kGreen+2;
  Color_t colPhi = kMagenta+1;
  Color_t colK0s = kCyan+1;
  Color_t colLambda = kOrange+1;

  // markers setting
  Int_t markCharged = kFullSquare;
  Int_t markPion = kFullCircle;
  Int_t markKaon = kFullTriangleUp;
  Int_t markProton = kFullCross;
  // Int_t markProton = kFullSquare;
  Int_t markPhi = kFullStar;
  Int_t markK0s = kFullTriangleDown;
  Int_t markLambda = kFullDiamond;
  // Int_t markLambda = kFullCross;

  Int_t markAlexV0A[] = {kOpenSquare, kOpenCircle, kOpenTriangleUp, kOpenCross};
  Int_t markAlexV0C[] = {kOpenSquare, kOpenCircle, kOpenTriangleUp, kOpenCross};
  // Int_t markAlexV0C[] = {kFullSquare, kFullCircle, kFullTriangleUp, kFullCross};

  Int_t markChargedHEP = kOpenSquare;
  Int_t markPionHEP = kOpenCircle;
  Int_t markKaonHEP = kOpenTriangleUp;
  Int_t markProtonHEP = kOpenCross;

  Double_t markSizeCharged = 1;
  Double_t markSizePion = 1;
  Double_t markSizeKaon = 1.2;
  Double_t markSizeProton = 1.5;
  Double_t markSizePhi = 1.5;
  Double_t markSizeK0s = 1.2;
  Double_t markSizeLambda = 1.6;
  // Double_t markSizeCharged = 2;
  // Double_t markSizePion = 2;
  // Double_t markSizeKaon = 2;
  // Double_t markSizeProton = 2;
  // Double_t markSizePhi = 2;
  // Double_t markSizeK0s = 2;
  // Double_t markSizeLambda = 2;


  // openning input file
  TFile* fInputFile = new TFile(sInputFile.Data(),"READ");
  if(!fInputFile->IsOpen()) { printf("No input file\n"); return; }
  fInputFile->ls();

  // Alex
  TFile* fInputFileAlexSDD = new TFile(sInputFileAlexSDD.Data(),"READ");
  if(!fInputFileAlexSDD->IsOpen()) { printf("No input file\n"); return; }
  fInputFileAlexSDD->ls();
  TFile* fInputFileAlexFAST = new TFile(sInputFileAlexFAST.Data(),"READ");
  if(!fInputFileAlexFAST->IsOpen()) { printf("No input file\n"); return; }
  fInputFileAlexFAST->ls();

  // LOADING HISTOGRAMS
  // mine
  std::vector<TH1D*> hMinePion;
  std::vector<TH1D*> hMineKaon;
  std::vector<TH1D*> hMineProton;
  std::vector<TH1D*> hMineIncl;

  // Alex
  std::vector<TH1D*> hAlexPionA;
  std::vector<TH1D*> hAlexKaonA;
  std::vector<TH1D*> hAlexProtonA;
  std::vector<TH1D*> hAlexInclA;
  std::vector<TH1D*> hAlexPionC;
  std::vector<TH1D*> hAlexKaonC;
  std::vector<TH1D*> hAlexProtonC;
  std::vector<TH1D*> hAlexInclC;

  TH1D* histTemp = 0x0;
  TH1D* histRatio = 0x0;

  std::vector<TH1D*> hMineInclRatioAlexA;
  std::vector<TH1D*> hMinePionRatioAlexA;
  std::vector<TH1D*> hMineKaonRatioAlexA;
  std::vector<TH1D*> hMineProtonRatioAlexA;

  for(Short_t cent(0); cent < iNumCentrality; cent++)
  {

    // mine
    histTemp = (TH1D*) fInputFile->Get(Form("hFlow2_Charged_harm2_gap08_cent%d",cent));
    histTemp->SetLineColor(colCharged);
    histTemp->SetMarkerColor(colCharged);
    histTemp->SetMarkerStyle(markCharged);
    hMineIncl.push_back(histTemp);

    histTemp = (TH1D*) fInputFile->Get(Form("hFlow2_Pion_harm2_gap08_cent%d",cent));
    histTemp->SetLineColor(colPion);
    histTemp->SetMarkerColor(colPion);
    histTemp->SetMarkerStyle(markPion);
    hMinePion.push_back(histTemp);

    histTemp = (TH1D*) fInputFile->Get(Form("hFlow2_Kaon_harm2_gap08_cent%d",cent));
    histTemp->SetLineColor(colKaon);
    histTemp->SetMarkerColor(colKaon);
    histTemp->SetMarkerStyle(markKaon);
    hMineKaon.push_back(histTemp);

    histTemp = (TH1D*) fInputFile->Get(Form("hFlow2_Proton_harm2_gap08_cent%d",cent));
    histTemp->SetLineColor(colProton);
    histTemp->SetMarkerColor(colProton);
    histTemp->SetMarkerStyle(markProton);
    hMineProton.push_back(histTemp);

    // Alex
    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2AllA_%d",cent));
    histTemp->SetLineColor(colCharged);
    histTemp->SetMarkerColor(colCharged);
    histTemp->SetMarkerStyle(markAlexV0A[0]);
    hAlexInclA.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2AllC_%d",cent));
    histTemp->SetLineColor(colCharged);
    histTemp->SetMarkerColor(colCharged);
    histTemp->SetMarkerStyle(markAlexV0C[0]);
    hAlexInclC.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2PiA_%d",cent));
    histTemp->SetLineColor(colPion);
    histTemp->SetMarkerColor(colPion);
    histTemp->SetMarkerStyle(markAlexV0A[1]);
    hAlexPionA.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2PiC_%d",cent));
    histTemp->SetLineColor(colPion);
    histTemp->SetMarkerColor(colPion);
    histTemp->SetMarkerStyle(markAlexV0C[1]);
    hAlexPionC.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2KA_%d",cent));
    histTemp->SetLineColor(colKaon);
    histTemp->SetMarkerColor(colKaon);
    histTemp->SetMarkerStyle(markAlexV0A[2]);
    hAlexKaonA.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2KC_%d",cent));
    histTemp->SetLineColor(colKaon);
    histTemp->SetMarkerColor(colKaon);
    histTemp->SetMarkerStyle(markAlexV0C[2]);
    hAlexKaonC.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2AntiPA_%d",cent));
    histTemp->SetLineColor(colProton);
    histTemp->SetMarkerColor(colProton);
    histTemp->SetMarkerStyle(markAlexV0A[3]);
    hAlexProtonA.push_back(histTemp);

    histTemp = (TH1D*) fInputFileAlexSDD->Get(Form("hv2AntiPC_%d",cent));
    histTemp->SetLineColor(colProton);
    histTemp->SetMarkerColor(colProton);
    histTemp->SetMarkerStyle(markAlexV0C[3]);
    hAlexProtonC.push_back(histTemp);
  }

  TLatex latex;
  latex.SetNDC();

  // making legend
  TLegend* legend = new TLegend(0.33, 0.62, 0.5, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0,0);
  legend->SetTextSize(gStyle->GetTextSize()*0.8);
  // legend->SetBorderSize(1);

  TLegend* legend2 = new TLegend(0.13, 0.62, 0.28, 0.88);
  legend2->SetBorderSize(0);
  legend2->SetFillColorAlpha(0,0);
  legend2->SetTextSize(gStyle->GetTextSize()*0.8);
  // legend2->SetBorderSize(1);

  legend->SetHeader("Mine\n( v2{2,|#eta| > 0.8})");
  legend->AddEntry(hMineIncl.at(0),"h^{#pm}","pl");
  legend->AddEntry(hMinePion.at(0),"#pi^{#pm}","pl");
  legend->AddEntry(hMineKaon.at(0),"K^{#pm}","pl");
  legend->AddEntry(hMineProton.at(0),"p/#bar{p}","pl");

  legend2->SetHeader("Alex (V0A-TPC)");
  legend2->AddEntry(hAlexInclA.at(0),"h^{#pm}","pl");
  legend2->AddEntry(hAlexPionA.at(0),"#pi^{#pm}","pl");
  legend2->AddEntry(hAlexKaonA.at(0),"K^{#pm}","pl");
  legend2->AddEntry(hAlexProtonA.at(0),"p/#bar{p}","pl");
  // legend->AddEntry(hFlowK0s,"K^{0}_{S}","pl");
  // legend2->AddEntry(hFlowPhi,"#phi","pl");
  // legend2->AddEntry(hFlowLambda,"#Lambda/#bar{#Lambda}","pl");

  // legend2->SetHeader("Run1 (SP)");
  // legend2->AddEntry(gHEP_Charged," ","fl");
  // legend2->AddEntry(gHEP_Pion," ","fl");
  // legend2->AddEntry(gHEP_Kaon," ","fl");
  // legend2->AddEntry(gHEP_Proton," ","fl");

  TH1D* hRatio = 0x0;

  TCanvas* canvas = new TCanvas("canvas","canvas");

  TCanvas* canRatio = new TCanvas("canRatio","canRatio");

  for(Short_t cent(0); cent < iNumCentrality; cent++)
  {
    canvas->cd();
    TH1* h = canvas->DrawFrame(0,0,4,0.5);
    // canvas->Clear();
    hAlexPionA.at(cent)->Draw("same hist p e1 x0 ");
    hAlexKaonA.at(cent)->Draw("same hist p e1 x0 ");
    hAlexProtonA.at(cent)->Draw("same hist p e1 x0 ");
    hAlexInclA.at(cent)->Draw("same hist p e1 x0 ");
    //
    // hAlexPionC.at(cent)->Draw("same hist p e1 x0 ");
    // hAlexKaonC.at(cent)->Draw("same hist p e1 x0 ");
    // hAlexProtonC.at(cent)->Draw("same hist p e1 x0  ");
    // hAlexInclC.at(cent)->Draw("same hist p e1 x0 ");

    hMineIncl.at(cent)->Draw("same hist p e1 x0 ");
    hMinePion.at(cent)->Draw("same hist p e1 x0 ");
    hMineKaon.at(cent)->Draw("same hist p e1 x0 ");
    hMineProton.at(cent)->Draw("same hist p e1 x0 ");
    legend->Draw();

    legend->Draw("same");
    legend2->Draw("same");
    latex.DrawLatex(0.15,0.53,Form("cent %s (V0A)",sCent[cent].Data()));
    // latex.DrawLatex(0.55,0.76,Form("Alex POIs V0C"));

    canvas->SaveAs(Form("%s/cent%d_A_low.pdf",sOutputFilePath.Data(),cent));

    // ratios
    histRatio = DoRatio(hMineIncl.at(cent),hAlexInclA.at(cent));
    hMineInclRatioAlexA.push_back(histRatio);
    histRatio = DoRatio(hMinePion.at(cent),hAlexPionA.at(cent));
    hMinePionRatioAlexA.push_back(histRatio);
    histRatio = DoRatio(hMineKaon.at(cent),hAlexKaonA.at(cent));
    hMineKaonRatioAlexA.push_back(histRatio);
    histRatio = DoRatio(hMineProton.at(cent),hAlexProtonA.at(cent));
    hMineProtonRatioAlexA.push_back(histRatio);

    canRatio->cd();
    TH1* hRat = canRatio->DrawFrame(0,0,4,10.);
    hMineInclRatioAlexA.at(cent)->Draw("same hist p e1 x0");
    canRatio->SaveAs(Form("%s/Ratio_Incl_cent%d_A_low.pdf",sOutputFilePath.Data(),cent));

  }
  return;
}
//_____________________________________________________________________________
TH1D* DoRatio(TH1D* up, TH1D* low)
{
  if(!up || !low) return 0x0;

  TH1D* hTemp = (TH1D*) up->Clone(Form("%s_ratio",up->GetName()));
  hTemp->Reset();

  const Short_t bins = hTemp->GetNbinsX();
  Double_t contUp = 0;
  Double_t contLow = 0;
  Double_t errUp = 0;
  Double_t errLow = 0;

  for(Short_t bin(1); bin < bins+1; bin++)
  {
    contUp = up->GetBinContent(bin);
    errLow = up->GetBinError(bin);
    contLow = low->GetBinContent(bin);
    errLow = low->GetBinError(bin);
    printf("low %g / up %g\n",contUp,contLow)
    hTemp->SetBinContent(bin, contUp/contLow);
  }

  return hTemp;
}
