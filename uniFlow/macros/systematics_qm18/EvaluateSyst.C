#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"

void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
void SetCustomPalette();
TH1D* ProcessList(TList* listDiff, TList* listBarlow, Double_t dCutBarlow);

Double_t dCutBarlow = 1.0;

TString sTopInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-test/";
TString sOutPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-test/figures_wBarlowCheck/";
// TString sOutPath = "./out/figures_noFB/";
// TString sOutPath = "./out/figures_wBarlowCheck_noFB/";

// TString sSpecies[] = {"Charged"}; TString sFiles[] = {"tracking/cls", "tracking/PV"};
// TString sSpecies[] = {"Pion","Kaon","Proton","Phi"}; TString sFiles[] = {"tracking/cls", "tracking/PV", "pid/3sigma", "pid/bayes90"};
// TString sSpecies[] = {"K0s","Lambda"}; TString sFiles[] = {"tracking/cls", "tracking/PV","v0s/3sigma", "v0s/CPA", "v0s/decayRad", "v0s/DCAdaughters"};

TString sSpecies[] = {"Charged"}; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB"};
// TString sSpecies[] = {"Pion","Kaon","Proton","Phi"}; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB", "pid/3sigma", "pid/bayes90"};
// TString sSpecies[] = {"K0s","Lambda"}; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB","v0s/3sigma", "v0s/CPA", "v0s/decayRad", "v0s/DCAdaughters"};

// // TString sSpecies[] = {"Charged"}; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB"};
// // TString sSpecies[] = {"Pion","Kaon","Proton","Phi"}; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB", "pid/3sigma", "pid/bayes90"};
// TString sSpecies = "Kaon"; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB", "pid/3sigma", "pid/bayes90"};
// // TString sSpecies = "Proton"; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB", "pid/3sigma", "pid/bayes90"};
// // TString sSpecies = "K0s"; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB","v0s/3sigma", "v0s/CPA", "v0s/decayRad", "v0s/DCAdaughters"};
// // TString sSpecies = "Lambda"; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB","v0s/3sigma", "v0s/CPA", "v0s/decayRad", "v0s/DCAdaughters"};
// // TString sSpecies = "Phi"; TString sFiles[] = {"tracking/cls", "tracking/PV", "tracking/FB", "pid/3sigma", "pid/bayes90"};
Int_t iNumFiles = sizeof(sFiles)/sizeof(sFiles[0]);
Int_t iNumSpecies = sizeof(sSpecies)/sizeof(sSpecies[0]);

TString sCentLabel[] = {"0-10%","10-20%","20-40%","40-60%"};
Int_t iNumCent = sizeof(sCentLabel)/sizeof(sCentLabel[0]);

Color_t colors[] = {kBlack,kRed, kBlue, kGreen+2, kViolet-1, kOrange+1, kMagenta-1, kOrange-4};

void EvaluateSyst()
{
  // SetCustomPalette();

  gSystem->mkdir(sOutPath.Data(),kTRUE);


  for(Int_t iSpec(0); iSpec < iNumSpecies; ++iSpec)
  {

    TFile* fileOut = TFile::Open(Form("%s/syst_%s.root",sOutPath.Data(),sSpecies[iSpec].Data()),"RECREATE");
    if(!fileOut) { printf("No fileOut \n"); return;}

    for(Int_t iCent(0); iCent < iNumCent; ++iCent)
    {

      TList* listDiff = new TList();
      TList* listBarlow = new TList();

      TLegend* leg = new TLegend(0.3,0.3,0.7,0.7);
      leg->SetBorderSize(0);
      leg->SetFillColorAlpha(0,0);
      leg->SetHeader(Form("%s (%s)",sSpecies[iSpec].Data(),sCentLabel[iCent].Data()));

      for(Int_t iFile(0); iFile < iNumFiles; ++iFile)
      {
        TString sInputFile = sFiles[iFile];
        TFile* fileInput = TFile::Open(Form("%s%s.root",sTopInputPath.Data(),sInputFile.Data()),"READ");
        if(!fileInput) { printf("ERROR : fileInput not found\n"); return; }
        // fileInput->ls();

        TList* list = (TList*) fileInput->Get(Form("%s_cent%d",sSpecies[iSpec].Data(),iCent));
        if(!list) { printf("ERROR : list not found\n"); return; }

        // list->ls();

        TH1D* histBarlow = (TH1D*) list->FindObject("histBarlow");
        if(!histBarlow) { printf("ERROR : histBarlow not found\n"); return; }
        TH1D* histDiff = (TH1D*) list->FindObject("histDiff");
        if(!histDiff) { printf("ERROR : histDiff not found\n"); return; }

        histBarlow->SetMaximum(10.0);

        StyleHist(histDiff, colors[iFile+1], kOpenCircle);
        StyleHist(histBarlow, colors[iFile+1], kOpenCircle);

        leg->AddEntry(histDiff,sInputFile.Data(),"pl");

        listDiff->Add(histDiff);
        listBarlow->Add(histBarlow);
      }

      // Int_t nPnt  = listDiff->GetEntries() ; // number of points
      // Int_t nnCol = gStyle->GetNumberOfColors();

      // processing lists
      TH1D* hSum = ProcessList(listDiff,listBarlow,dCutBarlow);
      StyleHist(hSum,colors[0],kDot);


      if(hSum->GetMaximum() > 1.0) { hSum->SetMaximum(1.0); }
      // hSum->SetMaximum(1.0);
      leg->AddEntry(hSum,"Sum in quad.","l");

      TCanvas* can = new TCanvas("can","can",1500,500);
      can->Divide(3,1);
      can->cd(1);
      leg->Draw();
      can->cd(2);
      hSum->Draw();
      for(Int_t i(0); i < listDiff->GetEntries(); ++i) { ((TH1D*) listDiff->At(i))->DrawCopy("same"); }

      can->cd(3);
      ((TH1D*) listBarlow->At(0))->Draw("hist");
      for(Int_t i(1); i < listBarlow->GetEntries(); ++i) { ((TH1D*) listBarlow->At(i))->DrawCopy("same  hist"); }


      can->SaveAs(Form("%s/syst_%s_cent%d.pdf",sOutPath.Data(),sSpecies[iSpec].Data(), iCent));
      fileOut->cd();
      hSum->Write(Form("hSum_cent%d",iCent));
    }
  }

  return;
}
// ==================================================================================================================
TH1D* ProcessList(TList* listDiff, TList* listBarlow, Double_t dCut)
{
  if(!listDiff) { printf("ERROR-ProcessList : listDiff not found!\n"); return 0x0; }
  if(!listBarlow) { printf("ERROR-ProcessList : listBarlow not found!\n"); return 0x0; }

  TH1D* hSum = (TH1D*) listDiff->At(0)->Clone("hSum");

  for(Int_t iEntry(1); iEntry < listDiff->GetEntries(); ++iEntry)
  {
    // 0 already done;
    TH1D* hNew = (TH1D*) listDiff->At(iEntry);

    TH1D* hBarlow = (TH1D*) listBarlow->At(iEntry);

    for(Int_t bin(1); bin < hSum->GetNbinsX()+1; ++bin)
    {
      Double_t dBarlow = hBarlow->GetBinContent(bin);

      if(dCut > 0.0 && dBarlow < dCut) continue;

      Double_t dOldCont = hSum->GetBinContent(bin);
      // Double_t dOldErr = hSum->GetBinError(bin);

      Double_t dNewCont = hNew->GetBinContent(bin);
      // Double_t dNewErr = 0.0;

      hSum->SetBinContent(bin,TMath::Sqrt(dOldCont*dOldCont + dNewCont*dNewCont));
      hSum->SetBinError(bin,0.0);
      // hSum->SetBinError(bin,TMath::Sqrt(dNewErrSq));

    }
  }

  return hSum;
}
// ==================================================================================================================
void StyleHist(TH1* hist, Color_t color, Style_t markerStyle, Bool_t showStats)
{
  if(!hist) { printf("ERROR-DrawHist: Hist does not found.\n"); return; }
  hist->SetStats(showStats);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  return;
};
// ==================================================================================================================
void SetCustomPalette()
{
  // # Usage for getting the color out of pallete
  // Int_t nPnt  = 2 ; // number of points
  // Int_t nnCol = gStyle->GetNumberOfColors();
  // Int_t iLegInx = 0;
  // for(Int_t iMult(iMultMin); iMult < 2; ++iMult)
  // {
  //   Int_t idx = loop * Float_t(nnCol-1) / (nPnt-1);
  //   StyleHist(histUnit, gStyle->GetColorPalette(idx), kOpenCircle);
  // }
  // ==================================================================================================================

  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};


  // Rain Bow
  // case 55:
  // Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
  // Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
  // Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};

  // Bird
  //case 57:
  // Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  // Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  // Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

  // Blue Green Yellow
  //case 71:
  // Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
  // Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
  // Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};

  // Solar
  // case 100:
  Double_t red[9]   = { 99./255., 116./255., 154./255., 174./255., 200./255., 196./255., 201./255., 201./255., 230./255.};
  Double_t green[9] = {  0./255.,   0./255.,   8./255.,  32./255.,  58./255.,  83./255., 119./255., 136./255., 173./255.};
  Double_t blue[9]  = {  5./255.,   6./255.,   7./255.,   9./255.,   9./255.,  14./255.,  17./255.,  19./255.,  24./255.};

  // Viridis
  // case 112:
  // Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
  // Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
  // Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};

  // Cividis
  // case 113:
  // Double_t red[9]   = {  0./255.,   5./255.,  65./255.,  97./255., 124./255., 156./255., 189./255., 224./255., 255./255.};
  // Double_t green[9] = { 32./255.,  54./255.,  77./255., 100./255., 123./255., 148./255., 175./255., 203./255., 234./255.};
  // Double_t blue[9]  = { 77./255., 110./255., 107./255., 111./255., 120./255., 119./255., 111./255.,  94./255.,  70./255.};

  Int_t pal = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
  const Int_t nCol = 255;
  Int_t colors[nCol];
  for (int i=0; i<nCol; i++) colors[i] =pal+i;

  gStyle->SetPalette(nCol,colors);
}
