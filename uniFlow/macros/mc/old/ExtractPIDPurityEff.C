// =================================================================================================
// ExtractPIDPurityEff.C
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
//
// Estimating purity & PID efficiency by making ratios of pt spectra as implemented in AliAnalysisTaskUniFlow.
// This is done for pions, kaons and protons.
//
// =================================================================================================

#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TColor.h"
#include "TLine.h"

void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);
void SetCustomPalette();

Color_t colRecoAll = kRed;
Color_t colRecoSelectedTrue = kBlue+2;
Color_t colRecoSelected = kGreen+1;
Color_t colGen = kViolet-1;

void ExtractPIDPurityEff()
{
  TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/mc-pid/pPb-LHC17f2a_cent_woSDD_fix";
  TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/mc-pid/pPb-LHC17f2a_cent_woSDD_fix";

  const Int_t iNumTasks = 4;
  TString sTasksName[iNumTasks] = {"UniFlow_2sigma", "UniFlow_3sigma", "UniFlow_5sigma", "UniFlow_bayes"};
  TString sTasksLabel[iNumTasks] = {"TPC+TOF (2#sigma)","TPC+TOF (3#sigma)","TPC+TOF (5#sigma)","Bayesian (>0.9)"};

  const Int_t iNumSpecies = 3;
  TString sSpeciesName[iNumSpecies] = {"Pion","Kaon","Proton"};
  Color_t colors[] = {kRed, kBlue+2, kGreen+1, kViolet-1};

  // ==================================================================================================================

  gSystem->mkdir(sOutputPath.Data(),kTRUE);
  TFile* fileOut = TFile::Open(Form("%s/outPIDPurityEff.root",sOutputPath.Data()), "RECREATE");
  if(!fileOut) { printf("ERROR: Input file not found!\n"); return; }

  TFile* fileIn = TFile::Open(Form("%s/AnalysisResults.root",sInputPath.Data()), "READ");
  if(!fileIn) { printf("ERROR: Input file not found!\n"); return; }


  for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
  {
    TString sTaskName = sTasksName[iTask];
    TString sTaskLabel = sTasksLabel[iTask];

    fileIn->cd(sTaskName.Data());
    TList* listInPID = (TList*) gDirectory->Get(Form("QA_PID_%s",sTaskName.Data()));
    if(!listInPID) { printf("ERROR: Input 'listInPID' not found"); gDirectory->ls(); return; }

    // preparing canvases
    TCanvas* can = new TCanvas("can","can",2400,800);
    can->Divide(3,1);

    TCanvas* canPurEff = new TCanvas("canPurEff","canPurEff",1600,800);
    canPurEff->Divide(2,1);
    canPurEff->cd(1);
    TH1* frame_canPurEff_1 = (TH1*) gPad->DrawFrame(0.0,0.0,10.0,1.5);
    frame_canPurEff_1->SetTitle("PID Purity; p_{T} (GeV/c); Purity");
    canPurEff->cd(2);
    TH1* frame_canPurEff_2 = (TH1*) gPad->DrawFrame(0.0,0.0,10.0,1.5);
    frame_canPurEff_2->SetTitle("PID Efficiency (+ Contamination); p_{T} (GeV/c); Efficiency");

    TLegend* leg_canPurEff = new TLegend(0.12,0.12,0.4,0.32);
    leg_canPurEff->SetBorderSize(0);
    leg_canPurEff->SetFillColorAlpha(0,0);
    leg_canPurEff->SetTextFont(42);
    leg_canPurEff->SetTextSize(0.04);
    leg_canPurEff->SetHeader(sTaskLabel.Data());


    // potential loop over species
    for(Int_t iSpecies(0); iSpecies < iNumSpecies; ++iSpecies)
    {
      TString sSpecies = sSpeciesName[iSpecies];

      // loading histos
      TH1D* fhMCRecoSelectedPt = (TH1D*) listInPID->FindObject(Form("fhMCRecoSelected%sPt",sSpecies.Data()));
      if(!fhMCRecoSelectedPt) { printf("ERROR:'fhMCRecoSelected%sPt' not found!\n",sSpecies.Data()); listInPID->ls(); return; }
      TH1D* fhMCRecoSelectedTruePt = (TH1D*) listInPID->FindObject(Form("fhMCRecoSelectedTrue%sPt",sSpecies.Data()));
      if(!fhMCRecoSelectedTruePt) { printf("ERROR:'fhMCRecoSelectedTrue%sPt' not found!\n",sSpecies.Data()); listInPID->ls(); return; }
      TH1D* fhMCRecoAllPt = (TH1D*) listInPID->FindObject(Form("fhMCRecoAll%sPt",sSpecies.Data()));
      if(!fhMCRecoAllPt) { printf("ERROR:'fhMCRecoAll%sPt' not found!\n",sSpecies.Data()); listInPID->ls(); return; }
      // TH1D* fhMCGenAllPt = 0x0;
      // if(!fhMCGenAllPt) { printf("ERROR:'fhMCGenAll%sPt' not found!\n",sSpecies.Data()); listInPID->ls(); return; }

      // estimating PID purity = fhMCRecoSelectedTruePt / fhMCRecoSelectedPt ( MC true selected / all selected )
      TH1D* hPurity = DivideHistos(fhMCRecoSelectedTruePt,fhMCRecoSelectedPt, kTRUE);

      // estimating PID efficiency = fhMCRecoSelectedPt / fhMCRecoAllPt ( all selected / all (reco, MC true) )
      TH1D* hPIDEff = DivideHistos(fhMCRecoSelectedPt,fhMCRecoAllPt, kTRUE);

      StyleHist(fhMCRecoSelectedPt,colRecoSelected,kDot);
      StyleHist(fhMCRecoSelectedTruePt,colRecoSelectedTrue,kDot);
      StyleHist(fhMCRecoAllPt,colRecoAll,kDot);
      // StyleHist(fhMCGenAllPt,colGen,kDot);

      StyleHist(hPurity,colors[iSpecies],kOpenCircle);
      StyleHist(hPIDEff,colors[iSpecies],kOpenCircle);

      TLegend* leg_can = new TLegend(0.55,0.7,0.89,0.89);
      leg_can->SetBorderSize(0);
      leg_can->SetFillColorAlpha(0,0);
      leg_can->SetTextFont(42);
      leg_can->SetTextSize(0.04);
      leg_can->SetHeader(Form("%s, %s", sSpecies.Data(),sTaskLabel.Data()));
      leg_can->AddEntry(fhMCRecoAllPt,"All","pel");
      leg_can->AddEntry(fhMCRecoSelectedPt,"All selected","pel");
      leg_can->AddEntry(fhMCRecoSelectedTruePt,"MC true selected","pel");

      can->cd(iSpecies+1);
      TH1* frame_can = (TH1*) gPad->DrawFrame(0.0,0.1,10.0,1e9);
      gPad->SetLogy();
      fhMCRecoAllPt->DrawCopy("same");
      fhMCRecoSelectedTruePt->DrawCopy("same");
      fhMCRecoSelectedPt->DrawCopy("same");
      leg_can->Draw();

      leg_canPurEff->AddEntry(hPurity,sSpecies.Data(),"pel");

      canPurEff->cd(1);
      hPurity->DrawCopy("same");
      canPurEff->cd(2);
      hPIDEff->DrawCopy("same");

      // saving to output ROOT file
      fileOut->cd();
      hPurity->Write(Form("hPIDPurity_%s_%s",sTaskName.Data(),sSpecies.Data()));
      hPIDEff->Write(Form("hPIDEfficiency_%s_%s",sTaskName.Data(),sSpecies.Data()));
    }

    canPurEff->cd(1);
    leg_canPurEff->Draw();
    canPurEff->cd(2);
    leg_canPurEff->Draw();

    can->SaveAs(Form("%s/distMCpt_%s.pdf",sOutputPath.Data(),sTaskName.Data()),"pdf");
    canPurEff->SaveAs(Form("%s/PurEff_%s.pdf",sOutputPath.Data(),sTaskName.Data()),"pdf");

  }

  return;
}
// ==================================================================================================================
void StyleHist(TH1* hist, Color_t color, Style_t markerStyle)
{
  if(!hist) { printf("ERROR-DrawHist: Hist does not found.\n"); return; }
  hist->SetLineColor(color);
  // hist->SetLineStyle(color);
  // hist->SetLineStyle(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  return;
};
// ==================================================================================================================
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  ratio->Reset();

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < binsDenom+1; iBin++)
  {
    if(iBin > binsNom) break;

    dContNom = nom->GetBinContent(iBin);
    dErrNom = nom->GetBinError(iBin);
    dContDenom = denom->GetBinContent(iBin);
    dErrDenom = denom->GetBinError(iBin);

    if(dContDenom == 0.0) continue;

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBin,dContRatio);
    ratio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
// ==================================================================================================================
void SetCustomPalette()
{
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};


  // Rain Bow
  // case 55:
  // Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
  // Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
  // Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};

  // Bird
  //case 57:
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

  // Blue Green Yellow
  //case 71:
  // Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
  // Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
  // Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};

  // Solar
  // case 100:
  // Double_t red[9]   = { 99./255., 116./255., 154./255., 174./255., 200./255., 196./255., 201./255., 201./255., 230./255.};
  // Double_t green[9] = {  0./255.,   0./255.,   8./255.,  32./255.,  58./255.,  83./255., 119./255., 136./255., 173./255.};
  // Double_t blue[9]  = {  5./255.,   6./255.,   7./255.,   9./255.,   9./255.,  14./255.,  17./255.,  19./255.,  24./255.};

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
