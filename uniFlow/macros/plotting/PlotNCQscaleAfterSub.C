/* Code for plotting NCQ scaling results of vn(pt)
*
*
*  Author: Vojtech Pacik, NBI (vojtech.pacik@cern.ch)
*/

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
TGraphErrors* ScaleByNCQ(TGraphErrors* points, Bool_t bScaleX, Bool_t bMeson, Double_t dMass = -1.0);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
void SetStyle(Bool_t graypalette=kFALSE);

TString sInputPathTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3/gap04/";
TString sOutputPath = Form("%s/results-ncq/",sInputPathTop.Data());

// // Original
TString sSpecies_labels[] = {"#pi^{#pm}","K^{#pm}","K^{0}_{S}","p(#bar{p})","#varphi","#Lambda(#bar{#Lambda})"};
TString sSpecies_list[] = {"Pion","Kaon","K0s","Proton","Phi","Lambda"};
Bool_t bMeson_list[] = {1,1,1,0,1,0};
Int_t iPDG_codes[] = {211, 321, 310, 2212, 333, 3122};

// W/O Phi
// TString sSpecies_labels[] = {"#pi^{#pm}","K^{#pm}","K^{0}_{S}","p(#bar{p})","#Lambda(#bar{#Lambda})"};
// TString sSpecies_list[] = {"Pion","Kaon","K0s","Proton","Lambda"};
// Bool_t bMeson_list[] = {1,1,1,0,0};
// Int_t iPDG_codes[] = {211, 321, 310, 2212, 3122};

TString sCent[] = {"0-10","10-20","20-40","40-60"};


Int_t iNumCent = sizeof(sCent) / sizeof(sCent[0]);
Int_t iNumSpecies = sizeof(sSpecies_list) / sizeof(sSpecies_list[0]);

void PlotNCQscaleAfterSub()
{
  SetStyle();

  TFile* fileIn = OpenFile(Form("%s/v2-syst-final.root",sInputPathTop.Data())); if(!fileIn) { printf("No fileIn\n"); return; }
  // fileIn->ls();

  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    TCanvas* canPlot = new TCanvas("canPlot","canPlot",1500,500);
    canPlot->Divide(3,1);
    canPlot->cd(1);
    TH1* frame_canPlot_1 = (TH1*) gPad->DrawFrame(0,-0.02,7.0,0.15);
    frame_canPlot_1->SetTitle(Form("; p_{T} (GeV/c); v_{2}^{sub} {2, |#Delta#eta| > 0.4} / NCQ" ));
    canPlot->cd(2);
    TH1* frame_canPlot_2 = (TH1*) gPad->DrawFrame(0,-0.02,4.0,0.15);
    frame_canPlot_2->SetTitle(Form("; p_{T} / NCQ (GeV/c); v_{2}^{sub} {2, |#Delta#eta| > 0.4} / NCQ"));
    canPlot->cd(3);
    TH1* frame_canPlot_3 = (TH1*) gPad->DrawFrame(0,-0.02,4.0,0.15);
    frame_canPlot_3->SetTitle(Form("; (m_{T}-m_{0}) / NCQ (GeV/c); v_{2}^{sub} {2, |#Delta#eta| > 0.4} / NCQ"));

    TLegend* leg = new TLegend(0.18,0.62,0.45,0.80);
    leg->SetFillColorAlpha(0,0);
    leg->SetBorderSize(0);


    for(Int_t i(0); i < iNumSpecies; ++i)
    {
      TString sSpecies = sSpecies_list[i];
      Bool_t bMeson = bMeson_list[i];
      Double_t dMass = TDatabasePDG::Instance()->GetParticle(iPDG_codes[i])->Mass();

      // loading
      TGraphErrors* histPoints = (TGraphErrors*) fileIn->Get(Form("graphPoints_%s_cent%d",sSpecies.Data(),iCent));
      if(!histPoints) { printf("No histPoints %s cent%d\n",sSpecies.Data(), iCent); return; }

      TGraphErrors* histPoints_NCQ_v2only = ScaleByNCQ(histPoints,kFALSE,bMeson);
      if(!histPoints_NCQ_v2only) { printf("No histPoints_NCQ_v2only %s cent%d\n",sSpecies.Data(), iCent); return; }

      TGraphErrors* histPoints_NCQ = ScaleByNCQ(histPoints,kTRUE,bMeson);
      if(!histPoints_NCQ) { printf("No histPoints_NCQ\n"); return; }

      TGraphErrors* histPoints_NCQ_kt = ScaleByNCQ(histPoints,kTRUE,bMeson,dMass);
      if(!histPoints_NCQ_kt) { printf("No histPoints_NCQ_kt\n"); return; }

      leg->AddEntry(histPoints,sSpecies.Data(),"p");

      canPlot->cd(1);
      histPoints_NCQ_v2only->Draw("same pZ");

      canPlot->cd(2);
      histPoints_NCQ->Draw("same Pz");

      canPlot->cd(3);
      histPoints_NCQ_kt->Draw("same Pz");
    }

    canPlot->cd(1);
    TLatex * text = new TLatex();
    TLatex * text2 = new TLatex();
    text->SetTextSizePixels(20);
    text->SetTextFont(42);
    text2->SetTextSizePixels(20);

    leg->Draw("same");
    text->DrawLatexNDC(0.186,0.85,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text2->DrawLatexNDC(0.186,0.81,Form("Multiplicity class %s%% (V0A)",sCent[iCent].Data()));

    gSystem->mkdir(sOutputPath.Data(),kTRUE);
    canPlot->SaveAs(Form("%s/PIDvn_cent%d.pdf",sOutputPath.Data(),iCent),"pdf");
  }


  return;
}
// ==================================================================================================================
TFile* OpenFile(TString sFileName, TString sMode)
{
  TFile* file = TFile::Open(sFileName.Data(),sMode.Data());
  if(!file) { printf("ERROR: Input file '%s' not found.\n",sFileName.Data()); return 0x0; }

  return file;
}
// ==================================================================================================================
TGraphErrors* ScaleByNCQ(TGraphErrors* points, Bool_t bScaleX, Bool_t bMeson, Double_t dMass)
{
  if(!points) { printf("ERROR-ScaleByNCQ : No points\n"); return 0x0; }

  Double_t dNCQ = 3.0;
  if(bMeson) dNCQ = 2.0;

  Bool_t bKtScale = (dMass > 0.0);
  // printf("bKt %d\n",bKtScale);

  TString sName = "NCQ";
  if(!bScaleX) { sName += "_v2only"; }
  if(bKtScale) { sName += "_kt"; }

  TGraphErrors* points_NCQ = (TGraphErrors*) points->Clone(Form("%s_%s",points->GetName(),sName.Data()));
  for(Int_t i(0); i < points_NCQ->GetN(); ++i)
  {
    Double_t dX = 0.0;
    Double_t dY = 0.0;
    points_NCQ->GetPoint(i,dX,dY);

    Double_t dXErr = points_NCQ->GetErrorX(i);
    Double_t dYErr = points_NCQ->GetErrorY(i);
    // printf("[%f +- %f;%f+-%f]\n",dX,dXErr,dY,dYErr);

    dY = dY / dNCQ;
    dYErr = dYErr / dNCQ;

    if(bKtScale)
    {
      Double_t dKt = TMath::Sqrt(dX*dX + dMass*dMass) - dMass;
      dX = dKt;
    }

    if(bScaleX) { dX = dX / dNCQ; }

    points_NCQ->SetPoint(i,dX, dY);
    points_NCQ->SetPointError(i,0.0,dYErr);
  }

  return points_NCQ;
}
// ==================================================================================================================
TH1D* LoadHisto(TString sHistName, TFile* file)
{
  if(!file) { printf("ERROR-LoadHisto: File does not found.\n"); return 0x0; }

  TH1D* hist = (TH1D*) file->Get(sHistName.Data());
  if(!hist) { printf("ERROR-LoadHisto: Histo '%s' not found\n",sHistName.Data()); file->ls(); return 0x0; }

  return hist;
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
//_____________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  // cout << "Setting style!" << std::end;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
//_____________________________________________________________________________
