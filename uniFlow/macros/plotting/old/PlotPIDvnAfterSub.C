#include "TFile.h"
#include "TH1.h"
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
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
void SetStyle(Bool_t graypalette=kFALSE);

// Color_t colCharged = kGray+1;
// Color_t colPion = kRed;
// Color_t colKaon = kBlue;
// Color_t colProton = kGreen+2;
// Color_t colPhi = kMagenta;
// Color_t colK0s = kCyan+1;
// Color_t colLambda = kOrange+1;
//
// Color_t filCharged = colCharged;
// Color_t filPhi = colPhi;
// Color_t filK0s = colK0s;
// Color_t filLambda = colLambda;
//
// Double_t dAlpha = 0.4;
//
// // markers setting
// Int_t markCharged = kFullCircle;
// Int_t markPion = kOpenCircle;
// Int_t markKaon = kOpenSquare;
// Int_t markProton = kOpenCross;
// // Int_t markProton = kFullSquare;
// Int_t markPhi = kFullStar;
// Int_t markK0s = kFullSquare;
// Int_t markLambda = kFullDiamond;
// // Int_t markLambda = kFullCross;

// TString sSpecies_list[] = {"Charged","Pion","Kaon","Proton"};
TString sSpecies_list[] = {"Charged","Pion","Kaon","K0s","Proton","Phi","Lambda"};
TString sSpecies_labels[] = {"h^{#pm}","#pi^{#pm}","K^{#pm}","K^{0}_{S}","p+#bar{p}","#phi","#Lambda+#bar{#Lambda}"};


// Color_t colors[] = {kGray+1,kRed,kBlue,kGreen+2,kMagenta,kCyan+1,kOrange+1};
// Int_t markers[] = {kOpenCircle, kFullCircle,kFullSquare,kFullCross, kFullStar,kFullSquare, kFullDiamond};
Color_t colors[] = {kBlack,kRed,kBlue,kCyan+1,kGreen+2,kMagenta,kOrange+1};
Int_t markers[] = {kOpenCircle, kFullStar,kFullTriangleUp,kFullTriangleDown,kFullSquare,kFullCross, kFullDiamond};

TString sGapVal= "0.0"; TString sGapName= "gap00";
// TString sGapVal= "0.4"; TString sGapName= "gap04";
// TString sGapVal= "0.8"; TString sGapName= "gap08";

void PlotPIDvnAfterSub()
{

  SetStyle();
  TString sInputPathTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3/"+sGapName+"/";
  // TString sInputSystTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out/figures_wBarlowCheck/";

  // TString sGap = "gap08";

  Int_t iNumSpecies = sizeof(sSpecies_list) / sizeof(sSpecies_list[0]);


  Int_t iNumCent = 5;
  TString sCent[] = {"0-10","10-20","20-40","40-60","60-100"};

  TFile* fileIn = OpenFile(Form("%s/Subtracted.root",sInputPathTop.Data())); if(!fileIn) {return;}

  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    TCanvas* canPlot = new TCanvas("canPlot","canPlot",800,800);
    canPlot->cd();
    TH1* frame_canPlot = (TH1*) gPad->DrawFrame(0,-0.02,7.0,0.30);
    frame_canPlot->SetTitle(Form("; p_{T} (GeV/c); v_{2}^{sub} {2, |#Delta#eta| > %s}",sGapVal.Data()));

    TLegend* leg = new TLegend(0.18,0.62,0.45,0.80);
    leg->SetFillColorAlpha(0,0);


    for(Int_t i(0); i < iNumSpecies; ++i)
    {
      // TList* list = (TList*) fileIn->Get("list_SubtPP_vn_ppint"); if(!list) { return; }
      // TH1D* htemp = (TH1D*) f

      // TFile* fileSyst = TFile::Open(Form("%s/syst_%s.root",sInputSystTop.Data(),sSpecies[iSpec].Data()),"RECREATE");
      // if(!fileSyst) { printf("nofilesyst\n"); return;}

      // fileSyst->ls();
      // TH1D* histSyst = (TH1D*) fileSyst->Get("hi");

      TH1D* htemp = (TH1D*) fileIn->Get(Form("hFlow2_%s_harm2_%s_cent%d",sSpecies_list[i].Data(),sGapName.Data(),iCent));
      if(!htemp) {return;}

      canPlot->cd(1);
      StyleHist(htemp,colors[i], markers[i]);
      htemp->SetFillColor(colors[i]);
      htemp->DrawCopy("hist p e same");
      leg->AddEntry(htemp,sSpecies_labels[i].Data(),"p");

    }
    canPlot->cd();
    TLatex * text = new TLatex();
    TLatex * text2 = new TLatex();
    text->SetTextSizePixels(20);
    text->SetTextFont(42);
    text2->SetTextSizePixels(20);

    leg->Draw("same");
    text->DrawLatexNDC(0.186,0.85,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text2->DrawLatexNDC(0.186,0.81,Form("Multiplicity class %s%% (V0A)",sCent[iCent].Data()));

    TString sOutputPath = Form("%s/results/",sInputPathTop.Data());
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
