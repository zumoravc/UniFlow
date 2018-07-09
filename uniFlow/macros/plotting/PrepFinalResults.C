/* PrepareFinalResults.case
*
* Making TGraphErrors with final results including systematics out of individual components.
*
*
*/

#include "TH1.h"
#include "TFile.h"
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
TH1D* PrepareSyst(TH1D* res, TH1D* syst);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
void SetStyle(Bool_t graypalette=kFALSE);

// regular order
TString sSpecies_labels[] = {"h^{#pm}","#pi^{#pm}","K^{#pm}","K^{0}_{S}","p(#bar{p})","#varphi","#Lambda(#bar{#Lambda})"};
TString sSpecies_list[] = {"Charged","Pion","Kaon","K0s","Proton","Phi","Lambda"};
// Color_t colors[] = {kBlack, kRed, kBlue, kCyan+1, kGreen+2, kMagenta, kOrange+1};
Color_t colors[] = {kBlack, kRed, kBlue, kCyan+2, kGreen+2, kMagenta+1, kOrange-1};
Color_t colorsFill[] = {kBlack, kRed, kBlue, kCyan+2, kGreen+2, kMagenta+1, kOrange-1};
// Color_t colorsFill[] = {kGray, kRed-9, kBlue-9, kCyan-8, kGreen-9, kMagenta-9, kOrange-9};
Int_t markers[] = {kOpenCircle, kFullStar, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullCross, kFullDiamond};
Double_t markersSize[] = {1.0,1.2,1.1,1.1,1.0,1.2,1.2};
Double_t dAlpha = 0.35;
// Double_t dAlpha = 1.0;

// // reverse order
// TString sSpecies_list[] = {"Phi","Proton","Lambda","K0s","Charged","Kaon","Pion"};
// Color_t colors[] = {kMagenta,kGreen+2,kOrange+1,kCyan+1,kBlack,kBlue,kRed};
// Int_t markers[] = {kFullCross,kFullSquare,kFullDiamond,kFullTriangleDown,kOpenCircle,kFullTriangleUp,kFullStar};


// TString sGapVal= "0.0"; TString sGapName= "gap00";
TString sGapVal= "0.4"; TString sGapName= "gap04";
// TString sGapVal= "0.8"; TString sGapName= "gap08";

TString sCent[] = {"0-10","10-20","20-40","40-60"};

Double_t dPtBinErr = 0.0; // vertical stat err
Double_t dPtBinSyst = 0.05; // width of syst err box

Int_t iNumCent = sizeof(sCent) / sizeof(sCent[0]);
Int_t iNumSpecies = sizeof(sSpecies_list) / sizeof(sSpecies_list[0]);

void PrepFinalResults()
{
  TString sInputPathTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3/"+sGapName+"/";
  TString sInputSystTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-binning-3/final/";

  TString sOutputPath = Form("%s/results-test-2/",sInputPathTop.Data());
  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  TFile* fileOut = TFile::Open(Form("%s/PIDv2-pPb-subt-syst.root",sInputPathTop.Data()),"RECREATE");
  if(!fileOut) { printf("No fileOut\n"); return; }

  TFile* fileIn = OpenFile(Form("%s/Subtracted.root",sInputPathTop.Data())); if(!fileIn) {return;}

  // Int_t iCent(0);
  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    // Int_t spec(0);
    for(Int_t spec(0); spec < iNumSpecies; ++spec)
    {
      TString sSpecies = sSpecies_list[spec];

      // loading systematics
      TFile* fileSyst = TFile::Open(Form("%s/syst_%s.root",sInputSystTop.Data(),sSpecies.Data()),"READ");
      if(!fileSyst) { printf("No fileSyst\n"); return; }

      TH1D* histSyst = (TH1D*) fileSyst->Get(Form("syst_cent%d",iCent));
      if(!histSyst) { printf("No histSyst \n"); fileSyst->ls(); return; }

      // loading central points
      TH1D* histPoints = (TH1D*) fileIn->Get(Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGapName.Data(),iCent));
      if(!histPoints) { printf("No histPoints\n"); fileIn->ls(); return; }

      // Preparing TGraphErrors
      Double_t dPt[30] = {0.0};
      Double_t dPtErr[30] = {0.0};
      Double_t dPtSyst[30] = {0.0};

      Double_t dV2[30] = {0.0};
      Double_t dV2Err[30] = {0.0};
      Double_t dV2Syst[30] = {0.0};

      Int_t iNumPoints = histPoints->GetNbinsX();
      if(iNumPoints > 30) { printf("iNumPoints: %d > 30\n",iNumPoints); return; }

      for(Int_t iBin(1); iBin < iNumPoints+1; ++iBin)
      {
        Double_t dBinPt = histPoints->GetBinCenter(iBin);
        // if(dBinPt > 3.6 && spec > 0) { dBinPt += 0.1*spec;}
        Double_t dBinPtErr = histPoints->GetBinWidth(iBin);

        Double_t dBinV2 = histPoints->GetBinContent(iBin);
        Double_t dBinV2Err = histPoints->GetBinError(iBin);
        Double_t dBinSystRel = histSyst->GetBinContent(iBin);
        Double_t dBinV2Syst = dBinSystRel*dBinV2;

        printf("v2(%f) = %f +- %f +- %f (%f%%)\n", dBinPt, dBinV2, dBinV2Err, dBinV2Syst, 100*dBinSystRel);

        dPt[iBin-1] = dBinPt;

        if(dPtBinErr < 0.0) { dPtErr[iBin-1] = dBinPtErr; }
        else { dPtErr[iBin-1] = dPtBinErr; }
        dPtSyst[iBin-1] = dPtBinSyst;

        dV2[iBin-1] = dBinV2;
        dV2Err[iBin-1] = dBinV2Err;
        dV2Syst[iBin-1] = dBinV2Syst;
      }

      TGraphErrors* graphPoints = new TGraphErrors(iNumPoints, dPt, dV2, dPtErr, dV2Err);
      graphPoints->SetName(Form("graphPoints_%s_cent%d",sSpecies.Data(),iCent));
      TGraphErrors* graphSyst = new TGraphErrors(iNumPoints, dPt, dV2, dPtSyst, dV2Syst);
      graphSyst->SetName(Form("graphSyst_%s_cent%d",sSpecies.Data(),iCent));

      graphPoints->SetMarkerStyle(markers[spec]);
      graphPoints->SetMarkerColor(colors[spec]);
      graphPoints->SetMarkerSize(markersSize[spec]);
      graphPoints->SetLineColor(colors[spec]);
      graphPoints->SetLineWidth(2);
      graphPoints->SetLineStyle(1);

      graphSyst->SetMarkerStyle(markers[spec]);
      graphSyst->SetMarkerColor(colors[spec]);
      graphSyst->SetMarkerSize(markersSize[spec]);
      graphSyst->SetLineColor(colors[spec]);
      graphSyst->SetLineWidth(2);
      graphSyst->SetLineStyle(1);
      graphSyst->SetFillStyle(1001);
      graphSyst->SetFillColorAlpha(colorsFill[spec],dAlpha);

      fileOut->cd();
      graphPoints->Write();
      graphSyst->Write();
    }
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
  // gStyle->SetHistLineWidth(1);
  // gStyle->SetHistLineColor(kRed);
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
  // gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  // gStyle->SetLegendTextSize(0.04);
}
//_____________________________________________________________________________
TH1D* PrepareSyst(TH1D* res, TH1D* syst)
{
  if(!res) {printf("nores\n"); return 0x0;}
  if(!syst) {printf("nosyst\n"); return 0x0;}

  TH1D* hResSyst = (TH1D*) res->Clone(Form("%s_sy",res->GetName()));

  for(Int_t bin(1); bin < hResSyst->GetNbinsX()+1; ++bin)
  {

    Double_t dCont = res->GetBinContent(bin);
    Double_t dSyst = syst->GetBinContent(bin);

    hResSyst->SetBinContent(bin,dCont);
    hResSyst->SetBinError(bin,dCont*dSyst);
  }

  return hResSyst;
}
//_____________________________________________________________________________
