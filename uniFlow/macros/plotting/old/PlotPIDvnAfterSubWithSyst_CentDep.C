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

// plotting order (mainly for syst)
Int_t plotOrder[] = {3,2,1,0}; // indexes in sSpecies_list

Bool_t bDrawSyst = kTRUE;

Int_t iPointsShiftSpec[] = {16,16,16,8,15,0,9};
// regular order
TString sSpecies_labels[] = {" h^{#pm}"," #pi^{#pm}"," K^{#pm}"," K^{0}_{S}"," p(#bar{p})"," #phi"," #Lambda(#bar{#Lambda})"};
TString sSpecies_list[] = {"Charged","Pion","Kaon","K0s","Proton","Phi","Lambda"};
// Color_t colors[] = {kBlack, kRed, kBlue, kCyan+1, kGreen+2, kMagenta, kOrange+1};
Color_t colors[] = {kGreen+2, kRed, kBlue, kOrange-1, kCyan+2, kMagenta+1, kOrange-1};
Color_t colorsFill[] = {kGreen+2, kRed, kBlue, kOrange-1, kCyan+2, kMagenta+1, kOrange-1};
// Color_t colorsFill[] = {kGray, kRed-9, kBlue-9, kCyan-8, kGreen-9, kMagenta-9, kOrange-9};
Int_t markers[] = {kOpenSquare, kOpenCircle, kOpenCross, kOpenDiamond};
Double_t markersSize[] = {1.0,1.0,1.2,1.35};
Double_t dAlpha = 0.25;
// Double_t dAlpha = 1.0;

// // reverse order
// TString sSpecies_list[] = {"Phi","Proton","Lambda","K0s","Charged","Kaon","Pion"};
// Color_t colors[] = {kMagenta,kGreen+2,kOrange+1,kCyan+1,kBlack,kBlue,kRed};
// Int_t markers[] = {kFullCross,kFullSquare,kFullDiamond,kFullTriangleDown,kOpenCircle,kFullTriangleUp,kFullStar};


// TString sGapVal= "0.0"; TString sGapName= "gap00";
TString sGapVal= "0.4"; TString sGapName= "gap04";
// TString sGapVal= "0.8"; TString sGapName= "gap08";

TString sCent[] = {"0-10","10-20","20-40","40-60"};


Int_t iNumCent = sizeof(sCent) / sizeof(sCent[0]);
Int_t iNumSpecies = sizeof(sSpecies_list) / sizeof(sSpecies_list[0]);

void PlotPIDvnAfterSubWithSyst_CentDep()
{
  SetStyle();

  TString sInputPathTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3/"+sGapName+"/";

  TString sOutputPath = Form("%s/results-centdep-4/",sInputPathTop.Data());
  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  TFile* fileIn = OpenFile(Form("%s/v2-syst-final.root",sInputPathTop.Data())); if(!fileIn) {return;}

  for(Int_t spec(0); spec < iNumSpecies; ++spec)
  {
    TString sSpecies = sSpecies_list[spec];

    // TLegend* leg = new TLegend(0.18,0.55,0.4,0.77);
    TLegend* leg = new TLegend(0.75,0.18,0.88,0.44);
    if(spec == 5) { leg = new TLegend(0.18,0.50,0.4,0.70); }
    leg->SetFillColorAlpha(0,0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    // leg->SetHeader(Form("%s |#eta| < 0.8 (V0A)",sSpecies_labels[spec].Data()));
    leg->SetHeader(Form("%s (V0A)",sSpecies_labels[spec].Data()));

    TCanvas* canPlot = new TCanvas("canPlot","canPlot",800,800);
    canPlot->cd();
    // TH1* frame_canPlot = (TH1*) gPad->DrawFrame(0.0,-0.02,7.3,0.30);
    TH1* frame_canPlot = (TH1*) gPad->DrawFrame(0.0,0.0,7.3,0.25);
    frame_canPlot->SetTitle(Form("; #it{p}_{T} (GeV/#it{c}); v_{2}^{sub} {2, |#Delta#eta| > %s}",sGapVal.Data()));

    TList* listSyst = new TList();
    TList* listPoints = new TList();

    for(Int_t iCent(0); iCent < iNumCent; ++iCent)
    {
      TGraphErrors* graphSyst = (TGraphErrors*) fileIn->Get(Form("graphSyst_%s_cent%d",sSpecies.Data(),iCent));
      if(!graphSyst) { printf("No graphSyst \n"); fileIn->ls(); return; }

      TGraphErrors* graphPoints = (TGraphErrors*) fileIn->Get(Form("graphPoints_%s_cent%d",sSpecies.Data(),iCent));
      if(!graphPoints) { printf("No graphPoints\n"); fileIn->ls(); return; }

      graphPoints->SetMarkerSize(markersSize[iCent]);
      graphPoints->SetMarkerStyle(markers[iCent]);
      graphPoints->SetMarkerColor(colors[iCent]);
      graphPoints->SetLineColor(colors[iCent]);

      graphSyst->SetMarkerSize(markersSize[iCent]);
      graphSyst->SetMarkerStyle(markers[iCent]);
      graphSyst->SetMarkerColor(colors[iCent]);
      graphSyst->SetLineColor(colors[iCent]);
      graphSyst->SetFillColorAlpha(colorsFill[iCent],dAlpha);

      // for(Int_t iPoint = 0; iPoint < graphPoints->GetN(); ++iPoint)
      for(Int_t iPoint = iPointsShiftSpec[spec]; iPoint < graphPoints->GetN(); ++iPoint)
      {
        Double_t dX;
        Double_t dY;
        graphPoints->GetPoint(iPoint,dX,dY);

        graphPoints->SetPoint(iPoint,dX+0.1*iCent, dY);
        graphSyst->SetPoint(iPoint,dX+0.1*iCent, dY);

      }

      listSyst->Add(graphSyst);
      listPoints->Add(graphPoints);
      if(iCent == iNumCent-1 && spec == 5) { continue; }
      leg->AddEntry(graphSyst,Form(" %s%%",sCent[iCent].Data()),"pf");
    }

    TLatex * text = new TLatex();
    TLatex * text2 = new TLatex();
    text->SetTextFont(42);
    text->SetTextSize(0.04);
    text2->SetTextFont(42);
    text2->SetTextSize(0.04);

    canPlot->cd();
    text2->DrawLatexNDC(0.186,0.84,Form("ALICE Preliminary"));
    text->DrawLatexNDC(0.186,0.79,Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    text->DrawLatexNDC(0.186,0.74,Form("|#eta| < 0.8"));
    leg->Draw();



    // Drawing

    if(bDrawSyst)
    {
      // systematics first
      for(Int_t i(0); i < iNumCent; ++i)
      {
        // Int_t index = i;
        Int_t index = plotOrder[i];
        TGraphErrors* graphSyst = (TGraphErrors*) listSyst->At(index);
        canPlot->cd();
        if(index == iNumCent-1 && spec == 5)
        {
          continue;
        }
        graphSyst->Draw("p2 same");
      }
    }

    // central points
    for(Int_t i(0); i < iNumCent; ++i)
    {
      Int_t index = i;
      // Int_t index = plotOrder[i];
      TGraphErrors* graphPoints = (TGraphErrors*) listPoints->At(index);
      canPlot->cd();
      if(index == iNumCent-1 && spec == 5)
      {
        continue;
      }
      graphPoints->Draw("pZ same");
    }

    canPlot->SaveAs(Form("%s/PIDvn_species%d.pdf",sOutputPath.Data(),spec),"pdf");
    canPlot->SaveAs(Form("%s/PIDvn_species%d.eps",sOutputPath.Data(),spec),"eps");
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
  gStyle->SetTitleOffset(1.4,"y");
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
