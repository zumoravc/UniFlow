/*
*  CompareKchK0s.C
*  ==============================
*  Macro for comparing flow results of Kch and K0s.
*
*  author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI
*/

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"

Bool_t SingleCase(TString sGap, TString sCent, TString sTag,TString sInputFile, TString sOutputDir);

void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
void SetCustomPalette();
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);
TList* PrepareCanvas(TCanvas* can,TString sGap);

void CompareKchK0s()
{
  TString sGaps[] = { "gap00","gap04","gap08"};
  // TString sGaps[] = {"gap04","gap08"};
  // TString sGaps[] = {"gap08"};
  Int_t iNumGaps = sizeof(sGaps) / sizeof(sGaps[0]);

  TString sCentLabel[] = {"0-20%","20-40%","40-60%"};
  // TString sCentLabel[] = {"0-10%","10-20%","20-40%","40-60%"};
  Int_t iNumCents = sizeof(sCentLabel) / sizeof(sCentLabel[0]);

  TString sPathTop = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/";
  TString sInputPathRel = "merged-pPb-16qt-nua/subt/output_binning-3-run1/";
  TString sInputFieName = "Subtracted.root";
  // ======================
  // looping over centrality and gaps
  for(Int_t cent(0); cent < iNumCents; ++cent)
  {
    for(Int_t gap(0); gap < iNumGaps; ++gap)
    {
      TString sGap = sGaps[gap];
      TString sCent = Form("cent%d",cent);
      TString sTag = "v_{2}^{sub} ("+sCentLabel[cent]+")";


      TString sInputFile = sPathTop+sInputPathRel+sGap+"/"+sInputFieName;
      TString sOutputDir = sPathTop+"outcome/comp-KchK0s/"+sInputPathRel+sGap;

      Bool_t bProcess = SingleCase(sGap,sCent,sTag,sInputFile,sOutputDir);
      if(!bProcess) { printf("ERROR: Single case (cent %d, gap %d) failed\n", cent, gap); return; }
    }
  }

  return;
}
// ===================================================================================================================
Bool_t SingleCase(TString sGap, TString sCent, TString sTag, TString sInputFile, TString sOutputDir)
{
  // ===================================================================================================================
  TFile* fileInput = TFile::Open(sInputFile.Data(),"READ");
  if(!fileInput) { printf("ERROR : file not found!\n"); return kFALSE; }

  TH1D* hFlowK0s = (TH1D*) fileInput->Get(Form("hFlow2_K0s_harm2_%s_%s",sGap.Data(),sCent.Data()));
  if(!hFlowK0s) { printf("ERROR : hFlowK0s not found!\n"); return kFALSE; }
  TH1D* hFlowKch = (TH1D*) fileInput->Get(Form("hFlow2_Kaon_harm2_%s_%s",sGap.Data(),sCent.Data()));
  if(!hFlowKch) { printf("ERROR : hFlowKch not found!\n"); return kFALSE; }
  TH1D* hRatio = DivideHistos(hFlowKch,hFlowK0s);
  if(!hRatio) { printf("ERROR : hRatio not found!\n"); return kFALSE; }

  StyleHist(hFlowKch, kRed, kOpenCircle);
  StyleHist(hFlowK0s, kGreen+2, kOpenSquare);
  StyleHist(hRatio, kBlue+2, kFullSquare);

  Double_t dXmin = hFlowKch->GetXaxis()->GetXmin();
  Double_t dXmax = hFlowKch->GetXaxis()->GetXmax();

  TLegend* leg = new TLegend(0.2,0.7,0.4,0.89);
  leg->SetFillColorAlpha(0,0);
  leg->SetBorderSize(0);
  leg->AddEntry(hFlowKch,"Kch","pl");
  leg->AddEntry(hFlowK0s,"K0s","pl");

  TLine* lUnity = new TLine();
  lUnity->SetLineColor(kGray+1);
  lUnity->SetLineStyle(kDashed);

  TCanvas* can_KchK0s = new TCanvas(Form("canDiff"),Form("canDiff"),500,1200);
  can_KchK0s->cd();
  TPad* padMain = new TPad("padMain","padMain", 0, 0.3, 1, 1.0);
  padMain->SetBottomMargin(0.0);
  padMain->SetRightMargin(0.03);
  padMain->SetLeftMargin(0.13);
  padMain->Draw();
  padMain->cd();
  TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(dXmin,-0.05,dXmax,0.5,Form("%s; p_{T} (GeV/c); v_{2}{2,%s}",sTag.Data(),sGap.Data()));
  frame_canDiff_1->SetTitleFont(43,"X");
  frame_canDiff_1->SetTitleSize(18,"X");
  frame_canDiff_1->SetTitleOffset(4.3,"X");
  frame_canDiff_1->SetLabelFont(43,"X");
  frame_canDiff_1->SetLabelSize(18,"X");
  // frame_canDiff_1->SetTitleFont(43,"Y");
  // frame_canDiff_1->SetTitleOffset(2.2,"Y");

  can_KchK0s->cd();
  TPad* padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.3);
  padRatio->SetTopMargin(0.0);
  padRatio->SetBottomMargin(0.25);
  padRatio->SetRightMargin(0.03);
  padRatio->SetLeftMargin(0.13);
  padRatio->Draw();
  padRatio->cd();
  TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(dXmin,0.7,dXmax,1.3);
  frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); Kch / K0s   "));
  frame_canDiff_2->SetNdivisions(505,"Y");
  frame_canDiff_2->SetTitleFont(43,"XY");
  frame_canDiff_2->SetTitleSize(18,"XY");
  frame_canDiff_2->SetTitleOffset(4.3,"X");
  frame_canDiff_2->SetTitleOffset(2.2,"Y");
  frame_canDiff_2->SetLabelFont(43,"XY");
  frame_canDiff_2->SetLabelSize(18,"XY");
  lUnity->DrawLine(dXmin,1.0,dXmax,1.0);

  padMain->cd();
  hFlowKch->DrawCopy("same");
  hFlowK0s->DrawCopy("same");
  leg->Draw();

  padRatio->cd();
  hRatio->DrawCopy("same");

  gSystem->mkdir(sOutputDir.Data(), kTRUE);
  can_KchK0s->SaveAs(Form("%s/KchK0s_%s.pdf",sOutputDir.Data(),sCent.Data()),"pdf");

  return kTRUE;
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
TList* PrepareCanvas(TCanvas* can,TString sGap)
{
  if(!can) { printf("ERROR-PrepareCanvas: no canvas found!\n"); return 0x0; }

  TLine* lUnity = new TLine();
  lUnity->SetLineColor(kGray+1);
  lUnity->SetLineStyle(kDashed);

  can->cd();
  TPad* padMain = new TPad("padMain","padMain", 0, 0.3, 1, 1.0);
  padMain->SetBottomMargin(0.0);
  padMain->SetRightMargin(0.03);
  padMain->SetLeftMargin(0.13);
  padMain->Draw();
  padMain->cd();
  TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(0.0,-0.05,7.0,1.0,Form("%s; p_{T} (GeV/c); v_{2}{2,%s}","",sGap.Data()));
  frame_canDiff_1->SetTitleFont(43,"XY");
  frame_canDiff_1->SetTitleSize(18,"XY");
  frame_canDiff_1->SetTitleOffset(4.3,"X");
  frame_canDiff_1->SetLabelFont(43,"X");
  frame_canDiff_1->SetLabelSize(18,"X");
  // frame_canDiff_1->SetTitleFont(43,"Y");
  frame_canDiff_1->SetTitleOffset(2.2,"Y");

  can->cd();
  TPad* padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.3);
  padRatio->SetTopMargin(0.0);
  padRatio->SetBottomMargin(0.25);
  padRatio->SetRightMargin(0.03);
  padRatio->SetLeftMargin(0.13);
  padRatio->Draw();
  padRatio->cd();
  TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(0.0,0.4,7.0,1.6);
  frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); Kch / K0s   "));
  frame_canDiff_2->SetNdivisions(505,"Y");
  frame_canDiff_2->SetTitleFont(43,"XY");
  frame_canDiff_2->SetTitleSize(18,"XY");
  frame_canDiff_2->SetTitleOffset(4.3,"X");
  frame_canDiff_2->SetTitleOffset(2.2,"Y");
  frame_canDiff_2->SetLabelFont(43,"XY");
  frame_canDiff_2->SetLabelSize(18,"XY");
  lUnity->DrawLine(0.0,1.0,4.0,1.0);

  TList* list = new TList();
  list->Add(padMain);
  list->Add(padRatio);
  list->Add(frame_canDiff_1);
  list->Add(frame_canDiff_2);


  return list;
}
