/* Comp_13TeV_5TeV
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

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

TList* PrepareCanvas(TCanvas* can);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);

TString sLegendHeader = "h^{#pm}  v_{2}^{sub} (0-10% V0A)"; TString sHistName = "hFlow2_Charged_harm2_gap04_cent0";
// TString sLegendHeader = "h^{#pm} v_{2}^{sub} (10-20% V0A)"; TString sHistName = "hFlow2_Charged_harm2_gap04_cent1";
// TString sLegendHeader = "h^{#pm} v_{2}^{sub} (20-40% V0A)"; TString sHistName = "hFlow2_Charged_harm2_gap04_cent2";
// TString sLegendHeader = "h^{#pm} v_{2}^{sub} (40-60% V0A)"; TString sHistName = "hFlow2_Charged_harm2_gap04_cent3";

TString sInput13TeV = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3/gap04/Subtracted.root";
TString sInput5TeV = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3-5TeVpp/gap04/Subtracted.root";
TString sOutput = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/Comp_13TeV_5TeV/";

void Comp_subt_13TeV_5TeV()
{
  TFile* file13TeV = TFile::Open(sInput13TeV.Data(),"READ");
  if(!file13TeV) { printf("no 13 TEV file\n"); return; }

  TFile* file5TeV = TFile::Open(sInput5TeV.Data(),"READ");
  if(!file5TeV) { printf("no 5 TEV file\n"); return; }

  TH1D* hist13TeV = (TH1D*) file13TeV->Get(sHistName.Data());
  if(!hist13TeV) { printf("no hist13TeV\n"); return; }

  TH1D* hist5TeV = (TH1D*) file5TeV->Get(sHistName.Data());
  if(!hist5TeV) { printf("no hist5TeV\n"); return; }


  StyleHist(hist5TeV,kGreen+2,kFullSquare);
  StyleHist(hist13TeV,kBlue+2,kFullCircle);
  // StyleHist(histRatio,kBlue+2,kFullCircle);

  TH1D* histRatio = DivideHistos(hist5TeV,hist13TeV,0);
  if(!histRatio) { printf("no histRatio\n"); return; }


  TCanvas* canCent = new TCanvas("canCent","canCent",500,1200);
  canCent->cd();
  TH1* frameMainCent = (TH1*) gPad->DrawFrame(0.0,0.0,7.0,0.4,Form("; p_{T} (GeV/c); v_{2}{2,|#Delta#eta|>0.4}"));

  // TCanvas* canCent = new TCanvas("canCent","canCent",500,1200);
  // TList* listCent = PrepareCanvas(canCent);
  // if(!listCent) { printf("ERROR : PrepareCanvas failed\n"); return; }
  // TPad* padMain = (TPad*) listCent->At(0);
  // TPad* padRatio = (TPad*) listCent->At(1);
  // TH1* frameMainCent = (TH1*) listCent->At(2);
  // TH1* frameRatioCent = (TH1*) listCent->At(3);



  frameMainCent->SetTitle(Form("; p_{T} (GeV/c); v^{sub}_{2}{2,|#Delta#eta|>0.4}"));
  // frameRatioCent->SetTitle(Form("; p_{T} (GeV/c); 5 TeV / 13 TeV   "));

  TLegend* legCent = new TLegend(0.2,0.7,0.5,0.89);
  legCent->SetFillColorAlpha(0,0);
  legCent->SetBorderSize(0);
  legCent->SetHeader(sLegendHeader.Data());
  legCent->AddEntry(hist13TeV,"pp 13 TeV (16kl)","p");
  legCent->AddEntry(hist5TeV,"pp 5 TeV (17p)","p");

  padMain->cd();
  hist5TeV->DrawCopy("same");
  hist13TeV->DrawCopy("same");
  legCent->Draw();

  padRatio->cd();
  histRatio->DrawCopy("same");

  gSystem->mkdir(sOutput.Data(),kTRUE);
  canCent->SaveAs(Form("%s/v2_sub_%s.pdf",sOutput.Data(),sHistName.Data()),"pdf");

  return;
}
// ==================================================================================================================
TList* PrepareCanvas(TCanvas* can)
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
  TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(0.0,0.0,7.0,0.4,Form("; p_{T} (GeV/c); v_{2}{2,|#Delta#eta|>0.4}"));
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
  TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(0.0,0.7,7.0,1.3);
  frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); Kch / K0s   "));
  frame_canDiff_2->SetNdivisions(505,"Y");
  frame_canDiff_2->SetTitleFont(43,"XY");
  frame_canDiff_2->SetTitleSize(18,"XY");
  frame_canDiff_2->SetTitleOffset(4.3,"X");
  frame_canDiff_2->SetTitleOffset(2.2,"Y");
  frame_canDiff_2->SetLabelFont(43,"XY");
  frame_canDiff_2->SetLabelSize(18,"XY");
  lUnity->DrawLine(0.0,1.0,10.0,1.0);

  TList* list = new TList();
  list->Add(padMain);
  list->Add(padRatio);
  list->Add(frame_canDiff_1);
  list->Add(frame_canDiff_2);

  return list;
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
TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor)
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
  for(Short_t iBinNom(1); iBinNom < binsNom+1; iBinNom++)
  {
    Double_t dBinCenter = nom->GetBinCenter(iBinNom);
    Int_t iBinDenom = denom->FindFixBin(dBinCenter);

    dContNom = nom->GetBinContent(iBinNom);
    dErrNom = nom->GetBinError(iBinNom);
    dContDenom = denom->GetBinContent(iBinDenom);
    dErrDenom = denom->GetBinError(iBinDenom);

    if(dContDenom == 0.0) continue;

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBinNom,dContRatio);
    ratio->SetBinError(iBinNom,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
// ==================================================================================================================
