#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLine.h"

TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TList* PrepareCanvas(TCanvas* can);
TH1D* MakeVn(TH1D* hDn, TH1D* hCn, Int_t iCent);



TString sCentLabels[] = {"0-10%","10-20%", "20-40%", "40-60%", "60-100%"};

Color_t colThis = kGreen+2;
Color_t colAlex = kBlue;

Int_t markAlex = kOpenCircle;
Int_t markThis = kOpenSquare;

Int_t markAlex_neg = kFullCircle;
Int_t markThis_neg = kFullSquare;

void PlotCompAlex()
{
  TString sFileAlex = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/comparisonAlex/histUQMVojtech/hist_uQ_gap04_Q03_V0Acent_pp16qSDD.root";
  // TString sPathThis = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pPb-run3-gaps-04-06-10-12/output_vn_compAlex/gap04/SP_nonscaled_noneventweighted/";
  // TString sPathOut = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/comparisonAlex/pPb-raw/";
  // Int_t iNumCent = 5;

  TString sPathOut = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/comparisonAlex/pp-5TeV-raw/";
  TString sPathThis = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/multiplicity-fluctuations/pp-NchRFP-gap04-17p/output_vn_compAlex/gap04/SP_nonscaled_noneventweighted/";
  // TString sPathThis = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pp-run3-gaps-04-06-10-12/output_vn_compAlex/gap04/SP_nonscaled_noneventweighted/";
  Int_t iNumCent = 1;


  gSystem->mkdir(sPathOut.Data(),1);

  TFile* fileAlex = TFile::Open(sFileAlex.Data(), "READ"); if(!fileAlex) { printf("No fileAlex\n"); return; }
  TFile* fileThisPos = TFile::Open(Form("%s/Processed_pos.root",sPathThis.Data()), "READ"); if(!fileThisPos) { printf("No fileThisPos\n"); return; }
  TFile* fileThisNeg = TFile::Open(Form("%s/Processed_neg.root",sPathThis.Data()), "READ"); if(!fileThisNeg) { printf("No fileThisNeg\n"); return; }

  fileAlex->ls();
  fileThisPos->ls();

  // reference
  TH1D* histAlex_cn = (TH1D*) fileAlex->Get("hQaQbpp"); if(!histAlex_cn) { printf("no histAlex_cn\n"); return; }
  // TH1D* histAlex_cn = (TH1D*) fileAlex->Get("hQaQbpPb"); if(!histAlex_cn) { printf("no histAlex_cn\n"); return; }
  StyleHist(histAlex_cn,colAlex,markAlex);
  TH1D* histThis_cn = (TH1D*) fileThisPos->Get("hCum2_Refs_harm2_gap04"); if(!histThis_cn) { printf("no histThis_cn\n"); return; }
  StyleHist(histThis_cn,colThis,markThis);

  TH1D* histRatio_cn = DivideHistos(histAlex_cn, histThis_cn,0);

  TLegend* leg_cn = new TLegend(0.7,0.6,0.88,0.88);
  leg_cn->SetFillColorAlpha(0,0);
  leg_cn->SetBorderSize(0);
  // leg_cn->SetHeader(Form("Unsubt. p-Pb"));
  leg_cn->SetHeader(Form("pp (0-100%)"));
  leg_cn->AddEntry(histThis_cn,"Method C","p");
  leg_cn->AddEntry(histAlex_cn,"Alex","p");

  TCanvas* can_cn = new TCanvas("can_cn","can_cn",500,700);
  can_cn->Divide(1,2);
  can_cn->cd(1);
  TH1* frame_cn_1 = (TH1*) gPad->DrawFrame(0.0,0.0,100.0,2.0);
  frame_cn_1->SetTitle("; p_{T} (GeV/c); Q*Q");
  histAlex_cn->DrawCopy("same");
  histThis_cn->DrawCopy("same");
  leg_cn->Draw();
  can_cn->cd(2);
  TH1* frame_cn_2 = (TH1*) gPad->DrawFrame(0.0,0.7,100.0,1.3);
  frame_cn_2->SetTitle("; p_{T} (GeV/c); ratio Alex/This");
  histRatio_cn->DrawCopy("same");
  can_cn->SaveAs(Form("%s/QQ.pdf",sPathOut.Data()),"pdf");

  // Int_t iCentInx = 0;
  for(Int_t iCentInx = 0; iCentInx < iNumCent; ++iCentInx)
  {
    // differenttial
    // TH1D* histAlex_dn_pos = (TH1D*) fileAlex->Get(Form("huQBpPb_%d",iCentInx)); if(!histAlex_dn_pos) { printf("no histAlex_dn_pos\n"); return; }
    TH1D* histAlex_dn_pos = (TH1D*) fileAlex->Get(Form("huQBpp")); if(!histAlex_dn_pos) { printf("no histAlex_dn_pos\n"); return; }
    StyleHist(histAlex_dn_pos,colAlex,markAlex);
    TH1D* histThis_dn_pos = (TH1D*) fileThisPos->Get(Form("hCum2_Charged_harm2_gap04_cent%d",iCentInx)); if(!histThis_dn_pos) { printf("no histThis_dn_pos\n"); return; }
    StyleHist(histThis_dn_pos,colThis,markThis);

    // TH1D* histAlex_dn_neg = (TH1D*) fileAlex->Get(Form("huQApPb_%d",iCentInx)); if(!histAlex_dn_neg) { printf("no histAlex_dn_neg\n"); return; }
    TH1D* histAlex_dn_neg = (TH1D*) fileAlex->Get(Form("huQApp")); if(!histAlex_dn_neg) { printf("no histAlex_dn_neg\n"); return; }
    StyleHist(histAlex_dn_neg,colAlex,markAlex_neg);
    TH1D* histThis_dn_neg = (TH1D*) fileThisNeg->Get(Form("hCum2_Charged_harm2_gap04_cent%d",iCentInx)); if(!histThis_dn_neg) { printf("no histThis_dn_neg\n"); return; }
    StyleHist(histThis_dn_neg,colThis,markThis_neg);

    TH1D* histRatio_dn_pos = DivideHistos(histAlex_dn_pos,histThis_dn_pos,0);
    TH1D* histRatio_dn_neg = DivideHistos(histAlex_dn_neg,histThis_dn_neg,0);

    TLegend* leg_dn = new TLegend(0.14,0.6,0.3,0.88);
    // leg_dn->SetHeader(Form("Unsubt. p-Pb (%s)",sCentLabels[iCentInx].Data()));
    leg_dn->SetHeader(Form("pp"));
    leg_dn->SetFillColorAlpha(0,0);
    leg_dn->SetBorderSize(0);
    leg_dn->AddEntry(histThis_dn_pos,"Method C (pos eta)","p");
    leg_dn->AddEntry(histThis_dn_neg,"Method C (neg eta)","p");
    leg_dn->AddEntry(histAlex_dn_pos,"Alex (Qb)","p");
    leg_dn->AddEntry(histAlex_dn_neg,"Alex (Qa)","p");


    TCanvas* can_diff = new TCanvas("can_diff","can_diff",500,700);
    can_diff->Divide(1,2);
    can_diff->cd(1);
    TH1* frame_diff_1 = (TH1*) gPad->DrawFrame(0.0,0.0,10.0,0.5);
    frame_diff_1->SetTitle("; p_{T} (GeV/c); u*Q");
    histThis_dn_pos->Draw("same");
    histThis_dn_neg->Draw("same");
    histAlex_dn_pos->Draw("same");
    histAlex_dn_neg->Draw("same");
    leg_dn->Draw();
    can_diff->cd(2);
    TH1* frame_diff_2 = (TH1*) gPad->DrawFrame(0.0,0.7,10.0,1.3);
    frame_diff_2->SetTitle("; p_{T} (GeV/c); ratio Alex / Method C");
    histRatio_dn_pos->Draw("same");
    histRatio_dn_neg->Draw("same");
    can_diff->SaveAs(Form("%s/uQ_cent%d.pdf",sPathOut.Data(),iCentInx),"pdf");

    // vns
    TH1D* histAlex_vn_pos = MakeVn(histAlex_dn_pos,histAlex_cn, iCentInx+1);
    TH1D* histAlex_vn_neg = MakeVn(histAlex_dn_neg,histAlex_cn, iCentInx+1);
    TH1D* histThis_vn_pos = MakeVn(histThis_dn_pos,histThis_cn, iCentInx+1);
    TH1D* histThis_vn_neg = MakeVn(histThis_dn_neg,histThis_cn, iCentInx+1);

    TH1D* histRatio_vn_pos = DivideHistos(histAlex_vn_pos,histThis_vn_pos,0);
    TH1D* histRatio_vn_neg = DivideHistos(histAlex_vn_neg,histThis_vn_neg,0);

    TCanvas* can_vn = new TCanvas("can_vn","can_vn",500,700);
    can_vn->Divide(1,2);
    can_vn->cd(1);
    TH1* frame_vn_1 = (TH1*) gPad->DrawFrame(0.0,0.0,10.0,1.2);
    frame_vn_1->SetTitle("; p_{T} (GeV/c); v_{2}{2,|#Delta#eta|> 0.4}");
    histThis_vn_pos->Draw("same");
    histThis_vn_neg->Draw("same");
    histAlex_vn_pos->Draw("same");
    histAlex_vn_neg->Draw("same");
    leg_dn->Draw();
    can_vn->cd(2);
    TH1* frame_vn_2 = (TH1*) gPad->DrawFrame(0.0,0.5,10.0,1.5);
    frame_vn_2->SetTitle("; p_{T} (GeV/c); ratio Alex / Method C");
    histRatio_vn_pos->Draw("same");
    histRatio_vn_neg->Draw("same");
    can_vn->SaveAs(Form("%s/vn_cent%d.pdf",sPathOut.Data(),iCentInx),"pdf");

  }

  return;
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
  TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(0.0,-0.05,4.0,0.4,Form("; p_{T} (GeV/c); v_{2}{2,|#Delta#eta|>0.8}"));
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
  TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(0.0,0.65,4.0,1.35);
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
// ==================================================================================================================
TH1D* MakeVn(TH1D* hDn, TH1D* hCn, Int_t iCent)
{
  if(!hDn) { printf("ERROR-MakeVn: Hist 'hDn' does not found.\n"); return 0x0; }
  if(!hCn) { printf("ERROR-MakeVn: Hist 'hCn' does not found.\n"); return 0x0; }

  TH1D* hVn = (TH1D*) hDn->Clone(Form("%s_vn",hDn->GetName()));
  if(!hVn) { printf("ERROR-MakeVn: Hist 'hVn' does not cloned properly.\n"); return 0x0; }
  hVn->Reset();

  Double_t dConCn = hCn->GetBinContent(iCent);
  Double_t dErrCn = hCn->GetBinError(iCent);
  if(dConCn <= 0.0) { printf("ERROR-MakeVn: Division by zero (cent %d : %f)\n",iCent, dConCn); hCn->Draw(); return hVn; }

  for(Int_t bin(0); bin < hDn->GetNbinsX()+1; ++bin)
  {
    Double_t dConDn = hDn->GetBinContent(bin);
    Double_t dErrDn = hDn->GetBinError(bin);

    Double_t dConVn = dConDn / TMath::Sqrt(dConCn);
    Double_t dErrVnSq = dErrDn*dErrDn / dConCn + 0.25*dConDn*dConDn*dErrCn*dErrCn/(dConCn*dConCn*dConCn);

    hVn->SetBinContent(bin, dConVn);
    hVn->SetBinError(bin, TMath::Sqrt(dErrVnSq));
  }

  return hVn;
}
// ==================================================================================================================
