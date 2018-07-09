#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TSystem.h"

TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
TH1D* BarlowTest(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
TH1D* ApplyBarlow(TH1* diff, TH1* barlow, Double_t dCut);
TH1D* Diff(TH1* ratio);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TList* PrepareCanvas(TCanvas* can);

Color_t colBase = kRed;
Color_t colSyst = kBlue+2;

Int_t markBase = kOpenSquare;
Int_t markSyst = kFullCircle;

Double_t dBarlowCut = -1.0;

TString sGap = "gap04";
TString sCentLabels[] = {"0-10%","10-20%","20-40%","40-60%"};
TString sType = "v0s"; TString sTags[] = {"3sigma","decayRad","DCAdaughters","CPA"}; TString sSpeciesList[] = {"K0s","Lambda"};
// TString sType = "tracking"; TString sTags[] = {"FB","PV","cls"}; TString sSpeciesList[] = {"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
// TString sType = "pid"; TString sTags[] = {"3sigma","bayes90"}; TString sSpeciesList[] = {"Pion","Kaon","Proton","Phi"};

void ProcessSyst()
{
  Int_t iTag = 2;
  Int_t iCent = 3;
  Int_t iSpecies = 1;
  Bool_t bCorrelated = 1;

  // X ranges = 0.3
  Double_t dFitXmin = 0.2;
  Double_t dFitXmax = 7.0;

  Double_t dRatioYmin = 0.7;
  Double_t dRatioYmax = 1.3;




  TString sTag = sTags[iTag];
  TString sSpecies = sSpeciesList[iSpecies];
  TString sCentLabel = sCentLabels[iCent];


  TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-binning-3/"+sType+"/";
  TString sInputBase = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/subt/output_binning-3/"+sGap+"/Subtracted.root";
  TString sInputSyst = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/"+sType+"/merged-16qt/subt/output_binning-3/"+sGap+"/"+sTag+"/Subtracted.root";

  //  --------------------------------------------------------
  TFile* fileInputBase = TFile::Open(sInputBase.Data(),"READ");
  if(!fileInputBase) { printf("ERROR : fileInputBase not found!\n"); return; }
  TFile* fileInputSyst = TFile::Open(sInputSyst.Data(),"READ");
  if(!fileInputSyst) { printf("ERROR : fileInputSyst not found!\n"); return; }

  gSystem->mkdir(sOutputPath.Data(),kTRUE);
  TFile* fileOutput = TFile::Open(Form("%s/%s.root",sOutputPath.Data(),sTag.Data()),"RECREATE");
  if(!fileOutput) { printf("ERROR : fileOutput not found!\n"); return; }

  TString sHistName = Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),iCent);

  TH1D* histBase = (TH1D*) fileInputBase->Get(sHistName.Data());
  if(!histBase) { printf("ERROR : Input histBase '%s' does not found\n",sHistName.Data()); fileInputBase->ls(); return; }
  StyleHist(histBase,colBase,markBase);

  TH1D* histSyst = (TH1D*) fileInputSyst->Get(sHistName.Data());
  if(!histSyst) { printf("ERROR : Input histSyst '%s' does not found\n",sHistName.Data()); fileInputSyst->ls(); return; }
  StyleHist(histSyst,colSyst,markSyst);

  TH1D* histRatio = DivideHistos(histSyst,histBase,bCorrelated);
  if(!histRatio) { printf("ERROR: Ratio failed\n"); return; }
  StyleHist(histRatio,colSyst,markSyst);

  TH1D* histBarlow = BarlowTest(histBase,histSyst,bCorrelated);
  if(!histBarlow) { printf("ERROR : Barlow failed\n"); return; }
  StyleHist(histBarlow,colSyst,markSyst);

  TH1D* histDiff = Diff(histRatio);
  if(!histDiff) { printf("ERROR: Diff failed\n"); return; }
  StyleHist(histRatio,colSyst,markSyst);

  // TH1D* histDiff_clone = (TH1D*) histDiff->Clone(Form("%s_clone",histDiff->GetName()));
  TH1D* histDiff_clone = ApplyBarlow(histDiff,histBarlow,dBarlowCut);
  // TF1* fitDiff = new TF1("fitDiff","[0]",histDiff->GetXaxis()->GetXmin(),histDiff->GetXaxis()->GetXmax());
  // histDiff_clone->Fit("fitDiff","MINR");

  TF1* fitRatio = new TF1("fitRatio","[0]",dFitXmin,dFitXmax);
  histRatio->Fit("fitRatio","MNR");


  TList* listOut = new TList();

  histBase->SetNameTitle("histBase",Form("%s (%s): %s; p_{T} (GeV/c); v_{2}{2,%s}",sSpecies.Data(),sCentLabels[iCent].Data(),sTag.Data(),sGap.Data()));
  histSyst->SetNameTitle("histSyst",Form("%s (%s): %s; p_{T} (GeV/c); v_{2}{2,%s}",sSpecies.Data(),sCentLabels[iCent].Data(),sTag.Data(),sGap.Data()));
  histRatio->SetNameTitle("histRatio",Form("Ratio x_{syst} / x_{default} ; p_{T} (GeV/c); syst/default"));
  histDiff->SetNameTitle("histDiff",Form("|1 - (x_{syst} / x_{default})|; p_{T} (GeV/c); Relative diff. (%%)"));
  histBarlow->SetNameTitle("histBarlow",Form("Barlow | x_{syst} - x_{default} | / #sqrt{#sigma^{2}_{syst} #pm #sigma^{2}_{default}}; p_{T} (GeV/c); Barlow"));
  listOut->Add(histBase);
  listOut->Add(histSyst);
  listOut->Add(histRatio);
  listOut->Add(histDiff);
  listOut->Add(histBarlow);
  listOut->Add(fitRatio);
  // listOut->Add(fitDiff);

  fileOutput->cd();
  listOut->Write(Form("%s_cent%d",sSpecies.Data(),iCent),TObject::kSingleKey);

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetNDC();

  TLegend* leg = new TLegend(0.55,0.15,0.88,0.38);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0,0);
  leg->SetHeader(Form("%s (%s)",sSpecies.Data(),sCentLabels[iCent].Data()));
  leg->AddEntry(histBase,"Default","pel");
  leg->AddEntry(histSyst,sTag.Data(),"pel");

  TLine* lineUnity = new TLine();
  lineUnity->SetLineStyle(2);
  lineUnity->SetLineWidth(2);
  lineUnity->SetLineColor(kGray+1);

  TLine* lineSigma = new TLine();
  lineSigma->SetLineStyle(2);
  lineSigma->SetLineWidth(2);
  lineSigma->SetLineColor(kRed);

  TCanvas* can = new TCanvas("can","can",1800,400);
  can->Divide(4,1);
  can->cd(1);
  // TH1* frame_main = (TH1*) gPad->DrawFrame(0.0,0.0,7.0,0.5);
  gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.15);
  histSyst->SetLabelFont(43,"XY");
  histSyst->SetLabelSize(14,"XY");
  histSyst->SetTitleFont(43,"XY");
  histSyst->SetTitleSize(14,"XY");
  histSyst->SetTitleOffset(2.1,"Y");
  histSyst->SetTitleOffset(1.5,"X");

  histSyst->DrawCopy("hist p e1");
  histBase->DrawCopy("same hist p e1");
  leg->Draw();

  can->cd(2);
  gPad->SetBottomMargin(0.13);
  gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.15);
  TH1* frame_ratio = (TH1*) gPad->DrawFrame(histRatio->GetXaxis()->GetXmin(),dRatioYmin,histRatio->GetXaxis()->GetXmax(),dRatioYmax);
  frame_ratio->SetTitle(Form("Ratio x_{syst} / x_{default} ; p_{T} (GeV/c); syst / default"));
  frame_ratio->SetLabelFont(43,"XY");
  frame_ratio->SetLabelSize(14,"XY");
  frame_ratio->SetTitleFont(43,"XY");
  frame_ratio->SetTitleSize(14,"XY");
  frame_ratio->SetTitleOffset(2.1,"Y");
  frame_ratio->SetTitleOffset(1.5,"X");

  histRatio->DrawCopy("same");
  fitRatio->DrawCopy("same");
  lineUnity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.0,histRatio->GetXaxis()->GetXmax(),1.0);
  latex->DrawLatex(0.2,0.85,Form("#color[2]{%s (%s): %s}", sSpecies.Data(), sCentLabel.Data(), sTag.Data()));
  latex->DrawLatex(0.2,0.80,Form("#color[2]{chi2/ndf = %.2f/%d = %.2f}",fitRatio->GetChisquare(), fitRatio->GetNDF(), fitRatio->GetChisquare() / fitRatio->GetNDF() ));
  latex->DrawLatex(0.2,0.75,Form("#color[2]{prob = %.6f }",fitRatio->GetProb()));
  latex->DrawLatex(0.2,0.70,Form("#color[2]{a = %.3f #pm %.3f}",fitRatio->GetParameter(0), fitRatio->GetParError(0)));

  can->cd(3);
  gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.15);
  histDiff->SetLabelFont(43,"XY");
  histDiff->SetLabelSize(14,"XY");
  histDiff->SetTitleFont(43,"XY");
  histDiff->SetTitleSize(14,"XY");
  histDiff->SetTitleOffset(2.1,"Y");
  histDiff->SetTitleOffset(1.5,"X");
  histDiff->DrawCopy();
  // fitDiff->DrawCopy("same");
  // latex->DrawLatex(0.2,0.8,Form("#color[2]{a = %.3f #pm %.3f}",fitDiff->GetParameter(0), fitDiff->GetParError(0)));
  Double_t dPar = TMath::Abs(1.0 - fitRatio->GetParameter(0));
  lineSigma->DrawLine(histDiff->GetXaxis()->GetXmin(),dPar,histDiff->GetXaxis()->GetXmax(),dPar);

  can->cd(4);
  gPad->SetBottomMargin(0.13);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.15);
  histBarlow->SetLabelFont(43,"XY");
  histBarlow->SetLabelSize(14,"XY");
  histBarlow->SetTitleFont(43,"XY");
  histBarlow->SetTitleSize(14,"XY");
  histBarlow->SetTitleOffset(2.1,"Y");
  histBarlow->SetTitleOffset(1.5,"X");
  histBarlow->DrawCopy();
  lineUnity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.0,histRatio->GetXaxis()->GetXmax(),1.0);


  TString sOut = Form("%s/plots_%s",sOutputPath.Data(),sTag.Data());
  gSystem->mkdir(sOut.Data(),kTRUE);
  can->SaveAs(Form("%s/%s_cent%d.pdf",sOut.Data(),sSpecies.Data(),iCent));

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

    if(bCor) { dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3)); }
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBinNom,dContRatio);
    ratio->SetBinError(iBinNom,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
// ==================================================================================================================
TH1D* Diff(TH1* ratio)
{
  if(!ratio) { printf("Diff-ERROR: ratio histo does not exists\n"); return 0x0; }

  TH1D* hDiff = (TH1D*) ratio->Clone(Form("Diff_%s",ratio->GetName()));
  hDiff->Reset();

  for(Int_t bin(1); bin < hDiff->GetNbinsX()+1; ++bin)
  {
    Double_t dCont = ratio->GetBinContent(bin);
    Double_t dErr = ratio->GetBinError(bin);

    hDiff->SetBinContent(bin, TMath::Abs(dCont-1.0));
    hDiff->SetBinError(bin,dErr);
  }

  hDiff->SetMinimum(0.0);

  return hDiff;
}
// ==================================================================================================================
TH1D* BarlowTest(TH1* nom, TH1* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* hBarlow = (TH1D*) nom->Clone(Form("Barlow_%s_%s",nom->GetName(),denom->GetName()));
  hBarlow->Reset();

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


    Double_t dBarlowNom = TMath::Abs(dContDenom-dContNom);
    Double_t dBarlowDenom = 0.0;

    if(bCor) { dBarlowDenom = TMath::Abs(dErrDenom*dErrDenom - dErrNom*dErrNom); }
    else { dBarlowDenom = TMath::Abs(dErrDenom*dErrDenom + dErrNom*dErrNom); }
    // Double_t dBarlowDenom = dErrDenom*dErrDenom - dErrNom*dErrNom;

    if(dBarlowDenom <= 0.0) continue;

    hBarlow->SetBinContent(iBinDenom,dBarlowNom/TMath::Sqrt(dBarlowDenom));
  }

  hBarlow->SetMinimum(0.0);

  return hBarlow;
}
// ==================================================================================================================
TH1D* ApplyBarlow(TH1* diff, TH1* barlow, Double_t dCut)
{
  if(!diff) { printf("ERROR-ApplyBarlow : no diff\n"); return 0x0; }
  if(!barlow) { printf("ERROR-ApplyBarlow : no barlow\n"); return 0x0; }

  TH1D* after = (TH1D*) diff->Clone(Form("%s_Barlowed",diff->GetName()));
  if(dCut < 0.0) { return after; }

  for(Int_t iBin(1); iBin < after->GetNbinsX()+1; ++iBin)
  {
    if(barlow->GetBinContent(iBin) <= dCut)
    {
      after->SetBinError(iBin, 999.9);
    }
  }

  return after;
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
