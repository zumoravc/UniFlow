#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TF1.h"

Bool_t bCorrelated = kTRUE;

Double_t dRatioToler = 0.2;
// Double_t dRatioYmin = 0.7;
// Double_t dRatioYmax = 1.3;

Color_t colBase = kBlue;
Color_t colSyst = kRed+1;

Int_t markBase = kFullSquare;
Int_t markSyst = kFullCircle;

// PlotTogether
// PlotRatio
// Plot Relative Difference
// Plot | relative difference|
// Barlow


TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
TH1D* BarlowTest(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
Double_t ApplyBarlow(TH1* diff, TH1* barlow);
TH1D* Diff(TH1* ratio);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TCanvas* PrepareCanvas(TString name);
TList* PrepareCanvasRatioSub(TCanvas* can);

Bool_t ProcessSingle(
    TString hist,
    TString path,
    TString syst,
    TString baseline = "default"
);


void Syst()
{

    std::vector<TString> vecSyst;

    // vecSyst.push_back("FB768");
    // vecSyst.push_back("PID3sigma");
    vecSyst.push_back("PVz8");
    vecSyst.push_back("TPCcls90");
    // vecSyst.push_back("V0sCPA099");
    vecSyst.push_back("V0sCrossFind1");
    // vecSyst.push_back("V0sDaugDCA3");
    // vecSyst.push_back("V0sDaugPt02");
    vecSyst.push_back("V0sDecRad10");
    vecSyst.push_back("V0sFinderOn");
    vecSyst.push_back("V0sPVDCA3");

    Int_t iNumCent = 7;
    TString species = "Lambda";

    TString path = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/%s/",species.Data());


    TString histoName[] = {
        "<<3>>(4,-2,-2)_2sub(0)",
        "<<3>>(5,-3,-2)_2sub(0)",
        "<<3>>(6,-3,-3)_2sub(0)"
    };

    Int_t iNumHist = sizeof(histoName) / sizeof(histoName[0]);

    for(Int_t iCent(0); iCent < iNumCent; ++iCent) {
        // Int_t iCent = 2;
        for(Int_t iHist(0); iHist < iNumHist; ++iHist) {
            // Int_t iHist = 0;
            for(Int_t iSyst(0); iSyst < (Int_t) vecSyst.size(); ++iSyst) {

                if(!ProcessSingle(Form("%s_%s_mult%d",species.Data(),histoName[iHist].Data(),iCent),path,vecSyst.at(iSyst))) { return; }

            }
        }
    }


    // ProcessSingle(histoName.Data(), path, syst);
}

Bool_t ProcessSingle(TString hist, TString path, TString syst, TString baseline)
{
    TFile* fileBase = TFile::Open(Form("%s/%s/Processed.root",path.Data(),baseline.Data()),"READ");
    if(!fileBase) { printf("E: Baseline file '%s/%s/Processed.root' not found!",path.Data(),baseline.Data()); return kFALSE; }

    TFile* fileSyst = TFile::Open(Form("%s/%s/Processed.root",path.Data(),syst.Data()),"READ");
    if(!fileSyst) { printf("E: Systematic file '%s/%s/Processed.root' not found!",path.Data(),syst.Data()); return kFALSE; }

    TH1D* histBase = (TH1D*) fileBase->Get(hist.Data());
    if(!histBase) { printf("E: Baseline histo '%s' not found!",hist.Data()); fileBase->ls(); return kFALSE; }
    StyleHist(histBase,colBase,markBase);
    histBase->SetName("Base");

    TH1D* histSyst = (TH1D*) fileSyst->Get(hist.Data());
    if(!histSyst) { printf("E: Systematic histo '%s' not found!",hist.Data()); fileSyst->ls(); return kFALSE; }
    StyleHist(histSyst,colSyst,markSyst);
    histBase->SetName("Syst");

    TH1D* histRatio = DivideHistos(histSyst,histBase,bCorrelated);
    if(!histRatio) { printf("ERROR: Ratio failed\n"); return kFALSE; }

    TH1D* histBarlow = BarlowTest(histBase,histSyst,bCorrelated);
    if(!histBarlow) { printf("ERROR : Barlow failed\n"); return kFALSE; }

    TH1D* histDiff = Diff(histRatio);
    if(!histDiff) { printf("ERROR: Diff failed\n"); return kFALSE; }

    Double_t xmin = histBase->GetXaxis()->GetXmin();
    Double_t xmax = histBase->GetXaxis()->GetXmax();

    // fit Diff
    TF1* fitDiff = new TF1("fitDiff","[0]",xmin,xmax);
    fitDiff->SetLineColor(kBlue+2);
    fitDiff->SetLineStyle(7);
    histDiff->Fit(fitDiff,"RN");

    // TH1D* histAfterBarlow = ApplyBarlow(histDiff, histBarlow);
    // TF1* fitAfterBarlow = new TF1("fitAfterBarlow","[0]",xmin,xmax);
    // fitAfterBarlow->SetLineColor(kGreen+2);
    // fitAfterBarlow->SetLineStyle(9);
    // histAfterBarlow->Fit(fitAfterBarlow,"RN");

    Double_t dAfterBarlow = ApplyBarlow(histDiff, histBarlow);
    TF1* fitAfterBarlow = new TF1("fitAfterBarlow","[0]",xmin,xmax);
    fitAfterBarlow->FixParameter(0,dAfterBarlow);
    fitAfterBarlow->SetLineStyle(9);
    fitAfterBarlow->SetLineColor(kGreen+2);

    // Plotting stuff


    // TLegend* leg = new TLegend(0.55,0.15,0.88,0.38);
    TLegend* leg = new TLegend(0.16,0.7,0.5,0.88);
    leg->SetTextFont(43);
    leg->SetTextSize(14);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0);
    leg->AddEntry(histBase,"Baseline","pel");
    leg->AddEntry(histSyst,syst,"pel");

    TLine* lineUnity = new TLine();
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(2);
    lineUnity->SetLineColor(kGray+1);

    TText* text = new TText();
    text->SetNDC();
    text->SetTextFont(43);
    text->SetTextSize(17);

    // TLatex* latex = new TLatex();
    // latex->SetTextSize(0.05);
    // latex->SetNDC();

    TCanvas* can = PrepareCanvas(Form("can_%s",syst.Data()));
    can->cd(1);
    histBase->DrawCopy("e1");
    histSyst->DrawCopy("same e1");
    leg->Draw();

    can->cd(2);
    histBarlow->DrawCopy("B");
    lineUnity->DrawLine(xmin,1.0,xmax,1.0);

    can->cd(3);
    histRatio->DrawCopy("e1");
    lineUnity->DrawLine(xmin,1.0,xmax,1.0);

    can->cd(4);
    histDiff->DrawCopy("e1");
    fitDiff->DrawCopy("same");
    fitAfterBarlow->DrawCopy("same");
    lineUnity->DrawLine(xmin,0.0,xmax,0.0);
    text->SetTextColor(kBlue+2);
    text->DrawText(0.34,0.78,Form("Pol0 = %.3f +- %.3f", fitDiff->GetParameter(0),fitDiff->GetParError(0)));
    text->SetTextColor(kGreen+2);
    text->DrawText(0.34,0.84,Form("Pol0 = %.3f +- %.3f", fitAfterBarlow->GetParameter(0),fitAfterBarlow->GetParError(0)));


    // saving output
    gSystem->mkdir(Form("%s/%s/plots_root/",path.Data(),syst.Data()),kTRUE);
    gSystem->mkdir(Form("%s/%s/plots_pdf/",path.Data(),syst.Data()),kTRUE);
    gSystem->mkdir(Form("%s/%s/plots_eps/",path.Data(),syst.Data()),kTRUE);

    TFile* fileOut = TFile::Open(Form("%s/%s/plots_root/%s.root",path.Data(),syst.Data(),hist.Data()),"RECREATE");
    if(!fileOut) { printf("E: Output file not created!"); return kFALSE; }

    TList* outList = new TList();
    outList->Add(histBase);
    outList->Add(histSyst);
    outList->Add(histRatio);
    outList->Add(histBarlow);
    // outList->Add(histAfterBarlow);
    outList->Add(histDiff);
    outList->Add(fitDiff);
    outList->Add(fitAfterBarlow);

    fileOut->cd();
    outList->Write("list",TObject::kSingleKey);

    can->SaveAs(Form("%s/%s/plots_pdf/%s.pdf",path.Data(),syst.Data(),hist.Data()),"pdf");
    can->SaveAs(Form("%s/%s/plots_eps/%s.eps",path.Data(),syst.Data(),hist.Data()),"eps");

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
  hist->SetLabelFont(43,"XY");
  hist->SetLabelSize(15,"XY");
  hist->SetTitleFont(43,"XY");
  hist->SetTitleSize(15,"XY");
  hist->SetTitleOffset(2.1,"Y");
  hist->SetTitleOffset(2.0,"X");
  return;
};
// ==================================================================================================================
TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  // TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  TH1D* ratio = (TH1D*) nom->Clone("Ratio");
  ratio->Reset();
  ratio->SetName("ratio");
  ratio->SetTitle("Ratio syst / default");
  ratio->SetMaximum(1.0+dRatioToler);
  ratio->SetMinimum(1.0-dRatioToler);
  StyleHist(ratio, colSyst, markSyst);

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

  // TH1D* hDiff = (TH1D*) ratio->Clone(Form("Diff_%s",ratio->GetName()));
  TH1D* hDiff = (TH1D*) ratio->Clone("Diff");
  hDiff->Reset();
  hDiff->SetTitle("Abs. deviation: |Ratio - 1|");
  hDiff->SetMinimum(-0.01);
  hDiff->SetMaximum(dRatioToler);
  StyleHist(hDiff, colSyst, markSyst);

  for(Int_t bin(1); bin < hDiff->GetNbinsX()+1; ++bin)
  {
    Double_t dCont = ratio->GetBinContent(bin);
    Double_t dErr = ratio->GetBinError(bin);

    hDiff->SetBinContent(bin, TMath::Abs(dCont-1.0));
    hDiff->SetBinError(bin,dErr);
  }

  return hDiff;
}
// ==================================================================================================================
TH1D* BarlowTest(TH1* nom, TH1* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  // TH1D* hBarlow = (TH1D*) nom->Clone(Form("Barlow_%s_%s",nom->GetName(),denom->GetName()));
  TH1D* hBarlow = (TH1D*) nom->Clone("Barlow");
  hBarlow->Reset();
  if(bCor) {
      hBarlow->SetTitle("Barlow: |x-y| / |#sigma^{2}_{x} - #sigma^{2}_{y}|");
  } else {
      hBarlow->SetTitle("Barlow: |x-y| / |#sigma^{2}_{x} + #sigma^{2}_{y}|");
  }
  StyleHist(hBarlow, Color_t(colSyst+1), markSyst);
  hBarlow->SetFillColorAlpha(colSyst,1.0);

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
  hBarlow->SetMaximum(10.0);

  return hBarlow;
}
// ==================================================================================================================
Double_t ApplyBarlow(TH1* diff, TH1* barlow)
{
  if(!diff) { printf("ERROR-ApplyBarlow : no diff\n"); return -1.0; }
  if(!barlow) { printf("ERROR-ApplyBarlow : no barlow\n"); return -1.0; }

  // TH1D* after = (TH1D*) diff->Clone(Form("%s_Barlowed",diff->GetName()));
  TH1D* after = (TH1D*) diff->Clone("AfterBarlow");
  // if(dCut < 0.0) { return after; }

  Double_t dSum = 0.0;
  Double_t dSumWeights = 0.0;

  for(Int_t iBin(1); iBin < after->GetNbinsX()+1; ++iBin) {
    Double_t dBarlow = barlow->GetBinContent(iBin);
    Double_t dWeight = 0.0;
    if(dBarlow > 0.0) { dWeight = 1.0 / dBarlow; }

    Double_t dContent = diff->GetBinContent(iBin);

    dSum += dWeight*dContent;
    dSumWeights += dWeight;

    // after->SetBinContent(iBin, diff->GetBinContent(iBin));
    // after->SetBinError(iBin, dWeight);
  }

  Double_t dFit = dSum / dSumWeights;

  return dFit;
}
// ==================================================================================================================
TCanvas* PrepareCanvas(TString name)
{
    TCanvas* can = new TCanvas(name.Data(), name.Data(), 700,700);
    can->Divide(2,2);

    can->cd(1);
    // TH1* frame_main = (TH1*) gPad->DrawFrame(0.0,0.0,7.0,0.5);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.15);
    // histSyst->SetLabelFont(43,"XY");
    // histSyst->SetLabelSize(14,"XY");
    // histSyst->SetTitleFont(43,"XY");
    // histSyst->SetTitleSize(14,"XY");
    // histSyst->SetTitleOffset(2.1,"Y");
    // histSyst->SetTitleOffset(1.5,"X");
    //
    // histSyst->DrawCopy("hist p e1");
    // histBase->DrawCopy("same hist p e1");
    // leg->DrawCopy();

    can->cd(2);
    gPad->SetBottomMargin(0.13);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.15);

    // TH1* frame_ratio = (TH1*) gPad->DrawFrame(histRatio->GetXaxis()->GetXmin(),dRatioYmin,histRatio->GetXaxis()->GetXmax(),dRatioYmax);
    // frame_ratio->SetTitle(Form("Ratio x_{syst} / x_{default} ; p_{T} (GeV/c); syst / default"));
    // frame_ratio->SetLabelFont(43,"XY");
    // frame_ratio->SetLabelSize(14,"XY");
    // frame_ratio->SetTitleFont(43,"XY");
    // frame_ratio->SetTitleSize(14,"XY");
    // frame_ratio->SetTitleOffset(2.1,"Y");
    // frame_ratio->SetTitleOffset(1.5,"X");

    // histRatio->DrawCopy("same");
    // fitRatio->DrawCopy("same");
    // lineUnity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.0,histRatio->GetXaxis()->GetXmax(),1.0);
    // latex->DrawLatex(0.2,0.85,Form("#color[2]{%s (%s): %s}", sSpecies.Data(), sCentLabel.Data(), sTag.Data()));
    // latex->DrawLatex(0.2,0.80,Form("#color[2]{chi2/ndf = %.2f/%d = %.2f}",fitRatio->GetChisquare(), fitRatio->GetNDF(), fitRatio->GetChisquare() / fitRatio->GetNDF() ));
    // latex->DrawLatex(0.2,0.75,Form("#color[2]{prob = %.6f }",fitRatio->GetProb()));
    // latex->DrawLatex(0.2,0.70,Form("#color[2]{a = %.3f #pm %.3f}",fitRatio->GetParameter(0), fitRatio->GetParError(0)));

    can->cd(3);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.15);
    // histDiff->SetLabelFont(43,"XY");
    // histDiff->SetLabelSize(14,"XY");
    // histDiff->SetTitleFont(43,"XY");
    // histDiff->SetTitleSize(14,"XY");
    // histDiff->SetTitleOffset(2.1,"Y");
    // histDiff->SetTitleOffset(1.5,"X");
    // histDiff->DrawCopy();
    // // fitDiff->DrawCopy("same");
    // // latex->DrawLatex(0.2,0.8,Form("#color[2]{a = %.3f #pm %.3f}",fitDiff->GetParameter(0), fitDiff->GetParError(0)));
    // Double_t dPar = TMath::Abs(1.0 - fitRatio->GetParameter(0));
    // lineSigma->DrawLine(histDiff->GetXaxis()->GetXmin(),dPar,histDiff->GetXaxis()->GetXmax(),dPar);

    can->cd(4);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.03);
    gPad->SetLeftMargin(0.15);
    // histBarlow->SetLabelFont(43,"XY");
    // histBarlow->SetLabelSize(14,"XY");
    // histBarlow->SetTitleFont(43,"XY");
    // histBarlow->SetTitleSize(14,"XY");
    // histBarlow->SetTitleOffset(2.1,"Y");
    // histBarlow->SetTitleOffset(1.5,"X");
    // histBarlow->DrawCopy();
    // lineUnity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.0,histRatio->GetXaxis()->GetXmax(),1.0);

    return can;
}
// ==================================================================================================================
TList* PrepareCanvasRatioSub(TCanvas* can)
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
