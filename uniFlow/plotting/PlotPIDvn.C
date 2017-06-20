/* ProcessUniFlow class
 *
 * Class implemented for plotting identified flow based on ProcessUniFlow class results
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();

void PlotPIDvn()
{
  LoadLibs();
  SetStyle();

  TString sInputFile = TString("/Users/vpacik/NBI/Flow/results/uniFlow_syst/NUA_cor/merged_16q/UniFlowTest.root");
  TString sInputFileRecon = TString("/Users/vpacik/NBI/Flow/results/uniFlow_syst/NUA_cor/merged_16q/UniFlowTest.root");
  TString sOutputFilePath = TString("/Users/vpacik/NBI/Flow/results/uniFlow_syst/NUA_cor/merged_16q/Run1comp");

  // TString sGap = "08";
  // Int_t iCent = 0;
  const Short_t iNumCent = 4;
  TString sCent[iNumCent] = {"0-20","20-40","40-60","60-100"};
  Double_t dYmin[iNumCent] = {0,0,0,0};
  Double_t dYmax[iNumCent] = {0.3,0.3,0.5,0.5};

  // ALICE Preferred colors and markers (from figure template)
  const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

  // colors setting
  Color_t colCharged = kBlack;
  Color_t colPion = kRed;
  Color_t colKaon = kBlue;
  Color_t colProton = kGreen+2;
  Color_t colPhi = kMagenta;
  Color_t colK0s = kCyan+1;
  Color_t colLambda = kOrange+1;

  // markers setting
  Int_t markCharged = kOpenSquare;
  Int_t markPion = kFullCircle;
  Int_t markKaon = kFullTriangleUp;
  Int_t markProton = kFullCross;
  // Int_t markProton = kFullSquare;
  Int_t markPhi = kFullStar;
  Int_t markK0s = kFullTriangleDown;
  Int_t markLambda = kFullDiamond;
  // Int_t markLambda = kFullCross;

  Double_t markSizeCharged = 1;
  Double_t markSizePion = 1;
  Double_t markSizeKaon = 1.2;
  Double_t markSizeProton = 1.5;
  Double_t markSizePhi = 1.5;
  Double_t markSizeK0s = 1.2;
  Double_t markSizeLambda = 1.6;
  // Double_t markSizeCharged = 2;
  // Double_t markSizePion = 2;
  // Double_t markSizeKaon = 2;
  // Double_t markSizeProton = 2;
  // Double_t markSizePhi = 2;
  // Double_t markSizeK0s = 2;
  // Double_t markSizeLambda = 2;

  // openning input file
  TFile* fInputFile = new TFile(sInputFile.Data(),"READ");
  if(!fInputFile->IsOpen()) return;
  fInputFile->ls();

  TFile* fInputFileRecon = new TFile(sInputFileRecon.Data(),"READ");
  if(!fInputFileRecon->IsOpen()) return;
  fInputFileRecon->ls();

  // preparing canvas
  //Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas* cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600);
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800);
  // cfig->SetLogy();
  // Set Titles etc..


  TLatex * text = new TLatex();
  TLatex * text2 = new TLatex();
  text2->SetTextSizePixels(22);


  // loop over centrality
  for(Short_t iCent(0); iCent < iNumCent; iCent++)
  {
    TH1* h = cfig->DrawFrame(0,dYmin[iCent],6,dYmax[iCent]);
    // Set titles
    h->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    h->SetYTitle("#it{v}_{2} {2,|#eta| > 0.8}");


    // loading histos

    // DIRECT
    TH1D* hFlowCharged = (TH1D*) fInputFile->Get(Form("hFlow2_Charged_harm2_gap08_cent%d",iCent));
    if(!hFlowCharged) { printf("No charged\n"); return; }
    TH1D* hFlowPion = (TH1D*) fInputFile->Get(Form("hFlow2_Pion_harm2_gap08_cent%d",iCent));
    if(!hFlowPion) { printf("No pion\n"); return; }
    TH1D* hFlowKaon = (TH1D*) fInputFile->Get(Form("hFlow2_Kaon_harm2_gap08_cent%d",iCent));
    if(!hFlowKaon) { printf("No kaon\n"); return; }
    TH1D* hFlowProton = (TH1D*) fInputFile->Get(Form("hFlow2_Proton_harm2_gap08_cent%d",iCent));
    if(!hFlowProton)  { printf("No proton\n"); return; }

    // reconstructed
    TH1D* hFlowK0s = (TH1D*) fInputFileRecon->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent));
    if(!hFlowK0s) { printf("No K0s\n"); return; }
    TH1D* hFlowLambda = (TH1D*) fInputFileRecon->Get(Form("hFlow2_Lambda_harm2_gap08_mult%d",iCent));
    if(!hFlowLambda) { printf("No Lambda\n"); return; }
    TH1D* hFlowPhi = (TH1D*) fInputFileRecon->Get(Form("hFlow2_Phi_harm2_gap08_mult%d",iCent));
    if(!hFlowPhi) { printf("No Phi\n"); return; }

    // setting histos
    // hFlowCharged->SetStats(kFALSE);
    // hFlowCharged->GetXaxis()->SetRangeUser(0.2,10.);
    // hFlowCharged->SetMinimum(0.);
    // hFlowCharged->SetMaximum(0.5);

    hFlowCharged->SetLineColor(colCharged);
    hFlowCharged->SetMarkerColor(colCharged);
    hFlowCharged->SetMarkerStyle(markCharged);
    hFlowCharged->SetMarkerSize(markSizeCharged);

    hFlowPion->SetLineColor(colPion);
    hFlowPion->SetMarkerColor(colPion);
    hFlowPion->SetMarkerStyle(markPion);
    hFlowPion->SetMarkerSize(markSizePion);

    hFlowKaon->SetLineColor(colKaon);
    hFlowKaon->SetMarkerColor(colKaon);
    hFlowKaon->SetMarkerStyle(markKaon);
    hFlowKaon->SetMarkerSize(markSizeKaon);

    hFlowProton->SetLineColor(colProton);
    hFlowProton->SetMarkerColor(colProton);
    hFlowProton->SetMarkerStyle(markProton);
    hFlowProton->SetMarkerSize(markSizeProton);

    hFlowPhi->SetLineColor(colPhi);
    hFlowPhi->SetMarkerColor(colPhi);
    hFlowPhi->SetMarkerStyle(markPhi);
    hFlowPhi->SetMarkerSize(markSizePhi);

    hFlowK0s->SetLineColor(colK0s);
    hFlowK0s->SetMarkerColor(colK0s);
    hFlowK0s->SetMarkerStyle(markK0s);
    hFlowK0s->SetMarkerSize(markSizeK0s);

    hFlowLambda->SetLineColor(colLambda);
    hFlowLambda->SetMarkerColor(colLambda);
    hFlowLambda->SetMarkerStyle(markLambda);
    hFlowLambda->SetMarkerSize(markSizeLambda);



    // drawing stuff
    // hFlowCharged->Draw("same");
    hFlowPion->Draw("hist p e1 x0 same");
    hFlowKaon->Draw("hist p e1 x0  same");
    hFlowProton->Draw("hist p e1 x0 same");
    hFlowK0s->Draw("hist p e1 x0 same");
    hFlowLambda->Draw("hist p e1 x0 same");
    hFlowCharged->Draw("hist p e1 x0 same");
    hFlowPhi->Draw("hist p e1 x0 same");

    // Draw the logo
    //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
    //  >0: ALICE Preliminary
    // DrawLogo(2, 0.186, 0.83);

    // You should always specify the colliding system
    // NOTATION: pp, p-Pb, Pb-Pb.
    // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
    // You can change the position of this with


    text->DrawLatexNDC(0.18,0.83,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    text2->DrawLatexNDC(0.18,0.78,Form("Centrality Class %s%% (V0A)",sCent[iCent].Data()));
    // TLatex * text = new TLatex(0.3,0.25,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    // text->Draw();
    // TLatex * text2 = new TLatex (0.3,0.23,Form("Multiplicity Class %s%% (V0A)",sCent.Data()));
    // text2->SetTextSizePixels(22);
    // text2->Draw();

    // TLatex * text3 = new TLatex (0.55,0.76,"|#eta| < 0.8");
    // text3->SetTextSizePixels(20);
    // text3->Draw();


    // making legend
    // TLegend * legend = new TLegend(0.62, 0.18, 0.8, 0.46);
    TLegend * legend = new TLegend(0.19, 0.48, 0.37, 0.76);
    legend->SetFillColorAlpha(0,0);
    legend->SetTextSize(gStyle->GetTextSize()*0.8);
    // legend->SetBorderSize(1);

    // TLegend * legend2 = new TLegend(0.72, 0.18, 0.91, 0.39);
    TLegend * legend2 = new TLegend(0.29, 0.55, 0.48, 0.76);
    legend2->SetFillColorAlpha(0,0);
    legend2->SetTextSize(gStyle->GetTextSize()*0.8);
    // legend2->SetBorderSize(1);

    legend->AddEntry(hFlowCharged,"h^{#pm}","pl");
    legend->AddEntry(hFlowPion,"#pi^{#pm}","pl");
    legend->AddEntry(hFlowKaon,"K^{#pm}","pl");
    legend->AddEntry(hFlowK0s,"K^{0}_{S}","pl");
    legend2->AddEntry(hFlowProton,"p/#bar{p}","pl");
    legend2->AddEntry(hFlowPhi,"#phi","pl");
    legend2->AddEntry(hFlowLambda,"#Lambda/#bar{#Lambda}","pl");

    legend->Draw();
    legend2->Draw();

    // cfig->SaveAs(Form("%s/plots/PID_%d.%s",sOutputFilePath.Data(),iCent,sOutputFormat.Data()));
    cfig->SaveAs(Form("%s/PID_%d.pdf",sOutputFilePath.Data(),iCent));
    cfig->SaveAs(Form("%s/PID_%d.png",sOutputFilePath.Data(),iCent));
    cfig->SaveAs(Form("%s/PID_%d.eps",sOutputFilePath.Data(),iCent));

  }


  return;
}
//_____________________________________________________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}
//_____________________________________________________________________________
void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"
  TString sText;

  switch (logo)
  {
    case 0:
      sText = TString("ALICE");
      break;

    case 1:
      sText = TString("ALICE Preliminary");
      break;

    case 2:
      sText = TString("ALICE This thesis");
      break;
  }


  TLatex *   tex = new TLatex(xmin,ymin, sText.Data());
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();

  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}
//_____________________________________________________________________________
void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}
//_____________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;

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
void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2();
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2();
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

}
//_____________________________________________________________________________
