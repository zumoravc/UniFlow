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

  TString sInputFile = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/testRefs/UniFlow_PID.root");
  TString sInputFileRecon = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test/UniFlow_Reconstructed.root");
  // TString sGap = "08";
  Int_t iCent = 5;

  // ALICE Preferred colors and markers (from figure template)
  const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

  // colors setting
  Color_t colCharged = kBlack;
  Color_t colPion = kRed+1;
  Color_t colKaon = kBlue+1;
  Color_t colProton = kGreen+3;
  Color_t colPhi = kMagenta+1;
  Color_t colK0s = kCyan+2;
  Color_t colLambda = kYellow+2;

  // markers setting
  Int_t markCharged = kOpenCircle;
  Int_t markPion = kFullCircle;
  Int_t markKaon = kFullTriangleDown;
  Int_t markProton = kFullSquare;
  Int_t markPhi = kFullStar;
  Int_t markK0s = kFullTriangleUp;
  Int_t markLambda = kFullDiamond;

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
  TH1* h = cfig->DrawFrame(0,0,10,1);

  // Set titles
  h->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  h->SetYTitle("#it{v}_{2}{2,|#eta| > 0.8}");

  // loading histos

  // DIRECT
  TH1D* hFlowCharged = (TH1D*) fInputFile->Get(Form("hFlow_Charged_harm2_gap08_cent%d_taskcharged",iCent));
  if(!hFlowCharged) return;
  TH1D* hFlowPion = (TH1D*) fInputFile->Get(Form("hFlow_Pion_harm2_gap08_cent%d_taskPi",iCent));
  if(!hFlowPion) return;
  TH1D* hFlowKaon = (TH1D*) fInputFile->Get(Form("hFlow_Kaon_harm2_gap08_cent%d_taskK",iCent));
  if(!hFlowKaon) return;
  TH1D* hFlowProton = (TH1D*) fInputFile->Get(Form("hFlow_Proton_harm2_gap08_cent%d_taskp",iCent));
  if(!hFlowProton) return;

  // reconstructed
  TH1D* hFlowK0s = (TH1D*) fInputFileRecon->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent));
  if(!hFlowK0s) return;
  TH1D* hFlowLambda = (TH1D*) fInputFileRecon->Get(Form("hFlow2_Lambda_harm2_gap08_mult%d",iCent));
  if(!hFlowLambda) return;
  TH1D* hFlowPhi = (TH1D*) fInputFileRecon->Get(Form("hFlow2_Phi_harm2_gap08_mult%d",iCent));
  if(!hFlowPhi) return;

  // setting histos
  hFlowCharged->SetStats(kFALSE);
  hFlowCharged->GetXaxis()->SetRangeUser(0.2,10.);
  hFlowCharged->SetMinimum(0.);
  hFlowCharged->SetMaximum(0.5);

  hFlowCharged->SetLineColor(colCharged);
  hFlowCharged->SetMarkerColor(colCharged);
  hFlowCharged->SetMarkerStyle(markCharged);

  hFlowPion->SetLineColor(colPion);
  hFlowPion->SetMarkerColor(colPion);
  hFlowPion->SetMarkerStyle(markPion);

  hFlowKaon->SetLineColor(colKaon);
  hFlowKaon->SetMarkerColor(colKaon);
  hFlowKaon->SetMarkerStyle(markKaon);

  hFlowProton->SetLineColor(colProton);
  hFlowProton->SetMarkerColor(colProton);
  hFlowProton->SetMarkerStyle(markProton);

  hFlowPhi->SetLineColor(colPhi);
  hFlowPhi->SetMarkerColor(colPhi);
  hFlowPhi->SetMarkerStyle(markPhi);

  hFlowK0s->SetLineColor(colK0s);
  hFlowK0s->SetMarkerColor(colK0s);
  hFlowK0s->SetMarkerStyle(markK0s);

  hFlowLambda->SetLineColor(colLambda);
  hFlowLambda->SetMarkerColor(colLambda);
  hFlowLambda->SetMarkerStyle(markLambda);



  // drawing stuff
  // hFlowCharged->Draw("same");
  hFlowCharged->Draw("hist p e1 x0 same");
  hFlowPion->Draw("hist p e1 x0 same");
  hFlowKaon->Draw("hist p e1 x0  same");
  hFlowProton->Draw("hist p e1 x0 same");
  hFlowPhi->Draw("hist p e1 x0 same");
  hFlowK0s->Draw("hist p e1 x0 same");
  hFlowLambda->Draw("hist p e1 x0 same");

  // Draw the logo
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  DrawLogo(2, 0.59, 0.81);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb.
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  TLatex * text = new TLatex (0.55,0.88,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  text->Draw();
  TLatex * text2 = new TLatex (0.55,0.81,"Multiplicity Class 0-10% (V0A)");
  text2->SetTextSizePixels(20);
  text2->Draw();
  TLatex * text3 = new TLatex (0.55,0.76,"|#eta| < 0.8");
  text3->SetTextSizePixels(20);
  text3->Draw();


  // making legend
  TLegend * legend = new TLegend(0.22, 0.40, 0.32, 0.7);
  legend->SetFillColor(0);
  legend->SetTextSize(gStyle->GetTextSize()*0.8);
  // legend->SetBorderSize(1);

  TLegend * legend2 = new TLegend(0.35, 0.49, 0.45, 0.7);
  legend2->SetFillColor(0);
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
