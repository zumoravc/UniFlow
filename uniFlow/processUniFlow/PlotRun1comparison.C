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
void RatioError(TH1* nominator = 0x0, TH1* nominator_errors = 0x0, TH1* denominator = 0x0, TH1* ratio = 0x0);

void PlotRun1comparison()
{
  LoadLibs();
  SetStyle();

  TString sInputFile = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/run1_comparison/UniFlow.root");
  TString sInputFileRecon = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/merged/UniFlow.root");
  TString sInputFilePublished = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/run1_comparison/HEPdata.root");

  TString sOutputFilePath = TString("/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/run1_comparison/withRatios");

  TString sOutputFormat = TString("png");

  // TString sGap = "08";
  Int_t iCent = 0;
  Int_t iTable = 5; // 5,9,13,17
  TString sCent = TString("0-20");

  // ALICE Preferred colors and markers (from figure template)
  const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
  const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

  // colors setting
  Color_t colCharged = kBlack;
  Color_t colPion = kRed;
  Color_t colKaon = kBlue;
  Color_t colProton = kGreen+2;
  Color_t colPhi = kMagenta+1;
  Color_t colK0s = kCyan+1;
  Color_t colLambda = kOrange+1;

  // markers setting
  Int_t markCharged = kFullSquare;
  Int_t markPion = kFullCircle;
  Int_t markKaon = kFullTriangleUp;
  Int_t markProton = kFullCross;
  // Int_t markProton = kFullSquare;
  Int_t markPhi = kFullStar;
  Int_t markK0s = kFullTriangleDown;
  Int_t markLambda = kFullDiamond;
  // Int_t markLambda = kFullCross;

  Int_t markChargedHEP = kOpenSquare;
  Int_t markPionHEP = kOpenCircle;
  Int_t markKaonHEP = kOpenTriangleUp;
  Int_t markProtonHEP = kOpenCross;

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

  TFile* fInputFilePublished = new TFile(sInputFilePublished.Data(),"READ");
  if(!fInputFilePublished->IsOpen()) return;
  fInputFilePublished->ls();

  // preparing canvas
  //Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
  // Rectangular
  TCanvas* cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600);
  // Square
  //TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800);
  // cfig->SetLogy();
  // Set Titles etc..
  TH1* h = cfig->DrawFrame(0,0,6,0.3);

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

  // published
  fInputFilePublished->cd(Form("Table %d",iTable)); // charged 0-20
  TH1F* hHEP_Charged = (TH1F*) gDirectory->Get("Hist1D_y2");
  if(!hHEP_Charged) { printf("No HEP charged\n"); return; }
  TH1F* hHEP_Charged_err = (TH1F*) gDirectory->Get("Hist1D_y2_e1");
  if(!hHEP_Charged_err) { printf("No HEP charged errors\n"); return; }
  TGraphAsymmErrors* gHEP_Charged = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y2");
  if(!gHEP_Charged) { printf("No HEP charged graph\n"); return; }

  fInputFilePublished->cd(Form("Table %d",iTable+1)); // pions 0-20
  TGraphAsymmErrors* gHEP_Pion = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y2");
  if(!gHEP_Pion) { printf("No HEP Pion graph\n"); return; }
  TH1F* hHEP_Pion = (TH1F*) gDirectory->Get("Hist1D_y2");
  if(!hHEP_Pion) { printf("No HEP Pion\n"); return; }
  TH1F* hHEP_Pion_err = (TH1F*) gDirectory->Get("Hist1D_y2_e1");
  if(!hHEP_Pion_err) { printf("No HEP Pion errors\n"); return; }

  fInputFilePublished->cd(Form("Table %d",iTable+2)); // kaons 0-20
  TGraphAsymmErrors* gHEP_Kaon = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y2");
  if(!gHEP_Kaon) { printf("No HEP Kaon graph\n"); return; }
  TH1F* hHEP_Kaon = (TH1F*) gDirectory->Get("Hist1D_y2");
  if(!hHEP_Kaon) { printf("No HEP Kaon\n"); return; }
  TH1F* hHEP_Kaon_err = (TH1F*) gDirectory->Get("Hist1D_y2_e1");
  if(!hHEP_Kaon_err) { printf("No HEP Kaon errors\n"); return; }

  fInputFilePublished->cd(Form("Table %d",iTable+3)); // protons 0-20
  TGraphAsymmErrors* gHEP_Proton = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y2");
  if(!gHEP_Proton) { printf("No HEP Proton graph\n"); return; }
  TH1F* hHEP_Proton = (TH1F*) gDirectory->Get("Hist1D_y2");
  if(!hHEP_Proton) { printf("No HEP Proton\n"); return; }
  TH1F* hHEP_Proton_err = (TH1F*) gDirectory->Get("Hist1D_y2_e1");
  if(!hHEP_Proton_err) { printf("No HEP Proton errors\n"); return; }



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

  gHEP_Charged->SetLineColor(colCharged);
  gHEP_Charged->SetFillColor(colCharged);
  gHEP_Charged->SetFillStyle(3001);
  gHEP_Charged->SetMarkerColor(colCharged);
  // gHEP_Charged->SetMarkerStyle(markChargedHEP);

  gHEP_Pion->SetLineColor(colPion);
  gHEP_Pion->SetFillColor(colPion);
  gHEP_Pion->SetFillStyle(3001);
  gHEP_Pion->SetMarkerColor(colPion);
  // gHEP_Pion->SetMarkerStyle(markPionHEP);

  gHEP_Kaon->SetLineColor(colKaon);
  gHEP_Kaon->SetFillColor(colKaon);
  gHEP_Kaon->SetFillStyle(3001);
  gHEP_Kaon->SetMarkerColor(colKaon);
  // gHEP_Kaon->SetMarkerStyle(markKaonHEP);

  gHEP_Proton->SetLineColor(colProton);
  gHEP_Proton->SetFillColor(colProton);
  gHEP_Proton->SetFillStyle(3001);
  gHEP_Proton->SetMarkerColor(colProton);
  // gHEP_Proton->SetMarkerStyle(markProtonHEP);

  // drawing stuff
  // drawing HEP
  gHEP_Pion->Draw("same p2");
  gHEP_Kaon->Draw("same p2");
  gHEP_Proton->Draw("same p2");
  gHEP_Charged->Draw("same p2");
  gHEP_Charged->Draw("same p");
  gHEP_Pion->Draw("same p");
  gHEP_Kaon->Draw("same p");
  gHEP_Proton->Draw("same p");
  // hHEP_Charged->Draw("hist p same");
  // hHEP_Charged_err->Draw("hist e1 same");


  // hFlowCharged->Draw("same");
  hFlowCharged->Draw("hist p e1 x0 same");
  hFlowPion->Draw("hist p e1 x0 same");
  hFlowKaon->Draw("hist p e1 x0  same");
  hFlowProton->Draw("hist p e1 x0 same");
  // hFlowK0s->Draw("hist p e1 x0 same");
  // hFlowLambda->Draw("hist p e1 x0 same");
  // hFlowCharged->Draw("hist p e1 x0 same");
  // hFlowPhi->Draw("hist p e1 x0 same");


  // Draw the logo
  //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
  //  >0: ALICE Preliminary
  DrawLogo(2, 0.61, 0.83);

  // You should always specify the colliding system
  // NOTATION: pp, p-Pb, Pb-Pb.
  // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
  // You can change the position of this with
  TLatex * text = new TLatex (0.3,0.27,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  text->Draw();
  TLatex * text2 = new TLatex (0.3,0.25,Form("Multiplicity Class %s%% (V0A)",sCent.Data()));
  text2->SetTextSizePixels(22);
  text2->Draw();
  // TLatex * text3 = new TLatex (0.55,0.76,"|#eta| < 0.8");
  // text3->SetTextSizePixels(20);
  // text3->Draw();


  // making legend
  TLegend * legend = new TLegend(0.62, 0.18, 0.8, 0.46);
  legend->SetFillColorAlpha(0,0);
  legend->SetTextSize(gStyle->GetTextSize()*0.8);
  // legend->SetBorderSize(1);

  TLegend * legend2 = new TLegend(0.72, 0.18, 0.91, 0.46);
  legend2->SetFillColorAlpha(0,0);
  legend2->SetTextSize(gStyle->GetTextSize()*0.8);
  // legend2->SetBorderSize(1);

  legend->SetHeader("Run2");
  legend->AddEntry(hFlowCharged,"h^{#pm}","pl");
  legend->AddEntry(hFlowPion,"#pi^{#pm}","pl");
  legend->AddEntry(hFlowKaon,"K^{#pm}","pl");
  legend->AddEntry(hFlowProton,"p/#bar{p}","pl");
  // legend->AddEntry(hFlowK0s,"K^{0}_{S}","pl");
  // legend2->AddEntry(hFlowPhi,"#phi","pl");
  // legend2->AddEntry(hFlowLambda,"#Lambda/#bar{#Lambda}","pl");

  legend2->SetHeader("Run1 (SP)");
  legend2->AddEntry(gHEP_Charged," ","fl");
  legend2->AddEntry(gHEP_Pion," ","fl");
  legend2->AddEntry(gHEP_Kaon," ","fl");
  legend2->AddEntry(gHEP_Proton," ","fl");

  legend->Draw();
  legend2->Draw();

  // making ratios
  TCanvas* canRatio = new TCanvas("canRatio","canRatio",600,800);
  TH1* hRatio = canRatio->DrawFrame(0,0.5,4,1.5);
  hRatio->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  hRatio->SetYTitle("Run1 / Run2");

  TH1F* hRatio_Charged = hHEP_Charged->Clone("hRatio_Charged");
  hRatio_Charged->Divide(hFlowCharged);
  RatioError(hHEP_Charged,hHEP_Charged_err,hFlowCharged,hRatio_Charged);
  hRatio_Charged->SetLineColor(colCharged);
  hRatio_Charged->SetMarkerColor(colCharged);
  hRatio_Charged->SetMarkerStyle(markCharged);

  TH1F* hRatio_Pion = hHEP_Pion->Clone("hRatio_Pion");
  hRatio_Pion->Divide(hFlowPion);
  RatioError(hHEP_Pion,hHEP_Pion_err,hFlowPion,hRatio_Pion);
  hRatio_Pion->SetLineColor(colPion);
  hRatio_Pion->SetMarkerColor(colPion);
  hRatio_Pion->SetMarkerStyle(markPion);

  TH1F* hRatio_Kaon = hHEP_Kaon->Clone("hRatio_Kaon");
  printf("Bins: kaon HEP: %d | my %d\n",hHEP_Kaon->GetNbinsX(),hFlowKaon->GetNbinsX());
  hRatio_Kaon->Divide(hFlowKaon);
  RatioError(hHEP_Kaon,hHEP_Kaon_err,hFlowKaon,hRatio_Kaon);
  hRatio_Kaon->SetLineColor(colKaon);
  hRatio_Kaon->SetMarkerColor(colKaon);
  hRatio_Kaon->SetMarkerStyle(markKaon);

  TH1F* hRatio_Proton = hHEP_Proton->Clone("hRatio_Proton");
  printf("Bins: proton HEP: %d | my %d\n",hHEP_Proton->GetNbinsX(),hFlowProton->GetNbinsX());
  hRatio_Proton->Divide(hFlowProton);
  RatioError(hHEP_Proton,hHEP_Proton_err,hFlowProton,hRatio_Proton);
  hRatio_Proton->SetLineColor(colProton);
  hRatio_Proton->SetMarkerColor(colProton);
  hRatio_Proton->SetMarkerStyle(markProton);


  hRatio_Charged->Draw("same hist p e");
  hRatio_Pion->Draw("same hist p e");
  hRatio_Kaon->Draw("same hist p e");
  hRatio_Proton->Draw("same hist p e");


  // cfig->SaveAs(Form("%s/PID_%d.%s",sOutputFilePath.Data(),iCent,sOutputFormat.Data()));
  cfig->SaveAs(Form("%s/PID_%d.pdf",sOutputFilePath.Data(),iCent));
  cfig->SaveAs(Form("%s/PID_%d.png",sOutputFilePath.Data(),iCent));
  cfig->SaveAs(Form("%s/PID_%d.eps",sOutputFilePath.Data(),iCent));

  canRatio->SaveAs(Form("%s/PID_ratio_%d.pdf",sOutputFilePath.Data(),iCent));
  canRatio->SaveAs(Form("%s/PID_ratio_%d.png",sOutputFilePath.Data(),iCent));
  canRatio->SaveAs(Form("%s/PID_ratio_%d.eps",sOutputFilePath.Data(),iCent));

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
void RatioError(TH1* nominator, TH1* nominator_error, TH1* denominator, TH1* ratio)
{
  if(!ratio || ! denominator || !nominator || !nominator_error) { printf("Either of input histos not found!\n"); return; }

  Int_t iNbins = ratio->GetNbinsX();
  // printf("points: %d | bins %d\n", iNpoints,iNbins);
  // return;
  if(nominator->GetNbinsX() != iNbins || denominator->GetNbinsX() != iNbins || nominator_error->GetNbinsX() != iNbins) { printf("Different binning!\n"); return; }


  Double_t content_deno = 0, error_deno = 0;
  Double_t content_no = 0;
  Double_t error_no = 0;
  Double_t final_error = 0;

  for(Int_t bin(1); bin < iNbins+1; bin++)
  {
    content_deno = denominator->GetBinContent(bin);
    error_deno = denominator->GetBinError(bin);
    //
    content_no = nominator->GetBinContent(bin);
    error_no = nominator_error->GetBinContent(bin);
    //
    final_error = TMath::Power(error_no/content_deno,2) + TMath::Power(error_deno*content_no/(content_deno*content_deno),2);
    ratio->SetBinError(bin,TMath::Sqrt(final_error));
  }
  printf("done\n");
  return;
}
