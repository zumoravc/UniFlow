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

void PlotPIDvn()
{
  TString sInputFile = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/testRefs/UniFlow_PID.root");
  TString sInputFileRecon = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/test/UniFlow_Reconstructed.root");
  // TString sGap = "08";
  Int_t iCent = 4;

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


  // making legend
  TLegend* legend = new TLegend(0.19,0.7,0.3,0.82);
  legend->SetBorderSize(0);
  legend->AddEntry(hFlowCharged,"h^{#pm}","pel");
  legend->AddEntry(hFlowPion,"#pi^{#pm}","pel");
  legend->AddEntry(hFlowKaon,"K^{#pm}","pel");
  legend->AddEntry(hFlowK0s,"K^{0}_{S}","pel");
  legend->AddEntry(hFlowProton,"p/#bar{p}","pel");
  legend->AddEntry(hFlowPhi,"#phi","pel");
  legend->AddEntry(hFlowLambda,"#Lambda/#bar{#Lambda}","pel");


  // drawing stuff
  hFlowCharged->Draw("hist p e1 x0");
  hFlowPion->Draw("hist p e1 x0 same");
  hFlowKaon->Draw("hist p e1 x0  same");
  hFlowProton->Draw("hist p e1 x0 same");
  hFlowPhi->Draw("hist p e1 x0 same");
  hFlowK0s->Draw("hist p e1 x0 same");
  hFlowLambda->Draw("hist p e1 x0 same");
  legend->Draw();

  return;
}
