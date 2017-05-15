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

  // colors setting
  Color_t colCharged = kBlack;
  Color_t colPion = kRed;
  Color_t colKaon = kBlue;
  Color_t colProton = kGreen+2;

  // markers setting
  Int_t markCharged = 20;
  Int_t markPion = 20;
  Int_t markKaon = 20;
  Int_t markProton = 20;

  // openning input file
  TFile* fInputFile = new TFile(sInputFile.Data(),"READ");
  if(!fInputFile->IsOpen()) return;
  fInputFile->ls();

  // loading histos
  TH1D* hFlowCharged = (TH1D*) fInputFile->Get("hFlow_Charged_harm2_gap08_cent0_taskcharged");
  if(!hFlowCharged) return;
  TH1D* hFlowPion = (TH1D*) fInputFile->Get("hFlow_Pion_harm2_gap08_cent0_taskPi");
  if(!hFlowPion) return;
  TH1D* hFlowKaon = (TH1D*) fInputFile->Get("hFlow_Kaon_harm2_gap08_cent0_taskK");
  if(!hFlowKaon) return;
  TH1D* hFlowProton = (TH1D*) fInputFile->Get("hFlow_Proton_harm2_gap08_cent0_taskp");
  if(!hFlowProton) return;

  // setting histos
  hFlowCharged->SetStats(kFALSE);
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

  // making legend
  TLegend* legend = new TLegend(0.19,0.7,0.5,0.82);
  legend->SetBorderSize(0);
  legend->AddEntry(hFlowCharged,"h^{#pm}","pel");
  legend->AddEntry(hFlowPion,"#pi","pel");
  legend->AddEntry(hFlowKaon,"K","pel");
  legend->AddEntry(hFlowProton,"p","pel");


  // drawing stuff
  hFlowCharged->Draw("hist p e1 x0");
  hFlowPion->Draw("hist p e1 x0 same");
  hFlowKaon->Draw("hist p e1 x0  same");
  hFlowProton->Draw("hist p e1 x0 same");
  legend->Draw();

  return;
}
