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

  hFlowPion->SetLineColor(colPion);
  hFlowPion->SetMarkerColor(colPion);

  hFlowKaon->SetLineColor(colKaon);
  hFlowKaon->SetMarkerColor(colKaon);

  hFlowProton->SetLineColor(colProton);
  hFlowProton->SetMarkerColor(colProton);

  // drawing stuff
  hFlowCharged->Draw();
  hFlowPion->Draw("same");
  hFlowKaon->Draw("same");
  hFlowProton->Draw("same");

  return;
}
