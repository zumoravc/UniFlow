/* ProcessUniFlow class
 *
 * Class implemented for plotting identified flow based on ProcessUniFlow class results
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void PlotPIDvn()
{
  TString sInputFile = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/testRefs/UniFlow_PID.root");
  TFile* fInputFile = new TFile(sInputFile.Data(),"READ");
  if(!fInputFile->IsOpen()) return;

  fInputFile->ls();

  TH1D* hFlowCharged = (TH1D*) fInputFile->Get("hFlow_Charged_harm2_gap08_cent0_taskcharged");
  if(!hFlowCharged) return;
  TH1D* hFlowPion = (TH1D*) fInputFile->Get("hFlow_Pion_harm2_gap08_cent0_taskPi");
  if(!hFlowPion) return;
  TH1D* hFlowKaon = (TH1D*) fInputFile->Get("hFlow_Kaon_harm2_gap08_cent0_taskK");
  if(!hFlowKaon) return;
  TH1D* hFlowProton = (TH1D*) fInputFile->Get("hFlow_Proton_harm2_gap08_cent0_taskp");
  if(!hFlowProton) return;

  Color_t colCharged = kBlack;
  Color_t colPion = kRed;
  Color_t colKaon = kBlue;
  Color_t colProton = kGreen+2;

  hFlowCharged->SetStats(kFALSE);
  hFlowCharged->SetLineColor(colCharged);
  hFlowCharged->SetMarkerColor(colCharged);

  hFlowPion->SetLineColor(colPion);
  hFlowPion->SetMarkerColor(colPion);

  hFlowKaon->SetLineColor(colKaon);
  hFlowKaon->SetMarkerColor(colKaon);

  hFlowProton->SetLineColor(colProton);
  hFlowProton->SetMarkerColor(colProton);


  hFlowCharged->Draw();
  hFlowPion->Draw("same");
  hFlowKaon->Draw("same");
  hFlowProton->Draw("same");


  return;
}
