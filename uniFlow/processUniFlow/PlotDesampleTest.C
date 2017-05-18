void PlotDesampleTest()
{
  TString sInputFile = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/sampleModTestNoRebin/UniFlow_DesampleTest.root");

  // colors setting
  Color_t colors[] = {kRed, kGreen, kBlue};

  TFile* fInput = new TFile(sInputFile.Data(),"READ");
  if(!fInput->IsOpen()) return;
  fInput->ls();

  TH1D* hRebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_Refs_Rebin");
  if(!hRebin) return;
  TH1D* hNoRebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_Refs_noRebin");
  if(!hNoRebin) return;
  TH1D* hNoRebin_rebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_testRebin_Refs_noRebin");
  if(!hNoRebin_rebin) return;

  TLegend* leg = new TLegend(0.1,0.6,0.4,0.8);
  leg->AddEntry(hNoRebin, "Desampled","pel");
  leg->AddEntry(hRebin, "Rebin->Desampled","pel");
  leg->AddEntry(hNoRebin_rebin, "Desampled->Rebin","pel");

  hRebin->SetLineColor(colors[0]);
  hNoRebin->SetLineColor(colors[1]);
  hNoRebin_rebin->SetLineColor(colors[2]);

  TH1D* hRatio = (TH1D*) hRebin->Clone("hRatio");
  hRatio->Divide(hNoRebin_rebin);

  // printing test output
  const Short_t testBin = 2;
  printf(" hRebin (%d): %g +- %g\n",testBin,hRebin->GetBinContent(testBin),hRebin->GetBinError(testBin));
  printf(" hNoRebin (%d): %g +- %g\n",testBin,hNoRebin->GetBinContent(testBin),hNoRebin->GetBinError(testBin));
  printf(" hNoRebin_rebin (%d): %g +- %g\n",testBin,hNoRebin_rebin->GetBinContent(testBin),hNoRebin_rebin->GetBinError(testBin));
  printf(" hRatio (%d): %g +- %g\n",testBin,hRatio->GetBinContent(testBin),hRatio->GetBinError(testBin));



  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,1);
  can->cd(1);
  TH1* h2 = can->DrawFrame(0,0.07,100,0.1);
  hRebin->Draw("same");
  hNoRebin->Draw("same");
  hNoRebin_rebin->Draw("same");
  leg->Draw();

  can->cd(2);
  TH1* h = can->DrawFrame(0,0.95,100,1.05);
  hRatio->Draw("same");

}
