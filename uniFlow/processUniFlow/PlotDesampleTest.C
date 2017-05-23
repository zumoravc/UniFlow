void PlotDesampleTest()
{
  TString sInputFile = TString("/Users/vpacik/NBI/Flow/uniFlow/processUniFlow/sampleModTestNoRebin/UniFlow_DesampleTest_2.root");

  // colors setting
  Color_t colors[] = {kRed, kGreen+2, kBlue, kRed-4,kOrange+1,kBlack,kOrange,kViolet};
  Color_t markers[] = {26,24,21,25,33,24,20,25};

  TFile* fInput = new TFile(sInputFile.Data(),"READ");
  if(!fInput->IsOpen()) return;
  fInput->ls();

  TH1D* hMerged = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_merged_Refs_noRebin");
  if(!hMerged) return;
  TH1D* hMergedRebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_merged_Refs_Rebin");
  if(!hMergedRebin) return;
  TH1D* hMergedRebinTest = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_merged_testRebin_Refs_noRebin");
  if(!hMergedRebinTest) return;
  TH1D* hRebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_Refs_Rebin");
  if(!hRebin) return;
  TH1D* hNoRebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_Refs_noRebin");
  if(!hNoRebin) return;
  TH1D* hNoRebin_unitRMS = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_Refs_noRebin_unitRMS");
  if(!hNoRebin_unitRMS) return;
  TH1D* hNoRebin_rebin = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_testRebin_Refs_noRebin");
  if(!hNoRebin_rebin) return;
  TH1D* hNoRebin_RMS = (TH1D*) fInput->Get("hFlow2_Refs_harm2_gap08_testRebin_Refs_noRebin_RMS");
  if(!hNoRebin_RMS) return;


  TLegend* leg = new TLegend(0.1,0.6,0.4,0.8);
  leg->AddEntry(hMerged, "Merged (true)","pel");
  leg->AddEntry(hNoRebin, "Desampled","pel");
  leg->AddEntry(hNoRebin_unitRMS, "Desampled (RMS)","pel");
  leg->AddEntry(hMergedRebin, "Rebin->Merged","pel");
  leg->AddEntry(hRebin, "Rebin->Desampled","pel");
  leg->AddEntry(hMergedRebinTest, "Merged->Rebin","pel");
  leg->AddEntry(hNoRebin_rebin, "Desampled->Rebin","pel");
  leg->AddEntry(hNoRebin_RMS, "Desampled->Rebin (RMS)","pel");

  hRebin->SetLineColor(colors[0]);
  hRebin->SetMarkerColor(colors[0]);
  hRebin->SetMarkerStyle(markers[0]);
  hNoRebin->SetLineColor(colors[1]);
  hNoRebin->SetMarkerColor(colors[1]);
  hNoRebin->SetMarkerStyle(markers[1]);
  hNoRebin_rebin->SetLineColor(colors[2]);
  hNoRebin_rebin->SetMarkerColor(colors[2]);
  hNoRebin_rebin->SetMarkerStyle(markers[2]);
  hNoRebin_RMS->SetLineColor(colors[3]);
  hNoRebin_RMS->SetMarkerColor(colors[3]);
  hNoRebin_RMS->SetMarkerStyle(markers[3]);
  hMerged->SetLineColor(colors[4]);
  hMerged->SetMarkerColor(colors[4]);
  hMerged->SetMarkerStyle(markers[4]);
  hMergedRebin->SetLineColor(colors[5]);
  hMergedRebin->SetMarkerColor(colors[5]);
  hMergedRebin->SetMarkerStyle(markers[5]);
  hMergedRebinTest->SetLineColor(colors[6]);
  hMergedRebinTest->SetMarkerColor(colors[6]);
  hMergedRebinTest->SetMarkerStyle(markers[6]);
  hNoRebin_unitRMS->SetLineColor(colors[7]);
  hNoRebin_unitRMS->SetMarkerColor(colors[7]);
  hNoRebin_unitRMS->SetMarkerStyle(markers[7]);

  // doing ratios
  TH1D* hRatio_Desampled = (TH1D*) hMerged->Clone("hRatio_Desampled");
  hRatio_Desampled->Divide(hNoRebin);
  hRatio_Desampled->SetLineColor(colors[1]);
  hRatio_Desampled->SetMarkerColor(colors[1]);
  hRatio_Desampled->SetMarkerStyle(markers[1]);

  TH1D* hRatio_unitRMS = (TH1D*) hMerged->Clone("hRatio_unitRMS");
  hRatio_unitRMS->Divide(hNoRebin_unitRMS);
  hRatio_unitRMS->SetLineColor(colors[7]);
  hRatio_unitRMS->SetMarkerColor(colors[7]);
  hRatio_unitRMS->SetMarkerStyle(markers[7]);

  TH1D* hRatio_hMergedRebin = (TH1D*) hMergedRebinTest->Clone("hRatio_hMergedRebin");
  hRatio_hMergedRebin->Divide(hMergedRebin);
  hRatio_hMergedRebin->SetLineColor(colors[5]);
  hRatio_hMergedRebin->SetMarkerColor(colors[5]);
  hRatio_hMergedRebin->SetMarkerStyle(markers[5]);

  TH1D* hRatio_noRebin_rebin = (TH1D*) hMergedRebinTest->Clone("hRatio_noRebin_rebin");
  hRatio_noRebin_rebin->Divide(hNoRebin_rebin);
  hRatio_noRebin_rebin->SetLineColor(colors[2]);
  hRatio_noRebin_rebin->SetMarkerColor(colors[2]);
  hRatio_noRebin_rebin->SetMarkerStyle(markers[2]);

  TH1D* hRatio_hRebin = (TH1D*) hMergedRebinTest->Clone("hRatio_hRebin");
  hRatio_hRebin->Divide(hRebin);
  hRatio_hRebin->SetLineColor(colors[0]);
  hRatio_hRebin->SetMarkerColor(colors[0]);
  hRatio_hRebin->SetMarkerStyle(markers[0]);

  TH1D* hRatio_hNoRebin_RMS = (TH1D*) hMergedRebinTest->Clone("hRatio_hNoRebin_RMS");
  hRatio_hNoRebin_RMS->Divide(hNoRebin_RMS);
  hRatio_hNoRebin_RMS->SetLineColor(colors[3]);
  hRatio_hNoRebin_RMS->SetMarkerColor(colors[3]);
  hRatio_hNoRebin_RMS->SetMarkerStyle(markers[3]);

  // preparing uncertainty comparison
  TH1D* hErr_Desampled = (TH1D*) hNoRebin->Clone("hErr_Desampled");
  hErr_Desampled->Reset();
  for(Short_t bin(1); bin < hNoRebin->GetNbinsX()+1; bin++) { hErr_Desampled->SetBinContent(bin, hNoRebin->GetBinError(bin)); }

  TH1D* hErr_unitRMS = (TH1D*) hNoRebin_unitRMS->Clone("hErr_Desampled");
  hErr_unitRMS->Reset();
  for(Short_t bin(1); bin < hNoRebin_unitRMS->GetNbinsX()+1; bin++) { hErr_unitRMS->SetBinContent(bin, hNoRebin_unitRMS->GetBinError(bin)); }

  TH1D* hErr_hMerged = (TH1D*) hMerged->Clone("hErr_hMerged");
  hErr_hMerged->Reset();
  for(Short_t bin(1); bin < hMerged->GetNbinsX()+1; bin++) { hErr_hMerged->SetBinContent(bin, hMerged->GetBinError(bin)); }

  TH1D* hErr_hMergedRebin = (TH1D*) hMergedRebin->Clone("hErr_hMergedRebin");
  hErr_hMergedRebin->Reset();
  for(Short_t bin(1); bin < hMergedRebin->GetNbinsX()+1; bin++) { hErr_hMergedRebin->SetBinContent(bin, hMergedRebin->GetBinError(bin)); }

  TH1D* hErr_hMergedRebinTest = (TH1D*) hMergedRebinTest->Clone("hErr_hMergedRebin");
  hErr_hMergedRebinTest->Reset();
  for(Short_t bin(1); bin < hMergedRebinTest->GetNbinsX()+1; bin++) { hErr_hMergedRebinTest->SetBinContent(bin, hMergedRebinTest->GetBinError(bin)); }

  TH1D* hErr_hRebin = (TH1D*) hRebin->Clone("hErr_hRebin");
  hErr_hRebin->Reset();
  for(Short_t bin(1); bin < hRebin->GetNbinsX()+1; bin++) { hErr_hRebin->SetBinContent(bin, hRebin->GetBinError(bin)); }

  TH1D* hErr_hNoRebin_rebin = (TH1D*) hNoRebin_rebin->Clone("hErr_hNoRebin_rebin");
  hErr_hNoRebin_rebin->Reset();
  for(Short_t bin(1); bin < hNoRebin_rebin->GetNbinsX()+1; bin++) { hErr_hNoRebin_rebin->SetBinContent(bin, hNoRebin_rebin->GetBinError(bin)); }

  TH1D* hErr_hNoRebin_RMS = (TH1D*) hNoRebin_RMS->Clone("hErr_hNoRebin_RMS");
  hErr_hNoRebin_RMS->Reset();
  for(Short_t bin(1); bin < hNoRebin_RMS->GetNbinsX()+1; bin++) { hErr_hNoRebin_RMS->SetBinContent(bin, hNoRebin_RMS->GetBinError(bin)); }


  // printing test output
  const Short_t testBin = 2;
  printf(" hRebin (%d): %g +- %g\n",testBin,hRebin->GetBinContent(testBin),hRebin->GetBinError(testBin));
  printf(" hNoRebin (%d): %g +- %g\n",testBin,hNoRebin->GetBinContent(testBin),hNoRebin->GetBinError(testBin));
  printf(" hNoRebin_rebin (%d): %g +- %g\n",testBin,hNoRebin_rebin->GetBinContent(testBin),hNoRebin_rebin->GetBinError(testBin));
  // printf(" hRatio (%d): %g +- %g\n",testBin,hRatio->GetBinContent(testBin),hRatio->GetBinError(testBin));



  TCanvas* can = new TCanvas("can","can",1200,400);
  can->Divide(3,1);
  can->cd(1);
  TH1* h2 = can->DrawFrame(0,0.07,100,0.1);
  h2->SetTitle("v_{2} {2}");
  hMerged->Draw("same");
  hNoRebin->Draw("same");
  hNoRebin_unitRMS->Draw("same");
  hMergedRebin->Draw("same");
  hMergedRebinTest->Draw("same");
  hRebin->Draw("same");
  hNoRebin_rebin->Draw("same");
  hNoRebin_RMS->Draw("same");
  leg->Draw();

  can->cd(2);
  TH1* h = can->DrawFrame(0,0.99,100,1.01);
  h->SetTitle("Ratio Merged (Merged->Rebin) / X");
  hRatio_Desampled->Draw("same hist p");
  hRatio_unitRMS->Draw("same hist p");
  hRatio_hRebin->Draw("same hist p");
  hRatio_hMergedRebin->Draw("same hist p");
  hRatio_noRebin_rebin->Draw("same hist p");
  hRatio_hNoRebin_RMS->Draw("same hist p");

  can->cd(3);
  TH1* h3 = can->DrawFrame(0,0.,100,0.005);
  h3->SetTitle("Uncertainty");
  hErr_Desampled->Draw("same hist p");
  hErr_unitRMS->Draw("same hist p");
  hErr_hMerged->Draw("same hist p");
  hErr_hRebin->Draw("same hist p");
  hErr_hNoRebin_rebin->Draw("same hist p");
  hErr_hMergedRebin->Draw("same hist p");
  hErr_hMergedRebinTest->Draw("same hist p");
  hErr_hNoRebin_RMS->Draw("same hist p");
}
