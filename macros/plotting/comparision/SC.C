void SC()
{
  TMultiGraph *mg = new TMultiGraph();

  Color_t colors[9] = {kBlack, kRed+2, kBlue, kGreen+2, kMagenta+2, kBlue+3, kYellow+2, kCyan +2, kGray+2};
  TFile* fileIn[9];
  TGraphErrors *FMC4_harm23[9];
  for(int i = 0; i < 9; i++)
  {
    fileIn[i] =  TFile::Open(Form("/home/alidock/ana/output/LHC15o/train_7364/stage1/processUniFlow/Processed%d.root",i+1));
    if(!fileIn[i]) {printf("File %d not opened! \n", i+1); return;}
    FMC4_harm23[i] = new TGraphErrors((TH1D*) fileIn[i]->Get("Refs_FMC4_harm23"));
    mg->Add(FMC4_harm23[i]);
    FMC4_harm23[i]->SetMarkerColor(colors[i]);
    // FMC4_harm23[i]->SetMarkerStyle(kOpenSquare);
    FMC4_harm23[i]->SetLineColor(colors[i]);
  }

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}

  TGraphErrors *FMC4_harm23_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC3232"));
  FMC4_harm23_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm23_YZ->SetMarkerColor(kRed);
  FMC4_harm23_YZ->SetMarkerSize(1.);
  FMC4_harm23_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm23_YZ);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("SC (2,3)");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-0.28);
  mg->SetMaximum(0.62);

  TLegend* leg = new TLegend(0.32,0.7,0.72,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  // leg->SetNColumns(2);
  leg->AddEntry(FMC4_harm23_YZ,"SC(2,3) - YZ","p");
  leg->Draw("same");

  can->SaveAs("SC23.pdf");

}
