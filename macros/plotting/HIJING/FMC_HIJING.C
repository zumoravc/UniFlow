// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void FMC_HIJING()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/MC/AMPT_train2380/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}


  TGraphErrors *FMC4_harm23 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm23"));
  FMC4_harm23->SetMarkerStyle(kOpenSquare);
  FMC4_harm23->SetMarkerColor(kRed);
  FMC4_harm23->SetMarkerSize(1.0);
  FMC4_harm23->SetLineColor(kRed);
  mg->Add(FMC4_harm23);

  TGraphErrors *FMC4_harm23_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC3232"));
  FMC4_harm23_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm23_YZ->SetMarkerColor(kRed);
  FMC4_harm23_YZ->SetMarkerSize(1.);
  FMC4_harm23_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm23_YZ);

  TGraphErrors *FMC4_harm24 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm24"));
  FMC4_harm24->SetMarkerStyle(kOpenCircle);
  FMC4_harm24->SetMarkerColor(kBlue+1);
  FMC4_harm24->SetMarkerSize(1.0);
  FMC4_harm24->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24);

  TGraphErrors *FMC4_harm24_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC4242"));
  FMC4_harm24_YZ->SetMarkerStyle(kFullCircle);
  FMC4_harm24_YZ->SetMarkerColor(kBlue+1);
  FMC4_harm24_YZ->SetMarkerSize(1.);
  FMC4_harm24_YZ->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24_YZ);

  TGraphErrors *FMC4_harm34 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm34"));
  FMC4_harm34->SetMarkerStyle(kOpenDiamond);
  FMC4_harm34->SetMarkerColor(kGreen+2);
  FMC4_harm34->SetMarkerSize(1.0);
  FMC4_harm34->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34);

  TGraphErrors *FMC4_harm34_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC4343"));
  FMC4_harm34_YZ->SetMarkerStyle(kFullDiamond);
  FMC4_harm34_YZ->SetMarkerColor(kGreen+2);
  FMC4_harm34_YZ->SetMarkerSize(1.);
  FMC4_harm34_YZ->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34_YZ);

  // TGraphErrors *fMC6_harm223 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm223"));
  // fMC6_harm223->SetMarkerStyle(kOpenCircle);
  // fMC6_harm223->SetMarkerColor(3);
  // fMC6_harm223->SetMarkerSize(1.0);
  // fMC6_harm223->SetLineColor(3);
  // mg->Add(fMC6_harm223);
  //
  // TGraphErrors *fMC6_harm233 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm233"));
  // fMC6_harm233->SetMarkerStyle(kOpenDiamond);
  // fMC6_harm233->SetMarkerColor(4);
  // fMC6_harm233->SetMarkerSize(1.0);
  // fMC6_harm233->SetLineColor(4);
  // mg->Add(fMC6_harm233);
  //
  // TGraphErrors *fMC8_harm2223 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2223"));
  // fMC8_harm2223->SetMarkerStyle(kOpenTriangleUp);
  // fMC8_harm2223->SetMarkerColor(1);
  // fMC8_harm2223->SetMarkerSize(1.0);
  // fMC8_harm2223->SetLineColor(1);
  // mg->Add(fMC8_harm2223);
  //
  // TGraphErrors *fMC8_harm2233 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2233"));
  // fMC8_harm2233->SetMarkerStyle(kOpenTriangleDown);
  // fMC8_harm2233->SetMarkerColor(6);
  // fMC8_harm2233->SetMarkerSize(1.0);
  // fMC8_harm2233->SetMarkerColor(6);
  // mg->Add(fMC8_harm2233);
  //
  // TGraphErrors *fMC8_harm2333 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2333"));
  // fMC8_harm2333->SetMarkerStyle(kOpenCross);
  // fMC8_harm2333->SetMarkerColor(7);
  // fMC8_harm2333->SetMarkerSize(1.0);
  // fMC8_harm2333->SetLineColor(7);
  // mg->Add(fMC8_harm2333);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("FMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-3E-6);
  mg->SetMaximum(5E-6);

  TLegend* leg = new TLegend(0.12,0.72,0.52,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(FMC4_harm23,"SC (2,3) - AMPT","p");
  leg->AddEntry(FMC4_harm23_YZ,"SC (2,3) - YZ","p");
  leg->AddEntry(FMC4_harm24,"SC (2,4) - AMPT","p");
  leg->AddEntry(FMC4_harm24_YZ,"SC (2,4) - YZ","p");
  leg->AddEntry(FMC4_harm34,"SC (3,4) - AMPT","p");
  leg->AddEntry(FMC4_harm34_YZ,"SC (3,4) - YZ","p");
  // leg->AddEntry(fMC6_harm223,"FMC (v_{2}^{4},v_{3}^{2})","p");
  // leg->AddEntry(fMC6_harm233,"FMC (v_{2}^{2},v_{3}^{4})","p");
  // leg->AddEntry(fMC8_harm2223,"FMC (v_{2}^{6},v_{3}^{2})","p");
  // leg->AddEntry(fMC8_harm2233,"FMC (v_{2}^{4},v_{3}^{4})","p");
  // leg->AddEntry(fMC8_harm2333,"FMC (v_{2}^{2},v_{3}^{6})","p");

  leg->Draw("same");

  can->SaveAs("SC_AMPT.pdf");



}
