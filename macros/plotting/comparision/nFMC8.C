// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void nFMC8()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}

  TGraphErrors *fMC8_harm2223 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2223"));
  fMC8_harm2223->SetMarkerStyle(kOpenSquare);
  fMC8_harm2223->SetMarkerColor(kRed);
  fMC8_harm2223->SetMarkerSize(1.);
  fMC8_harm2223->SetLineColor(kRed);
  mg->Add(fMC8_harm2223);

  TGraphErrors *fMC8_harm2223_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC32223222"));
  fMC8_harm2223_YZ->SetMarkerStyle(kFullSquare);
  fMC8_harm2223_YZ->SetMarkerColor(kRed);
  fMC8_harm2223_YZ->SetMarkerSize(1.);
  fMC8_harm2223_YZ->SetLineColor(kRed);
  mg->Add(fMC8_harm2223_YZ);

  TGraphErrors *fMC8_harm2233 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2233"));
  fMC8_harm2233->SetMarkerStyle(kOpenCircle);
  fMC8_harm2233->SetMarkerColor(kBlue);
  fMC8_harm2233->SetMarkerSize(1.1);
  fMC8_harm2233->SetMarkerColor(kBlue);
  mg->Add(fMC8_harm2233);

  TGraphErrors *fMC8_harm2233_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC33223322"));
  fMC8_harm2233_YZ->SetMarkerStyle(kFullCircle);
  fMC8_harm2233_YZ->SetMarkerColor(kBlue);
  fMC8_harm2233_YZ->SetMarkerSize(1.1);
  fMC8_harm2233_YZ->SetMarkerColor(kBlue);
  mg->Add(fMC8_harm2233_YZ);

  TGraphErrors *fMC8_harm2333 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2333"));
  fMC8_harm2333->SetMarkerStyle(kOpenDiamond);
  fMC8_harm2333->SetMarkerColor(kGreen+2);
  fMC8_harm2333->SetMarkerSize(1.6);
  fMC8_harm2333->SetLineColor(kGreen+2);
  mg->Add(fMC8_harm2333);

  TGraphErrors *fMC8_harm2333_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC33323332"));
  fMC8_harm2333_YZ->SetMarkerStyle(kFullDiamond);
  fMC8_harm2333_YZ->SetMarkerColor(kGreen+2);
  fMC8_harm2333_YZ->SetMarkerSize(1.6);
  fMC8_harm2333_YZ->SetLineColor(kGreen+2);
  mg->Add(fMC8_harm2333_YZ);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("nFMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-0.4);
  mg->SetMaximum(0.2);

  TLegend* leg = new TLegend(0.32,0.72,0.72,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, LHC15o (sample), additional pile up rejection, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(fMC8_harm2223_YZ,"nFMC (v_{2}^{6},v_{3}^{2}) - YZ","p");
  leg->AddEntry(fMC8_harm2223,"nFMC (v_{2}^{6},v_{3}^{2})","p");
  leg->AddEntry(fMC8_harm2233_YZ,"nFMC (v_{2}^{4},v_{3}^{4}) - YZ","p");
  leg->AddEntry(fMC8_harm2233,"nFMC (v_{2}^{4},v_{3}^{4})","p");
  leg->AddEntry(fMC8_harm2333_YZ,"nFMC (v_{2}^{2},v_{3}^{6}) - YZ","p");
  leg->AddEntry(fMC8_harm2333,"nFMC (v_{2}^{2},v_{3}^{6})","p");
  leg->Draw("same");

  can->SaveAs("FMC8_comp+new.pdf");

  TCanvas* can2 = new TCanvas("can2", "can2", 600, 400);
  gStyle->SetOptStat(kFALSE);
  TH1D *ratio = (TH1D*) fileIn->Get("Refs_fMC8_harm2223");
  TH1D *denom = (TH1D*) fileYZ->Get("his_nFMC32223222");
  ratio->Divide(denom);
  ratio->SetMarkerStyle(kFullSquare);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);
  ratio->SetTitle(";Centrality; nFMC (v_{2}^{6},v_{3}^{2}) YZ / nFMC (v_{2}^{6},v_{3}^{2})");
  ratio->Draw("h");
  TF1 *fu1 = new TF1("fu1", "1", 0, 100);
  fu1->SetLineColor(kRed);
  fu1->Draw("same");
  can2->SaveAs("ratio_nFMC2223.pdf");


}
