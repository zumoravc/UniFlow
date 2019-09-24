// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void nFMC()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}

  TGraphErrors *FMC4_harm23 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm23"));
  FMC4_harm23->SetMarkerStyle(kOpenSquare);
  FMC4_harm23->SetMarkerColor(kRed);
  FMC4_harm23->SetMarkerSize(1.);
  FMC4_harm23->SetLineColor(kRed);
  mg->Add(FMC4_harm23);

  TGraphErrors *FMC4_harm23_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC3232"));
  FMC4_harm23_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm23_YZ->SetMarkerColor(kRed);
  FMC4_harm23_YZ->SetMarkerSize(1.);
  FMC4_harm23_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm23_YZ);

  TGraphErrors *fMC6_harm223 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm223"));
  fMC6_harm223->SetMarkerStyle(kOpenCircle);
  fMC6_harm223->SetMarkerColor(kBlue);
  fMC6_harm223->SetMarkerSize(1.1);
  fMC6_harm223->SetMarkerColor(kBlue);
  mg->Add(fMC6_harm223);

  TGraphErrors *fMC6_harm223_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC322322"));
  fMC6_harm223_YZ->SetMarkerStyle(kFullCircle);
  fMC6_harm223_YZ->SetMarkerColor(kBlue);
  fMC6_harm223_YZ->SetMarkerSize(1.1);
  fMC6_harm223_YZ->SetMarkerColor(kBlue);
  mg->Add(fMC6_harm223_YZ);

  TGraphErrors *fMC6_harm233 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm233"));
  fMC6_harm233->SetMarkerStyle(kOpenDiamond);
  fMC6_harm233->SetMarkerColor(kGreen+2);
  fMC6_harm233->SetMarkerSize(1.6);
  fMC6_harm233->SetLineColor(kGreen+2);
  mg->Add(fMC6_harm233);

  TGraphErrors *fMC6_harm233_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_nFMC332332"));
  fMC6_harm233_YZ->SetMarkerStyle(kFullDiamond);
  fMC6_harm233_YZ->SetMarkerColor(kGreen+2);
  fMC6_harm233_YZ->SetMarkerSize(1.6);
  fMC6_harm233_YZ->SetLineColor(kGreen+2);
  mg->Add(fMC6_harm233_YZ);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("nFMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-0.25);
  mg->SetMaximum(0.3);

  TLegend* leg = new TLegend(0.32,0.72,0.72,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, LHC15o (sample), additional pile up rejection, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(FMC4_harm23_YZ,"nFMC (v_{2}^{2},v_{3}^{2}) - YZ","p");
  leg->AddEntry(FMC4_harm23,"nFMC (v_{2}^{2},v_{3}^{2})","p");
  leg->AddEntry(fMC6_harm223_YZ,"nFMC (v_{2}^{4},v_{3}^{2}) - YZ","p");
  leg->AddEntry(fMC6_harm223,"nFMC (v_{2}^{4},v_{3}^{2})","p");
  leg->AddEntry(fMC6_harm233_YZ,"nFMC (v_{2}^{2},v_{3}^{4}) - YZ","p");
  leg->AddEntry(fMC6_harm233,"nFMC (v_{2}^{2},v_{3}^{4})","p");
  leg->Draw("same");

  can->SaveAs("FMC_comp_new.pdf");

  // TCanvas* can2 = new TCanvas("can2", "can2", 600, 400);
  // gStyle->SetOptStat(kFALSE);
  // TH1D *data_histov22 = (TH1D*) fileIn->Get("Refs_hFlow4_harm2_gap00");
  // TH1D *published_v22 = (TH1D*)data_histov22->Clone("published_v22");
  // for(int i = 1; i < 10; i++)
  // {
  //   published_v22->SetBinContent(i,v24Run2[i-1]);
  //   published_v22->SetBinError(i,v24Run2CombErr[i-1]);
  // }
  // published_v22->Divide(data_histov22);
  // published_v22->SetMarkerStyle(kFullSquare);
  // published_v22->SetMarkerColor(kBlack);
  // published_v22->SetLineColor(kBlack);
  // published_v22->SetTitle(";Centrality; v_{2}{4} published / v_{2}{4} reconstructed");
  // published_v22->Draw("h");
  // TF1 *fu1 = new TF1("fu1", "1", 0, 100);
  // fu1->SetLineColor(kRed);
  // fu1->Draw("same");
  // can2->SaveAs("ratio_v24.pdf");


}
