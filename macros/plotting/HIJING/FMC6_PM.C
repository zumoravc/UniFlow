// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void FMC6_PM()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn_p = TFile::Open("/home/alidock/ana/output/MC/HIJING_grid_full_RBR/processUniFlow/Plus.root","READ");
  if(!fileIn_p) {printf("File plus not opened! \n"); return;}

  TFile* fileIn_m = TFile::Open("/home/alidock/ana/output/MC/HIJING_grid_full_RBR/processUniFlow/Minus.root","READ");
  if(!fileIn_m) {printf("File minus not opened! \n"); return;}

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/MC/HIJING_grid_full_RBR/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File processed not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}


  TGraphErrors *fMC6_harm223 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm223"));
  fMC6_harm223->SetMarkerStyle(kOpenSquare);
  fMC6_harm223->SetMarkerColor(kRed);
  fMC6_harm223->SetMarkerSize(1.0);
  fMC6_harm223->SetLineColor(kRed);
  mg->Add(fMC6_harm223);

  TGraphErrors *FMC4_harm23_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC322322"));
  FMC4_harm23_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm23_YZ->SetMarkerColor(kRed);
  FMC4_harm23_YZ->SetMarkerSize(1.);
  FMC4_harm23_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm23_YZ);

  TGraphErrors *fMC6_harm233 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm233"));
  fMC6_harm233->SetMarkerStyle(kOpenCircle);
  fMC6_harm233->SetMarkerColor(kBlue+1);
  fMC6_harm233->SetMarkerSize(1.0);
  fMC6_harm233->SetLineColor(kBlue+1);
  mg->Add(fMC6_harm233);

  TGraphErrors *fMC6_harm234 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC6_harm234"));
  fMC6_harm234->SetMarkerStyle(kOpenDiamond);
  fMC6_harm234->SetMarkerColor(kGreen+2);
  fMC6_harm234->SetMarkerSize(1.0);
  fMC6_harm234->SetLineColor(kGreen+2);
  mg->Add(fMC6_harm234);


  TGraphErrors *FMC4_harm24_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC332332"));
  FMC4_harm24_YZ->SetMarkerStyle(kFullCircle);
  FMC4_harm24_YZ->SetMarkerColor(kBlue+1);
  FMC4_harm24_YZ->SetMarkerSize(1.);
  FMC4_harm24_YZ->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24_YZ);


  TGraphErrors *FMC4_harm34_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC432432"));
  FMC4_harm34_YZ->SetMarkerStyle(kFullDiamond);
  FMC4_harm34_YZ->SetMarkerColor(kGreen+2);
  FMC4_harm34_YZ->SetMarkerSize(1.);
  FMC4_harm34_YZ->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34_YZ);

  // TGraphErrors *fMC6_harm223_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC6_harm223"));
  // fMC6_harm223_p->SetMarkerStyle(kOpenCircle);
  // fMC6_harm223_p->SetMarkerColor(kBlue);
  // fMC6_harm223_p->SetMarkerSize(1.1);
  // fMC6_harm223_p->SetMarkerColor(kBlue);
  // mg->Add(fMC6_harm223_p);
  //
  // TGraphErrors *fMC6_harm223_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC6_harm223"));
  // fMC6_harm223_m->SetMarkerStyle(kFullCircle);
  // fMC6_harm223_m->SetMarkerColor(kBlue);
  // fMC6_harm223_m->SetMarkerSize(1.1);
  // fMC6_harm223_m->SetMarkerColor(kBlue);
  // mg->Add(fMC6_harm223_m);
  //
  // TGraphErrors *fMC6_harm233_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC6_harm233"));
  // fMC6_harm233_p->SetMarkerStyle(kOpenDiamond);
  // fMC6_harm233_p->SetMarkerColor(kGreen+2);
  // fMC6_harm233_p->SetMarkerSize(1.6);
  // fMC6_harm233_p->SetLineColor(kGreen+2);
  // mg->Add(fMC6_harm233_p);
  //
  // TGraphErrors *fMC6_harm233_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC6_harm233"));
  // fMC6_harm233_m->SetMarkerStyle(kFullDiamond);
  // fMC6_harm233_m->SetMarkerColor(kGreen+2);
  // fMC6_harm233_m->SetMarkerSize(1.6);
  // fMC6_harm233_m->SetLineColor(kGreen+2);
  // mg->Add(fMC6_harm233_m);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("FMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-3E-9);
  mg->SetMaximum(30E-9);

  TLegend* leg = new TLegend(0.12,0.72,0.52,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("HIJING, LHC18e1, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(fMC6_harm223,"FMC (v_{2}^{4},v_{3}^{2}) - HIJING","p");
  leg->AddEntry(FMC4_harm23_YZ,"FMC (v_{2}^{4},v_{3}^{2}) - YZ","p");
  leg->AddEntry(fMC6_harm233,"FMC (v_{2}^{2},v_{3}^{4}) - HIJING","p");
  leg->AddEntry(FMC4_harm24_YZ,"FMC (v_{2}^{2},v_{3}^{4}) - YZ","p");
  leg->AddEntry(fMC6_harm234,"FMC (v_{2}^{2},v_{3}^{2},v_{4}^{2}) - HIJING","p");
  leg->AddEntry(FMC4_harm34_YZ,"FMC (v_{2}^{2},v_{3}^{2},v_{4}^{2}) - YZ","p");

  leg->Draw("same");

  can->SaveAs("FMC6.pdf");

  TCanvas* can2[6];
  TH1D* ratio;
  TH1D* denom;
  TString name[6] = {"harm223", "harm233", "harm234", "322322", "332332", "432432"};
  for(int i = 0; i < 3; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    ratio = (TH1D*) fileIn->Get("Refs_fMC6_" + name[i]);
    denom = (TH1D*) fileYZ->Get("his_FMC" + name[i+3]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality; FMC6 HIJING / data");
    // ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetRangeUser(0.,60.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0, 60);
    ratio->Fit(fu1);
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    can2[i]->SaveAs("ratio_" + name[i]+ ".pdf");
  }


}
