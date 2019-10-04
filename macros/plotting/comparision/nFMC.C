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

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/train_7387/processUniFlow/Processed_NoNorm.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_HighOrder.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}

  TGraphErrors *FMC4_harm34 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm34"));
  FMC4_harm34->SetMarkerStyle(kOpenSquare);
  FMC4_harm34->SetMarkerColor(kRed);
  FMC4_harm34->SetMarkerSize(1.);
  FMC4_harm34->SetLineColor(kRed);
  mg->Add(FMC4_harm34);

  TGraphErrors *FMC4_harm34_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC4343"));
  FMC4_harm34_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm34_YZ->SetMarkerColor(kRed);
  FMC4_harm34_YZ->SetMarkerSize(1.);
  FMC4_harm34_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm34_YZ);

  TGraphErrors *FMC4_harm35 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm35"));
  FMC4_harm35->SetMarkerStyle(kOpenCircle);
  FMC4_harm35->SetMarkerColor(kBlue);
  FMC4_harm35->SetMarkerSize(1.1);
  FMC4_harm35->SetMarkerColor(kBlue);
  mg->Add(FMC4_harm35);

  TGraphErrors *FMC4_harm35_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC5353"));
  FMC4_harm35_YZ->SetMarkerStyle(kFullCircle);
  FMC4_harm35_YZ->SetMarkerColor(kBlue);
  FMC4_harm35_YZ->SetMarkerSize(1.1);
  FMC4_harm35_YZ->SetMarkerColor(kBlue);
  mg->Add(FMC4_harm35_YZ);

  TGraphErrors *FMC4_harm45 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm45"));
  FMC4_harm45->SetMarkerStyle(kOpenDiamond);
  FMC4_harm45->SetMarkerColor(kGreen+2);
  FMC4_harm45->SetMarkerSize(1.6);
  FMC4_harm45->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm45);

  TGraphErrors *FMC4_harm45_YZ = new TGraphErrors((TH1D*) fileYZ->Get("his_FMC5454"));
  FMC4_harm45_YZ->SetMarkerStyle(kFullDiamond);
  FMC4_harm45_YZ->SetMarkerColor(kGreen+2);
  FMC4_harm45_YZ->SetMarkerSize(1.6);
  FMC4_harm45_YZ->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm45_YZ);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("SC(m,n)");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-6.5E-8);
  mg->SetMaximum(8.5E-8);

  TLegend* leg = new TLegend(0.12,0.72,0.52,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(FMC4_harm34,"SC(3,4) ","p");
  leg->AddEntry(FMC4_harm34_YZ,"SC(3,4) - YZ","p");
  leg->AddEntry(FMC4_harm35,"SC(3,5)","p");
  leg->AddEntry(FMC4_harm35_YZ,"SC(3,5) - YZ","p");
  leg->AddEntry(FMC4_harm45_YZ,"SC(4,5)","p");
  leg->AddEntry(FMC4_harm45,"SC(4,5) - YZ","p");
  leg->Draw("same");

  can->SaveAs("SC_for345.pdf");

  TCanvas* can2[3];
  TH1D* ratio;
  TH1D* denom;
  TString namesZM[3] = {"Refs_FMC4_harm34", "Refs_FMC4_harm35", "Refs_FMC4_harm45"};
  TString namesYZ[3] = {"his_FMC4343", "his_FMC5353", "his_FMC5454"};
  TString names[3] = {"FMC4343", "FMC5353", "FMC5454"};
  for(int i = 0; i < 3; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    ratio = (TH1D*) fileIn->Get(namesZM[i]);
    denom = (TH1D*) fileYZ->Get(namesYZ[i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality; nSC: ZM / YZ");
    ratio->SetTitle(names[i]);
    // ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetRangeUser(0.,60.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 5, 60);
    ratio->Fit(fu1,"R");
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    can2[i]->SaveAs("ratio_" + names[i]+ ".pdf");
  }


}
