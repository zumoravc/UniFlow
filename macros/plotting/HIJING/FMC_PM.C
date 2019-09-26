// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void FMC_PM()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn_p = TFile::Open("/home/alidock/ana/output/MC/HIJING_grid_full_RBR/processUniFlow/Plus.root","READ");
  if(!fileIn_p) {printf("File plus not opened! \n"); return;}

  TFile* fileIn_m = TFile::Open("/home/alidock/ana/output/MC/HIJING_grid_full_RBR/processUniFlow/Minus.root","READ");
  if(!fileIn_m) {printf("File minus not opened! \n"); return;}

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/MC/HIJING_grid_full_RBR/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File processed not opened! \n"); return;}

  TGraphErrors *FMC4_harm23 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm23"));
  FMC4_harm23->SetMarkerStyle(kOpenSquare);
  FMC4_harm23->SetMarkerColor(kRed);
  FMC4_harm23->SetMarkerSize(1.0);
  FMC4_harm23->SetLineColor(kRed);
  mg->Add(FMC4_harm23);

  TGraphErrors *FMC4_harm23_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_FMC4_harm23"));
  FMC4_harm23_p->SetMarkerStyle(kFullSquare);
  FMC4_harm23_p->SetMarkerColor(kGreen+2);
  FMC4_harm23_p->SetMarkerSize(1.);
  FMC4_harm23_p->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm23_p);

  TGraphErrors *FMC4_harm23_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_FMC4_harm23"));
  FMC4_harm23_m->SetMarkerStyle(kFullSquare);
  FMC4_harm23_m->SetMarkerColor(kBlue+1);
  FMC4_harm23_m->SetMarkerSize(1.);
  FMC4_harm23_m->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm23_m);

  TGraphErrors *FMC4_harm24 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm24"));
  FMC4_harm24->SetMarkerStyle(kOpenCircle);
  FMC4_harm24->SetMarkerColor(kRed);
  FMC4_harm24->SetMarkerSize(1.0);
  FMC4_harm24->SetLineColor(kRed);
  mg->Add(FMC4_harm24);

  TGraphErrors *FMC4_harm24_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_FMC4_harm24"));
  FMC4_harm24_p->SetMarkerStyle(kFullCircle);
  FMC4_harm24_p->SetMarkerColor(kGreen+2);
  FMC4_harm24_p->SetMarkerSize(1.);
  FMC4_harm24_p->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm24_p);

  TGraphErrors *FMC4_harm24_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_FMC4_harm24"));
  FMC4_harm24_m->SetMarkerStyle(kFullCircle);
  FMC4_harm24_m->SetMarkerColor(kBlue+1);
  FMC4_harm24_m->SetMarkerSize(1.);
  FMC4_harm24_m->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24_m);

  TGraphErrors *FMC4_harm34 = new TGraphErrors((TH1D*) fileIn->Get("Refs_FMC4_harm34"));
  FMC4_harm34->SetMarkerStyle(kOpenDiamond);
  FMC4_harm34->SetMarkerColor(kRed);
  FMC4_harm34->SetMarkerSize(1.0);
  FMC4_harm34->SetLineColor(kRed);
  mg->Add(FMC4_harm34);

  TGraphErrors *FMC4_harm34_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_FMC4_harm34"));
  FMC4_harm34_p->SetMarkerStyle(kFullDiamond);
  FMC4_harm34_p->SetMarkerColor(kGreen+2);
  FMC4_harm34_p->SetMarkerSize(1.);
  FMC4_harm34_p->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34_p);

  TGraphErrors *FMC4_harm34_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_FMC4_harm34"));
  FMC4_harm34_m->SetMarkerStyle(kFullDiamond);
  FMC4_harm34_m->SetMarkerColor(kBlue+1);
  FMC4_harm34_m->SetMarkerSize(1.);
  FMC4_harm34_m->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm34_m);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("SC(m,n)");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-0.15E-6);
  mg->SetMaximum(0.25E-6);

  TLegend* leg = new TLegend(0.32,0.72,0.72,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("HIJING, LHC18e1, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(3);
  leg->AddEntry(FMC4_harm23,"SC(2,3)","p");
  leg->AddEntry(FMC4_harm23_p,"SC(2,3) - Plus","p");
  leg->AddEntry(FMC4_harm23_m,"SC(2,3) - Minus","p");
  leg->AddEntry(FMC4_harm24,"SC(2,4)","p");
  leg->AddEntry(FMC4_harm24_p,"SC(2,4) - Plus","p");
  leg->AddEntry(FMC4_harm24_m,"SC(2,4) - Minus","p");
  leg->AddEntry(FMC4_harm34,"SC(3,4)","p");
  leg->AddEntry(FMC4_harm34_p,"SC(3,4) - Plus","p");
  leg->AddEntry(FMC4_harm34_m,"SC(3,4) - Minus","p");
  leg->Draw("same");

  can->SaveAs("FMC4_comp_PM.pdf");

  TCanvas* can2[6];
  TH1D* ratio;
  TH1D* denom;
  TString name[6] = {"harm23", "harm24", "harm34", "harm23", "harm24", "harm34"};
  for(int i = 0; i < 6; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    if(i < 3) ratio = (TH1D*) fileIn_p->Get("Refs_FMC4_" + name[i]);
    else ratio = (TH1D*) fileIn_m->Get("Refs_FMC4_" + name[i]);
    denom = (TH1D*) fileIn->Get("Refs_FMC4_" + name[i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality");
    if(i < 3) ratio->SetTitle("pos " + name[i]);
    else ratio->SetTitle("neg " + name[i]);
    ratio->GetYaxis()->SetRangeUser(-1.,3.);
    ratio->GetXaxis()->SetRangeUser(0,60);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0, 60);
    ratio->Fit(fu1);
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    if(i<3) can2[i]->SaveAs("ratio_p_" + name[i]+ ".pdf");
    else can2[i]->SaveAs("ratio_m_" + name[i]+ ".pdf");
  }


}
