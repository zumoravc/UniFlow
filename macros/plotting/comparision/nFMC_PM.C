// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void nFMC_PM()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn_p = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Plus.root","READ");
  if(!fileIn_p) {printf("File plus not opened! \n"); return;}

  TFile* fileIn_m = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Minus.root","READ");
  if(!fileIn_m) {printf("File minus not opened! \n"); return;}

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File processed not opened! \n"); return;}

  TGraphErrors *FMC4_harm23_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_FMC4_harm23"));
  FMC4_harm23_p->SetMarkerStyle(kOpenSquare);
  FMC4_harm23_p->SetMarkerColor(kRed);
  FMC4_harm23_p->SetMarkerSize(1.);
  FMC4_harm23_p->SetLineColor(kRed);
  mg->Add(FMC4_harm23_p);

  TGraphErrors *FMC4_harm23_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_FMC4_harm23"));
  FMC4_harm23_m->SetMarkerStyle(kFullSquare);
  FMC4_harm23_m->SetMarkerColor(kRed);
  FMC4_harm23_m->SetMarkerSize(1.);
  FMC4_harm23_m->SetLineColor(kRed);
  mg->Add(FMC4_harm23_m);

  TGraphErrors *fMC6_harm223_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC6_harm223"));
  fMC6_harm223_p->SetMarkerStyle(kOpenCircle);
  fMC6_harm223_p->SetMarkerColor(kBlue);
  fMC6_harm223_p->SetMarkerSize(1.1);
  fMC6_harm223_p->SetMarkerColor(kBlue);
  mg->Add(fMC6_harm223_p);

  TGraphErrors *fMC6_harm223_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC6_harm223"));
  fMC6_harm223_m->SetMarkerStyle(kFullCircle);
  fMC6_harm223_m->SetMarkerColor(kBlue);
  fMC6_harm223_m->SetMarkerSize(1.1);
  fMC6_harm223_m->SetMarkerColor(kBlue);
  mg->Add(fMC6_harm223_m);

  TGraphErrors *fMC6_harm233_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC6_harm233"));
  fMC6_harm233_p->SetMarkerStyle(kOpenDiamond);
  fMC6_harm233_p->SetMarkerColor(kGreen+2);
  fMC6_harm233_p->SetMarkerSize(1.6);
  fMC6_harm233_p->SetLineColor(kGreen+2);
  mg->Add(fMC6_harm233_p);

  TGraphErrors *fMC6_harm233_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC6_harm233"));
  fMC6_harm233_m->SetMarkerStyle(kFullDiamond);
  fMC6_harm233_m->SetMarkerColor(kGreen+2);
  fMC6_harm233_m->SetMarkerSize(1.6);
  fMC6_harm233_m->SetLineColor(kGreen+2);
  mg->Add(fMC6_harm233_m);

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
  leg->SetHeader("ALICE experiment, LHC15o, additional pile up rejection, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(FMC4_harm23_p,"nFMC (v_{2}^{2},v_{3}^{2}) - Plus","p");
  leg->AddEntry(FMC4_harm23_m,"nFMC (v_{2}^{2},v_{3}^{2}) - Minus","p");
  leg->AddEntry(fMC6_harm223_p,"nFMC (v_{2}^{4},v_{3}^{2}) - Plus","p");
  leg->AddEntry(fMC6_harm223_m,"nFMC (v_{2}^{4},v_{3}^{2}) - Minus","p");
  leg->AddEntry(fMC6_harm233_p,"nFMC (v_{2}^{2},v_{3}^{4}) - Plus","p");
  leg->AddEntry(fMC6_harm233_m,"nFMC (v_{2}^{2},v_{3}^{4}) - Minus","p");
  leg->Draw("same");

  can->SaveAs("FMC_comp_PM.pdf");

  TCanvas* can2[6];
  TH1D* ratio;
  TH1D* denom;
  TString name[6] = {"harm23", "harm223", "harm233", "harm23", "harm223", "harm233"};
  for(int i = 0; i < 6; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    if(i == 0) ratio = (TH1D*) fileIn_p->Get("Refs_FMC4_" + name[i]);
    else if(i == 3) ratio = (TH1D*) fileIn_m->Get("Refs_FMC4_" + name[i]);
    else if(i < 3) ratio = (TH1D*) fileIn_p->Get("Refs_fMC6_" + name[i]);
    else ratio = (TH1D*) fileIn_m->Get("Refs_fMC6_" + name[i]);
    if(i == 0 || i == 3) denom = (TH1D*) fileIn->Get("Refs_FMC4_" + name[i]);
    else denom = (TH1D*) fileIn->Get("Refs_fMC6_" + name[i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality; nFMC (v_{2}^{k},v_{3}^{l}) / nFMC (v_{2}^{k},v_{3}^{l})");
    if(i < 3) ratio->SetTitle("pos " + name[i]);
    else ratio->SetTitle("neg " + name[i]);
    ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0, 100);
    ratio->Fit(fu1);
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    if(i<3) can2[i]->SaveAs("ratio_p_" + name[i]+ ".pdf");
    else can2[i]->SaveAs("ratio_m_" + name[i]+ ".pdf");
  }


}
