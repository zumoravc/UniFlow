// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void nFMC8_PM()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn_p = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Plus.root","READ");
  if(!fileIn_p) {printf("File plus not opened! \n"); return;}

  TFile* fileIn_m = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Minus.root","READ");
  if(!fileIn_m) {printf("File minus not opened! \n"); return;}

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File processed not opened! \n"); return;}

  TGraphErrors *fMC8_harm2223_p= new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC8_harm2223"));
  fMC8_harm2223_p->SetMarkerStyle(kOpenSquare);
  fMC8_harm2223_p->SetMarkerColor(kRed);
  fMC8_harm2223_p->SetMarkerSize(1.);
  fMC8_harm2223_p->SetLineColor(kRed);
  mg->Add(fMC8_harm2223_p);

  TGraphErrors *fMC8_harm2223_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC8_harm2223"));
  fMC8_harm2223_m->SetMarkerStyle(kFullSquare);
  fMC8_harm2223_m->SetMarkerColor(kRed);
  fMC8_harm2223_m->SetMarkerSize(1.);
  fMC8_harm2223_m->SetLineColor(kRed);
  mg->Add(fMC8_harm2223_m);

  TGraphErrors *fMC8_harm2233_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC8_harm2233"));
  fMC8_harm2233_p->SetMarkerStyle(kOpenCircle);
  fMC8_harm2233_p->SetMarkerColor(kBlue);
  fMC8_harm2233_p->SetMarkerSize(1.1);
  fMC8_harm2233_p->SetMarkerColor(kBlue);
  mg->Add(fMC8_harm2233_p);

  TGraphErrors *fMC8_harm2233_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC8_harm2233"));
  fMC8_harm2233_m->SetMarkerStyle(kFullCircle);
  fMC8_harm2233_m->SetMarkerColor(kBlue);
  fMC8_harm2233_m->SetMarkerSize(1.1);
  fMC8_harm2233_m->SetMarkerColor(kBlue);
  mg->Add(fMC8_harm2233_m);

  TGraphErrors *fMC8_harm2333_p = new TGraphErrors((TH1D*) fileIn_p->Get("Refs_fMC8_harm2333"));
  fMC8_harm2333_p->SetMarkerStyle(kOpenDiamond);
  fMC8_harm2333_p->SetMarkerColor(kGreen+2);
  fMC8_harm2333_p->SetMarkerSize(1.6);
  fMC8_harm2333_p->SetLineColor(kGreen+2);
  mg->Add(fMC8_harm2333_p);

  TGraphErrors *fMC8_harm2333_m = new TGraphErrors((TH1D*) fileIn_m->Get("Refs_fMC8_harm2333"));
  fMC8_harm2333_m->SetMarkerStyle(kFullDiamond);
  fMC8_harm2333_m->SetMarkerColor(kGreen+2);
  fMC8_harm2333_m->SetMarkerSize(1.6);
  fMC8_harm2333_m->SetLineColor(kGreen+2);
  mg->Add(fMC8_harm2333_m);

  // TGraphErrors *fMC8_harm2223 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2223"));
  // fMC8_harm2223->SetMarkerStyle(kOpenSquare);
  // fMC8_harm2223->SetMarkerColor(kRed);
  // fMC8_harm2223->SetMarkerSize(1.);
  // fMC8_harm2223->SetLineColor(kRed);
  // mg->Add(fMC8_harm2223);
  //
  // TGraphErrors *fMC8_harm2233 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2233"));
  // fMC8_harm2233->SetMarkerStyle(kOpenCircle);
  // fMC8_harm2233->SetMarkerColor(kBlue);
  // fMC8_harm2233->SetMarkerSize(1.1);
  // fMC8_harm2233->SetMarkerColor(kBlue);
  // mg->Add(fMC8_harm2233);
  //
  // TGraphErrors *fMC8_harm2333 = new TGraphErrors((TH1D*) fileIn->Get("Refs_fMC8_harm2333"));
  // fMC8_harm2333->SetMarkerStyle(kOpenDiamond);
  // fMC8_harm2333->SetMarkerColor(kGreen+2);
  // fMC8_harm2333->SetMarkerSize(1.6);
  // fMC8_harm2333->SetLineColor(kGreen+2);
  // mg->Add(fMC8_harm2333);

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
  leg->SetHeader("ALICE experiment, LHC15o, additional pile up rejection, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(fMC8_harm2223_p,"nFMC (v_{2}^{6},v_{3}^{2}) - Plus","p");
  leg->AddEntry(fMC8_harm2223_m,"nFMC (v_{2}^{6},v_{3}^{2}) - Minus","p");
  leg->AddEntry(fMC8_harm2233_p,"nFMC (v_{2}^{4},v_{3}^{4}) - Plus","p");
  leg->AddEntry(fMC8_harm2233_m,"nFMC (v_{2}^{4},v_{3}^{4}) - Minus","p");
  leg->AddEntry(fMC8_harm2333_p,"nFMC (v_{2}^{2},v_{3}^{6}) - Plus","p");
  leg->AddEntry(fMC8_harm2333_m,"nFMC (v_{2}^{2},v_{3}^{6}) - Minus","p");
  leg->Draw("same");

  can->SaveAs("FMC8_comp_PM.pdf");

  TCanvas* can2[6];
  TH1D* ratio;
  TH1D* denom;
  TString name[6] = {"harm2223", "harm2233", "harm2333", "harm2223", "harm2233", "harm2333"};
  for(int i = 0; i < 6; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    if(i<3) ratio = (TH1D*) fileIn_p->Get("Refs_fMC8_" + name[i]);
    else ratio = (TH1D*) fileIn_m->Get("Refs_fMC8_" + name[i]);
    denom = (TH1D*) fileIn->Get("Refs_fMC8_" + name[i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality; nFMC (v_{2}^{k},v_{3}^{l}) / nFMC (v_{2}^{k},v_{3}^{l})");
    if(i < 3) ratio->SetTitle("pos " + name[i]);
    else ratio->SetTitle("neg " + name[i]);
    if(i == 0 || i == 3) ratio->GetYaxis()->SetRangeUser(0.,2.);
    else ratio->GetYaxis()->SetRangeUser(-1.,3.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0, 100);
    ratio->Fit(fu1);
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    if(i<3) can2[i]->SaveAs("ratio_p_" + name[i]+ ".pdf");
    else can2[i]->SaveAs("ratio_m_" + name[i]+ ".pdf");
  }




}
