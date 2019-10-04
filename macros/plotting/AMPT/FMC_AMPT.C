TString namesZM[9] = {"Refs_FMC4_harm23", "Refs_FMC4_harm24", "Refs_FMC4_harm34", "Refs_fMC6_harm223", "Refs_fMC6_harm233", "Refs_fMC6_harm234", "Refs_fMC8_harm2223", "Refs_fMC8_harm2233", "Refs_fMC8_harm2333"};
TString namesYZ[9] = {"his_nFMC3232", "his_nFMC4242", "his_nFMC4343", "his_nFMC322322", "his_nFMC332332", "his_nFMC432432", "his_nFMC32223222", "his_nFMC33223322", "his_nFMC33323332"};
TString names[9] = {"FMC3232", "FMC4242", "FMC4343", "FMC322322", "FMC332332", "FMC432432", "FMC32223222", "FMC33223322", "FMC33323332"};

void FMC_AMPT()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/MC/AMPT_train2388/processUniFlow/Processed_NOGAP.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}

  // int id = 0; //FMC4
  // int id = 3; //FMC6
  int id = 6; //FMC8


  TGraphErrors *FMC4_harm23 = new TGraphErrors((TH1D*) fileIn->Get(namesZM[id]));
  FMC4_harm23->SetMarkerStyle(kOpenSquare);
  FMC4_harm23->SetMarkerColor(kRed);
  FMC4_harm23->SetMarkerSize(1.0);
  FMC4_harm23->SetLineColor(kRed);
  mg->Add(FMC4_harm23);

  TGraphErrors *FMC4_harm23_YZ = new TGraphErrors((TH1D*) fileYZ->Get(namesYZ[id]));
  FMC4_harm23_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm23_YZ->SetMarkerColor(kRed);
  FMC4_harm23_YZ->SetMarkerSize(1.);
  FMC4_harm23_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm23_YZ);

  TGraphErrors *FMC4_harm24 = new TGraphErrors((TH1D*) fileIn->Get(namesZM[id+1]));
  FMC4_harm24->SetMarkerStyle(kOpenCircle);
  FMC4_harm24->SetMarkerColor(kBlue+1);
  FMC4_harm24->SetMarkerSize(1.0);
  FMC4_harm24->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24);

  TGraphErrors *FMC4_harm24_YZ = new TGraphErrors((TH1D*) fileYZ->Get(namesYZ[id+1]));
  FMC4_harm24_YZ->SetMarkerStyle(kFullCircle);
  FMC4_harm24_YZ->SetMarkerColor(kBlue+1);
  FMC4_harm24_YZ->SetMarkerSize(1.);
  FMC4_harm24_YZ->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24_YZ);

  TGraphErrors *FMC4_harm34 = new TGraphErrors((TH1D*) fileIn->Get(namesZM[id+2]));
  FMC4_harm34->SetMarkerStyle(kOpenDiamond);
  FMC4_harm34->SetMarkerColor(kGreen+2);
  FMC4_harm34->SetMarkerSize(1.0);
  FMC4_harm34->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34);

  TGraphErrors *FMC4_harm34_YZ = new TGraphErrors((TH1D*) fileYZ->Get(namesYZ[id+2]));
  FMC4_harm34_YZ->SetMarkerStyle(kFullDiamond);
  FMC4_harm34_YZ->SetMarkerColor(kGreen+2);
  FMC4_harm34_YZ->SetMarkerSize(1.);
  FMC4_harm34_YZ->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34_YZ);


  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("nFMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  TLegend* leg;
  mg->Draw("ap");
  double min = 0;
  double max = 0;
  if(id == 0){
    min = -0.3;
    max = 1.1;
    leg = new TLegend(0.12,0.72,0.52,0.88);
  }
  else if(id == 3){
    min = -0.9;
    max = 0.2;
    leg = new TLegend(0.12,0.22,0.52,0.38);
  }
  else{
    min = -0.35;
    max = 0.2;
    leg = new TLegend(0.12,0.72,0.52,0.88);
  }
  mg->SetMinimum(min);
  mg->SetMaximum(max);

  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  if(id == 0) {
    leg->AddEntry(FMC4_harm23,"nSC (2,3) - AMPT","p");
    leg->AddEntry(FMC4_harm23_YZ,"nSC (2,3) - YZ","p");
    leg->AddEntry(FMC4_harm24,"nSC (2,4) - AMPT","p");
    leg->AddEntry(FMC4_harm24_YZ,"nSC (2,4) - YZ","p");
    leg->AddEntry(FMC4_harm34,"nSC (3,4) - AMPT","p");
    leg->AddEntry(FMC4_harm34_YZ,"nSC (3,4) - YZ","p");
  }
  else if(id == 3) {
    leg->AddEntry(FMC4_harm23,"nFMC(v_{2}^{4},v_{3}^{2}) - AMPT","p");
    leg->AddEntry(FMC4_harm23_YZ,"nFMC(v_{2}^{4},v_{3}^{2}) - YZ","p");
    leg->AddEntry(FMC4_harm24,"nFMC(v_{2}^{2},v_{3}^{4}) - AMPT","p");
    leg->AddEntry(FMC4_harm24_YZ,"nFMC(v_{2}^{2},v_{3}^{4}) - YZ","p");
    leg->AddEntry(FMC4_harm34,"nFMC(v_{2}^{2},v_{3}^{2},v_{4}^{2}) - AMPT","p");
    leg->AddEntry(FMC4_harm34_YZ,"nFMC(v_{2}^{2},v_{3}^{2},v_{4}^{2}) - YZ","p");
  }
  else {
    leg->AddEntry(FMC4_harm23,"nFMC(v_{2}^{6},v_{3}^{2}) - AMPT","p");
    leg->AddEntry(FMC4_harm23_YZ,"nFMC(v_{2}^{6},v_{3}^{2}) - YZ","p");
    leg->AddEntry(FMC4_harm24,"nFMC(v_{2}^{4},v_{3}^{4}) - AMPT","p");
    leg->AddEntry(FMC4_harm24_YZ,"nFMC(v_{2}^{4},v_{3}^{4}) - YZ","p");
    leg->AddEntry(FMC4_harm34,"nFMC(v_{2}^{2},v_{3}^{6}) - AMPT","p");
    leg->AddEntry(FMC4_harm34_YZ,"nFMC(v_{2}^{2},v_{3}^{6}) - YZ","p");
  }
  leg->Draw("same");

  if(id == 0) can->SaveAs("SC_AMPT.pdf");
  else if(id == 3) can->SaveAs("FMC6_AMPT.pdf");
  else can->SaveAs("FMC8_AMPT.pdf");

  TCanvas* can2[3];
  TH1D* ratio;
  TH1D* denom;
  for(int i = 0; i < 3; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    ratio = (TH1D*) fileIn->Get(namesZM[id+i]);
    denom = (TH1D*) fileYZ->Get(namesYZ[id+i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality; nFMC: AMPT / data");
    ratio->SetTitle(names[id+i]);
    // ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetRangeUser(0.,60.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0, 60);
    ratio->Fit(fu1);
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    can2[i]->SaveAs("ratio_" + names[id+i]+ ".pdf");
  }



}
