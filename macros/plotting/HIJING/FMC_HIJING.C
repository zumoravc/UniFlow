TString namesZM[9] = {"Refs_FMC4_harm23", "Refs_FMC4_harm24", "Refs_FMC4_harm34", "Refs_fMC6_harm223", "Refs_fMC6_harm233", "Refs_fMC6_harm234", "Refs_fMC8_harm2223", "Refs_fMC8_harm2233", "Refs_fMC8_harm2333"};
TString namesYZ[9] = {"his_FMC3232", "his_FMC4242", "his_FMC4343", "his_FMC322322", "his_FMC332332", "his_FMC432432", "his_FMC32223222", "his_FMC33223322", "his_FMC33323332"};
TString names[9] = {"FMC3232", "FMC4242", "FMC4343", "FMC322322", "FMC332332", "FMC432432", "FMC32223222", "FMC33223322", "FMC33323332"};

void FMC_HIJING()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/MC/HIJING_train_2381/processUniFlow/Processed.root","READ");
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
  mg->GetYaxis()->SetTitle("FMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-1.3E-9);
  mg->SetMaximum(0.5E-9);

  TLegend* leg = new TLegend(0.12,0.72,0.52,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  if(id == 0) {
    leg->AddEntry(FMC4_harm23,"SC (2,3) - HIJING","p");
    leg->AddEntry(FMC4_harm23_YZ,"SC (2,3) - YZ","p");
    leg->AddEntry(FMC4_harm24,"SC (2,4) - HIJING","p");
    leg->AddEntry(FMC4_harm24_YZ,"SC (2,4) - YZ","p");
    leg->AddEntry(FMC4_harm34,"SC (3,4) - HIJING","p");
    leg->AddEntry(FMC4_harm34_YZ,"SC (3,4) - YZ","p");
  }
  else if(id == 3) {
    leg->AddEntry(FMC4_harm23,"FMC(v_{2}^{4},v_{3}^{2}) - HIJING","p");
    leg->AddEntry(FMC4_harm23_YZ,"FMC(v_{2}^{4},v_{3}^{2}) - YZ","p");
    leg->AddEntry(FMC4_harm24,"FMC(v_{2}^{2},v_{3}^{4}) - HIJING","p");
    leg->AddEntry(FMC4_harm24_YZ,"FMC(v_{2}^{2},v_{3}^{4}) - YZ","p");
    leg->AddEntry(FMC4_harm34,"FMC(v_{2}^{2},v_{3}^{2},v_{4}^{2}) - HIJING","p");
    leg->AddEntry(FMC4_harm34_YZ,"FMC(v_{2}^{2},v_{3}^{2},v_{4}^{2}) - YZ","p");
  }
  else {
    leg->AddEntry(FMC4_harm23,"FMC(v_{2}^{6},v_{3}^{2}) - HIJING","p");
    leg->AddEntry(FMC4_harm23_YZ,"FMC(v_{2}^{6},v_{3}^{2}) - YZ","p");
    leg->AddEntry(FMC4_harm24,"FMC(v_{2}^{4},v_{3}^{4}) - HIJING","p");
    leg->AddEntry(FMC4_harm24_YZ,"FMC(v_{2}^{4},v_{3}^{4}) - YZ","p");
    leg->AddEntry(FMC4_harm34,"FMC(v_{2}^{2},v_{3}^{6}) - HIJING","p");
    leg->AddEntry(FMC4_harm34_YZ,"FMC(v_{2}^{2},v_{3}^{6}) - YZ","p");
  }

  leg->Draw("same");

  if(id == 0) can->SaveAs("SC_HIJING.pdf");
  else if(id == 3) can->SaveAs("FMC6_HIJING.pdf");
  else can->SaveAs("FMC8_HIJING.pdf");

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
    ratio->SetTitle(";Centrality;HIJING / data");
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
