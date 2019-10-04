TString namesZM[10] = {"Refs_FMC4_harm23", "Refs_FMC4_harm24", "Refs_FMC4_harm34", "Refs_fMC6_harm223", "Refs_fMC6_harm233", "Refs_fMC6_harm234", "Refs_fMC6_harm345", "Refs_fMC8_harm2223", "Refs_fMC8_harm2233", "Refs_fMC8_harm2333"};
TString namesGap[10] = {"Refs_FMC4_harm23_2sub(0)", "Refs_FMC4_harm24_2sub(0)", "Refs_FMC4_harm34_2sub(0)", "Refs_fMC6_harm223_2sub(0)", "Refs_fMC6_harm233_2sub(0)", "Refs_fMC6_harm234_2sub(0)", "Refs_fMC6_harm345_2sub(0)", "Refs_fMC8_harm2223_2sub(0)", "Refs_fMC8_harm2233_2sub(0)", "Refs_fMC8_harm2333_2sub(0)"};
TString names[10] = {"FMC3232", "FMC4242", "FMC4343", "FMC322322", "FMC332332", "FMC432432", "FMC543543", "FMC32223222", "FMC33223322", "FMC33323332"};
TString namesYZ[10] = {"his_nFMC3232", "his_nFMC4242", "his_nFMC4343", "his_nFMC322322", "his_nFMC332332", "his_nFMC432432", "his_nFMC543543", "his_nFMC32223222", "his_nFMC33223322", "his_nFMC33323332"};


void FMC_comp()
{
  TMultiGraph *mg = new TMultiGraph();

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/train_7387/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileNoGap = TFile::Open("/home/alidock/ana/output/LHC15o/train_7387/processUniFlow/Processed_NOGAP.root","READ");
  if(!fileNoGap) {printf("File NoGap not opened! \n"); return;}

  TFile* fileYZ = TFile::Open("/home/alidock/ana/output/YZ_results/Output_all.root","READ");
  if(!fileYZ) {printf("File YZ not opened! \n"); return;}

  // int id = 0; //FMC4
  // int id = 3; //FMC6
  int id = 7; //FMC8

  int set = 0; //gap vs. no gap
  // int set = 1; //YZ
  if(set == 1){
    fileIn = fileYZ;
    for(int i = 0; i < 10; i ++) namesGap[i] = namesYZ[i];
  }

  TGraphErrors *FMC4_harm23 = new TGraphErrors((TH1D*) fileNoGap->Get(namesZM[id]));
  FMC4_harm23->SetMarkerStyle(kOpenSquare);
  FMC4_harm23->SetMarkerColor(kRed);
  FMC4_harm23->SetMarkerSize(1.0);
  FMC4_harm23->SetLineColor(kRed);
  mg->Add(FMC4_harm23);

  TGraphErrors *FMC4_harm23_YZ = new TGraphErrors((TH1D*) fileIn->Get(namesGap[id]));
  FMC4_harm23_YZ->SetMarkerStyle(kFullSquare);
  FMC4_harm23_YZ->SetMarkerColor(kRed);
  FMC4_harm23_YZ->SetMarkerSize(1.);
  FMC4_harm23_YZ->SetLineColor(kRed);
  mg->Add(FMC4_harm23_YZ);

  TGraphErrors *FMC4_harm24 = new TGraphErrors((TH1D*) fileNoGap->Get(namesZM[id+1]));
  FMC4_harm24->SetMarkerStyle(kOpenCircle);
  FMC4_harm24->SetMarkerColor(kBlue+1);
  FMC4_harm24->SetMarkerSize(1.0);
  FMC4_harm24->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24);

  TGraphErrors *FMC4_harm24_YZ = new TGraphErrors((TH1D*) fileIn->Get(namesGap[id+1]));
  FMC4_harm24_YZ->SetMarkerStyle(kFullCircle);
  FMC4_harm24_YZ->SetMarkerColor(kBlue+1);
  FMC4_harm24_YZ->SetMarkerSize(1.);
  FMC4_harm24_YZ->SetLineColor(kBlue+1);
  mg->Add(FMC4_harm24_YZ);

  TGraphErrors *FMC4_harm34 = new TGraphErrors((TH1D*) fileNoGap->Get(namesZM[id+2]));
  FMC4_harm34->SetMarkerStyle(kOpenDiamond);
  FMC4_harm34->SetMarkerColor(kGreen+2);
  FMC4_harm34->SetMarkerSize(1.0);
  FMC4_harm34->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34);

  TGraphErrors *FMC4_harm34_YZ = new TGraphErrors((TH1D*) fileIn->Get(namesGap[id+2]));
  FMC4_harm34_YZ->SetMarkerStyle(kFullDiamond);
  FMC4_harm34_YZ->SetMarkerColor(kGreen+2);
  FMC4_harm34_YZ->SetMarkerSize(1.);
  FMC4_harm34_YZ->SetLineColor(kGreen+2);
  mg->Add(FMC4_harm34_YZ);

  TGraphErrors *FMC4_harm345 = nullptr;
  TGraphErrors *FMC4_harm345_YZ = nullptr;
  // if(id == 3)
  // {
  //   FMC4_harm345 = new TGraphErrors((TH1D*) fileIn->Get(namesZM[id+3]));
  //   FMC4_harm345->SetMarkerStyle(kOpenStar);
  //   FMC4_harm345->SetMarkerColor(kMagenta+2);
  //   FMC4_harm345->SetMarkerSize(1.0);
  //   FMC4_harm345->SetLineColor(kMagenta+2);
  //   mg->Add(FMC4_harm345);
  //
  //   FMC4_harm345_YZ = new TGraphErrors((TH1D*) fileNoGap->Get(namesGap[id+3]));
  //   FMC4_harm345_YZ->SetMarkerStyle(kFullStar);
  //   FMC4_harm345_YZ->SetMarkerColor(kMagenta+2);
  //   FMC4_harm345_YZ->SetMarkerSize(1.);
  //   FMC4_harm345_YZ->SetLineColor(kMagenta+2);
  //   mg->Add(FMC4_harm345_YZ);
  // }


  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("nFMC (v_{2}^{k},v_{3}^{l})");
  mg->GetXaxis()->SetRangeUser(0, 60);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-0.32);
  mg->SetMaximum(.13);

  TLegend* leg = new TLegend(0.32,0.72,0.72,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, LHC15o, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  if(id == 0) {
    if(set == 0){
      leg->AddEntry(FMC4_harm23,"nSC (2,3) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm23_YZ,"nSC (2,3)","p");
      leg->AddEntry(FMC4_harm24,"nSC (2,4) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm24_YZ,"nSC (2,4)","p");
      leg->AddEntry(FMC4_harm34,"nSC (3,4) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm34_YZ,"nSC (3,4)","p");
    }
    else {
      leg->AddEntry(FMC4_harm23,"nSC (2,3)","p");
      leg->AddEntry(FMC4_harm23_YZ,"nSC (2,3) - YZ","p");
      leg->AddEntry(FMC4_harm24,"nSC (2,4)","p");
      leg->AddEntry(FMC4_harm24_YZ,"nSC (2,4) - YZ","p");
      leg->AddEntry(FMC4_harm34,"nSC (3,4)","p");
      leg->AddEntry(FMC4_harm34_YZ,"nSC (3,4) - YZ","p");
    }
  }
  else if(id == 3) {
    if(set == 0) {
      leg->AddEntry(FMC4_harm23,"FMC(v_{2}^{4},v_{3}^{2}) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm23_YZ,"FMC(v_{2}^{4},v_{3}^{2})","p");
      leg->AddEntry(FMC4_harm24,"FMC(v_{2}^{2},v_{3}^{4})(|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm24_YZ,"FMC(v_{2}^{2},v_{3}^{4})","p");
      leg->AddEntry(FMC4_harm34,"FMC(v_{2}^{2},v_{3}^{2},v_{4}^{2}) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm34_YZ,"FMC(v_{2}^{2},v_{3}^{2},v_{4}^{2})","p");
      leg->AddEntry(FMC4_harm345,"FMC(v_{3}^{2},v_{4}^{2},v_{5}^{2}) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm345_YZ,"FMC(v_{3}^{2},v_{4}^{2},v_{5}^{2})","p");
    }
    else {
      leg->AddEntry(FMC4_harm23,"nFMC(v_{2}^{4},v_{3}^{2})","p");
      leg->AddEntry(FMC4_harm23_YZ,"nFMC(v_{2}^{4},v_{3}^{2}) - YZ","p");
      leg->AddEntry(FMC4_harm24,"nFMC(v_{2}^{2},v_{3}^{4})","p");
      leg->AddEntry(FMC4_harm24_YZ,"nFMC(v_{2}^{2},v_{3}^{4}) - YZ","p");
      leg->AddEntry(FMC4_harm34,"nFMC(v_{2}^{2},v_{3}^{2},v_{4}^{2})","p");
      leg->AddEntry(FMC4_harm34_YZ,"nFMC(v_{2}^{2},v_{3}^{2},v_{4}^{2}) - YZ","p");

    }
  }
  else {
    if(set == 0) {
      leg->AddEntry(FMC4_harm23,"FMC(v_{2}^{6},v_{3}^{2}) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm23_YZ,"FMC(v_{2}^{6},v_{3}^{2})","p");
      leg->AddEntry(FMC4_harm24,"FMC(v_{2}^{4},v_{3}^{4}) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm24_YZ,"FMC(v_{2}^{4},v_{3}^{4})","p");
      leg->AddEntry(FMC4_harm34,"FMC(v_{2}^{2},v_{3}^{6}) (|#Delta#eta|>0.0)","p");
      leg->AddEntry(FMC4_harm34_YZ,"FMC(v_{2}^{2},v_{3}^{6}) ","p");
    }
    else {
      leg->AddEntry(FMC4_harm23,"nFMC(v_{2}^{6},v_{3}^{2})","p");
      leg->AddEntry(FMC4_harm23_YZ,"nFMC(v_{2}^{6},v_{3}^{2}) - YZ","p");
      leg->AddEntry(FMC4_harm24,"nFMC(v_{2}^{4},v_{3}^{4})","p");
      leg->AddEntry(FMC4_harm24_YZ,"nMC(v_{2}^{4},v_{3}^{4}) - YZ","p");
      leg->AddEntry(FMC4_harm34,"nFMC(v_{2}^{2},v_{3}^{6})","p");
      leg->AddEntry(FMC4_harm34_YZ,"nFMC(v_{2}^{2},v_{3}^{6}) - YZ","p");
    }
  }

  leg->Draw("same");

  if(id == 0) can->SaveAs("YZ/SC.pdf");
  else if(id == 3) can->SaveAs("YZ/FMC6.pdf");
  else can->SaveAs("YZ/FMC8.pdf");

  TCanvas* can2[3];
  TH1D* ratio;
  TH1D* denom;
  int max = 3;
  if(id == 3 && set == 0) max = 4;
  for(int i = 0; i < max; i++)
  {
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    ratio = (TH1D*) fileNoGap->Get(namesZM[id+i]);
    denom = (TH1D*) fileYZ->Get(namesYZ[id+i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";Centrality; ZM / YZ");
    ratio->SetTitle(names[id+i]);
    // ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetRangeUser(5.,60.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 5, 60);
    ratio->Fit(fu1);
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    can2[i]->SaveAs("YZ/ratio_" + names[id+i]+ ".pdf");
  }


}
