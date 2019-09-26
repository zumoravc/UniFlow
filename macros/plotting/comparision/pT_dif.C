void pT_dif()
{
  TMultiGraph *mg[9];
  TCanvas* can[9];

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/grid_fullSt/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TGraphErrors* v22[9];
  TGraphErrors* v32[9];
  TGraphErrors* v42[9];
  for(int i = 0; i < 9; i++)
  {
    v22[i] = new TGraphErrors((TH1D*) fileIn->Get(Form("Charged_hFlow2_harm2_gap10_cent%d",i)));
    v32[i] = new TGraphErrors((TH1D*) fileIn->Get(Form("Charged_hFlow2_harm3_gap10_cent%d",i)));
    v42[i] = new TGraphErrors((TH1D*) fileIn->Get(Form("Charged_hFlow2_harm4_gap10_cent%d",i)));

    v22[i]->SetMarkerStyle(kOpenSquare);
    v22[i]->SetMarkerColor(kRed);
    v22[i]->SetMarkerSize(1.);
    v22[i]->SetLineColor(kRed);

    v32[i]->SetMarkerStyle(kOpenCircle);
    v32[i]->SetMarkerColor(kBlue);
    v32[i]->SetMarkerSize(1.1);
    v32[i]->SetMarkerColor(kBlue);

    v42[i]->SetMarkerStyle(kOpenDiamond);
    v42[i]->SetMarkerColor(kGreen+2);
    v42[i]->SetMarkerSize(1.6);
    v42[i]->SetLineColor(kGreen+2);

    mg[i] = new TMultiGraph();
    mg[i]->Add(v22[i]);
    mg[i]->Add(v32[i]);
    mg[i]->Add(v42[i]);

    mg[i]->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
    mg[i]->GetYaxis()->SetTitle("v_{n}{2}");
    mg[i]->GetXaxis()->SetRangeUser(0, 5);



  }

  TFile* filePUB[24];
  TGraphErrors* pub[24];
  TDirectory* dir;
  //0-5, 5-10, 30-40
  int num[24] = {15, 17, 19, 23, 25, 27, 33, 35, 37, 43, 45, 47, 53, 55, 57, 63, 65, 67, 73, 75, 77, 82, 84, 86};
  for(int i = 0; i < 24; i++)
  {
  filePUB[i] = TFile::Open("/home/alidock/ana/output/published/published.root","READ");
  if(!filePUB[i]) {printf("File YZ not opened! \n"); return;}

  dir = gFile->GetDirectory(Form("Table %d",num[i]));
  // dir->GetObject("Graph1D_y1", pub[i]);

  pub[i] = new TGraphErrors((TH1D*) dir->Get("Hist1D_y1"));

  }

  Color_t col[3] = {kRed, kBlue, kGreen+2};
  for(int i = 0; i < 3; i++)
  {
    if(i == 0){
      for(int j = 0; j < 8; j++){
        pub[j*3+i]->SetMarkerStyle(kFullSquare);
        pub[j*3+i]->SetMarkerColor(col[i]);
        pub[j*3+i]->SetLineColor(col[i]);
      }
    }
    if(i == 1){
      for(int j = 0; j < 8; j++){
        pub[j*3+i]->SetMarkerStyle(kFullCircle);
        pub[j*3+i]->SetMarkerColor(col[i]);
        pub[j*3+i]->SetLineColor(col[i]);
      }
    }
    if(i == 2){
      for(int j = 0; j < 8; j++)
      {
        pub[j*3+i]->SetMarkerStyle(kFullDiamond);
        pub[j*3+i]->SetMarkerColor(col[i]);
        pub[j*3+i]->SetLineColor(col[i]);
      }
    }
  }


  TLegend* leg[9];


  for(int i = 0; i < 9; i++)
  {
    if(i<8) {
      mg[i]->Add(pub[3*i]);
      mg[i]->Add(pub[3*i+1]);
      mg[i]->Add(pub[3*i+1]);
    }
    can[i] = new TCanvas(Form("can%d",i), "can", 600, 400);
    mg[i]->Draw("ap");
    mg[i]->SetMinimum(0.0);
    if(i < 4) mg[i]->SetMaximum(0.2);
    else mg[i]->SetMaximum(0.5);

   leg[i]  = new TLegend(0.12,0.65,0.52,0.88);
   leg[i]->SetBorderSize(0);
   leg[i]->SetFillColorAlpha(0.0,0.0);
   gStyle->SetLegendTextSize(0.026);
   leg[i]->SetHeader("ALICE experiment, LHC15o (sample), additional pile up rejection, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
   leg[i]->SetNColumns(2);

    leg[i]->AddEntry(v22[i], "v_{2}{2}(|#Delta#eta|>1.0)", "p");
    if(i < 8) leg[i]->AddEntry(pub[i*3], "v_{2}{2}(|#Delta#eta|>1.0) - published", "p");
    leg[i]->AddEntry(v32[i], "v_{3}{2}(|#Delta#eta|>1.0)", "p");
    if(i < 8) leg[i]->AddEntry(pub[i*3+1], "v_{3}{2}(|#Delta#eta|>1.0)  - published", "p");
    leg[i]->AddEntry(v42[i], "v_{4}{2}(|#Delta#eta|>1.0)", "p");
    if(i < 8) leg[i]->AddEntry(pub[i*3+2], "v_{4}{2}(|#Delta#eta|>1.0) - published ", "p");
    leg[i]->Draw("same");

    can[i]->SaveAs(Form("v234_h_cent%d.pdf",i));
  }



  // can->SaveAs("FMC_comp_new.pdf");

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
