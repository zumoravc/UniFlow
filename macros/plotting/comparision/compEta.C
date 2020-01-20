void compEta()
{
  TMultiGraph *mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.12,0.72,0.52,0.88);

  TFile* filePos = TFile::Open("/home/alidock/output/LHC15o/train_7387/processUniFlow/Processed_pos.root","READ");
  if(!filePos) {printf("File pos not opened! \n"); return;}

  TFile* fileNeg = TFile::Open("/home/alidock/output/LHC15o/train_7387/processUniFlow/Processed_neg.root","READ");
  if(!fileNeg) {printf("File neg not opened! \n"); return;}

  TGraphErrors *gr[3][2];

  Color_t col[3] = {kRed, kBlue, kGreen+2};

  for(Int_t i(0); i < 3; i++){
    gr[i][0] = new TGraphErrors((TH1D*) filePos->Get(Form("Refs_hFlow2_harm%d_gap00",i+2)));
    gr[i][1] = new TGraphErrors((TH1D*) fileNeg->Get(Form("Refs_hFlow2_harm%d_gap00",i+2)));

    gr[i][0]->SetMarkerStyle(kOpenSquare);
    gr[i][1]->SetMarkerStyle(kOpenCircle);

    for(Int_t j(0); j < 2; j++){
      gr[i][j]->SetMarkerColor(col[i]);
      gr[i][j]->SetLineColor(col[i]);

      mg->Add(gr[i][j]);
    }
    leg->AddEntry(gr[i][0], Form("v_{%d}{2}, #eta > 0.0",i+2), "p");
    leg->AddEntry(gr[i][1], Form("v_{%d}{2}, #eta < 0.0",i+2), "p");
  }



  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("v_{n}{2}");
  mg->GetXaxis()->SetRangeUser(0, 80);


  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(0.0);
  mg->SetMaximum(0.15);

  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, RFPs, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->Draw("same");

  can->SaveAs("comp_vn2.pdf");

  TCanvas* can2[3];
  TH1D* ratio;
  TH1D* denom;
  for(Int_t i(0); i < 3; i++){
    can2[i] = new TCanvas(Form("can2_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    ratio = (TH1D*) filePos->Get(Form("Refs_hFlow2_harm%d_gap00",i+2));
    denom = (TH1D*) fileNeg->Get(Form("Refs_hFlow2_harm%d_gap00",i+2));
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(Form(";Centrality; v_{%d}{2}: pos / neg #eta",i+2));
    ratio->GetXaxis()->SetRangeUser(0.,80.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0, 80);
    ratio->Fit(fu1,"R");
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    can2[i]->SaveAs(Form("ratio_v%d2.pdf",i+2));
  }


}
