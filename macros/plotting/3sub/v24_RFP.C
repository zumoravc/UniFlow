char sides[4] = "LMR";
Color_t color[] = {kRed, kBlue, kGreen+3};
Marker_t marker[] = {kOpenCircle, kOpenSquare, kOpenDiamond};

void v24_RFP(){

  TFile* fileIn = TFile::Open("/home/alidock/output/LHC15o/train_7456/processUniFlow/Processed_binningVP.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileV = TFile::Open("/home/alidock/output/published/Vojta/Processed.root","READ");
  if(!fileV) {printf("File VP not opened! \n"); return;}

  TMultiGraph *mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.12,0.75,0.72,0.88);

  TGraphErrors *zm = new TGraphErrors((TH1D*) fileIn->Get("Refs_hFlow4_harm2_gap(08,08)_3sub"));
  zm->SetMarkerStyle(kOpenSquare);
  zm->SetMarkerColor(kRed);
  zm->SetMarkerSize(1.);
  zm->SetLineColor(kRed);
  mg->Add(zm);

  TGraphErrors *combi[3] = {nullptr};
  for(Int_t i(0); i < 3; i++){
    combi[i] = new TGraphErrors((TH1D*) fileIn->Get(Form("Refs_hFlow4_harm2_gap(08,08)_3sub_two_%c",sides[i])));
    combi[i]->SetMarkerStyle(marker[i]);
    combi[i]->SetMarkerColor(color[i]);
    combi[i]->SetLineColor(color[i]);
    // mg->Add(combi[i]);
    // leg->AddEntry(combi[i], Form("v_{2}{4}_{3 sub}, 2p. in %c",sides[i]), "p");
  }

  TGraphErrors *vp = new TGraphErrors((TH1D*) fileV->Get("Refs_hFlow4_harm2_gap00"));
  vp->SetMarkerStyle(kOpenCircle);
  vp->SetMarkerColor(kBlue);
  vp->SetMarkerSize(1.);
  vp->SetLineColor(kBlue);
  mg->Add(vp);

  mg->GetXaxis()->SetTitle("Centtrality class (V0M)");
  mg->GetYaxis()->SetTitle("v_{2}{4}");
  // mg->GetXaxis()->SetRangeUser(0, 60);

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(0.03);
  mg->SetMaximum(0.11);

  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  TString legString = "ALICE experiment, PbPb @ 5.02TeV, LHC15o";

  leg->SetHeader(legString.Data());
  gStyle->SetLegendTextSize(0.025);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->AddEntry(zm,"v_{2}{4}_{3sub} - ZM","p");
  leg->AddEntry(vp,"v_{2}{4, |#Delta#eta|>0.0} - VP","p");
  leg->Draw("same");

  can->SaveAs("v24_RFP_comp.pdf");

  TCanvas* cRat = new TCanvas("cRat", "cRat", 600, 400);
  TH1D* hisZ = (TH1D*) fileIn->Get("Refs_hFlow4_harm2_gap(08,08)_3sub");
  TH1D* hisV = (TH1D*) fileV->Get("Refs_hFlow4_harm2_gap00");

  TLegend* leg2 = new TLegend(0.12,0.75,0.55,0.88);
  TF1 *fu1 = new TF1("fu1", "pol0", 0., 60.);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(1);

  hisV->Divide(hisZ);
  hisV->SetMarkerStyle(kFullSquare);
  hisV->SetMarkerColor(kBlack);

  hisV->Draw("h");
  hisV->GetYaxis()->SetRangeUser(0.98,1.01);
  hisV->GetXaxis()->SetTitle("Centtrality class (V0M)");
  hisV->GetYaxis()->SetTitle("VP / ZM ");

  leg2->SetFillColorAlpha(0.0,0.0);
  leg2->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.02);
  leg2->SetHeader(legString.Data());
  leg2->Draw("same");

  hisV->Fit(fu1,"R");
  fu1->SetLineColor(kRed);
  fu1->Draw("same");

  cRat->SaveAs("v24_RFP_ratio.pdf");

}
