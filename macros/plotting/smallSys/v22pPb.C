void v22pPb(){
  TFile* fileIn = TFile::Open("/home/alidock/ana/output/pPb_LHC16q/HADD/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TMultiGraph *mg = new TMultiGraph();

  TGraphErrors *nogap = new TGraphErrors((TH1D*) fileIn->Get("Kaon_hFlow2_harm2_gap-10_cent0"));
  nogap->SetMarkerStyle(kOpenSquare);
  nogap->SetMarkerColor(kRed);
  nogap->SetMarkerSize(1.);
  nogap->SetLineColor(kRed);
  mg->Add(nogap);

  TGraphErrors *gap00 = new TGraphErrors((TH1D*) fileIn->Get("Kaon_hFlow2_harm2_gap00_cent0"));
  gap00->SetMarkerStyle(kOpenCircle);
  gap00->SetMarkerColor(kGreen+3);
  gap00->SetMarkerSize(1.);
  gap00->SetLineColor(kGreen+3);
  mg->Add(gap00);

  TGraphErrors *gap08 = new TGraphErrors((TH1D*) fileIn->Get("Kaon_hFlow2_harm2_gap08_cent0"));
  gap08->SetMarkerStyle(kOpenDiamond);
  gap08->SetMarkerColor(kBlue);
  gap08->SetMarkerSize(1.);
  gap08->SetLineColor(kBlue);
  mg->Add(gap08);

  TGraphErrors *gap10 = new TGraphErrors((TH1D*) fileIn->Get("Kaon_hFlow2_harm2_gap10_cent0"));
  gap10->SetMarkerStyle(kOpenStar);
  gap10->SetMarkerColor(kMagenta+2);
  gap10->SetMarkerSize(1.);
  gap10->SetLineColor(kMagenta+2);
  mg->Add(gap10);

  mg->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
  mg->GetYaxis()->SetTitle("v_{2}{2}");
  // mg->GetXaxis()->SetRangeUser(0, 60);

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  // mg->SetMinimum(0.0);
  // mg->SetMaximum(0.3);

  TLegend* leg = new TLegend(0.12,0.65,0.52,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, pPb @ 5.02TeV, LHC16q, K^{#pm}, |#eta| < 0.8, V0A cent: 0-10%");
  // leg->SetHeader("ALICE experiment, pPb @ 5.02TeV, LHC16q, p/#bar{p}, |#eta| < 0.8, V0A cent: 0-1%");
  gStyle->SetLegendTextSize(0.025);
  leg->SetFillStyle(0);
  // leg->SetNColumns(2);
  // leg->AddEntry((TObject*) 0, "0.2 < p_{T}(RFPs) < 5 GeV/c, RFPs: |#eta| < 0.8", "");
  // leg->AddEntry((TObject*) 0, "RFPs: |#eta| < 0.8", "");
  leg->AddEntry(nogap,"v_{2}{2}","p");
  leg->AddEntry(gap00,"v_{2}{2,|#Delta#eta|>0.0}","p");
  leg->AddEntry(gap08,"v_{2}{2,|#Delta#eta|>0.8}","p");
  leg->AddEntry(gap10,"v_{2}{2,|#Delta#eta|>1.0}","p");
  leg->Draw("same");

  can->SaveAs("v22_pPb_K.pdf");

}
