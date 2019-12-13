void v24_nch(Int_t part, Int_t centClass){

  TString specLong = "";
  TString spec = "";
  switch(part){
    case 0: spec = "h"; specLong = "Charged"; break;
    case 1: spec = "pi"; specLong = "Pion"; break;
    case 2: spec = "K"; specLong = "Kaon"; break;
    case 3: spec = "p"; specLong = "Proton"; break;
    default:
      printf("There is a problem... \n"); return;
  }


  TFile* fileIn = TFile::Open("/home/alidock/output/LHC15o/train_7456/processUniFlow/Processed_Nch1.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.12,0.8,0.75,0.88);

  // TGraphErrors* v24nogap = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap-10_cent%d",specLong.Data(),centClass)));
  // v24nogap->SetMarkerStyle(kOpenCircle);
  // v24nogap->SetMarkerColor(kGreen+2);
  // v24nogap->SetLineColor(kGreen+2);
  // mg->Add(v24nogap);
  // leg->AddEntry(v24nogap, "v_{2}{4}", "p");
  //
  // TGraphErrors* v24gap10 = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap10_cent%d",specLong.Data(),centClass)));
  // v24gap10->SetMarkerStyle(kOpenDiamond);
  // v24gap10->SetMarkerColor(kBlue);
  // v24gap10->SetLineColor(kBlue);
  // mg->Add(v24gap10);
  // leg->AddEntry(v24gap10, "v_{2}{4,|#Delta#eta| > 1.0}", "p");
  //
  // TGraphErrors* v24gap00 = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap00_cent%d",specLong.Data(),centClass)));
  // v24gap00->SetMarkerStyle(kOpenSquare);
  // v24gap00->SetMarkerColor(kRed);
  // v24gap00->SetLineColor(kRed);
  // mg->Add(v24gap00);
  // leg->AddEntry(v24gap00, "v_{2}{4,|#Delta#eta| > 0.0}", "p");

  // TGraphErrors* gap08 = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap08_cent%d",specLong.Data(),centClass)));
  // gap08->SetMarkerStyle(kOpenDiamond);
  // gap08->SetMarkerColor(kBlue);
  // gap08->SetLineColor(kBlue);
  // // mg->Add(gap08);
  // // leg->AddEntry(gap08, "v_{2}{4,|#Delta#eta| > 0.8}", "p");
  //
  TGraphErrors* allCombi = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub",specLong.Data(),centClass)));
  allCombi->SetMarkerStyle(kOpenSquare);
  allCombi->SetMarkerColor(kRed);
  allCombi->SetLineColor(kRed);
  mg->Add(allCombi);
  leg->AddEntry(allCombi, "v_{2}{4}_{3sub}", "p");


  mg->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
  mg->GetYaxis()->SetTitle("v_{2}{4}_{3 sub}");
  mg->GetXaxis()->SetRangeUser(0, 5);

  mg->Draw("ap");
  mg->SetMinimum(0.0);
  mg->SetMaximum(0.4);

  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  TString legString = "ALICE experiment, PbPb @ 5.02TeV, LHC15o, ";
  switch(part){
    case 0: legString += "h^{#pm}, |#eta| < 0.8, N_{ch}: "; break;
    case 1: legString += "#pi^{#pm}, |#eta| < 0.8, N_{ch}: "; break;
    case 2: legString += "K^{#pm}, |#eta| < 0.8, N_{ch}: "; break;
    case 3: legString += "p/#bar{p}, |#eta| < 0.8, N_{ch}: "; break;
  }

  // switch(centClass){
  //   case 0: legString += "0-25"; break;
  //   case 1: legString += "25-50"; break;
  //   case 2: legString += "50-75"; break;
  //   case 3: legString += "75-100"; break;
  //   case 4: legString += "100-200"; break;
  // }

  switch(centClass){
    case 0: legString += "50-100"; break;
    case 1: legString += "100-150"; break;
  }

  leg->SetHeader(legString.Data());
  gStyle->SetLegendTextSize(0.025);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->Draw("same");

  can->SaveAs(Form("ptdif/v24_PbPb_%s_cent%d_Nch1.pdf",spec.Data(),centClass));

  // TCanvas* cRat = new TCanvas("cRat", "cRat", 600, 400);
  //
  // TH1D* his[3][3] = {nullptr};
  // TH1D* hisAllCombi = (TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub",specLong.Data(),centClass));
  // TH1D* hisV = (TH1D*) fileV->Get(Form("%s_hFlow4_harm2_gap00_cent%d",specLong.Data(),centClass));
  //
  // TLegend* leg2 = new TLegend(0.12,0.75,0.55,0.88);
  // TF1 *fu1 = new TF1("fu1", "pol0", 0.4, 5.);
  //
  // gStyle->SetOptStat(kFALSE);
  // gStyle->SetOptFit(1);
  //
  // hisV->Divide(hisAllCombi);
  // hisV->SetMarkerStyle(kFullSquare);
  // hisV->SetMarkerColor(kBlack);
  //
  // hisV->Draw("h");
  // hisV->GetYaxis()->SetRangeUser(0.75,1.25);
  // hisV->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
  // hisV->GetYaxis()->SetTitle("VP / ZM ");
  //
  // leg2->SetFillColorAlpha(0.0,0.0);
  // leg2->SetBorderSize(0);
  // gStyle->SetLegendTextSize(0.02);
  // leg2->SetHeader(legString.Data());
  // leg2->Draw("same");
  //
  // hisV->Fit(fu1,"R");
  // fu1->SetLineColor(kRed);
  // fu1->Draw("same");
  //
  // cRat->SaveAs(Form("ptdif/ratios/v24_PbPb_%s_cent%d_VV.pdf",spec.Data(),centClass));


}
