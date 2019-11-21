char sides[4] = "LMR";

void combinations(Int_t part, Int_t centClass){

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

  TFile* fileIn = TFile::Open("/home/alidock/output/LHC15o/train_7456/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TH1D* hisAllCombi = (TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub",specLong.Data(),centClass));

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  can->Divide(3,3);

  TH1D* his[3][3] = {nullptr};
  TF1 *fu1 = new TF1("fu1", "pol0", 0.4, 5.);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(1);

  TLegend* leg = new TLegend(0.1,0.11,0.85,0.18);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  // leg->SetFillStyle(0);
  TString legString = "ALICE experiment, PbPb @ 5.02TeV, LHC15o, ";
  switch(part){
    case 0: legString += "h^{#pm}, |#eta| < 0.8, V0M cent: "; break;
    case 1: legString += "#pi^{#pm}, |#eta| < 0.8, V0M cent: "; break;
    case 2: legString += "K^{#pm}, |#eta| < 0.8, V0M cent: "; break;
    case 3: legString += "p/#bar{p}, |#eta| < 0.8, V0M cent: "; break;
  }

  switch(centClass){
    case 0: legString += "5-10%"; break;
    case 1: legString += "10-20%"; break;
    case 2: legString += "20-30%"; break;
    case 3: legString += "30-40%"; break;
    case 4: legString += "40-50%"; break;
    case 5: legString += "50-60%"; break;
  }

  leg->SetHeader(legString.Data());
  gStyle->SetLegendTextSize(0.036);

  for(Int_t poiPos(0); poiPos < 3; poiPos++){
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      can->cd(1+twoPos+3*poiPos);
      his[poiPos][twoPos] = (TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub_poi_%c_two_%c",specLong.Data(),centClass,sides[poiPos],sides[twoPos]));
      if(!his[poiPos][twoPos]) {printf("Problem in poi %c two %c \n", sides[poiPos],sides[twoPos]); fileIn->ls(); return;}

      his[poiPos][twoPos]->Divide(hisAllCombi);
      his[poiPos][twoPos]->SetMarkerStyle(kFullSquare);
      his[poiPos][twoPos]->SetMarkerColor(kBlack);
      his[poiPos][twoPos]->SetMarkerSize(0.1);
      his[poiPos][twoPos]->SetLineColor(kBlack);

      his[poiPos][twoPos]->Draw("h");
      his[poiPos][twoPos]->GetYaxis()->SetRangeUser(0.9,1.1);
      his[poiPos][twoPos]->GetXaxis()->SetRangeUser(0.2,5.);
      his[poiPos][twoPos]->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
      his[poiPos][twoPos]->GetYaxis()->SetTitle(Form("v_{2}{4} 3 sub poi %c two %c / v_{2}{4} 3 sub ", sides[poiPos],sides[twoPos]));
      his[poiPos][twoPos]->SetTitle(Form("h^{#pm}: v_{2}{4} "));

      his[poiPos][twoPos]->Fit(fu1,"R");
      fu1->SetLineColor(kRed);
      leg->Draw("same");
      // fu1->Draw("same");
    }
  }

  can->SaveAs(Form("ptdif/v24_PbPb_%s_cent%d_combinations.pdf",spec.Data(),centClass));

}
