char sides[4] = "LMR";
Color_t color[] = {kRed, kBlue, kGreen, kOrange, kPink-4, kCyan+1, kMagenta+2, kViolet-6, kAzure+6, kSpring+9, kPink-9, kGreen+3};

void v24(Int_t part, Int_t centClass){

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

  TFile* fileV = TFile::Open("/home/alidock/output/published/Vojta/Processed.root","READ");
  if(!fileV) {printf("File VV not opened! \n"); return;}

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.12,0.7,0.85,0.88);

  TGraphErrors* gap08 = new TGraphErrors((TH1D*) fileV->Get(Form("%s_hFlow4_harm2_gap00_cent%d",specLong.Data(),centClass)));
  gap08->SetMarkerStyle(kOpenDiamond);
  gap08->SetMarkerColor(kBlue);
  gap08->SetLineColor(kBlue);
  mg->Add(gap08);
  leg->AddEntry(gap08, "v_{2}{4,|#Delta#eta| > 0.0} - VP", "p");

  TGraphErrors* allCombi = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub",specLong.Data(),centClass)));
  allCombi->SetMarkerStyle(kOpenSquare);
  allCombi->SetMarkerColor(kRed);
  allCombi->SetLineColor(kRed);
  mg->Add(allCombi);
  leg->AddEntry(allCombi, "v_{2}{4}_{3sub} - ZM", "p");


  mg->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
  mg->GetYaxis()->SetTitle("v_{2}{4}");
  mg->GetXaxis()->SetRangeUser(0, 10);

  mg->Draw("ap");
  mg->SetMinimum(0.0);
  mg->SetMaximum(0.35);

  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
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
  gStyle->SetLegendTextSize(0.025);
  // leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->Draw("same");

  can->SaveAs(Form("ptdif/v24_PbPb_%s_cent%d_VV.pdf",spec.Data(),centClass));

  TCanvas* cRat = new TCanvas("cRat", "cRat", 600, 400);

  TH1D* his[3][3] = {nullptr};
  TH1D* hisAllCombi = (TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub",specLong.Data(),centClass));
  TH1D* hisV = (TH1D*) fileV->Get(Form("%s_hFlow4_harm2_gap00_cent%d",specLong.Data(),centClass));

  TLegend* leg2 = new TLegend(0.12,0.75,0.55,0.88);
  TF1 *fu1 = new TF1("fu1", "pol0", 0.4, 5.);

  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(1);

  hisV->Divide(hisAllCombi);
  hisV->SetMarkerStyle(kFullSquare);
  hisV->SetMarkerColor(kBlack);

  hisV->Draw("h");
  hisV->GetYaxis()->SetRangeUser(0.75,1.25);
  hisV->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
  hisV->GetYaxis()->SetTitle("VP / ZM ");

  leg2->SetFillColorAlpha(0.0,0.0);
  leg2->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.02);
  leg2->SetHeader(legString.Data());
  leg2->Draw("same");

  hisV->Fit(fu1,"R");
  fu1->SetLineColor(kRed);
  fu1->Draw("same");

  cRat->SaveAs(Form("ptdif/ratios/v24_PbPb_%s_cent%d_VV.pdf",spec.Data(),centClass));

  TMultiGraph* comp = new TMultiGraph();
  TGraphErrors* graph[3][3] = {nullptr};
  TLegend* legGr = new TLegend(0.12,0.75,0.75,0.88);

  legGr->SetFillColorAlpha(0.0,0.0);
  legGr->SetBorderSize(0);
  gStyle->SetLegendTextSize(0.02);
  legGr->SetHeader(legString.Data());


  for(Int_t poiPos(0); poiPos < 3; poiPos++){
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      TCanvas* can2 = new TCanvas("can2", "can2", 600, 400);

      his[poiPos][twoPos] = (TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub_poi_%c_two_%c",specLong.Data(),centClass,sides[poiPos],sides[twoPos]));
      if(!his[poiPos][twoPos]) {printf("Problem in poi %c two %c \n", sides[poiPos],sides[twoPos]); fileIn->ls(); return;}

      graph[poiPos][twoPos] = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub_poi_%c_two_%c",specLong.Data(),centClass,sides[poiPos],sides[twoPos])));

      if(!graph[poiPos][twoPos]) {printf("Problem in poi %c two %c \n", sides[poiPos],sides[twoPos]); fileIn->ls(); return;}

      graph[poiPos][twoPos]->SetMarkerStyle(kOpenCircle);
      graph[poiPos][twoPos]->SetMarkerColor(color[poiPos + 3*twoPos]);
      graph[poiPos][twoPos]->SetMarkerSize(1.);
      graph[poiPos][twoPos]->SetLineColor(color[poiPos + 3*twoPos]);
      comp->Add(graph[poiPos][twoPos]);

      legGr->AddEntry(graph[poiPos][twoPos],Form("v_{2}{4} (poi: %c, two: %c, central #eta: -0.4, 0.4)",sides[poiPos],sides[twoPos]),"p");


      gStyle->SetOptStat(kFALSE);
      gStyle->SetOptFit(1);

      his[poiPos][twoPos]->Divide(hisAllCombi);
      his[poiPos][twoPos]->SetMarkerStyle(kFullSquare);
      his[poiPos][twoPos]->SetMarkerColor(kBlack);

      his[poiPos][twoPos]->Draw("h");
      his[poiPos][twoPos]->GetYaxis()->SetRangeUser(0.85,1.15);
      his[poiPos][twoPos]->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
      his[poiPos][twoPos]->GetYaxis()->SetTitle(Form("v_{2}{4} 3 sub poi %c two %c / v_{2}{4} 3 sub ", sides[poiPos],sides[twoPos]));

      leg2->Draw("same");

      his[poiPos][twoPos]->Fit(fu1,"R");
      fu1->SetLineColor(kRed);
      fu1->Draw("same");

      can2->SaveAs(Form("ptdif/ratios/v24_PbPb_%s_cent%d_poi_%c_two_%c.pdf",spec.Data(),centClass,sides[poiPos],sides[twoPos]));

      delete can2;
    }
  }

 TCanvas* cGr = new TCanvas("cGr", "cGr", 600, 400);
 comp->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
 comp->GetYaxis()->SetTitle("v_{2}{4}_{3 sub}");
 comp->GetXaxis()->SetRangeUser(0, 10);

 comp->Draw("ap");
 comp->SetMinimum(0.0);
 comp->SetMaximum(0.15);

 legGr->SetNColumns(2);
 legGr->Draw("same");

 cGr->SaveAs(Form("ptdif/v24_PbPb_%s_cent%d_combi.pdf",spec.Data(),centClass));

}
