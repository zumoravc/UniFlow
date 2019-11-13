char sides[4] = "LMR";
Color_t color[] = {kRed, kBlue, kGreen, kOrange, kPink-4, kCyan+1, kMagenta+2, kViolet-6, kAzure+6, kSpring+9, kPink-9, kGreen+3};
Int_t mult[] = {0,5,10,20,30,40,50,60,70,80};


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
  TGraphErrors* graph[3][3] = {nullptr};
  TLegend* leg = new TLegend(0.12,0.7,0.85,0.88);


  for(Int_t poiPos(0); poiPos < 3; poiPos++){
    for(Int_t twoPos(0); twoPos < 3; twoPos++){
      graph[poiPos][twoPos] = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub_poi_%c_two_%c",specLong.Data(),centClass,sides[poiPos],sides[twoPos])));
      if(!graph[poiPos][twoPos]) {printf("Problem in poi %c two %c \n", sides[poiPos],sides[twoPos]); fileIn->ls(); return;}

      graph[poiPos][twoPos]->SetMarkerStyle(kOpenCircle);
      graph[poiPos][twoPos]->SetMarkerColor(color[poiPos + 3*twoPos]);
      graph[poiPos][twoPos]->SetMarkerSize(1.);
      graph[poiPos][twoPos]->SetLineColor(color[poiPos + 3*twoPos]);
      mg->Add(graph[poiPos][twoPos]);

      leg->AddEntry(graph[poiPos][twoPos],Form("v_{2}{4} (poi: %c, two: %c, central #eta: -0.4, 0.4)",sides[poiPos],sides[twoPos]),"p");
    }
  }

  TGraphErrors* gap0 = new TGraphErrors((TH1D*) fileVV->Get(Form("v24_Gap0_%d_%d_Syst",mult[centClass],mult[centClass+1])));
  TGraphErrors* gap1 = new TGraphErrors((TH1D*) fileVV->Get("v24_Gap1_5_10_Syst"));
  gap0->SetMarkerStyle(kOpenDiamond);
  gap0->SetMarkerColor(kBlack);
  gap0->SetLineColor(kBlack);
  gap1->SetMarkerStyle(kOpenSquare);
  gap1->SetMarkerColor(kGray+2);
  gap1->SetLineColor(kGray+2);
  if(part == 0)
  {
    mg->Add(gap0);
    mg->Add(gap1);

    leg->AddEntry(gap0, "v_{2}{4, |#Delta#eta| > 0.0} - VV");
    leg->AddEntry(gap1, "v_{2}{4, |#Delta#eta| > 1.0} - VV");
  }

  mg->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
  mg->GetYaxis()->SetTitle("v_{2}{4}_{3 sub}");
  mg->GetXaxis()->SetRangeUser(0, 10);

  mg->Draw("ap");
  mg->SetMinimum(0.0);
  mg->SetMaximum(0.25);

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
    case 0: legString += "0-5%"; break;
    case 1: legString += "5-10%"; break;
    case 2: legString += "10-20%"; break;
    case 3: legString += "20-30%"; break;
    case 4: legString += "30-40%"; break;
    case 5: legString += "40-50%"; break;
    case 6: legString += "50-60%"; break;
    case 7: legString += "60-70%"; break;
  }

  leg->SetHeader(legString.Data());
  gStyle->SetLegendTextSize(0.025);
  // leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->Draw("same");

  can->SaveAs(Form("ptdif/ratios/v24_PbPb_%s_cent%d.pdf",spec.Data(),centClass));

  TGraphErrors merged = nullptr;
  for(Int_t poiPos(0); poiPos < 3; poiPos++){
    for(Int_t twoPos(0); twoPos < 3; twoPos++){

    }
  }

}
