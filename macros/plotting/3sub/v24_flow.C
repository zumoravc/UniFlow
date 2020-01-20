char sides[4] = "LMR";
Color_t color[] = {kRed, kBlue, kGreen, kOrange, kPink-4, kCyan+1, kMagenta+2, kViolet-6, kAzure+6, kSpring+9, kPink-9, kGreen+3};

void v24_flow(Int_t part){

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

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  TLegend* leg = new TLegend(0.1,0.82,0.85,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetFillStyle(0);

  TString legString = "ALICE experiment, PbPb @ 5.02TeV, LHC15o, ";
  switch(part){
    case 0: legString += "h^{#pm}"; break;
    case 1: legString += "#pi^{#pm}"; break;
    case 2: legString += "K^{#pm}"; break;
    case 3: legString += "p/#bar{p}"; break;
  }
  legString += "\n |#eta| < 0.8";

  leg->SetHeader(legString.Data());
  gStyle->SetLegendTextSize(0.03);

  can->Divide(4,2);
  can->cd(1);
  TLatex* l = new TLatex();
  // l->SetTextAlign(23);
  l->SetTextSize(0.073);
  l->DrawLatex(0.1,0.6,"ALICE experiment");
  l->DrawLatex(0.1,0.5,"PbPb @ 5.02TeV (LHC15o)");
  TString lat = "";
  switch(part){
    case 0: lat += "h^{#pm}"; break;
    case 1: lat += "#pi^{#pm}"; break;
    case 2: lat += "K^{#pm}"; break;
    case 3: lat += "p/#bar{p}"; break;
  }
  lat += ", |#eta| < 0.8";
  l->DrawLatex(0.1,0.4,Form("%s",lat.Data()));



  // l->SetTextSize(0.05);
  l->Draw();
  for(Int_t i(0); i < 7; i++){
    can->cd(i+2);

    TGraphErrors* allCombi = new TGraphErrors((TH1D*) fileIn->Get(Form("%s_hFlow4_harm2_gap(08,08)_cent%d_3sub",specLong.Data(),i)));
    allCombi->SetMarkerStyle(kOpenSquare);
    allCombi->SetMarkerColor(kRed);
    allCombi->SetLineColor(kRed);

    allCombi->SetMinimum(0.0);
    allCombi->SetMaximum(0.35);

    TString cent = "VOM cent: ";
    switch(i){
      case 0: cent += "0-5%"; break;
      case 1: cent += "5-10%"; break;
      case 2: cent += "10-20%"; break;
      case 3: cent += "20-30%"; break;
      case 4: cent += "30-40%"; break;
      case 5: cent += "40-50%"; break;
      case 6: cent += "50-60%"; break;
    }

    allCombi->GetXaxis()->SetTitle("p_{#it{T}} [GeV/#it{c}]");
    allCombi->GetYaxis()->SetTitle("v_{2}{4}");
    allCombi->SetTitle(Form("%s",cent.Data()));
    allCombi->GetXaxis()->SetRangeUser(0, 10);

    allCombi->Draw("ap");
    // leg->Draw("same");
  }
  can->SaveAs(Form("v24_%s.pdf",specLong.Data()));



}
