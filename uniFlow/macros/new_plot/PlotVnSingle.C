void PlotVnSingle()
{
  TString sOutDir = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/output/plots_new_2/";
  gSystem->mkdir(sOutDir,1);

  TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/output/Processed.root","READ");
  fileIn->ls();

  Int_t iMultCent = 1;
  TString sCent = "20-30% V0M";

  TString sSpecies = "Charged";  TString sSpesLabel = "h^{#pm}"; Color_t color = kGray+2;
  // TString sSpecies = "Pion";  TString sSpesLabel = "#pi^{#pm}"; Color_t color = kRed+1;
  // TString sSpecies = "K0s";  TString sSpesLabel = "K^{0}_{S}"; Color_t color = kGreen+2;
  // TString sSpecies = "Lambda";  TString sSpesLabel = "#Lambda(#bar{#Lambda})"; Color_t color = kOrange+1;
  // TString sSpecies = "Phi"; TString  sSpesLabel = "#phi";   Color_t color = kMagenta;

  Int_t markers[] = {
    kFullCircle,
    kFullSquare,
    kFullDiamond,
    kOpenCircle,
    kOpenSquare,
    kOpenDiamond
  };


  TString sNames[] = {
    Form("hFlow2_%s_harm2_gap-10_cent%d",sSpecies.Data(), iMultCent),
    Form("hFlow2_%s_harm2_gap00_cent%d",sSpecies.Data(), iMultCent),
    Form("hFlow2_%s_harm2_gap04_cent%d",sSpecies.Data(), iMultCent),
    Form("hFlow4_%s_harm2_gap-10_cent%d",sSpecies.Data(), iMultCent),
    Form("hFlow4_%s_harm2_gap00_cent%d",sSpecies.Data(), iMultCent),
    Form("hFlow4_%s_harm2_gap04_cent%d",sSpecies.Data(), iMultCent),
  };

  TString sLabels[] = {
    "v_{2}{2}",
    "v_{2}{2,|#Delta#eta| > 0}",
    "v_{2}{2,|#Delta#eta| > 0.4}",
    "v_{2}{4}",
    "v_{2}{4,|#Delta#eta| > 0}",
    "v_{2}{4,|#Delta#eta| > 0.4}"
  };


  TCanvas* can = new TCanvas("can","can",800,1000);
  can->cd();
  TH1* frame = (TH1*) gPad->DrawFrame(0.0,0.0,6.0,0.3);
  frame->SetTitle("; #it{p}_{T} (GeV/#it{c}); v_{2}");

  TLegend* leg = new TLegend(0.14,0.65,0.4,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader(Form("%s (%s)",sSpesLabel.Data(),sCent.Data()));
  leg->Draw();

  Int_t iNumPlots = sizeof(sNames) / sizeof(sNames[0]);
  for(Int_t i(0); i < iNumPlots; ++i) {
    TH1D* hist = (TH1D*) fileIn->Get(sNames[i].Data());
    if(!hist) continue;

    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markers[i]);
    hist->Draw("hist ape1x0 same");

    leg->AddEntry(hist,sLabels[i].Data(),"p");
    can->SaveAs(Form("%s/can_%s_cent%d_stage%d.pdf",sOutDir.Data(),sSpecies.Data(),iMultCent,i),"pdf");
  }

  return;
}
