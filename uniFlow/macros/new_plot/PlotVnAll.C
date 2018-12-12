void PlotVnAll()
{
  TString sOutDir = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/output_bk/plots_new/";
  gSystem->mkdir(sOutDir,1);

  TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/output_bk/Processed.root","READ");
  // TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/output_v3/Processed.root","READ");
  fileIn->ls();


  // Int_t iMultCent = 0; TString sCent = "10-20% V0M";
  // // Int_t iMultCent = 1; TString sCent = "20-30% V0M";
  // // Int_t iMultCent = 2; TString sCent = "30-40% V0M";
  // // Int_t iMultCent = 3; TString sCent = "40-50% V0M";

  TString sCent[] = {
    // "0-5% V0M",
    // "5-10% V0M",
    "10-20% V0M",
    "20-30% V0M",
    "30-40% V0M",
    "40-50% V0M",
    // "50-60% V0M",
  };

  Int_t iNumMult = sizeof(sCent) / sizeof(sCent[0]);

  TString sSpecies[] = {
    "Charged",
    "Pion",
    "Kaon",
    "K0s",
    "Proton",
    "Lambda",
    // "Phi",
  };

  TString sSpesLabels[] = {
    "h^{#pm}",
    "#pi^{#pm}",
    "K^{#pm}",
    "K^{0}_{S}",
    "p^{#pm}",
    "#Lambda(#bar{#Lambda})",
    "#phi"
  };

  Color_t colors[] = {
    kGray+2,
    kRed+1,
    kGreen+3,
    kGreen+1,
    kBlue,
    kOrange+1,
    kMagenta
  };

  Int_t markers[] = {
    kOpenCircle,
    kFullCircle,
    kFullTriangleDown,
    kFullTriangleUp,
    kFullSquare,
    kFullDiamond,
    kFullCross,
  };

  Double_t markerSizes[] = {
    1.0,
    1.0,
    1.1,
    1.1,
    0.8,
    1.4,
    1.2
  };


  // TString sNames[] = {
  //   Form("hFlow2_%s_harm2_gap-10_cent%d_sample0",sSpecies.Data(), iMultCent),
  //   Form("hFlow2_%s_harm2_gap00_cent%d_sample0",sSpecies.Data(), iMultCent),
  //   Form("hFlow2_%s_harm2_gap04_cent%d_sample0",sSpecies.Data(), iMultCent),
  //   Form("hFlow4_%s_harm2_gap-10_cent%d_sample0",sSpecies.Data(), iMultCent),
  //   Form("hFlow4_%s_harm2_gap00_cent%d_sample0",sSpecies.Data(), iMultCent),
  //   Form("hFlow4_%s_harm2_gap04_cent%d_sample0",sSpecies.Data(), iMultCent),
  // };

  // TString sLabels[] = {
  //   "v_{2}{2}",
  //   "v_{2}{2,|#Delta#eta| > 0}",
  //   "v_{2}{2,|#Delta#eta| > 0.4}",
  //   "v_{2}{4}",
  //   "v_{2}{4,|#Delta#eta| > 0}",
  //   "v_{2}{4,|#Delta#eta| > 0.4}"
  // };

  for(Int_t iMult(0); iMult < iNumMult; ++iMult) {

    TCanvas* can = new TCanvas("can","can",800,600);
    can->cd();
    TH1* frame = (TH1*) gPad->DrawFrame(0.0,0.0,6.0,0.3);
    frame->SetTitle("; #it{p}_{T} (GeV/#it{c}); v_{2}");

    TLegend* leg = new TLegend(0.14,0.65,0.4,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.0,0.0);

    Int_t iNumPlots = sizeof(sSpecies) / sizeof(sSpecies[0]);
    for(Int_t i(0); i < iNumPlots; ++i) {

      TString sName = Form("hFlow4_%s_harm2_gap-10_cent%d",sSpecies[i].Data(), iMult);

      TH1D* hist = (TH1D*) fileIn->Get(sName.Data());
      if(!hist) continue;

      hist->SetLineColor(colors[i]);
      hist->SetMarkerColor(colors[i]);
      hist->SetMarkerStyle(markers[i]);
      hist->SetMarkerSize(markerSizes[i]);
      hist->Draw("hist ape1x0 same");

      leg->AddEntry(hist,sSpesLabels[i].Data(),"p");
    }

    TString sLabel = "v_{2}{4}";
    leg->SetHeader(Form("%s (%s)",sLabel.Data(),sCent[iMult].Data()));
    leg->Draw();
    can->SaveAs(Form("%s/can_v24_all_cent%d.pdf",sOutDir.Data(),iMult),"pdf");

    delete leg;
    delete frame;
    delete can;
  }


  return;
}
