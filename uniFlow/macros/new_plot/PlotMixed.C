void PlotMixed()
{

  // TString sCorr = "<<3>>(4,-2,-2)"; TString sLabel = "v_{4,-2,-2}"; TString sGap = "gap-10";
  // TString sCorr = "<<3>>(4,-2,-2)_2sub(0)"; TString sLabel = "v_{4,-2,-2} {|#Delta#eta| > 0} "; TString sGap = "gap0";
  // TString sCorr = "<<3>>(4,-2,-2)_2sub(0.4)"; TString sLabel = "v_{4,-2,-2} {|#Delta#eta| > 0.4} "; TString sGap = "gap4";

  // TString sCorr = "<<3>>(5,-3,-2)"; TString sLabel = "v_{5,-3,-2}"; TString sGap = "gap-10";
  // TString sCorr = "<<3>>(5,-3,-2)_2sub(0)"; TString sLabel = "v_{5,-3,-2} {|#Delta#eta| > 0} "; TString sGap = "gap0";
  // TString sCorr = "<<3>>(5,-3,-2)_2sub(0.4)"; TString sLabel = "v_{5,-3,-2} {|#Delta#eta| > 0.4} "; TString sGap = "gap4";

  // TString sCorr = "<<3>>(6,-3,-3)"; TString sLabel = "v_{6,-3,-3}"; TString sGap = "gap-10";
  // TString sCorr = "<<3>>(6,-3,-3)_2sub(0)"; TString sLabel = "v_{6,-3,-3} {|#Delta#eta| > 0} "; TString sGap = "gap0";
  TString sCorr = "<<3>>(6,-3,-3)_2sub(0.4)"; TString sLabel = "v_{6,-3,-3} {|#Delta#eta| > 0.4} "; TString sGap = "gap4";

  TString sCent[] = {
    "0-5% V0M",
    "5-10% V0M",
    "10-20% V0M",
    "20-30% V0M",
    "30-40% V0M",
    "40-50% V0M",
    "50-60% V0M",
  };

  Int_t iNumMult = sizeof(sCent) / sizeof(sCent[0]);

  TString sInputPath = Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/lhc15o/pass1_AOD194/afterMemLeak/flow/output_mixed/%s",sGap.Data());
  // TString sOutDir = sInputPath+"/plots";
  TString sOutDir = sInputPath+"/plots_all";


  // Int_t iMultCent = 0; TString sCent = "10-20% V0M";
  // Int_t iMultCent = 1; TString sCent = "20-30% V0M";
  // Int_t iMultCent = 2; TString sCent = "30-40% V0M";
  // Int_t iMultCent = 3; TString sCent = "40-50% V0M";

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
    kGreen+1,
    kGreen+2,
    kBlue,
    kOrange+1,
    kMagenta
  };

  Int_t markers[] = {
    kOpenCircle,
    kFullCircle,
    kFullTriangleDown,
    kFullSquare,
    kFullTriangleUp,
    kFullCross,
    kFullDiamond,
  };

  Double_t markerSizes[] = {
    1.0,
    1.0,
    1.1,
    1.1,
    1.0,
    1.2,
    1.4
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

  gSystem->mkdir(sOutDir.Data(),1);

  TFile* fileIn = TFile::Open(Form("%s/Processed.root",sInputPath.Data()),"READ");
  fileIn->ls();


  for(Int_t iMult(0); iMult < iNumMult; ++iMult) {
    TCanvas* can = new TCanvas("can","can",800,1000);
    can->cd();
    TH1* frame = (TH1*) gPad->DrawFrame(0.0,-0.02,6.0,0.15);
    frame->SetTitle("; #it{p}_{T} (GeV/#it{c});");

    TLegend* leg = new TLegend(0.14,0.65,0.4,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.0,0.0);

    Int_t iNumPlots = sizeof(sSpecies) / sizeof(sSpecies[0]);
    for(Int_t i(0); i < iNumPlots; ++i) {

      TString sName = Form("%s_%s_mult%d",sSpecies[i].Data(), sCorr.Data(), iMult);

      TH1D* hist = (TH1D*) fileIn->Get(sName.Data());
      if(!hist) { printf("Hist '%s' not found!\n",sName.Data()); continue; }

      hist->SetLineColor(colors[i]);
      hist->SetMarkerColor(colors[i]);
      hist->SetMarkerStyle(markers[i]);
      hist->SetMarkerSize(markerSizes[i]);
      hist->Draw("hist ape1x0 same");

      leg->AddEntry(hist,sSpesLabels[i].Data(),"p");
    }

    leg->SetHeader(Form("%s (%s)",sLabel.Data(),sCent[iMult].Data()));
    leg->Draw();
    can->SaveAs(Form("%s/can_%s_cent%d.pdf",sOutDir.Data(),sLabel.Data(),iMult),"pdf");

    delete frame;
    delete leg;
    delete can;
  }



  return;
}
