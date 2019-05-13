void SetEx(TGraphErrors* gae, Double_t Ex);

void PlotMixedGraph()
{

  // TString sCorr = "<<3>>(4,-2,-2)"; TString sLabel = "v_{4,-2,-2}";
  // TString sCorr = "<<3>>(5,-3,-2)"; TString sLabel = "v_{5,-3,-2}";
  // TString sCorr = "<<3>>(6,-3,-3)"; TString sLabel = "v_{6,-3,-3}";


  TString sGap = "_gap00";
  TString sCorr = "Mergedv422"; TString sLabel = "v_{4,-2,-2} {|#Delta#eta| > 0} ";  TString sCorrSyst = "MergedSystv422"; Double_t dYlow = -0.01; Double_t dYhigh = 0.1;
  // TString sCorr = "Mergedv523"; TString sLabel = "v_{5,-3,-2} {|#Delta#eta| > 0} "; TString sCorrSyst = "MergedSystv523"; Double_t dYlow = -0.01; Double_t dYhigh = 0.1;
  // TString sCorr = "Mergedv633"; TString sLabel = "v_{6,-3,-3} {|#Delta#eta| > 0} ";TString sCorrSyst = "MergedSystv633"; Double_t dYlow = -0.01; Double_t dYhigh = 0.06;

  // TString sCorr = "<<3>>(4,-2,-2)_2sub(0.4)"; TString sLabel = "v_{4,-2,-2} {|#Delta#eta| > 0.4} ";
  // TString sCorr = "<<3>>(5,-3,-2)_2sub(0.4)"; TString sLabel = "v_{5,-3,-2} {|#Delta#eta| > 0.4} ";
  // TString sCorr = "<<3>>(6,-3,-3)_2sub(0.4)"; TString sLabel = "v_{6,-3,-3} {|#Delta#eta| > 0.4} ";

  TString sCent[] = {
    "0-5% V0M",
    "5-10% V0M",
    "10-20% V0M",
    "20-30% V0M",
    "30-40% V0M",
    "40-50% V0M",
    "50-60% V0M",
  };

  TString sCentLab[] = {
    "0-5",
    "5-10",
    "10-20",
    "20-30",
    "30-40",
    "40-50",
    "50-60"
  };

  Int_t iNumMult = sizeof(sCent) / sizeof(sCent[0]);

  TString sInputPath = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/All/graphs_bk_2/");
  TString sOutDir = sInputPath+"/../plots/";


  TString sSpecies[] = {
    // "Charged",
    "Pion",
    "Kaon",
    "K0s",
    "Proton",
    "Lambda",
    "Phi",
  };

  TString sSpesLabels[] = {
    // "h^{#pm}",
    "#pi^{#pm}",
    "K^{#pm}",
    "K^{0}_{S}",
    "p^{#pm}",
    "#Lambda(#bar{#Lambda})",
    "#phi"
  };

  Color_t colors[] = {
    // kGray+2,
    kRed+1,
    kGreen+1,
    kGreen+2,
    kBlue,
    kOrange+1,
    kMagenta
  };

  Int_t markers[] = {
    // kOpenCircle,
    kFullCircle,
    kFullSquare,
    kFullTriangleDown,
    kFullTriangleUp,
    kFullCross,
    kFullDiamond
  };

  Double_t markerSizes[] = {
    // 1.0,
    1.0,
    1.0,
    1.1,
    1.1,
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

  for(Int_t iMult(0); iMult < iNumMult; ++iMult) {

    TCanvas* can = new TCanvas("can","can",600,600);
    can->cd();
    TH1* frame = (TH1*) gPad->DrawFrame(0.0,dYlow,6.0,dYhigh);
    frame->SetTitle("; #it{p}_{T} (GeV/#it{c});");

    TLegend* leg = new TLegend(0.14,0.5,0.4,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.0,0.0);
    leg->SetTextFont(43);
    leg->SetTextSize(24);

    Int_t iNumPlots = sizeof(sSpecies) / sizeof(sSpecies[0]);
    for(Int_t i(0); i < iNumPlots; ++i) {


      TGraphErrors* histSyst = nullptr;
      TFile* fileInSyst = TFile::Open(Form("%s/%s%s.root",sInputPath.Data(),sCorrSyst.Data(),sSpecies[i].Data()),"READ");
      if(fileInSyst) {
          TString sName = Form("%spT_%s%s",sCorrSyst.Data(),sCentLab[iMult].Data(),sGap.Data());
          histSyst = (TGraphErrors*) fileInSyst->Get(sName.Data());
          if(!histSyst) { printf("HistSyst '%s' not found!\n",sName.Data()); fileInSyst->ls(); return; }

          histSyst->SetLineColor(colors[i]);
          histSyst->SetFillColorAlpha(0.0,0.0);
          // hist->SetMarkerColor(colors[i]);
          // hist->SetMarkerStyle(markers[i]);
          // hist->SetMarkerSize(markerSizes[i]);
          histSyst->Draw("same5");
      }

      TFile* fileIn = TFile::Open(Form("%s/%s%s.root",sInputPath.Data(),sCorr.Data(),sSpecies[i].Data()),"READ");
      if(!fileIn) {continue; }
      // TString sName = Form("%s_%s_mult%d",sSpecies[i].Data(), sCorr.Data(), iMult);
      TString sName = Form("%spT_%s%s",sCorr.Data(),sCentLab[iMult].Data(),sGap.Data());

      TGraphErrors* hist = (TGraphErrors*) fileIn->Get(sName.Data());
      if(!hist) { printf("Hist '%s' not found!\n",sName.Data()); fileIn->ls(); return; }

      hist->SetLineColor(colors[i]);
      hist->SetMarkerColor(colors[i]);
      hist->SetMarkerStyle(markers[i]);
      hist->SetMarkerSize(markerSizes[i]);
      SetEx(hist,0.0);
      hist->Draw("same p");

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


void SetEx(TGraphErrors* gae, Double_t Ex)
{
   Int_t np = gae->GetN();
   for (Int_t i=0; i<np; i++) {
      gae->SetPointError(i,Ex,gae->GetErrorY(i));
      // gae->SetPointError(i,Ex);
   }
}
