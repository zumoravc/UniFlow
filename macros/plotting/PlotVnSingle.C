void PlotVnSingle()
{
  TString sOutDir = "/home/alidock/ana/output/test/";
  gSystem->mkdir(sOutDir,1);

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/test/Processed.root","READ");
  fileIn->ls();

  Int_t iMultCent = 1;
  //TString sCent = "20-30% V0M";

  TString sSpecies = "Refs";  TString sSpesLabel = "RFP"; //Color_t color = kBlue+2;
  //TString sSpecies = "Charged";  TString sSpesLabel = "h^{#pm}"; Color_t color = kGray+2;
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

  Color_t color[] = {kRed, kBlue, kGreen, kRed, kBlue, kGreen};


  TString sNames[] = {
    "Refs_hFlow2_harm2_gap10","Refs_hFlow2_harm3_gap10","Refs_hFlow2_harm4_gap10"
    // Form("hFlow2_%s_harm2_gap-10_cent%d",sSpecies.Data(), iMultCent),
    // Form("hFlow2_%s_harm2_gap00_cent%d",sSpecies.Data(), iMultCent),
    // Form("hFlow2_%s_harm2_gap04_cent%d",sSpecies.Data(), iMultCent),
    // Form("hFlow4_%s_harm2_gap-10_cent%d",sSpecies.Data(), iMultCent),
    // Form("hFlow4_%s_harm2_gap00_cent%d",sSpecies.Data(), iMultCent),
    // Form("hFlow4_%s_harm2_gap04_cent%d",sSpecies.Data(), iMultCent),
  };

  TString sLabels[] = {
    "v_{2}{2,|#Delta#eta| > 1.0}",
    "v_{3}{2,|#Delta#eta| > 1.0}",
    "v_{4}{2,|#Delta#eta| > 1.0}",
    "v_{2}{2,|#Delta#eta| > 1.0} PUBLISHED",
    "v_{3}{2,|#Delta#eta| > 1.0} PUBLISHED",
    "v_{4}{2,|#Delta#eta| > 1.0} PUBLISHED"
    // "v_{2}{2,|#Delta#eta| > 0}",
    // "v_{2}{2,|#Delta#eta| > 0.4}",
    // "v_{2}{4}",
    // "v_{2}{4,|#Delta#eta| > 0}",
    // "v_{2}{4,|#Delta#eta| > 0.4}"
  };


  // Data
  Double_t v22Gap10Run2[10]={0.0283859, 0.0456604, 0.0655068, 0.0870721, 0.099105, 0.104143, 0.10286, 0.0974591, 0.0888104};
  Double_t v22Gap10Run2CombErr[10]={0.000711585, 0.000940029, 0.00105135, 0.00137951, 0.00158637, 0.00172523, 0.00187963, 0.00236759, 0.0045766};
  Double_t v32Gap10Run2[10]={0.0206723, 0.0231991, 0.0279915, 0.0309315, 0.0331377, 0.0323443, 0.0270494, 0.0276913};
  Double_t v32Gap10Run2CombErr[10]={0.00114039, 0.0013322, 0.00136947, 0.00155171, 0.00174155, 0.00201624, 0.00309182, 0.00561635};
  Double_t v42Gap10Run2[10]={0.0114652, 0.0129068, 0.0138326, 0.0156838, 0.017018, 0.0162893, 0.0177113, 0.0198688, 0.00916289};
  Double_t v42Gap10Run2CombErr[10]={-1, 0.0366982, 0.0586746, 0.0779774, 0.0867288, 0.0890264, 0.0835125, 0.0806163, 0.0464366};

  TH1D* v22Gap10 = (TH1D*) fileIn->Get(sNames[1].Data());
  v22Gap10->Clear();
  TH1D* v32Gap10 = (TH1D*) fileIn->Get(sNames[1].Data());
  v32Gap10->Clear();
  TH1D* v42Gap10 = (TH1D*) fileIn->Get(sNames[1].Data());
  v42Gap10->Clear();

  for(int i = 1; i < 10; i++)
  {
    v22Gap10->SetBinContent(i, v22Gap10Run2[i]);
    v22Gap10->SetBinError(i, v22Gap10Run2CombErr[i]);
    v32Gap10->SetBinContent(i, v32Gap10Run2[i]);
    v32Gap10->SetBinError(i, v32Gap10Run2CombErr[i]);
    v42Gap10->SetBinContent(i, v42Gap10Run2[i]);
    v42Gap10->SetBinError(i, v42Gap10Run2CombErr[i]);
  }

  TCanvas* can = new TCanvas("can","can",800,1000);
  can->cd();
  TH1* frame = (TH1*) gPad->DrawFrame(0.0,0.0,80.0,0.15);
  //frame->SetTitle("; #it{p}_{T} (GeV/#it{c}); v_{2}");
  frame->SetTitle("; centrality; v_{n}{2}");


  TLegend* leg = new TLegend(0.14,0.65,0.4,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  //leg->SetHeader(Form("%s (%s)",sSpesLabel.Data(),sCent.Data()));
  leg->SetHeader(Form("%s",sSpesLabel.Data()));
  leg->Draw();

  Int_t iNumPlots = sizeof(sNames) / sizeof(sNames[0]);
  for(Int_t i(0); i < iNumPlots; ++i) {
  for(Int_t i(0); i < iNumPlots; ++i) {
    if(i < 3)
      TH1D* hist = (TH1D*) fileIn->Get(sNames[i].Data());
    if (i==3) TH1D* hist = v22Gap10;
    if (i==4) TH1D* hist = v32Gap10;
    if (i==5) TH1D* hist = v42Gap10;
    //if(!hist) continue;


    hist->SetLineColor(color[i]);
    hist->SetMarkerColor(color[i]);
    hist->SetMarkerStyle(markers[i]);
    hist->Draw("hist ape1x0 same");

    leg->AddEntry(hist,sLabels[i].Data(),"p");
    can->SaveAs(Form("%s/can_%s_cent%d_stage%d.pdf",sOutDir.Data(),sSpecies.Data(),iMultCent,i),"pdf");
  }

  return;
}
