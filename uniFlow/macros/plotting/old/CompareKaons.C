TH1F* Divide(TH1* Nom = 0x0, TH1* Denom = 0x0);

void CompareKaons()
{
  const char* sInputKchRun1 = "/Users/vpacik/NBI/Flow/results/uniFlow_ver4_V0A/run1_comparison/HEPdata_extracted.root";

  const char* sInputKchDefault = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTestKaon/UniFlow.root";
  const char* sInputK0sDefault = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTestKaon/UniFlow.root";
  const char* sInputK0sDefault_app = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/Alex/UniFlow.root";
  // const char* sInputK0sDefault_app = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTest/UniFlow.root";
  const char* sInputK0sPbPb = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/Alex/UniFlow.root";
  const char* sInputK0sTight = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/loosePhi_tightV0s/kaon/UniFlow.root";
  const char* sInputK0sDCA = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/dcaZ/UniFlow.root";
  const char* sInputK0sFB = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit_withNUA/UniFlow_fb768_kaons/UniFlow.root";
  // const char* sInputK0sNoArm = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/noArm/UniFlow.root";

  // const char* sInputKch = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit_withNUA/UniFlow_fb768_kaons/UniFlow_fb768_kaons.root";



  // const char* sOutput = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTestKaon/";
  const char* sOutput = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTestKaon/";
  // const char* sOutput = "/Users/vpacik/NBI/Flow/results/NUA_cor/merged_16qt_noFASTt/fitTest/";



  const Int_t iNumCent = 4;
  Double_t dK0sSyst[iNumCent] = {0.07,0.1,0.12,0.18};
  TString sCent[iNumCent] = {"0-20","20-40","40-60","60-100"};

  Int_t markKaonRun1 = kFullSquare;
  Int_t markKaon = kFullCircle;
  Int_t markK0s = kOpenCircle;
  Int_t markK0sTight = kOpenSquare;
  Int_t markK0sDCA = kOpenSquare;
  Int_t markK0sFB = kOpenSquare;
  // Int_t markK0sNoArm = kOpenCircle;

  Color_t colKaon = kRed;
  Color_t colKaonRun1 = kGray+1;
  Color_t colK0s = kBlue;
  Color_t colK0sTight = kGreen+2;
  Color_t colK0sDCA = kBlack;
  Color_t colK0sFB =  kGreen+2;

  TFile* fileKchRun1 = new TFile(sInputKchRun1,"READ"); if(!fileKchRun1) { printf("Kch Run1 file not found!\n"); return; }
  TFile* fileK0sDefault = new TFile(sInputK0sDefault,"READ"); if(!fileK0sDefault) { printf("K0s file not found!\n"); return; }
  TFile* fileK0sDefault_app = new TFile(sInputK0sDefault_app,"READ"); if(!fileK0sDefault_app) { printf("K0s file not found!\n"); return; }
  TFile* fileKchDefault = new TFile(sInputKchDefault,"READ"); if(!fileKchDefault) { printf("K0s file not found!\n"); return; }
  TFile* fileK0sPbPb = new TFile(sInputK0sPbPb,"READ"); if(!fileK0sPbPb) { printf("K0s file not found!\n"); return; }
  TFile* fileK0sTight = new TFile(sInputK0sTight,"READ"); if(!fileK0sTight) { printf("K0s file not found!\n"); return; }
  TFile* fileK0sDCA = new TFile(sInputK0sDCA,"READ"); if(!fileK0sDCA) { printf("K0s file not found!\n"); return; }
  TFile* fileK0sFB = new TFile(sInputK0sFB,"READ"); if(!fileK0sFB) { printf("K0s file not found!\n"); return; }

  // fileKch->ls();
  // fileK0s->ls();

  TH1* hFlow_KchRun1 = 0x0;
  TH1* hFlow_KchDefault = 0x0;
  TH1* hFlow_K0sAlex = 0x0;
  TH1* hFlow_K0sDefault = 0x0;
  TH1* hFlow_K0sDefault_app = 0x0;
  TH1* hFlow_K0sTight = 0x0;
  TH1* hFlow_K0sDCA = 0x0;
  TH1* hFlow_K0sFB = 0x0;
  TH1F* hRatioDefault = 0x0;
  TH1F* hRatioPbPb = 0x0;
  TH1F* hRatioTight = 0x0;
  TH1F* hRatioDCA = 0x0;
  TH1F* hRatioFB = 0x0;

  TH1F* hRatioLooseTight = 0x0;


  TH1D* hK0sDefaultSyst = 0x0;
  TH1D* hK0sDefaultSyst_app = 0x0;

  TLegend* legend = new TLegend(0.2,0.65,0.7,0.84);
  legend->SetFillColorAlpha(0,0);
  legend->SetBorderSize(0);
  legend->SetTextSize(gStyle->GetTextSize()*0.8);

  TLegend* legendFinal = new TLegend(0.2,0.65,0.7,0.84);
  legendFinal->SetFillColorAlpha(0,0);
  legendFinal->SetBorderSize(0);
  legendFinal->SetTextSize(gStyle->GetTextSize()*0.8);

  TCanvas* cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600);
  cfig->Divide(2,1);
  TH1* h1 = 0x0;
  TH1* h2 = 0x0;

  TCanvas* cK0s = new TCanvas("cK0s","cK0s",800,800);
  TH1* hK0s = 0x0;

  TCanvas* cComp = new TCanvas("cComp","cComp",800,800);
  TH1* hComp = 0x0;

  TCanvas* cComp2 = new TCanvas("cComp2","cComp2",800,800);
  TH1* hComp2 = 0x0;



  TF1* fitDefault = new TF1("fitDefault","[0]",0.75,3);
  fitDefault->SetLineColor(colKaon);
  TF1* fitPbPb = new TF1("fitPbPb","[0]",0.75,3);
  fitPbPb->SetLineColor(colK0s);
  TF1* fitTight = new TF1("fitTight","[0]",0.75,3);
  fitTight->SetLineColor(colK0sTight);
  TF1* fitDCA = new TF1("fitDCA","[0]",0.75,3);
  fitDCA->SetLineColor(colK0sDCA);
  TF1* fitFB = new TF1("fitFB","[0]",0.75,3);
  fitFB->SetLineColor(colK0sFB);

  TF1* fitLooseTight = new TF1("fitLooseTight","[0]",0.75,3);



  TLine* unity = new TLine();
  unity->SetLineStyle(7);
  TText* text = new TText();

  for(Int_t iCent(0); iCent < iNumCent; iCent++)
  {
    h1 = cfig->cd(1)->DrawFrame(0,0,4,0.4);
    h2 = cfig->cd(2)->DrawFrame(0,0,4,2);

    hK0s = cK0s->cd()->DrawFrame(0,0,4,2);

    h1->SetTitle("; #it{p}_{T} (GeV/#it{c}); v2");
    h2->SetTitle("K0s/Kch; #it{p}_{T} (GeV/#it{c})");

    hFlow_KchRun1 = (TH1D*) fileKchRun1->Get(Form("hFlow2_Kaon_harm2_gap08_cent%d",iCent)); if(!hFlow_KchRun1) { printf("Hist Kch not found\n"); return; }
    hFlow_KchDefault = (TH1D*) fileKchDefault->Get(Form("hFlow2_Kaon_harm2_gap08_cent%d",iCent)); if(!hFlow_KchRun1) { printf("Hist Kch not found\n"); return; }
    // hFlow_Kch = (TH1D*) fileKch->Get(Form("hFlow2_Kaon_harm2_gap08_cent%d",iCent)); if(!hFlow_Kch) { printf("Hist Kch not found\n"); return; }
    hFlow_K0sDefault = (TH1D*) fileK0sDefault->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent)); if(!hFlow_K0sDefault) { printf("Hist K0s not found\n"); return; }
    hFlow_K0sDefault_app = (TH1D*) fileK0sDefault_app->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent)); if(!hFlow_K0sDefault) { printf("Hist K0s not found\n"); return; }
    hFlow_K0sPbPb = (TH1D*) fileK0sPbPb->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent)); if(!hFlow_K0sPbPb) { printf("Hist K0s not found\n"); return; }
    hFlow_K0sTight = (TH1D*) fileK0sTight->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent)); if(!hFlow_K0sTight) { printf("Hist K0s not found\n"); return; }
    hFlow_K0sDCA = (TH1D*) fileK0sDCA->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent)); if(!hFlow_K0sDCA) { printf("Hist K0s not found\n"); return; }
    hFlow_K0sFB = (TH1D*) fileK0sFB->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent)); if(!hFlow_K0sFB) { printf("Hist K0s not found\n"); return; }

    hRatioDefault = Divide(hFlow_K0sDefault,hFlow_KchRun1);
    hRatioPbPb = Divide(hFlow_K0sPbPb,hFlow_KchRun1);
    hRatioTight = Divide(hFlow_K0sTight,hFlow_KchRun1);
    hRatioDCA = Divide(hFlow_K0sDCA,hFlow_KchRun1);
    hRatioFB = Divide(hFlow_K0sFB,hFlow_KchRun1);
    hRatioLooseTight = Divide(hFlow_K0sDefault,hFlow_K0sPbPb);


    hK0sDefaultSyst = (TH1D*) hFlow_K0sDefault->Clone("syst");
    for(Short_t bin(1); bin < hK0sDefaultSyst->GetNbinsX()+1; bin++) hK0sDefaultSyst->SetBinError(bin,dK0sSyst[iCent]*hFlow_K0sDefault->GetBinContent(bin));

    hK0sDefaultSyst_app = (TH1D*) hFlow_K0sDefault_app->Clone("syst");
    for(Short_t bin(1); bin < hK0sDefaultSyst_app->GetNbinsX()+1; bin++) hK0sDefaultSyst_app->SetBinError(bin,dK0sSyst[iCent]*hFlow_K0sDefault_app->GetBinContent(bin));

    hFlow_K0sDefault->SetStats(0);
    hFlow_K0sPbPb->SetStats(0);
    hRatioDefault->SetStats(0);
    hRatioPbPb->SetStats(0);

    hFlow_K0sPbPb->SetMarkerStyle(markK0s);
    hFlow_K0sPbPb->SetMarkerColor(colK0s);
    hFlow_K0sPbPb->SetLineColor(colK0s);

    hFlow_K0sDefault->SetMarkerStyle(markKaon);
    hFlow_K0sDefault->SetMarkerColor(colKaon);
    hFlow_K0sDefault->SetLineColor(colKaon);

    hFlow_K0sDefault_app->SetMarkerStyle(markKaon);
    hFlow_K0sDefault_app->SetMarkerColor(colKaon);
    hFlow_K0sDefault_app->SetLineColor(colKaon);

    hFlow_KchRun1->SetMarkerStyle(markKaonRun1);
    hFlow_KchRun1->SetMarkerColor(colKaonRun1);
    hFlow_KchRun1->SetLineColor(colKaonRun1);

    hFlow_K0sTight->SetMarkerStyle(markK0sTight);
    hFlow_K0sTight->SetMarkerColor(colK0sTight);
    hFlow_K0sTight->SetLineColor(colK0sTight);

    hFlow_K0sDCA->SetMarkerStyle(markK0sDCA);
    hFlow_K0sDCA->SetMarkerColor(colK0sDCA);
    hFlow_K0sDCA->SetLineColor(colK0sDCA);

    hFlow_K0sFB->SetMarkerStyle(markK0sFB);
    hFlow_K0sFB->SetMarkerColor(colK0sFB);
    hFlow_K0sFB->SetLineColor(colK0sFB);

    hRatioDefault->SetMarkerStyle(markKaon);
    hRatioDefault->SetMarkerColor(colKaon);
    hRatioDefault->SetLineColor(colKaon);

    hRatioPbPb->SetMarkerStyle(markK0s);
    hRatioPbPb->SetMarkerColor(colK0s);
    hRatioPbPb->SetLineColor(colK0s);

    hRatioTight->SetMarkerStyle(markK0sTight);
    hRatioTight->SetMarkerColor(colK0sTight);
    hRatioTight->SetLineColor(colK0sTight);

    hRatioDCA->SetMarkerStyle(markK0sDCA);
    hRatioDCA->SetMarkerColor(colK0sDCA);
    hRatioDCA->SetLineColor(colK0sDCA);

    hRatioFB->SetMarkerStyle(markK0sFB);
    hRatioFB->SetMarkerColor(colK0sFB);
    hRatioFB->SetLineColor(colK0sFB);

    if(iCent == 0)
    {
      legend->AddEntry(hFlow_KchRun1,"Kch (Run1)","pel");
      legend->AddEntry(hFlow_K0sDefault,"K0s (default)","pel");
      // legend->AddEntry(hFlow_K0sTight,"K0s (tight)","pel");
      legend->AddEntry(hFlow_K0sDCA,"K0s (daughters dca-z)","pel");
      legend->AddEntry(hFlow_K0sPbPb,"K0s (Pb-Pb like)","pel");
      legend->AddEntry(hFlow_K0sFB,"K0s (Charged FB768)","pel");

      legendFinal->AddEntry(hFlow_KchRun1,"Kch (Run1)","pel");
      legendFinal->AddEntry(hFlow_K0sDefault,"K0s (default)","pel");
    }

    cK0s->cd();
    hRatioLooseTight->Draw("same p e1");
    hRatioLooseTight->Fit("fitLooseTight","RL");
    cK0s->SaveAs(Form("%s/CompKaons_Loose_Tight_cent%d.pdf",sOutput,iCent));

    cfig->cd(1);
    hFlow_KchRun1->Draw("same p e1");
    hFlow_K0sPbPb->Draw("same p e1");
    hFlow_K0sDefault->Draw("same p e1");
    // hFlow_K0sTight->Draw("same p e1");
    hFlow_K0sDCA->Draw("same p e1");
    hFlow_K0sFB->Draw("same p e1");
    legend->SetHeader(Form("Cent %s",sCent[iCent].Data()));
    legend->Draw();


    cfig->cd(2);
    hRatioDefault->Draw("same p e1");
    hRatioPbPb->Draw("same p e1");
    // hRatioTight->Draw("same p e1");
    hRatioDCA->Draw("same p e1");
    hRatioFB->Draw("same p e1");

    hRatioDefault->Fit("fitDefault","RL");
    hRatioPbPb->Fit("fitPbPb","RL");
    // hRatioTight->Fit("fitTight","RL");
    hRatioDCA->Fit("fitDCA","RL");
    hRatioFB->Fit("fitFB","RL");
    unity->DrawLine(0,1.,4.,1.);

    // text->SetTextColor(colKaon);
    // text->DrawTextNDC(0.3,0.2,Form("%g+-%g",fitDefault->GetParameter(0),fitDefault->GetParError(0)));
    // text->SetTextColor(colKaonRun1);
    // text->DrawTextNDC(0.3,0.15,Form("%g+-%g",fitPbPb->GetParameter(0),fitPbPb->GetParError(0)));

    cfig->SaveAs(Form("%s/CompKaons_cent%d.pdf",sOutput,iCent));


    hComp = cComp->cd()->DrawFrame(0,0,4,0.4);
    cComp->cd();
    hK0sDefaultSyst->SetFillColorAlpha(colKaon,0.5);
    hK0sDefaultSyst->Draw("same p e2");
    hFlow_K0sDefault->Draw("same p e1 x0");
    hFlow_KchRun1->SetLineColor(kBlack);
    hFlow_KchRun1->SetMarkerColor(kBlack);
    hFlow_KchRun1->Draw("same p e1 x0");
    legendFinal->SetHeader(Form("Cent %s",sCent[iCent].Data()));
    legendFinal->Draw();
    cComp->SaveAs(Form("%s/CompKaons_final_cent%d.pdf",sOutput,iCent));

    hComp2 = cComp2->cd()->DrawFrame(0,0,4,0.4);
    cComp2->cd();
    hK0sDefaultSyst_app->SetFillColorAlpha(colKaon,0.5);
    hK0sDefaultSyst_app->Draw("same p e2");
    hFlow_K0sDefault_app->Draw("same p e1 x0");
    hFlow_KchDefault->Draw("same p e1 x0");
    hFlow_KchRun1->SetLineColor(kBlack);
    hFlow_KchRun1->SetMarkerColor(kBlack);
    hFlow_KchRun1->Draw("same p e1 x0");
    legendFinal->SetHeader(Form("Cent %s",sCent[iCent].Data()));
    legendFinal->Draw();
    cComp2->SaveAs(Form("%s/CompKaons_final_Bin_cent%d.pdf",sOutput,iCent));
  }






  return;
}

TH1F* Divide(TH1* Nom, TH1* Denom)
{
  if(!Nom || !Denom) { printf("At least one of input not found\n"); return 0x0; }

  Int_t iNumBins = Nom->GetNbinsX();
  Int_t iNumBinsErr = Denom->GetNbinsX();
  // printf("%d | %d\n",hist->GetNbinsX(),histErr->GetNbinsX());
  printf("%d | %d\n",iNumBins,iNumBinsErr);

  if(iNumBins != iNumBinsErr) { printf("Not the same number of bins\n"); return 0x0; }

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dErr = 0;

  TH1F* Ratio = (TH1F*) Nom->Clone("Ratio");
  for(Int_t iBin(1); iBin < iNumBins+1; iBin++)
  {
    dContNom = Nom->GetBinContent(iBin);
    dErrNom = Nom->GetBinError(iBin);
    dContDenom = Denom->GetBinContent(iBin);
    dErrDenom = Denom->GetBinError(iBin);

    dErr = TMath::Power(dErrNom/dContDenom,2) + TMath::Power((dContNom*dErrDenom)/(dContDenom*dContDenom),2);
    // dErr -= 2*(dContNom*TMath::Power(dContDenom,-3)*dErrNom*dErrDenom);

    Ratio->SetBinContent(iBin,dContNom/dContDenom);
    Ratio->SetBinError(iBin,TMath::Sqrt(dErr));
    // hist->SetBinError(iBin,histErr->GetBinError(iBin));
  }

  return Ratio;
  // return hist;

}
