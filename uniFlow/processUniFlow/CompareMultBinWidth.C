TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

void CompareMultBinWidth()
{
  TFile* unitBin = TFile::Open("~/NBI/Flow/uniFlow/results/KchK0s-rerun/preliminary/output_fixedMultBins/Processed.root");
  TFile* fixedBin = TFile::Open("~/NBI/Flow/uniFlow/results/KchK0s-rerun/preliminary_fixedMultBins/plots/Processed.root");

  unitBin->ls();

  const char* histNames[] = {
    "hFlow2_Refs_harm2_gap08",
    "hFlow2_K0s_harm2_gap08_mult0",
    "hFlow2_K0s_harm2_gap08_mult1",
    "hFlow2_K0s_harm2_gap08_mult2",
    "hFlow2_K0s_harm2_gap08_mult3",
    "hFlow2_K0s_harm2_gap08_mult4",
    "hFlow2_K0s_harm2_gap08_mult5",
    "hFlow2_Kaon_harm2_gap08_cent0",
    "hFlow2_Kaon_harm2_gap08_cent1",
    "hFlow2_Kaon_harm2_gap08_cent2",
    "hFlow2_Kaon_harm2_gap08_cent3",
    "hFlow2_Kaon_harm2_gap08_cent4",
    "hFlow2_Kaon_harm2_gap08_cent5",
  };

  const Int_t numHists = sizeof(histNames)/sizeof(histNames[0]);


  gSystem->mkdir("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s-rerun/multBinsComparison");

  for(Int_t hist(0); hist < numHists; hist++)
  {
    const char* histName = histNames[hist];

    TH1D* histUnit = (TH1D*) unitBin->Get(histName);
    TH1D* histFixed_old = (TH1D*) fixedBin->Get(histName);

    TH1D* histFixed = (TH1D*) histUnit->Clone("histFixed");
    histFixed->Reset();


    for(Int_t bin(1); bin < histFixed->GetNbinsX()+1; bin++)
    {
      histFixed->SetBinContent(bin,histFixed_old->GetBinContent(bin));
      histFixed->SetBinError(bin,histFixed_old->GetBinError(bin));
    }

    hRatio = DivideHistos(histUnit,histFixed,0);
    hRatio->SetStats(0);
    hRatio->SetMinimum(0.95);
    hRatio->SetMaximum(1.05);

    histUnit->SetStats(0);
    // histUnit->SetMaximum(0.5);
    histUnit->SetMinimum(0.);
    histFixed->SetLineColor(kRed);

    TCanvas* can = new TCanvas("can","can");
    can->Divide(2,1);
    can->cd(1);
    histUnit->DrawCopy("e1");
    histFixed->DrawCopy("same e1");
    can->cd(2);
    // histFixed->DrawCopy("e1");
    // can->cd(3);
    hRatio->DrawCopy();
    // can->cd(4);
    // histFixed_old->DrawCopy("e1");
    can->SaveAs(Form("~/NBI/Flow/uniFlow/results/KchK0s-rerun/multBinsComparison_2/%s.pdf",histName),"pdf");

  }




}
// =============================================================================
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  ratio->Reset();

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < binsDenom+1; iBin++)
  {
    if(iBin > binsNom) break;

    dContNom = nom->GetBinContent(iBin);
    dErrNom = nom->GetBinError(iBin);
    dContDenom = denom->GetBinContent(iBin);
    dErrDenom = denom->GetBinError(iBin);

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBin,dContRatio);
    ratio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
