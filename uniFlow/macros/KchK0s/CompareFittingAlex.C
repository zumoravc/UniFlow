/*
Macro for comparison of results of my fitting procedure with the one provided by Alex
and implemented in method ExtractFlowK0sAlex()
- commit ac92ad6a67511bac6652fbc93d409e6f0d175f35

Author: Vojtech Pacik, NBI (vojtech.pacik@cern.ch)
*/

TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor);


void CompareFittingAlex()
{
  TString inputDirectory = "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA";
  TString outDirectory = "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA";
  const nBinsMult = 4;


  TFile* fileAlex = TFile::Open(Form("%s/output/Processed.root",inputDirectory.Data()),"READ");
  TFile* fileMine = TFile::Open(Form("%s/output_seqFit/Processed.root",inputDirectory.Data()),"READ");

  Color_t colMine = kRed;
  Color_t colAlex = kBlue;

  if(!fileAlex || !fileMine) return;

  fileAlex->ls();

  TH1D* hMine = 0x0;
  TH1D* hAlex = 0x0;
  TH1D* hRatio = 0x0;

  TCanvas* can = new TCanvas("can","Can",400,800);
  TLine* lUnity = new TLine(0.2,1.,10.,1.);
  lUnity->SetLineColor(kBlack);
  lUnity->SetLineStyle(kDashed);
  lUnity->SetLineWidth(1);

  TLegend* legend = new TLegend(0.13,0.88,0.3,0.68);
  legend->SetBorderSize(0);

  for(Short_t mult(0); mult < nBinsMult; mult++)
  {

    hMine = (TH1D*) fileMine->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));
    hAlex = (TH1D*) fileAlex->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));

    if(!hMine || !hAlex) return;

    if(mult == 0)
    {
      legend->AddEntry(hMine,"sequential fit","pel");
      legend->AddEntry(hAlex,"one-step fit","pel");
    }

    hMine->SetLineColor(colMine);
    hMine->SetMarkerColor(colMine);
    hAlex->SetLineColor(colAlex);
    hAlex->SetMarkerColor(colAlex);


    hRatio = (TH1D*) hMine->Clone(Form("hRatio_%d",mult));
    hRatio->SetTitle("Ratio sequential / one-step fit");
    // hRatio->SetLineColor(kBlack);
    hRatio->Divide(hMine,hAlex,1,1,"b");

    // hRatio = DivideHistos(hMine,hAlex,1);

    hRatio->SetMinimum(0.8);
    hRatio->SetMaximum(1.2);

    can->Clear();
    can->Divide(1,2);

    can->cd(1);
    hMine->DrawCopy("e1");
    hAlex->DrawCopy("same e1");
    legend->Draw();

    can->cd(2);
    hRatio->DrawCopy("e1");
    lUnity->Draw("same");

    can->SaveAs(Form("%s/mult%d.pdf",outDirectory.Data(),mult));

  }

}

// ===========================================================================
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  ratio->Reset();

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < binsNom+1; iBin++)
  {
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
