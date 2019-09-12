


TString sInputFileRecon = TString("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTest/UniFlow.root");
TString sOutputFilePath = TString("/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/");


void PrepareSystUnc()
{
  // loading flow
  TFile* fInputFileRecon = new TFile(sInputFileRecon.Data(),"READ");
  if(!fInputFileRecon->IsOpen()) return;
  // fInputFileRecon->ls();

  TCanvas* canFinal = new TCanvas("canFinal","canFinal",1200,1200);
  canFinal->Divide(2,2);

  const Int_t iNumCent = 4;

  Double_t dSyst_K0s[iNumCent] = { 0.02, 0.02, 0.11, 10 };


  Short_t iCent = 0;

  TH1D* hFlowK0s = (TH1D*) fInputFileRecon->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",iCent));
  if(!hFlowK0s) { printf("No K0s\n"); return; }
  TH1D* hFlowLambda = (TH1D*) fInputFileRecon->Get(Form("hFlow2_Lambda_harm2_gap08_mult%d",iCent));
  if(!hFlowLambda) { printf("No Lambda\n"); return; }
  TH1D* hFlowPhi = (TH1D*) fInputFileRecon->Get(Form("hFlow2_Phi_harm2_gap08_mult%d",iCent));
  if(!hFlowPhi) { printf("No Phi\n"); return; }

  TH1D* hSyst_K0s = (TH1D*) hFlowK0s->Clone(Form("hSyst_K0s_harm2_gap08_mult%d",iCent));
  hSyst_K0s->Reset();

  printf("bins: %d",hSyst_K0s->GetNbinsX());
  for(Int_t iBin(1); iBin < hSyst_K0s->GetNbinsX()+1; iBin++)
  {
    // hSyst_K0s->SetBinContent(iBin,dSyst_K0s[iCent]*hFlowK0s->GetBinContent(iBin));
    hSyst_K0s->SetBinContent(iBin,hFlowK0s->GetBinContent(iBin));
    hSyst_K0s->SetBinError(iBin,dSyst_K0s[iCent]*hFlowK0s->GetBinContent(iBin));
    hSyst_K0s->SetFillColor(kRed);
  }
  hSyst_K0s->Draw();


  canFinal->cd(iCent+1);
  // hSyst_K0s->Draw("same p2");
  hSyst_K0s->Draw("hist p2 e2");
  hFlowK0s->Draw("same e p");




  return;
}
