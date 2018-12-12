void CompareYou()
{
  TFile* fileYou = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/you/Analysis_vnPsim.root","READ");
  if(!fileYou) { printf("ERROR: fileYou not found!\n"); return; }

  // TFile* fileMe = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/flow-modes/cross-charged-you/15o_hi_pass1/mixed/Processed.root","READ");
  TFile* fileMe = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/processUniFlow/mixed_test/Processed.root","READ");
  if(!fileMe) { printf("ERROR: fileMe not found!\n"); return; }

  TString sOutput = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/processUniFlow/mixed_test/compYou/";

  // preparing canvas
  TCanvas* can_v422 = new TCanvas("can_v422","can_v422", 1000,1000);
  for(Int_t cent(4); cent > -1; --cent)
  {
    const char* sNameYou = Form("hisv4psi2_cent%d",cent);
    const char* sNameMe = Form("Charged_<<3>>(4,-2,-2)_Pos_sample0_mult%d_rebin_px",cent);

    TH1D* histYou = (TH1D*) fileYou->Get(sNameYou);
    if(!histYou) { printf("ERROR: histYou '%s' not found!\n",sNameYou); fileYou->ls(); return; }

    TH1D* histMe = (TH1D*) fileMe->Get(sNameMe);
    if(!histMe) { printf("ERROR: histMe '%s' not found!\n", sNameMe); fileMe->ls(); return; }

    histMe->SetLineColor(kOrange+7);
    histMe->SetMarkerColor(kOrange+7);
    histMe->SetMarkerStyle(kOpenSquare);

    can_v422->cd();
    histYou->Draw("same");
    histMe->Draw("same");
    can_v422->SaveAs(Form("%s/comp_v422.pdf",sOutput.Data()));

  }

  TCanvas* can_v422gap10 = new TCanvas("can_v422gap10","can_v422gap10", 1000,1000);
  for(Int_t cent(4); cent > -1; --cent)
  {
    const char* sNameYou = Form("hisv4psi2Gap10_cent%d",cent);
    const char* sNameMe = Form("Charged_<<3>>(4,-2,-2)_2sub(1)_Pos_sample0_mult%d_rebin_px",cent);

    TH1D* histYou = (TH1D*) fileYou->Get(sNameYou);
    if(!histYou) { printf("ERROR: histYou '%s' not found!\n",sNameYou); fileYou->ls(); return; }

    TH1D* histMe = (TH1D*) fileMe->Get(sNameMe);
    if(!histMe) { printf("ERROR: histMe '%s' not found!\n", sNameMe); fileMe->ls(); return; }

    histMe->SetLineColor(kOrange+7);
    histMe->SetMarkerColor(kOrange+7);
    histMe->SetMarkerStyle(kOpenSquare);

    can_v422gap10->cd();
    histYou->Draw("same");
    histMe->Draw("same");
    can_v422gap10->SaveAs(Form("%s/comp_v422gap10.pdf",sOutput.Data()));
  }

  TCanvas* can_v532gap10 = new TCanvas("can_v532gap10","can_v532gap10", 1000,1000);
  for(Int_t cent(4); cent > -1; --cent)
  {
    const char* sNameYou = Form("hisv5psi32Gap10_cent%d",cent);
    const char* sNameMe = Form("Charged_<<3>>(5,-3,-2)_2sub(1)_Pos_sample0_mult%d_rebin_px",cent);

    TH1D* histYou = (TH1D*) fileYou->Get(sNameYou);
    if(!histYou) { printf("ERROR: histYou '%s' not found!\n",sNameYou); fileYou->ls(); return; }

    TH1D* histMe = (TH1D*) fileMe->Get(sNameMe);
    if(!histMe) { printf("ERROR: histMe '%s' not found!\n", sNameMe); fileMe->ls(); return; }

    histMe->SetLineColor(kOrange+7);
    histMe->SetMarkerColor(kOrange+7);
    histMe->SetMarkerStyle(kOpenSquare);

    can_v532gap10->cd();
    histYou->Draw("same");
    histMe->Draw("same");
    can_v532gap10->SaveAs(Form("%s/comp_v532gap10.pdf",sOutput.Data()));
  }




  TCanvas* can_v633gap10 = new TCanvas("can_v633gap10","can_v633gap10", 1000,1000);
  for(Int_t cent(4); cent > -1; --cent)
  {
    const char* sNameYou = Form("hisv6psi3Gap10_cent%d",cent);
    const char* sNameMe = Form("Charged_<<3>>(6,-3,-3)_2sub(1)_Pos_sample0_mult%d_rebin_px",cent);

    TH1D* histYou = (TH1D*) fileYou->Get(sNameYou);
    if(!histYou) { printf("ERROR: histYou '%s' not found!\n",sNameYou); fileYou->ls(); return; }

    TH1D* histMe = (TH1D*) fileMe->Get(sNameMe);
    if(!histMe) { printf("ERROR: histMe '%s' not found!\n", sNameMe); fileMe->ls(); return; }

    histMe->SetLineColor(kOrange+7);
    histMe->SetMarkerColor(kOrange+7);
    histMe->SetMarkerStyle(kOpenSquare);

    can_v633gap10->cd();
    histYou->Draw("same");
    histMe->Draw("same");
    can_v633gap10->SaveAs(Form("%s/comp_v633gap10.pdf",sOutput.Data()));
  }

  return;
}
