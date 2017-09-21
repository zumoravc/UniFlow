/*
Macro for comparison of results of my fitting procedure with the one provided by Alex
and implemented in method ExtractFlowK0sAlex()
- commit ac92ad6a67511bac6652fbc93d409e6f0d175f35

Author: Vojtech Pacik, NBI (vojtech.pacik@cern.ch)
*/

void CompareK0sSelectionAlex()
{
  TString outDirectory = "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/comparison";
  const nBinsMult = 4;

  TFile* filePrel = new TFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA/plots/Processed.root","READ");
  TFile* fileAlex = new TFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_Alex/plots/Processed.root","READ");
  TFile* fileAN = new TFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_AN/plots/Processed.root","READ");
  TFile* filePID = new TFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_pidAlone/plots/Processed.root","READ");
  TFile* fileDecay = new TFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_decay5/plots/Processed.root","READ");
  TFile* filePIDtrackingDecay = new TFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay/plots/Processed.root","READ");

  Color_t colAlex = kBlue;
  Color_t colAN = kRed;
  Color_t colPID = kGreen+2;
  Color_t colDecay = kOrange;
  Color_t colPrel = 38;
  Color_t colPIDtrackingDecay = kMagenta+2;

  if(!fileAlex || !filePID || !fileAN || !filePIDtrackingDecay ||! filePrel) return;


  gSystem->mkdir(outDirectory.Data(),kTRUE);

  fileAlex->ls();

  TH1D* hAlex = 0x0;
  TH1D* hAN = 0x0;
  TH1D* hPID = 0x0;
  TH1D* hDecay = 0x0;
  TH1D* hPIDtrackingDecay = 0x0;
  TH1D* hPrel = 0x0;

  TCanvas* can = new TCanvas("can","Can");
  TLine* lUnity = new TLine(0.,1.,20.,1.);
  lUnity->SetLineColor(kBlack);
  lUnity->SetLineStyle(kDashed);
  lUnity->SetLineWidth(1);

  TLegend* legend = new TLegend(0.13,0.88,0.3,0.68);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0,0);

  for(Short_t mult(0); mult < nBinsMult; mult++)
  {

    hAlex = (TH1D*) fileAlex->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));
    hAN = (TH1D*) fileAN->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));
    hPID = (TH1D*) filePID->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));
    hDecay = (TH1D*) fileDecay->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));
    hPrel = (TH1D*) filePrel->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));
    hPIDtrackingDecay = (TH1D*) filePIDtrackingDecay->Get(Form("hFlow2_K0s_harm2_gap08_mult%d",mult));


    if(!hAN || !hPID || !hAlex || !hPIDtrackingDecay) return;

    if(mult == 0)
    {
      legend->AddEntry(hAlex,"Alex","pel");
      legend->AddEntry(hAN,"AN","pel");
      legend->AddEntry(hPrel,"Prel","pel");
      legend->AddEntry(hPID,"PID","pel");
      legend->AddEntry(hDecay,"Decay5","pel");
      legend->AddEntry(hPIDtrackingDecay,"PIDTrackingDecay","pel");
    }

    hAlex->SetLineColor(colAlex);
    hAlex->SetMarkerColor(colAlex);
    hAN->SetLineColor(colAN);
    hAN->SetMarkerColor(colAN);
    hPID->SetLineColor(colPID);
    hPID->SetMarkerColor(colPID);
    hDecay->SetLineColor(colDecay);
    hDecay->SetMarkerColor(colDecay);
    hPIDtrackingDecay->SetLineColor(colPIDtrackingDecay);
    hPIDtrackingDecay->SetMarkerColor(colPIDtrackingDecay);
    hPrel->SetMarkerColor(colPrel);
    hPrel->SetLineColor(colPrel);

    hRatio = (TH1D*) hAN->Clone(Form("hRatio_%d",mult));
    hRatio->SetTitle("Ratio X / Alex");
    // hRatio->SetLineColor(kBlack);
    hRatio->Divide(hAlex);
    hRatio->SetMinimum(0.);
    hRatio->SetMaximum(2.);

    hRatio2 = (TH1D*) hPID->Clone(Form("hRatio2_%d",mult));
    hRatio2->SetTitle("Ratio X / Alex");
    // hRatio2->SetLineColor(kBlack);
    hRatio2->Divide(hAlex);
    // hRatio2->SetMinimum(0.7);
    // hRatio2->SetMaximum(1.3);

    hRatio3 = (TH1D*) hPIDtrackingDecay->Clone(Form("hRatio3_%d",mult));
    hRatio3->SetTitle("Ratio X / Alex");
    // hRatio2->SetLineColor(kBlack);
    hRatio3->Divide(hAlex);
    // hRatio3->SetMinimum(0.7);
    // hRatio3->SetMaximum(1.3);

    hRatio4 = (TH1D*) hDecay->Clone(Form("hRatio4_%d",mult));
    hRatio4->SetTitle("Ratio X / Alex");
    // hRatio4->SetLineColor(kBlack);
    hRatio4->Divide(hAlex);
    // hRatio4->SetMinimum(0.7);
    // hRatio4->SetMaximum(1.3);

    hRatio5 = (TH1D*) hPrel->Clone(Form("hRatio5_%d",mult));
    hRatio5->SetTitle("Ratio X / Alex");
    // hRatio5->SetLineColor(kBlack);
    hRatio5->Divide(hAlex);
    // hRatio5->SetMinimum(0.7);
    // hRatio5->SetMaximum(1.3);

    can->Clear();
    can->Divide(1,2);

    can->cd(1);
    hAN->DrawCopy("e1");
    hAlex->DrawCopy("same e1");
    hPID->DrawCopy("same e1");
    hDecay->DrawCopy("same e1");
    hPIDtrackingDecay->DrawCopy("same e1");
    hPrel->DrawCopy("same e1");
    legend->Draw();

    can->cd(2);
    hRatio->DrawCopy("e1");
    hRatio2->DrawCopy("e1 same");
    hRatio3->DrawCopy("e1 same");
    hRatio4->DrawCopy("e1 same");
    hRatio5->DrawCopy("e1 same");
    lUnity->Draw("same");

    can->SaveAs(Form("%s/mult%d.pdf",outDirectory.Data(),mult));

  }

}
