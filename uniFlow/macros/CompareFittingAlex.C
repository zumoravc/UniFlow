/*
Macro for comparison of results of my fitting procedure with the one provided by Alex
and implemented in method ExtractFlowK0sAlex()
- commit ac92ad6a67511bac6652fbc93d409e6f0d175f35

Author: Vojtech Pacik, NBI (vojtech.pacik@cern.ch)
*/

void CompareFittingAlex()
{
  TString outDirectory = "/Users/vpacik/NBI/Flow/uniFlow/testing";
  const nBinsMult = 4;


  TFile* fileAlex = new TFile("fitAlex/fittingAlex.root","READ");
  TFile* fileMine = new TFile("fitMine/fittingMine.root","READ");

  Color_t colMine = kRed;
  Color_t colAlex = kBlue;

  if(!fileAlex || !fileMine) return;

  fileAlex->ls();

  TH1D* hMine = 0x0;
  TH1D* hAlex = 0x0;
  TH1D* hRatio = 0x0;

  TCanvas* can = new TCanvas("can","Can");
  TLine* lUnity = new TLine(0.,1.,20.,1.);
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
      legend->AddEntry(hMine,"Vojtech","pel");
      legend->AddEntry(hAlex,"Alex","pel");
    }

    hMine->SetLineColor(colMine);
    hMine->SetMarkerColor(colMine);
    hAlex->SetLineColor(colAlex);
    hAlex->SetMarkerColor(colAlex);

    hRatio = (TH1D*) hMine->Clone(Form("hRatio_%d",mult));
    hRatio->SetTitle("Ratio Vojtech / Alex");
    // hRatio->SetLineColor(kBlack);
    hRatio->Divide(hAlex);
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
