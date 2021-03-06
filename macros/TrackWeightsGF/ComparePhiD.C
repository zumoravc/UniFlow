Int_t iNumRuns = 32; Int_t iRunList[] = {
  265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265335, 265334, 265332, 265309
};

TString sPath = "/home/alidock/ana/output/pPb_LHC16q/CENT_woSSD/train_2867";

Color_t color[] = {kRed, kBlue, kGreen, kOrange, kPink-4, kCyan+1, kMagenta+2, kViolet-6, kAzure+6, kGray+1, kBlack, kSpring+9, kPink-9, kGreen+3};

void ComparePhiD(){
  TFile* firstRun = TFile::Open(Form("%s/%d/AnalysisResults.root",sPath.Data(),iRunList[0]),"READ");
  if(!firstRun) { printf("Run %d | Input file with weights not found\n",iRunList[0]); return; }
  firstRun->cd("UniFlow");
  TList* listFR = (TList*) gDirectory->Get("Weights_UniFlow");
  if(!listFR) { printf("Run %d | TList with weights not found\n",iRunList[0]); gDirectory->ls() ;return; }
  TH2D* h2WeightsInputFR = (TH2D*) listFR->FindObject("fh2WeightsRefs");
  TH1D* hPhiFR = (TH1D*) h2WeightsInputFR->ProjectionX("hPhiFR");
  hPhiFR->Rebin(2);
  hPhiFR->Scale(1/hPhiFR->GetEntries());
  printf("Number of entries of Run %d: %f\n", iRunList[0], hPhiFR->GetEntries());

  TCanvas* can = new TCanvas("can", "", 600, 400);
  TMultiGraph *mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.42,0.22,0.82,0.38);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment | LHC16q - CENT_woSSD | pPb @ 5.02 TeV");
  leg->SetNColumns(3);
  Int_t counter = 0;

  for(Int_t iRun = 1; iRun < iNumRuns; ++iRun) {
      can->Clear();
      printf(" ### Run %d (out of %d)\n",iRun+1,iNumRuns);
      TFile* fileTemp = TFile::Open(Form("%s/%d/AnalysisResults.root",sPath.Data(),iRunList[iRun]),"READ");
      if(!fileTemp) { printf("Run %d | Input file with weights not found\n",iRunList[iRun]); continue; }
      fileTemp->cd("UniFlow");
      TList* listTemp = (TList*) gDirectory->Get("Weights_UniFlow");
      if(!listTemp) { printf("Run %d | TList with weights not found\n",iRunList[iRun]); gDirectory->ls() ;return; }
      TH2D* h2WeightsInput = (TH2D*) listTemp->FindObject("fh2WeightsRefs");
      TH1D* hPhi = (TH1D*) h2WeightsInput->ProjectionX("hPhi",0, -1);
      Double_t nOfEntries = hPhi->GetEntries() / 10000000;
      if(nOfEntries < 1) counter++;
      printf("Number of entries of Run %d (divided by 10M): %f \n\n", iRunList[iRun], nOfEntries );
      hPhi->Rebin(2);
      hPhi->Scale(1/hPhi->GetEntries());
      hPhi->SetTitle(Form("Run %d | LHC16q | pPb @ 5.02 TeV; #varphi; First run / this run",iRunList[iRun]));
      // hPhi->Draw("h");
      // can->SaveAs(Form("RBR_pPb/%d.pdf",iRunList[iRun]));

      hPhi->Divide(hPhiFR);

      gStyle->SetPalette(kRainBow);

      TGraph *thisRun = new TGraph((TH1D*) hPhi);
      // thisRun->SetLineColor(iRun+33);
      if(nOfEntries < 1) thisRun->SetLineWidth(2);
      mg->Add(thisRun);
      leg->AddEntry(thisRun,Form("%d",iRunList[iRun]),"l");
    }

    TCanvas* canMg = new TCanvas("canMg", "can", 600, 400);
    mg->GetXaxis()->SetRangeUser(0,TMath::TwoPi());
    mg->SetMinimum(0.65);
    mg->SetMaximum(1.08);
    mg->GetXaxis()->SetTitle("#varphi");
    mg->GetYaxis()->SetTitle("First run / this run");
    mg->Draw("AL PLC PFC");
    leg->Draw("same");
    canMg->SaveAs("RBR_pPb/CENT_woSSD.pdf");

    printf("\n\n COUNTER: %d \n", counter);

    return;
}
