void v22pPb_ratios(){
  TFile* fileIn = TFile::Open("/home/alidock/ana/output/pPb_LHC16q/HADD/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileInNW = TFile::Open("/home/alidock/ana/output/pPb_LHC16q/HADD/processUniFlow/ProcessedNW.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  // TGraphErrors *nogap = new TGraphErrors((TH1D*) fileIn->Get("Charged_hFlow2_harm2_gap-10_cent0"));

  // leg->AddEntry(nogap,"v_{2}","p");
  // leg->AddEntry(gap00,"v_{2}(|#Delta#eta|>0.0)","p");
  // leg->AddEntry(gap08,"v_{2}(|#Delta#eta|>0.8)","p");
  // leg->AddEntry(gap10,"v_{2}(|#Delta#eta|>1.0)","p");
  // leg->Draw("same");

  TString particles[4] = {"Charged", "Pion", "Kaon", "Proton"};

  TCanvas* can2[3];
  TH1D* ratio;
  TH1D* denom;
  for(int i = 0; i < 3; i++)
  {
    particles[i] += "_hFlow2_harm2_gap-10_cent0";
    can2[i] = new TCanvas(Form("can_%d",i), "can2", 600, 400);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    ratio = (TH1D*) fileIn->Get(particles[i]);
    denom = (TH1D*) fileInNW->Get(particles[i]);
    ratio->Divide(denom);
    ratio->SetMarkerStyle(kFullSquare);
    ratio->SetMarkerColor(kBlack);
    ratio->SetLineColor(kBlack);
    ratio->SetTitle(";p_{#it{T}} [GeV/#it{c}]; v_{2}: with weights / without weights");
    ratio->SetTitle(particles[i]);
    // ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetRangeUser(0.2,5.);
    ratio->Draw("h");
    TF1 *fu1 = new TF1("fu1", "pol0", 0.4, 5.);
    ratio->Fit(fu1,"R");
    fu1->SetLineColor(kRed);
    fu1->Draw("same");
    can2[i]->SaveAs("ratio_v22_gap-1_" + particles[i]+ ".pdf");
  }

}
