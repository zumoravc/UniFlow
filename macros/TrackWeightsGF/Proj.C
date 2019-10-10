Int_t iNumRuns = 80; Int_t iRunList[] = {
  246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851,
  246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804,
  246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495,
  246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272,
  246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153,
  246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049,
  246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952,
  245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683
};

void Proj(){
  TFile* file = TFile::Open("/home/alidock/ana/output/YZ_results/TrackingEfficiency_PbPb5TeV_AMPT.root","READ");
  // TFile* file = TFile::Open("/home/alidock/ana/output/YZ_results/PhiWeight_LHC15o_AMPT_HIR.root","READ");
  if(!file) { printf("File not opened!"); return; }

  // TCanvas *can = new TCanvas("can", "", 600, 600);
  // TH3F* h3Weights = (TH3F*) file->Get(Form("fPhiWeight_%d",iRunList[0]));
  // // TH3F* h3Weights = (TH3F*) file->Get("fPhiWeight_245683");
  // // if(!h3Weights) { printf("Histgram not found!"); return; }
  // TH1D* phi = h3Weights->ProjectionX("phi", 0, -1, 19, 21);
  // // TH1D* phi = h3Weights->ProjectionZ("phi", 0, -1, 0, -1);
  // // printf("8 - 10 : bins: %d %d \n", phi->FindBin(8), phi->FindBin(10));
  // // printf("-10 - -8 : bins: %d %d \n", phi->FindBin(-10), phi->FindBin(-8));
  // // printf("-1 - 1 : bins: %d %d \n", phi->FindBin(-1), phi->FindBin(1));
  // gStyle->SetOptStat(kFALSE);
  // phi->SetTitle("#varphi distribution, |#eta| < 0.8, 8 < v_{z} < 10; #varphi; Entries");
  // phi->Draw("h");
  // can->SaveAs("PhiDis_Vz_8_10.pdf");
  // // h3Weights->Draw();

  TCanvas *can = new TCanvas("can", "", 600, 600);
  TH3F* h3Weights = (TH3F*) file->Get(Form("eff_LHC15o_AMPT_%d",iRunList[0]));
    TH1D* phi = h3Weights->ProjectionX("phi", 0, -1, 0, -1);
  // TH1D* phi = h3Weights->ProjectionZ("phi", 0, -1, 0, -1);
  // printf("8 - 10 : bins: %d %d \n", phi->FindBin(8), phi->FindBin(10));
  // printf("-10 - -8 : bins: %d %d \n", phi->FindBin(-10), phi->FindBin(-8));
  // printf("-1 - 1 : bins: %d %d \n", phi->FindBin(-1), phi->FindBin(1));
  gStyle->SetOptStat(kFALSE);
  phi->SetTitle("NUE AMPT, |#eta| < 1.0, |v_{z}| < 10; p_{T}; Entries");
  phi->Draw("h");
  can->SaveAs("NUE_AMPT.pdf");
  // h3Weights->Draw();
  return;
}
