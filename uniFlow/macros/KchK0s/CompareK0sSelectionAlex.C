/*
Macro for comparison of results of my fitting procedure with the one provided by Alex
and implemented in method ExtractFlowK0sAlex()
- commit ac92ad6a67511bac6652fbc93d409e6f0d175f35

Author: Vojtech Pacik, NBI (vojtech.pacik@cern.ch)
*/

const Short_t iNumFiles = 13;
const Short_t nBinsMult = 4;
TString sMultLabel[nBinsMult] =
{
  "0-10",
  "10-20",
  "20-40",
  "40-60"
};

TString outDirectory = "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/comparison/compKchK0s";

Double_t dXlow = 0.; Double_t dXupp = 8.;
Double_t dYlow = 0.; Double_t dYupp = 0.4;


TString sInputFile[iNumFiles] = {
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/run1/HEPdata_extracted.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA/output/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA/output/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_NUAcor/output_run1/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_Alex/output_compAlex/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_AN/output_compAlex/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_decay5/output_run1/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_pidAlone/output_run1/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtracking/output_compAlex/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay/output_compAlex/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/Alex_results/v2_K0_Kch_pPb_16qSDD.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/Alex_results/v2_K0_Kch_pPb_16qSDD.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/Alex_results/v2_K0_Kch_pPb_16qSDD.root"
};

TString sHistName[iNumFiles] = {
  "hFlow2_Kaon_harm2_gap08_cent",
  "hFlow2_Kaon_harm2_gap08_cent",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hFlow2_K0s_harm2_gap08_mult",
  "hv2K0rad05_",
  "hv2K0rad5_",
  "hv2K_"
};

TString sLabel[iNumFiles] = {
  "Kch (Run1)",
  // "Kch",
  "Original (Kch)",
  "Original (K0s)",
  "Prel(NUAcor)",
  "Alex-like",
  "DCA_1",
  "Radius_5",
  "PID",
  "PID+Tracking",
  "PID+Tracking+Radius_5",
  "Alex(rad 0.5)",
  "Alex(rad 5)",
  "Alex(Kch)"
};

Short_t iReference = 2;
Bool_t bPlot[iNumFiles] = {
  0,
  1,
  1,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0
};

Color_t color[] = {
  kBlue,
  kBlack,
  kGreen,
  kGreen,// kMagenta,
  kOrange+1,
  kGreen+2,
  kRed,
  kBlue+2,
  kRed,
  kBlue-2,
  kRed,
  kGreen-2

  // 1,
  // 2,
  // 3,
  // 4,
  // 5,
  // 6,
  // 7,
  // 8,
  // 9,
  // 28
};


void ProcessFile(Short_t fileIndex, Short_t fileIndexBase, Short_t mult, TCanvas* canvas, TLegend* legend);
TH1D* LoadHistFromFile(const char* fileName, const char* histName);
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

void CompareK0sSelectionAlex()
{

  // ==============
  gSystem->mkdir(outDirectory.Data(),kTRUE);


  TCanvas* can = new TCanvas("can","can",400,800);

  TLine* lUnity = new TLine(dXlow,1.,dXupp,1.);
  lUnity->SetLineColor(kBlack);
  lUnity->SetLineStyle(kDashed);
  lUnity->SetLineWidth(1);

  TLegend* legend = new TLegend(0.13,0.88,0.4,0.48);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0,0);



  for(Short_t mult(0); mult < nBinsMult; mult++)
  {
    // TH1D* hKch = LoadHistFromFile("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA/plots/Processed.root",Form("hFlow2_Kaon_harm2_gap08_cent%d",mult));


    can->Divide(1,2);
    can->cd(1);
    TH1* frame_upper = gPad->DrawFrame(dXlow,dYlow,dXupp,dYupp,Form("Multiplicity %s%%; #it{p}_{T} (GeV/#it{c}); K0s #it{v}_{2}{2, #Delta#eta < |0.8|}",sMultLabel[mult].Data()));
    // hKch->Draw("same");
    can->cd(2);
    TH1* frame_lower = gPad->DrawFrame(dXlow,0.5,dXupp,1.5,Form("Ratio wrt %s; #it{p}_{T} (GeV/#it{c}); ratio",sLabel[iReference].Data()));

    for(Short_t file(0); file < iNumFiles; file++)
    {
      if(bPlot[file]) ProcessFile(file,iReference,mult,can,legend);
    }

    can->cd(1);
    legend->Draw("same");
    can->cd(2);
    lUnity->Draw("same");
    can->SaveAs(Form("%s/mult%d.pdf",outDirectory.Data(),mult));
    can->Clear();
    legend->Clear();
  }

  return;
}
// ===========================================================================
void ProcessFile(
    Short_t fileIndex,
    Short_t fileIndexBase,
    Short_t mult,
    TCanvas* canvas,
    TLegend* legend
  )
{
  if(!canvas) { printf("ERROR: Canvas does not exists!\n"); return; }
  if(!legend) { printf("ERROR: Legend does not exists!\n"); return; }

  const char* histName = Form("%s%d",sHistName[fileIndex].Data(),mult);
  const char* histNameBase = Form("%s%d",sHistName[fileIndexBase].Data(),mult);

  TH1D* hist = LoadHistFromFile(sInputFile[fileIndex].Data(),histName);
  TH1D* histBase = LoadHistFromFile(sInputFile[fileIndexBase].Data(),histNameBase);

  hist->SetLineColor(color[fileIndex]);
  hist->SetMarkerColor(color[fileIndex]);
  hist->SetMarkerStyle(kOpenCircle);
  legend->AddEntry(hist,sLabel[fileIndex].Data(),"pel");

  TH1D* hRatio = (TH1D*) hist->Clone(Form("%s_ratio",hist->GetName()));
  // hRatio->Divide(hist,histBase);

  TH1D* hRatio = DivideHistos(hist,histBase,0);

  TF1* fit = new TF1("fit","pol0",dXlow,dXupp);
  fit->SetLineColor(color[fileIndex]);
  hRatio->Fit("fit","NR");
  // hRatio->SetTitle("Ratio X / Alex");

  canvas->cd(1);
  hist->DrawCopy("same e1 p0");

  canvas->cd(2);
  if(fileIndex != fileIndexBase)
  hRatio->DrawCopy("same e1 p0");
  fit->Draw("same");

  return;
}
// ===========================================================================
TH1D* LoadHistFromFile(const char* fileName, const char* histName)
{
  if(!fileName) { printf("Filename '%s' does not exists!\n",fileName); return 0x0; }
  if(!histName) { printf("Histname '%s' does not exists!\n",fileName); return 0x0; }

  TFile* file = TFile::Open(fileName,"READ");
  if(!file) { printf("File '%s' does not exists!\n",fileName); return 0x0; }

  TH1D* hist = (TH1D*) file->Get(histName);
  if(!hist) { printf("Hist '%s' does not exists!\n",histName); return 0x0; }

  return hist;
}
// ===========================================================================
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
