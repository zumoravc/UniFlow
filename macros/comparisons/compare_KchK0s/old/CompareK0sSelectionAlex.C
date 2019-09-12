/*
Macro for comparison of results of my fitting procedure with the one provided by Alex
and implemented in method ExtractFlowK0sAlex()
- commit ac92ad6a67511bac6652fbc93d409e6f0d175f35

Author: Vojtech Pacik, NBI (vojtech.pacik@cern.ch)
*/

const Short_t iNumFiles = 6;
const Short_t nBinsMult = 4;
TString sMultLabel[nBinsMult] =
{
  "0-20",
  "20-40",
  "40-60",
  "60-100"
};

TString outDirectory = "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/test";
TString sHistName = "hFlow2_K0s_harm2_gap08_mult";

TString sInputFile[iNumFiles] = {
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/prel_noNUA/plots/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_Alex/plots/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_AN/plots/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_pidAlone/plots/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_decay5/plots/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0sComp4_PIDtrackingDecay/plots/Processed.root"
};

TString sLabel[iNumFiles] = {
  "Prel",
  "Alex",
  "AN",
  "PID",
  "Decay",
  "PID+Tracking+Decay"
};

Color_t color[iNumFiles] = {
  kBlue,
  kRed,
  kGreen+2,
  kOrange,
  38,
  kMagenta+2
};

void ProcessFile(TString fileName, TString label, Color_t color, TString histName, TString fileNameBase, TString histNameBase, TCanvas* canvas, TLegend* legend);

void CompareK0sSelectionAlex()
{

  // ==============
  gSystem->mkdir(outDirectory.Data(),kTRUE);

  TCanvas* can = new TCanvas("can","can");

  TLine* lUnity = new TLine(0.,1.,20.,1.);
  lUnity->SetLineColor(kBlack);
  lUnity->SetLineStyle(kDashed);
  lUnity->SetLineWidth(1);

  TLegend* legend = new TLegend(0.13,0.88,0.3,0.68);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0,0);

  for(Short_t mult(0); mult < nBinsMult; mult++)
  {
    can->Divide(1,2);
    can->cd(1);
    TH1* frame_upper = gPad->DrawFrame(0.,0.,20.,1.,Form("Multiplicity %s%%; #it{p}_{T} (GeV/#it{c}); #it{v}_{2}{2}",sMultLabel[mult].Data()));
    can->cd(2);
    TH1* frame_lower = gPad->DrawFrame(0.,0.,20.,2.,Form("; #it{p}_{T} (GeV/#it{c}); ratio X /%s",sLabel[1].Data()));

    for(Short_t file(0); file < iNumFiles; file++)
    {
      ProcessFile(sInputFile[file],sLabel[file],color[file],sHistName.Format("%s%d",sHistName.Data(),mult),sInputFile[1],sHistName.Format("%s0",sHistName.Data()),can,legend);
    }

    can->cd(1);
    legend->Draw();
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
    TString fileName,
    TString label,
    Color_t color,
    TString histName,
    TString fileNameBase,
    TString histNameBase,
    TCanvas* canvas,
    TLegend* legend
  )
{
  if(!canvas) { printf("ERROR: Canvas does not exists!\n"); return; }
  if(!legend) { printf("ERROR: Legend does not exists!\n"); return; }

  TFile* fileBase = TFile::Open(fileNameBase.Data(),"READ");
  if(!fileBase) { printf("ERROR: File (base) '%s' does not exists!\n",fileNameBase.Data()); return; }

  TH1D* histBase = (TH1D*) fileBase->Get(histNameBase.Data());
  if(!histBase) { printf("ERROR: Hist (base) '%s' does not exists!\n",histNameBase.Data()); return; }

  TFile* file = TFile::Open(fileName.Data(),"READ");
  if(!file) { printf("ERROR: File '%s' does not exists!\n",fileName.Data()); return; }

  TH1D* hist = (TH1D*) file->Get(histName.Data());
  if(!hist) { printf("ERROR: Hist '%s' does not exists!\n",histName.Data()); return; }

  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  legend->AddEntry(hist,label.Data(),"pel");

  TH1D* hRatio = (TH1D*) hist->Clone(Form("%s_ratio",hist->GetName()));
  hRatio->Divide(hist,histBase);
  // hRatio->SetTitle("Ratio X / Alex");

  canvas->cd(1);
  hist->DrawCopy("same p0 e1");

  canvas->cd(2);
  hRatio->DrawCopy("same e1 p0");

  return;
}
