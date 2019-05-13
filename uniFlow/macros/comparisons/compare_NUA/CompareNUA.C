#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"

Color_t colors[] = {kBlue, kRed, kGreen};

TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

void CompareNUA()
{
  const char* dirs[] = {
    "preliminary",
    "preliminary_NUA-runbyrun",
  };

  const char* species = "Lambda";

  Int_t numDirs = sizeof(dirs)/sizeof(dirs[0]);

  TList* listPions[4] = {
    new TList(),
    new TList(),
    new TList(),
    new TList()
  };

  TCanvas* canPions = new TCanvas("canPions","canPions",1600,800);
  canPions->Divide(4,2);

  TLegend* leg = new TLegend(0.1,0.5,0.6,0.89);

  for(Int_t dir(0); dir < numDirs; dir++)
  {
    const char* path = Form("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s-rerun/%s/output/Processed.root",dirs[dir]);

    TFile* file = TFile::Open(path,"READ");
    if(!file) return;
    file->ls();

    TH1D* hPion[4] = {
      (TH1D*) file->Get(Form("hFlow2_%s_harm2_gap08_mult0",species)),
      (TH1D*) file->Get(Form("hFlow2_%s_harm2_gap08_mult1",species)),
      (TH1D*) file->Get(Form("hFlow2_%s_harm2_gap08_mult2",species)),
      (TH1D*) file->Get(Form("hFlow2_%s_harm2_gap08_mult3",species))
    };

    leg->AddEntry(hPion[0],dirs[dir],"pel");

    for(Int_t hist(0); hist < 4; hist++)
    {
      listPions[hist]->Add(hPion[hist]);

      canPions->cd(hist+1);
      hPion[hist]->SetMinimum(0.);
      hPion[hist]->SetMaximum(1.);
      hPion[hist]->SetLineColor(colors[dir]);
      hPion[hist]->SetMarkerColor(colors[dir]);
      hPion[hist]->SetStats(0);
      hPion[hist]->Draw("same");
    }
  }

  for(Int_t hist(0); hist < 4; hist++)
  {
    canPions->cd(5+hist);
    TH1D* hRatio = DivideHistos((TH1D*)listPions[hist]->At(0),(TH1D*)listPions[hist]->At(1),kTRUE);
    hRatio->SetMinimum(0.8);
    hRatio->SetMaximum(1.2);
    hRatio->Draw();
  }

  canPions->cd(1);
  leg->Draw();



}

// ===========================================================================
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
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
