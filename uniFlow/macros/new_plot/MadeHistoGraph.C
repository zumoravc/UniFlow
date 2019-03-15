TGraphErrors* Transfer(TH1D* histo);
void SetEx(TGraphErrors* gae, Double_t Ex);

void MadeHistoGraph()
{
    TString sFilePath = "/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/All/";

    Int_t iNumCent = 7;

    TString sFileOutName;
    TString sHistoName;
    TString sHistoOutName;

    TString sCent[] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60"};

    // sHistoName = "K0s_<<3>>(4,-2,-2)_2sub(0)"; sFileOutName = "Mergedv422K0s"; sHistoOutName = "Mergedv422pT";
    // sHistoName = "K0s_<<3>>(5,-3,-2)_2sub(0)"; sFileOutName = "Mergedv523K0s"; sHistoOutName = "Mergedv523pT";
    // sHistoName = "K0s_<<3>>(6,-3,-3)_2sub(0)"; sFileOutName = "Mergedv633K0s"; sHistoOutName = "Mergedv633pT";

    // sHistoName = "Lambda_<<3>>(4,-2,-2)_2sub(0)"; sFileOutName = "Mergedv422Lambda"; sHistoOutName = "Mergedv422pT";
    // sHistoName = "Lambda_<<3>>(5,-3,-2)_2sub(0)"; sFileOutName = "Mergedv523Lambda"; sHistoOutName = "Mergedv523pT";
    // sHistoName = "Lambda_<<3>>(6,-3,-3)_2sub(0)"; sFileOutName = "Mergedv633Lambda"; sHistoOutName = "Mergedv633pT";

    // sHistoName = "Phi_<<3>>(4,-2,-2)_2sub(0)"; sFileOutName = "Mergedv422Phi"; sHistoOutName = "Mergedv422pT";
    // sHistoName = "Phi_<<3>>(5,-3,-2)_2sub(0)"; sFileOutName = "Mergedv523Phi"; sHistoOutName = "Mergedv523pT";
    sHistoName = "Phi_<<3>>(6,-3,-3)_2sub(0)"; sFileOutName = "Mergedv633Phi"; sHistoOutName = "Mergedv633pT";

    TFile* fileIn = TFile::Open(Form("%s/Processed.root",sFilePath.Data()),"READ");
    if(!fileIn) { printf("E: file input not open!\n"); return; }

    gSystem->mkdir(Form("%s/graphs/",sFilePath.Data()),kTRUE);
    TFile* fileOut = TFile::Open(Form("%s/graphs/%s.root",sFilePath.Data(),sFileOutName.Data()),"RECREATE");
    if(!fileOut)  { printf("E: file output not open!\n"); return; }

    for(Int_t iCent(0); iCent < iNumCent; ++iCent) {

        TH1D* hist = (TH1D*) fileIn->Get(Form("%s_mult%d",sHistoName.Data(),iCent));
        if(!hist) { printf("E: histo not found!\n"); fileIn->ls(); return; }

        TGraphErrors* graph = Transfer(hist);
        graph->SetName(Form("%s_%s_gap00",sHistoOutName.Data(),sCent[iCent].Data()));
        fileOut->cd();
        graph->Write();
    }

}

TGraphErrors* Transfer(TH1D* histo)
{
    if(!histo) { printf("E::Transfer : Input not found!\n"); return nullptr; }

    TGraphErrors* err = new TGraphErrors(histo);
    SetEx(err,0.0);

    // for(Int_t i(0); i < err->GetN(); ++i) {

    // }
    return err;

}

void SetEx(TGraphErrors* gae, Double_t Ex)
{
   Int_t np = gae->GetN();
   for (Int_t i=0; i<np; i++) {
      gae->SetPointError(i,Ex,gae->GetErrorY(i));
      // gae->SetPointError(i,Ex);
   }
}
