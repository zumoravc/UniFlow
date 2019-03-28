TGraphErrors* Transfer(TH1D* histo, TF1* syst);
void SetEx(TGraphErrors* gae, Double_t Ex);

void MakeSyst(TH1* hist, TF1* syst);

void MadeHistoGraphSyst()
{
    TString sFilePath = "/mnt/CodesALICE/Flow/uniFlow/results/nlf/output/Phi/";
    TString sFilePathSyst = "/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/Phi/final/values/";

    Int_t iNumCent = 7;

    TString sFileInSyst;
    TString sFileOutName;
    TString sHistoName;
    TString sHistoOutName;

    TString sCent[] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60"};

    // sHistoName = "K0s_<<3>>(4,-2,-2)_2sub(0)"; sFileOutName = "MergedSystv422K0s"; sHistoOutName = "MergedSystv422pT"; sFileInSyst = "K0s_v422.root";
    // sHistoName = "K0s_<<3>>(5,-3,-2)_2sub(0)"; sFileOutName = "MergedSystv523K0s"; sHistoOutName = "MergedSystv523pT"; sFileInSyst = "K0s_v523.root";
    // sHistoName = "K0s_<<3>>(6,-3,-3)_2sub(0)"; sFileOutName = "MergedSystv633K0s"; sHistoOutName = "MergedSystv633pT"; sFileInSyst = "K0s_v633.root";

    // sHistoName = "Lambda_<<3>>(4,-2,-2)_2sub(0)"; sFileOutName = "MergedSystv422Lambda"; sHistoOutName = "MergedSystv422pT"; sFileInSyst = "Lambda_v422.root";
    // sHistoName = "Lambda_<<3>>(5,-3,-2)_2sub(0)"; sFileOutName = "MergedSystv523Lambda"; sHistoOutName = "MergedSystv523pT"; sFileInSyst = "Lambda_v523.root";
    // sHistoName = "Lambda_<<3>>(6,-3,-3)_2sub(0)"; sFileOutName = "MergedSystv633Lambda"; sHistoOutName = "MergedSystv633pT"; sFileInSyst = "Lambda_v633.root";

    sHistoName = "Phi_<<3>>(4,-2,-2)_2sub(0)"; sFileOutName = "MergedSystv422Phi"; sHistoOutName = "MergedSystv422pT";  sFileInSyst = "Phi_v422.root";
    // sHistoName = "Phi_<<3>>(5,-3,-2)_2sub(0)"; sFileOutName = "Mergedv523Phi"; sHistoOutName = "Mergedv523pT";
    // sHistoName = "Phi_<<3>>(6,-3,-3)_2sub(0)"; sFileOutName = "Mergedv633Phi"; sHistoOutName = "Mergedv633pT";

    TFile* fileIn = TFile::Open(Form("%s/Processed.root",sFilePath.Data()),"READ");
    if(!fileIn) { printf("E: file input not open!\n"); return; }

    TFile* fileInSyst = TFile::Open(Form("%s/%s",sFilePathSyst.Data(),sFileInSyst.Data()),"READ");
    if(!fileInSyst) { printf("E: file syst input not open!\n"); return; }

    gSystem->mkdir(Form("%s/graphs/",sFilePath.Data()),kTRUE);
    TFile* fileOut = TFile::Open(Form("%s/graphs/%s.root",sFilePath.Data(),sFileOutName.Data()),"RECREATE");
    if(!fileOut)  { printf("E: file output not open!\n"); return; }

    for(Int_t iCent(0); iCent < iNumCent; ++iCent) {

        TH1D* hist = (TH1D*) fileIn->Get(Form("%s_mult%d",sHistoName.Data(),iCent));
        if(!hist) { printf("E: histo not found!\n"); fileIn->ls(); return; }

        TList* list = (TList*) fileInSyst->Get(Form("mult%d",iCent));
        if(!list){ printf("E: fit list not found!\n"); return; }
        TF1* fit = (TF1*) list->FindObject("total");
        if(!fit){ printf("E: fit funcionf not found!\n"); return; }

        TGraphErrors* graph = Transfer(hist,fit);
        graph->SetName(Form("%s_%s_gap00",sHistoOutName.Data(),sCent[iCent].Data()));
        fileOut->cd();
        graph->Write();
    }

}

TGraphErrors* Transfer(TH1D* histo, TF1* syst)
{
    if(!histo) { printf("E::Transfer : Input not found!\n"); return nullptr; }

    Double_t dSystRel = syst->GetParameter(0);

    Double_t dX = 0.0;
    Double_t dY = 0.0;

    TGraphErrors* err = new TGraphErrors(histo);
    // SetEx(err,0.05);
    Int_t np = err->GetN();
    for (Int_t i=0; i<np; i++) {
        err->GetPoint(i,dX,dY);
        err->SetPointError(i,0.05, dSystRel*dY);
    }

    return err;

}

void SetEx(TGraphErrors* gae, Double_t Ex, Double_t)
{
   Int_t np = gae->GetN();
   for (Int_t i=0; i<np; i++) {
      gae->SetPointError(i,Ex,gae->GetErrorY(i));
      // gae->SetPointError(i,Ex);
   }
}

void MakeSyst(TH1* hist, TF1* syst)
{



}
