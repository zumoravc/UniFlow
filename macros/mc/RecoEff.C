enum eMC {kGen = 0, kReco, kRecoTrue, kNum};

TString sHistName[eMC::kNum] = {
    "fh2MCPtEtaGen",
    "fh2MCPtEtaReco",
    "fh2MCPtEtaRecoTrue"
};

void AllRecoEff(
    Bool_t bOver = kFALSE,
    TString sTag = "UniFlow",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/running/grid/PbPb_MC/LHC18e1/purity_eff_v2/"
);

void RecoEff(
    TString sSpecies = "Charged",
    Bool_t bOverwrite = kFALSE,
    std::vector<Double_t> dPtBins = {},
    TString sTag = "UniFlow",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/running/grid/PbPb_MC/LHC18e1/purity_eff/"
)
{
    // Openning files
    TFile* fileOut = TFile::Open(Form("%s/McEffPurity.root",sPath.Data()),"UPDATE");
    if(!fileOut) { printf("fileOut not opened!\n"); return; }

    TFile* fileIn = TFile::Open(Form("%s/AnalysisResults.root",sPath.Data()),"READ");
    if(!fileIn) { printf("fileIn not opened!\n"); return; }

    fileIn->cd(sTag.Data());
    TList* listMC = (TList*) gDirectory->Get(Form("MC_%s",sTag.Data()));
    if(!listMC) { printf("ListMC not found!\n"); gDirectory->ls(); return; }

    // Loading histos (2D) & making 1D projections
    TH2D* h2[eMC::kNum] = {};
    TH1D* h1[eMC::kNum] = {};

    for(Int_t i(0); i < eMC::kNum; ++i) {
        TString sName = Form("%s%s",sHistName[i].Data(),sSpecies.Data());
        TH2D* hist2 = (TH2D*) listMC->FindObject(sName.Data());
        if(!hist2) { printf("'%s' not found!\n",sName.Data()); listMC->ls(); return; }
        hist2->Sumw2();
        h2[i] = hist2;


        TH1D* hist_projPt = (TH1D*) hist2->ProjectionX(Form("%s_pt",sName.Data()));
        if(!hist_projPt) { printf("ProjectionX (pt) of '%s' failed!\n",sName.Data()); return; }
        hist_projPt->Sumw2();
        h1[i] =  hist_projPt;

        if(dPtBins.size()) {
            h1[i] = (TH1D*) hist_projPt->Rebin(dPtBins.size()-1, Form("%s_pt_rebin",sName.Data()), dPtBins.data());
        }

    }

    // Making RecoEff
    TString sName_RecoEff_Pt = Form("%s_RecoEff_Pt",sSpecies.Data());
    TH1D* hRecoEff_Pt = (TH1D*) fileOut->Get(sName_RecoEff_Pt.Data());
    if(hRecoEff_Pt && !bOverwrite) {
        printf("RecoEff '%s' already exists in the file and overwrite is turn off!\n",sName_RecoEff_Pt.Data());
        // fileOut->ls();
    } else {
        hRecoEff_Pt = (TH1D*) h1[kReco]->Clone(sName_RecoEff_Pt.Data());
        hRecoEff_Pt->Divide(h1[kGen]);
        if(!hRecoEff_Pt) { printf("RecoEff_Pt division failed!\n"); return; }
    }

    TString sName_RecoEff_PtEta = Form("%s_RecoEff_PtEta",sSpecies.Data());
    TH2D* h2RecoEff_PtEta = (TH2D*) fileOut->Get(sName_RecoEff_PtEta.Data());
    if(h2RecoEff_PtEta && !bOverwrite) {
        printf("hRecoEff '%s' already exists in the file and overwrite is turn off!\n",sName_RecoEff_PtEta.Data());
        // fileOut->ls();
    } else {
        h2RecoEff_PtEta = (TH2D*) h2[kReco]->Clone(sName_RecoEff_PtEta.Data());
        h2RecoEff_PtEta->Divide(h2[kGen]);
        if(!h2RecoEff_PtEta) { printf("h2RecoEff_PtEta division failed!\n"); return; }
    }

    // Making Purity
    TString sName_Purity_Pt = Form("%s_Purity_Pt",sSpecies.Data());
    TH1D* hPurity_Pt = (TH1D*) fileOut->Get(sName_Purity_Pt.Data());
    if(hPurity_Pt && !bOverwrite) {
        printf("hPurity '%s' already exists in the file and overwrite is turn off!\n",sName_Purity_Pt.Data());
        // fileOut->ls();
    } else {
        hPurity_Pt = (TH1D*) h1[kRecoTrue]->Clone(sName_Purity_Pt.Data());
        hPurity_Pt->Divide(h1[kReco]);
        if(!hPurity_Pt) { printf("hPurity_Pt division failed!\n"); return; }
    }

    TString sName_Purity_PtEta = Form("%s_Purity_PtEta",sSpecies.Data());
    TH2D* h2Purity_PtEta = (TH2D*) fileOut->Get(sName_Purity_PtEta.Data());
    if(h2Purity_PtEta && !bOverwrite) {
        printf("Purity_PtEta '%s' already exists in the file and overwrite is turn off!\n",sName_Purity_PtEta.Data());
        // fileOut->ls();
    } else {
        h2Purity_PtEta = (TH2D*) h2[kRecoTrue]->Clone(sName_Purity_PtEta.Data());
        h2Purity_PtEta->Divide(h2[kReco]);
        if(!h2Purity_PtEta) { printf("h2Purity_PtEta division failed!\n"); return; }
    }

    // writing stuff to output file

    fileOut->cd();

    Int_t opt = 0;
    if(bOverwrite) { opt = TObject::kOverwrite; }


    h2[kGen]->Write(Form("%s_Gen_PtEta",sSpecies.Data()), opt);
    h2[kReco]->Write(Form("%s_Reco_PtEta",sSpecies.Data()), opt);
    h2[kRecoTrue]->Write(Form("%s_RecoTrue_PtEta",sSpecies.Data()), opt);
    h2RecoEff_PtEta->Write(Form("%s_RecoEff_PtEta",sSpecies.Data()), opt);
    h2Purity_PtEta->Write(Form("%s_Purity_PtEta",sSpecies.Data()), opt);

    h1[kGen]->Write(Form("%s_Gen_Pt",sSpecies.Data()), opt);
    h1[kReco]->Write(Form("%s_Reco_Pt",sSpecies.Data()), opt);
    h1[kRecoTrue]->Write(Form("%s_RecoTrue_Pt",sSpecies.Data()), opt);
    hRecoEff_Pt->Write(Form("%s_RecoEff_Pt",sSpecies.Data()), opt);
    hPurity_Pt->Write(Form("%s_Purity_Pt",sSpecies.Data()), opt);

    fileOut->Close();
    fileIn->Close();

    return;
}

void AllRecoEff(Bool_t bOver, TString sTag, TString sPath)
{
    std::vector<Double_t> bins = {
        0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
    RecoEff("Charged", bOver, bins, sTag, sPath);
    RecoEff("Pion", bOver, bins, sTag, sPath);
    RecoEff("Kaon", bOver, bins, sTag, sPath);
    RecoEff("Proton", bOver, bins, sTag, sPath);
    RecoEff("K0s", bOver, bins, sTag, sPath);
    RecoEff("Lambda", bOver, bins, sTag, sPath);
    RecoEff("Phi", bOver, bins, sTag, sPath);
}
