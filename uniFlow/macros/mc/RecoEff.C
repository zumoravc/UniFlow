void RecoEff()
{
    TString sPath = "~/Codes/Flow/uniFlow/class/";
    TString sTag = "UniFlow";
    TString sSpecies = "Refs";

    TFile* fileIn = TFile::Open(Form("%s/AnalysisResults.root",sPath.Data()),"READ");
    if(!fileIn) { printf("File not opened!\n"); return; }

    fileIn->cd(sTag.Data());

    TList* listMC = (TList*) gDirectory->Get(Form("MC_%s",sTag.Data()));
    if(!listMC) { printf("ListMC not found!\n"); gDirectory->ls(); return; }

    TH2D* h2Gen = (TH2D*) listMC->FindObject(Form("fh2MCPtEtaGen%s",sSpecies.Data()));
    if(!h2Gen) { printf("h2Gen not found!\n"); listMC->ls(); return; }

    TH2D* h2Reco = (TH2D*) listMC->FindObject(Form("fh2MCPtEtaReco%s",sSpecies.Data()));
    if(!h2Reco) { printf("h2Reco not found!\n"); listMC->ls(); return; }

    TH2D* h2True = (TH2D*) listMC->FindObject(Form("fh2MCPtEtaRecoTrue%s",sSpecies.Data()));
    if(!h2True) { printf("True not found!\n"); listMC->ls(); return; }

    TH2D* h2RecoEff = (TH2D*) h2Reco->Clone(Form("h2RecoEff%s",sSpecies.Data()));
    h2RecoEff->Divide(h2Gen);
    if(!h2RecoEff) { printf("h2RecoEff division failed!\n"); return; }

    h2RecoEff->DrawCopy("colz");

    delete h2RecoEff;
    return;
}
