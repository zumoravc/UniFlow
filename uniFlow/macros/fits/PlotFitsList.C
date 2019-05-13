// Macro for plotting fits within fits.root file from ProcessUniFlow
// March 2018 (for NLF modes systematics run)

void PlotFitsList(
    TString sSpecies,
    Int_t iNumPt,
    Int_t iNumCent,
    TString sPath,
    std::vector<const char*> vTags,
    std::vector<const char*> vCorr
);

void PlotFitsList(
    TString sSpecies,
    Int_t iNumPt,
    Int_t iNumCent,
    std::vector<const char*> vPath,
    std::vector<const char*> vCorr
);

void PlotFitsList(
    TString sSpecies,
    Int_t iNumPt,
    Int_t iNumCent,
    TString sPath,
    std::vector<const char*> vCorr
);


void PlotFitsList(
    TString sSpecies = "K0s",
    Int_t iNumPt = 14,
    Int_t iNumCent = 6,
    TString sPath = "../results/cums/PbPb/syst_6815/K0s/CL1/",
    TString sCorr = "<<2>>(2,-2)"
)
{
    // TString sPath = "../results/nlf/output/K0s/", TString corrName = "<<2>>(2,-2)", TString sSpecies = "K0s", Int_t iNumPt = 8, Int_t iNumCent = 5)
    TMacro plot = TMacro(gSystem->ExpandPathName("~/Codes/ALICE/Flow/uniFlow/processUniFlow/PlotFits.C"));
    plot.Exec(Form("\"%s\",\"%s\",\"%s\",%d,%d",sPath.Data(), sCorr.Data(), sSpecies.Data(), iNumPt, iNumCent));
}


void PlotFitsList(TString sSpecies, Int_t iNumPt, Int_t iNumCent,std::vector<const char*> vPath, std::vector<const char*> vCorr)
{
    Int_t iNumPaths = vPath.size();
    Int_t iNumCorrs = vCorr.size();
    for(Int_t iP(0); iP < iNumPaths; ++iP) {
        TString sPath = vPath.at(iP);

        for(Int_t iC(0); iC < iNumCorrs; ++iC) {
            TString sCorr = vCorr.at(iC);
            PlotFitsList(sSpecies,iNumPt,iNumCent,sPath,sCorr);
        }
    }

    return;
}

void PlotFitsList(TString sSpecies,Int_t iNumPt,Int_t iNumCent,TString sPath,std::vector<const char*> vTags,std::vector<const char*> vCorr)
{
    Int_t iNumTags = vTags.size();
    Int_t iNumCorrs = vCorr.size();
    for(Int_t iT(0); iT < iNumTags; ++iT) {

        TString sPathThis = Form("%s/%s/",sPath.Data(), vTags.at(iT));

        for(Int_t iC(0); iC < iNumCorrs; ++iC) {
            TString sCorr = vCorr.at(iC);
            PlotFitsList(sSpecies,iNumPt,iNumCent,sPathThis,sCorr);
        }
    }
}

void PlotFitsList(
    TString sSpecies,
    Int_t iNumPt,
    Int_t iNumCent,
    TString sPath,
    std::vector<const char*> vCorr
)
{
    Int_t iNumCorrs = vCorr.size();
    for(Int_t iC(0); iC < iNumCorrs; ++iC) {
        TString sCorr = vCorr.at(iC);
        PlotFitsList(sSpecies,iNumPt,iNumCent,sPath,sCorr);
    }
}
