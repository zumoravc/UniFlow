void PlotSingle(
        const char* fileName = "./default",
        const char* histoName = "Lambda_<<3>>(4,-2,-2)_2sub(0)"
);

void PlotAll(const char* fileName = "./default")
{
    PlotSingle(fileName,"Lambda_<<3>>(4,-2,-2)_2sub(0)");
    PlotSingle(fileName,"Lambda_<<3>>(5,-3,-2)_2sub(0)");
    PlotSingle(fileName,"Lambda_<<3>>(6,-3,-3)_2sub(0)");
}

void PlotSingle( const char* fileName, const char* histoName)
{
    // const char* fileName = "~/Codes/Flow/uniFlow/results/trains/CF_PbPb/6527_20190218-2140/output_mixed/gap00/Lambda/Processed.root";
    // const char* histoName = "Lambda_<<3>>(4,-2,-2)_2sub(0)";

    Int_t iNumCent = 7;

    TFile* file = TFile::Open(Form("%s/Processed.root",fileName),"READ");
    if(!file) return;

    TCanvas* can = new TCanvas(Form("can_%s",histoName),Form("can_%s",histoName),1200,700);
    can->Divide(4,2);

    for(Int_t cent(0); cent < iNumCent; ++cent) {
        TH1D* hist = (TH1D*) gDirectory->Get(Form("%s_mult%d",histoName,cent));
        if(!hist) gDirectory->ls();

        hist->SetStats(0);
        can->cd(cent+1);
        hist->Draw("e1");
    }
}
