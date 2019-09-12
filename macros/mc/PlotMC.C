enum PosLegend {kLegTopLeft = 1, kLegTopRight, kLegBotLeft, kLegBotRight};
TLegend* MakeLegend(PosLegend pos);

Color_t colors[] = {kRed, kGreen+2, kBlue};

void PlotMC(
    TString sHeader = "Pb-Pb (LHC18e1)",
    std::vector<TString> vecSpecies = {"Pion","Kaon","Proton"},
    TString sPathIn = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/running/grid/PbPb_MC/LHC18e1/purity_eff_v2/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/running/grid/PbPb_MC/LHC18e1/purity_eff_v2/plots/"
)
{
    TFile* fileMC = TFile::Open(Form("%s/McEffPurity.root",sPathIn.Data()),"READ");
    if(!fileMC) { printf("File not found!\n"); return; }

    TLegend* legPur = MakeLegend(kLegBotLeft);
    legPur->SetHeader(sHeader.Data());
    legPur->SetBorderSize(0);
    legPur->SetFillColorAlpha(0.0,0.0);

    TLegend* legRecoEff = MakeLegend(kLegTopRight);
    legRecoEff->SetHeader(sHeader.Data());
    legRecoEff->SetBorderSize(0);
    legRecoEff->SetFillColorAlpha(0.0,0.0);

    TLine* lUnity = new TLine();
    lUnity->SetLineColor(kBlack);
    lUnity->SetLineStyle(kDashed);
    lUnity->SetLineWidth(2);

    TCanvas* canPur = new TCanvas("canPur","Purity",600,600);
    TCanvas* canRecoEff = new TCanvas("canRecoEff","RecoEff",600,600);

    Int_t iNumS = vecSpecies.size();
    for(Int_t iS(0); iS < iNumS; ++iS) {

        TString sSpecies = vecSpecies.at(iS);

        TH1D* hPurity = (TH1D*) fileMC->Get(Form("%s_Purity_Pt",sSpecies.Data()));
        if(!hPurity) { printf("hPurity not found!\n"); fileMC->ls(); return; }

        TH1D* hRecoEff = (TH1D*) fileMC->Get(Form("%s_RecoEff_Pt",sSpecies.Data()));
        if(!hRecoEff) { printf("hRecoEff not found!\n"); fileMC->ls(); return; }

        if(iS == 0) {
            canPur->cd();
            gPad->SetLeftMargin(0.15);
            TH1* framePur = (TH1*) gPad->DrawFrame(0.0,0.6,hPurity->GetXaxis()->GetXmax(),1.05);
            framePur->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            framePur->GetYaxis()->SetTitle("Purity");
            lUnity->DrawLine(0.0,1.0,hPurity->GetXaxis()->GetXmax(),1.0);

            canRecoEff->cd();
            gPad->SetLeftMargin(0.15);
            TH1* frameRecoEff = (TH1*) gPad->DrawFrame(0.0,0.0,hPurity->GetXaxis()->GetXmax(),1.0);
            frameRecoEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            frameRecoEff->GetYaxis()->SetTitle("Efficiency");
            lUnity->DrawLine(0.0,1.0,hPurity->GetXaxis()->GetXmax(),1.0);

        }

        legPur->AddEntry(hPurity,sSpecies.Data(),"pl");
        legRecoEff->AddEntry(hPurity,sSpecies.Data(),"pl");

        Color_t col = colors[iS];

        canPur->cd();
        hPurity->SetMarkerColor(col);
        hPurity->SetMarkerStyle(kFullCircle);
        hPurity->SetMarkerSize(0.8);
        hPurity->SetLineColor(col);
        hPurity->DrawCopy("same");

        canRecoEff->cd();
        hRecoEff->SetMarkerColor(col);
        hRecoEff->SetMarkerStyle(kFullCircle);
        hRecoEff->SetMarkerSize(0.8);
        hRecoEff->SetLineColor(col);
        hRecoEff->DrawCopy("same");
    }

    canPur->cd();
    legPur->Draw();

    canRecoEff->cd();
    legRecoEff->Draw();

    gSystem->mkdir(sPathOut.Data(),1);

    canPur->SaveAs(Form("%s/Purity_%s.pdf",sPathOut.Data(),sHeader.Data()),"pdf");
    canRecoEff->SaveAs(Form("%s/RecoEff_%s.pdf",sPathOut.Data(),sHeader.Data()),"pdf");

}

TLegend* MakeLegend(PosLegend pos)
{
    switch(pos) {
        case kLegTopLeft: return new TLegend(0.26,0.6,0.5,0.88); break;
        case kLegTopRight: return new TLegend(0.62,0.6,0.88,0.88); break;
        case kLegBotLeft: return new TLegend(0.26,0.12,0.5,0.4); break;
        case kLegBotRight: return new TLegend(0.62,0.12,0.88,0.4); break;
        default: return nullptr;
    }
}
