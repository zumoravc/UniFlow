enum PosLegend {kLegTopLeft = 1, kLegTopRight, kLegBotLeft, kLegBotRight};
TLegend* MakeLegend(PosLegend pos);

void PlotFluctPubl(
    TString sInDir = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/6815/V0M/gap_all/",
    TString sFileName = "FlowFluctCor",
    TString sOutDir = ""
)
{

    if(sOutDir.Length() < 1) { sOutDir = Form("%s/plots_%s/",sInDir.Data(),sFileName.Data()); }

    // published
    Bool_t bPubl = 1;
    TString sInFilePubl = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/published/HEPData-ins1116150-v1-root.root";
    TString sTable = "Table 7";

    Double_t xMin = 0.0; Double_t xMax = 8.0;
    Double_t yMin = 0.0; Double_t yMax = 0.75;
    PosLegend pos = kLegBotLeft;

    // TString sLabel = "<v_{2}>"; TString sHistName = "mean"; yMin = 0.0; yMax = 0.25;  pos = kLegTopRight; bPubl = 0;
    TString sLabel = "#sigma(v_{2})"; TString sHistName = "std"; yMin = 0.0; yMax = 0.12;  pos = kLegTopRight; bPubl = 0;
    // TString sLabel = "F(v_{2})"; TString sHistName = "rel"; yMin = 0.2; yMax = 0.65;  pos = kLegTopLeft; bPubl = 1;

    std::vector<TString> sSpecies = {};
    std::vector<TString> sSpesLabels = {};
    std::vector<Color_t> colors = {};
    std::vector<Int_t> markers = {};
    std::vector<Double_t> markerSizes = {};

    sSpecies.push_back("Charged"); sSpesLabels.push_back("h^{#pm}");
    // sSpecies.push_back("Pion"); sSpesLabels.push_back("#pi^{#pm}"); bPubl = 0;
    // sSpecies.push_back("Kaon"); sSpesLabels.push_back("K^{#pm}"); bPubl = 0; pos = kLegTopRight;
    // sSpecies.push_back("K0s"); sSpesLabels.push_back("K^{0}_{S}");bPubl = 0;
    // sSpecies.push_back("Proton"); sSpesLabels.push_back( "p(#bar{p})");bPubl = 0; pos = kLegTopRight;
    // sSpecies.push_back("Lambda"); sSpesLabels.push_back("#Lambda(#bar{#Lambda})"); bPubl = 0;
    // sSpecies.push_back("Phi"); sSpesLabels.push_back("#phi"); bPubl = 0;

    colors.push_back(kBlack);
    colors.push_back(kRed+1);
    colors.push_back(kGreen-1);
    // colors.push_back(kGreen-1);
    colors.push_back(kBlue);
    colors.push_back(kOrange+1);
    colors.push_back(kMagenta);

    markers.push_back(kFullCircle);
    markers.push_back(kFullCircle);
    markers.push_back(kFullTriangleDown);
    markers.push_back(kFullTriangleUp);
    markers.push_back(kFullSquare);
    markers.push_back(kFullDiamond);
    markers.push_back(kFullCross);

     markerSizes.push_back(1.0);
     markerSizes.push_back(1.0);
     markerSizes.push_back(1.1);
     markerSizes.push_back(1.1);
     markerSizes.push_back(0.8);
     markerSizes.push_back(1.4);
     markerSizes.push_back(1.2);

    TString sCent[] = {
        // "0-5% V0M",
        "5-10% V0M",
        // "0-10% V0M",
        "10-20% V0M",
        "20-30% V0M",
        "30-40% V0M",
        "40-50% V0M",
        "50-60% V0M",
        // "60-70% V0M",
        // "70-80% V0M",
        // "80-90% V0M",
        // "90-100% V0M"
    };


    TFile* fileIn = TFile::Open(Form("%s/%s.root",sInDir.Data(),sFileName.Data()),"READ");

    TFile* filePubl = TFile::Open(Form("%s",sInFilePubl.Data()),"READ");

    gSystem->Exec(Form("mkdir -p %s", sOutDir.Data()));

    TCanvas* canCent = new TCanvas("canCent","canCent",900,600);
    canCent->cd();
    TH1* frameCent = (TH1*) gPad->DrawFrame(xMin, yMin, xMax, yMax);
    frameCent->SetTitle("; #it{p}_{T} (GeV/#it{c});");
    TLegend* legCent = MakeLegend(pos);
    legCent->SetBorderSize(0);
    legCent->SetFillColorAlpha(0.0,0.0);

    Int_t iNumMult = sizeof(sCent) / sizeof(sCent[0]);
    for(Int_t iMult(0); iMult < iNumMult; ++iMult) {

        TString sName = Form("%s_%s_cent%d",sSpecies[0].Data(), sHistName.Data(), iMult);

        TH1D* histCL = (TH1D*) fileIn->Get(sName.Data());

        canCent->cd();
        histCL->SetLineColor(colors[iMult]);
        histCL->SetMarkerColor(colors[iMult]);
        // histCL->SetMarkerStyle(markers[iMult]);
        // histCL->SetMarkerSize(markerSizes[iMult]);
        histCL->SetMarkerStyle(kFullCircle);
        histCL->SetMarkerSize(1.0);
        histCL->DrawCopy("hist ape1x0 same");

        // published
        if(bPubl && iMult < 5) {
            filePubl->cd(sTable.Data());
            TString sNamePubl = Form("Graph1D_y%d",iMult+2);
            TGraphAsymmErrors* gr = (TGraphAsymmErrors*) gDirectory->Get(sNamePubl.Data());
            if(!gr) { printf("Error- not found publ %s\n",sNamePubl.Data()); gDirectory->ls(); return; }

            gr->SetFillColorAlpha(colors[iMult],0.2);
            gr->SetLineColor(colors[iMult]);
            gr->SetLineWidth(1);
            gr->SetMarkerColor(colors[iMult]);
            gr->SetMarkerStyle(kOpenSquare);
            gr->SetMarkerSize(1.0);
            gr->Draw("p5 same");
            // gr->Draw("|| same");
            if(iMult == 0) {
                legCent->AddEntry(gr,"ALICE PLB 719 (2013)","p");
            }
        }
        legCent->AddEntry(histCL,sCent[iMult].Data(),"p");

    }

    canCent->cd();
    legCent->SetHeader(Form("%s %s",sSpesLabels[0].Data(),sLabel.Data()));
    legCent->Draw();

    if(bPubl) {
        canCent->SaveAs(Form("%s/Publ_%s_%s.pdf",sOutDir.Data(),sSpecies[0].Data(),sLabel.Data()),"pdf");
    } else {
        canCent->SaveAs(Form("%s/%s_%s.pdf",sOutDir.Data(),sSpecies[0].Data(),sLabel.Data()),"pdf");
    }


    return;
}



TLegend* MakeLegend(PosLegend pos)
{
    switch(pos) {
        case kLegTopLeft: return new TLegend(0.14,0.65,0.5,0.88); break;
        case kLegTopRight: return new TLegend(0.68,0.65,0.88,0.88); break;
        case kLegBotLeft: return new TLegend(0.14,0.12,0.5,0.4); break;
        case kLegBotRight: return new TLegend(0.68,0.12,0.88,0.4); break;
        default: return nullptr;
    }
}
