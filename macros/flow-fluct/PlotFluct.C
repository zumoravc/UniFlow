enum PosLegend {kLegTopLeft = 1, kLegTopRight, kLegBotLeft, kLegBotRight};
TLegend* MakeLegend(PosLegend pos);

void PlotFluctAll(TString sInDir = "~/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/MC_2210/fluct_2050/", TString sFileName = "FlowFluctCor");

void PlotFluct(
    TString sInDir = "~/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/MC_2210/species_direct/",
    TString sFileName = "FlowFluct",
    TString sLabel = "F(v_{2})",
    TString sHistName = "rel",
    Double_t yMin = 0.2,
    Double_t yMax = 0.650,
    PosLegend pos = kLegTopRight,
    Double_t xMin = 0.0,
    Double_t xMax = 4.0
)
{
    // TString sInDir = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/6815/V0M/gap_all/";
    TString sOutDir = Form("%s/plots_%s/",sInDir.Data(),sFileName.Data());

    // TString sLabel = "<v_{2}>"; TString sHistName = "mean"; yMin = 0.0; yMax = 0.3;  pos = kLegTopRight;
    // TString sLabel = "#sigma(v_{2})"; TString sHistName = "std"; yMin = 0.0; yMax = 0.15;  pos = kLegTopRight;
    // TString sLabel = "F(v_{2})"; TString sHistName = "rel"; yMin = 0.2; yMax = 0.650;  pos = kLegTopRight;

    TString sCent[] = {
        // "0-5% V0M",
        // "5-10% V0M",
        // "0-10% V0M",
        // "10-20% V0M",
        // "20-30% V0M",
        // "30-40% V0M",
        // "40-50% V0M"
        "20-50% V0M",
        // "50-60% V0M",
        // "60-70% V0M",
        // "70-80% V0M",
        // "80-90% V0M",
        // "90-100% V0M"
    };

    std::vector<TString> sSpecies = {};
    std::vector<TString> sSpesLabels = {};
    std::vector<Color_t> colors = {};
    std::vector<Int_t> markers = {};
    std::vector<Double_t> markerSizes = {};

    sSpecies.push_back("Charged"); sSpesLabels.push_back("h^{#pm}");
    sSpecies.push_back("Pion"); sSpesLabels.push_back("#pi^{#pm}");
    sSpecies.push_back("Kaon"); sSpesLabels.push_back("K^{#pm}");
    sSpecies.push_back("Proton"); sSpesLabels.push_back("p(#bar{p})");
    // sSpecies.push_back("K0s"); sSpesLabels.push_back("K^{0}_{S}");
    // sSpecies.push_back("Lambda"); sSpesLabels.push_back("#Lambda(#bar{#Lambda})");
    // sSpecies.push_back("Phi"); sSpesLabels.push_back("#phi");



    colors.push_back(kGray+2); markers.push_back(kFullCircle); markerSizes.push_back(1.0);
    colors.push_back(kRed+1); markers.push_back(kFullCircle); markerSizes.push_back(1.0);
    colors.push_back(kGreen+3); markers.push_back(kFullTriangleDown); markerSizes.push_back(1.1);
    colors.push_back(kBlue); markers.push_back(kFullSquare); markerSizes.push_back(0.8);
    colors.push_back(kGreen+1); markers.push_back(kFullTriangleUp); markerSizes.push_back(1.1);
    colors.push_back(kOrange+1); markers.push_back(kFullDiamond); markerSizes.push_back(1.4);
    colors.push_back(kMagenta); markers.push_back(kFullCross); markerSizes.push_back(1.2);




    TFile* fileIn = TFile::Open(Form("%s/%s.root",sInDir.Data(),sFileName.Data()),"READ");

    // TFile* fileIn = TFile::Open(Form("%s/FlowFluct.root",sInDir.Data()),"READ");
    // fileIn->ls();

    // TFile* fileInPhi = TFile::Open(Form("%s/Processed.root",sInDir.Data()),"READ");
    // TFile* fileInPhi = TFile::Open(Form("%s/ProcessedPhi.root",sInDir.Data()),"READ");
    // fileIn->ls();

    gSystem->Exec(Form("mkdir %s",sOutDir.Data()));

    TCanvas* canCent = new TCanvas("canCent","canCent",600,600);
    canCent->cd();
    TH1* frameCent = (TH1*) gPad->DrawFrame(xMin, yMin, xMax, yMax);
    frameCent->SetTitle("; #it{p}_{T} (GeV/#it{c});");
    TLegend* legCent = MakeLegend(pos);
    legCent->SetBorderSize(0);
    legCent->SetFillColorAlpha(0.0,0.0);

    Int_t iNumMult = sizeof(sCent) / sizeof(sCent[0]);
    for(Int_t iMult(0); iMult < iNumMult; ++iMult) {

        TCanvas* can = new TCanvas("can","can",900,600);
        can->cd();
        TH1* frame = (TH1*) gPad->DrawFrame(xMin, yMin, xMax, yMax);
        frame->SetTitle("; #it{p}_{T} (GeV/#it{c});");

        TLegend* leg = MakeLegend(pos);
        // TLegend* leg = new TLegend(0.14,0.12,0.5,0.4);
        // TLegend* leg = new TLegend(0.14,0.65,0.5,0.88);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(0.0,0.0);

        Int_t iNumPlots = sSpecies.size();
        for(Int_t i(0); i < iNumPlots; ++i) {

            // if(sSpecies[i].Contains("Charged")) {

            // }

            Bool_t bPhi = (sSpecies[i].Contains("Phi"));

            TString sName = Form("%s_%s_cent%d",sSpecies[i].Data(), sHistName.Data(), iMult);

            TH1D* hist = (TH1D*) fileIn->Get(sName.Data());

            if(bPhi) {
                // if(iMult == 0 || iMult == iNumMult-1) continue;
                // hist = (TH1D*) fileInPhi->Get(sName.Data());
            }

            if(!hist) { printf("ERROR: histo '%s' not found!\n",sName.Data()); fileIn->ls(); return; }
            can->cd();
            hist->SetLineColor(colors[i]);
            hist->SetMarkerColor(colors[i]);
            hist->SetMarkerStyle(markers[i]);
            hist->SetMarkerSize(markerSizes[i]);
            hist->DrawCopy("hist ape1x0 same");

            leg->AddEntry(hist,sSpesLabels[i].Data(),"p");

            if(sSpecies[i].Contains("Charged")) {
                canCent->cd();
                TH1D* histCL = (TH1D*) hist->Clone("histCL");
                histCL->SetLineColor(colors[iMult]);
                histCL->SetMarkerColor(colors[iMult]);
                histCL->SetMarkerStyle(markers[iMult]);
                histCL->SetMarkerSize(markerSizes[iMult]);
                histCL->DrawCopy("hist ape1x0 same");
                legCent->AddEntry(histCL,sCent[iMult].Data(),"p");
            }
        }

        leg->SetHeader(Form("%s (%s)",sLabel.Data(),sCent[iMult].Data()));
        leg->Draw();

        TLine* line_0 = new TLine(xMin, 0.0, xMax, 0.0);
        line_0->SetLineColor(kGray+2);
        line_0->SetLineStyle(7);
        line_0->Draw();

        can->SaveAs(Form("%s/All_%s_cent%d.pdf",sOutDir.Data(),sLabel.Data(),iMult),"pdf");

        delete leg;
        delete frame;
        delete can;



    }

    canCent->cd();
    legCent->SetHeader(Form("%s",sLabel.Data()));
    legCent->Draw();
    canCent->SaveAs(Form("%s/Charged_%s.pdf",sOutDir.Data(),sLabel.Data()),"pdf");

    return;
}


void PlotFluctAll(TString sInDir, TString sFileName)
{
    PlotFluct(sInDir, sFileName, "<v_{2}>","mean",0.0, 0.3, kLegTopRight);
    PlotFluct(sInDir, sFileName, "#sigma(v_{2})","std", 0.0, 0.15, kLegTopRight);
    PlotFluct(sInDir, sFileName, "F(v_{2})", "rel", 0.2, 0.650, kLegTopRight,0.0,4.0);
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
