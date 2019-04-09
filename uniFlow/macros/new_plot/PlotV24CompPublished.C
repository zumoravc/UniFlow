void PlotV24CompPublished()
{
    // Comparison with https://www.hepdata.net/record/ins1666817
    // JHEP 1807 (2018) 103, 2018

    TString sOutputPath = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/plots/comp_published/";
    TString sPathPubl = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/published/HEPData-ins1666817-v1-root.root";

    TString sPathFile, sTable, sHeader, sHistoName;
    TString sHistLabel = "this";

    // sTable = "Table 1"; sHeader = "v_{2}{2,|#Delta#eta|>1.0}"; sHistLabel = "this (|#Delta#eta| > 0.8)"; sHistoName = "Refs_hFlow2_harm2_gap08"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_norebin_refs/gap_all/Processed.root";
    // sTable = "Table 2"; sHeader = "v_{2}{4}"; sHistoName = "Refs_hFlow4_harm2_gap-10"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_norebin_refs/gap_all/Processed.root";
    // sTable = "Table 31"; sHeader = "v_{2}{4} (5-10%)"; sHistoName = "Charged_hFlow4_harm2_gap-10_cent1"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_publ/gap_all/Processed.root";
    // sTable = "Table 41"; sHeader = "v_{2}{4} (10-20%)"; sHistoName = "Charged_hFlow4_harm2_gap-10_cent2"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_publ/gap_all/Processed.root";
    // sTable = "Table 51"; sHeader = "v_{2}{4} (20-30%)"; sHistoName = "Charged_hFlow4_harm2_gap-10_cent3"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_publ/gap_all/Processed.root";
    // sTable = "Table 61"; sHeader = "v_{2}{4} (30-40%)"; sHistoName = "Charged_hFlow4_harm2_gap-10_cent4"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_publ/gap_all/Processed.root";
    // sTable = "Table 71"; sHeader = "v_{2}{4} (40-50%)"; sHistoName = "Charged_hFlow4_harm2_gap-10_cent5"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_publ/gap_all/Processed.root";
    sTable = "Table 81"; sHeader = "v_{2}{4} (50-60%)"; sHistoName = "Charged_hFlow4_harm2_gap-10_cent6"; sPathFile = "/mnt/CodesALICE/Flow/uniFlow/results/cums/PbPb/6710/V0M_publ/gap_all/Processed.root";




    TFile* fileIn = TFile::Open(sPathFile.Data(),"READ");
    if(!fileIn) { printf("ERROR: fileIn not found!\n"); return; }

    TFile* filePubl = TFile::Open(sPathPubl.Data(), "READ");
    if(!filePubl) { printf("ERROR: filePubl not found!\n"); return; }

    TH1D* histo = (TH1D*) fileIn->Get(sHistoName.Data());
    if(!histo) { printf("ERROR: histo not found!\n"); fileIn->ls(); return; }

    filePubl->cd(sTable.Data());
    // gDirectory->ls();

    // TH1F* histo_publ = (TH1F*) gDirectory->Get("Hist1D_y1");
    TGraphAsymmErrors* histo_publ = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1");
    if(!histo_publ) { printf("ERROR: histo_publ nto found!\n"); return; }

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.0,0.0);
    leg->SetHeader(sHeader.Data());
    leg->AddEntry(histo,sHistLabel.Data(),"pel");
    leg->AddEntry(histo_publ,"published","pf");

    TCanvas* can = new TCanvas("can","can",800,800);
    can->cd();

    Double_t xMin = 0.0;
    Double_t xMax = histo->GetXaxis()->GetXmax();
    Double_t yMin = 0.0;
    Double_t yMax = 0.3;

    TH1* frame = (TH1*) gPad->DrawFrame(xMin,yMin,xMax,yMax);

    histo_publ->SetFillColor(kGreen-7);
    histo_publ->SetLineColor(kGreen+1);
    histo_publ->SetMarkerColor(kGreen+1);
    histo_publ->SetMarkerStyle(7);
    histo_publ->Draw("same pe2");

    // histo->SetMinimum(0);
    // histo->SetMaximum(0.3);
    histo->DrawCopy("same");

    leg->Draw();

    gSystem->mkdir(sOutputPath.Data(),1);
    can->SaveAs(Form("%s/%s.pdf",sOutputPath.Data(), sHistoName.Data()),"pdf");
}
