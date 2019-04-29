// Macro for plotting two panel figure with inv. mass distribution & corr-mass figure
// Author: Vojtech Pacik (2019)

void PlotFlowMass() {

    TString sPathIn = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/non_linear_modes/output/K0s/fits.root";
    TString sListName = "K0s_<<3>>(4,-2,-2)_2sub(0)_cent2_pt5;1";
    TString sHistMass = "histMass";
    TString sHistCorr = "histCorr";

    TFile* fileIn = TFile::Open(sPathIn.Data(),"READ");
    if(!fileIn) { printf("FileIn not open!\n"); return; }

    TList* listIn = (TList*) fileIn->Get(sListName.Data());
    if(!listIn) { printf("List '%s' not found!\n",sListName.Data()); fileIn->ls(); return; }


    TH1D* hMass = (TH1D*) listIn->FindObject(sHistMass.Data());
    if(!hMass) { printf("hMass '%s' not found!\n",sHistMass.Data()); listIn->ls(); return; }

    TH1D* hCorr = (TH1D*) listIn->FindObject(sHistCorr.Data());
    if(!hCorr) { printf("hCorr '%s' not found!\n",sHistCorr.Data()); listIn->ls(); return; }

    // temp
    TH1D* h1 = hMass;
    // TH1D* h2 = hMass;
    TH1D* h3 = hCorr;


    // h1->SetMaximum(0);
        // h1->SetMaximum(350000);
        // Define the ratio plot
        // TH1F *h3 = (TH1F*)h1->Clone("h3");
        h3->SetLineColor(kBlack);
        // h3->SetMinimum(0.8);  // Define Y ..
        // h3->SetMaximum(1.35); // .. range
        // h3->Sumw2();
        // h3->SetStats(0);      // No statistics on lower plot
        h3->SetMarkerStyle(21);

        // h1 settings
        h1->SetLineColor(kBlue+1);
        h1->SetLineWidth(2);



    Double_t dXmin = h1->GetXaxis()->GetXmin();
    Double_t dXmax = h1->GetXaxis()->GetXmax();

    Double_t dYmin = 0.0;
    Double_t dYmax = 300000;

    printf("X: min : %f | max %f\n",dXmin,dXmax);
    printf("Y: min : %f | max %f\n",dYmin,dYmax);


    Double_t dXmin2 = dXmin;
    Double_t dXmax2 = dXmax;
    Double_t dYmin2 = h3->GetMinimum();
    Double_t dYmax2 = h3->GetMaximum();


    // Define the Canvas
    TCanvas* c = new TCanvas("c", "canvas", 600, 600);

    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
    pad1->SetBottomMargin(0.01); // Upper and lower plot are joined
    pad1->SetRightMargin(0.1);
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1

    // c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.4);
    pad2->SetTopMargin(0.01);
    // pad2->SetLeftMargin(0.01);
    pad2->SetRightMargin(0.1);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();

    pad1->cd();
    TH1* frame1 = (TH1*) gPad->DrawFrame(dXmin, dYmin, dXmax, dYmax);
    h1->Draw("same ep");

    pad2->cd();
    TH1* frame2 = (TH1*) gPad->DrawFrame(dXmin2, dYmin2, dXmax2, dYmax2);
    h3->Draw("same ep");

    // TGaxis *axis = new TGaxis( dXmin, dYmin, dXmin, dYmax, dYmin,dYmax,510,"");
    // axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    // axis->SetLabelSize(15);
    // axis->Draw();


    // Main frame (frame1) setting
    frame1->GetYaxis()->SetTitle("Counts");
    // hiding X axis label (due to margin)
    frame1->GetXaxis()->SetLabelSize(0.0);
    // Y axis h1 plot settings
    frame1->GetYaxis()->SetTitleSize(20);
    frame1->GetYaxis()->SetTitleFont(43);
    frame1->GetYaxis()->SetTitleOffset(1.55);

    // Ratio plot (frame2) settings
    frame2->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    frame2->GetYaxis()->SetTitle("dn");
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetYaxis()->SetTitleSize(20);
    frame2->GetYaxis()->SetTitleFont(43);
    frame2->GetYaxis()->SetTitleOffset(1.55);
    frame2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame2->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    frame2->GetXaxis()->SetTitle("M_{inv} (GeV/#it{c}^{2})");
    frame2->GetXaxis()->SetTitleSize(20);
    frame2->GetXaxis()->SetTitleFont(43);
    frame2->GetXaxis()->SetTitleOffset(4.);
    frame2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame2->GetXaxis()->SetLabelSize(15);

    c->SaveAs("can.pdf","pdf");

}
}
