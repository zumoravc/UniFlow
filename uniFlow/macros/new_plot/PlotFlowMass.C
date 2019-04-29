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

    // Define two gaussian histograms. Note the X and Y title are defined
    // at booking time using the convention "Hist_title ; X_title ; Y_title"
    // Define two gaussian histograms. Note the X and Y title are defined
    // at booking time using the convention "Hist_title ; X_title ; Y_title"

    // TH1F *h1 = new TH1F("h1", "Two gaussian plots and their ratio;x title; h1 and h2 gaussian histograms", 100, -5, 5);
    // h1->FillRandom("gaus");

    // Define the Canvas
    TCanvas *c = new TCanvas("c", "canvas", 800, 800);

    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.01); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    h1->SetStats(0);          // No statistics on upper plot
    h1->Draw();               // Draw h1

    // Do not draw the Y axis label on the upper plot and redraw a small
    // axis instead, in order to avoid the first label (0) to be clipped.
    // h1->GetYaxis()->SetLabelSize(0.0);
    // removing original ticks
    // h1->GetYaxis()->SetTickLength(0.0);

    h1->SetMaximum(0);
    h1->SetMaximum(350000);

    // hiding X axis label (due to margin)
    h1->GetXaxis()->SetLabelSize(0.0);


    Double_t dXmin = h1->GetXaxis()->GetXmin();
    Double_t dXmax = h1->GetXaxis()->GetXmax();
    // Double_t dYmin = h1->GetMinimum();
    // Double_t dYmax = h1->GetMaximum();
    Double_t dYmin = 0.0;
    Double_t dYmax = 300000;

    printf("X: min : %f | max %f\n",dXmin,dXmax);
    printf("Y: min : %f | max %f\n",dYmin,dYmax);

    // TGaxis *axis = new TGaxis( dXmin, dYmin, dXmin, dYmax, dYmin,dYmax,510,"");
    // axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    // axis->SetLabelSize(15);
    // axis->Draw();

    // lower plot will be in pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    // Define the ratio plot
    // TH1F *h3 = (TH1F*)h1->Clone("h3");
    h3->SetLineColor(kBlack);
    // h3->SetMinimum(0.8);  // Define Y ..
    // h3->SetMaximum(1.35); // .. range
    h3->Sumw2();
    h3->SetStats(0);      // No statistics on lower plot
    h3->SetMarkerStyle(21);
    h3->Draw("ep");       // Draw the ratio plot

    // h1 settings
    h1->SetLineColor(kBlue+1);
    h1->SetLineWidth(2);

    // Y axis h1 plot settings
    h1->GetYaxis()->SetTitleSize(20);
    h1->GetYaxis()->SetTitleFont(43);
    h1->GetYaxis()->SetTitleOffset(1.55);

    // Ratio plot (h3) settings
    h3->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    h3->GetYaxis()->SetTitle("ratio h1/h2 ");
    h3->GetYaxis()->SetNdivisions(505);
    h3->GetYaxis()->SetTitleSize(20);
    h3->GetYaxis()->SetTitleFont(43);
    h3->GetYaxis()->SetTitleOffset(1.55);
    h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h3->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    h3->GetXaxis()->SetTitleSize(20);
    h3->GetXaxis()->SetTitleFont(43);
    h3->GetXaxis()->SetTitleOffset(4.);
    h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h3->GetXaxis()->SetLabelSize(15);
}
