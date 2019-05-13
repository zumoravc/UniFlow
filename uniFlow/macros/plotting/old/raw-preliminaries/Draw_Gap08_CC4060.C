void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};


void Draw_Gap08_CC4060() {
   
    
    const Double_t sys_Charged = 0.02;
    const Double_t sys_K0 = 0.11;
    const Double_t sys_Lambda = 0.13;
    const Double_t sys_Phi = 0.07;
    
    // Load necessary libraries
    LoadLibs();
    // Set the default style
    SetStyle();
    

    TFile* Root1 = TFile::Open("/Users/youzhou/Dropbox/Research/ALICE/Flow_ch/Code_502TeVpPb/VP/Plots/UniFlow.root");
    TH1D* fv2_Charged_gap08_cent2 = (TH1D*)Root1->Get("hFlow2_Charged_harm2_gap08_cent2");
    
    TH1D* fv2_K0_gap08_cent2 = (TH1D*)Root1->Get("hFlow2_K0s_harm2_gap08_mult2");
    
    
    TH1D* fv2_Lambda_gap08_cent2 = (TH1D*)Root1->Get("hFlow2_Lambda_harm2_gap08_mult2");
    
    TH1D* fv2_Phi_gap08_cent2 = (TH1D*)Root1->Get("hFlow2_Phi_harm2_gap08_mult2");
    
    Double_t pt_Charged[]={0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.25, 3.75, 4.5, 5.5};
    Double_t pt_K0[]={0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.25, 3.75, 4.5, 5.5};
    Double_t pt_Lambda[]={0.95, 1.25, 1.5, 1.7, 1.9, 2.1, 2.3, 2.55, 2.95, 3.4, 3.95, 4.7975, 5.48229, 5.85687, 6.23146};
    Double_t pt_Phi[] = {0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 5};
    Double_t ptErr[]={0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.02};
    Double_t ptBin[20]={0.};
  
    Double_t v22Gap08_Charged_cc4060[20]={0.};
    Double_t v22Gap08Err_Charged_cc4060[20]={0.};
    Double_t v22Gap08Sys_Charged_cc4060[20]={0.};
    
    Double_t v22Gap08_K0_cc4060[20]={0.};
    Double_t v22Gap08Err_K0_cc4060[20]={0.};
    Double_t v22Gap08Sys_K0_cc4060[20]={0.};

    Double_t v22Gap08_Lambda_cc4060[20]={0.};
    Double_t v22Gap08Err_Lambda_cc4060[20]={0.};
    Double_t v22Gap08Sys_Lambda_cc4060[20]={0.};
    
    Double_t v22Gap08_Phi_cc4060[20]={0.};
    Double_t v22Gap08Err_Phi_cc4060[20]={0.};
    Double_t v22Gap08Sys_Phi_cc4060[20]={0.};
    
    for(Int_t k=0;k<20;k++)
    {
        //cout << fv2_Phi_gap08_cent2->GetBinCenter(k+1)<<", ";
        v22Gap08_Charged_cc4060[k] = fv2_Charged_gap08_cent2->GetBinContent(k+1);
        v22Gap08Err_Charged_cc4060[k] = fv2_Charged_gap08_cent2->GetBinError(k+1);
        v22Gap08Sys_Charged_cc4060[k] = sys_Charged * fv2_Charged_gap08_cent2->GetBinContent(k+1);
        
    }
    
    for(Int_t k=0;k<17;k++)
    {
        //cout << fv2_Phi_gap08_cent2->GetBinCenter(k+1)<<", ";
        v22Gap08_K0_cc4060[k] = fv2_K0_gap08_cent2->GetBinContent(k+1);
        v22Gap08Err_K0_cc4060[k] = fv2_K0_gap08_cent2->GetBinError(k+1);
        v22Gap08Sys_K0_cc4060[k] = sys_K0 * fv2_K0_gap08_cent2->GetBinContent(k+1);
    }
    
    for(Int_t k=0;k<15;k++)
    {
        //cout << fv2_Phi_gap08_cent2->GetBinCenter(k+1)<<", ";
        v22Gap08_Lambda_cc4060[k] = fv2_Lambda_gap08_cent2->GetBinContent(k+1);
        v22Gap08Err_Lambda_cc4060[k] = fv2_Lambda_gap08_cent2->GetBinError(k+1);
        v22Gap08Sys_Lambda_cc4060[k] = sys_Lambda * fv2_Lambda_gap08_cent2->GetBinContent(k+1);
        
    }
    
    for(Int_t k=0;k<15;k++)
    {
        //cout << fv2_Phi_gap08_cent2->GetBinCenter(k+1)<<", ";
        v22Gap08_Phi_cc4060[k] = fv2_Phi_gap08_cent2->GetBinContent(k+1);
        v22Gap08Err_Phi_cc4060[k] = fv2_Phi_gap08_cent2->GetBinError(k+1);
        v22Gap08Sys_Phi_cc4060[k] = sys_Phi * fv2_Phi_gap08_cent2->GetBinContent(k+1);
        
    }
    
    
    TGraphErrors *figv2_Charged_gap08_cent2 = new TGraphErrors(20,pt_Charged,v22Gap08_Charged_cc4060,ptBin, v22Gap08Err_Charged_cc4060);
    figv2_Charged_gap08_cent2->SetMarkerStyle(20);
    figv2_Charged_gap08_cent2->SetMarkerColor(1);
    figv2_Charged_gap08_cent2->SetMarkerSize(1.);
    figv2_Charged_gap08_cent2->SetLineColor(1);
    figv2_Charged_gap08_cent2->SetLineWidth(2);
    figv2_Charged_gap08_cent2->SetLineStyle(1);
    
    TGraphErrors *gv2_Charged_gap08_cent2 = new TGraphErrors(20,pt_Charged,v22Gap08_Charged_cc4060,ptErr, v22Gap08Sys_Charged_cc4060);
    gv2_Charged_gap08_cent2->SetMarkerStyle(20);
    gv2_Charged_gap08_cent2->SetMarkerColor(1);
    gv2_Charged_gap08_cent2->SetMarkerSize(1.);
    gv2_Charged_gap08_cent2->SetLineColor(kGray);
    gv2_Charged_gap08_cent2->SetLineWidth(2);
    gv2_Charged_gap08_cent2->SetLineStyle(1);
    gv2_Charged_gap08_cent2->SetFillStyle(1001);
    gv2_Charged_gap08_cent2->SetFillColor(kGray);

    
    
    TGraphErrors *figv2_K0_gap08_cent2 = new TGraphErrors(15,pt_K0,v22Gap08_K0_cc4060,ptBin, v22Gap08Err_K0_cc4060);
    figv2_K0_gap08_cent2->SetMarkerStyle(21);
    figv2_K0_gap08_cent2->SetMarkerColor(kCyan+2);
    figv2_K0_gap08_cent2->SetMarkerSize(1.);
    figv2_K0_gap08_cent2->SetLineColor(kCyan+2);
    figv2_K0_gap08_cent2->SetLineWidth(2);
    figv2_K0_gap08_cent2->SetLineStyle(1);
    
    TGraphErrors *gv2_K0_gap08_cent2 = new TGraphErrors(15,pt_K0,v22Gap08_K0_cc4060,ptErr, v22Gap08Sys_K0_cc4060);
    gv2_K0_gap08_cent2->SetMarkerStyle(21);
    gv2_K0_gap08_cent2->SetMarkerColor(kCyan+2);
    gv2_K0_gap08_cent2->SetMarkerSize(1.);
    gv2_K0_gap08_cent2->SetLineColor(kCyan-8);
    gv2_K0_gap08_cent2->SetLineWidth(2);
    gv2_K0_gap08_cent2->SetLineStyle(1);
    gv2_K0_gap08_cent2->SetFillStyle(1001);
    gv2_K0_gap08_cent2->SetFillColor(kCyan-8);

    
    TGraphErrors *figv2_Lambda_gap08_cent2 = new TGraphErrors(18,pt_Lambda,v22Gap08_Lambda_cc4060,ptBin, v22Gap08Err_Lambda_cc4060);
    figv2_Lambda_gap08_cent2->SetMarkerStyle(33);
    figv2_Lambda_gap08_cent2->SetMarkerColor(kOrange-1);
    figv2_Lambda_gap08_cent2->SetMarkerSize(1.5);
    figv2_Lambda_gap08_cent2->SetLineColor(kOrange-1);
    figv2_Lambda_gap08_cent2->SetLineWidth(2);
    figv2_Lambda_gap08_cent2->SetLineStyle(1);
    
    TGraphErrors *gv2_Lambda_gap08_cent2 = new TGraphErrors(18,pt_Lambda,v22Gap08_Lambda_cc4060,ptErr, v22Gap08Sys_Lambda_cc4060);
    gv2_Lambda_gap08_cent2->SetMarkerStyle(33);
    gv2_Lambda_gap08_cent2->SetMarkerColor(kOrange-1);
    gv2_Lambda_gap08_cent2->SetMarkerSize(1.5);
    gv2_Lambda_gap08_cent2->SetLineColor(kOrange-9);
    gv2_Lambda_gap08_cent2->SetLineWidth(2);
    gv2_Lambda_gap08_cent2->SetLineStyle(1);
    gv2_Lambda_gap08_cent2->SetFillStyle(1001);
    gv2_Lambda_gap08_cent2->SetFillColor(kOrange-9);

   
    TGraphErrors *figv2_Phi_gap08_cent2 = new TGraphErrors(18,pt_Phi,v22Gap08_Phi_cc4060,ptBin, v22Gap08Err_Phi_cc4060);
    figv2_Phi_gap08_cent2->SetMarkerStyle(34);
    figv2_Phi_gap08_cent2->SetMarkerColor(kMagenta+1);
    figv2_Phi_gap08_cent2->SetMarkerSize(1.3);
    figv2_Phi_gap08_cent2->SetLineColor(kMagenta+1);
    figv2_Phi_gap08_cent2->SetLineWidth(2);
    figv2_Phi_gap08_cent2->SetLineStyle(1);
    
    TGraphErrors *gv2_Phi_gap08_cent2 = new TGraphErrors(18,pt_Phi,v22Gap08_Phi_cc4060,ptErr, v22Gap08Sys_Phi_cc4060);
    gv2_Phi_gap08_cent2->SetMarkerStyle(34);
    gv2_Phi_gap08_cent2->SetMarkerColor(kMagenta+1);
    gv2_Phi_gap08_cent2->SetMarkerSize(1.3);
    gv2_Phi_gap08_cent2->SetLineColor(kMagenta-9);
    gv2_Phi_gap08_cent2->SetLineWidth(2);
    gv2_Phi_gap08_cent2->SetLineStyle(1);
    gv2_Phi_gap08_cent2->SetFillStyle(1001);
    gv2_Phi_gap08_cent2->SetFillColor(kMagenta-9);
    
    
   /*
    fv2_K0s_gap08_cent2->SetMarkerStyle(kFullDiamond);
    fv2_K0s_gap08_cent2->SetMarkerStyle(kFullDiamond);
    fv2_K0s_gap08_cent2->SetMarkerStyle(kFullDiamond);
    fv2_K0s_gap08_cent3->SetMarkerStyle(kFullDiamond);
    fv2_K0s_gap08_cent4->SetMarkerStyle(kFullDiamond);
    fv2_K0s_gap08_cent2->SetMarkerColor(kMagenta+1);
    fv2_K0s_gap08_cent2->SetMarkerColor(kMagenta+1);
    fv2_K0s_gap08_cent2->SetMarkerColor(kMagenta+1);
    fv2_K0s_gap08_cent3->SetMarkerColor(kMagenta+1);
    fv2_K0s_gap08_cent4->SetMarkerColor(kMagenta+1);
    fv2_K0s_gap08_cent2->SetLineColor(kMagenta+1);
    fv2_K0s_gap08_cent2->SetLineColor(kMagenta+1);
    fv2_K0s_gap08_cent2->SetLineColor(kMagenta+1);
    fv2_K0s_gap08_cent3->SetLineColor(kMagenta+1);
    fv2_K0s_gap08_cent4->SetLineColor(kMagenta+1);


    fv2_Lambda_gap08_cent2->SetMarkerStyle(kFullStar);
    fv2_Lambda_gap08_cent2->SetMarkerStyle(kFullStar);
    fv2_Lambda_gap08_cent2->SetMarkerStyle(kFullStar);
    fv2_Lambda_gap08_cent3->SetMarkerStyle(kFullStar);
    fv2_Lambda_gap08_cent4->SetMarkerStyle(kFullStar);
    fv2_Lambda_gap08_cent2->SetMarkerColor(kOrange-1);
    fv2_Lambda_gap08_cent2->SetMarkerColor(kOrange-1);
    fv2_Lambda_gap08_cent2->SetMarkerColor(kOrange-1);
    fv2_Lambda_gap08_cent3->SetMarkerColor(kOrange-1);
    fv2_Lambda_gap08_cent4->SetMarkerColor(kOrange-1);
    fv2_Lambda_gap08_cent2->SetLineColor(kOrange-1);
    fv2_Lambda_gap08_cent2->SetLineColor(kOrange-1);
    fv2_Lambda_gap08_cent2->SetLineColor(kOrange-1);
    fv2_Lambda_gap08_cent3->SetLineColor(kOrange-1);
    fv2_Lambda_gap08_cent4->SetLineColor(kOrange-1);
    */
    
    
   
    TLegend *lALICE = new TLegend(0.2,0.9,0.4,0.96);
    lALICE->SetFillColor(0);
    lALICE->SetBorderSize(1);
    lALICE->SetTextSize(0.04);
    lALICE->SetTextFont(42);
    lALICE->SetLineColor(0);
    lALICE->SetHeader("ALICE Preliminary");
   
    
    TLatex *TEXT1 = new TLatex(0.35, 0.275, "ALICE Preliminary");
    TEXT1->SetTextFont(42);
    TEXT1->SetTextSize(0.04);
    
    TLatex *TEXT2 = new TLatex(4., 0.275, "40-60% (V0-A)");
    TEXT2->SetTextFont(42);
    TEXT2->SetTextSize(0.04);
    
    // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
    // Rectangular
    TCanvas *cfig_cent2 = new TCanvas("cfig_cent2", "Alice Figure Template", 800, 600);
    gStyle->SetOptStat(0);
    TH1 * h_cent2 = cfig_cent2->DrawFrame(0,0.01,6,0.35);
    const char *  texX="#it{p}_{T} (GeV/#it{c})";
    const char *  texY="#it{v}_{2}{2,|#Delta#eta|>0.8}";
    h_cent2->SetXTitle(texX);
    h_cent2->SetYTitle(texY);
    
    // Legend for ALICE:
    TLegend *legend_pidv2_cent2 = new TLegend(0.19,0.56,0.35,0.82);
    legend_pidv2_cent2->SetFillStyle(0); // white legend background
    legend_pidv2_cent2->SetTextSize(0.04);
    legend_pidv2_cent2->SetTextFont(42);
    legend_pidv2_cent2->SetLineColor(0);
    legend_pidv2_cent2->SetMargin(0.3);
    legend_pidv2_cent2->SetHeader("p-Pb #sqrt{s_{NN}} = 5.02 TeV");
    legend_pidv2_cent2->AddEntry(gv2_Charged_gap08_cent2,"h^{ #pm} ","pf");
    //legend_pidv2_cent2->AddEntry(fv2_Kaon_gap08_cent2,"#it{K}^{ #pm}, #it{v}_{2}{2,|#Delta#eta|>0.8}","p");
    legend_pidv2_cent2->AddEntry(gv2_K0_gap08_cent2,"#it{K}^{ 0}_{S} ","pf");
    //legend_pidv2_cent2->AddEntry(fv2_Proton_gap08_cent2,"#it{p} (#it{#bar{p}}), #it{v}_{2}{2,|#Delta#eta|>0.8}","p");
    //legend_pidv2_cent2->AddEntry(gv2_Phi_gap08_cent2,"#it{#phi} ","pf");
    legend_pidv2_cent2->AddEntry(gv2_Lambda_gap08_cent2,"#it{#Lambda} (#it{#bar{#Lambda}}) ","pf");
    
    legend_pidv2_cent2->Draw();
    
    TEXT1->Draw();
     TEXT2->Draw();
    
   // gv2_Phi_gap08_cent2->Draw("2psame");
    gv2_Lambda_gap08_cent2->Draw("2psame");
    gv2_K0_gap08_cent2->Draw("2psame");
    gv2_Charged_gap08_cent2->Draw("2psame");
    

   // figv2_Phi_gap08_cent2->Draw("pZsame");
    figv2_Lambda_gap08_cent2->Draw("pZsame");
    figv2_K0_gap08_cent2->Draw("pZsame");
    figv2_Charged_gap08_cent2->Draw("pzsame");
    
    
    
    //fv2_Kaon_gap08_cent2->Draw("psame");
    //fv2_Proton_gap08_cent2->Draw("psame");
    //fv2_K0s_gap08_cent2->Draw("psame");
    //fv2_Lambda_gap08_cent2->Draw("psame");
    

  
    
    
    
}

//________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}


void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);


}

void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  

}
