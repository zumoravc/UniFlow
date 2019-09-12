#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>

void SetCustomPalette();


void PlotTinyBins()
{
    TString sFileIn = "./../../../results/pPb/unitMultBin/V0A0_10_Nch00_140/bins1/Processed.root"; Int_t iCent = 140;
    // TString sFileIn = "./V0A0_10_Nch00_140/bins5/Processed.root"; Int_t iCent = 28;
    // TString sFileIn = "./V0A0_10_Nch00_140/bins10/Processed.root"; Int_t iCent = 14;
    TString sHistoName = "Charged_hFlow4_harm2_gap-10";

    TFile* file = TFile::Open(sFileIn.Data(),"READ");

    SetCustomPalette();

    TCanvas* can = new TCanvas();
    TLegend* leg = Utils::MakeLegend(Utils::kLegTopLeft);

    Int_t nnCol = gStyle->GetNumberOfColors();
    // printf("%d\n",nnCol);

    for (Int_t i = 0; i < iCent; i++) {
      Int_t idx = i * Float_t(nnCol-1) / (iCent-1);
      Color_t col = gStyle->GetColorPalette(idx);

      TH1D* hist = (TH1D*) file->Get(Form("%s_cent%d",sHistoName.Data(),i));
      if(!hist) return;

      hist->SetLineColor(col);
      hist->SetMarkerColor(col);

      can->cd();
      hist->SetStats(0);
      hist->SetMinimum(-1.0);
      hist->SetMaximum(1.0);
      hist->Draw("same hist l");
      leg->AddEntry(hist,Form("%d",i),"l");
    }


    can->cd();
    leg->Draw();


    return;
}

void SetCustomPalette()
{
  // Setting custom (ala ROOT6) color palette
  // See https://root.cern.ch/doc/master/TColor_8cxx_source.html#l02400
  // for definition of color setting (array bellow) for ROOT 6 defined palettes


    Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};

    // // // Rain Bow
    // Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
    // Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
    // Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};

    // Visible Spectrum
    Double_t red[9]   = { 18./255.,  72./255.,   5./255.,  23./255.,  29./255., 201./255., 200./255., 98./255., 29./255.};
    Double_t green[9] = {  0./255.,   0./255.,  43./255., 167./255., 211./255., 117./255.,   0./255.,  0./255.,  0./255.};
    Double_t blue[9]  = { 51./255., 203./255., 177./255.,  26./255.,  10./255.,   9./255.,   8./255.,  3./255.,  0./255.};


    // Bird
    //case 57:
    // Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
    // Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
    // Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

    // Blue Green Yellow
    // //case 71:
    // Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
    // Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
    // Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};

    // Viridis
    // case 112:
    // Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
    // Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
    // Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};

    // Cividis
    // case 113:
    // Double_t red[9]   = {  0./255.,   5./255.,  65./255.,  97./255., 124./255., 156./255., 189./255., 224./255., 255./255.};
    // Double_t green[9] = { 32./255.,  54./255.,  77./255., 100./255., 123./255., 148./255., 175./255., 203./255., 234./255.};
    // Double_t blue[9]  = { 77./255., 110./255., 107./255., 111./255., 120./255., 119./255., 111./255.,  94./255.,  70./255.};


  Int_t pal = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
  const Int_t nCol = 255;
  Int_t colors[nCol];
  for (int i=0; i<nCol; i++) colors[i] = pal+i;

  gStyle->SetPalette(nCol,colors);
}
