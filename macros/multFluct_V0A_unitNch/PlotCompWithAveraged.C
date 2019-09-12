// Macro for plotting v2 in double differential event selection (i.e. 0-10 V0A and tiny bins of Nch) with nCH integrated
// To studt mult. fluctuations of v2{4}

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>

void SetCustomPalette();
void PlotList(TList* list, TString sDraw = "", TLegend* leg = nullptr, std::vector<TString> labels = {});

void PlotCompWithAveraged()
{
    TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch40_80/";
    TString sSpecies = "Charged";
    TString sFileIn = Form("Averaged_%s.root",sSpecies.Data());
    TString sHistoName = Form("%s_hFlow4_harm2_gap-10",sSpecies.Data());
    Int_t mark = kFullCircle;
    TString sPlotSingleFrom = "bins5"; Int_t iNumCent = 8;


    std::vector<TString> binnings = {
        "bins20",
        "bins10",
        "bins5",
        "bins1"
    };

    std::vector<TString> labels = {};
    TList* list = new TList();

    // adding v2{4} integrated
    TFile* fileInt = TFile::Open(Form("%s/binsInt/Processed.root",sPath.Data()),"READ");
    if(!fileInt) {
        Utils::Error(Form("File '%s' not found!","binsInt/Processed.root"));
        delete list;
        return;
    }

    TH1* hist = (TH1*) fileInt->Get(Form("%s_cent0",sHistoName.Data()));
    if(!hist) {
        Utils::Error(Form("Histo '%s' not found!",sHistoName.Data()));
        fileInt->ls();
        delete list;
        return;
    }

    hist->SetLineWidth(2);
    hist->SetMarkerStyle(mark);
    list->Add(hist);
    labels.push_back("integ");

    // adding different binnins (averaged)
    Int_t iNumHistos = binnings.size();
    for(Int_t iH(0); iH < iNumHistos; ++iH) {
        TString sFile = Form("%s/%s/%s", sPath.Data(), binnings.at(iH).Data(), sFileIn.Data());
        TFile* file = TFile::Open(sFile.Data(),"READ");
        if(!file) {
            Utils::Error(Form("File '%s' not found!",sFile.Data()));
            delete list;
            return;
        }

        TH1* hist = (TH1*) file->Get(sHistoName.Data());
        if(!hist) {
            Utils::Error(Form("Histo '%s' not found!",sHistoName.Data()));
            file->ls();
            delete list;
            return;
        }

        hist->SetLineWidth(2);
        hist->SetMarkerStyle(mark);
        list->Add(hist);
        labels.push_back(binnings.at(iH));
    }

    TCanvas* can = new TCanvas();
    TLegend* leg = Utils::MakeLegend(Utils::kLegBotLeft);

    can->cd();
    TH1* frame = (TH1*) gPad->DrawFrame(0.0,-0.5,4.0,0.5);
    frame->SetTitle(sHistoName.Data());

    // additional single bins from unit
    if(sPlotSingleFrom.Length()) {
        TList* listUnit = new TList();

        TFile* fileUnit = TFile::Open(Form("%s/%s/Processed.root",sPath.Data(),sPlotSingleFrom.Data()),"READ");
        if(!fileUnit) {
            Utils::Error(Form("File '%s' not found!","bins1/Processed.root"));
            delete listUnit;
            return;
        }

        for (Int_t i = 0; i < iNumCent; i++) {
            TH1D* hist = (TH1D*) fileUnit->Get(Form("%s_cent%d",sHistoName.Data(),i));
            if(!hist) {
                Utils::Error(Form("Histo unit '%s' index %d not found!",sHistoName.Data(),i));
                fileUnit->ls();
                delete listUnit;
                return;
            }


            listUnit->Add(hist);
        }

        TLegend* legUnit = new TLegend(0.9,0.03,1.0,0.98);
        legUnit->SetBorderSize(0);

        can->cd();
        PlotList(listUnit,"hist l", legUnit);
        legUnit->Draw();
    }

    PlotList(list,"e1",leg,labels);
    leg->SetHeader(sPath.Data());
    leg->Draw();

    can->SaveAs(Form("%s/CompAvgs_%s.pdf",sPath.Data(), sHistoName.Data()), "pdf");

    return;
}

void PlotList(TList* list, TString sDraw, TLegend* leg, std::vector<TString> labels)
{
    if(!list) { Utils::Error("Input list not found!","PlotList"); return; }

    SetCustomPalette();

    Int_t iNum = list->GetEntries();

    Int_t nnCol = gStyle->GetNumberOfColors();
    // printf("%d\n",nnCol);

    for (Int_t i = 0; i < iNum; i++) {
      Int_t idx = i * Float_t(nnCol-1) / (iNum);
      // Int_t idx = i * Float_t(nnCol-1) / (iNum-1);
      Color_t col = gStyle->GetColorPalette(idx);

      TH1* hist = (TH1*) list->At(i);
      if(!hist) { Utils::Error(Form("Histo at index %d not found!",i),"PlotList"); return; }

      hist->SetLineColor(col);
      hist->SetMarkerColor(col);
      hist->Draw(Form("same %s",sDraw.Data()));
      if(leg) {
          TString sLabel = "";
          labels.size() ? sLabel = labels.at(i) : sLabel = Form("%d",i);
          leg->AddEntry(hist,sLabel.Data(),"l");
      }

    }

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
    // //case 57:
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
