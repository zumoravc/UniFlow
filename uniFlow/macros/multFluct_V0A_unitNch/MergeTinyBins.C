// Macro for averaging results from tiny bins (~unit Nch) into integrated one

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>
#include <vector>

void SetCustomPalette();

void GetWeightedMean(std::vector<Double_t> cont, std::vector<Double_t> err, Double_t& dAverage, Double_t& dError);
TH1* AverageHistos(TList* list);

void MergeTinyBins()
{
    TString sFileIn = "Processed.root";
    TString sSpecies = "Proton";
    TString sFileOut = Form("Averaged_%s.root",sSpecies.Data());

    std::vector<TString> histos = {
        "pCor2_harm2_gap-10",
        "hCum2_harm2_gap-10",
        "hFlow2_harm2_gap-10",
        "pCor2_harm2_gap08",
        "hCum2_harm2_gap08",
        "hFlow2_harm2_gap08",
        "pCor4_harm2_gap-10",
        "hCum4_harm2_gap-10",
        "hFlow4_harm2_gap-10"
    };

    TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch00_140/bins1/"; Int_t iCent = 140;
    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch00_140/bins5/"; Int_t iCent = 28;
    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch00_140/bins10/"; Int_t iCent = 14;
    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch00_140/bins20/"; Int_t iCent = 7;

    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch40_80/bins1/"; Int_t iCent = 40;
    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch40_80/bins5/"; Int_t iCent = 8;
    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch40_80/bins10/"; Int_t iCent = 4;
    // TString sPath = "./../../results/pPb/unitMultBin/V0A0_10_Nch40_80/bins20/"; Int_t iCent = 2;

    TFile* file = TFile::Open(Form("%s/%s",sPath.Data(),sFileIn.Data()),"READ");
    TFile* fileOut = TFile::Open(Form("%s/%s",sPath.Data(),sFileOut.Data()),"RECREATE");

    Int_t iNumHistos = histos.size();

    for(Int_t iH(0); iH < iNumHistos; ++iH) {
        TString sHistoName = Form("%s_%s",sSpecies.Data(), histos.at(iH).Data());

        TList* list = new TList();
        list->SetOwner(kTRUE);

        for (Int_t i = 0; i < iCent; i++) {
            TH1D* hist = (TH1D*) file->Get(Form("%s_cent%d",sHistoName.Data(),i));
            if(!hist) {
                Utils::Error(Form("Histo %s (cent %d) not found!",sHistoName.Data(), i));
                file->ls();
                delete list;
                return;
            }

            list->Add(hist);
        }

        TH1D* average = (TH1D*) AverageHistos(list);
        average->SetName(sHistoName.Data());
        average->Draw();

        fileOut->cd();
        average->Write(sHistoName.Data());

        delete list;
    }

    file->Close();
    fileOut->Close();

    return;
}

TH1* AverageHistos(TList* list)
{
    TH1* hist = (TH1*) list->At(0)->Clone("temp");
    if(!hist) { Utils::Error("0-th histo not found!","AverageHistos"); return nullptr; }

    Int_t iNumBins = hist->GetNbinsX();
    Int_t iNumHistos = list->GetEntries();

    for(Int_t iB(1); iB < iNumBins+1; ++iB) {
        std::vector<Double_t> val = {};
        std::vector<Double_t> err = {};

        for(Int_t iH(0); iH < iNumHistos; ++iH) {
            TH1* h = (TH1*) list->At(iH);

            Double_t dCont = h->GetBinContent(iB);
            Double_t dErr = h->GetBinError(iB);

            if(dErr > 100) { continue; }

            val.push_back(dCont);
            err.push_back(dErr);
        }

        Double_t dVal = 0.0;
        Double_t dErr = 0.0;

        GetWeightedMean(val,err,dVal,dErr);
        hist->SetBinContent(iB,dVal);
        hist->SetBinError(iB,dErr);
    }

    return hist;
}

void GetWeightedMean(std::vector<Double_t> cont, std::vector<Double_t> err, Double_t& dAverage, Double_t& dError)
{
    Int_t iNumValues = cont.size();

    Double_t dSum = 0.0;
    Double_t dSumErr = 0.0;

    for(Int_t i(0); i < iNumValues; ++i) {
        Double_t dCon = cont.at(i);
        Double_t dErr = err.at(i);

        if(!(dErr > 0.0)) { continue; }

        // Double_t dW = TMath::Power(dErr,-2);
        Double_t dW = 1.0;

        dSum += dCon * dW;
        dSumErr += dW;
    }

    dAverage = 0.0;
    dError = 0.0;

    if(!(dSumErr > 0.0)) return;

    dAverage = dSum / dSumErr;
    dError = 1.0 / TMath::Sqrt(dSumErr);
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
