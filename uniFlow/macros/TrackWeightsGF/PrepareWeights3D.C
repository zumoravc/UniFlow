/* PrepareWeights
 *
 * Macro for extracting the (eta,phi,z) distribution and make them into 3D (eta,phi,z) and 2D (eta,phi) weights
 * for Generic Framework calculations and preparing source ROOT file used for running the analysis.
 *
 * This can be done on run-integrated results and by run-by-run basis (if bRunByRun = kTRUE)
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

 #include "TFile.h"
 #include "TObject.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "TH3.h"
 #include "TCanvas.h"
 #include "TSystem.h"
 #include "TPad.h"

TString sTag = "_bayes90";
TString sTaskTag = "UniFlow" + sTag;
TString sPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/pid/pPb-16q-FAST";
TString sOutputPath = sPath+"/weights"+sTag+"/";
TString sOutFileName = "weights"+sTag+"_16q_FAST.root";

const Short_t iNumPart = 8;
const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
// const Short_t iNumPart = 1;
// const TString species[iNumPart] = {"Refs"};

const Bool_t bRunByRun = kTRUE;
// Run Lists

// pPb 5.02 TeV
// RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20171129_v2.txt [31 runs]
Int_t iNumRuns = 31; Int_t iRunList[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
// RunList_LHC16t_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [4 runs]
// Int_t iNumRuns = 4; Int_t iRunList[] = {267166, 267165, 267164, 267163};

// pp 13 TeV
// RunList_LHC16l_pass1_CentralBarrelTracking_hadronPID_20170509_v2.txt [70]
// Int_t iNumRuns = 70; Int_t iRunList[] = {259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259713, 259711, 259705, 259704, 259703, 259700, 259697, 259668, 259650, 259649, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962, 258923, 258919};

// RunList_LHC16k_pass1_CentralBarrelTracking_hadronPID_20170516_v3.txt [194]
// Int_t iNumRuns = 194; Int_t iRunList[] = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};



TH2D* ProcessWeights2D(const TH3D* dist);
TH3D* ProcessWeights3D(const TH3D* dist);
TH1D* CalculateWeight(TH1D* proj);
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh);

void PrepareWeights3D()
{
  gSystem->mkdir(sOutputPath.Data(),kTRUE);
  TFile* fOutput = new TFile(Form("%s/%s",sPath.Data(),sOutFileName.Data()),"RECREATE");

  TFile* fInput = TFile::Open(Form("%s/AnalysisResults.root",sPath.Data()),"READ");
  if(!fInput) return;

  fInput->cd(sTaskTag.Data());
  TList* list =  (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));
  if(!list) return;

  TList* listRun = new TList();
  listRun->SetOwner(kTRUE);

  for(Short_t part = 0; part < iNumPart; part++)
  {
    printf(" -part %d (out of %d)\n",part+1,iNumPart);

    TH3D* h3Weights = (TH3D*) list->FindObject(Form("fh3Weights%s",species[part].Data()));
    if(!h3Weights) { printf("Hist 'fh3Weights%s' not found\n",species[part].Data()); return; }

    // Plotting overview plots with dist slices in z
    // TCanvas* canDist = new TCanvas("canDist","canDist",1600,2000);
    // canDist->Divide(5,4);
    // for(Int_t iBinZ(1); iBinZ < h3Weights->GetNbinsZ()+1; ++iBinZ)
    // {
    //   TH2D* slice = SliceWeights(h3Weights,iBinZ,iBinZ);
    //   slice->SetStats(0);
    //   // if(!slice) continue;
    //
    //   canDist->cd(iBinZ);
    //   // gPad->SetLogz();
    //   slice->DrawCopy("colz");
    // }
    // canDist->SaveAs(Form("%s/dist_%s_all.pdf",sOutputPath.Data(), species[part].Data()), "pdf");

    TH2D* weights = ProcessWeights2D(h3Weights);
    if(!weights) { printf("Hist with 2D weights not processed!\n"); return; }

    // CheckWeights(h3Weights, iRunList[iRun]);
    weights->SetName(Form("%s",species[part].Data()));
    weights->SetTitle(Form("%s",species[part].Data()));
    listRun->Add(weights);

    TH3D* weights3D = ProcessWeights3D(h3Weights);
    if(!weights3D) { printf("Hist with 3D weights not processed!\n"); return; }

    // CheckWeights(h3Weights, iRunList[iRun]);
    weights3D->SetName(Form("%s3D",species[part].Data()));
    weights3D->SetTitle(Form("%s3D",species[part].Data()));
    listRun->Add(weights3D);

    // Plotting overview plots with dist slices in z
    // TCanvas* canSlices = new TCanvas("canSlices","canSlices",1600,2000);
    // canSlices->Divide(5,4);
    // for(Int_t iBinZ(1); iBinZ < weights->GetNbinsZ()+1; ++iBinZ)
    // {
    //   TH2D* slice = SliceWeights(weights3D,iBinZ,iBinZ);
    //   slice->SetStats(0);
    //   // if(!slice) continue;
    //
    //   canSlices->cd(iBinZ);
    //   // gPad->SetLogz();
    //   slice->DrawCopy("colz");
    // }
    // canSlices->SaveAs(Form("%s/weights_%s_all.pdf",sOutputPath.Data(), species[part].Data()), "pdf");

  } // end-for{species}

  fOutput->cd();
  listRun->Write("weights",TObject::kSingleKey);
  delete listRun;

  if(bRunByRun)
  {
    for(Short_t iRun = 0; iRun < iNumRuns; iRun++)
    {
      printf(" ### Run %d (out of %d)\n",iRun+1,iNumRuns);
      TFile* fileTemp = new TFile(Form("%s/merge/merge_%d/AnalysisResults.root",sPath.Data(),iRunList[iRun]),"READ");
      if(!fileTemp) { printf("Run %d | Input file with weights not found\n",iRunList[iRun]); return; }

      fileTemp->cd(sTaskTag.Data());

      TList* listTemp = (TList*) gDirectory->Get(Form("Flow_Weights_%s",sTaskTag.Data()));
      if(!listTemp) { printf("Run %d | TList with weights not found\n",iRunList[iRun]); gDirectory->ls() ;return; }

      TList* listRun = new TList();
      listRun->SetOwner(kTRUE);

      for(Short_t part = 0; part < iNumPart; part++)
      {
        printf(" -part %d (out of %d)\n",part+1,iNumPart);
        TH3D* h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
        if(!h3Weights) { printf("Run %d | Hist 'fh3Weights%s' not found\n",iRunList[iRun],species[part].Data()); listTemp->ls(); return; }

        TH2D* h2Weights = ProcessWeights2D(h3Weights);

        h2Weights->SetName(Form("%s",species[part].Data()));
        h2Weights->SetTitle(Form("%s | Run %d",species[part].Data(),iRunList[iRun]));
        listRun->Add(h2Weights);
      }

      fOutput->cd();
      listRun->Write(Form("%d",iRunList[iRun]),TObject::kSingleKey);
      delete listRun;
    }
  }

  return;
}
//_____________________________________________________________________________
TH2D* ProcessWeights2D(const TH3D* dist)
{
  if(!dist) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH2D* dist_2d = (TH2D*) dist->Project3D("xy");
  TH2D* weights_2d = (TH2D*) dist_2d->Clone();
  weights_2d->Reset();

  for(Short_t binEta(1); binEta < dist_2d->GetNbinsX()+1; ++binEta)
  {
    // projection onto phi (in eta bin)
    TH1D* proj = (TH1D*) dist_2d->ProjectionY(Form("y_%d",binEta), binEta,binEta);
    if(!proj) { printf("No projection! Something went wrong!\n"); return 0x0; }

    TH1D* hWeight = CalculateWeight(proj);

    for(Short_t binPhi(1); binPhi < weights_2d->GetNbinsY()+1; ++binPhi)
    {
      weights_2d->SetBinContent(binEta,binPhi, hWeight->GetBinContent(binPhi));
      weights_2d->SetBinError(binEta,binPhi, 0.0);
    }

  }

  return weights_2d;
}
//_____________________________________________________________________________
TH3D* ProcessWeights3D(const TH3D* dist)
{
  if(!dist) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH3D* distTemp = (TH3D*) dist->Clone("distTemp");
  // phi dist is on y-axis
  TH3D* weights = (TH3D*) dist->Clone("weights");
  weights->Reset();

  // preparing weigths:  bin content = max / bincontent

  // z-dimension
  for(Short_t binZ(1); binZ < dist->GetNbinsZ()+1; ++binZ)
  {
    // y-dimension : eta
    for(Short_t binEta(1); binEta < dist->GetNbinsX()+1; ++binEta)
    {
      distTemp->GetXaxis()->SetRange(binEta,binEta);
      distTemp->GetZaxis()->SetRange(binZ,binZ);

      // projection into phi
      TH1D* proj = (TH1D*) distTemp->Project3D(Form("y_%d_%d",binEta,binZ));
      if(!proj) { printf("No projection! Something went wrong!\n"); return 0x0; }

      // calculating weights in this 1D projection
      TH1D* hWeight = CalculateWeight(proj);
      for(Short_t binPhi(1); binPhi < hWeight->GetNbinsX()+1; ++binPhi)
      {
        weights->SetBinContent(binEta,binPhi,binZ, hWeight->GetBinContent(binPhi));
        weights->SetBinError(binEta,binPhi,binZ, 0.0);
      }

    }
  }

  return weights;
}
//_____________________________________________________________________________
TH1D* CalculateWeight(TH1D* proj)
{
  if(!proj) { printf("CalculateWeight-ERROR: Input proj not found!\n"); return 0x0; }

  TH1D* hWeight = (TH1D*) proj->Clone("weight");
  if(!hWeight) { printf("CalculateWeight-ERROR: hWeight not cloned!\n"); return 0x0; }

  Double_t dMax = proj->GetMaximum();
  // printf("dMax %f\n",dMax);

  for(Int_t iBin(0); iBin < proj->GetNbinsX()+2; ++iBin)
  {
    Double_t dContent = proj->GetBinContent(iBin);
    if(dContent > 0.0) { hWeight->SetBinContent(iBin, dMax / dContent); }
    else { hWeight->SetBinContent(iBin, 1.0); }
    hWeight->SetBinError(iBin, 0.0);
  }

  return hWeight;
}
//_____________________________________________________________________________
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh)
{
  if(!weights) { printf("PrepareWeights: no input histo\n"); return 0x0;}

  TH3D* weightsTemp = (TH3D*) weights->Clone();

  weightsTemp->GetZaxis()->SetRange(iBinLow,iBinHigh);

  TH2D* weights_2d = (TH2D*) weightsTemp->Project3D("xy");
  return weights_2d;
}
//_____________________________________________________________________________
