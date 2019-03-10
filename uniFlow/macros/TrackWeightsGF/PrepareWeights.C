/* PrepareWeights
 *
 * Macro for extracting the (eta,phi,z) distribution and make them into 3D (eta,phi,z) and 2D (eta,phi) weights
 * for Generic Framework calculations and preparing source ROOT file used for running the analysis.
 *
 * This can be done on run-integrated results and by run-by-run basis (if bRunByRun = kTRUE)
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

#include <iostream>
#include "TFile.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"

Bool_t fbDebug = kFALSE;
Bool_t fbUse3Dweights = kFALSE;
TString sTaskName = "UniFlow";
// TString sTag = "TPCcls90";
TString sTag = "FB768";
TString sPath = "/mnt/CodesALICE/";
// TString sOutputPath = "./" + sTag + "/";
TString sOutputPath = "./weights/";

TString sOutFileName = "weights.root";
// TString sOutputPath = sPath+"/weights"+sTag+"/";

// const Short_t iNumPart = 1; const TString species[iNumPart] = {"Refs"};
// const Short_t iNumPart = 5; const TString species[iNumPart] = {"Refs","Charged","K0s","Lambda","Phi"};
// const Short_t iNumPart = 5; const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton"};
const Short_t iNumPart = 6; const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton","Phi"};
// const Short_t iNumPart = 8; const TString species[iNumPart] = {"Refs","Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};

const Bool_t bRunByRun = kFALSE;
// Run Lists

// pPb 5.02 TeV
// RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20171129_v2.txt [31 runs]
// Int_t iNumRuns = 31; Int_t iRunList[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
// RunList_LHC16t_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [4 runs]
// Int_t iNumRuns = 4; Int_t iRunList[] = {267166, 267165, 267164, 267163};
// Int_t iNumRuns = 49; Int_t iRunList[] = {246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847,
// 246846, 246845, 246844, 246810, 246809, 246808, 246805, 246804, 246766,
// 246765 ,246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493,
// 246488, 246487, 246434, 246431,
//
// 246424, 246276, 246275, 246272, 246271, 246225,
// 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151};

Int_t iNumRuns = 77; Int_t iRunList[] = {
  246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847,
  246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766,
  246765 ,246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493,
  246488, 246487, 246434, 246431, 246424, 246276, 246275, 246272, 246271, 246225,
  246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151,
  246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042,
  246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923,
  245833, 245831, 245829, 245705, 245702, 245692, 245683
};


// pp 13 TeV
// RunList_LHC16l_pass1_CentralBarrelTracking_hadronPID_20170509_v2.txt [70]
// Int_t iNumRuns = 70; Int_t iRunList[] = {259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259713, 259711, 259705, 259704, 259703, 259700, 259697, 259668, 259650, 259649, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962, 258923, 258919};

// RunList_LHC16k_pass1_CentralBarrelTracking_hadronPID_20170516_v3.txt [194]
// Int_t iNumRuns = 194; Int_t iRunList[] = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};


TH2D* Process2DWeights(const TH2* dist);
TH2D* Process3DWeightsTo2D(const TH3D* dist);
TH3D* Process3DWeightsTo3D(const TH3D* dist);
TH1D* CalculateWeight(TH1D* proj);
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh);
TH1* MergeListProfiles(TList* list);

void Fatal(TString sMessage, TString sMethod = "");
void Error(TString sMessage, TString sMethod = "");
void Warning(TString sMessage, TString sMethod = "");
void Info(TString sMessage, TString sMethod = "");
void Debug(TString sMessage, TString sMethod = "");

void PrepareWeights()
{
    gSystem->mkdir(sOutputPath.Data(),kTRUE);

    TFile* fOutput = TFile::Open(Form("%s/%s",sOutputPath.Data(),sOutFileName.Data()),"UPDATE");
    if(!fOutput) { Error("Output opening failed!"); return; }

    if(!sTag.IsNull()) { sTaskName += sTag; }

    // check if "weights" list already exists
    TList* outList = nullptr;
    outList = (TList*) fOutput->Get("weights");
    if(outList) {
        Warning("Outlist already exists, using it instead of creating it now!");
    } else {
        outList = new TList();
        outList->SetOwner(kTRUE);
    }

    if(!bRunByRun) {

        TFile* fInput = TFile::Open(Form("%s/AnalysisResults.root",sPath.Data()),"READ");
        if(!fInput) { Error("Input file not opened!"); return; }
        if(!fInput->cd(sTaskName.Data())) { Error(Form("Directory named '%s' not found!",sTaskName.Data())); fInput->ls(); return; }

        TList* list =  (TList*) gDirectory->Get(Form("Weights_%s",sTaskName.Data()));
        if(!list) { Error("Input TList with weights not opened!"); return; }

        TString sListName = "averaged";
        if(!sTag.IsNull()) { sListName = sTag; }

        if(outList->FindObject(sListName.Data())) {
            Error(Form("This listRun '%s' already exits within output file. Terminating!",sListName.Data()));
            outList->ls();
            return;
        }

        TList* listRun = new TList();
        listRun->SetOwner(kTRUE);
        listRun->SetName(sListName.Data());

        for(Int_t part = 0; part < iNumPart; ++part) {
            Info(Form(" -part %d (out of %d)",part+1,iNumPart));

            TH1* weights = nullptr;

            if(fbUse3Dweights) {
                TH3D* h3Weights = (TH3D*) list->FindObject(Form("fh3Weights%s",species[part].Data()));
                if(!h3Weights) { Error(Form("Hist 'fh3Weights%s' not found",species[part].Data())); list->ls(); return; }

                weights = Process3DWeightsTo2D(h3Weights);
                if(!weights) { Error("Weights processing (3D->2D) failed!"); return; }
            } else {
                TH2D* h2Weights = (TH2D*) list->FindObject(Form("fh2Weights%s",species[part].Data()));
                if(!h2Weights) { Error(Form("Hist 'fh2Weights%s' not found",species[part].Data())); list->ls(); return; }

                weights = Process2DWeights(h2Weights);
                if(!weights) { Error("Weights processing (2D->2D) failed!"); return; }
            }

            // CheckWeights(h3Weights, iRunList[iRun]);
            weights->SetName(Form("%s",species[part].Data()));
            weights->SetTitle(Form("%s",species[part].Data()));
            listRun->Add(weights);
        } // end-for{species}

        outList->Add(listRun);

    } else {

        Warning("Not re-implemented!");

        TList* listGlob[iNumPart]{};

        for(Int_t part(0); part < iNumPart; ++part) {
            listGlob[part] = new TList();
        }

        for(Int_t iRun = 0; iRun < iNumRuns; ++iRun) {
            printf(" ### Run %d (out of %d)\n",iRun+1,iNumRuns);
            TFile* fileTemp = new TFile(Form("%s/merge/merge_%d/AnalysisResults.root",sPath.Data(),iRunList[iRun]),"READ");
            if(!fileTemp) { printf("Run %d | Input file with weights not found\n",iRunList[iRun]); return; }

            fileTemp->cd(sTaskName.Data());

            TList* listTemp = (TList*) gDirectory->Get(Form("Weights_%s",sTaskName.Data()));
            if(!listTemp) { printf("Run %d | TList with weights not found\n",iRunList[iRun]); gDirectory->ls() ;return; }

            TList* listRun = new TList();
            listRun->SetOwner(kTRUE);
            listRun->SetName(Form("%d",iRunList[iRun]));

            for(Int_t part = 0; part < iNumPart; ++part) {
                printf(" -part %d (out of %d)\n",part+1,iNumPart);
                TH3D* h3Weights = (TH3D*) listTemp->FindObject(Form("fh3Weights%s",species[part].Data()));
                if(!h3Weights) { printf("Run %d | Hist 'fh3Weights%s' not found\n",iRunList[iRun],species[part].Data()); listTemp->ls(); return; }

                TH2D* h2Weights = Process3DWeightsTo2D(h3Weights);
                h2Weights->SetName(Form("%s",species[part].Data()));
                h2Weights->SetTitle(Form("%s | Run %d",species[part].Data(),iRunList[iRun]));
                listRun->Add(h2Weights);

                listGlob[part]->Add(h3Weights);
            }

            outList->Add(listRun);
        }

        // making run-averaged weights per-species out of global lists
        TList* listRunAver = new TList();
        listRunAver->SetName("averaged");

        printf("=== Making run-averaged weights by merging individual runs. ===\n");
        for(Int_t part(0); part < iNumPart; ++part) {
            printf(" -part %d (out of %d)\n",part+1,iNumPart);

            TH3D* merged = (TH3D*) MergeListProfiles(listGlob[part]);
            if(!merged) { printf("Merge failed!"); return; }

            TH2D* h2Weights = Process3DWeightsTo2D(merged);
            h2Weights->SetName(Form("%s",species[part].Data()));
            h2Weights->SetTitle(Form("%s | Run-averaged",species[part].Data()));
            listRunAver->Add(h2Weights);

            delete merged;
            delete listGlob[part];
        }

        outList->Add(listRunAver);
    }

    fOutput->cd();
    outList->Write("weights",TObject::kSingleKey+TObject::kOverwrite);

    return;
}
//_____________________________________________________________________________
TH2D* Process2DWeights(const TH2* dist)
{
    if(!dist) { Error("No input histo","Process2DWeights"); return nullptr; }

    TH2D* dist_2d = (TH2D*) dist->Clone();

    TH2D* weights_2d = (TH2D*) dist->Clone();
    weights_2d->Reset();

    for(Int_t binEta(1); binEta < dist_2d->GetNbinsX()+1; ++binEta) {
        // projection onto phi (in eta bin)
        TH1D* proj = (TH1D*) dist_2d->ProjectionY(Form("y_%d",binEta), binEta,binEta);
        if(!proj) { printf("No projection! Something went wrong!\n"); return 0x0; }

        TH1D* hWeight = CalculateWeight(proj);

        for(Int_t binPhi(1); binPhi < weights_2d->GetNbinsY()+1; ++binPhi) {
            weights_2d->SetBinContent(binEta,binPhi, hWeight->GetBinContent(binPhi));
            weights_2d->SetBinError(binEta,binPhi, 0.0);
        }
    }

    return weights_2d;
}
//_____________________________________________________________________________
TH2D* Process3DWeightsTo2D(const TH3D* dist)
{
    if(!dist) { Error("No input histo","Process3DWeightsTo2D"); return nullptr; }

    TH2D* dist_2d = (TH2D*) dist->Project3D("xy");
    if(!dist_2d) { Error("Projection of 3D->2D failed!","Process3DWeightsTo2D"); return nullptr; }

    return Process2DWeights(dist_2d);
}
//_____________________________________________________________________________
TH3D* Process3DWeightsTo3D(const TH3D* dist)
{
    if(!dist) { Error("No input histo","Process3DWeightsTo3D"); return nullptr;}

    TH3D* distTemp = (TH3D*) dist->Clone("distTemp");
    // phi dist is on y-axis
    TH3D* weights = (TH3D*) dist->Clone("weights");
    weights->Reset();

    // preparing weigths:  bin content = max / bincontent

    // z-dimension
    for(Int_t binZ(1); binZ < dist->GetNbinsZ()+1; ++binZ) {
        // y-dimension : eta
        for(Int_t binEta(1); binEta < dist->GetNbinsX()+1; ++binEta) {
            distTemp->GetXaxis()->SetRange(binEta,binEta);
            distTemp->GetZaxis()->SetRange(binZ,binZ);

            // projection into phi
            TH1D* proj = (TH1D*) distTemp->Project3D(Form("y_%d_%d",binEta,binZ));
            if(!proj) { Error("No projection! Something went wrong!","Process3DWeightsTo3D"); return nullptr; }

            // calculating weights in this 1D projection
            TH1D* hWeight = CalculateWeight(proj);

            for(Int_t binPhi(1); binPhi < hWeight->GetNbinsX()+1; ++binPhi) {
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
    if(!proj) { Error("Input proj not found!","CalculateWeight"); return nullptr; }

    TH1D* hWeight = (TH1D*) proj->Clone("weight");
    if(!hWeight) { Error("hWeight not cloned!","CalculateWeight"); return nullptr; }

    Double_t dMax = proj->GetMaximum();
    // printf("dMax %f\n",dMax);

    for(Int_t iBin(0); iBin < proj->GetNbinsX()+2; ++iBin) {
        Double_t dContent = proj->GetBinContent(iBin);
        if(dContent > 0.0) {
            hWeight->SetBinContent(iBin, dMax / dContent);
        } else {
            hWeight->SetBinContent(iBin, 1.0);
        }
        hWeight->SetBinError(iBin, 0.0);
    }

    return hWeight;
}
//_____________________________________________________________________________
TH2D* SliceWeights(const TH3D* weights, Int_t iBinLow, Int_t iBinHigh)
{
    if(!weights) { Error("No input histo","SliceWeights"); return nullptr; }

    TH3D* weightsTemp = (TH3D*) weights->Clone();

    weightsTemp->GetZaxis()->SetRange(iBinLow,iBinHigh);

    TH2D* weights_2d = (TH2D*) weightsTemp->Project3D("xy");
    return weights_2d;
}
//_____________________________________________________________________________
TH1* MergeListProfiles(TList* list)
{
    // merge list of TProfiles into single TProfile and return it
    if(!list || list->IsEmpty()) { Error("List not valid or empty","MergeListProfiles"); return nullptr; }

    TH1* merged = (TH1*) list->At(0)->Clone();
    // merged->SetName(Form("%s_merged",merged->GetName()));

    // only 1 entry
    if(list->GetEntries() < 2) {
        Warning("Only one entry for merging; returning it directly instead!","MergeListProfiles");
        return merged;
    }

    merged->Reset();
    Double_t mergeStatus = merged->Merge(list);
    if(mergeStatus == -1) { Error("Merging failed!","MergeListProfiles"); return nullptr; }

    return merged;
}
//_____________________________________________________________________________
void Fatal(TString sMsg, TString sMethod)
{
    if(!sMethod.IsNull()) { sMethod = Form("::%s",sMethod.Data()); }
	printf("\033[91mFatal%s  %s. Terminating!\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void Error(TString sMsg, TString sMethod)
{
    if(!sMethod.IsNull()) { sMethod = Form("::%s",sMethod.Data()); }
	printf("\033[91mError%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void Info(TString sMsg, TString sMethod)
{
    if(!sMethod.IsNull()) { sMethod = Form("::%s",sMethod.Data()); }
	printf("\033[96mInfo%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void Warning(TString sMsg, TString sMethod)
{
    if(!sMethod.IsNull()) { sMethod = Form("::%s",sMethod.Data()); }
	printf("\033[93mWarning%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void Debug(TString sMsg, TString sMethod)
{
    if(!fbDebug) { return; }
    if(!sMethod.IsNull()) { sMethod = Form("::%s",sMethod.Data()); }
	printf("\033[95mDebug%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
