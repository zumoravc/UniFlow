#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TMath.h"

TList* flFlowRefs;
TList* flFlowCharged;
TList* flFlowPID;
TList* flFlowPhi;
TList* flFlowK0s;
TList* flFlowLambda;

// Double_t fdMultBinning[] = {0,100};
Double_t fdMultBinning[] = {0,10,20,40,60,100};
Int_t fiNumMultBins = sizeof(fdMultBinning) / sizeof(fdMultBinning[0]) - 1;
Double_t fdPtBinning[] = {0.2,0.4,0.6,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.0};
Int_t fiNumPtBins = sizeof(fdPtBinning) / sizeof(fdPtBinning[0]) - 1;

Bool_t bUsePtBinning = kTRUE;
Bool_t bMergePosNeg = kTRUE;

TString sInputFile = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/AnalysisResults.root";
TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/study-slices/rebin/";
TString sTaskName = "UniFlow";
TString sGap = "gap04";
TString sSpecies = "K0s";
// TString sSpecies = "Lambda";
// TString sSpecies = "Phi";

Bool_t LoadLists(TFile* ffInputFile, TString fsTaskName);
Bool_t PrepareSlices(TH3D* h3Entries, TProfile3D* p3Flow);
TProfile2D* Project3DProfile(const TProfile3D* prof3dorig = 0x0); // making projection out of TProfile3D
TProfile2D* DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,bool originalRange, bool useUF, bool useOF);
TH2D*       DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY, bool computeErrors, bool originalRange, bool useUF, bool useOF);

void StudySlices()
{

  sOutputPath = sOutputPath + "/" + sGap + "/" + sSpecies + "/";

  TFile* fileInput = TFile::Open(sInputFile.Data(),"READ");
  if(!LoadLists(fileInput,sTaskName)) return;

  TList* list = 0x0;
  TString sNameMass;
  TString sNameFlow;

  if(sSpecies.EqualTo("K0s"))
  {
    list = flFlowK0s;
    sNameMass = "fh3V0sEntries"+sSpecies;
    sNameFlow = "fp3V0sCorr"+sSpecies;
  }
  else if(sSpecies.EqualTo("Lambda"))
  {
    list = flFlowLambda;
    sNameMass = "fh3V0sEntries"+sSpecies;
    sNameFlow = "fp3V0sCorr"+sSpecies;
  }
  else if(sSpecies.EqualTo("Phi"))
  {
    list = flFlowPhi;
    sNameMass = "fh3PhiEntriesSignal";
    sNameFlow = "fp3PhiCorr";
  }

  list->ls();

  TH3D* h3EntriesPos = (TH3D*) list->FindObject(Form("%s_%s_Pos",sNameMass.Data(),sGap.Data()));
  if(!h3EntriesPos) { printf("E-StudySlices: h3EntriesPos not found!\n"); return; }

  TH3D* h3EntriesNeg = (TH3D*) list->FindObject(Form("%s_%s_Neg",sNameMass.Data(),sGap.Data()));
  if(!h3EntriesNeg) { printf("E-StudySlices: h3EntriesNeg not found!\n"); return; }

  TProfile3D* p3FlowPos = (TProfile3D*) list->FindObject(Form("%s_<2>_harm2_%s_Pos",sNameFlow.Data(),sGap.Data()));
  if(!p3FlowPos) { printf("E-StudySlices: p3FlowPos not found!\n"); return; }

  TProfile3D* p3FlowNeg = (TProfile3D*) list->FindObject(Form("%s_<2>_harm2_%s_Neg",sNameFlow.Data(),sGap.Data()));
  if(!p3FlowNeg) { printf("E-StudySlices: p3FlowNeg not found!\n"); return; }

  TH3D* h3Entries = 0x0;
  TProfile3D* p3Flow = 0x0;

  if(bMergePosNeg)
  {
    h3Entries = (TH3D*) h3EntriesPos->Clone("h3Entries");
    h3Entries->Add(h3EntriesNeg);

    p3Flow = (TProfile3D*) p3FlowPos->Clone("p3Flow");
    p3Flow->Reset();
    TList* listMerge = new TList();
    listMerge->Add(p3FlowPos);
    listMerge->Add(p3FlowNeg);
    Double_t dStatus = p3Flow->Merge(listMerge);
    if(dStatus == -1.0) { printf("E-PrepareSlices : Merged failed \n"); return; }
  }
  else
  {
    h3Entries = h3EntriesPos;
    p3Flow = p3FlowPos;
  }

  // h3Entries->Draw();
  p3Flow->Draw();

  if(!PrepareSlices(h3Entries,p3Flow)) { printf("E-StudySlices: Slices not prepared"); }

  return;
}
//_____________________________________________________________________________
Bool_t LoadLists(TFile* ffInputFile, TString fsTaskName)
{
  // loading TLists into task
  if(!ffInputFile) { printf("Input file does not exists!\n"); return kFALSE; }
  ffInputFile->cd(fsTaskName.Data());

  flFlowRefs = (TList*) gDirectory->Get(Form("Flow_Refs_%s",fsTaskName.Data()));
  if(!flFlowRefs) { printf("flFlow_Refs list does not exists!\n"); ffInputFile->ls(); return kFALSE; }
  flFlowCharged = (TList*) gDirectory->Get(Form("Flow_Charged_%s",fsTaskName.Data()));
  if(!flFlowCharged) { printf("flFlow_Charged list does not exists!\n"); return kFALSE; }
  flFlowPID = (TList*) gDirectory->Get(Form("Flow_PID_%s",fsTaskName.Data()));
  if(!flFlowPID) { printf("flFlow_PID list does not exists!\n"); return kFALSE; }
  flFlowPhi = (TList*) gDirectory->Get(Form("Flow_Phi_%s",fsTaskName.Data()));
  if(!flFlowPhi) { printf("flFlow_Phi list does not exists!\n"); return kFALSE; }
  flFlowK0s = (TList*) gDirectory->Get(Form("Flow_K0s_%s",fsTaskName.Data()));
  if(!flFlowK0s) { printf("flFlow_K0s list does not exists!\n"); return kFALSE; }
  flFlowLambda = (TList*) gDirectory->Get(Form("Flow_Lambda_%s",fsTaskName.Data()));
  if(!flFlowLambda) { printf("flFlow_Lambda list does not exists!\n"); return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t PrepareSlices(TH3D* h3Entries, TProfile3D* p3Flow)
{
  if(!h3Entries) { printf("E-PrepareSlices : h3Entries does not exists!\n"); return kFALSE; }
  if(!p3Flow) { printf("E-PrepareSlices : p3Flow does not exists!\n"); return kFALSE; }

  Int_t iNumPtBins = h3Entries->GetNbinsY();
  if(bUsePtBinning)
  {
    if(fiNumPtBins < 1) { printf("E-PrepareSlices: Low number of pt bins\n"); return kFALSE; }
    iNumPtBins = fiNumPtBins;
  }

  TList* listFlow[fiNumMultBins];
  TList* listMass[fiNumMultBins];

  for(Int_t multBin(0); multBin < fiNumMultBins; ++multBin)
  {
    gSystem->mkdir(Form("%s/slices_mult%d/",sOutputPath.Data(), multBin),kTRUE);

    Int_t binMultLow = h3Entries->GetXaxis()->FindFixBin(fdMultBinning[multBin]);
    Int_t binMultHigh = h3Entries->GetXaxis()->FindFixBin(fdMultBinning[multBin+1]) - 1;

    TH1D* hNumEntriesMult = (TH1D*) h3Entries->ProjectionY(Form("hNumEntriesMult%d",multBin),binMultLow,binMultHigh,1,1);
    hNumEntriesMult->Reset();

    listFlow[multBin] = new TList();
    listMass[multBin] = new TList();

    for(Int_t ptBin(0); ptBin < iNumPtBins; ++ptBin)
    {
      Int_t binPtLow = ptBin+1;
      Int_t binPtHigh = ptBin+1;

      if(bUsePtBinning)
      {
        binPtLow = h3Entries->GetYaxis()->FindFixBin(fdPtBinning[ptBin]);
        binPtHigh = h3Entries->GetYaxis()->FindFixBin(fdPtBinning[ptBin+1]) - 1;
      }

      TH1D* hEntries = (TH1D*) h3Entries->ProjectionZ(Form("%s_pz_mult%d_pt%d",h3Entries->GetName(),multBin,ptBin), binMultLow,binMultHigh, binPtLow,binPtHigh,"e");
      if(!hEntries) { printf("E-PrepareSlices: projection failed\n"); return kFALSE; }

      p3Flow->GetXaxis()->SetRange(binMultLow,binMultHigh);
      TProfile2D* p2Flow = Project3DProfile(p3Flow);
      TProfile* p1Flow = (TProfile*) p2Flow->ProfileX(Form("profFlowMass_cent%d_pt%d",multBin,ptBin),binPtLow,binPtHigh);
      TH1D* hFlow = (TH1D*) p1Flow->ProjectionX(Form("hFlowMass_cent%d_pt%d",multBin,ptBin));
      // TH1D* hFlow = (TH1D*) p3Flow->ProjectionZ(Form("%s_pz_mult%d_pt%d",p3Flow->GetName(),multBin,ptBin), binMultLow,binMultHigh, binPtLow,binPtHigh);
      if(!hFlow) { printf("E-PrepareSlices: projection failed\n"); return kFALSE; }


      hNumEntriesMult->SetBinContent(ptBin+1, hEntries->GetEntries());

      listFlow[multBin]->Add(hFlow);
      listMass[multBin]->Add(hEntries);

      TCanvas* canTemp = new TCanvas("canTemp","canTemp",1600,800);
      canTemp->Divide(2,1);
      canTemp->cd(1);
      canTemp->SetLogy();
      hEntries->DrawCopy();
      canTemp->cd(2);
      hFlow->DrawCopy();

      canTemp->SaveAs(Form("%s/slices_mult%d/slice_pt%d.pdf",sOutputPath.Data(),multBin,ptBin),"pdf");
    }

    TCanvas* can = new TCanvas("can","can",800,800);
    can->cd();
    gPad->SetLogy();
    hNumEntriesMult->Draw();
    can->SaveAs(Form("%s/hEntriesMult%d.pdf",sOutputPath.Data(),multBin),"pdf");
  }

  gSystem->mkdir(Form("%s/slices_multall/",sOutputPath.Data()),kTRUE);

  TCanvas* canOver = new TCanvas("canOver","canOver",fiNumMultBins*400,800);
  canOver->Divide(fiNumMultBins,2);

  for(Int_t ptBin(0); ptBin < iNumPtBins; ++ptBin)
  {
    for(Int_t multBin(0); multBin < fiNumMultBins; ++multBin)
    {
      canOver->cd(multBin+1);
      ( (TH1D*)listMass[multBin]->At(ptBin))->DrawCopy();

      canOver->cd(multBin+fiNumMultBins+1);
      ( (TH1D*)listFlow[multBin]->At(ptBin))->DrawCopy();

    }
    canOver->SaveAs(Form("%s/slices_multall/slice_pt%d.pdf",sOutputPath.Data(),ptBin),"pdf");
  }
  return kTRUE;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
TProfile2D* Project3DProfile(const TProfile3D* prof3dorig)
{
  if(!prof3dorig) return 0x0;
  TProfile3D* prof3d = (TProfile3D*) prof3dorig->Clone();

  Int_t iBinFirst = prof3d->GetXaxis()->GetFirst();
  Int_t iBinLast = prof3d->GetXaxis()->GetLast();
  Int_t iNumBins = prof3d->GetNbinsX();
  Int_t iNumBinsAxis = prof3d->GetXaxis()->GetNbins();
  // printf("Bins:  %d - %d (%d | %d) \n", iBinFirst,iBinLast,iNumBins,iNumBinsAxis);

  // // making 3d hist from 3d profile
  // TH3D* hist3d = prof3d->ProjectionXYZ();   //NOTE do not care about range !!!
  // TH3D* hist3d_entry = prof3d->ProjectionXYZ("hist3d_entry","B");   //NOTE do not care about range !!!
  // TH3D* hist3d_weight = prof3d->ProjectionXYZ("hist3d_weight","W");   //NOTE do not care about range !!!

  TProfile2D* prof2d_test = DoProjectProfile2D(prof3d,"prof2d_test","",prof3d->GetYaxis(),prof3d->GetZaxis(),1,0,0);
  // prof2d_test->Draw("colz");

  return prof2d_test;

  // resulting profile
  // TProfile* result = new TProfile("result","result",100,0,10);
  // for(Int_t i(0); i < 10; i++) result->Fill(i,1);

  //
  // TCanvas* canTest = new TCanvas("canTest");
  // canTest->Divide(2,2);
  // canTest->cd(1);
  // prof3d->Draw("box");
  // canTest->cd(2);
  // hist3d->Draw("box");
  // canTest->cd(3);
  // hist3d_entry->Draw("box");
  // canTest->cd(4);
  // hist3d_weight->Draw("box");


  // return result;
}
//_____________________________________________________________________________
TProfile2D * DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,
                                           bool originalRange, bool useUF, bool useOF)
{
// internal method to project to a 2D Profile
 // called from TH3::Project3DProfile but re-implemented in case of the TPRofile3D since what is done is different

 // projX, projY: axes of the orifinal histogram to which the projection is done (e.g. xy)

 // Get the ranges where we will work.
 Int_t ixmin = projX->GetFirst();
 Int_t ixmax = projX->GetLast();
 Int_t iymin = projY->GetFirst();
 Int_t iymax = projY->GetLast();
 if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
 if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
 Int_t nx = ixmax-ixmin+1;
 Int_t ny = iymax-iymin+1;

 // Create the projected profiles
 TProfile2D *p2 = 0;
 // Create always a new TProfile2D (not as in the case of TH3 projection)

 const TArrayD *xbins = projX->GetXbins();
 const TArrayD *ybins = projY->GetXbins();
 // assume all axis have variable bins or have fixed bins
 if ( originalRange ) {
    if (xbins->fN == 0 && ybins->fN == 0) {
       p2 = new TProfile2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                           ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
    } else {
       p2 = new TProfile2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
    }
 } else {
    if (xbins->fN == 0 && ybins->fN == 0) {
       p2 = new TProfile2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                           ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
    } else {
       p2 = new TProfile2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
    }
 }

 // new profile p2 is set according to axis ranges (keeping originals or not)

 // weights
 bool useWeights = (h3->fBinSumw2.fN != 0); //array elements
 if (useWeights) p2->Sumw2();

 // make projection in a 3D first // from 3D profile -> TH3
 TH3D * h3dW = h3->ProjectionXYZ("h3temp-W","W"); // getbincontent*getBinentries
 TH3D * h3dN = h3->ProjectionXYZ("h3temp-N","B"); // bin content is original profile = GetEntriesBin

 // fix ???
 h3dW->GetXaxis()->SetRange(h3->GetXaxis()->GetFirst(),h3->GetXaxis()->GetLast());
 h3dW->GetYaxis()->SetRange(h3->GetYaxis()->GetFirst(),h3->GetYaxis()->GetLast());
 h3dW->GetZaxis()->SetRange(h3->GetZaxis()->GetFirst(),h3->GetZaxis()->GetLast());
 h3dN->GetXaxis()->SetRange(h3->GetXaxis()->GetFirst(),h3->GetXaxis()->GetLast());
 h3dN->GetYaxis()->SetRange(h3->GetYaxis()->GetFirst(),h3->GetYaxis()->GetLast());
 h3dN->GetZaxis()->SetRange(h3->GetZaxis()->GetFirst(),h3->GetZaxis()->GetLast());




 h3dW->SetDirectory(0); h3dN->SetDirectory(0); // istograms does not bellow to any directorz ???

 // note that h3dW is always a weighted histogram - so we need to compute error in the projection
 TAxis * projX_hW = h3dW->GetXaxis();
 TAxis * projX_hN = h3dN->GetXaxis();
 if (projX == h3->GetYaxis() ) {  projX_hW =  h3dW->GetYaxis();  projX_hN =  h3dN->GetYaxis(); }
 if (projX == h3->GetZaxis() ) {  projX_hW =  h3dW->GetZaxis();  projX_hN =  h3dN->GetZaxis(); }
 TAxis * projY_hW = h3dW->GetYaxis();
 TAxis * projY_hN = h3dN->GetYaxis();
 if (projY == h3->GetXaxis() ) {  projY_hW =  h3dW->GetXaxis();  projY_hN =  h3dN->GetXaxis(); }
 if (projY == h3->GetZaxis() ) {  projY_hW =  h3dW->GetZaxis();  projY_hN =  h3dN->GetZaxis(); }
 // checking the axes

 // TH3 -> TH2
 TH2D * h2W = DoProject2D(h3dW,"htemp-W","",projX_hW, projY_hW, true, originalRange, useUF, useOF);
 TH2D * h2N = DoProject2D(h3dN,"htemp-N","",projX_hN, projY_hN, useWeights, originalRange, useUF, useOF);
 h2W->SetDirectory(0); h2N->SetDirectory(0);


 // fill the bin content
 R__ASSERT( h2W->fN == p2->fN );
 R__ASSERT( h2N->fN == p2->fN );
 R__ASSERT( h2W->GetSumw2()->fN != 0); // h2W should always be a weighted histogram since h3dW is weighted


 // filling the new tprofile2D
 for (int i = 0; i < p2->fN ; ++i) {
    //std::cout << " proj bin " << i << "  " <<  h2W->fArray[i] << "  " << h2N->fArray[i] << std::endl;
    p2->fArray[i] = h2W->fArray[i];   // array of profile is sum of all values
    p2->GetSumw2()->fArray[i]  = h2W->GetSumw2()->fArray[i];   // array of content square of profile is weight square of the W projected histogram
    p2->SetBinEntries(i, h2N->fArray[i] );
    if (useWeights) p2->GetBinSumw2()->fArray[i] = h2N->GetSumw2()->fArray[i];    // sum of weight squares are stored to compute errors in h1N histogram
 }
 // delete the created histograms
 delete h3dW;
 delete h3dN;
 delete h2W;
 delete h2N;

 // Also we need to set the entries since they have not been correctly calculated during the projection
 // we can only set them to the effective entries
 p2->SetEntries( p2->GetEffectiveEntries() );

 return p2;
}
//_____________________________________________________________________________
TH2D* DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY,
                    bool computeErrors, bool originalRange,
                    bool useUF, bool useOF)
{
  // internal method performing the projection to a 2D histogram
     // called from TH3::Project3D

     TH2D *h2 = 0;

     // Get range to use as well as bin limits
     Int_t ixmin = projX->GetFirst();
     Int_t ixmax = projX->GetLast();
     Int_t iymin = projY->GetFirst();
     Int_t iymax = projY->GetLast();
     if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
     if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
     Int_t nx = ixmax-ixmin+1;
     Int_t ny = iymax-iymin+1;

      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange )
      {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1]
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
         }
      } else {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1]
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
         }
      }

    //  // Copy the axis attributes and the axis labels if needed.
    //  THashList* labels1 = 0;
    //  THashList* labels2 = 0;
    //  // "xy"
    //  h2->GetXaxis()->ImportAttributes(projY);
    //  h2->GetYaxis()->ImportAttributes(projX);
    //  labels1 = projY->GetLabels();
    //  labels2 = projX->GetLabels();
    //  if (labels1) {
    //     TIter iL(labels1);
    //     TObjString* lb;
    //     Int_t i = 1;
    //     while ((lb=(TObjString*)iL())) {
    //        h2->GetXaxis()->SetBinLabel(i,lb->String().Data());
    //        i++;
    //     }
    //  }
    //  if (labels2) {
    //     TIter iL(labels2);
    //     TObjString* lb;
    //     Int_t i = 1;
    //     while ((lb=(TObjString*)iL())) {
    //        h2->GetYaxis()->SetBinLabel(i,lb->String().Data());
    //        i++;
    //     }
    //  }
    //  h2->SetLineColor(this->GetLineColor());
    //  h2->SetFillColor(this->GetFillColor());
    //  h2->SetMarkerColor(this->GetMarkerColor());
    //  h2->SetMarkerStyle(this->GetMarkerStyle());

     // Activate errors
     if ( computeErrors) h2->Sumw2();

     // Set references to the axis, so that the bucle has no branches.
     TAxis* out = 0;
     if ( projX != h3->GetXaxis() && projY != h3->GetXaxis() ) {
        out = h3->GetXaxis();
     } else if ( projX != h3->GetYaxis() && projY != h3->GetYaxis() ) {
        out = h3->GetYaxis();
     } else {
        out = h3->GetZaxis();
     }

     Int_t *refX = 0, *refY = 0, *refZ = 0;
     Int_t ixbin, iybin, outbin;
     if ( projX == h3->GetXaxis() && projY == h3->GetYaxis() ) { refX = &ixbin;  refY = &iybin;  refZ = &outbin; }
     if ( projX == h3->GetYaxis() && projY == h3->GetXaxis() ) { refX = &iybin;  refY = &ixbin;  refZ = &outbin; }
     if ( projX == h3->GetXaxis() && projY == h3->GetZaxis() ) { refX = &ixbin;  refY = &outbin; refZ = &iybin;  }
     if ( projX == h3->GetZaxis() && projY == h3->GetXaxis() ) { refX = &iybin;  refY = &outbin; refZ = &ixbin;  }
     if ( projX == h3->GetYaxis() && projY == h3->GetZaxis() ) { refX = &outbin; refY = &ixbin;  refZ = &iybin;  }
     if ( projX == h3->GetZaxis() && projY == h3->GetYaxis() ) { refX = &outbin; refY = &iybin;  refZ = &ixbin;  }
     R__ASSERT (refX != 0 && refY != 0 && refZ != 0);

     // Fill the projected histogram excluding underflow/overflows if considered in the option
     // if specified in the option (by default they considered)
     Double_t totcont  = 0;

     Int_t outmin = out->GetFirst();
     Int_t outmax = out->GetLast();
     // GetFirst(), GetLast() can return (0,0) when the range bit is set artifically (see TAxis::SetRange)
     if (outmin == 0 && outmax == 0) { outmin = 1; outmax = out->GetNbins(); }
     // correct for underflow/overflows
     if (useUF && !out->TestBit(TAxis::kAxisRange) )  outmin -= 1;
     if (useOF && !out->TestBit(TAxis::kAxisRange) )  outmax += 1;

     for (ixbin=0;ixbin<=1+projX->GetNbins();ixbin++){
        if ( projX->TestBit(TAxis::kAxisRange) && ( ixbin < ixmin || ixbin > ixmax )) continue;
        Int_t ix = h2->GetYaxis()->FindBin( projX->GetBinCenter(ixbin) );

        for (iybin=0;iybin<=1+projY->GetNbins();iybin++){
           if ( projY->TestBit(TAxis::kAxisRange) && ( iybin < iymin || iybin > iymax )) continue;
           Int_t iy = h2->GetXaxis()->FindBin( projY->GetBinCenter(iybin) );

           Double_t cont = 0;
           Double_t err2 = 0;

           // loop on the bins to be integrated (outbin should be called inbin)
           for (outbin = outmin; outbin <= outmax; outbin++){

              Int_t bin = h3->GetBin(*refX,*refY,*refZ);

              // sum the bin contents and errors if needed
              cont += h3->GetBinContent(bin);
              if (computeErrors) {
                 Double_t exyz = h3->GetBinError(bin);
                 err2 += exyz*exyz;
              }

           }

           // remember axis are inverted
           h2->SetBinContent(iy , ix, cont);
           if (computeErrors) h2->SetBinError(iy, ix, TMath::Sqrt(err2) );
           // sum all content
           totcont += cont;

        }
     }

     // since we use fill we need to reset and recalculate the statistics (see comment in DoProject1D )
     // or keep original statistics if consistent sumw2
     bool resetStats = true;
     double eps = 1.E-12;
     if (h3->IsA() == TH3F::Class() ) eps = 1.E-6;
     if (h3->fTsumw != 0 && TMath::Abs( h3->fTsumw - totcont) <  TMath::Abs(h3->fTsumw) * eps) resetStats = false;

     bool resetEntries = resetStats;
     // entries are calculated using underflow/overflow. If excluded entries must be reset
     resetEntries |= !useUF || !useOF;

     if (!resetStats) {
        Double_t stats[TH1::kNstat];
        Double_t oldst[TH1::kNstat]; // old statistics
        for (Int_t i = 0; i < TH1::kNstat; ++i) { oldst[i] = 0; }
        h3->GetStats(oldst);
        std::copy(oldst,oldst+TH1::kNstat,stats);
        // not that projX refer to Y axis and projX refer to the X axis of projected histogram
        // nothing to do for projection in Y vs X
        if ( projY == h3->GetXaxis() && projX == h3->GetZaxis() ) {  // case XZ
           stats[4] = oldst[7];
           stats[5] = oldst[8];
           stats[6] = oldst[9];
        }
        if ( projY == h3->GetYaxis() ) {
           stats[2] = oldst[4];
           stats[3] = oldst[5];
           if ( projX == h3->GetXaxis() )  { // case YX
              stats[4] = oldst[2];
              stats[5] = oldst[3];
           }
           if ( projX == h3->GetZaxis() )  { // case YZ
              stats[4] = oldst[7];
              stats[5] = oldst[8];
              stats[6] = oldst[10];
           }
        }
        else if  ( projY == h3->GetZaxis() ) {
           stats[2] = oldst[7];
           stats[3] = oldst[8];
           if ( projX == h3->GetXaxis() )  { // case ZX
              stats[4] = oldst[2];
              stats[5] = oldst[3];
              stats[6] = oldst[9];
           }
           if ( projX == h3->GetYaxis() )  { // case ZY
              stats[4] = oldst[4];
              stats[5] = oldst[5];
              stats[6] = oldst[10];
           }
        }
        // set the new statistics
        h2->PutStats(stats);
     }
     else {
        // recalculate the statistics
        h2->ResetStats();
     }

     if (resetEntries) {
        // use the effective entries for the entries
        // since this  is the only way to estimate them
        Double_t entries =  h2->GetEffectiveEntries();
        if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
        h2->SetEntries( entries );
     }
     else {
        h2->SetEntries( h3->fEntries );
     }


     return h2;
}
//_____________________________________________________________________________
