// plotting macro
// v_2 from two- and multi-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void elFlow()
{
  TMultiGraph *mg = new TMultiGraph();


  Double_t xCross[10] = {2.5,7.5,15,25,35,45,55,65,75};
  Double_t xCrossErr[10] = {2.5,2.5,5,5,5,5,5,5,5};
  // RUN2 data
  Double_t v22Gap10Run2[10]={0.0283859, 0.0456604, 0.0655068, 0.0870721, 0.099105, 0.104143, 0.10286, 0.0974591, 0.0888104};
  Double_t v22Gap10Run2Err[10]={0.000570139, 0.000643862, 0.000373923, 0.000444072, 0.000553776, 0.000732208, 0.00107352, 0.00186236, 0.00437842};
  Double_t v22Gap10Run2Sys[10]={0.000425788, 0.000684906, 0.000982603, 0.00130608, 0.00148658, 0.00156214, 0.00154291, 0.00146189, 0.00133216};
  Double_t v22Gap10Run2CombErr[10]={0.000711585, 0.000940029, 0.00105135, 0.00137951, 0.00158637, 0.00172523, 0.00187963, 0.00236759, 0.0045766};

  Double_t v24Run2[10]={-1, 0.0366982, 0.0586746, 0.0779774, 0.0867288, 0.0890264, 0.0835125, 0.0806163, 0.0464366};
  Double_t v24Run2Err[10]={0, 0.00258022, 0.000331448, 0.00138072, 0.00257394, 0.00275397, 0.00207766, 0.00961209, 0.0835755};
  Double_t v24Run2Sys[10]={0, 0.000623869, 0.000997468, 0.00132562, 0.00147439, 0.00151345, 0.00141971, 0.00137048, 0.000789422};
  Double_t v24Run2CombErr[10]={0, 0.00265457, 0.00105109, 0.00191407, 0.00296631, 0.00314243, 0.0025164, 0.0097093, 0.0835792};

  Double_t v26Run2[10]={-1, 0.0361417, 0.0581481, 0.0774404, 0.086501, 0.0882722, 0.0827871, 0.0685936};
  Double_t v26Run2Err[10]={0, 0.00312685, 0.00031263, 0.00125004, 0.00244129, 0.00239164, 0.00455578, 0.0183562};
  Double_t v26Run2Sys[10]={0, 0.000614409, 0.000988518, 0.00131649, 0.00147052, 0.00150063, 0.00140738, 0.00116609};
  Double_t v26Run2CombErr[10]={0, 0.00318664, 0.00103678, 0.00181542, 0.00284997, 0.00282345, 0.00476821, 0.0183932};

  Double_t v28Run2[10]={-1, 0.0364682, 0.0580561, 0.0773885, 0.0866066, 0.0885098, 0.0826735, 0.0646657};
  Double_t v28Run2Err[10]={0, 0.00297485, 0.000345201, 0.00127119, 0.00241142, 0.00223079, 0.00610433, 0.039863};
  Double_t v28Run2Sys[10]={0, 0.00182341, 0.0029028, 0.00386943, 0.00433033, 0.00442549, 0.00413367, 0.00323328};
  Double_t v28Run2CombErr[10]={0, 0.0034892, 0.00292326, 0.00407288, 0.00495648, 0.00495595, 0.00737225, 0.0399939};

  TGraphErrors *graphv22Gap10 = new TGraphErrors(9,xCross,v22Gap10Run2,xCrossErr,v22Gap10Run2CombErr);
  graphv22Gap10->SetMarkerStyle(kFullSquare);
  graphv22Gap10->SetMarkerColor(kRed+1);
  graphv22Gap10->SetMarkerSize(1.);
  graphv22Gap10->SetLineColor(kRed+1);
  mg->Add(graphv22Gap10);
  //ShiftAlongXaxis(graphv22Gap10,1.0);

  TGraphErrors *graph_v24 = new TGraphErrors(8,xCross,v24Run2,xCrossErr,v24Run2CombErr);
  graph_v24->SetMarkerStyle(kFullCircle);
  graph_v24->SetMarkerColor(kBlue);
  graph_v24->SetLineColor(kBlue);
  graph_v24->SetMarkerSize(1.2);
  mg->Add(graph_v24);
  //ShiftAlongXaxis(graph_v24,-1);

  TGraphErrors *graph_v26 = new TGraphErrors(8,xCross,v26Run2,xCrossErr,v26Run2CombErr);
  graph_v26->SetMarkerStyle(33);
  graph_v26->SetMarkerColor(kGreen+1);
  graph_v26->SetMarkerSize(1.6);
  graph_v26->SetLineColor(kGreen+1);
  mg->Add(graph_v26);
  //ShiftAlongXaxis(graph_v26,-1);

  TGraphErrors *graph_v28 = new TGraphErrors(8,xCross,v28Run2,xCrossErr,v28Run2CombErr);
  graph_v28->SetMarkerStyle(3);
  graph_v28->SetMarkerColor(kBlack);
  graph_v28->SetLineColor(kBlack);
  graph_v28->SetMarkerSize(1.5);
  mg->Add(graph_v28);
  //ShiftAlongXaxis(graph_v28,1);


  TFile* fileIn = TFile::Open("/home/alidock/ana/output/LHC15o/train_7302/processUniFlow/Processed.root","READ");

  TGraphErrors *data_v22Gap10 = new TGraphErrors((TH1D*) fileIn->Get("Refs_hFlow2_harm2_gap10"));
  data_v22Gap10->SetMarkerStyle(kOpenSquare);
  data_v22Gap10->SetMarkerColor(kRed);
  data_v22Gap10->SetMarkerSize(1.);
  data_v22Gap10->SetLineColor(kRed);
  mg->Add(data_v22Gap10);

  TGraphErrors *data_v24 = new TGraphErrors((TH1D*) fileIn->Get("Refs_hFlow4_harm2_gap-10"));
  data_v24->SetMarkerStyle(kOpenCircle);
  data_v24->SetMarkerColor(kBlue);
  data_v24->SetMarkerSize(1.);
  data_v24->SetLineColor(kBlue);
  mg->Add(data_v24);


  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("v_{2}{n}");

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(0.00);
  mg->SetMaximum(0.18);

  TLegend* leg = new TLegend(0.12,0.6,0.62,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("ALICE experiment, LHC15o (full dataset), RFPs flow, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(graphv22Gap10,"v_{2}{2},|#Delta#eta| > 1.0} PUBLISHED","p");
  leg->AddEntry(data_v22Gap10,"v_{2}{2},|#Delta#eta| > 1.0}","p");
  leg->AddEntry(graph_v24,"v_{2}{4} PUBLISHED","p");
  leg->AddEntry(data_v24,"v_{2}{4}","p");
  leg->AddEntry(graph_v26,"v_{2}{6} PUBLISHED","p");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(graph_v28,"v_{2}{8} PUBLISHED","p");

  leg->Draw("same");

  can->SaveAs("v2n_comp.pdf");


}
