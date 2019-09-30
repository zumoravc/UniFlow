// plotting macro
// v_2, v_3, and v_4 from two-particle correlations
// input:
// published data as numbers
// data produced in analysis as histograms
//
// ZM
// 29/08/2019

void MCvsPublished()
{
  TMultiGraph *mg = new TMultiGraph();


  Double_t xCross[10] = {2.5,7.5,15,25,35,45,55,65,75};
  Double_t xCrossErr[10] = {2.5,2.5,5,5,5,5,5,5,5};
  // RUN2 data
  Double_t v22Gap10Run2[10]={0.0283859, 0.0456604, 0.0655068, 0.0870721, 0.099105, 0.104143, 0.10286, 0.0974591, 0.0888104};
  Double_t v22Gap10Run2Err[10]={0.000570139, 0.000643862, 0.000373923, 0.000444072, 0.000553776, 0.000732208, 0.00107352, 0.00186236, 0.00437842};
  Double_t v22Gap10Run2Sys[10]={0.000425788, 0.000684906, 0.000982603, 0.00130608, 0.00148658, 0.00156214, 0.00154291, 0.00146189, 0.00133216};
  Double_t v22Gap10Run2CombErr[10]={0.000711585, 0.000940029, 0.00105135, 0.00137951, 0.00158637, 0.00172523, 0.00187963, 0.00236759, 0.0045766};

  Double_t v32Gap10Run2[10]={0.0206723, 0.0231991, 0.0279915, 0.0309315, 0.0331377, 0.0323443, 0.0270494, 0.0276913};
  Double_t v32Gap10Run2Err[10]={0.000629474, 0.000797443, 0.000466369, 0.000619115, 0.000842254, 0.00136072, 0.00283039, 0.00546999};
  Double_t v32Gap10Run2Sys[10]={0.000950925, 0.00106716, 0.00128761, 0.00142285, 0.00152434, 0.00148784, 0.00124427, 0.0012738};
  Double_t v32Gap10Run2CombErr[10]={0.00114039, 0.0013322, 0.00136947, 0.00155171, 0.00174155, 0.00201624, 0.00309182, 0.00561635};

  Double_t v42Gap10Run2[10]={0.0114652, 0.0129068, 0.0138326, 0.0156838, 0.017018, 0.0162893, 0.0177113, 0.0198688, 0.00916289};
  Double_t v42Gap10Run2Err[10]={0.000949648, 0.00120259, 0.000788619, 0.00104461, 0.00147443, 0.00245469, 0.00406993, 0.0073516, 0.0383903};
  Double_t v42Gap10Run2Sys[10]={0.000550329, 0.000619527, 0.000663966, 0.000752823, 0.000816866, 0.000781885, 0.000850141, 0.000953704, 0.000439819};
  Double_t v42Gap10Run2CombErr[10]={0.00109759, 0.00135279, 0.00103091, 0.00128761, 0.00168559, 0.00257621, 0.00415777, 0.0074132, 0.0383928};

  Double_t v24Run2[10]={-1, 0.0366982, 0.0586746, 0.0779774, 0.0867288, 0.0890264, 0.0835125, 0.0806163, 0.0464366};
  Double_t v24Run2Err[10]={0, 0.00258022, 0.000331448, 0.00138072, 0.00257394, 0.00275397, 0.00207766, 0.00961209, 0.0835755};
  Double_t v24Run2Sys[10]={0, 0.000623869, 0.000997468, 0.00132562, 0.00147439, 0.00151345, 0.00141971, 0.00137048, 0.000789422};
  Double_t v24Run2CombErr[10]={0, 0.00265457, 0.00105109, 0.00191407, 0.00296631, 0.00314243, 0.0025164, 0.0097093, 0.0835792};


  TGraphErrors *graphv22Gap10 = new TGraphErrors(9,xCross,v22Gap10Run2,xCrossErr,v22Gap10Run2CombErr);
  graphv22Gap10->SetMarkerStyle(kFullSquare);
  graphv22Gap10->SetMarkerColor(kRed+1);
  graphv22Gap10->SetMarkerSize(1.);
  graphv22Gap10->SetLineColor(kRed+1);
  mg->Add(graphv22Gap10);
  //ShiftAlongXaxis(graphv22Gap10,1.0);

  TGraphErrors *graphv32Gap10 = new TGraphErrors(8,xCross,v32Gap10Run2,xCrossErr,v32Gap10Run2CombErr);
  graphv32Gap10->SetMarkerStyle(kFullCircle);
  graphv32Gap10->SetMarkerColor(kBlue+1);
  graphv32Gap10->SetMarkerSize(1.1);
  graphv32Gap10->SetMarkerColor(kBlue+1);
  mg->Add(graphv32Gap10);
  //ShiftAlongXaxis(graphv32Gap10,-1.0);

  TGraphErrors *graphv42Gap10 = new TGraphErrors(8,xCross,v42Gap10Run2,xCrossErr,v42Gap10Run2CombErr);
  graphv42Gap10->SetMarkerStyle(kFullDiamond);
  graphv42Gap10->SetMarkerColor(kGreen+3);
  graphv42Gap10->SetMarkerSize(1.6);
  graphv42Gap10->SetLineColor(kGreen+3);
  mg->Add(graphv42Gap10);
  //ShiftAlongXaxis(graphv42Gap10,1.5);

  TFile* fileIn = TFile::Open("/home/alidock/ana/output/MC/AMPT_train2380/processUniFlow/Processed.root","READ");
  if(!fileIn) {printf("File not opened! \n"); return;}

  TFile* fileInH = TFile::Open("/home/alidock/ana/output/MC/HIJING_train_2381/processUniFlow/Processed.root","READ");
  if(!fileInH) {printf("File not opened! \n"); return;}

  TGraphErrors *data_v22Gap10 = new TGraphErrors((TH1D*) fileIn->Get("Refs_hFlow2_harm2_gap10"));
  data_v22Gap10->SetMarkerStyle(kOpenSquare);
  data_v22Gap10->SetMarkerColor(kRed);
  data_v22Gap10->SetMarkerSize(1.);
  data_v22Gap10->SetLineColor(kRed);
  mg->Add(data_v22Gap10);

  TGraphErrors *data_v32Gap10 = new TGraphErrors((TH1D*) fileIn->Get("Refs_hFlow2_harm3_gap10"));
  data_v32Gap10->SetMarkerStyle(kOpenCircle);
  data_v32Gap10->SetMarkerColor(kBlue);
  data_v32Gap10->SetMarkerSize(1.1);
  data_v32Gap10->SetMarkerColor(kBlue);
  mg->Add(data_v32Gap10);

  TGraphErrors *data_v42Gap10 = new TGraphErrors((TH1D*) fileIn->Get("Refs_hFlow2_harm4_gap10"));
  data_v42Gap10->SetMarkerStyle(kOpenDiamond);
  data_v42Gap10->SetMarkerColor(kGreen+2);
  data_v42Gap10->SetMarkerSize(1.6);
  data_v42Gap10->SetLineColor(kGreen+2);
  mg->Add(data_v42Gap10);

  TGraphErrors *HIJING_v22Gap10 = new TGraphErrors((TH1D*) fileInH->Get("Refs_hFlow2_harm2_gap10"));
  HIJING_v22Gap10->SetMarkerStyle(kOpenStar);
  HIJING_v22Gap10->SetMarkerColor(kBlack);
  HIJING_v22Gap10->SetMarkerSize(1.);
  HIJING_v22Gap10->SetLineColor(kBlack);
  mg->Add(HIJING_v22Gap10);

  mg->GetXaxis()->SetTitle("Centrality class (V0M)");
  mg->GetYaxis()->SetTitle("v_{n}{2}");

  TCanvas* can = new TCanvas("can", "can", 600, 400);
  mg->Draw("ap");
  mg->SetMinimum(-0.02);
  mg->SetMaximum(0.18);

  TLegend* leg = new TLegend(0.12,0.65,0.62,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.0,0.0);
  leg->SetHeader("RFPs flow, 0.2 < p_{T}(RFPs) < 5 GeV/c, |#eta| < 0.8");
  leg->SetNColumns(2);
  leg->AddEntry(graphv22Gap10,"v_{2}{2},|#Delta#eta| > 1.0} PUBLISHED","p");
  leg->AddEntry(data_v22Gap10,"v_{2}{2},|#Delta#eta| > 1.0} AMPT","p");
  leg->AddEntry(graphv32Gap10,"v_{3}{2},|#Delta#eta| > 1.0} PUBLISHED","p");
  leg->AddEntry(data_v32Gap10,"v_{3}{2},|#Delta#eta| > 1.0} AMPT","p");
  leg->AddEntry(graphv42Gap10,"v_{4}{2},|#Delta#eta| > 1.0} PUBLISHED","p");
  leg->AddEntry(data_v42Gap10,"v_{4}{2},|#Delta#eta| > 1.0} AMPT","p");
  leg->AddEntry(HIJING_v22Gap10,"v_{2}{2},|#Delta#eta| > 1.0} HIJING","p");
  leg->Draw("same");

  can->SaveAs("vn2_models.pdf");

  // TCanvas* can2 = new TCanvas("can2", "can2", 600, 400);
  // gStyle->SetOptStat(kFALSE);
  // TH1D *data_histov22 = (TH1D*) fileIn->Get("Refs_hFlow4_harm2_gap00");
  // TH1D *published_v22 = (TH1D*)data_histov22->Clone("published_v22");
  // for(int i = 1; i < 10; i++)
  // {
  //   published_v22->SetBinContent(i,v24Run2[i-1]);
  //   published_v22->SetBinError(i,v24Run2CombErr[i-1]);
  // }
  // published_v22->Divide(data_histov22);
  // published_v22->SetMarkerStyle(kFullSquare);
  // published_v22->SetMarkerColor(kBlack);
  // published_v22->SetLineColor(kBlack);
  // published_v22->SetTitle(";Centrality; v_{2}{4} published / v_{2}{4} reconstructed");
  // published_v22->Draw("h");
  // TF1 *fu1 = new TF1("fu1", "1", 0, 100);
  // fu1->SetLineColor(kRed);
  // fu1->Draw("same");
  // can2->SaveAs("ratio_v24.pdf");


}
