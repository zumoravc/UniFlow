TCanvas* PlotCent(const Short_t iCent, TString sTag);
TGraphErrors* Kv2_1020_QC2(Int_t color=1, Int_t marker=20);

TGraphAsymmErrors* v2Pion1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1);
TGraphAsymmErrors* v2Kaon1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1);
TGraphAsymmErrors* v2Antiproton1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1);
TGraphErrors* Kv2_1020_QC4(Int_t color=1, Int_t marker=20);

void PlotPID()
{
	TFile* fInput = new TFile("~/NBI/Flow/classProcess/ProcessedFlow.root","READ");
	fInput->cd();
	//fInput->ls();

	//TGraphErrors* Kv2_1020_QC2 = Kv2_1020_QC2();
	TGraphAsymmErrors* v2Pion1020 = v2Pion1020();	
	TGraphAsymmErrors* v2Kaon1020 = v2Kaon1020();	
	TGraphAsymmErrors* v2aProton1020 = v2Antiproton1020();
  TGraphErrors* Kv2_1020_QC4 = Kv2_1020_QC4();

	TCanvas* Diff22_cent3 = PlotCent(2,"Diffv22");
	Diff22_cent3->cd();
	v2Pion1020->Draw("same");
	v2Kaon1020->Draw("same");
	v2aProton1020->Draw("same");

	//Kv2_1020_QC2->Draw("same");
	TCanvas* Diff24_cent3 = PlotCent(2,"Diffv24");
  Diff24_cent3->cd();
  Kv2_1020_QC4->Draw("same");
/*
	const Short_t iCent = 3;

	TH1D* hPion22 = (TH1D*) gDirectory->Get(Form("hPionsDiffv22_cent%d",iCent));
	TH1D* hKaon22 = (TH1D*) gDirectory->Get(Form("hKaonsDiffv22_cent%d",iCent));
	TH1D* hProton22 = (TH1D*) gDirectory->Get(Form("hProtonsDiffv22_cent%d",iCent));
	

	Color_t colors[3] = {kRed,kBlue,kGreen+2};
	TCanvas* cCent = new TCanvas();
	cCent->cd();
	hPion22->SetLineColor(colors[0]);
	hPion22->SetMarkerColor(colors[0]);
	hKaon22->SetLineColor(colors[1]);
	hKaon22->SetMarkerColor(colors[1]);
	hProton22->SetLineColor(colors[2]);
	hProton22->SetMarkerColor(colors[2]);
	
	hPion22->Draw();
	hKaon22->Draw("same");
	hProton22->Draw("same");
*/


	//PlotCent(3);

	return;
}

TCanvas* PlotCent(const Short_t iCent, TString sTag)
{
	TH1D* hPion22 = (TH1D*) gDirectory->Get(Form("hPions%s_cent%d",sTag.Data(),iCent));
	TH1D* hKaon22 = (TH1D*) gDirectory->Get(Form("hKaons%s_cent%d",sTag.Data(),iCent));
	TH1D* hProton22 = (TH1D*) gDirectory->Get(Form("hProtons%s_cent%d",sTag.Data(),iCent));
	

	Color_t colors[3] = {kRed,kBlue,kGreen+2};
	TCanvas* cCent = new TCanvas();
	cCent->cd();
	hPion22->SetMinimum(-0.05);
	hPion22->SetMaximum(0.3);
	hPion22->SetLineColor(colors[0]);
	hPion22->SetMarkerColor(colors[0]);
	hKaon22->SetLineColor(colors[1]);
	hKaon22->SetMarkerColor(colors[1]);
	hProton22->SetLineColor(colors[2]);
	hProton22->SetMarkerColor(colors[2]);
	
	hPion22->Draw();
	hKaon22->Draw("same");
	hProton22->Draw("same");

	return cCent;
}

TGraphErrors* Kv2_1020_QC2(Int_t color, Int_t marker) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.007051, 0.016100, 0.037009, 0.057360, 0.076636, 0.096443, 0.110868, 0.124309, 0.134059, 0.143032, 0.148407, 0.156759, 0.161398, 0.156367, 0.148881, 0.142463, 0.135274, 0.108699, 0.109107};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.003485, 0.001355, 0.000849, 0.000704, 0.000690, 0.000741, 0.000824, 0.000949, 0.001078, 0.001346, 0.001635, 0.001515, 0.002217, 0.003181, 0.004494, 0.004505, 0.008140, 0.011124, 0.020627};
  Double_t _ysys[] = {0.001480, 0.000588, 0.001159, 0.001734, 0.002305, 0.002927, 0.003347, 0.003831, 0.004036, 0.004330, 0.004453, 0.004743, 0.005010, 0.005319, 0.004618, 0.004824, 0.004091, 0.003268, 0.013869};
  //if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  //if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  //if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Pion1020(Int_t color, Int_t marker, Int_t first,Int_t last)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0.0201224, 0.0252062, 0.0311323, 0.0369631, 0.0428155, 0.0483589, 0.0550696, 0.0651319, 0.0746561, 0.0834061, 0.0919292, 0.0993882, 0.106362, 0.11341, 0.119441, 0.125327, 0.130489, 0.135768, 0.140111, 0.1436, 0.147624, 0.15197, 0.157217, 0.15934, 0.159981, 0.157703, 0.15448, 0.151198, 0.142915, 0.132599, 0.113408, 0.10271, 0.0943079, 0.0879974, 0.0790374, 0.0795734, 0.0614621, 0.067106, 0.0609654, 0.0515049, 0.0223675};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 7.16967e-05, 7.04768e-05, 7.14158e-05, 7.34838e-05, 7.64088e-05, 8.00149e-05, 8.9464e-05, 9.55187e-05, 0.000105152, 0.000116857, 0.000130562, 0.000146276, 0.000163809, 0.000183749, 0.00020652, 0.000232603, 0.000262503, 0.000296109, 0.000334244, 0.000377637, 0.00042637, 0.00033071, 0.000441212, 0.000580034, 0.000749139, 0.000958424, 0.00120529, 0.00149142, 0.00182397, 0.0016818, 0.00233162, 0.00308243, 0.00393626, 0.00487286, 0.0058437, 0.00539979, 0.00725315, 0.00952612, 0.0093786, 0.0123107, 0.0176851};
  Double_t _yerr2[] = {0, 0, 0, 0, 7.16967e-05, 7.04768e-05, 7.14158e-05, 7.34838e-05, 7.64088e-05, 8.00149e-05, 8.9464e-05, 9.55187e-05, 0.000105152, 0.000116857, 0.000130562, 0.000146276, 0.000163809, 0.000183749, 0.00020652, 0.000232603, 0.000262503, 0.000296109, 0.000334244, 0.000377637, 0.00042637, 0.00033071, 0.000441212, 0.000580034, 0.000749139, 0.000958424, 0.00120529, 0.00149142, 0.00182397, 0.0016818, 0.00233162, 0.00308243, 0.00393626, 0.00487286, 0.0058437, 0.00539979, 0.00725315, 0.00952612, 0.0093786, 0.0123107, 0.0176851};

  /*
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.36284000000000027e-04;
	Float_t pol1 = 8.43209000000000037e-04;
	
	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPi(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }
  */

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

TGraphAsymmErrors* v2Kaon1020(Int_t color, Int_t marker, Int_t first,Int_t last)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0.00640075, 0.00890777, 0.0114231, 0.0121152, 0.0193763, 0.0262718, 0.0373136, 0.04773, 0.0585364, 0.0685442, 0.0781197, 0.0884599, 0.0959829, 0.104197, 0.110977, 0.118551, 0.124407, 0.130268, 0.137113, 0.141952, 0.149588, 0.158874, 0.165006, 0.168625, 0.168325, 0.16697, 0.173543, 0.161645, 0.15523, 0.14532, 0.128752, 0.142204, 0.0843617, 0.269387, 0.0231167, 0.431226, 0.0167267, 0.315711, 0.235569, -0.000340595};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0.000525284, 0.000429881, 0.000408989, 0.000491132, 0.000713549, 0.000405664, 0.000346767, 0.000322585, 0.000314851, 0.000318631, 0.00032842, 0.000344056, 0.000364314, 0.00039006, 0.000423204, 0.000459161, 0.000504151, 0.000555309, 0.000614542, 0.000683822, 0.000519796, 0.000688975, 0.00091857, 0.00122243, 0.00160751, 0.00211499, 0.00275193, 0.0035775, 0.00370145, 0.00631253, 0.0106874, 0.0182738, 0.0313279, 0.0476483, 0.0611653, 0.106198, 0.207285, 0.182173, 0.191574, 0.224581};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.000525284, 0.000429881, 0.000408989, 0.000491132, 0.000713549, 0.000405664, 0.000346767, 0.000322585, 0.000314851, 0.000318631, 0.00032842, 0.000344056, 0.000364314, 0.00039006, 0.000423204, 0.000459161, 0.000504151, 0.000555309, 0.000614542, 0.000683822, 0.000519796, 0.000688975, 0.00091857, 0.00122243, 0.00160751, 0.00211499, 0.00275193, 0.0035775, 0.00370145, 0.00631253, 0.0106874, 0.0182738, 0.0313279, 0.0476483, 0.0611653, 0.106198, 0.207285, 0.182173, 0.191574, 0.224581};

  /*
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 7.36284000000000027e-04;
	Float_t pol1 = 8.43209000000000037e-04;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systKa(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }
	*/

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}


TGraphAsymmErrors* v2Antiproton1020(Int_t color, Int_t marker, Int_t first,Int_t last)
{
  //commentme
  Int_t _nPoints = 45;
  if (last>_nPoints-1) last=_nPoints-1;
  if (last<0 && first<0) last=_nPoints-1;
  if (last<0) last=_nPoints-1+last;
  if (first<0) first=0;
  Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
  Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00514464, 0.00344225, 0.00468868, 0.00546, 0.0082079, 0.012222, 0.0165934, 0.0233094, 0.0288221, 0.0374016, 0.0454486, 0.0563712, 0.064759, 0.075207, 0.0867709, 0.0963642, 0.106439, 0.116716, 0.12656, 0.138558, 0.157754, 0.173931, 0.185026, 0.199048, 0.206145, 0.208087, 0.213577, 0.213699, 0.21574, 0.192142, 0.152677, 0.172787, 0.122158, 0.11594, 0.170479, 0.143031, 0.0113144, 0.188084, 0.909081};
  Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00163153, 0.00131053, 0.00108878, 0.000950215, 0.000582294, 0.000516432, 0.000703097, 0.000643592, 0.000611608, 0.000594501, 0.000588832, 0.000590229, 0.000600458, 0.000619565, 0.00064412, 0.00067877, 0.000720346, 0.000770869, 0.000828332, 0.000600635, 0.000741676, 0.000928701, 0.00116982, 0.00102807, 0.00131593, 0.00169927, 0.00218807, 0.00222766, 0.00363962, 0.00585315, 0.00929629, 0.013941, 0.020436, 0.0231194, 0.046922, 0.0895645, 0.111115, 0.118019, 0.100696};
  Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00163153, 0.00131053, 0.00108878, 0.000950215, 0.000582294, 0.000516432, 0.000703097, 0.000643592, 0.000611608, 0.000594501, 0.000588832, 0.000590229, 0.000600458, 0.000619565, 0.00064412, 0.00067877, 0.000720346, 0.000770869, 0.000828332, 0.000600635, 0.000741676, 0.000928701, 0.00116982, 0.00102807, 0.00131593, 0.00169927, 0.00218807, 0.00222766, 0.00363962, 0.00585315, 0.00929629, 0.013941, 0.020436, 0.0231194, 0.046922, 0.0895645, 0.111115, 0.118019, 0.100696};

  /*
  if(!kStat){
    for(Int_t i=0;i<_nPoints;i++){
      _yerr[i] = 0;
      _yerr2[i] = 0;
      _xerr[i] = 0.05;
    }
  }


  if(kSyst){
    for(Int_t i=0;i<_nPoints;i++){
	Float_t pol0 = 2.13364693747489361e-03;
	Float_t pol1 = 1.52812297971200361e-03;

	Float_t nonflow = (pol0 + pol1*_x[i]);
	Float_t systerr = _y[i]*systPr(2, _x[i]);
	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
    }
  }
	*/
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}
TGraphErrors* Kv2_1020_QC4(Int_t color, Int_t marker) {
  Int_t _nPoints = 19;
  Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
  Double_t _y[] = {0.010304, 0.016204, 0.033282, 0.050568, 0.068130, 0.082275, 0.098983, 0.107410, 0.118314, 0.122915, 0.127885, 0.135673, 0.133027, 0.125030, 0.119996, 0.114883, 0.103795, 0.056606, 0.073532};
  Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
  Double_t _yerr[] = {0.008246, 0.003348, 0.002248, 0.001928, 0.001887, 0.001998, 0.002171, 0.002489, 0.002762, 0.003426, 0.004034, 0.003694, 0.005368, 0.007697, 0.010785, 0.010791, 0.019094, 0.026001, 0.046208};
  Double_t _ysys[] = {0.003572, 0.000550, 0.001082, 0.001520, 0.002068, 0.003303, 0.003411, 0.003818, 0.003650, 0.003819, 0.003868, 0.004158, 0.005506, 0.003772, 0.003747, 0.003794, 0.003375, 0.008648, 0.025117};
  //if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
  //if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
  //if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
  TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  return graph;
}

