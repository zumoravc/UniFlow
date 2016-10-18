TCanvas* ExtractFlow(TH1D* hInvMass, TH1D* hFlowMass, Double_t* dConV2,const Short_t iPartSpecies = 0);

void V0sExtractFlow(
		const TString sInput = "~/NBI/Codes/results/V0s/8/plusplus_part1/plots/V0sFlow.root",
		const TString sOutput = "~/NBI/Codes/results/V0s/8/plusplus/plots",
		const TString sTag = "_JHEP",
		const TString sEtaGap = "Gap10",
		const TString sOutputFormat = "png",
		const Bool_t bLoadHEP = "kFALSE"
	)
{


	const Int_t iNumPtBins = 10; // pT bins
	const Int_t iNumCentBins = 9; // centrality bins
	// bins edges
	Double_t fPtBinEdges[] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; 
	Double_t fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

	// =======================================
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	
	fInput->cd();

	// ===== Loading input ===== 
	TList* lInputInvMassK0s = (TList*) gDirectory->Get("lInvMass_K0s");
	TList* lInputInvMassLambda = (TList*) gDirectory->Get("lInvMass_Lambda");
	TList* lInputFlowMassK0s = (TList*) gDirectory->Get("lFlowMass_K0s");
	TList* lInputFlowMassLambda = (TList*) gDirectory->Get("lFlowMass_Lambda");
	

	TH1D* hInvMassK0s[iNumCentBins][iNumPtBins];
	TH1D* hInvMassLambda[iNumCentBins][iNumPtBins];
	TH1D* hFlowMassK0s[iNumCentBins][iNumPtBins];
	TH1D* hFlowMassLambda[iNumCentBins][iNumPtBins];
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			hInvMassK0s[i][j] = (TH1D*) lInputInvMassK0s->FindObject(Form("hInvMass_K0s_Cent%d_pt%d",i,j))->Clone(Form("hInvMassK0s_Cent%d_pt%d",i,j));
			hInvMassLambda[i][j] = (TH1D*) lInputInvMassLambda->FindObject(Form("hInvMass_Lambda_Cent%d_pt%d",i,j))->Clone(Form("hInvMassLambda_Cent%d_pt%d",i,j));
			hFlowMassK0s[i][j] = (TH1D*) lInputFlowMassK0s->FindObject(Form("hFlowMass_K0s_Cent%d_pt%d",i,j))->Clone(Form("hFlowMassK0s_Cent%d_pt%d",i,j));
			hFlowMassLambda[i][j] = (TH1D*) lInputFlowMassLambda->FindObject(Form("hFlowMass_Lambda_Cent%d_pt%d",i,j))->Clone(Form("hFlowMassLambda_Cent%d_pt%d",i,j));
		}
	}

	// ===== Extracting signal & flow ===== 
	TCanvas* cCan;

	//cCan = ExtractFlow(hInvMassK0s[4][0],hFlowMassK0s[4][0]);
	//cCan->Print(Form("%s/fitK0s/test.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

	TCanvas* cPtFlow = new TCanvas("cPtFlow","PtFlow",1500,1500);
	TH1D* hFlowPt_K0s[iNumCentBins]; 
	TH1D* hFlowPt_Lambda[iNumCentBins]; 
	Double_t dExtractedV2_K0s[iNumPtBins] = {0};
	Double_t dExtractedV2_Lambda[iNumPtBins] = {0};
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		hFlowPt_K0s[i] = new TH1D(Form("hFlowPt_K0s_%s_Cent%d",sEtaGap.Data(),i),Form("K0s: v2 %s Cent %d; #it{p}_{T} (GeV/#it{c}); #it{v_2}",sEtaGap.Data(),i),iNumPtBins,fPtBinEdges);
		hFlowPt_Lambda[i] = new TH1D(Form("hFlowPt_Lambda_%s_Cent%d",sEtaGap.Data(),i),Form("#Lambda: v2 %s Cent %d; #it{p}_{T} (GeV/#it{c}); #it{v_2}",sEtaGap.Data(),i),iNumPtBins,fPtBinEdges);
		
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			cCan = ExtractFlow(hInvMassK0s[i][j],hFlowMassK0s[i][j],&dExtractedV2_K0s[j],0);
			cCan->Print(Form("%s/fitK0s/fit_K0s_Cent_%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			hFlowPt_K0s[i]->SetBinContent(j+1,dExtractedV2_K0s[j]);
			
			cCan = ExtractFlow(hInvMassLambda[i][j],hFlowMassLambda[i][j],&dExtractedV2_Lambda[j],1);
			cCan->Print(Form("%s/fitLambda/fit_Lambda_Cent_%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			hFlowPt_Lambda[i]->SetBinContent(j+1,dExtractedV2_Lambda[j]);
		}

		cPtFlow->cd();
		hFlowPt_K0s[i]->Draw();
		cPtFlow->Print(Form("%s/finalFlow/flowPt_K0s_%s_Cent_%d.%s",sOutput.Data(),sEtaGap.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());

		cPtFlow->cd();
		hFlowPt_Lambda[i]->Draw();
		cPtFlow->Print(Form("%s/finalFlow/flowPt_Lambda_%s_Cent_%d.%s",sOutput.Data(),sEtaGap.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
	}



	// loading JHEP HEP data

	if(bLoadHEP)
	{
		TFile* fJHEP = new TFile("/Users/vpacik/NBI/Flow/hepData/PID_v2_PbPb/HEPData-ins1297103-1-root.root","READ");

		const Short_t iNumCentBinsHEP = 7;
		TGraphAsymmErrors* graphJHEP_K0s[iNumCentBinsHEP];
		TGraphAsymmErrors* graphJHEP_Lambda[iNumCentBinsHEP];

		for(Int_t i(0); i < 7; i++)
		{
			fJHEP->cd(Form("Table %d",15+i));
			graphJHEP_K0s[i] = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1")->Clone(Form("graphJHEP_K0s_Cent_%d",i));	
			cPtFlow->cd();
			hFlowPt_K0s[i]->Draw();
			graphJHEP_K0s[i]->Draw("same");
			cPtFlow->Print(Form("%s/finalFlow/flowPt_K0s_withJHEP_%s_Cent_%d.%s",sOutput.Data(),sEtaGap.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());

			fJHEP->cd(Form("Table %d",41+i));
			graphJHEP_Lambda[i] = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1")->Clone(Form("graphJHEP_Lambda_Cent_%d",i));	
			cPtFlow->cd();
			hFlowPt_Lambda[i]->Draw();
			graphJHEP_Lambda[i]->Draw("same");
			cPtFlow->Print(Form("%s/finalFlow/flowPt_Lambda_withJHEP_%s_Cent_%d.%s",sOutput.Data(),sEtaGap.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
		} 

		
	}
}

TCanvas* ExtractFlow(TH1D* hInvMass, TH1D* hFlowMass, Double_t* dConV2, const Short_t iPartSpecies)
{
	const Short_t dNumSigma = 5;

	// ==========================================

	if(!hInvMass || !hFlowMass)
	{
		printf("Invalid pointer to input histos!\n");
		return 0x0;
	}
	
	Double_t dMassPeak = 0., dMassPeakLimit = 0.;
	Double_t dFitRange[2] = {0};
	TString sPartName;

	switch(iPartSpecies)
	{
		case 0: //K0s
			sPartName = TString("K0s");
			dMassPeak = 0.49;
			dMassPeakLimit = 0.03;
			dFitRange[0] = 0.4;
			dFitRange[1] = 0.6;

		break;

		case 1: //Lambda
			sPartName = TString("Lambda");
			dMassPeak = 1.116;
			dMassPeakLimit = 0.01;
			dFitRange[0] = 1.08;
			dFitRange[1] = 1.16;
		break;

		default:
			return;
		break;
	}

	printf("=== Fitting %s ===\n", sPartName.Data());

	TH1D* hInvMass_temp = (TH1D*) hInvMass->Clone("hInvMass_temp"); // cloning inv mass dist for peak window fitting
	TH1D* hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side"); // cloning inv mass dist for sidebands fitting
	TH1D* hInvMass_subs = (TH1D*) hInvMass->Clone("hInvMass_subs"); // cloning inv mass dist for BG subtracttion
	TH1D* hInvMassRebin = (TH1D*) hInvMass->Clone("hInvMassRebin"); // cloning inv mass dist for rebining

	TCanvas* cCan = new TCanvas("cCan","cFit",2000,1000);
	cCan->Divide(2,1);
	cCan->cd(1);
	hInvMass->Draw();
	cCan->cd(2);
	hFlowMass->Draw();

	// fitting inv mass dist
	TCanvas* cCanTemp = new TCanvas("cCanTemp","Temp",2000,2000);
	cCanTemp->Divide(3,2);
	cCanTemp->cd(1);

	TF1* fFitInvMass = new TF1("fFitInvMass","gaus(0)+pol3(3)",dFitRange[0],dFitRange[1]); 
	fFitInvMass->SetParNames("Amp","Mean","Sigma");
	fFitInvMass->SetNpx(1000000);
	fFitInvMass->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMass->SetParameter(1,dMassPeak);
	fFitInvMass->SetParameter(2,0.005);
	fFitInvMass->SetParLimits(2,0.,0.01);
	//fFitInvMass->SetLineColor(kGreen);
	fFitInvMass->SetLineWidth(2);
	//hInvMass->SetMinimum(0);
	//hInvMass->SetMaximum(4000);	

	hInvMass->Fit("fFitInvMass","R");

	// * TODO - fit verification procedure * //

	// extracting mean & sigma
	Double_t dMean = fFitInvMass->GetParameter(1);
	Double_t dSigma = fFitInvMass->GetParameter(2);
	
	Double_t dMassWindow[2] = {dMean - dNumSigma*dSigma, dMean + dNumSigma*dSigma}; // setting inv. mass peak
	printf("dMassWindow: %g - %g\n", dMassWindow[0],dMassWindow[1]);


	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hInvMass->GetBinCenter(i) > dMassWindow[0] && hInvMass->GetBinCenter(i) < dMassWindow[1])
		{
			hInvMass_side->SetBinError(i,9999999999999); 
		}
	}

	/*
	TF1* fFitInvMassPeak = new TF1("fFitInvMassPeak","gaus(0)+pol3(3)",dMassWindow[0],dMassWindow[1]);
	fFitInvMassPeak->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMassPeak->SetParameter(1,0.49);
	fFitInvMassPeak->SetParameter(2,0.005);
	fFitInvMassPeak->SetParLimits(2,0.,0.01);
	fFitInvMassPeak->SetLineWidth(2);
	fFitInvMassPeak->SetLineStyle(2);
	fFitInvMassPeak->SetLineColor(kPink+2);
	fFitInvMassPeak->SetNpx(1000000);
	*/


	TF1* fFitInvMassBg = new TF1("fFitInvMassBg","pol3",dFitRange[0],dFitRange[1]);
	fFitInvMassBg->SetLineWidth(2);
	fFitInvMassBg->SetLineStyle(2);
	fFitInvMassBg->SetLineColor(kGreen+2);
	
	cCanTemp->cd(2);
	//hInvMass_temp->SetMarkerColor(kRed);
	//hInvMass_temp->SetLineColor(kRed);
	hInvMass_side->SetMinimum(0);
	hInvMass_side->SetMaximum(100000);
	hInvMass_side->Fit("fFitInvMassBg","R");
	hInvMass_side->Draw("");
	/*
	cCanTemp->cd(3);
	hInvMass_temp->SetMinimum(0);
	hInvMass_temp->SetMaximum(100000);
	hInvMass_temp->Fit("fFitInvMassPeak","R");
	hInvMass_temp->Draw();

	TF1* fFitInvMassFixed = new TF1("fFitInvMassFixed","gaus(0)+pol3(3)",0.4,0.6);
	fFitInvMassFixed->SetParameter(1,0.49);
	fFitInvMassFixed->SetParameter(2,0.005);
	fFitInvMassFixed->SetParLimits(2,0.,0.01);
	fFitInvMassFixed->FixParameter(3, fFitInvMassBg->GetParameter(0));
	fFitInvMassFixed->FixParameter(4, fFitInvMassBg->GetParameter(1));	
	fFitInvMassFixed->FixParameter(5, fFitInvMassBg->GetParameter(2));	
	fFitInvMassFixed->FixParameter(6, fFitInvMassBg->GetParameter(3));	
	fFitInvMassFixed->SetLineWidth(2);
	fFitInvMassFixed->SetLineStyle(2);
	fFitInvMassFixed->SetLineColor(kYellow+2);
	fFitInvMassFixed->SetNpx(1000000);
	hInvMass_temp->Fit("fFitInvMassFixed","R");

	fFitInvMassGauss = new TF1("fFitInvMassGauss","gaus(0)",0.4,0.6);
	fFitInvMassGauss->SetLineWidth(2);
	fFitInvMassGauss->SetLineStyle(3);
	fFitInvMassGauss->SetLineColor(kRed);
	fFitInvMassGauss->SetParameters(fFitInvMassFixed->GetParameter(0),fFitInvMassFixed->GetParameter(1),fFitInvMassFixed->GetParameter(2));
	*/

	cCanTemp->cd(4);
	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		// subtracting the BG
		hInvMass_subs->SetBinContent(i,hInvMass->GetBinContent(i) - fFitInvMassBg->Eval(hInvMass_subs->GetBinCenter(i)));
	}

	TH1D* hInvMass_RatioSigTot = (TH1D*) hInvMass_subs->Clone("hInvMass_RatioSigTot");

	TF1* fFitInvMassSubs = new TF1("fFitInvMassSubs","gaus(0)+pol3(3)",dFitRange[0],dFitRange[1]); 
	fFitInvMassSubs->SetParameter(1,dMassPeak);
	fFitInvMassSubs->SetParameter(2,0.005);
	fFitInvMassSubs->SetParLimits(2,0.,0.02);

	fFitInvMassSubs->SetLineWidth(2);
	fFitInvMassSubs->SetLineStyle(2);
	fFitInvMassSubs->SetLineColor(kRed+2);
	fFitInvMassSubs->SetNpx(1000000);
	hInvMass_subs->Fit("fFitInvMassSubs","R");
	hInvMass_subs->SetMaximum(2000);
	hInvMass_subs->Draw();

	cCanTemp->cd(1);
	fFitInvMassSubs->Draw("same");

	cCanTemp->cd(5);
	hInvMass_RatioSigTot->Divide(hInvMass);
	TF1* fFitInvMassRatioSigTot = new TF1("fFitInvMassRatioSigTot","gaus(0)+gaus(3)+pol2(6)",dFitRange[0],dFitRange[1]);
	fFitInvMassRatioSigTot->SetParameter(0,0.7);
	fFitInvMassRatioSigTot->SetParameter(1,dMassPeak);
	fFitInvMassRatioSigTot->SetParLimits(1,dMassPeak-dMassPeakLimit,dMassPeak+dMassPeakLimit);
	fFitInvMassRatioSigTot->SetParameter(2,0.05);
	fFitInvMassRatioSigTot->SetParLimits(2,0.,0.03);
	fFitInvMassRatioSigTot->SetParameter(3,0.7);
	fFitInvMassRatioSigTot->SetParameter(4,dMassPeak);
	fFitInvMassRatioSigTot->SetParLimits(4,dMassPeak-dMassPeakLimit,dMassPeak+dMassPeakLimit);
	fFitInvMassRatioSigTot->SetParameter(5,0.05);
	fFitInvMassRatioSigTot->SetParLimits(5,0.,0.03);
	//fFitInvMassRatioSigTot->FixParameter(6,fFitFlowMass_side->GetParameter(0));
	//fFitInvMassRatioSigTot->FixParameter(7,fFitFlowMass_side->GetParameter(1));
	//fFitInvMassRatioSigTot->FixParameter(8,fFitFlowMass_side->GetParameter(2));
	hInvMass_RatioSigTot->Fit("fFitInvMassRatioSigTot","R");
	hInvMass_RatioSigTot->Draw();


	//cCanTemp->Print("~/NBI/Flow/results/V0s/10/plots_JHEP/fitK0s/temp.png","png");

	// Drawing over the original plot

	cCan->cd(1);
	fFitInvMassBg->Draw("same");
	//fFitInvMassPeak->Draw("same");
	//fFitInvMassFixed->Draw("same");
	//fFitInvMassGauss->Draw("same");

	/*

	// estimating number of candidates
	const Int_t iNumMethod = 7; // 0 from approx fit / 1 signal from peak fit / 2 signal bin counting 
	TString sMethod[] = {"AproxFit","PeakFit","BinCounting","AproxFitIntegral","PeakFitIntegral","FixedFitIntegral","GaussFitIntegral"};
	Double_t dCandTot[iNumMethod] = {0};
	Double_t dCandSig[iNumMethod] = {0};
	Double_t dCandBg[iNumMethod] = {0};
	
	Short_t iIndexMethod = 0;
	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		if(hInvMass->GetBinCenter(i) > dMassWindow[0] && hInvMass->GetBinCenter(i) < dMassWindow[1])
		{
			iIndexMethod = 0;
			dCandTot[iIndexMethod] += hInvMass->GetBinContent(i);
			dCandBg[iIndexMethod] += fFitInvMassBg->Eval(hInvMass->GetBinCenter(i));
			dCandSig[iIndexMethod] += fFitInvMass->Eval(hInvMass->GetBinCenter(i)) - fFitInvMassBg->Eval(hInvMass->GetBinCenter(i) );

			iIndexMethod++;

			dCandTot[iIndexMethod] += hInvMass->GetBinContent(i);
			dCandBg[iIndexMethod] += fFitInvMassBg->Eval(hInvMass->GetBinCenter(i));
			dCandSig[iIndexMethod] += fFitInvMassPeak->Eval(hInvMass->GetBinCenter(i)) - fFitInvMassBg->Eval(hInvMass->GetBinCenter(i));

			iIndexMethod++;

			dCandTot[iIndexMethod] += hInvMass->GetBinContent(i);
			dCandBg[iIndexMethod] += fFitInvMassBg->Eval(hInvMass->GetBinCenter(i));
			dCandSig[iIndexMethod] += hInvMass->GetBinContent(i) - fFitInvMassBg->Eval(hInvMass->GetBinCenter(i)); 
	
			iIndexMethod++;
		}
	}

	Int_t nPoints = 1000;
	Double_t* dX = new Double_t[nPoints];
	Double_t* dW = new Double_t[nPoints];
	fFitInvMass->CalcGaussLegendreSamplingPoints(nPoints,dX,dW,1e-15);
	fFitInvMassBg->CalcGaussLegendreSamplingPoints(nPoints,dX,dW,1e-15);
	fFitInvMassFixed->CalcGaussLegendreSamplingPoints(nPoints,dX,dW,1e-15);
	fFitInvMassGauss->CalcGaussLegendreSamplingPoints(nPoints,dX,dW,1e-15);

	Double_t dIntTot = fFitInvMass->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) / hInvMass->GetBinWidth(5);
	Double_t dIntBg = fFitInvMassBg->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) /  hInvMass->GetBinWidth(5);
	Double_t dIntSig = (dIntTot - dIntBg);

	dCandTot[iIndexMethod] = dIntTot;
	dCandBg[iIndexMethod] = dIntBg;
	dCandSig[iIndexMethod] = dIntSig;

	iIndexMethod++;

	dIntTot = fFitInvMassPeak->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) / hInvMass->GetBinWidth(5);
 	dIntBg = fFitInvMassBg->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) /  hInvMass->GetBinWidth(5);
	dIntSig = (dIntTot - dIntBg);

	dCandTot[iIndexMethod] = dIntTot;
	dCandBg[iIndexMethod] = dIntBg;
	dCandSig[iIndexMethod] = dIntSig;

	iIndexMethod++;

	dIntTot = fFitInvMassFixed->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) / hInvMass->GetBinWidth(5);
 	dIntBg = fFitInvMassBg->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) /  hInvMass->GetBinWidth(5);
	dIntSig = (dIntTot - dIntBg);

	dCandTot[iIndexMethod] = dIntTot;
	dCandBg[iIndexMethod] = dIntBg;
	dCandSig[iIndexMethod] = dIntSig;

	iIndexMethod++;

	dIntTot = fFitInvMassGauss->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]) / hInvMass->GetBinWidth(5);
 	dIntBg = 0; // bg subtracted
	dIntSig = dIntTot;

	dCandTot[iIndexMethod] = dIntTot;
	dCandBg[iIndexMethod] = dIntBg;
	dCandSig[iIndexMethod] = dIntSig;

	iIndexMethod++;


	for(Int_t i(0); i < iNumMethod; i++)
	{
		printf("Cand::Method %s: Sig %g / Bg %g / Total %g (Sig+BG %g) / Sig/Bg %g\n",sMethod[i].Data(),dCandSig[i],dCandBg[i], dCandTot[i], dCandSig[i]+dCandBg[i],dCandSig[i]/dCandBg[i]);
	}
	
	*/

	// fitting flow mass dist
	TCanvas* cCanFlowTemp = new TCanvas("cCanFlowTemp","FlowTemp",2000,2000);
	cCanFlowTemp->Divide(2,2);
	cCanFlowTemp->cd(1);
	
	TH1D* hFlowMass_temp = (TH1D*) hFlowMass->Clone("hFlowMass_temp");
	TH1D* hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");

	// side bands fitting

	for(Int_t i(1); i < hFlowMass->GetNbinsX()+1; i++)
	{
		if(hFlowMass->GetBinCenter(i) > dMassWindow[0] && hFlowMass->GetBinCenter(i) < dMassWindow[1])
		{
			hFlowMass_side->SetBinError(i,9999999999999); 
		}
	}
	
	TF1* fFitFlowMass_side = new TF1("fFitFlowMass_side","pol2(0)",dFitRange[0],dFitRange[1]);
	hFlowMass_side->Fit("fFitFlowMass_side","R");

	// Fitting by vTot (LH side)

	TF1* fFitFlowMassTot = new TF1("fFitFlowMassTot","[12]*(gaus(0)+gaus(3)+pol2(6))+(1-(gaus(0)+gaus(3)+pol2(6)))*pol2(9)",dFitRange[0],dFitRange[1]);
	// inv mass fit parameters
	for(Int_t i(0); i < 9; i++)
	{
		fFitFlowMassTot->FixParameter(i,fFitInvMassRatioSigTot->GetParameter(i));
	}
	// vn bg fit parameters
	fFitFlowMassTot->FixParameter(9,fFitFlowMass_side->GetParameter(0));
	fFitFlowMassTot->FixParameter(10,fFitFlowMass_side->GetParameter(1));
	fFitFlowMassTot->FixParameter(11,fFitFlowMass_side->GetParameter(2));
	hFlowMass_temp->Fit("fFitFlowMassTot","R");

	*dConV2 = fFitFlowMassTot->GetParameter(12);

	// drawing

	cCanFlowTemp->cd(1);
	hFlowMass_temp->SetMinimum(0.00);
	hFlowMass_temp->SetMaximum(0.11);
	hFlowMass_temp->Draw();

	cCanFlowTemp->cd(2);
	hFlowMass_side->SetMinimum(0.00);
	hFlowMass_side->SetMaximum(0.11);
	hFlowMass_side->Draw();

	cCan->cd(2);
	fFitFlowMassTot->Draw("same");

	return cCan;	
}
