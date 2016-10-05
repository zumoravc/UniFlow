TCanvas* ExtractFlow(TH1D* hInvMass, TH1D* hFlowMass);

void V0sExtractFlow(
		const TString sInput = "~/NBI/Codes/results/V0s/8/plusplus_part1/plots/V0sFlow.root",
		const TString sOutput = "~/NBI/Codes/results/V0s/8/plusplus/plots",
		const TString sTag = "_JHEP",
		const TString sEtaGap = "Gap09",
		const TString sOutputFormat = "png"
	)
{


	const Int_t iNumPtBins = 10; // pT bins
	const Int_t iNumCentBins = 9; // centrality bins
	// bins edges
	Double_t fPtBinEdges[] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; 
	Double_t fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

	// =======================================
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareHistos.C");

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

	cCan = ExtractFlow(hInvMassK0s[0][2],hFlowMassK0s[0][2]);
	cCan->Print(Form("%s/fitK0s/test.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

	/*
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			cCan = ExtractFlow(hInvMassK0s[i][j],hFlowMassK0s[i][j]);
			cCan->Print(Form("%s/fitK0s/fit_K0s_Cent_%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
		}
	}
	*/

}

TCanvas* ExtractFlow(TH1D* hInvMass, TH1D* hFlowMass)
{
	const Short_t dNumSigma = 5.;

	// ==========================================

	if(!hInvMass || !hFlowMass)
	{
		printf("Invalid pointer to input histos!\n");
		return 0x0;
	}
	
	TH1D* hInvMass_temp = (TH1D*) hInvMass->Clone("hInvMass_temp"); // cloning inv mass dist for peak window fitting
	TH1D* hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side"); // cloning inv mass dist for sidebands fitting

	TCanvas* cCan = new TCanvas("cCan","temp",2000,2000);
	cCan->Divide(2,2);
	cCan->cd(1);
	hInvMass->Draw();
	cCan->cd(2);
	hFlowMass->Draw();

	// ===== Fitting K0s =====
	// fitting inv mass dist
	cCan->cd(1);
	TF1* fFitInvMass = new TF1("fFitInvMass","gaus(0)+pol3(3)",0.4,0.6); 
	fFitInvMass->SetParNames("Amp","Mean","Sigma");
	fFitInvMass->SetNpx(1000000);
	fFitInvMass->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMass->SetParameter(1,0.49);
	fFitInvMass->SetParameter(2,0.005);
	fFitInvMass->SetParLimits(2,0.,0.01);
	//fFitInvMass->SetLineColor(kGreen);
	fFitInvMass->SetLineWidth(2);
	hInvMass->SetMinimum(0);
	hInvMass->SetMaximum(4000);
	hInvMass->Fit("fFitInvMass","R");

	// * TODO - fit verification procedure * //

	// extracting mean & sigma
	Double_t dMean = fFitInvMass->GetParameter(1);
	Double_t dSigma = fFitInvMass->GetParameter(2);
	
	Double_t dMassWindow[2] = {dMean - dNumSigma*dSigma, dMean + dNumSigma*dSigma}; // setting inv. mass peak
	printf("dMassWindow: %g - %g\n", dMassWindow[0],dMassWindow[1]);


	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		// Excluding sidebands window (setting errors to inf)
		if(hInvMass->GetBinCenter(i) < dMassWindow[0] || hInvMass->GetBinCenter(i) > dMassWindow[1])
		{
			//hInvMass_temp->SetBinError(i,9999999999999); // setting huge errors outside the window peak
		}

		// Excluding mass peak window (setting errors to inf)
		if(hInvMass->GetBinCenter(i) > dMassWindow[0] && hInvMass->GetBinCenter(i) < dMassWindow[1])
		{
			hInvMass_side->SetBinError(i,9999999999999); 
		}
	}

	TF1* fFitInvMassPeak = new TF1("fFitInvMassPeak","gaus(0)+pol3(3)",dMassWindow[0],dMassWindow[1]);
	fFitInvMassPeak->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMassPeak->SetParameter(1,0.49);
	fFitInvMassPeak->SetParameter(2,0.005);
	fFitInvMassPeak->SetParLimits(2,0.,0.01);
	fFitInvMassPeak->SetLineWidth(2);
	fFitInvMassPeak->SetLineStyle(2);
	fFitInvMassPeak->SetLineColor(kPink+2);
	fFitInvMassPeak->SetNpx(1000000);

	TF1* fFitInvMassBg = new TF1("fFitInvMassBg","pol3",0.4,0.6);
	fFitInvMassBg->SetLineWidth(2);
	fFitInvMassBg->SetLineStyle(2);
	fFitInvMassBg->SetLineColor(kGreen+2);
	
	cCan->cd(3);
	//hInvMass_temp->SetMarkerColor(kRed);
	//hInvMass_temp->SetLineColor(kRed);
	hInvMass_side->SetMinimum(0);
	hInvMass_side->SetMaximum(100000);
	hInvMass_side->Fit("fFitInvMassBg","R");
	hInvMass_side->Draw("");

	cCan->cd(4);
	hInvMass_temp->SetMinimum(0);
	hInvMass_temp->SetMaximum(100000);
	hInvMass_temp->Fit("fFitInvMassPeak","R");
	hInvMass_temp->Draw();

	cCan->cd(1);
	fFitInvMassBg->Draw("same");
	fFitInvMassPeak->Draw("same");


	

	// estimating number of candidates
	const Int_t iNumMethod = 5; // 0 from approx fit / 1 signal from peak fit / 2 signal bin counting 
	TString sMethod[] = {"AproxFit","PeakFit","BinCounting","AproxFitIntegral","PeakFitIntegral"};
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

	dCandTot[iIndexMethod] = dCandTot[0];
	dCandBg[iIndexMethod] = dCandBg[0];
	dCandSig[iIndexMethod] = fFitInvMass->IntegralFast(nPoints,dX,dW,dMassWindow[0],dMassWindow[1]);// - dCandBg[iIndexMethod]; 

	iIndexMethod++;

	dCandTot[iIndexMethod] = dCandTot[0];
	dCandBg[iIndexMethod] = dCandBg[0];
	dCandSig[iIndexMethod] = fFitInvMassPeak->Integral(dMassWindow[0],dMassWindow[1]);// - dCandBg[iIndexMethod]; 

	iIndexMethod++;


	for(Int_t i(0); i < iNumMethod; i++)
	{
		printf("Cand::Method %s: Sig %g / Bg %g / Total %g (Sig+BG %g) / Sig/Bg %g\n",sMethod[i].Data(),dCandSig[i],dCandBg[i], dCandTot[i], dCandSig[i]+dCandBg[i],dCandSig[i]/dCandBg[i]);
	}
		// fitting flow mass dist
	cCan->cd(2);
	//TF1* fFitFlowMass = new TF1("fFitFlowMass","pol1(0)",0.4,0.6); 
	//hFlowMass->Fit("fFitFlowMass","R");


	return cCan;
	
}
