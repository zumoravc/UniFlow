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

	TCanvas* cCan = new TCanvas();
	cCan->Divide(2,2);
	cCan->cd(1);
	hInvMass->Draw();
	cCan->cd(2);
	hFlowMass->Draw();

	// ===== Fitting K0s =====
	// fitting inv mass dist
	cCan->cd(1);
	TF1* fFitInvMass = new TF1("fFitInvMass","gaus(0)+pol3(3)",0.4,0.6); 
	fFitInvMass->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMass->SetParameter(1,0.49);
	fFitInvMass->SetParameter(2,0.01);
	fFitInvMass->SetLineWidth(1);
	hInvMass->Fit("fFitInvMass","R");

	// * TODO - fit verification procedure * //

	// extracting mean & sigma
	Double_t dMean = fFitInvMass->GetParameter(1);
	Double_t dSigma = fFitInvMass->GetParameter(2);

	Double_t dMassWindow[2] = {dMean - dNumSigma*dSigma, dMean + dNumSigma*dSigma}; // setting inv. mass peak
	printf("dMassWindow: %g - %g\n", dMassWindow[0],dMassWindow[1]);

	// Fitting only in the mass peak window

	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		if(hInvMass->GetBinCenter(i) < dMassWindow[0] || hInvMass->GetBinCenter(i) > dMassWindow[1])
		{
			hInvMass_temp->SetBinError(i,9999); // setting huge errors outside the window peak
		}

		if(hInvMass->GetBinCenter(i) > dMassWindow[0] && hInvMass->GetBinCenter(i) < dMassWindow[1])
		{
			hInvMass_side->SetBinError(i,999999999999); // setting huge errors outside the window peak
		}
	}

	
	cCan->cd(3);
	//hInvMass_temp->SetMarkerColor(kRed);
	//hInvMass_temp->SetLineColor(kRed);
	hInvMass_temp->SetMinimum(0);
	hInvMass_temp->Draw();



	// fitting clone of inv mass dist in peak window
	TF1* fFitInvMassPeak = (TF1*) fFitInvMass->Clone("fFitInvMassPeak");
	fFitInvMassPeak->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMassPeak->SetParameter(1,0.49);
	fFitInvMassPeak->SetParameter(2,0.01);
	fFitInvMassPeak->SetLineColor(kGreen);
	fFitInvMassPeak->SetLineWidth(1);
	hInvMass_temp->Fit("fFitInvMassPeak","LR");

	TF1* fFitInvMassBg = new TF1("fFitInvMassBg","pol3(3)",0.4,0.6);
	fFitInvMassBg->SetLineWidth(1);
	fFitInvMassBg->SetLineStyle(3);
	fFitInvMassBg->SetLineColor(kGreen);

	for(Int_t i(3); i < 7; i++)
	{
		fFitInvMassBg->SetParameter(i,fFitInvMassPeak->GetParameter(i));
	}

	cCan->cd(4);
	hInvMass_side->SetMaximum(5000);
	hInvMass_side->SetMinimum(0);
	hInvMass_side->Draw();
	// fitting clone of inv mass dist in sidebands
	TF1* fFitInvMassSide = new TF1("fFitInvMassSide","pol3(0)",0.4,0.6);
	fFitInvMassSide->SetLineColor(kBlack);
	fFitInvMassSide->SetLineWidth(1);
	fFitInvMassSide->SetLineStyle(2);
	hInvMass_side->Fit("fFitInvMassSide","LR");
	fFitInvMassSide->Draw("same");

	// estimating number of candidates
	const Int_t iNumMethod = 3; // 0 from fits / 1 signal bin counting // 2 sidebands BG fit
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
			dCandSig[iIndexMethod] += fFitInvMassPeak->Eval(hInvMass->GetBinCenter(i));
			dCandBg[iIndexMethod] += fFitInvMassBg->Eval(hInvMass->GetBinCenter(i));

			iIndexMethod++;

			dCandTot[iIndexMethod] += hInvMass->GetBinContent(i);
			dCandBg[iIndexMethod] += fFitInvMassBg->Eval(hInvMass->GetBinCenter(i));
			dCandSig[iIndexMethod] += hInvMass->GetBinContent(i) - fFitInvMassBg->Eval(hInvMass->GetBinCenter(i)); 
			
		}
	}

	for(Int_t i(0); i < iNumMethod; i++)
	{
		printf("Cand method %d: Sig %g / Bg %g / Total %g (Sig+BG %g)\n",i,dCandSig[i],dCandBg[i], dCandTot[i], dCandSig[i]+dCandBg[i]);
	}



	cCan->cd(1);
	fFitInvMassPeak->Draw("same");
	fFitInvMassBg->Draw("same");
	fFitInvMass->Draw("same");



	// fitting flow mass dist
	cCan->cd(2);
	TF1* fFitFlowMass = new TF1("fFitFlowMass","pol1(0)",0.4,0.6); 
	//hFlowMass->Fit("fFitFlowMass","R");


	return cCan;
}