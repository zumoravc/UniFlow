TCanvas* ExtractFlow(TH1D* hInvMass, TH1D* hFlowMass, Double_t* dConV2,const Short_t iPartSpecies = 0);

void V0sExtractFlow(
		const TString sInput = "~/NBI/Codes/results/V0s/8/plusplus_part1/plots/V0sFlow.root",
		const TString sOutput = "~/NBI/Codes/results/V0s/8/plusplus/plots",
		const TString sTag = "_JHEP",
		const TString sEtaGap = "Gap10",
		const TString sOutputFormat = "png"
	)
{


	const Int_t iNumPtBins = 22; // pT bins
	const Int_t iNumCentBins = 9; // centrality bins
	// bins edges
	Double_t fPtBinEdges[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4}; 
	Double_t fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

	// =======================================
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	TFile* fOutput = new TFile(Form("%s/../PtFlow_V0s_%s.root",sOutput.Data(),sEtaGap.Data()),"RECREATE");
	
	if(!fInput->IsOpen())
		return;

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

	/*
	// testing in for one distribution
	Double_t dVal = 0;
	cCan = ExtractFlow(hInvMassK0s[4][0],hFlowMassK0s[4][0],&dVal,0);
	cCan->Print(Form("%s/fitK0s/test.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());
	return;
	*/
	
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
			
			//cCan = ExtractFlow(hInvMassLambda[i][j],hFlowMassLambda[i][j],&dExtractedV2_Lambda[j],1);
			//cCan->Print(Form("%s/fitLambda/fit_Lambda_Cent_%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			//hFlowPt_Lambda[i]->SetBinContent(j+1,dExtractedV2_Lambda[j]);
		}

		// writing to output root file
		fOutput->cd();
		hFlowPt_K0s[i]->Write();
		hFlowPt_Lambda[i]->Write();

		cPtFlow->cd();
		hFlowPt_K0s[i]->Draw();
		cPtFlow->Print(Form("%s/PtFlow/flowPt_K0s_%s_Cent_%d.%s",sOutput.Data(),sEtaGap.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());

		cPtFlow->cd();
		hFlowPt_Lambda[i]->Draw();
		cPtFlow->Print(Form("%s/PtFlow/flowPt_Lambda_%s_Cent_%d.%s",sOutput.Data(),sEtaGap.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
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
			dFitRange[0] = 1.09;
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

	// ==== Fitting inv mass dist =====
	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass first shot\n");
	printf("-----------------------------------------\n");
	
	TF1* fFitInvMass = new TF1("fFitInvMass","gaus(0)+pol3(3)",dFitRange[0],dFitRange[1]); 
	fFitInvMass->SetParNames("Amp","Mean","Sigma");
	fFitInvMass->SetNpx(100000);
	fFitInvMass->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMass->SetParameter(1,dMassPeak);
	fFitInvMass->SetParameter(2,0.005);
	fFitInvMass->SetParLimits(2,0.,0.01);
	//fFitInvMass->SetLineColor(kGreen);
	fFitInvMass->SetLineWidth(2);
	//hInvMass->SetMinimum(0);
	//hInvMass->SetMaximum(4000);	
	hInvMass->Fit("fFitInvMass","R0");

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

	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass in sidebins\n");
	printf("-----------------------------------------\n");
	
	TF1* fFitInvMassBg = new TF1("fFitInvMassBg","pol3",dFitRange[0],dFitRange[1]);
	fFitInvMassBg->SetLineWidth(2);
	fFitInvMassBg->SetLineStyle(2);
	fFitInvMassBg->SetLineColor(kGreen+2);
	hInvMass_side->SetMinimum(0);
	hInvMass_side->SetMaximum(100000);
  hInvMass_side->Fit("fFitInvMassBg","R");
  	
	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		// subtracting the BG
		hInvMass_subs->SetBinContent(i,hInvMass->GetBinContent(i) - fFitInvMassBg->Eval(hInvMass_subs->GetBinCenter(i)));
	}

	TH1D* hInvMass_RatioSigTot = (TH1D*) hInvMass_subs->Clone("hInvMass_RatioSigTot");
	hInvMass_RatioSigTot->Divide(hInvMass);
	
	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass BG subtracted\n");
	printf("-----------------------------------------\n");
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
	
	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass Ratio Sig/Tot\n");
	printf("-----------------------------------------\n");
	

	TF1* fFitInvMassRatioSigTot = new TF1("fFitInvMassRatioSigTot","gaus(0)+gaus(3)+pol2(6)",dFitRange[0],dFitRange[1]);
	fFitInvMassRatioSigTot->SetParameter(0,0.7);
	fFitInvMassRatioSigTot->SetParameter(1,dMassPeak);
	fFitInvMassRatioSigTot->SetParLimits(1,dMassPeak-dMassPeakLimit,dMassPeak+dMassPeakLimit);
	fFitInvMassRatioSigTot->SetParameter(2,0.02);
	fFitInvMassRatioSigTot->SetParLimits(2,0.01,0.03);
	fFitInvMassRatioSigTot->SetParameter(3,0.7);
	fFitInvMassRatioSigTot->SetParameter(4,dMassPeak);
	fFitInvMassRatioSigTot->SetParLimits(4,dMassPeak-dMassPeakLimit,dMassPeak+dMassPeakLimit);
	fFitInvMassRatioSigTot->SetParameter(5,0.02);
	fFitInvMassRatioSigTot->SetParLimits(5,0.01,0.03);
	hInvMass_RatioSigTot->Fit("fFitInvMassRatioSigTot","R");

	// ==== Fitting flow mass dist ====
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
	printf("-----------------------------------------\n");
	printf("Fitting Flow Mass sidebands\n");
	printf("-----------------------------------------\n");
	
	TF1* fFitFlowMass_side = new TF1("fFitFlowMass_side","pol1(0)",dFitRange[0],dFitRange[1]);
	hFlowMass_side->Fit("fFitFlowMass_side","R0");

	// Fitting by vTot (LH side)

	printf("-----------------------------------------\n");
	printf("Fitting Flow Mass total \n");
	printf("-----------------------------------------\n");

	TF1* fFitFlowMassTot = new TF1("fFitFlowMassTot","[11]*(gaus(0)+gaus(3)+pol2(6))+(1-(gaus(0)+gaus(3)+pol2(6)))*pol1(9)",dFitRange[0],dFitRange[1]);
	// inv mass fit parameters
	for(Int_t i(0); i < 9; i++)
	{
		fFitFlowMassTot->FixParameter(i,fFitInvMassRatioSigTot->GetParameter(i));
	}
	// vn bg fit parameters
	fFitFlowMassTot->FixParameter(9,fFitFlowMass_side->GetParameter(0));
	fFitFlowMassTot->FixParameter(10,fFitFlowMass_side->GetParameter(1));
	//fFitFlowMassTot->FixParameter(11,fFitFlowMass_side->GetParameter(2));
	hFlowMass_temp->Fit("fFitFlowMassTot","R0");

	*dConV2 = fFitFlowMassTot->GetParameter(11);
	
	TLatex* latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);
	latex->SetTextSize(0.05);

	// drawing
	TCanvas* cCan = new TCanvas("cCan","cFit",2000,1000);
	cCan->Divide(3,1);
	cCan->cd(1);
	hInvMass->Draw();
	fFitInvMass->Draw("same");
	fFitInvMassBg->Draw("same");
	fFitInvMassSubs->Draw("same");
	cCan->cd(2);
	hInvMass_RatioSigTot->SetStats(0);
	hInvMass_RatioSigTot->Draw();
	fFitInvMassRatioSigTot->Draw("same");
	cCan->cd(3);
	hFlowMass->SetStats(0);
	hFlowMass->Draw();
	fFitFlowMassTot->Draw("same");
	latex->DrawLatex(0.2,0.2,Form("#it{v}_{2}: %g #pm %g",fFitFlowMassTot->GetParameter(11),fFitFlowMassTot->GetParError(11)));
	latex->DrawLatex(0.2,0.11,Form("Chi2/ndf: %g / %d",fFitFlowMassTot->GetChisquare(),fFitFlowMassTot->GetNDF()));
	
	return cCan;	
}
