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
			hInvMassK0s[i][j] = (TH1D*) lInputInvMassK0s->FindObject(Form("hInvMass_K0s_Cent%d_pt%d",i,j));
			hInvMassLambda[i][j] = (TH1D*) lInputInvMassLambda->FindObject(Form("hInvMass_Lambda_Cent%d_pt%d",i,j));
			hFlowMassK0s[i][j] = (TH1D*) lInputFlowMassK0s->FindObject(Form("hFlowMass_K0s_Cent%d_pt%d",i,j));
			hFlowMassLambda[i][j] = (TH1D*) lInputFlowMassLambda->FindObject(Form("hFlowMass_Lambda_Cent%d_pt%d",i,j));
		}
	}

	// ===== Extracting signal & flow ===== 
	TCanvas* cCan;

	//cCan = ExtractFlow(hInvMassK0s[0][2],hFlowMassK0s[0][2]);
	//cCan->Print(Form("%s/fitK0s/test.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			cCan = ExtractFlow(hInvMassK0s[i][j],hFlowMassK0s[i][j]);
			cCan->Print(Form("%s/fitK0s/fit_K0s_Cent_%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
		}
	}

}

TCanvas* ExtractFlow(TH1D* hInvMass, TH1D* hFlowMass)
{
	if(!hInvMass || !hFlowMass)
	{
		printf("Invalid pointer to input histos!\n");
		return 0x0;
	}

	TCanvas* cCan = new TCanvas();
	cCan->Divide(2,1);
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
	hInvMass->Fit("fFitInvMass","R");

	// fitting flow mass dist
	cCan->cd(2);
	TF1* fFitFlowMass = new TF1("fFitFlowMass","pol1(0)",0.4,0.6); 
	hFlowMass->Fit("fFitFlowMass","R");


	return cCan;
}