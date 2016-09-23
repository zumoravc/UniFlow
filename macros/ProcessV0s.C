void ProcessV0s()
{
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareHistos.C");
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareRatio.C");
	
	TString sInput = "~/NBI/Codes/flow/AnalysisResults.root";

	TFile* fInput = new TFile(sInput.Data(),"READ");
	fInput->cd("FlowPID");

	const Int_t iNumCentBins = 9; // centrality bins
	const Int_t iNumPtBins = 9; // pT bins

	// ===== Loading input ===== 
	TList* lInputTracks = (TList*) gDirectory->Get("Tracks");
	TList* lInputV0s = (TList*) gDirectory->Get("V0s");

	// reference 
	TProfile* pRefCorTwo2_Gap09 = (TProfile*) lInputTracks->FindObject("fRefCorTwo2_Gap09"); 
	
	// V0s 
	TProfile2D* p2V0sDiffTwo2_Gap09P_K0s[iNumCentBins];
	TProfile2D* p2V0sDiffTwo2_Gap09N_K0s[iNumCentBins];
	TProfile2D* p2V0sDiffTwo2_Gap09P_Lambda[iNumCentBins];
	TProfile2D* p2V0sDiffTwo2_Gap09N_Lambda[iNumCentBins];

	TH2D* h2InvMass_Gap09_K0s[iNumCentBins];
	TH2D* h2InvMass_Gap09_Lambda[iNumCentBins];
	
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		h2InvMass_Gap09_K0s[i] = (TH2D*) lInputV0s->FindObject(Form("fV0sK0s_Gap09_Cent%d",i))->Clone(Form("h2InvMass_Gap09_K0s_Cent%d",i));
		h2InvMass_Gap09_Lambda[i] = (TH2D*) lInputV0s->FindObject(Form("fV0sLambda_Gap09_Cent%d",i))->Clone(Form("h2InvMass_Gap09_Lambda_Cent%d",i));

		p2V0sDiffTwo2_Gap09P_K0s[i] = (TProfile2D*) lInputV0s->FindObject(Form("fV0sDiffTwo2_Gap09P_K0s_Cent%d",i))->Clone(Form("p2V0sDiffTwo2_Gap09P_K0s_Cent%d",i));
		p2V0sDiffTwo2_Gap09N_K0s[i] = (TProfile2D*) lInputV0s->FindObject(Form("fV0sDiffTwo2_Gap09N_K0s_Cent%d",i))->Clone(Form("p2V0sDiffTwo2_Gap09N_K0s_Cent%d",i));
	
		p2V0sDiffTwo2_Gap09P_Lambda[i] = (TProfile2D*) lInputV0s->FindObject(Form("fV0sDiffTwo2_Gap09P_Lambda_Cent%d",i))->Clone(Form("p2V0sDiffTwo2_Gap09P_Lambda_Cent%d",i));
		p2V0sDiffTwo2_Gap09P_Lambda[i] = (TProfile2D*) lInputV0s->FindObject(Form("fV0sDiffTwo2_Gap09N_Lambda_Cent%d",i))->Clone(Form("p2V0sDiffTwo2_Gap09N_Lambda_Cent%d",i));
	}

	// ===== Making projections ======
	TH1D* hInvMass_K0s[iNumCentBins][iNumPtBins];
	TH1D* hInvMass_Lambda[iNumCentBins][iNumPtBins];

	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			hInvMass_K0s[i][j] = 
		}
	}



}