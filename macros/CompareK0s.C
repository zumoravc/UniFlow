void CompareK0s()
{	
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	// loading TGraphs QM12 
	TFile* fInK0s1020 = new TFile("~/NBI/Flow/hepData/KsQM_1020.root","READ");
	TFile* fInK0s4050 = new TFile("~/NBI/Flow/hepData/KsQM_4050.root","READ");

	TH1D* hQMK0s[2];
	fInK0s1020->cd();
	hQMK0s[0] = (TH1D*) gDirectory->Get("hist")->Clone("hQMK0s_1020");

	fInK0s4050->cd();
	hQMK0s[1] = (TH1D*) gDirectory->Get("hist")->Clone("hQMK0s_4050");
	
	// loading my input
	TFile* fInMine = new TFile("~/NBI/Flow/results/V0s/13-QM12-check-2/PtFlow_V0s_Gap10.root","READ");
	fInMine->cd();

	TH1D* hK0s[2];
	hK0s[0] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent2")->Clone("hK0s_1020");
	hK0s[1] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent5")->Clone("hK0s_4050");


	TList* lCompK0s[2];
	TString sLabel[] = {"Mine","QM"};

	for(Int_t i(0); i < 2; i++)
	{
		hK0s[i]->SetStats(0);
		hK0s[i]->SetMinimum(0);
		hK0s[i]->SetMaximum(0.3);

		lCompK0s[i] = new TList();
		lCompK0s[i]->Add(hK0s[i]);
		lCompK0s[i]->Add(hQMK0s[i]);
	}

	/*
	TCanvas* cCompare = new TCanvas("cCompare","K0s QM12 comparison");
	cCompare->Divide(2,1);
	cCompare->cd(1);
	hK0s_1020->Draw();
	graphQMK0s_1020->Draw("same");
	cCompare->cd(2);
	hK0s_4050->Draw();
	graphQMK0s_4050->Draw("same");
	*/

	TCanvas* cCompRatio[2];
		//cCompRatio[1] = CompareRatio(hK0s[1],hQMK0s[1]);

	for(Int_t i(0); i < 2; i++)
	{
		//cCompRatio[i] = CompareHistos(lCompK0s[i],sLabel,0.,0.3,kFALSE);
		cCompRatio[i] = CompareRatio(hK0s[i],hQMK0s[i]);
	}
}