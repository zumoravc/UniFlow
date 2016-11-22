/*
TCanvas* CompareHistos(
	TList* lList,
	TString* sLabel,
	Double_t dMinY,
	Double_t dMaxY,
	Bool_t bLog
)
*/

void Compare(TH1D* hOne, TH1D* hTwo, TPad* pad);



void PlotGFK()
{
	
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C++");


	TString sOutPath = "~/NBI/Flow/temp/comp/";

	TFile* fInput = new TFile("~/NBI/Flow/temp/Flow-10samples.root","READ");
	fInput->cd();
	//fInput->ls();




	// ============================================================================

	// loading all lists in input file
	
	const Short_t iNumEtaGaps = 3;
	TString sEtaGaps[iNumEtaGaps] = {"-10","00","10"};

	const Short_t iNumHarmonics = 3;
	Short_t iHarmonics[iNumHarmonics] = {2,3,4};

	TH1D* hRef4[iNumHarmonics]; 
	TH1D* hRef2[iNumHarmonics][iNumEtaGaps];
	TList* lsTracks[iNumHarmonics][iNumEtaGaps];
	TList* lsPions[iNumHarmonics][iNumEtaGaps];
	TList* lsKaons[iNumHarmonics][iNumEtaGaps];
	TList* lsProtons[iNumHarmonics][iNumEtaGaps];

	for(Short_t iHarm(0); iHarm < iNumHarmonics; iHarm++)
	{
		hRef4[iHarm] = (TH1D*) gDirectory->Get(Form("hRef_n%d4_Gap-10",iHarmonics[iHarm]));
		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			hRef2[iHarm][iGap] = (TH1D*) gDirectory->Get(Form("hRef_n%d2_Gap%s",iHarmonics[iHarm],sEtaGaps[iGap].Data()));
			lsTracks[iHarm][iGap] = (TList*) gDirectory->Get(Form("Tracks_n%d_Gap%s",iHarmonics[iHarm],sEtaGaps[iGap].Data()));
			lsKaons[iHarm][iGap] = (TList*) gDirectory->Get(Form("Kaons_n%d_Gap%s",iHarmonics[iHarm],sEtaGaps[iGap].Data()));
			lsPions[iHarm][iGap] = (TList*) gDirectory->Get(Form("Pions_n%d_Gap%s",iHarmonics[iHarm],sEtaGaps[iGap].Data()));
			lsProtons[iHarm][iGap] = (TList*) gDirectory->Get(Form("Protons_n%d_Gap%s",iHarmonics[iHarm],sEtaGaps[iGap].Data()));
		}
	}

	// ######## Tlists loaded ##################

	// ============================================================================
	
	// pT diff v22
	// Loading You

	TFile* fYouPtV2 = new TFile("~/NBI/Flow/results/you/Analysis_v2pt.root","READ");
	TFile* fYouPtV3 = new TFile("~/NBI/Flow/results/you/Analysis_v3pt.root","READ");
	TFile* fYouPtV4 = new TFile("~/NBI/Flow/results/you/Analysis_v4pt.root","READ");
	//fYouPtV2->cd();
	//fYouPtV2->ls();
	//fYouPtV3->ls();

	const Short_t iNumCent = 9;
	TString sYouLabel[] = {"NoGap","Gap00","Gap10"};
	TH1D* hYouPtV22[iNumEtaGaps][iNumCent];
	TH1D* hYouPtV32[iNumEtaGaps][iNumCent];
	TH1D* hYouPtV42[iNumEtaGaps][iNumCent];
	TH1D* hYouPtV24[iNumCent];
	TH1D* hYouPtV34[iNumCent];
	TH1D* hYouPtV44[iNumCent];
	
	TH1D* hMinePtV22[iNumEtaGaps][iNumCent];
	TH1D* hMinePtV32[iNumEtaGaps][iNumCent];
	TH1D* hMinePtV42[iNumEtaGaps][iNumCent];
	TH1D* hMinePtV24[iNumCent];
	TH1D* hMinePtV34[iNumCent];
	TH1D* hMinePtV44[iNumCent];

	//lsTracks[0][0]->ls();

	// loading
	for(Short_t iCent(0); iCent < iNumCent; iCent++)
	{
		printf("You: %s\n",Form("hisv24_cent%d",iCent));
		hYouPtV24[iCent] = (TH1D*) fYouPtV2->Get(Form("hisv24_cent%d",iCent));
		if(!hYouPtV24[iCent]) return;
		printf("34\n");
		hYouPtV34[iCent] = (TH1D*) fYouPtV3->Get(Form("hisv34_cent%d",iCent));
		if(!hYouPtV34[iCent]) return;
		printf("44\n");
		hYouPtV44[iCent] = (TH1D*) fYouPtV4->Get(Form("hisv44_cent%d",iCent));
		if(!hYouPtV44[iCent]) return;

		printf("Mine: %s\n", Form("fTracks_n24_gap-10_cent%d_number0_px_desampled",iCent));
		hMinePtV24[iCent] = (TH1D*) lsTracks[0][0]->FindObject(Form("fTracks_n24_gap-10_cent%d_number0_px_desampled",iCent));
		if(!hMinePtV24[iCent]) return;
		printf("34\n");
		hMinePtV34[iCent] = (TH1D*) lsTracks[1][0]->FindObject(Form("fTracks_n34_gap-10_cent%d_number0_px_desampled",iCent));
		if(!hMinePtV34[iCent]) return;
		printf("44\n");
		hMinePtV44[iCent] = (TH1D*) lsTracks[2][0]->FindObject(Form("fTracks_n44_gap-10_cent%d_number0_px_desampled",iCent));
		if(!hMinePtV44[iCent]) return;

		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			printf("You: %s\n",Form("hisv22%s_cent%d",sYouLabel[iGap].Data(),iCent));
			hYouPtV22[iGap][iCent] = (TH1D*) fYouPtV2->Get(Form("hisv22%s_cent%d",sYouLabel[iGap].Data(),iCent));
			if(!hYouPtV22[iGap][iCent]) return;
			hYouPtV32[iGap][iCent] = (TH1D*) fYouPtV3->Get(Form("hisv32%s_cent%d",sYouLabel[iGap].Data(),iCent));
			if(!hYouPtV32[iGap][iCent]) return;
			hYouPtV42[iGap][iCent] = (TH1D*) fYouPtV4->Get(Form("hisv42%s_cent%d",sYouLabel[iGap].Data(),iCent));
			if(!hYouPtV42[iGap][iCent]) return;	
			
			printf("Mine: %s\n", Form("fTracks_n22_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
			hMinePtV22[iGap][iCent] = (TH1D*) lsTracks[0][iGap]->FindObject(Form("fTracks_n22_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
			if(!hMinePtV22[iGap][iCent]) return;
			hMinePtV32[iGap][iCent] = (TH1D*) lsTracks[1][iGap]->FindObject(Form("fTracks_n32_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
			if(!hMinePtV32[iGap][iCent]) return;
			hMinePtV42[iGap][iCent] = (TH1D*) lsTracks[2][iGap]->FindObject(Form("fTracks_n42_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
			if(!hMinePtV42[iGap][iCent]) return;
		}
	}

	// plotting v22, v24
	TCanvas* cTemp = new TCanvas("cTemp");
	cTemp->cd();

	for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
	{
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			hMinePtV22[iGap][iCent]->Draw();
			hYouPtV22[iGap][iCent]->Draw("same");
			cTemp->SaveAs(Form("%s/PtV22_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
		}

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			hMinePtV32[iGap][iCent]->Draw();
			hYouPtV32[iGap][iCent]->Draw("same");
			cTemp->SaveAs(Form("%s/PtV32_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
		}

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			hMinePtV42[iGap][iCent]->Draw();
			hYouPtV42[iGap][iCent]->Draw("same");
			cTemp->SaveAs(Form("%s/PtV42_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
		}
	}

	for(Short_t iCent(0); iCent < iNumCent; iCent++)
	{
		hMinePtV24[iCent]->Draw();
		hYouPtV24[iCent]->Draw("same");
		cTemp->SaveAs(Form("%s/PtV24_cent%d_comp.pdf",sOutPath.Data(),iCent));
	}

	for(Short_t iCent(0); iCent < iNumCent; iCent++)
	{
		hMinePtV34[iCent]->Draw();
		hYouPtV34[iCent]->Draw("same");
		cTemp->SaveAs(Form("%s/PtV34_cent%d_comp.pdf",sOutPath.Data(),iCent));
	}

	
	for(Short_t iCent(0); iCent < iNumCent; iCent++)
	{
		hMinePtV44[iCent]->Draw();
		//hYouPtV44[iCent]->Draw("same");
		cTemp->SaveAs(Form("%s/PtV44_cent%d_comp.pdf",sOutPath.Data(),iCent));
	}










	return;
	// ============================================================================

	// Plotting stuff

	// reference flow
	TList* lsTempRef = new TList();
	TString sLabelRef[iNumHarmonics][iNumEtaGaps];


	for(Short_t iHarm(0); iHarm < iNumHarmonics; iHarm++)
	{
		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			lsTempRef->Add(hRef2[iHarm][iGap]);
			sLabelRef[iHarm][iGap] = Form("v%d2 Gap%s",iHarmonics[iHarm], sEtaGaps[iGap].Data());

			// plotting

		}
	}

	//lsTempRef->ls();

	//TCanvas* cRef = CompareHistos(lsTempRef,sLabelRef,-1,0.5,kFALSE);

	// ============================================================================

	// Comparing with published / preliminary / You results

	TFile* fYouRef = new TFile("~/NBI/Flow/results/you/Analysis_vn_2468.root","READ");
	fYouRef->cd();
	//fYouRef->ls();

	// reference flow
	TH1D* fv22NoGap = (TH1D*) gDirectory->Get("fv22NoGap");
	TH1D* fv22Gap00 = (TH1D*) gDirectory->Get("fv22Gap00");
	TH1D* fv22Gap04 = (TH1D*) gDirectory->Get("fv22Gap04");
	TH1D* fv22Gap08 = (TH1D*) gDirectory->Get("fv22Gap08");
	TH1D* fv22Gap10 = (TH1D*) gDirectory->Get("fv22Gap10");

	//=========================================================================================================================
	//  v2:
	// QC2_v2 = v2{2}:
	Double_t xQC2_v2[] = {2.5,7.5,15.,25.,35.,45.,55.,65.,75.};
	Double_t yQC2_v2[] = {0.03057544, 0.04782869, 0.06730953, 0.08693288, 0.09905272, 
	    0.10507712, 0.10546323, 0.10293218, 0.10292639};
	Double_t xErrQC2_v2[9] = {0.};
	Double_t yErrQC2_v2[] = {0.00006919, 0.00008033, 0.00006915, 0.00008314, 0.00010003,
	    0.00012251, 0.00015668, 0.00022385, 0.00039746};
	Int_t nPointsQC2_v2 = sizeof(xQC2_v2)/sizeof(Double_t);         
	TGraphErrors *QC2_v2 = new TGraphErrors(nPointsQC2_v2,xQC2_v2,yQC2_v2,xErrQC2_v2,yErrQC2_v2);

 	//  v2{SP}(pt) for 30-40%, rapidity gap = 1.0:
    const Int_t nPointsSP_3040ALICE_v2_etaGap10 = 22;
    Double_t xSP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.250000,0.350000,0.450000,0.550000,0.650000,
        0.750000,0.850000,0.950000,1.100000,1.300000,1.500000,1.700000,1.900000,2.100000,2.300000,2.500000,2.800000,3.200000,
        3.600000,4.000000,4.400000,4.800000};
    Double_t ySP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.039107,0.053955,0.068788,0.083037,0.095749,
        0.107080,0.119547,0.127823,0.143002,0.160203,0.180128,0.191935,0.201014,0.210955,0.220726,0.227580,0.244640,0.230650,
        0.234024,0.220566,0.217360,0.208241};
    Double_t xErrSP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.};
    Double_t yErrSP_3040ALICE_v2_etaGap10[nPointsSP_3040ALICE_v2_etaGap10] = {0.000746,0.000744,0.000790,0.000867,0.000960,
        0.001069,0.001187,0.001321,0.001091,0.001339,0.001642,0.002016,0.002483,0.003050,0.003736,0.004558,0.004267,0.006108,
        0.008517,0.011491,0.017248,0.017496};
    TGraphErrors *GrSP_3040ALICE_v2_etaGap10 = new TGraphErrors(nPointsSP_3040ALICE_v2_etaGap10,xSP_3040ALICE_v2_etaGap10,
                                                              ySP_3040ALICE_v2_etaGap10,xErrSP_3040ALICE_v2_etaGap10,
                                                              yErrSP_3040ALICE_v2_etaGap10);

	TCanvas* cYouRef = new TCanvas("cYouRef","YouRef");
	cYouRef->cd();
	//fv22NoGap->Draw();
	//fv22Gap00->Draw("same");
	hRef2[0][3]->Draw();
	fv22Gap10->Draw("same");
	//fv22Gap00->Draw();
	//hRef2[0][0]->Draw("same");
	//hRef2[0][1]->Draw("same");
	QC2_v2->Draw("same");
	//hRef2[0][0]->Draw("same");
	//Compare(fv22Gap00,hRef2[0][1],gPad);

/*
	TCanvas* cYouRefRatio = new TCanvas();
	cYouRefRatio->cd();
	
	TH1D* hRatioNoGap = (TH1D*) fv22NoGap->Clone("RatioNoGap");
	hRatioNoGap->Divide(hRef2[0][0]);
	hRatioNoGap->Draw();
*/

	// pt diff
	TFile* fYouDiff2 = new TFile("~/NBI/Flow/results/you/Analysis_v2pt.root","READ");
	fYouDiff2->cd();
	fYouDiff2->ls();
	

	TH1D* hisv22Gap10_cent4 = (TH1D*) gDirectory->Get("hisv22Gap10_cent4");
	lsTracks[0][0]->ls();
	TH1D* hTracks22Gap10_cent4 = (TH1D*) lsTracks[0][3]->FindObject("fTracks_n22_gap10_cent4_number0_0_px_desampled");

	TCanvas* cYouDiff2 = new TCanvas("cYouDiff2","YouDiff2");
	cYouDiff2->cd();
	hTracks22Gap10_cent4->Draw("l");
	hisv22Gap10_cent4->Draw("same");
	GrSP_3040ALICE_v2_etaGap10->Draw("same");


	return;
}

