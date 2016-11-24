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


/*
TGraphErrors* Kv2_1020_QC2(Int_t color=1, Int_t marker=20);
TGraphAsymmErrors* v2Pion1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1);
TGraphAsymmErrors* v2Kaon1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1);
TGraphAsymmErrors* v2Antiproton1020(Int_t color=1, Int_t marker=20, Int_t first=-1,Int_t last=-1);
TGraphErrors* Kv2_1020_QC4(Int_t color=1, Int_t marker=20);
*/

void PlotGFK(TString sOutPath = "~/NBI/Flow/temp/comp/", TString sInputFile = "~/NBI/Flow/temp/Flow-10samples.root")
{
	
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C++");

	Bool_t bYou = kTRUE;
	Bool_t bYouRatio = kTRUE;
	Bool_t bPID = kFALSE;


	//TString sOutPath = "~/NBI/Flow/temp/comp/";

	TFile* fInput = new TFile(sInputFile.Data(),"READ");
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

		const Short_t iNumCent = 9;
	// ============================================================================
	// Comparison with YOU 

	if(bYou)
	{

		TString sYouLabel[] = {"NoGap","Gap00","Gap10"};
		// reference 
		
		TFile* fYouRef= new TFile("~/NBI/Flow/results/you/Analysis_vn_2468.root","READ");
		fYouRef->ls();

		TH1D* hYouRef22[iNumEtaGaps];
		TH1D* hYouRef32[iNumEtaGaps];
		TH1D* hYouRef42[iNumEtaGaps];
		TH1D* hYouRef24 = (TH1D*) fYouRef->Get("fv24");

		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			hYouRef22[iGap] = (TH1D*) fYouRef->Get(Form("fv22%s",sYouLabel[iGap].Data()));
			if(!hYouRef22[iGap]) return;
			hYouRef32[iGap] = (TH1D*) fYouRef->Get(Form("fv32%s",sYouLabel[iGap].Data()));
			if(!hYouRef32[iGap]) return;
			hYouRef42[iGap] = (TH1D*) fYouRef->Get(Form("fv42%s",sYouLabel[iGap].Data()));
			if(!hYouRef42[iGap]) return;
		}

		TLine* unity = new TLine(0.,1.,80.,1.);
		TCanvas* cRatio = new TCanvas("cRatio");
		cRatio->Divide(2,1);
		TH1D* hTempRatio = 0x0;
		TH1D* hDummy = new TH1D();

		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			cRatio->cd(1);
			hYouRef22[iGap]->Draw();
			hYouRef22[iGap]->SetMinimum(0.);
			hRef2[0][iGap]->Draw("same");


			cRatio->cd(2);
			hTempRatio = (TH1D*) hYouRef22[iGap]->Clone("hTempRatio");
			hTempRatio->Divide(hRef2[0][iGap]);
			hTempRatio->SetMinimum(0.8);
			hTempRatio->SetMaximum(1.2);
			hTempRatio->Draw();
			unity->Draw("same");
			cRatio->SaveAs(Form("%s/Ref/Ref_v22_Gap%s_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data()));
		}


		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			cRatio->cd(1);
			hYouRef32[iGap]->Draw();
			hYouRef32[iGap]->SetMinimum(0.);
			hRef2[1][iGap]->Draw("same");

			
			cRatio->cd(2);
			if(iGap == 2) 
			{
				hDummy->Draw();
				continue; 
			}
			hTempRatio = (TH1D*) hYouRef32[iGap]->Clone("hTempRatio");
			hTempRatio->Divide(hRef2[1][iGap]);
			hTempRatio->SetMinimum(0.8);
			hTempRatio->SetMaximum(1.2);
			hTempRatio->Draw();
			unity->Draw("same");
			
			cRatio->SaveAs(Form("%s/Ref/Ref_v32_Gap%s_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data()));
		}


		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			cRatio->cd(1);
			hYouRef42[iGap]->Draw();
			hYouRef42[iGap]->SetMinimum(0.);
			hRef2[2][iGap]->Draw("same");

			
			cRatio->cd(2);
			if(iGap == 2) 
			{
				hDummy->Draw();
				continue; 
			}
			hTempRatio = (TH1D*) hYouRef42[iGap]->Clone("hTempRatio");
			hTempRatio->Divide(hRef2[2][iGap]);
			hTempRatio->SetMinimum(0.8);
			hTempRatio->SetMaximum(1.2);
			hTempRatio->Draw();
			unity->Draw("same");

			cRatio->SaveAs(Form("%s/Ref/Ref_v42_Gap%s_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data()));
		}

		cRatio->cd(1);
		hYouRef24->Draw();
		hYouRef24->SetMinimum(0.);
		hRef4[0]->Draw("same");

		cRatio->cd(2);
		hTempRatio = (TH1D*) hYouRef24->Clone("hTempRatio");
		hTempRatio->Divide(hRef4[0]);
		hTempRatio->SetMinimum(0.8);
		hTempRatio->SetMaximum(1.2);
		hTempRatio->Draw();
		unity->Draw("same");
		cRatio->SaveAs(Form("%s/Ref/Ref_v24_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data()));

		// pT diff v22
		// Loading You

		TCanvas* cTemp = new TCanvas("cTemp");
		TFile* fYouPtV2 = new TFile("~/NBI/Flow/results/you/Analysis_v2pt.root","READ");
		TFile* fYouPtV3 = new TFile("~/NBI/Flow/results/you/Analysis_v3pt.root","READ");
		TFile* fYouPtV4 = new TFile("~/NBI/Flow/results/you/Analysis_v4pt.root","READ");
		//fYouPtV2->cd();
		//fYouPtV2->ls();
		//fYouPtV3->ls();

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
			//printf("You: %s\n",Form("hisv24_cent%d",iCent));
			hYouPtV24[iCent] = (TH1D*) fYouPtV2->Get(Form("hisv24_cent%d",iCent));
			if(!hYouPtV24[iCent]) return;
			//printf("34\n");
			hYouPtV34[iCent] = (TH1D*) fYouPtV3->Get(Form("hisv34_cent%d",iCent));
			if(!hYouPtV34[iCent]) return;
			//printf("44\n");
			hYouPtV44[iCent] = (TH1D*) fYouPtV4->Get(Form("hisv44_cent%d",iCent));
			if(!hYouPtV44[iCent]) return;

			//printf("Mine: %s\n", Form("fTracks_n24_gap-10_cent%d_number0_px_desampled",iCent));
			hMinePtV24[iCent] = (TH1D*) lsTracks[0][0]->FindObject(Form("fTracks_n24_gap-10_cent%d_number0_px_desampled",iCent));
			if(!hMinePtV24[iCent]) return;
			//printf("34\n");
			hMinePtV34[iCent] = (TH1D*) lsTracks[1][0]->FindObject(Form("fTracks_n34_gap-10_cent%d_number0_px_desampled",iCent));
			if(!hMinePtV34[iCent]) return;
			//printf("44\n");
			hMinePtV44[iCent] = (TH1D*) lsTracks[2][0]->FindObject(Form("fTracks_n44_gap-10_cent%d_number0_px_desampled",iCent));
			if(!hMinePtV44[iCent]) return;

			for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
			{
				//printf("You: %s\n",Form("hisv22%s_cent%d",sYouLabel[iGap].Data(),iCent));
				hYouPtV22[iGap][iCent] = (TH1D*) fYouPtV2->Get(Form("hisv22%s_cent%d",sYouLabel[iGap].Data(),iCent));
				if(!hYouPtV22[iGap][iCent]) return;
				hYouPtV32[iGap][iCent] = (TH1D*) fYouPtV3->Get(Form("hisv32%s_cent%d",sYouLabel[iGap].Data(),iCent));
				if(!hYouPtV32[iGap][iCent]) return;
				hYouPtV42[iGap][iCent] = (TH1D*) fYouPtV4->Get(Form("hisv42%s_cent%d",sYouLabel[iGap].Data(),iCent));
				if(!hYouPtV42[iGap][iCent]) return;	
				
				//printf("Mine: %s\n", Form("fTracks_n22_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
				hMinePtV22[iGap][iCent] = (TH1D*) lsTracks[0][iGap]->FindObject(Form("fTracks_n22_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
				if(!hMinePtV22[iGap][iCent]) return;
				hMinePtV32[iGap][iCent] = (TH1D*) lsTracks[1][iGap]->FindObject(Form("fTracks_n32_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
				if(!hMinePtV32[iGap][iCent]) return;
				hMinePtV42[iGap][iCent] = (TH1D*) lsTracks[2][iGap]->FindObject(Form("fTracks_n42_gap%s_cent%d_number0_0_px_desampled",sEtaGaps[iGap].Data(),iCent));
				if(!hMinePtV42[iGap][iCent]) return;
			}
		}

		printf("loaded\n");
		// plotting v22, v24

		TLegend* legTracks = new TLegend(0.12,0.6,0.3,0.89);
		legTracks->SetBorderSize(0);
		legTracks->AddEntry(hMinePtV22[0][0], "Mine","pel");
		//legTracks->AddEntry(hYouPtV22[0][0], "You","pel");
		
		TLatex latexTracks;
		latexTracks.SetNDC();
		
		for(Short_t iGap(0); iGap < iNumEtaGaps; iGap++)
		{
			for(Short_t iCent(0); iCent < iNumCent; iCent++)
			{
				cTemp->cd();
				hMinePtV22[iGap][iCent]->Draw();
				hYouPtV22[iGap][iCent]->Draw("same");
				legTracks->Draw();
				cTemp->SaveAs(Form("%s/Tracks/PtV22_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
				
				if(bYouRatio)
				{
					hTempRatio = (TH1D*) hYouPtV22[iGap][iCent]->Clone("hTempRatio");
					hTempRatio->Divide(hMinePtV22[iGap][iCent]);
					hTempRatio->SetMinimum(0.8);
					hTempRatio->SetMaximum(1.2);

					cRatio->cd(1);
					hMinePtV22[iGap][iCent]->Draw();
					hYouPtV22[iGap][iCent]->Draw("same");
					legTracks->Draw();

					cRatio->cd(2);
					hTempRatio->Draw();
					unity->Draw("same");
					cRatio->SaveAs(Form("%s/Tracks/Ratios/RatioPtV22_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
				}
			}

			for(Short_t iCent(0); iCent < iNumCent; iCent++)
			{
				cTemp->cd();
				hMinePtV32[iGap][iCent]->Draw();
				hYouPtV32[iGap][iCent]->Draw("same");
				legTracks->Draw();
				cTemp->SaveAs(Form("%s/Tracks/PtV32_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));

				if(bYouRatio)
				{
					hTempRatio = (TH1D*) hYouPtV32[iGap][iCent]->Clone("hTempRatio");
					hTempRatio->Divide(hMinePtV32[iGap][iCent]);
					hTempRatio->SetMinimum(0.8);
					hTempRatio->SetMaximum(1.2);

					cRatio->cd(1);
					hMinePtV32[iGap][iCent]->Draw();
					hYouPtV32[iGap][iCent]->Draw("same");
					legTracks->Draw();

					cRatio->cd(2);
					hTempRatio->Draw();
					unity->Draw("same");
					cRatio->SaveAs(Form("%s/Tracks/Ratios/RatioPtV32_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
				}
			}

			for(Short_t iCent(0); iCent < iNumCent; iCent++)
			{
				cTemp->cd();
				hMinePtV42[iGap][iCent]->Draw();
				hYouPtV42[iGap][iCent]->Draw("same");
				legTracks->Draw();
				cTemp->SaveAs(Form("%s/Tracks/PtV42_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));

				if(bYouRatio)
				{
					hTempRatio = (TH1D*) hYouPtV42[iGap][iCent]->Clone("hTempRatio");
					hTempRatio->Divide(hMinePtV42[iGap][iCent]);
					hTempRatio->SetMinimum(0.8);
					hTempRatio->SetMaximum(1.2);

					cRatio->cd(1);
					hMinePtV42[iGap][iCent]->Draw();
					hYouPtV42[iGap][iCent]->Draw("same");
					legTracks->Draw();

					cRatio->cd(2);
					hTempRatio->Draw();
					unity->Draw("same");
					cRatio->SaveAs(Form("%s/Tracks/Ratios/RatioPtV42_Gap%s_cent%d_comp.pdf",sOutPath.Data(),sEtaGaps[iGap].Data(),iCent));
				}
			}
		}

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			cTemp->cd();
			hMinePtV24[iCent]->Draw();
			hYouPtV24[iCent]->Draw("same");
			legTracks->Draw();
			cTemp->SaveAs(Form("%s/Tracks/PtV24_cent%d_comp.pdf",sOutPath.Data(),iCent));
		}

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			cTemp->cd();
			hMinePtV34[iCent]->Draw();
			hYouPtV34[iCent]->Draw("same");
			legTracks->Draw();
			cTemp->SaveAs(Form("%s/Tracks/PtV34_cent%d_comp.pdf",sOutPath.Data(),iCent));
		}

		
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			cTemp->cd();
			hMinePtV44[iCent]->Draw();
			//hYouPtV44[iCent]->Draw("same");
			legTracks->Draw();
			cTemp->SaveAs(Form("%s/Tracks/PtV44_cent%d_comp.pdf",sOutPath.Data(),iCent));
		}
	}

	delete cTemp;
	delete cRatio;

	// ============================================================================
	// Comparison of pi, K, p



	if(bPID)
	{

		gROOT->LoadMacro("~/NBI/Flow/results/pid/v2All.C");
		
		// loading PID 
		Int_t iCentBins[10] = {0,5,10,20,30,40,50,60,70,80};

		TGraph* PionV2[iNumCent];
		TGraph* KaonV2[iNumCent];
		TGraph* ProtonV2[iNumCent];
		
		TGraph* PionV3[iNumCent];
		TGraph* KaonV3[iNumCent];
		TGraph* ProtonV3[iNumCent];
		
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			PionV2[iCent] = 0x0;
			PionV3[iCent] = 0x0;
			KaonV2[iCent] = 0x0;
			KaonV3[iCent] = 0x0;
			ProtonV2[iCent] = 0x0;
			ProtonV3[iCent] = 0x0;
		}


		PionV3[0] = v3Pion0005(1,20,4,36);
		PionV3[1] = v3Pion0510(1,20,4,36);
		PionV3[2] = v3Pion1020(1,20,4,36);
		PionV3[3] = v3Pion2030(1,20,4,36);
		PionV3[4] = v3Pion3040(1,20,4,36);
		PionV3[5] = v3Pion4050(1,20,4,36);
		PionV3[6] = v3Pion5060(1,20,4,36);
		PionV3[7] = v3Pion6070(1,20,4,36);
		//PionV3[8] = v3Pion7080(1,20,4,36);
		
		KaonV3[0] = v3Kaon0005(1,20,4,36);
		KaonV3[1] = v3Kaon0510(1,20,4,36);
		KaonV3[2] = v3Kaon1020(1,20,4,36);
		KaonV3[3] = v3Kaon2030(1,20,4,36);
		KaonV3[4] = v3Kaon3040(1,20,4,36);
		KaonV3[5] = v3Kaon4050(1,20,4,36);
		KaonV3[6] = v3Kaon5060(1,20,4,36);
		KaonV3[7] = v3Kaon6070(1,20,4,36);
		//KaonV3[8] = v3Kaon7080(1,20,4,36);

		ProtonV3[0] = v3Antiproton0005(1,20,4,36);
		ProtonV3[1] = v3Antiproton0510(1,20,4,36);
		ProtonV3[2] = v3Antiproton1020(1,20,4,36);
		ProtonV3[3] = v3Antiproton2030(1,20,4,36);
		ProtonV3[4] = v3Antiproton3040(1,20,4,36);
		ProtonV3[5] = v3Antiproton4050(1,20,4,36);
		ProtonV3[6] = v3Antiproton5060(1,20,4,36);
		ProtonV3[7] = v3Antiproton6070(1,20,4,36);
		//ProtonV3[8] = v3Antiproton7080(1,20,4,36);
		
		TH1D* hMinePionV2[iNumCent];
		TH1D* hMineKaonV2[iNumCent];
		TH1D* hMineProtonV2[iNumCent];
		TH1D* hMinePionV32[iNumCent];
		TH1D* hMineKaonV32[iNumCent];
		TH1D* hMineProtonV32[iNumCent];
		
		// mikolaj QM11

		TGraph* PionV22[iNumCent];
		TGraph* KaonV22[iNumCent];
		TGraph* ProtonV22[iNumCent];

		TGraph* PionV32[iNumCent];
		TGraph* KaonV32[iNumCent];
		TGraph* ProtonV32[iNumCent];

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			PionV22[iCent] = 0x0;
			PionV32[iCent] = 0x0;
			KaonV22[iCent] = 0x0;
			KaonV32[iCent] = 0x0;
			ProtonV22[iCent] = 0x0;
			ProtonV32[iCent] = 0x0;
		}

		PionV22[2] = v22_etagap10_1020_pion(2,20,3,31);
    PionV22[3] = v22_etagap10_2030_pion(2,20,3,31);
    PionV22[4] = v22_etagap10_3040_pion(2,20,3,31);
    PionV22[5] = v22_etagap10_4050_pion(2,20,3,31);
    PionV22[6] = v22_etagap10_5060_pion(2,20,3,31);

    
		KaonV22[2] = v22_etagap10_1020_kaon(2,20,4,23);
    KaonV22[3] = v22_etagap10_2030_kaon(2,20,4,23);
    KaonV22[4] = v22_etagap10_3040_kaon(2,20,4,23);
    KaonV22[5] = v22_etagap10_4050_kaon(2,20,4,23);
    KaonV22[6] = v22_etagap10_5060_kaon(2,20,4,23);

		ProtonV22[2] = v22_etagap10_1020_antiproton(2,20,5,35);
    ProtonV22[3] = v22_etagap10_2030_antiproton(2,20,5,35);
    ProtonV22[4] = v22_etagap10_3040_antiproton(2,20,5,35);
    ProtonV22[5] = v22_etagap10_4050_antiproton(2,20,5,35);
    ProtonV22[6] = v22_etagap10_5060_antiproton(2,20,5,35);

    PionV32[2] = v32_1020_pion(2,20);
    PionV32[5] = v32_4050_pion(2,20);
    KaonV32[2] = v32_1020_kaon(2,20);
    KaonV32[5] = v32_4050_kaon(2,20);
    ProtonV32[2] = v32_1020_antiproton(2,20);
    ProtonV32[5] = v32_4050_antiproton(2,20);

    TString sErrorOption = TString("P");


		//lsPions[0][3]->ls();
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			PionV2[iCent] = GetV2(1,iCentBins[iCent],iCentBins[iCent+1]);
			KaonV2[iCent] = GetV2(2,iCentBins[iCent],iCentBins[iCent+1]);
			ProtonV2[iCent] = GetV2(3,iCentBins[iCent],iCentBins[iCent+1]);

			hMinePionV2[iCent] = (TH1D*) lsPions[0][2]->FindObject(Form("fPion_n22_gap10_cent%d_number0_0_px_desampled",iCent));
			hMinePionV2[iCent]->SetMarkerStyle(20);
			hMinePionV2[iCent]->SetMarkerColor(kBlue);
			hMineKaonV2[iCent] = (TH1D*) lsKaons[0][2]->FindObject(Form("fKaon_n22_gap10_cent%d_number0_0_px_desampled",iCent));
			hMineKaonV2[iCent]->SetMarkerStyle(20);
			hMineKaonV2[iCent]->SetMarkerColor(kBlue);
			hMineProtonV2[iCent] = (TH1D*) lsProtons[0][2]->FindObject(Form("fProton_n22_gap10_cent%d_number0_0_px_desampled",iCent));
			hMineProtonV2[iCent]->SetMarkerStyle(20);
			hMineProtonV2[iCent]->SetMarkerColor(kBlue);
			
			hMinePionV32[iCent] = (TH1D*) lsPions[1][2]->FindObject(Form("fPion_n32_gap10_cent%d_number0_0_px_desampled",iCent));
			hMinePionV32[iCent]->SetMarkerStyle(20);
			hMinePionV32[iCent]->SetMarkerColor(kBlue);
			hMineKaonV32[iCent] = (TH1D*) lsKaons[1][2]->FindObject(Form("fKaon_n32_gap10_cent%d_number0_0_px_desampled",iCent));
			hMineKaonV32[iCent]->SetMarkerStyle(20);
			hMineKaonV32[iCent]->SetMarkerColor(kBlue);
			hMineProtonV32[iCent] = (TH1D*) lsProtons[1][2]->FindObject(Form("fProton_n32_gap10_cent%d_number0_0_px_desampled",iCent));
			hMineProtonV32[iCent]->SetMarkerStyle(20);
			hMineProtonV32[iCent]->SetMarkerColor(kBlue);
		}
		
		printf("Loaded\n");

		TLegend* leg = new TLegend(0.12,0.6,0.3,0.89);
		leg->SetBorderSize(0);
		leg->AddEntry(hMinePionV2[2], "Mine QC2 #Delta#eta > 0","pel");
		leg->AddEntry(ProtonV22[2], "QC2 QM11 #Delta#eta > 1","pel");
		leg->AddEntry(ProtonV2[2], "SP #Delta#eta > 0.9","pel");
		
		TLatex latex;
		latex.SetNDC();


		// plotting
		TCanvas* cPID = new TCanvas("cPID");

		cPID->cd();
		
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			// v22
			if(!hMinePionV2[iCent])
			{
				printf(hMinePionV2[iCent]->GetName());
			  continue;
			}
			hMinePionV2[iCent]->SetMinimum(0.);
			hMinePionV2[iCent]->SetMaximum(0.5);
			hMinePionV2[iCent]->Draw();
			
			if(PionV2[iCent])
				PionV2[iCent]->Draw(sErrorOption.Data());

			if(PionV22[iCent])
				PionV22[iCent]->Draw(sErrorOption.Data());
			
			leg->Draw();
			latex.DrawLatex(0.65,0.8, Form("#pi: %d%% - %d%%",iCentBins[iCent],iCentBins[iCent+1]));
			cPID->SaveAs(Form("%s/Pions/PionV22_cent%d_comp.pdf",sOutPath.Data(),iCent));	
		}

	
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			if(!hMineKaonV2[iCent])
			{
				printf(hMineKaonV2[iCent]->GetName());
			  continue;
			}
			hMineKaonV2[iCent]->SetMinimum(0.);
			hMineKaonV2[iCent]->SetMaximum(0.5);
			hMineKaonV2[iCent]->Draw();
			
			if(KaonV2[iCent])
				KaonV2[iCent]->Draw(sErrorOption.Data());

			if(KaonV22[iCent])
				KaonV22[iCent]->Draw(sErrorOption.Data());
			
			leg->Draw();
			latex.DrawLatex(0.65,0.8, Form("K: %d%% - %d%%",iCentBins[iCent],iCentBins[iCent+1]));
			cPID->SaveAs(Form("%s/Kaons/KaonV22_cent%d_comp.pdf",sOutPath.Data(),iCent));
		}

		
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			if(!hMineProtonV2[iCent])
			{
				printf(hMineProtonV2[iCent]->GetName());
			  continue;
			}
			hMineProtonV2[iCent]->SetMinimum(0.);
			hMineProtonV2[iCent]->SetMaximum(0.5);
			hMineProtonV2[iCent]->Draw();

			if(ProtonV2[iCent])
				ProtonV2[iCent]->Draw(sErrorOption.Data());

			if(ProtonV22[iCent])
				ProtonV22[iCent]->Draw(sErrorOption.Data());

			leg->Draw();
			latex.DrawLatex(0.65,0.8, Form("p: %d%% - %d%%",iCentBins[iCent],iCentBins[iCent+1]));
			cPID->SaveAs(Form("%s/Protons/ProtonV22_cent%d_comp.pdf",sOutPath.Data(),iCent));
		}
 
		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			// v32
			if(!hMinePionV32[iCent])
			{
				printf(hMinePionV32[iCent]->GetName());
			  continue;
			}
			hMinePionV32[iCent]->SetMinimum(0.);
			hMinePionV32[iCent]->SetMaximum(0.5);
			hMinePionV32[iCent]->Draw();
			
			if(PionV3[iCent])
				PionV3[iCent]->Draw(sErrorOption.Data());

			if(PionV32[iCent])
				PionV32[iCent]->Draw(sErrorOption.Data());

			leg->Draw();
			latex.DrawLatex(0.65,0.8, Form("#pi: %d%% - %d%%",iCentBins[iCent],iCentBins[iCent+1]));
			cPID->SaveAs(Form("%s/Pions/PionV32_cent%d_comp.pdf",sOutPath.Data(),iCent));
		}

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			if(!hMineKaonV32[iCent])
			{
				printf(hMineKaonV32[iCent]->GetName());
			  continue;
			}
			hMineKaonV32[iCent]->SetMinimum(0.);
			hMineKaonV32[iCent]->SetMaximum(0.5);
			hMineKaonV32[iCent]->Draw();
			
			if(KaonV3[iCent])
				KaonV3[iCent]->Draw(sErrorOption.Data());

			if(KaonV32[iCent])
				KaonV32[iCent]->Draw(sErrorOption.Data());

			leg->Draw();
			latex.DrawLatex(0.65,0.8, Form("K: %d%% - %d%%",iCentBins[iCent],iCentBins[iCent+1]));
			cPID->SaveAs(Form("%s/Kaons/KaonV32_cent%d_comp.pdf",sOutPath.Data(),iCent));	
		}

		for(Short_t iCent(0); iCent < iNumCent; iCent++)
		{
			if(!hMineProtonV32[iCent])
			{
				printf(hMineProtonV3[iCent]->GetName());
			  continue;
			}
			hMineProtonV32[iCent]->SetMinimum(0.);
			hMineProtonV32[iCent]->SetMaximum(0.5);
			hMineProtonV32[iCent]->Draw();
			
			if(ProtonV3[iCent])
				ProtonV3[iCent]->Draw(sErrorOption.Data());

			if(ProtonV32[iCent])
				ProtonV32[iCent]->Draw(sErrorOption.Data());

			leg->Draw();
			latex.DrawLatex(0.65,0.8, Form("p: %d%% - %d%%",iCentBins[iCent],iCentBins[iCent+1]));
			cPID->SaveAs(Form("%s/Protons/ProtonV32_cent%d_comp.pdf",sOutPath.Data(),iCent));	
		}	
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



// TGraphErrors* Kv2_1020_QC2(Int_t color, Int_t marker) {
//   Int_t _nPoints = 19;
//   Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
//   Double_t _y[] = {0.007051, 0.016100, 0.037009, 0.057360, 0.076636, 0.096443, 0.110868, 0.124309, 0.134059, 0.143032, 0.148407, 0.156759, 0.161398, 0.156367, 0.148881, 0.142463, 0.135274, 0.108699, 0.109107};
//   Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
//   Double_t _yerr[] = {0.003485, 0.001355, 0.000849, 0.000704, 0.000690, 0.000741, 0.000824, 0.000949, 0.001078, 0.001346, 0.001635, 0.001515, 0.002217, 0.003181, 0.004494, 0.004505, 0.008140, 0.011124, 0.020627};
//   Double_t _ysys[] = {0.001480, 0.000588, 0.001159, 0.001734, 0.002305, 0.002927, 0.003347, 0.003831, 0.004036, 0.004330, 0.004453, 0.004743, 0.005010, 0.005319, 0.004618, 0.004824, 0.004091, 0.003268, 0.013869};
//   //if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
//   //if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
//   //if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
//   TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
//   graph->SetLineColor(color);
//   graph->SetMarkerColor(color);
//   graph->SetMarkerStyle(marker);
//   return graph;
// }

// TGraphAsymmErrors* v2Pion1020(Int_t color, Int_t marker, Int_t first,Int_t last)
// {
//   //commentme
//   Int_t _nPoints = 45;
//   if (last>_nPoints-1) last=_nPoints-1;
//   if (last<0 && first<0) last=_nPoints-1;
//   if (last<0) last=_nPoints-1+last;
//   if (first<0) first=0;
//   Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
//   Double_t _y[] = {0, 0, 0, 0, 0.0201224, 0.0252062, 0.0311323, 0.0369631, 0.0428155, 0.0483589, 0.0550696, 0.0651319, 0.0746561, 0.0834061, 0.0919292, 0.0993882, 0.106362, 0.11341, 0.119441, 0.125327, 0.130489, 0.135768, 0.140111, 0.1436, 0.147624, 0.15197, 0.157217, 0.15934, 0.159981, 0.157703, 0.15448, 0.151198, 0.142915, 0.132599, 0.113408, 0.10271, 0.0943079, 0.0879974, 0.0790374, 0.0795734, 0.0614621, 0.067106, 0.0609654, 0.0515049, 0.0223675};
//   Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//   Double_t _yerr[] = {0, 0, 0, 0, 7.16967e-05, 7.04768e-05, 7.14158e-05, 7.34838e-05, 7.64088e-05, 8.00149e-05, 8.9464e-05, 9.55187e-05, 0.000105152, 0.000116857, 0.000130562, 0.000146276, 0.000163809, 0.000183749, 0.00020652, 0.000232603, 0.000262503, 0.000296109, 0.000334244, 0.000377637, 0.00042637, 0.00033071, 0.000441212, 0.000580034, 0.000749139, 0.000958424, 0.00120529, 0.00149142, 0.00182397, 0.0016818, 0.00233162, 0.00308243, 0.00393626, 0.00487286, 0.0058437, 0.00539979, 0.00725315, 0.00952612, 0.0093786, 0.0123107, 0.0176851};
//   Double_t _yerr2[] = {0, 0, 0, 0, 7.16967e-05, 7.04768e-05, 7.14158e-05, 7.34838e-05, 7.64088e-05, 8.00149e-05, 8.9464e-05, 9.55187e-05, 0.000105152, 0.000116857, 0.000130562, 0.000146276, 0.000163809, 0.000183749, 0.00020652, 0.000232603, 0.000262503, 0.000296109, 0.000334244, 0.000377637, 0.00042637, 0.00033071, 0.000441212, 0.000580034, 0.000749139, 0.000958424, 0.00120529, 0.00149142, 0.00182397, 0.0016818, 0.00233162, 0.00308243, 0.00393626, 0.00487286, 0.0058437, 0.00539979, 0.00725315, 0.00952612, 0.0093786, 0.0123107, 0.0176851};

//   /*
//   if(!kStat){
//     for(Int_t i=0;i<_nPoints;i++){
//       _yerr[i] = 0;
//       _yerr2[i] = 0;
//       _xerr[i] = 0.05;
//     }
//   }


//   if(kSyst){
//     for(Int_t i=0;i<_nPoints;i++){
// 	Float_t pol0 = 7.36284000000000027e-04;
// 	Float_t pol1 = 8.43209000000000037e-04;
	
// 	Float_t nonflow = (pol0 + pol1*_x[i]);
// 	Float_t systerr = _y[i]*systPi(2, _x[i]);
// 	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
// 	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
//     }
//   }
//   */

//   TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
//   graph->SetLineColor(color);
//   graph->SetMarkerColor(color);
//   graph->SetMarkerStyle(marker);
//   return graph;
// }

// TGraphAsymmErrors* v2Kaon1020(Int_t color, Int_t marker, Int_t first,Int_t last)
// {
//   //commentme
//   Int_t _nPoints = 45;
//   if (last>_nPoints-1) last=_nPoints-1;
//   if (last<0 && first<0) last=_nPoints-1;
//   if (last<0) last=_nPoints-1+last;
//   if (first<0) first=0;
//   Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
//   Double_t _y[] = {0, 0, 0, 0, 0, 0.00640075, 0.00890777, 0.0114231, 0.0121152, 0.0193763, 0.0262718, 0.0373136, 0.04773, 0.0585364, 0.0685442, 0.0781197, 0.0884599, 0.0959829, 0.104197, 0.110977, 0.118551, 0.124407, 0.130268, 0.137113, 0.141952, 0.149588, 0.158874, 0.165006, 0.168625, 0.168325, 0.16697, 0.173543, 0.161645, 0.15523, 0.14532, 0.128752, 0.142204, 0.0843617, 0.269387, 0.0231167, 0.431226, 0.0167267, 0.315711, 0.235569, -0.000340595};
//   Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//   Double_t _yerr[] = {0, 0, 0, 0, 0, 0.000525284, 0.000429881, 0.000408989, 0.000491132, 0.000713549, 0.000405664, 0.000346767, 0.000322585, 0.000314851, 0.000318631, 0.00032842, 0.000344056, 0.000364314, 0.00039006, 0.000423204, 0.000459161, 0.000504151, 0.000555309, 0.000614542, 0.000683822, 0.000519796, 0.000688975, 0.00091857, 0.00122243, 0.00160751, 0.00211499, 0.00275193, 0.0035775, 0.00370145, 0.00631253, 0.0106874, 0.0182738, 0.0313279, 0.0476483, 0.0611653, 0.106198, 0.207285, 0.182173, 0.191574, 0.224581};
//   Double_t _yerr2[] = {0, 0, 0, 0, 0, 0.000525284, 0.000429881, 0.000408989, 0.000491132, 0.000713549, 0.000405664, 0.000346767, 0.000322585, 0.000314851, 0.000318631, 0.00032842, 0.000344056, 0.000364314, 0.00039006, 0.000423204, 0.000459161, 0.000504151, 0.000555309, 0.000614542, 0.000683822, 0.000519796, 0.000688975, 0.00091857, 0.00122243, 0.00160751, 0.00211499, 0.00275193, 0.0035775, 0.00370145, 0.00631253, 0.0106874, 0.0182738, 0.0313279, 0.0476483, 0.0611653, 0.106198, 0.207285, 0.182173, 0.191574, 0.224581};

//   /*
//   if(!kStat){
//     for(Int_t i=0;i<_nPoints;i++){
//       _yerr[i] = 0;
//       _yerr2[i] = 0;
//       _xerr[i] = 0.05;
//     }
//   }


//   if(kSyst){
//     for(Int_t i=0;i<_nPoints;i++){
// 	Float_t pol0 = 7.36284000000000027e-04;
// 	Float_t pol1 = 8.43209000000000037e-04;

// 	Float_t nonflow = (pol0 + pol1*_x[i]);
// 	Float_t systerr = _y[i]*systKa(2, _x[i]);
// 	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
// 	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
//     }
//   }
// 	*/

//   TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
//   graph->SetLineColor(color);
//   graph->SetMarkerColor(color);
//   graph->SetMarkerStyle(marker);
//   return graph;
// }


// TGraphAsymmErrors* v2Antiproton1020(Int_t color, Int_t marker, Int_t first,Int_t last)
// {
//   //commentme
//   Int_t _nPoints = 45;
//   if (last>_nPoints-1) last=_nPoints-1;
//   if (last<0 && first<0) last=_nPoints-1;
//   if (last<0) last=_nPoints-1+last;
//   if (first<0) first=0;
//   Double_t _x[] = {0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 11, 13.5, 17.5};
//   Double_t _y[] = {0, 0, 0, 0, 0, 0, 0.00514464, 0.00344225, 0.00468868, 0.00546, 0.0082079, 0.012222, 0.0165934, 0.0233094, 0.0288221, 0.0374016, 0.0454486, 0.0563712, 0.064759, 0.075207, 0.0867709, 0.0963642, 0.106439, 0.116716, 0.12656, 0.138558, 0.157754, 0.173931, 0.185026, 0.199048, 0.206145, 0.208087, 0.213577, 0.213699, 0.21574, 0.192142, 0.152677, 0.172787, 0.122158, 0.11594, 0.170479, 0.143031, 0.0113144, 0.188084, 0.909081};
//   Double_t _xerr[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//   Double_t _yerr[] = {0, 0, 0, 0, 0, 0, 0.00163153, 0.00131053, 0.00108878, 0.000950215, 0.000582294, 0.000516432, 0.000703097, 0.000643592, 0.000611608, 0.000594501, 0.000588832, 0.000590229, 0.000600458, 0.000619565, 0.00064412, 0.00067877, 0.000720346, 0.000770869, 0.000828332, 0.000600635, 0.000741676, 0.000928701, 0.00116982, 0.00102807, 0.00131593, 0.00169927, 0.00218807, 0.00222766, 0.00363962, 0.00585315, 0.00929629, 0.013941, 0.020436, 0.0231194, 0.046922, 0.0895645, 0.111115, 0.118019, 0.100696};
//   Double_t _yerr2[] = {0, 0, 0, 0, 0, 0, 0.00163153, 0.00131053, 0.00108878, 0.000950215, 0.000582294, 0.000516432, 0.000703097, 0.000643592, 0.000611608, 0.000594501, 0.000588832, 0.000590229, 0.000600458, 0.000619565, 0.00064412, 0.00067877, 0.000720346, 0.000770869, 0.000828332, 0.000600635, 0.000741676, 0.000928701, 0.00116982, 0.00102807, 0.00131593, 0.00169927, 0.00218807, 0.00222766, 0.00363962, 0.00585315, 0.00929629, 0.013941, 0.020436, 0.0231194, 0.046922, 0.0895645, 0.111115, 0.118019, 0.100696};

//   /*
//   if(!kStat){
//     for(Int_t i=0;i<_nPoints;i++){
//       _yerr[i] = 0;
//       _yerr2[i] = 0;
//       _xerr[i] = 0.05;
//     }
//   }


//   if(kSyst){
//     for(Int_t i=0;i<_nPoints;i++){
// 	Float_t pol0 = 2.13364693747489361e-03;
// 	Float_t pol1 = 1.52812297971200361e-03;

// 	Float_t nonflow = (pol0 + pol1*_x[i]);
// 	Float_t systerr = _y[i]*systPr(2, _x[i]);
// 	_yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + systerr*systerr);
// 	_yerr2[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + nonflow*nonflow);
//     }
//   }
// 	*/
//   TGraphAsymmErrors* graph = new TGraphAsymmErrors(last-first+1, &_x[first], &_y[first], &_xerr[first],&_xerr[first], &_yerr2[first],&_yerr[first]);
//   graph->SetLineColor(color);
//   graph->SetMarkerColor(color);
//   graph->SetMarkerStyle(marker);
//   return graph;
// }
// TGraphErrors* Kv2_1020_QC4(Int_t color, Int_t marker) {
//   Int_t _nPoints = 19;
//   Double_t _x[] = {0.300000, 0.500000, 0.700000, 0.900000, 1.100000, 1.300000, 1.500000, 1.700000, 1.900000, 2.100000, 2.300000, 2.600000, 3.000000, 3.400000, 3.800000, 4.500000, 5.500000, 7.000000, 10.000000};
//   Double_t _y[] = {0.010304, 0.016204, 0.033282, 0.050568, 0.068130, 0.082275, 0.098983, 0.107410, 0.118314, 0.122915, 0.127885, 0.135673, 0.133027, 0.125030, 0.119996, 0.114883, 0.103795, 0.056606, 0.073532};
//   Double_t _xerr[] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
//   Double_t _yerr[] = {0.008246, 0.003348, 0.002248, 0.001928, 0.001887, 0.001998, 0.002171, 0.002489, 0.002762, 0.003426, 0.004034, 0.003694, 0.005368, 0.007697, 0.010785, 0.010791, 0.019094, 0.026001, 0.046208};
//   Double_t _ysys[] = {0.003572, 0.000550, 0.001082, 0.001520, 0.002068, 0.003303, 0.003411, 0.003818, 0.003650, 0.003819, 0.003868, 0.004158, 0.005506, 0.003772, 0.003747, 0.003794, 0.003375, 0.008648, 0.025117};
//   //if(!kStat) for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = 0;
//   //if(!kStat) for(Int_t i=0;i!=_nPoints;++i){ _yerr[i] = 0; _xerr[i] = 0.05;}
//   //if(kSyst)  for(Int_t i=0;i!=_nPoints;++i) _yerr[i] = TMath::Sqrt(_yerr[i]*_yerr[i] + _ysys[i]*_ysys[i]);
//   TGraphErrors* graph = new TGraphErrors(_nPoints, _x, _y, _xerr, _yerr);
//   graph->SetLineColor(color);
//   graph->SetMarkerColor(color);
//   graph->SetMarkerStyle(marker);
//   return graph;
// }


