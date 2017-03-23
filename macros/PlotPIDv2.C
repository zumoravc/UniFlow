void PlotPIDv2()
{
	TFile* fInput = new TFile("/Users/vpacik/NBI/Flow/results/v2-full-ver2-pPbRun2-lhc16q/Flow_CENT_wSDD_ver2-2.root","READ");
	//TFile* fInput = new TFile("/Users/vpacik/NBI/Flow/results/v2-fullPID-pp-lhc16-kl-all/Flow_kINT7.root","READ");
	if(!fInput)
	{
		printf("File not found!\n");
		return;
	}


	TString outPath = "/Users/vpacik/NBI/Flow/results/v2-full-ver2-pPbRun2-lhc16q/plots3";
	TString sOutputFormat = "eps";

	Double_t dCentEdges[] = {0.,50.,100};
	const Short_t iCent = 2;

	TString sGaps[] = {"-10","00","08"};
	TString sGapsNice[] = {"No Gap","Gap > 0","Gap > 0.8"};
	const Short_t iGaps = 3;

	Int_t GapStyle[iGaps+1] = {20,21,22};
	Int_t PIDStyle[] = {22,29,33,20,34};
	Double_t MarkerSize[] = {1.2,1.5,1.5,1.2,1.5};
	// Color_t PIDCol[] = {kRed,kOrange+4,kMagenta,kGreen-2,kBlue};
	Color_t PIDCol[] = {kCyan+2,kAzure+1,kMagenta-4,kBlue,kRed};


	TH1D* hRef[iGaps];
	TH1D* hTracks[iGaps][iCent];
	TH1D* hPions[iGaps][iCent];
	TH1D* hKaons[iGaps][iCent];
	TH1D* hProtons[iGaps][iCent];
	TH1D* hK0s[iGaps][iCent];
	TH1D* hLambda[iGaps][iCent];
	TH1D* hRef4;
	TH1D* hTracks4[iCent];
	TH1D* hPions4[iCent];
	TH1D* hKaons4[iCent];
	TH1D* hProtons4[iCent];
	TList* listTracks = 0x0;
	TList* listPions = 0x0;
	TList* listKaons = 0x0;
	TList* listProtons = 0x0;
	TList* listTracks4 = 0x0;
	TList* listPions4 = 0x0;
	TList* listKaons4 = 0x0;
	TList* listProtons4 = 0x0;

	fInput->ls();

	gStyle->SetErrorX(0);

	fInput->cd();
	for(Short_t gap(0); gap < iGaps; gap++)
	{
		hRef[gap] = (TH1D*) gDirectory->Get(Form("hRef_n22_Gap%s",sGaps[gap].Data()));
		// if(gap == 0)
		// 	hRef4 = (TH1D*) gDirectory->Get(Form("hRef_n24_Gap%s",sGaps[gap].Data()));

		listTracks = (TList*) gDirectory->Get(Form("Tracks_n2_Gap%s",sGaps[gap].Data()));
		listPions = (TList*) gDirectory->Get(Form("Pions_n2_Gap%s",sGaps[gap].Data()));
		listKaons = (TList*) gDirectory->Get(Form("Kaons_n2_Gap%s",sGaps[gap].Data()));
		listProtons = (TList*) gDirectory->Get(Form("Protons_n2_Gap%s",sGaps[gap].Data()));

		// if(gap == 0)
		// {
		// 	listTracks4 = (TList*) gDirectory->Get(Form("Tracks_n24_Gap%s",sGaps[gap].Data()));
		// 	listPions4 = (TList*) gDirectory->Get(Form("Pions_n24_Gap%s",sGaps[gap].Data()));
		// 	listKaons4 = (TList*) gDirectory->Get(Form("Kaons_n24_Gap%s",sGaps[gap].Data()));
		// 	listProtons4 = (TList*) gDirectory->Get(Form("Protons_n24_Gap%s",sGaps[gap].Data()));
		// }

		for(Short_t cent(0); cent < iCent; cent++)
		{
			hTracks[gap][cent] = (TH1D*) listTracks->FindObject(Form("fTracks_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			hPions[gap][cent] = (TH1D*) listPions->FindObject(Form("fPion_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			hKaons[gap][cent] = (TH1D*) listKaons->FindObject(Form("fKaon_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			hProtons[gap][cent] = (TH1D*) listProtons->FindObject(Form("fProton_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			hK0s[gap][cent] = (TH1D*) fInput->Get(Form("hK0s_n22_gap%s_cent%d",sGaps[gap].Data(),cent));
			hLambda[gap][cent] = (TH1D*) fInput->Get(Form("hLambda_n22_gap%s_cent%d",sGaps[gap].Data(),cent));

			// if(gap == 0)
			// {
			// 	hTracks4[cent] = (TH1D*) listTracks->FindObject(Form("fTracks_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			// 	hPions4[cent] = (TH1D*) listPions->FindObject(Form("fPion_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			// 	hKaons4[cent] = (TH1D*) listKaons->FindObject(Form("fKaon_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			// 	hProtons4[cent] = (TH1D*) listProtons->FindObject(Form("fProton_n22_gap%s_cent%d_number0_0_px_desampled",sGaps[gap].Data(),cent));
			// }
		}
	}

	// // loading katarinas cum's
	// TFile* fKatar = new TFile("/Users/vpacik/NBI/Flow/results/v2-PID-pp-lhc16-kl/pp13TeVLHC16kl_pt0230_FB768_default.root","READ");
	// if(!fKatar)
	// {
	// 	printf("File not found!\n");
	// 	return;
	// }
	//
	// fKatar->ls();
	// TH1D* hKatarC22 = (TH1D*) fKatar->Get("Hcum22Ntrks");
	// Double_t rebin[] = {0,50,150};
	// const Int_t irebin = 2;
	// Short_t leng[irebin] = {50,100};
	//
	// TH1D* hKatarC22_rebin = (TH1D*) hKatarC22->Rebin(irebin,"hKatarC22_rebin",rebin);
	// hKatarC22_rebin->SetMarkerStyle(24);
	//
	// Double_t dCon;
	// for(Short_t i(1); i < irebin+1; i++)
	// {
	// 	dCon = hKatarC22_rebin->GetBinContent(i);
	// 	hKatarC22_rebin->SetBinContent(i,TMath::Sqrt(dCon/leng[i-1]));
	// }
	//

	// setting errors
	const Int_t iNumBinsRef = iCent;
	for(Short_t gap(0); gap < iGaps; gap++)
		for(Int_t bin(1); bin < iNumBinsRef+1; bin++)
		{
			// hRef[gap]->SetBinError(bin,0.001+(bin-1)*0.02);
			//hRef4->SetBinError(bin,0.001+(bin-1)*0.03);
		}

	const Int_t iNumBinsDiff = hTracks[0][0]->GetNbinsX();
	for(Short_t cent(0); cent < iCent; cent++)
		for(Int_t bin(1); bin < iNumBinsDiff+1; bin++)
		{
			// hTracks[0][cent]->SetBinError(bin,bin*0.0001+cent*0.002);
			// hPions[0][cent]->SetBinError(bin,bin*0.0001+cent*0.0005);
			// hKaons[0][cent]->SetBinError(bin,bin*0.0001+cent*0.0005);
			// hProtons[0][cent]->SetBinError(bin,bin*0.0001+cent*0.0005);
			//
			// hTracks4[cent]->SetBinError(bin,bin*0.002+cent*0.002);
			// hPions4[cent]->SetBinError(bin,bin*0.004+cent*0.004);
			// hKaons4[cent]->SetBinError(bin,bin*0.004+cent*0.004);
			// hProtons4[cent]->SetBinError(bin,bin*0.004+cent*0.004);
		}


	// TCanvas* ctemp = new TCanvas("ctem","temp");
	// ctemp->cd();
	// hKatarC22->Draw();
	// hKatarC22_rebin->Draw("same");

	gStyle->SetLegendBorderSize(0);
	TLegend* legGap = new TLegend(0.15,0.7,0.4,0.8);


	TLegend* leg = new TLegend(0.15,0.7,0.5,0.88);
	// leg->SetBorderSize(1.);
	//leg->AddEntry(hTracks[0][0],"charged","pel");
	leg->AddEntry(hPions[0][0],"pions","p");
	leg->AddEntry(hKaons[0][0],"kaons","p");
	leg->AddEntry(hK0s[0][0],"K0s","p");
	leg->AddEntry(hProtons[0][0],"protons","p");
	leg->AddEntry(hLambda[0][0],"Lambda","p");


	legGap->AddEntry(hRef[0],Form("No Gap","p"));
	legGap->AddEntry(hRef[1],Form("Gap 0","p"));
	legGap->AddEntry(hRef[2],Form("Gap 0.8","p"));
	//legGap->AddEntry(hRef4,Form("v2{4}","pel"));

	TCanvas* canRef = new TCanvas("canRef","Ref");
	for(Short_t gap(0); gap < iGaps; gap++)
	{
		canRef->cd();
		hRef[gap]->SetMarkerStyle(GapStyle[gap]);
		hRef[gap]->SetMarkerColor(kBlack);
		hRef[gap]->SetLineColor(kBlack);
		hRef[gap]->SetMinimum(-0.05);
		hRef[gap]->SetMaximum(0.2);
		hRef[gap]->Draw("ep same");
	}
	// //hKatarC22_rebin->Draw("same");
	// hRef4->SetLineColor(kBlack);
	// hRef4->SetMarkerColor(kBlack);
	// hRef4->SetMarkerStyle(GapStyle[3]);
	//hRef4->Draw("ep same");
	legGap->Draw("same");
	canRef->SaveAs(Form("%s/v22-Ref.%s",outPath.Data(),sOutputFormat.Data()));

	TCanvas* canCharged = new TCanvas("canCharged","Charged");
	for(Short_t cent(0); cent < iCent; cent++)
	{
		for(Short_t gap(0); gap < iGaps; gap++)
		{
			//legGap->AddEntry(hRef[gap],Form("Gap %s",sGaps[gap].Data()),"pel");
			canCharged->cd(1);
			hTracks[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hTracks[gap][cent]->SetMarkerColor(kBlack);
			hTracks[gap][cent]->SetLineColor(kBlack);
			hTracks[gap][cent]->SetMinimum(0);
			hTracks[gap][cent]->SetMaximum(0.6);
			if(gap == 0)
				hTracks[gap][cent]->Draw("ep");
			else
				hTracks[gap][cent]->Draw("ep same");
		}
		//hKatarC22_rebin->Draw("same");
		//hTracks4[cent]->SetMarkerStyle(GapStyle[3]);
		//hTracks4[cent]->SetMarkerColor(kBlack);
		//hTracks4[cent]->SetLineColor(kBlack);
		//hTracks4[cent]->Draw("ep same");
		legGap->Draw("same");
		canCharged->SaveAs(Form("%s/v22-charged-cent%d.%s",outPath.Data(),cent,sOutputFormat.Data()));
	}

	TCanvas* canPID = new TCanvas("canPID","PID",600,600);
	//canPID->Divide(3);


	for(Short_t cent(0); cent < iCent; cent++)
	{
		for(Short_t gap(0); gap < iGaps; gap++)
		{
			canPID->cd(1);
			leg->SetHeader(Form("%s cent %d",sGapsNice[gap].Data(),cent));
			hPions[gap][cent]->SetTitle(";#it{p}_{T} (GeV/c);");
			// hPions[gap][cent]->SetLabelOffset(0.01,"Y");
			hPions[gap][cent]->SetMinimum(0.);
			hPions[gap][cent]->SetMaximum(0.45);
			hPions[gap][cent]->SetLineColor(PIDCol[0]);
			hPions[gap][cent]->SetMarkerColor(PIDCol[0]);
			// hPions[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hPions[gap][cent]->SetMarkerStyle(PIDStyle[0]);
			hPions[gap][cent]->SetMarkerSize(MarkerSize[0]);

			hKaons[gap][cent]->SetLineColor(PIDCol[1]);
			// hKaons[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hKaons[gap][cent]->SetMarkerStyle(PIDStyle[1]);
			hKaons[gap][cent]->SetMarkerColor(PIDCol[1]);
			hKaons[gap][cent]->SetMarkerSize(MarkerSize[1]);

			hProtons[gap][cent]->SetLineColor(PIDCol[2]);
			// hProtons[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hProtons[gap][cent]->SetMarkerStyle(PIDStyle[2]);
			hProtons[gap][cent]->SetMarkerColor(PIDCol[2]);
			hProtons[gap][cent]->SetMarkerSize(MarkerSize[2]);

			hK0s[gap][cent]->SetLineColor(PIDCol[3]);
			// hK0s[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hK0s[gap][cent]->SetMarkerStyle(PIDStyle[3]);
			hK0s[gap][cent]->SetMarkerColor(PIDCol[3]);
			hK0s[gap][cent]->SetMarkerSize(MarkerSize[3]);
			//hLambda[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hLambda[gap][cent]->SetLineColor(PIDCol[4]);
			hLambda[gap][cent]->SetMarkerStyle(PIDStyle[4]);
			hLambda[gap][cent]->SetMarkerColor(PIDCol[4]);
			hLambda[gap][cent]->SetMarkerSize(MarkerSize[4]);

			hTracks[gap][cent]->SetLineColor(PIDCol[3]);
			hTracks[gap][cent]->SetMarkerStyle(GapStyle[gap]);
			hTracks[gap][cent]->SetMarkerColor(PIDCol[3]);

			hPions[gap][cent]->Draw("p");
			hKaons[gap][cent]->Draw("p same");
			hProtons[gap][cent]->Draw("p same");
			hK0s[gap][cent]->Draw("pe0 same");
			hLambda[gap][cent]->Draw("p same");
			//hTracks[gap][cent]->Draw("ep same");
			leg->Draw("same");
			canPID->SaveAs(Form("%s/v22-cent%d-Gap%s.%s",outPath.Data(),cent,sGaps[gap].Data(),sOutputFormat.Data()));
		}

		// canPID->cd(1);
		// leg->SetHeader(Form("gap %s","v2{4}"));
		// hPions4[cent]->SetMinimum(0.);
		// hPions4[cent]->SetMaximum(0.6);
		// hPions4[cent]->SetLineColor(PIDCol[0]);
		// hPions4[cent]->SetMarkerColor(PIDCol[0]);
		// hPions4[cent]->SetMarkerStyle(GapStyle[3]);
		// hKaons4[cent]->SetLineColor(PIDCol[1]);
		// hKaons4[cent]->SetMarkerStyle(GapStyle[3]);
		// hKaons4[cent]->SetMarkerColor(PIDCol[1]);
		// hProtons4[cent]->SetLineColor(PIDCol[2]);
		// hProtons4[cent]->SetMarkerStyle(GapStyle[3]);
		// hProtons4[cent]->SetMarkerColor(PIDCol[2]);
		// hTracks4[cent]->SetLineColor(PIDCol[3]);
		// hTracks4[cent]->SetMarkerStyle(GapStyle[3]);
		// hTracks4[cent]->SetMarkerColor(PIDCol[3]);
		//
		// hPions4[cent]->Draw("ep");
		// hKaons4[cent]->Draw("ep same");
		// hProtons4[cent]->Draw("ep same");
		// // //hTracks[gap][cent]->Draw("ep same");
		// leg->Draw("same");
		// canPID->SaveAs(Form("%s/v24-cent%d-Gap%s.pdf",outPath.Data(),cent,"NoGap"));
	}




}
