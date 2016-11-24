TCanvas* PlotCent(const Short_t iCent, TString sTag);

void PlotFullPID()
{
	TString sOutput = TString("~/NBI/Flow/temp");

	TFile* fInput = new TFile("~/NBI/Flow/temp/Flow-10samples.root","READ");
	fInput->cd();
	fInput->ls();


	TList* lTracks_n22 = (TList*) fInput->Get("Tracks_n2_Gap-10");
	TList* lDiffv22_Pions = (TList*) fInput->Get("n2_Pions");
	TList* lDiffv24_Pions = (TList*) fInput->Get("Diffv24_Pions");
	TList* lDiffv22_Kaons = (TList*) fInput->Get("Diffv22_Kaons");
	TList* lDiffv24_Kaons = (TList*) fInput->Get("Diffv24_Kaons");
	TList* lDiffv22_Protons = (TList*) fInput->Get("Diffv22_Protons");
	TList* lDiffv24_Protons = (TList*) fInput->Get("Diffv24_Protons");
	return;

	//TH1D* hRefv22_Tracks = (TH1D*) fInput->Get("Refv22_Tracks");
	//TH1D* hRefv24_Tracks = (TH1D*) fInput->Get("Refv24_Tracks");

	TFile* fInputRef = new TFile("~/NBI/Flow/temp/ProcessFlow_FB768_GFK.root","READ");
	TH1D* hRefv22 = (TH1D*) fInputRef->Get("Refv22_Tracks");
	TH1D* hRefv22Gap00 = (TH1D*) fInputRef->Get("Refv22_Gap00_Tracks");
	TH1D* hRefv22Gap10 = (TH1D*) fInputRef->Get("Refv22_Gap10_Tracks");
	TH1D* hRefv32 = (TH1D*) fInputRef->Get("Refv32_Tracks");
	TH1D* hRefv42 = (TH1D*) fInputRef->Get("Refv42_Tracks");
	TH1D* hRefv24 = (TH1D*) fInputRef->Get("Refv24_Tracks");

	TFile* fInputV0s = new TFile("~/NBI/Flow/results/V0s/15-sampling-5bins/sampling/PtFlow_V0s_Gap00.root","READ");


	Color_t colorPID[] = {kRed, kGreen, kBlue, kBlack, kMagenta, kGreen+2, kRed+2, kBlue+2,kMagenta+2};
	TString sLegend[] = {"#pi (no gap)","K (no gap)" ,"p (no gap)","#Lambda (v22 Gap 0)","K0s (v22 Gap 0)"};
	
	TH1D* histPion = 0x0;
	TH1D* histKaon = 0x0;
	TH1D* histProton = 0x0;
	TH1D* histK0s = 0x0;
	TH1D* histLambda = 0x0;

	const Short_t iNumCent = 9;

	TCanvas* cCanPions = new TCanvas();
	TLegend* legPions = new TLegend(0.15,0.6,0.26,0.9);
	legPions->SetBorderSize(0);
	TCanvas* cCan = new TCanvas();

	// reference
	hRefv22->SetStats(0);
	hRefv22->SetMinimum(0.);
	hRefv22->SetMaximum(0.2);
	hRefv22->SetLineColor(colorPID[0]);
	hRefv22Gap00->SetLineColor(colorPID[1]);
	hRefv22Gap10->SetLineColor(colorPID[2]);
	hRefv32->SetLineColor(colorPID[3]);
	hRefv42->SetLineColor(colorPID[4]);
	hRefv24->SetLineColor(colorPID[5]);
	hRefv22->Draw();
	hRefv22Gap00->Draw("same");
	hRefv22Gap10->Draw("same");
	hRefv32->Draw("same");
	hRefv42->Draw("same");
	hRefv24->Draw("same");
	TLegend* legRef = new TLegend(0.15,0.6,0.35,0.9);
	legRef->SetBorderSize(0);
	legRef->AddEntry(hRefv22,"v22 NoGap","pel");
	legRef->AddEntry(hRefv22Gap00,"v22 Gap00","pel");
	legRef->AddEntry(hRefv22Gap10,"v22 Gap10","pel");
	legRef->AddEntry(hRefv24,"v24 NoGap","pel");
	legRef->AddEntry(hRefv32,"v32 NoGap","pel");
	legRef->AddEntry(hRefv42,"v42 NoGap","pel");
	legRef->Draw("same");
	cCan->Print(Form("%s/Ref_Tracks.png",sOutput.Data()));

	// PID v22
	for(Int_t i(0); i < iNumCent; i++)
	{
		TLegend* leg = new TLegend(0.15,0.6,0.35,0.9);
		leg->SetBorderSize(0);
		
		cCan->cd();

		//v22 
		histPion = (TH1D*) lDiffv22_Pions->At(i);
		histKaon = (TH1D*) lDiffv22_Kaons->At(i);
		histProton = (TH1D*) lDiffv22_Protons->At(i);
		histLambda = (TH1D*) fInputV0s->Get(Form("hFlowPt_Lambda_Gap00_Cent%d",i));
		histK0s = (TH1D*) fInputV0s->Get(Form("hFlowPt_K0s_Gap00_Cent%d",i));
		
		histProton->SetStats(0);
		histProton->SetMinimum(-0.05);
		histProton->SetMaximum(0.45);
		histPion->SetLineColor(colorPID[0]);
		histPion->SetMarkerColor(colorPID[0]);
		histKaon->SetLineColor(colorPID[1]);
		histKaon->SetMarkerColor(colorPID[1]);
		histProton->SetLineColor(colorPID[2]);
		histProton->SetMarkerColor(colorPID[2]);
		histLambda->SetLineColor(colorPID[3]);
		histLambda->SetMarkerColor(colorPID[3]);
		histK0s->SetLineColor(colorPID[4]);
		histK0s->SetMarkerColor(colorPID[4]);

		leg->AddEntry(histPion,sLegend[0].Data(),"pel");
		leg->AddEntry(histKaon,sLegend[1].Data(),"pel");
		leg->AddEntry(histProton,sLegend[2].Data(),"pel");
		leg->AddEntry(histLambda,sLegend[3].Data(),"pel");
		leg->AddEntry(histK0s,sLegend[4].Data(),"pel");

		histProton->Draw();
		histPion->Draw("same");
		histKaon->Draw("same");
		histLambda->Draw("same");
		histK0s->Draw("same");
		leg->Draw("same");
		cCan->Print(Form("%s/PID_v22_cent%d.png",sOutput.Data(),i));


		// v24
		histPion = (TH1D*) lDiffv24_Pions->At(i);
		histKaon = (TH1D*) lDiffv24_Kaons->At(i);
		histProton = (TH1D*) lDiffv24_Protons->At(i);
		
		histProton->SetStats(0);
		histProton->SetMinimum(-0.05);
		histProton->SetMaximum(0.45);
		histPion->SetLineColor(colorPID[0]);
		histPion->SetMarkerColor(colorPID[0]);
		histKaon->SetLineColor(colorPID[1]);
		histKaon->SetMarkerColor(colorPID[1]);
		histProton->SetLineColor(colorPID[2]);
		histProton->SetMarkerColor(colorPID[2]);

		histProton->Draw();
		histPion->Draw("same");
		histKaon->Draw("same");
		histLambda->Draw("same");
		histK0s->Draw("same");
		leg->Draw("same");
		cCan->Print(Form("%s/PID_v24_cent%d.png",sOutput.Data(),i));

		// v22 cent for pions
		
		cCanPions->cd();
		histPion = (TH1D*) lDiffv24_Pions->At(i);
		histPion->SetMinimum(-0.05);
		histPion->SetMaximum(0.3);
		if(i == 8 || i == 7 || i == 0)
			continue; 
		histPion->SetLineColor(colorPID[i]);
		histPion->SetMarkerColor(colorPID[i]);
		legPions->AddEntry(histPion,Form("Pion cent%d",i),"pel");
		histPion->Draw("same");


	}
	legPions->Draw("same");
	cCanPions->Print(Form("%s/PID_v24_Pions.png",sOutput.Data()));



	return;
}

TCanvas* PlotCent(const Short_t iCent, TString sTag)
{
	TH1D* hPion22 = (TH1D*) gDirectory->Get(Form("hPions%s_cent%d",sTag.Data(),iCent));
	TH1D* hKaon22 = (TH1D*) gDirectory->Get(Form("hKaons%s_cent%d",sTag.Data(),iCent));
	TH1D* hProton22 = (TH1D*) gDirectory->Get(Form("hProtons%s_cent%d",sTag.Data(),iCent));
	

	Color_t colors[3] = {kRed,kBlue,kGreen+2};
	TCanvas* cCent = new TCanvas();
	cCent->cd();
	hPion22->SetMinimum(-0.05);
	hPion22->SetMaximum(0.3);
	hPion22->SetLineColor(colors[0]);
	hPion22->SetMarkerColor(colors[0]);
	hKaon22->SetLineColor(colors[1]);
	hKaon22->SetMarkerColor(colors[1]);
	hProton22->SetLineColor(colors[2]);
	hProton22->SetMarkerColor(colors[2]);
	
	hPion22->Draw();
	hKaon22->Draw("same");
	hProton22->Draw("same");

	return cCent;
}
