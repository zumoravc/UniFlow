TH1D* RedoFlow(TH1D* fInput, const char* name);
TH1D* MakeRatio(TH1D* fNom, TH1D* fDenom, const char* name);
//TCanvas* CompareRatio(TH1D* fNom, TH1D* fDenom);


void ProcessFlow(
			TString sInput = "~/NBI/Codes/flow/results/V0s/9/merge/AnalysisResults.root",
			TString sOutput = "~/NBI/Codes/results/TPConly",
			TString sOutputFormat = "png",
			Bool_t bKatarinaDiff = kTRUE,
			Bool_t bYouRef = kTRUE
	)
{
	/*
	TString sInput = "~/NBI/Codes/flow/results/V0s/9/merge/AnalysisResults.root";
	TString sOutput = "~/NBI/Codes/results/TPConly";
	TString sOutputFormat = "png";
	Bool_t bKatarinaDiff = kTRUE;
	Bool_t bYouRef = kTRUE;
	*/

	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareHistos.C");
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareRatio.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	fInput->cd("flowPID_JHEP");
	//fInput->ls();

	// ===== Loading input ===== 
	TList* fInputList = (TList*) gDirectory->Get("Tracks_flowPID_JHEP");

	// reference
	TProfile* fRefCorTwo2 = (TProfile*) fInputList->FindObject("fRefCorTwo2"); 
	TProfile* fRefCorTwo3 = (TProfile*) fInputList->FindObject("fRefCorTwo3"); 
	TProfile* fRefCorTwo4 = (TProfile*) fInputList->FindObject("fRefCorTwo4"); 
	TProfile* fRefCorTwo5 = (TProfile*) fInputList->FindObject("fRefCorTwo5"); 
	
	TProfile* fRefCorTwo2_Gap00 = (TProfile*) fInputList->FindObject("fRefCorTwo2_Gap00"); 
	TProfile* fRefCorTwo2_Gap04 = (TProfile*) fInputList->FindObject("fRefCorTwo2_Gap04"); 
	TProfile* fRefCorTwo2_Gap08 = (TProfile*) fInputList->FindObject("fRefCorTwo2_Gap08"); 
	TProfile* fRefCorTwo2_Gap10 = (TProfile*) fInputList->FindObject("fRefCorTwo2_Gap10"); 

	TH1D* hRefFlowTwo2 = (TH1D*) fRefCorTwo2->ProjectionX()->Clone("hRefFlowTwo2");
	TH1D* hRefFlowTwo3 = (TH1D*) fRefCorTwo3->ProjectionX()->Clone("hRefFlowTwo3");
	TH1D* hRefFlowTwo4 = (TH1D*) fRefCorTwo4->ProjectionX()->Clone("hRefFlowTwo4");
	TH1D* hRefFlowTwo5 = (TH1D*) fRefCorTwo5->ProjectionX()->Clone("hRefFlowTwo5");
	
	TH1D* hRefFlowTwo2_Gap00 = (TH1D*) fRefCorTwo2_Gap00->ProjectionX()->Clone("hRefFlowTwo2_Gap00");
	TH1D* hRefFlowTwo2_Gap04 = (TH1D*) fRefCorTwo2_Gap04->ProjectionX()->Clone("hRefFlowTwo2_Gap04");
	TH1D* hRefFlowTwo2_Gap08 = (TH1D*) fRefCorTwo2_Gap08->ProjectionX()->Clone("hRefFlowTwo2_Gap08");
	TH1D* hRefFlowTwo2_Gap10 = (TH1D*) fRefCorTwo2_Gap10->ProjectionX()->Clone("hRefFlowTwo2_Gap10");


	// differential 
	const Int_t iNumCentBins = 9;
	TProfile* fDiffCorTwo2[iNumCentBins];
	TProfile* fDiffCorTwo2_Gap00[iNumCentBins];
	TProfile* fDiffCorTwo2_Gap04[iNumCentBins];
	TProfile* fDiffCorTwo2_Gap08[iNumCentBins];
	TProfile* fDiffCorTwo2_Gap10[iNumCentBins];
	TProfile* fDiffCorTwo3[iNumCentBins];

	TH1D* hDiffFlowTwo2[iNumCentBins];
	TH1D* hDiffFlowTwo2_Gap00[iNumCentBins];
	TH1D* hDiffFlowTwo2_Gap04[iNumCentBins];
	TH1D* hDiffFlowTwo2_Gap08[iNumCentBins];
	TH1D* hDiffFlowTwo2_Gap10[iNumCentBins];
	TH1D* hDiffFlowTwo3[iNumCentBins];

	TList* lCent[iNumCentBins];

	for(Int_t i(0); i < iNumCentBins; i++)
	{
		lCent[i] = new TList();

		fDiffCorTwo2[i] = (TProfile*) fInputList->FindObject(Form("fDiffCorTwo2_Cent%d",i)); 
		fDiffCorTwo2_Gap00[i] = (TProfile*) fInputList->FindObject(Form("fDiffCorTwo2_Gap00_Cent%d",i)); 
		fDiffCorTwo2_Gap04[i] = (TProfile*) fInputList->FindObject(Form("fDiffCorTwo2_Gap04_Cent%d",i)); 
		fDiffCorTwo2_Gap08[i] = (TProfile*) fInputList->FindObject(Form("fDiffCorTwo2_Gap08_Cent%d",i)); 
		fDiffCorTwo2_Gap10[i] = (TProfile*) fInputList->FindObject(Form("fDiffCorTwo2_Gap10_Cent%d",i)); 
		fDiffCorTwo3[i] = (TProfile*) fInputList->FindObject(Form("fDiffCorTwo3_Cent%d",i)); 

		hDiffFlowTwo2[i] = (TH1D*) fDiffCorTwo2[i]->ProjectionX()->Clone(Form("hDiffFlowTwo2_Cent%d",i));
		hDiffFlowTwo2_Gap00[i] = (TH1D*) fDiffCorTwo2_Gap00[i]->ProjectionX()->Clone(Form("hDiffFlowTwo2_Gap00_Cent%d",i));
		hDiffFlowTwo2_Gap04[i] = (TH1D*) fDiffCorTwo2_Gap04[i]->ProjectionX()->Clone(Form("hDiffFlowTwo2_Gap04_Cent%d",i));
		hDiffFlowTwo2_Gap08[i] = (TH1D*) fDiffCorTwo2_Gap08[i]->ProjectionX()->Clone(Form("hDiffFlowTwo2_Gap08_Cent%d",i));
		hDiffFlowTwo2_Gap10[i] = (TH1D*) fDiffCorTwo2_Gap10[i]->ProjectionX()->Clone(Form("hDiffFlowTwo2_Gap10_Cent%d",i));
		hDiffFlowTwo3[i] = (TH1D*) fDiffCorTwo3[i]->ProjectionX()->Clone(Form("hDiffFlowTwo3_Cent%d",i));

		lCent[i]->Add(hDiffFlowTwo2[i]);
		lCent[i]->Add(hDiffFlowTwo2_Gap00[i]);
		lCent[i]->Add(hDiffFlowTwo2_Gap04[i]);
		lCent[i]->Add(hDiffFlowTwo2_Gap08[i]);
		lCent[i]->Add(hDiffFlowTwo2_Gap10[i]);
		lCent[i]->Add(hDiffFlowTwo3[i]);
	}

	// loaded



	// ===== Making flow ===== 
	const Int_t iNumBins = hRefFlowTwo2->GetNbinsX();
	for(Int_t j(1); j < iNumBins+1; j++)
	{
		hRefFlowTwo2->SetBinContent(j,TMath::Sqrt(hRefFlowTwo2->GetBinContent(j)));
		hRefFlowTwo3->SetBinContent(j,TMath::Sqrt(hRefFlowTwo3->GetBinContent(j)));
		hRefFlowTwo4->SetBinContent(j,TMath::Sqrt(hRefFlowTwo4->GetBinContent(j)));
		hRefFlowTwo5->SetBinContent(j,TMath::Sqrt(hRefFlowTwo5->GetBinContent(j)));
		
		hRefFlowTwo2_Gap00->SetBinContent(j,TMath::Sqrt(hRefFlowTwo2_Gap00->GetBinContent(j)));
		hRefFlowTwo2_Gap04->SetBinContent(j,TMath::Sqrt(hRefFlowTwo2_Gap04->GetBinContent(j)));
		hRefFlowTwo2_Gap08->SetBinContent(j,TMath::Sqrt(hRefFlowTwo2_Gap08->GetBinContent(j)));
		hRefFlowTwo2_Gap10->SetBinContent(j,TMath::Sqrt(hRefFlowTwo2_Gap10->GetBinContent(j)));
	}

	// hRefFlow = vn{2}
	Double_t dDiff = 0, dRef = 0;


	const Int_t iNumPtBins = hDiffFlowTwo2[0]->GetNbinsX();
	
	for(Int_t j(0); j < iNumCentBins; j++)
	{
		for(Int_t i(1); i < iNumPtBins+1; i++)
		{
			dRef = hRefFlowTwo2->GetBinContent(j+1);
			dDiff = hDiffFlowTwo2[j]->GetBinContent(i);
			hDiffFlowTwo2[j]->SetBinContent(i,dDiff/dRef);
			
			dRef = hRefFlowTwo2_Gap00->GetBinContent(j+1);
			dDiff = hDiffFlowTwo2_Gap00[j]->GetBinContent(i);
			hDiffFlowTwo2_Gap00[j]->SetBinContent(i,dDiff/dRef);
			
			dRef = hRefFlowTwo2_Gap04->GetBinContent(j+1);
			dDiff = hDiffFlowTwo2_Gap04[j]->GetBinContent(i);
			hDiffFlowTwo2_Gap04[j]->SetBinContent(i,dDiff/dRef);

			dRef = hRefFlowTwo2_Gap08->GetBinContent(j+1);
			dDiff = hDiffFlowTwo2_Gap08[j]->GetBinContent(i);
			hDiffFlowTwo2_Gap08[j]->SetBinContent(i,dDiff/dRef);

			dRef = hRefFlowTwo2_Gap10->GetBinContent(j+1);
			dDiff = hDiffFlowTwo2_Gap10[j]->GetBinContent(i);
			hDiffFlowTwo2_Gap10[j]->SetBinContent(i,dDiff/dRef);

			dRef = hRefFlowTwo3->GetBinContent(j+1);
			dDiff = hDiffFlowTwo3[j]->GetBinContent(i);
			hDiffFlowTwo3[j]->SetBinContent(i,dDiff/dRef);			
		}
	}

	// hDiff = vn'{2}


	TString sLabel[] = {hDiffFlowTwo2[0]->GetTitle(), hDiffFlowTwo2_Gap00[0]->GetTitle(), hDiffFlowTwo2_Gap04[0]->GetTitle(), hDiffFlowTwo2_Gap08[0]->GetTitle(), hDiffFlowTwo2_Gap10[0]->GetTitle(),hDiffFlowTwo3[0]->GetTitle()};
	TCanvas* cTemp;
	for (Int_t i = 0; i < iNumCentBins; ++i)
	{
		cTemp = CompareHistos(lCent[i],sLabel,0.,0.3,0);
		cTemp->Print(Form("%s/FlowTwo2_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
	}

	// ===== My v2 done !!! =====

	// ===== You ref flow comparison =====
	if(bYouRef)
	{
		TString sInputYou = "~/NBI/Codes/results/You";

		TFile* fYouRef = new TFile(Form("%s/Analysis_vn_2468.root",sInputYou.Data()),"READ");
		fYouRef->cd();

		// loading
		TH1D* hYouRefTwo2 = (TH1D*) gDirectory->Get("fv22NoGap")->Clone("fYouRefTwo2");
		TH1D* hYouRefTwo3 = (TH1D*) gDirectory->Get("fv32NoGap")->Clone("fYouRefTwo3");
		TH1D* hYouRefTwo4 = (TH1D*) gDirectory->Get("fv42NoGap")->Clone("fYouRefTwo4");

		TH1D* hYouRefTwo2_Gap00 = (TH1D*) gDirectory->Get("fv22Gap00")->Clone("fYouRefTwo2_Gap00");
		TH1D* hYouRefTwo2_Gap04 = (TH1D*) gDirectory->Get("fv22Gap04")->Clone("fYouRefTwo2_Gap04");
		TH1D* hYouRefTwo2_Gap08 = (TH1D*) gDirectory->Get("fv22Gap08")->Clone("fYouRefTwo2_Gap08");
		TH1D* hYouRefTwo2_Gap10 = (TH1D*) gDirectory->Get("fv22Gap10")->Clone("fYouRefTwo2_Gap10");

		// making ratios
		TCanvas* cTemp;
		cTemp = CompareRatio(hRefFlowTwo2,hYouRefTwo2);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo2.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

		cTemp = CompareRatio(hRefFlowTwo3,hYouRefTwo3);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo3.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());
		
		cTemp = CompareRatio(hRefFlowTwo4,hYouRefTwo4);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo4.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

		cTemp = CompareRatio(hRefFlowTwo2_Gap00,hYouRefTwo2_Gap00);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo2_Gap00.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());
		
		cTemp = CompareRatio(hRefFlowTwo2_Gap04,hYouRefTwo2_Gap04);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo2_Gap04.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

		cTemp = CompareRatio(hRefFlowTwo2_Gap08,hYouRefTwo2_Gap08);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo2_Gap08.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());

		cTemp = CompareRatio(hRefFlowTwo2_Gap10,hYouRefTwo2_Gap10);
		cTemp->Print(Form("%s/CompYouRef/RatioMeYou_RefTwo2_Gap10.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());
		
		// You pt diff flow comparison
		TFile* fInputYouDiff2 = new TFile("~/NBI/Codes/results/You/Analysis_v2pt.root","READ");
		fInputYouDiff2->cd();

		TH1D* hYouDiffFlowTwo2[iNumCentBins];
		TH1D* hYouDiffFlowTwo2_Gap00[iNumCentBins];
		TH1D* hYouDiffFlowTwo2_Gap04[iNumCentBins];
		TH1D* hYouDiffFlowTwo2_Gap08[iNumCentBins];
		TH1D* hYouDiffFlowTwo2_Gap10[iNumCentBins];
		for(Int_t i(0); i < iNumCentBins; i++)
		{
			hYouDiffFlowTwo2[i] = (TH1D*) gDirectory->Get(Form("hisv22NoGap_cent%d",i))->Clone(Form("hYouDiffFlowTwo2_Cent%d",i));
			hYouDiffFlowTwo2_Gap00[i] = (TH1D*) gDirectory->Get(Form("hisv22Gap00_cent%d",i))->Clone(Form("hYouDiffFlowTwo2_Gap00_Cent%d",i));
			hYouDiffFlowTwo2_Gap04[i] = (TH1D*) gDirectory->Get(Form("hisv22Gap04_cent%d",i))->Clone(Form("hYouDiffFlowTwo2_Gap04_Cent%d",i));
			hYouDiffFlowTwo2_Gap08[i] = (TH1D*) gDirectory->Get(Form("hisv22Gap08_cent%d",i))->Clone(Form("hYouDiffFlowTwo2_Gap08_Cent%d",i));
			hYouDiffFlowTwo2_Gap10[i] = (TH1D*) gDirectory->Get(Form("hisv22Gap10_cent%d",i))->Clone(Form("hYouDiffFlowTwo2_Gap10_Cent%d",i));
		}

		TCanvas* cTemp;
		for(Int_t i(0); i < iNumCentBins; i++)
		{
			cTemp = CompareRatio(hDiffFlowTwo2[i], hYouDiffFlowTwo2[i]);
			cTemp->Print(Form("%s/CompYouDiff/RatioMeYou_DiffTwo2_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap00[i], hYouDiffFlowTwo2_Gap00[i]);
			cTemp->Print(Form("%s/CompYouDiff/RatioMeYou_DiffTwo2_Gap00_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap04[i], hYouDiffFlowTwo2_Gap04[i]);
			cTemp->Print(Form("%s/CompYouDiff/RatioMeYou_DiffTwo2_Gap04_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap08[i], hYouDiffFlowTwo2_Gap08[i]);
			cTemp->Print(Form("%s/CompYouDiff/RatioMeYou_DiffTwo2_Gap08_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap10[i], hYouDiffFlowTwo2_Gap10[i]);
			cTemp->Print(Form("%s/CompYouDiff/RatioMeYou_DiffTwo2_Gap10_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
		}

	}

	// ===== Katarina comparison =====
	if(bKatarinaDiff)
	{
		TString sInputKatarina = "~/NBI/Codes/results/Katarina";

		TFile* fInputKatarinaDiff = new TFile(Form("%s/OutFileV2_TPConly_hist.root",sInputKatarina.Data()),"READ");
		fInputKatarinaDiff->cd();

		// loading
		TH1D* hKatarinaDiffFlowTwo2[iNumCentBins];
		TH1D* hKatarinaDiffFlowTwo2_Gap00[iNumCentBins];
		TH1D* hKatarinaDiffFlowTwo2_Gap04[iNumCentBins];
		TH1D* hKatarinaDiffFlowTwo2_Gap08[iNumCentBins];
		TH1D* hKatarinaDiffFlowTwo2_Gap10[iNumCentBins];
		for(Int_t i(0); i < iNumCentBins; i++)
		{
			hKatarinaDiffFlowTwo2[i] = (TH1D*) gDirectory->Get(Form("HChd22ReX_cent%d",i))->Clone(Form("hKatarinaDiffFlowTwo2_Cent%d",i));
			hKatarinaDiffFlowTwo2_Gap00[i] = (TH1D*) gDirectory->Get(Form("HChd22Gap0QMReX_cent%d",i))->Clone(Form("hKatarinaDiffFlowTwo2_Gap00_Cent%d",i));
			hKatarinaDiffFlowTwo2_Gap04[i] = (TH1D*) gDirectory->Get(Form("HChd22Gap4QMReX_cent%d",i))->Clone(Form("hKatarinaDiffFlowTwo2_Gap04_Cent%d",i));
			hKatarinaDiffFlowTwo2_Gap08[i] = (TH1D*) gDirectory->Get(Form("HChd22Gap8QMReX_cent%d",i))->Clone(Form("hKatarinaDiffFlowTwo2_Gap08_Cent%d",i));
			hKatarinaDiffFlowTwo2_Gap10[i] = (TH1D*) gDirectory->Get(Form("HChd22GapQMReX_cent%d",i))->Clone(Form("hKatarinaDiffFlowTwo2_Gap10_Cent%d",i));
		}

		TCanvas* cTemp;
		for(Int_t i(0); i < iNumCentBins; i++)
		{
			cTemp = CompareRatio(hDiffFlowTwo2[i],hKatarinaDiffFlowTwo2[i]);
			cTemp->Print(Form("%s/CompKatarinaDiff/RatioMeKatarina_DiffTwo2_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap00[i],hKatarinaDiffFlowTwo2_Gap00[i]);
			cTemp->Print(Form("%s/CompKatarinaDiff/RatioMeKatarina_DiffTwo2_Gap00_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap04[i],hKatarinaDiffFlowTwo2_Gap04[i]);
			cTemp->Print(Form("%s/CompKatarinaDiff/RatioMeKatarina_DiffTwo2_Gap04_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap08[i],hKatarinaDiffFlowTwo2_Gap08[i]);
			cTemp->Print(Form("%s/CompKatarinaDiff/RatioMeKatarina_DiffTwo2_Gap08_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
			
			cTemp = CompareRatio(hDiffFlowTwo2_Gap10[i],hKatarinaDiffFlowTwo2_Gap10[i]);
			cTemp->Print(Form("%s/CompKatarinaDiff/RatioMeKatarina_DiffTwo2_Gap10_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
		}
	}	
	
	return;
}

TH1D* RedoFlow(TH1D* fInput, const char* name)
{
	Double_t dBinsX[] = {0.,5.,10.,20,30,40,50,60,70,80}; 
	TH1D* hOut = new TH1D(Form("h%s",name),name,9,dBinsX);
	

	for(Int_t i = 1; i < 10; i++)
	{
		hOut->SetBinContent(i,fInput->GetBinContent(i));
		hOut->SetBinError(i,fInput->GetBinError(i));
	}

	hOut->SetMinimum(0.);
	hOut->SetMaximum(0.25);
	return hOut;

}

TH1D* MakeRatio(TH1D* fNom, TH1D* fDenom, const char* name)
{
	TH1D* fOut = (TH1D*) fNom->Clone(Form("h%s",name));
	fOut->SetTitle(name);
	fOut->Divide(fDenom);

	fOut->SetMaximum(1.05);
	fOut->SetMinimum(0.95);
	return fOut;
}

