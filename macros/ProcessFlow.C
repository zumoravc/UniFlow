TH1D* RedoFlow(TH1D* fInput, const char* name);
TH1D* MakeRatio(TH1D* fNom, TH1D* fDenom, const char* name);


void ProcessFlow()
{
	TString sInput = "~/NBI/Codes/flow/AnalysisResults.root";
	

	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	fInput->cd("FlowPID");
	//fInput->ls();

	// ===== Loading input ===== 
	TList* fInputList = (TList*) gDirectory->Get("FlowPID");

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
	CompareHistos(lCent[0],sLabel,-0.3,0.3,0);

	return;
	
	TLine* fLine = new TLine(-0.5,0,9.5,0);
	
	TList* lCum = new TList();
	lCum->Add(hProj2->Clone("hProj2Cum"));
	lCum->Add(hProj2Gap04->Clone("hProj2Gap04Cum"));
	lCum->Add(hProj2Gap08->Clone("hProj2Gap08Cum"));
	lCum->Add(hProj2Gap10->Clone("hProj2Gap10Cum"));
	lCum->Add(hProj3->Clone("hProj3Cum"));
	lCum->Add(hProj4->Clone("hProj4Cum"));
	lCum->Add(hProj5->Clone("hProj5Cum"));
	

	TString sCum[] = {"c2","c2","c2Gap04","c2Gap08","c2Gap10","c3","c4","c5"};
	//CompareHistos(lCum,sCum,-0.01,0.05,0);
	//fLine->Draw("same");

	Int_t iNbins = hProj2->GetNbinsX();
	for(Int_t i = 1; i < iNbins+1; i++)
	{
		hProj2->SetBinContent(i,TMath::Sqrt(hProj2->GetBinContent(i)));
		hProj2Gap00->SetBinContent(i,TMath::Sqrt(hProj2Gap00->GetBinContent(i)));
		hProj2Gap04->SetBinContent(i,TMath::Sqrt(hProj2Gap04->GetBinContent(i)));
		hProj2Gap08->SetBinContent(i,TMath::Sqrt(hProj2Gap08->GetBinContent(i)));
		hProj2Gap10->SetBinContent(i,TMath::Sqrt(hProj2Gap10->GetBinContent(i)));
	
		hProj3->SetBinContent(i,TMath::Sqrt(hProj3->GetBinContent(i)));
		//hProj4->SetBinContent(i,TMath::Sqrt(hProj4->GetBinContent(i)));
		//hProj5->SetBinContent(i,TMath::Sqrt(hProj5->GetBinContent(i)));
	}
	
	/*
	TCanvas* cTest = new TCanvas("cTest","TE|ST");
	cTest->cd();
	hProj4->Draw();
	*/

	TList* lCompare = new TList();
	lCompare->Add(hProj2);
	lCompare->Add(hProj2Gap04);
	lCompare->Add(hProj2Gap08);
	lCompare->Add(hProj2Gap10);
	lCompare->Add(hProj3);
	//lCompare->Add(hProj4);
	//lCompare->Add(hProj5);

	TString sData[] = {"v2{2}","v2{2,|#Delta#eta|>0.4}","v2{2,|#Delta#eta|>0.8}","v2{2,|#Delta#eta|>1}","v3{2}"/*,"c4{2}","c5{2}"*/};

	//CompareHistos(lCompare,sData,-0.00,0.25,0);
	//fLine->Draw("same");


	TH1D* rebinProj2 = RedoFlow(hProj2,"v2{2} Me (rebin)");
	TH1D* rebinProj2Gap00 = RedoFlow(hProj2Gap00,"v2{2,|#Delta#eta|>0} Me (rebin)");
	TH1D* rebinProj2Gap04 = RedoFlow(hProj2Gap04,"v2{2,|#Delta#eta|>0.4} Me (rebin)");
	TH1D* rebinProj2Gap08 = RedoFlow(hProj2Gap08,"v2{2,|#Delta#eta|>0.8} Me (rebin)");
	TH1D* rebinProj2Gap10 = RedoFlow(hProj2Gap10,"v2{2,|#Delta#eta|>1} Me (rebin)");
	TH1D* rebinProj3 = RedoFlow(hProj3,"v3{2} (rebin)");
	
	TString sRebin[] = {"v2{2}","v2{2,|#Delta#eta|>0}","v2{2,|#Delta#eta|>0.4}","v2{2,|#Delta#eta|>0.8}","v2{2,|#Delta#eta|>1}","v3{2}"/*,"c4{2}","c5{2}"*/};

	TList* lRebin = new TList();
	lRebin->Add(rebinProj2);
	lRebin->Add(rebinProj2Gap00);
	lRebin->Add(rebinProj2Gap04);
	lRebin->Add(rebinProj2Gap08);
	lRebin->Add(rebinProj2Gap10);
	lRebin->Add(rebinProj3);

	CompareHistos(lRebin,sRebin,0.,0.25,0);
	
	TFile* fInputYou = new TFile("~/NBI/Codes/macros/you/Analysis_vn_2468.root","READ");
	fInputYou->cd();

	TH1D* hYouV22 = (TH1D*) gDirectory->Get("fv22NoGap")->Clone("hYouV22");
	TH1D* hYouV22Gap04 = (TH1D*) gDirectory->Get("fv22Gap04")->Clone("hYouV22Gap04");
	TH1D* hYouV22Gap08 = (TH1D*) gDirectory->Get("fv22Gap08")->Clone("hYouV22Gap08");
	TH1D* hYouV22Gap1 = (TH1D*) gDirectory->Get("fv22Gap10")->Clone("hYouV22Gap1");
	TH1D* hYouV32 = (TH1D*) gDirectory->Get("fv32NoGap")->Clone("hYouV32");



	TList* lYou = new TList();
	lYou->Add(hYouV22);
	lYou->Add(hYouV22Gap04);
	lYou->Add(hYouV22Gap08);
	lYou->Add(hYouV22Gap1);
	lYou->Add(hYouV32);

	TString sYou[] = {"v2{2}","v2{2,|#Delta#eta|>0.4}","v2{2,|#Delta#eta|>0.8}","v2{2,|#Delta#eta|>1}","v3{2}"};
	CompareHistos(lYou,sYou,0,0.25,0);
	fLine->Draw("same");

/*
	hProj2->Draw();
	hProj3->Draw("same");
	hProj4->Draw("same");
	hProj5->Draw("same");
*/

	// ratio plots 

	return;
	TLine* lUnity = new TLine(0.,1.,80.,1);

	TH1D* hRatioProj2 = (TH1D*) rebinProj2->Clone("hRatioProj2");
	//hRatioProj2->Sumw2();
	hRatioProj2->SetTitle("v2{2} Me/You");
	hRatioProj2->SetMaximum(1.05);
	hRatioProj2->SetMinimum(0.95);
	hRatioProj2->Divide(hYouV22);

	TH1D* hRatioProj3 = (TH1D*) rebinProj3->Clone("hRatioProj3");
	//hRatioProj2->Sumw2();
	hRatioProj3->SetTitle("v3{2} Me/You");
	hRatioProj3->SetMaximum(1.05);
	hRatioProj3->SetMinimum(0.95);
	hRatioProj3->Divide(hYouV32);

	TH1D* hRatioProj2Gap04 = (TH1D*) rebinProj2Gap04->Clone("hRatioProj2Gap04");
	//hRatioProj2->Sumw2();
	hRatioProj2Gap04->SetTitle("v2{2} Gap 04 Me/You");
	hRatioProj2Gap04->SetMaximum(1.05);
	hRatioProj2Gap04->SetMinimum(0.95);
	hRatioProj2Gap04->Divide(hYouV22Gap04);


	TH1D* hRatioProj2Gap08 = (TH1D*) rebinProj2Gap08->Clone("hRatioProj2Gap08");
	//hRatioProj2->Sumw2();
	hRatioProj2Gap08->SetTitle("v2{2} Gap 08 Me/You");
	hRatioProj2Gap08->SetMaximum(1.05);
	hRatioProj2Gap08->SetMinimum(0.95);
	hRatioProj2Gap08->Divide(hYouV22Gap08);

	TH1D* hRatioProj2Gap10 = (TH1D*) rebinProj2Gap10->Clone("hRatioProj2Gap10");
	//hRatioProj2->Sumw2();
	hRatioProj2Gap10->SetTitle("v2{2} Gap 1 Me/You");
	hRatioProj2Gap10->SetMaximum(1.05);
	hRatioProj2Gap10->SetMinimum(0.95);
	hRatioProj2Gap10->Divide(hYouV22Gap1);
	
	TCanvas* Ratio = new TCanvas("Ratio","Ratio");
	Ratio->Divide(1,2);
	Ratio->cd(1);

	rebinProj2->SetLineColor(kRed);
	rebinProj2->SetMarkerColor(kRed);
	rebinProj2->Draw();
	hYouV22->Draw("same");
	
	Ratio->cd(2);
	hRatioProj2->Draw();
	lUnity->Draw("same");
	
	TCanvas* Ratiov3 = new TCanvas("Ratiov3","RatioV3");
	Ratiov3->Divide(1,2);
	Ratiov3->cd(1);

	rebinProj3->SetLineColor(kRed);
	rebinProj3->SetMarkerColor(kRed);
	rebinProj3->Draw();
	hYouV32->Draw("same");
	
	Ratiov3->cd(2);
	hRatioProj3->Draw();
	lUnity->Draw("same");
	
	TCanvas* RatioV2Gap04 = new TCanvas("Ratiov2Gap04","RatioV2Gap04");
	RatioV2Gap04->Divide(1,2);
	RatioV2Gap04->cd(1);

	rebinProj2Gap04->SetLineColor(kRed);
	rebinProj2Gap04->SetMarkerColor(kRed);
	rebinProj2Gap04->Draw();
	hYouV22Gap04->Draw("same");
	
	RatioV2Gap04->cd(2);
	hRatioProj2Gap04->Draw();
	lUnity->Draw("same");

		TCanvas* RatioV2Gap08 = new TCanvas("Ratiov2Gap08","RatioV2Gap08");
	RatioV2Gap08->Divide(1,2);
	RatioV2Gap08->cd(1);

	rebinProj2Gap08->SetLineColor(kRed);
	rebinProj2Gap08->SetMarkerColor(kRed);
	rebinProj2Gap08->Draw();
	hYouV22Gap08->Draw("same");
	
	RatioV2Gap08->cd(2);
	hRatioProj2Gap08->Draw();
	lUnity->Draw("same");

	TCanvas* RatioV2Gap10 = new TCanvas("Ratiov2Gap10","RatioV2Gap10");
	RatioV2Gap10->Divide(1,2);
	RatioV2Gap10->cd(1);

	rebinProj2Gap10->SetLineColor(kRed);
	rebinProj2Gap10->SetMarkerColor(kRed);
	rebinProj2Gap10->Draw();
	hYouV22Gap1->Draw("same");
	
	RatioV2Gap10->cd(2);
	hRatioProj2Gap10->Draw();
	lUnity->Draw("same");

	Ratio->Print("~/NBI/Codes/flow/results/t5/plots/RatioV2.pdf","pdf");
	Ratiov3->Print("~/NBI/Codes/flow/results/t5/plots/RatioV3.pdf","pdf");
	RatioV2Gap04->Print("~/NBI/Codes/flow/results/t5/plots/RatioV2Gap04.pdf","pdf");
	RatioV2Gap08->Print("~/NBI/Codes/flow/results/t5/plots/RatioV2Gap08.pdf","pdf");
	RatioV2Gap10->Print("~/NBI/Codes/flow/results/t5/plots/RatioV2Gap10.pdf","pdf");



	/*
	TH1D* hRatioV2 = MakeRatio(rebinProj2,hYouV22,"ratioV2");
	TH1D* hRatioV2Gap04 = MakeRatio(rebinProj2Gap04,hYouV22Gap04,"ratioV2Gap04");
	TH1D* hRatioV2Gap08 = MakeRatio(rebinProj2Gap08,hYouV22Gap08,"ratioV2Gap08");
	TH1D* hRatioV2Gap10 = MakeRatio(rebinProj2Gap10,hYouV22Gap1,"ratioV2Gap10");
	TH1D* hRatioV3 = MakeRatio(rebinProj3,hYouV32,"ratioV3");

	TCanvas* cRatio = new TCanvas("cRatio","cRatio");
	cRatio->cd();
	hRatioV3->Draw();


	TList* lRatio = new TList();
	lRatio->Add(hRatioV2);
	lRatio->Add(hRatioV2Gap04);
	lRatio->Add(hRatioV2Gap08);
	lRatio->Add(hRatioV2Gap10);
	lRatio->Add(hRatioV3);
	
	TString sRatio[] = {"v2{2}","v2{2,|#Delta#eta|>0.4}","v2{2,|#Delta#eta|>0.8}","v2{2,|#Delta#eta|>1}","v3{2}"};

	CompareHistos(lRatio,sRatio,0.5,1.5,0);
	*/
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
	fOut->SetTitle(Form("%s (Me/You)",name));
	fOut->Divide(hYouV22);

	fOut->SetMaximum(1.05);
	fOut->SetMinimum(0.95);
	return fOut;
}