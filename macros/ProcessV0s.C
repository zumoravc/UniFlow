void ProcessV0s(
		const TString sInput = "~/NBI/Flow/results/V0s/5/plusplus/merge/AnalysisResults_merged.root",
		const TString sOutput = "~/NBI/Flow/results/V0s/5/plusplus/plots",
		const TString sTag = "_JHEP",
		const TString sEtaGap = "Gap09",
		const TString sOutputFormat = "png"
	)
{
	//const TString sInput = "~/NBI/Codes/results/V0s/5/plusplus/merge/AnalysisResults_merged.root";
	//const TString sOutput = "~/NBI/Codes/results/V0s/5/plusplus/plots";
	//const TString sOutputFormat = "png";
	//const TString sEtaGap = "Gap09";
 
	const Int_t iNumPtBins = 10; // pT bins
	const Int_t iNumCentBins = 9; // centrality bins
	// bins edges
	Double_t fPtBinEdges[] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6}; 
	Double_t fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

	// =======================================
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	TFile* fOutput = new TFile(Form("%s/V0sFlow.root",sOutput.Data()),"RECREATE");

	fInput->cd(Form("%s",sTag.Data()));

	// ===== Loading input ===== 
	TList* lInputTracks = (TList*) gDirectory->Get(Form("Tracks_%s",sTag.Data()));
	TList* lInputV0s = (TList*) gDirectory->Get(Form("V0s_%s",sTag.Data()));
	//TList* lInputTracks = (TList*) gDirectory->Get(Form("Tracks",sTag.Data()));
	//TList* lInputV0s = (TList*) gDirectory->Get(Form("V0s"));

	// reference 
	TProfile* pRefCorTwo2_Gap09 = (TProfile*) (lInputTracks->FindObject(Form("fRefCorTwo2_%s",sEtaGap.Data())) )->Clone(Form("pRefCorTwo2_%s",sEtaGap.Data())); 
	TH1D* hRefCorTwo2_Gap09 = (TH1D*) pRefCorTwo2_Gap09->ProjectionX()->Clone(Form("hRefCorTwo2_%s",sEtaGap.Data()));
	
	// V0s 
	TProfile2D* p2V0sDiffTwo2_Gap09P_K0s[iNumCentBins];
	TProfile2D* p2V0sDiffTwo2_Gap09N_K0s[iNumCentBins];
	TProfile2D* p2V0sDiffTwo2_Gap09P_Lambda[iNumCentBins];
	TProfile2D* p2V0sDiffTwo2_Gap09N_Lambda[iNumCentBins];

	TH2D* h2InvMass_Gap09_K0s[iNumCentBins];
	TH2D* h2InvMass_Gap09_Lambda[iNumCentBins];
	
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		h2InvMass_Gap09_K0s[i] = (TH2D*) (lInputV0s->FindObject(Form("fV0sK0s_%s_Cent%d",sEtaGap.Data(),i)))->Clone(Form("h2InvMass_%s_K0s_Cent%d",sEtaGap.Data(),i));
		h2InvMass_Gap09_Lambda[i] = (TH2D*) (lInputV0s->FindObject(Form("fV0sLambda_%s_Cent%d",sEtaGap.Data(),i)))->Clone(Form("h2InvMass_%s_Lambda_Cent%d",sEtaGap.Data(),i));

		p2V0sDiffTwo2_Gap09P_K0s[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwo2_%sP_K0s_Cent%d",sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwo2_%sP_K0s_Cent%d",sEtaGap.Data(),i));
		p2V0sDiffTwo2_Gap09N_K0s[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwo2_%sN_K0s_Cent%d",sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwo2_%sN_K0s_Cent%d",sEtaGap.Data(),i));
	
		p2V0sDiffTwo2_Gap09P_Lambda[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwo2_%sP_Lambda_Cent%d",sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwo2_%sP_Lambda_Cent%d",sEtaGap.Data(),i));
		p2V0sDiffTwo2_Gap09N_Lambda[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwo2_%sN_Lambda_Cent%d",sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwo2_%sN_Lambda_Cent%d",sEtaGap.Data(),i));
	}

	// ===== Making projections ======
	TH1D* hInvMass_K0s[iNumCentBins][iNumPtBins];
	TH1D* hInvMass_Lambda[iNumCentBins][iNumPtBins];

	TH1D* hCummMass_Pos_K0s[iNumCentBins][iNumPtBins];
	TH1D* hCummMass_Neg_K0s[iNumCentBins][iNumPtBins];
	TH1D* hCummMass_Pos_Lambda[iNumCentBins][iNumPtBins];
	TH1D* hCummMass_Neg_Lambda[iNumCentBins][iNumPtBins];

	TList* lInvMass_K0s = new TList();
	TList* lInvMass_Lambda = new TList();
	TList* lCummMass_Pos_K0s = new TList();
	TList* lCummMass_Neg_K0s = new TList();
	TList* lCummMass_Pos_Lambda = new TList();
	TList* lCummMass_Neg_Lambda = new TList();
	
	TCanvas* cTemp = new TCanvas();
	cTemp->cd();
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			// inv mass plots
			hInvMass_K0s[i][j] = (TH1D*) h2InvMass_Gap09_K0s[i]->ProjectionY(Form("hInvMass_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hInvMass_K0s[i][j]->SetTitle(Form("K_{S}^{0}: InvMass |#it{#eta}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hInvMass_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/InvMassK0s/InvMass_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lInvMass_K0s->Add(hInvMass_K0s[i][j]);

			hInvMass_Lambda[i][j] = (TH1D*) h2InvMass_Gap09_Lambda[i]->ProjectionY(Form("hInvMass_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hInvMass_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: InvMass |#it{#eta}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hInvMass_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/InvMassLambda/InvMass_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lInvMass_Lambda->Add(hInvMass_Lambda[i][j]);
			
			// cum x inv mass plots
			hCummMass_Pos_K0s[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09P_K0s[i]->ProjectionY(Form("hCummMass_Pos_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Pos_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCummMass_Pos_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassK0s/CummMass_Pos_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCummMass_Pos_K0s->Add(hCummMass_Pos_K0s[i][j]);
			
			hCummMass_Neg_K0s[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09N_K0s[i]->ProjectionY(Form("hCummMass_Neg_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Neg_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|<-0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCummMass_Neg_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassK0s/CummMass_Neg_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCummMass_Neg_K0s->Add(hCummMass_Neg_K0s[i][j]);

			hCummMass_Pos_Lambda[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09P_Lambda[i]->ProjectionY(Form("hCummMass_Pos_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Pos_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCummMass_Pos_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassLambda/CummMass_Pos_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCummMass_Pos_Lambda->Add(hCummMass_Pos_Lambda[i][j]);
			
			hCummMass_Neg_Lambda[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09N_Lambda[i]->ProjectionY(Form("hCummMass_Neg_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Neg_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|<-0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCummMass_Neg_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassLambda/CummMass_Neg_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCummMass_Neg_Lambda->Add(hCummMass_Neg_Lambda[i][j]);
		}
	}

	// ===== Making Flow ======
	TF1* fUnity = new TF1("fUnity","1",0.,2.);

	TH1D* hFlowMass_Pos_K0s[iNumCentBins][iNumPtBins];
	TH1D* hFlowMass_Neg_K0s[iNumCentBins][iNumPtBins];
	TH1D* hFlowMass_Pos_Lambda[iNumCentBins][iNumPtBins];
	TH1D* hFlowMass_Neg_Lambda[iNumCentBins][iNumPtBins];

	TH1D* hFlowMass_K0s[iNumCentBins][iNumPtBins];
	TH1D* hFlowMass_Lambda[iNumCentBins][iNumPtBins];

	TList* lFlowMass_K0s = new TList();
	TList* lFlowMass_Pos_K0s = new TList();
	TList* lFlowMass_Neg_K0s = new TList();
	TList* lFlowMass_Lambda = new TList();
	TList* lFlowMass_Pos_Lambda = new TList();
	TList* lFlowMass_Neg_Lambda = new TList();
	
	Double_t dRefFlow = 0;

	for(Int_t i(0); i < iNumCentBins; i++)
	{
		dRefFlow = TMath::Sqrt(hRefCorTwo2_Gap09->GetBinContent(i+1));

		for(Int_t j(0); j < iNumPtBins; j++)
		{
			hFlowMass_Pos_K0s[i][j] = (TH1D*) hCummMass_Pos_K0s[i][j]->Clone(Form("hFlowMass_Pos_K0s_Cent%d_pt%d",i,j));
			hFlowMass_Neg_K0s[i][j] = (TH1D*) hCummMass_Neg_K0s[i][j]->Clone(Form("hFlowMass_Neg_K0s_Cent%d_pt%d",i,j));
			hFlowMass_Pos_Lambda[i][j] = (TH1D*) hCummMass_Pos_Lambda[i][j]->Clone(Form("hFlowMass_Pos_Lambda_Cent%d_pt%d",i,j));
			hFlowMass_Neg_Lambda[i][j] = (TH1D*) hCummMass_Neg_Lambda[i][j]->Clone(Form("hFlowMass_Neg_Lambda_Cent%d_pt%d",i,j));

			hFlowMass_Pos_K0s[i][j]->Divide(fUnity,dRefFlow);
			hFlowMass_Neg_K0s[i][j]->Divide(fUnity,dRefFlow);
			hFlowMass_Pos_Lambda[i][j]->Divide(fUnity,dRefFlow);
			hFlowMass_Neg_Lambda[i][j]->Divide(fUnity,dRefFlow);

			hFlowMass_Pos_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #it{v}_{2} %s |#it{#eta}^{POI}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",sEtaGap.Data(),fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hFlowMass_Neg_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #it{v}_{2} %s |#it{#eta}^{POI}|<-0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",sEtaGap.Data(),fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hFlowMass_Pos_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #it{v}_{2} %s |#it{#eta}^{POI}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",sEtaGap.Data(),fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hFlowMass_Neg_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #it{v}_{2} %s |#it{#eta}^{POI}|<-0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",sEtaGap.Data(),fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			
			hFlowMass_Pos_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassK0s/FlowMass_Pos_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lFlowMass_Pos_K0s->Add(hFlowMass_Pos_K0s[i][j]);

			hFlowMass_Neg_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassK0s/FlowMass_Neg_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lFlowMass_Neg_K0s->Add(hFlowMass_Neg_K0s[i][j]);
			
			hFlowMass_Pos_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassLambda/FlowMass_Pos_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lFlowMass_Pos_Lambda->Add(hFlowMass_Pos_Lambda[i][j]);
			
			hFlowMass_Neg_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassLambda/FlowMass_Neg_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lFlowMass_Neg_Lambda->Add(hFlowMass_Neg_Lambda[i][j]);

			// Making average of positive & negative eta POI
			hFlowMass_K0s[i][j] = (TH1D*) hFlowMass_Pos_K0s[i][j]->Clone(Form("hFlowMass_K0s_Cent%d_pt%d",i,j));
			hFlowMass_K0s[i][j]->Add(hFlowMass_Neg_K0s[i][j]);
			hFlowMass_K0s[i][j]->Divide(fUnity,2);
			hFlowMass_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #it{v}_{2} %s %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",sEtaGap.Data(),fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hFlowMass_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassK0s/FlowMass_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lFlowMass_K0s->Add(hFlowMass_K0s[i][j]);

			hFlowMass_Lambda[i][j] = (TH1D*) hFlowMass_Pos_Lambda[i][j]->Clone(Form("hFlowMass_Lambda_Cent%d_pt%d",i,j));
			hFlowMass_Lambda[i][j]->Add(hFlowMass_Neg_Lambda[i][j]);
			hFlowMass_Lambda[i][j]->Divide(fUnity,2);
			hFlowMass_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #it{v}_{2} %s %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",sEtaGap.Data(),fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hFlowMass_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassLambda/FlowMass_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lFlowMass_Lambda->Add(hFlowMass_Lambda[i][j]);
		}
	}

	// ===== Making comparison plots ======
	// inv mass plots
	TCanvas* cInvMassK0sPt = new TCanvas("cInvMassK0sPt","InvMassK0sPt",1800,1800);
	cInvMassK0sPt->Divide(3,3);
	TCanvas* cInvMassLambdaPt = new TCanvas("cInvMassLambdaPt","InvMassLambdaPt",1800,1800);
	cInvMassLambdaPt->Divide(3,3);
	TCanvas* cInvMassK0sCent = new TCanvas("cInvMassK0sCent","InvMassK0sCent",1800,1800);
	cInvMassK0sCent->Divide(3,3);
	TCanvas* cInvMassLambdaCent = new TCanvas("cInvMassLambdaCent","InvMassLambdaCent",1800,1800);
	cInvMassLambdaCent->Divide(3,3);

	// flow mass plots
	TCanvas* cFlowMassK0sPt = new TCanvas("cFlowMassK0sPt","FlowMassK0sPt",1800,1800);
	cFlowMassK0sPt->Divide(3,3);
	TCanvas* cFlowMassLambdaPt = new TCanvas("cFlowMassLambdaPt","FlowMassLambdaPt",1800,1800);
	cFlowMassLambdaPt->Divide(3,3);
	TCanvas* cFlowMassK0sCent = new TCanvas("cFlowMassK0sCent","FlowMassK0sCent",1800,1800);
	cFlowMassK0sCent->Divide(3,3);
	TCanvas* cFlowMassLambdaCent = new TCanvas("cFlowMassLambdaCent","FlowMassLambdaCent",1800,1800);
	cFlowMassLambdaCent->Divide(3,3);

	// loop for all cent and fixed pT
	for(Int_t j(0); j < iNumPtBins; j++) // pT
	{
		for(Int_t i(0); i < 9; i++) // cent
		{
			cInvMassK0sPt->cd(i+1);
			hInvMass_K0s[i][j]->Draw();

			cInvMassLambdaPt->cd(i+1);
			hInvMass_Lambda[i][j]->Draw();

			cFlowMassK0sPt->cd(i+1);
			hFlowMass_K0s[i][j]->Draw();
		
			cFlowMassLambdaPt->cd(i+1);
			hFlowMass_Lambda[i][j]->Draw();
		}
		cInvMassK0sPt->Print(Form("%s/compInvMass/InvMass_K0s_allCent_pt%d.%s",sOutput.Data(),j,sOutputFormat.Data()),sOutputFormat.Data());
		cInvMassLambdaPt->Print(Form("%s/compInvMass/InvMass_Lambda_allCent_pt%d.%s",sOutput.Data(),j,sOutputFormat.Data()),sOutputFormat.Data());
		cFlowMassK0sPt->Print(Form("%s/compFlowMass/FlowMass_K0s_allCent_pt%d.%s",sOutput.Data(),j,sOutputFormat.Data()),sOutputFormat.Data());
		cFlowMassLambdaPt->Print(Form("%s/compFlowMass/FlowMass_Lambda_allCent_pt%d.%s",sOutput.Data(),j,sOutputFormat.Data()),sOutputFormat.Data());
	}

	// loop for all pT and fixed cent
	for(Int_t i(0); i < iNumCentBins; i++) // cent
	{
	for(Int_t j(0); j < 9; j++) // pT
		{
			cInvMassK0sCent->cd(j+1);
			hInvMass_K0s[i][j]->Draw();

			cInvMassLambdaCent->cd(j+1);
			hInvMass_Lambda[i][j]->Draw();

			cFlowMassK0sCent->cd(j+1);
			hFlowMass_K0s[i][j]->Draw();

			cFlowMassLambdaCent->cd(j+1);
			hFlowMass_Lambda[i][j]->Draw();
		}
		cInvMassK0sCent->Print(Form("%s/compInvMass/InvMass_K0s_allPt_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
		cInvMassLambdaCent->Print(Form("%s/compInvMass/InvMass_Lambda_allPt_Cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
		cFlowMassK0sCent->Print(Form("%s/compFlowMass/FlowMass_K0s_allPt_cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
		cFlowMassLambdaCent->Print(Form("%s/compFlowMass/FlowMass_Lambda_allPt_Cent%d.%s",sOutput.Data(),i,sOutputFormat.Data()),sOutputFormat.Data());
	}


	

	// ===== Saving output ======
	fOutput->cd();
	lInvMass_K0s->Write("lInvMass_K0s",TObject::kSingleKey);
	lInvMass_Lambda->Write("lInvMass_Lambda",TObject::kSingleKey);
	lFlowMass_K0s->Write("lFlowMass_K0s",TObject::kSingleKey);
	lFlowMass_Lambda->Write("lFlowMass_Lambda",TObject::kSingleKey);
	lCummMass_Pos_K0s->Write("lCummMass_Pos_K0s",TObject::kSingleKey);
	lCummMass_Neg_K0s->Write("lCummMass_Neg_K0s",TObject::kSingleKey);
	lCummMass_Pos_Lambda->Write("lCummMass_Pos_Lambda",TObject::kSingleKey);
	lCummMass_Neg_Lambda->Write("lCummMass_Neg_Lambda",TObject::kSingleKey);
	lFlowMass_Pos_K0s->Write("lFlowMass_Pos_K0s",TObject::kSingleKey);
	lFlowMass_Neg_K0s->Write("lFlowMass_Neg_K0s",TObject::kSingleKey);
	lFlowMass_Pos_Lambda->Write("lFlowMass_Pos_Lambda",TObject::kSingleKey);
	lFlowMass_Neg_Lambda->Write("lFlowMass_Neg_Lambda",TObject::kSingleKey);
}