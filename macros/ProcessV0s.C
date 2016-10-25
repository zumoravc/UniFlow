void ProcessV0s(
		const TString sInput = "~/NBI/Flow/flow/AnalysisResults.root",
		const TString sOutput = "~/NBI/Flow/results/test",
		const TString sTag = "_JHEP",
		const TString sEtaGap = "Gap00",
		const TString sOutputFormat = "png"
	)
{
	//const TString sInput = "~/NBI/Codes/results/V0s/5/plusplus/merge/AnalysisResults_merged.root";
	//const TString sOutput = "~/NBI/Codes/results/V0s/5/plusplus/plots";
	//const TString sOutputFormat = "png";
	//const TString sEtaGap = "Gap09";
 
	const Int_t iNumPtBins = 22; // pT bins
	const Int_t iNumCentBins = 9; // centrality bins
	const Int_t iHarm = 2; // harmonics
	// bins edges
	Double_t fPtBinEdges[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4}; 
	Double_t fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

	// =======================================
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	TFile* fOutput = new TFile(Form("%s/../MassDist_V0s_%s.root",sOutput.Data(),sEtaGap.Data()),"RECREATE");

	fInput->cd(Form("%s",sTag.Data()));

	// ===== Loading input ===== 
	TList* lInputTracks = (TList*) gDirectory->Get(Form("Tracks_%s",sTag.Data()));
	TList* lInputV0s = (TList*) gDirectory->Get(Form("V0s_%s",sTag.Data()));
	//TList* lInputTracks = (TList*) gDirectory->Get(Form("Tracks",sTag.Data()));
	//TList* lInputV0s = (TList*) gDirectory->Get(Form("V0s"));

	// reference 
	TProfile* pTracksRefTwo = (TProfile*) (lInputTracks->FindObject(Form("fTracksRefTwo_n%d_%s_sample0",iHarm,sEtaGap.Data())) )->Clone(Form("pTracksRefTwo_%s",sEtaGap.Data())); 
	TH1D* hTracksRefTwo = (TH1D*) pTracksRefTwo->ProjectionX()->Clone(Form("hTracksRefTwo_%s",sEtaGap.Data()));
	
	// V0s 
	TProfile2D* p2V0sDiffTwoPos_K0s[iNumCentBins];
	TProfile2D* p2V0sDiffTwoNeg_K0s[iNumCentBins];
	TProfile2D* p2V0sDiffTwoPos_Lambda[iNumCentBins];
	TProfile2D* p2V0sDiffTwoNeg_Lambda[iNumCentBins];

	TH2D* h2InvMass_K0s[iNumCentBins];
	TH2D* h2InvMass_Lambda[iNumCentBins];
	
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		h2InvMass_K0s[i] = (TH2D*) (lInputV0s->FindObject(Form("fV0sPtInvMass_K0s_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i)))->Clone(Form("h2InvMass_K0s_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i));
		h2InvMass_Lambda[i] = (TH2D*) (lInputV0s->FindObject(Form("fV0sPtInvMass_Lambda_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i)))->Clone(Form("h2InvMass_Lambda_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i));

		p2V0sDiffTwoPos_K0s[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwoPos_K0s_n%d_%s_Cent%d_sample0",iHarm,sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwoPos_K0s_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i));
		p2V0sDiffTwoNeg_K0s[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwoNeg_K0s_n%d_%s_Cent%d_sample0",iHarm,sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwoNeg_K0s_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i));
	
		p2V0sDiffTwoPos_Lambda[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwoPos_Lambda_n%d_%s_Cent%d_sample0",iHarm,sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwoPos_Lambda_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i));
		p2V0sDiffTwoNeg_Lambda[i] = (TProfile2D*) (lInputV0s->FindObject(Form("fV0sDiffTwoNeg_Lambda_n%d_%s_Cent%d_sample0",iHarm,sEtaGap.Data(),i)))->Clone(Form("p2V0sDiffTwoNeg_Lambda_n%d_%s_Cent%d",iHarm,sEtaGap.Data(),i));
	}

	// ===== Making projections ======
	TH1D* hInvMass_K0s[iNumCentBins][iNumPtBins];
	TH1D* hInvMass_Lambda[iNumCentBins][iNumPtBins];

	TH1D* hCumMass_Pos_K0s[iNumCentBins][iNumPtBins];
	TH1D* hCumMass_Neg_K0s[iNumCentBins][iNumPtBins];
	TH1D* hCumMass_Pos_Lambda[iNumCentBins][iNumPtBins];
	TH1D* hCumMass_Neg_Lambda[iNumCentBins][iNumPtBins];

	TList* lInvMass_K0s = new TList();
	TList* lInvMass_Lambda = new TList();
	TList* lCumMass_Pos_K0s = new TList();
	TList* lCumMass_Neg_K0s = new TList();
	TList* lCumMass_Pos_Lambda = new TList();
	TList* lCumMass_Neg_Lambda = new TList();
	

	TCanvas* cTemp = new TCanvas();
	cTemp->cd();
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			printf("Cent %d pT %d\n",i,j);

			// inv mass plots
			hInvMass_K0s[i][j] = (TH1D*) h2InvMass_K0s[i]->ProjectionY(Form("hInvMass_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hInvMass_K0s[i][j]->SetTitle(Form("K_{S}^{0}: InvMass |#it{#eta}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hInvMass_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/InvMassK0s/InvMass_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lInvMass_K0s->Add(hInvMass_K0s[i][j]);

			hInvMass_Lambda[i][j] = (TH1D*) h2InvMass_Lambda[i]->ProjectionY(Form("hInvMass_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hInvMass_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: InvMass |#it{#eta}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hInvMass_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/InvMassLambda/InvMass_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lInvMass_Lambda->Add(hInvMass_Lambda[i][j]);
			
			// cum x inv mass plots
			hCumMass_Pos_K0s[i][j] = (TH1D*) p2V0sDiffTwoPos_K0s[i]->ProjectionY(Form("hCumMass_Pos_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCumMass_Pos_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCumMass_Pos_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/CumMassK0s/CumMass_Pos_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCumMass_Pos_K0s->Add(hCumMass_Pos_K0s[i][j]);
			
			hCumMass_Neg_K0s[i][j] = (TH1D*) p2V0sDiffTwoNeg_K0s[i]->ProjectionY(Form("hCumMass_Neg_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCumMass_Neg_K0s[i][j]->SetTitle(Form("K_{S}^{0}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|<-0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCumMass_Neg_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/CumMassK0s/CumMass_Neg_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCumMass_Neg_K0s->Add(hCumMass_Neg_K0s[i][j]);

			hCumMass_Pos_Lambda[i][j] = (TH1D*) p2V0sDiffTwoPos_Lambda[i]->ProjectionY(Form("hCumMass_Pos_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCumMass_Pos_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|>0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCumMass_Pos_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/CumMassLambda/CumMass_Pos_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCumMass_Pos_Lambda->Add(hCumMass_Pos_Lambda[i][j]);
			
			hCumMass_Neg_Lambda[i][j] = (TH1D*) p2V0sDiffTwoNeg_Lambda[i]->ProjectionY(Form("hCumMass_Neg_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCumMass_Neg_Lambda[i][j]->SetTitle(Form("#Lambda+#bar{#Lambda}: #LT#LT2'#GT#GT |#it{#eta}^{POI}|<-0.45 %g<#it{p}_{T}<%g GeV/#it{c} Cent %g-%g%%",fPtBinEdges[j],fPtBinEdges[j+1],fCentBinEdges[i],fCentBinEdges[i+1]));
			hCumMass_Neg_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/CumMassLambda/CumMass_Neg_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			lCumMass_Neg_Lambda->Add(hCumMass_Neg_Lambda[i][j]);

			
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

	TH1D* hRefFlow = (TH1D*) hTracksRefTwo->Clone("hRefFlowTwo");

	TList* lFlowMass_K0s = new TList();
	TList* lFlowMass_Pos_K0s = new TList();
	TList* lFlowMass_Neg_K0s = new TList();
	TList* lFlowMass_Lambda = new TList();
	TList* lFlowMass_Pos_Lambda = new TList();
	TList* lFlowMass_Neg_Lambda = new TList();
	
	Double_t dRefFlow = 0;

	for(Int_t i(0); i < iNumCentBins; i++)
	{
		dRefFlow = TMath::Sqrt(hTracksRefTwo->GetBinContent(i+1));
		hRefFlow->SetBinContent(i+1, dRefFlow);

		for(Int_t j(0); j < iNumPtBins; j++)
		{
			hFlowMass_Pos_K0s[i][j] = (TH1D*) hCumMass_Pos_K0s[i][j]->Clone(Form("hFlowMass_Pos_K0s_Cent%d_pt%d",i,j));
			hFlowMass_Neg_K0s[i][j] = (TH1D*) hCumMass_Neg_K0s[i][j]->Clone(Form("hFlowMass_Neg_K0s_Cent%d_pt%d",i,j));
			hFlowMass_Pos_Lambda[i][j] = (TH1D*) hCumMass_Pos_Lambda[i][j]->Clone(Form("hFlowMass_Pos_Lambda_Cent%d_pt%d",i,j));
			hFlowMass_Neg_Lambda[i][j] = (TH1D*) hCumMass_Neg_Lambda[i][j]->Clone(Form("hFlowMass_Neg_Lambda_Cent%d_pt%d",i,j));

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
		hRefFlow->Draw();
		cTemp->Print(Form("%s/RefFlow/RefFlowTwo.%s",sOutput.Data(),sOutputFormat.Data()),sOutputFormat.Data());	
	}

	// ===== Saving output ======
	fOutput->cd();
	hRefFlow->Write();
	lInvMass_K0s->Write("lInvMass_K0s",TObject::kSingleKey);
	lInvMass_Lambda->Write("lInvMass_Lambda",TObject::kSingleKey);
	lFlowMass_K0s->Write("lFlowMass_K0s",TObject::kSingleKey);
	lFlowMass_Lambda->Write("lFlowMass_Lambda",TObject::kSingleKey);
	lCumMass_Pos_K0s->Write("lCumMass_Pos_K0s",TObject::kSingleKey);
	lCumMass_Neg_K0s->Write("lCumMass_Neg_K0s",TObject::kSingleKey);
	lCumMass_Pos_Lambda->Write("lCumMass_Pos_Lambda",TObject::kSingleKey);
	lCumMass_Neg_Lambda->Write("lCumMass_Neg_Lambda",TObject::kSingleKey);
	lFlowMass_Pos_K0s->Write("lFlowMass_Pos_K0s",TObject::kSingleKey);
	lFlowMass_Neg_K0s->Write("lFlowMass_Neg_K0s",TObject::kSingleKey);
	lFlowMass_Pos_Lambda->Write("lFlowMass_Pos_Lambda",TObject::kSingleKey);
	lFlowMass_Neg_Lambda->Write("lFlowMass_Neg_Lambda",TObject::kSingleKey);
}