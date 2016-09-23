void ProcessV0s()
{
	TString sInput = "~/NBI/Codes/flow/AnalysisResults.root";
	TString sOutput = "~/NBI/Codes/results/V0s/3";
	TString sOutputFormat = "png";
 
  	// =======================================
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareRatio.C");
	gROOT->LoadMacro("~/NBI/Codes/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	TFile* fOutput = new TFile(Form("%s/V0sFlow.root",sOutput.Data()),"RECREATE");

	fInput->cd("FlowPID");

	const Int_t iNumCentBins = 9; // centrality bins
	const Int_t iNumPtBins = 10; // pT bins

	// ===== Loading input ===== 
	TList* lInputTracks = (TList*) gDirectory->Get("Tracks");
	TList* lInputV0s = (TList*) gDirectory->Get("V0s");

	// reference 
	TProfile* pRefCorTwo2_Gap09 = (TProfile*) lInputTracks->FindObject("fRefCorTwo2_Gap09")->Clone("pRefCorTwo2_Gap09"); 
	TH1D* hRefCorTwo2_Gap09 = (TH1D*) pRefCorTwo2_Gap09->ProjectionX()->Clone("hRefCorTwo2_Gap09");
	
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
		p2V0sDiffTwo2_Gap09N_Lambda[i] = (TProfile2D*) lInputV0s->FindObject(Form("fV0sDiffTwo2_Gap09N_Lambda_Cent%d",i))->Clone(Form("p2V0sDiffTwo2_Gap09N_Lambda_Cent%d",i));
	}

	// ===== Making projections ======
	TH1D* hInvMass_K0s[iNumCentBins][iNumPtBins];
	TH1D* hInvMass_Lambda[iNumCentBins][iNumPtBins];

	TH1D* hCummMass_Pos_K0s[iNumCentBins][iNumPtBins];
	TH1D* hCummMass_Neg_K0s[iNumCentBins][iNumPtBins];
	TH1D* hCummMass_Pos_Lambda[iNumCentBins][iNumPtBins];
	TH1D* hCummMass_Neg_Lambda[iNumCentBins][iNumPtBins];

	TCanvas* cTemp = new TCanvas();
	cTemp->cd();
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			// inv mass plots
			hInvMass_K0s[i][j] = (TH1D*) h2InvMass_Gap09_K0s[i]->ProjectionY(Form("hInvMass_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hInvMass_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/InvMassK0s/InvMass_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());

			hInvMass_Lambda[i][j] = (TH1D*) h2InvMass_Gap09_Lambda[i]->ProjectionY(Form("hInvMass_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hInvMass_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/InvMassLambda/InvMass_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			
			// cum x inv mass plots
			hCummMass_Pos_K0s[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09P_K0s[i]->ProjectionY(Form("hCummMass_Pos_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Pos_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassK0s/CummMass_Pos_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			
			hCummMass_Neg_K0s[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09N_K0s[i]->ProjectionY(Form("hCummMass_Neg_K0s_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Neg_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassK0s/CummMass_Neg_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());

			hCummMass_Pos_Lambda[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09P_Lambda[i]->ProjectionY(Form("hCummMass_Pos_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Pos_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassLambda/CummMass_Pos_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			
			hCummMass_Neg_Lambda[i][j] = (TH1D*) p2V0sDiffTwo2_Gap09N_Lambda[i]->ProjectionY(Form("hCummMass_Neg_Lambda_Cent%d_pt%d",i,j),j+1,j+1,"e");
			hCummMass_Neg_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/CummMassLambda/CummMass_Neg_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
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

			hFlowMass_Pos_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassK0s/FlowMass_Pos_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			
			hFlowMass_Neg_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassK0s/FlowMass_Neg_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			
			hFlowMass_Pos_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassLambda/FlowMass_Pos_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			
			hFlowMass_Neg_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassLambda/FlowMass_Neg_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());

			hFlowMass_K0s[i][j] = (TH1D*) hFlowMass_Pos_K0s[i][j]->Clone(Form("hFlowMass_K0s_Cent%d_pt%d",i,j));
			hFlowMass_K0s[i][j]->Add(hFlowMass_Neg_K0s[i][j]);
			hFlowMass_K0s[i][j]->Divide(fUnity,2);
			hFlowMass_K0s[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassK0s/FlowMass_K0s_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());

			hFlowMass_Lambda[i][j] = (TH1D*) hFlowMass_Pos_Lambda[i][j]->Clone(Form("hFlowMass_Lambda_Cent%d_pt%d",i,j));
			hFlowMass_Lambda[i][j]->Add(hFlowMass_Neg_Lambda[i][j]);
			hFlowMass_Lambda[i][j]->Divide(fUnity,2);
			hFlowMass_Lambda[i][j]->Draw();
			cTemp->Print(Form("%s/FlowMassLambda/FlowMass_Lambda_Cent%d_pt%d.%s",sOutput.Data(),i,j,sOutputFormat.Data()),sOutputFormat.Data());
			

		}
	}

	// ===== Saving output ======
	fOutput->cd();

	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			hInvMass_K0s[i][j]->Write();
			hInvMass_Lambda[i][j]->Write();

			hCummMass_Pos_K0s[i][j]->Write();
			hCummMass_Neg_K0s[i][j]->Write();
			hCummMass_Pos_Lambda[i][j]->Write();
			hCummMass_Neg_Lambda[i][j]->Write();

			hFlowMass_Pos_K0s[i][j]->Write();
			hFlowMass_Neg_K0s[i][j]->Write();
			hFlowMass_Pos_Lambda[i][j]->Write();
			hFlowMass_Neg_Lambda[i][j]->Write();

			hFlowMass_K0s[i][j]->Write();
			hFlowMass_Lambda[i][j]->Write();
		}
	}


}