void CompareK0s()
{	
	TString sPath = "~/NBI/Flow/results/V0s/15-sampling/sampling_test2/";
	TString sGap = "Gap10";

	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatios.C");
	gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	// loading TGraphs QM12 
	TFile* fInK0s1020 = new TFile("~/NBI/Flow/hepData/KsQM_1020.root","READ");
	TFile* fInK0s4050 = new TFile("~/NBI/Flow/hepData/KsQM_4050.root","READ");

	TH1D* hQMK0s[2];
	fInK0s1020->cd();
	hQMK0s[0] = (TH1D*) gDirectory->Get("hist")->Clone("hQMK0s_1020");

	fInK0s4050->cd();
	hQMK0s[1] = (TH1D*) gDirectory->Get("hist")->Clone("hQMK0s_4050");
	
	TFile* fInLambda1020 = new TFile("~/NBI/Flow/hepData/LambdaQM_1020.root","READ");
	TFile* fInLambda4050 = new TFile("~/NBI/Flow/hepData/LambdaQM_4050.root","READ");

	TH1D* hQMLambda[2];
	fInLambda1020->cd();
	hQMLambda[0] = (TH1D*) gDirectory->Get("hist")->Clone("hQMLambda_1020");

	fInLambda4050->cd();
	hQMLambda[1] = (TH1D*) gDirectory->Get("hist")->Clone("hQMLambda_4050");
	printf("===== QM prel loaded =====\n");

	// loading my input
	printf("5 bins\n");
	TFile* fInMine5bins = new TFile(Form("~/NBI/Flow/results/V0s/15-sampling-5bins/sampling_test/PtFlow_V0s_%s.root",sGap.Data()),"READ");
	fInMine5bins->cd();
	gDirectory->ls();

	TH1D* hK0s_5bins[2];
	TH1D* hLambda_5bins[2];
	hK0s_5bins[0] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent2")->Clone("hK0s_5bins_1020");
	hK0s_5bins[1] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent5")->Clone("hK0s_5bins_4050");
	hLambda_5bins[0] = (TH1D*) gDirectory->Get("hFlowPt_Lambda_Gap10_Cent2")->Clone("hLambda_5bins_1020");
	hLambda_5bins[1] = (TH1D*) gDirectory->Get("hFlowPt_Lambda_Gap10_Cent5")->Clone("hLambda_5bins_4050");
	
	
	printf("10 bins\n");
	TFile* fInMine10bins = new TFile(Form("~/NBI/Flow/results/V0s/15-sampling/testGap10/sampling/PtFlow_V0s_%s.root",sGap.Data()),"READ");
	fInMine10bins->cd();
	gDirectory->ls();
	
	TH1D* hK0s_10bins[2];
	TH1D* hLambda_10bins[2];
	hK0s_10bins[0] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent2")->Clone("hK0s_10bins_1020");
	hK0s_10bins[1] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent5")->Clone("hK0s_10bins_4050");
	hLambda_10bins[0] = (TH1D*) gDirectory->Get("hFlowPt_Lambda_Gap10_Cent2")->Clone("hLambda_10bins_1020");
	hLambda_10bins[1] = (TH1D*) gDirectory->Get("hFlowPt_Lambda_Gap10_Cent5")->Clone("hLambda_10bins_4050");
	
	printf("No bins\n");
	TFile* fInMineNoBins = new TFile(Form("~/NBI/Flow/results/V0s/15-sampling-5bins/noSampling_test/PtFlow_V0s_%s.root",sGap.Data()),"READ");
	fInMineNoBins->cd();
	gDirectory->ls();

	TH1D* hK0s_Nobins[2];
	TH1D* hLambda_Nobins[2];
	hK0s_Nobins[0] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent2")->Clone("hK0s_Nobins_1020");
	hK0s_Nobins[1] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent5")->Clone("hK0s_NObins_4050");
	hLambda_Nobins[0] = (TH1D*) gDirectory->Get("hFlowPt_Lambda_Gap10_Cent2")->Clone("hLambda_Nobins_1020");
	hLambda_Nobins[1] = (TH1D*) gDirectory->Get("hFlowPt_Lambda_Gap10_Cent5")->Clone("hLambda_NObins_4050");

	printf("No bins2\n");
	TFile* fInMineNoBins2 = new TFile(Form("~/NBI/Flow/results/V0s/15-sampling/noSampling_test/PtFlow_V0s_%s.root",sGap.Data()),"READ");
	fInMineNoBins2->cd();
	gDirectory->ls();

	TH1D* hK0s_Nobins2[2];
	hK0s_Nobins2[0] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent2")->Clone("hK0s_Nobins2_1020");
	hK0s_Nobins2[1] = (TH1D*) gDirectory->Get("hFlowPt_K0s_Gap10_Cent5")->Clone("hK0s_NObins2_4050");

	
	printf("===== My results loaded =====\n");
	Color_t colQM = kRed;
	Color_t colMine = kBlue;
	Short_t styleK0s = 20;
	Short_t styleK0sOpen = 24;
	Short_t styleLambda = 20;
	Short_t styleLambdaOpen = 24;

	TLatex latex = TLatex();
	latex.SetNDC();
	TLegend* leg = new TLegend(0.33,0.7,0.6,0.89);
	leg->SetBorderSize(0);
	leg->AddEntry(hQMK0s[0],"K^{0}_{S} (QM12 Prel)","pel"); 
	leg->AddEntry(hQMLambda[0],"#Lambda (QM12 Prel)","pel"); 

	TLegend* legMine = new TLegend(0.15,0.7,0.3,0.89);
	legMine->SetBorderSize(0);
	legMine->AddEntry(hK0s_Nobins[0],"K^{0}_{S} (this)","pel"); 
	legMine->AddEntry(hLambda_Nobins[0],"#Lambda (this)","pel"); 

	TCanvas* canV0s = new TCanvas("canV0s","canV0s",1000,500);
	canV0s->Divide(2,1);
	for(Short_t i(0); i < 2; i++)
	{
	canV0s->cd(i+1);
	hQMK0s[i]->SetStats(0);
	hQMK0s[i]->SetTitle("; #it{p}_{T} (GeV/#it{c}); #it{v}_{2} {2, |#Delta#eta| > 1.0}");
	hQMK0s[i]->SetTitleOffset(1.5,"y");
	hQMK0s[i]->SetMinimum(-0.0);
	hQMK0s[i]->SetMaximum(0.35);
	hQMK0s[i]->SetLineColor(kRed);
	hQMK0s[i]->SetMarkerColor(kRed);
	hQMK0s[i]->SetMarkerStyle(styleK0sOpen);
	hQMK0s[i]->Draw();
	hK0s_Nobins[i]->SetLineColor(kRed);
	hK0s_Nobins[i]->SetMarkerColor(kRed);
	hK0s_Nobins[i]->SetMarkerStyle(styleK0s);
	hK0s_Nobins[i]->Draw("same");


	hQMLambda[i]->SetLineColor(kBlue);
	hQMLambda[i]->SetMarkerColor(kBlue);
	hQMLambda[i]->SetMarkerStyle(styleK0sOpen);
	hQMLambda[i]->Draw("same");
	hLambda_Nobins[i]->SetLineColor(kBlue);
	hLambda_Nobins[i]->SetMarkerColor(kBlue);
	hLambda_Nobins[i]->SetMarkerStyle(styleK0s);
	hLambda_Nobins[i]->Draw("same");
}
	canV0s->cd(1);
	leg->Draw();
	legMine->Draw("same");
	latex.DrawLatex(0.6,0.15,"Cent 10-20%");
	canV0s->cd(2);
	latex.DrawLatex(0.6,0.15,"Cent 40-50%");
	
	canV0s->SaveAs(Form("%s/V0sComp.pdf",sPath.Data()));
	return;

	//hQMK0s[0]->SetMinimum(0.);
	//hQMK0s[0]->SetMaximum(0.3);
	hQMK0s[1]->SetLineColor(kBlue);
	hQMK0s[1]->SetMarkerColor(kBlue);
	hQMK0s[1]->SetMarkerStyle(styleK0sOpen);
	hQMK0s[1]->Draw("same");
	hK0s_Nobins[1]->SetLineColor(kBlue);
	hK0s_Nobins[1]->SetMarkerColor(kBlue);
	hK0s_Nobins[1]->SetMarkerStyle(styleK0s);
	hK0s_Nobins[1]->Draw("same");

	//hK0s_Nobins2[0]->SetLineColor(kGreen);
	//hK0s_Nobins2[0]->Draw("same");
	//hK0s_5bins[0]->SetLineColor(kMagenta);
	//hK0s_5bins[0]->Draw("same");
	//hK0s_10bins[0]->SetLineColor(kPink);
	//hK0s_10bins[0]->Draw("same");

	canV0s->cd(2);
	hQMLambda[0]->SetLineColor(colQM);
	hQMLambda[0]->SetMarkerColor(colQM);
	hQMLambda[0]->SetMarkerStyle(styleLambda);
	hQMLambda[0]->Draw("");
	hLambda_Nobins[0]->SetLineColor(colMine);
	hLambda_Nobins[0]->SetMarkerColor(colMine);
	hLambda_Nobins[0]->SetMarkerStyle(styleLambdaOpen);
	hLambda_Nobins[0]->Draw("same");

	//hLambda_10bins[0]->Draw("same");
	//hLambda_Nobins[0]->Draw("same");

  return;

	TList* lCompK0s[2];
	TList* lCompLambda[2];
	TString sLabel[] = {"QM","10 samples","5 samples","no samples"};
	TString sLabelLambda[] = {"10 samples","5 samples","no samples","QM"};
	TString sCent[2] = {"Cent 10-20%", "Cent 40-50%"}
	for(Int_t i(0); i < 2; i++)
	{
		hQMK0s[i]->SetTitle(Form("v_{2} {2,|#Delta#it{#eta}| > 1} %s; #it{p}_{T} (GeV/#it{c}); v_{2} {2,|#Delta#it{#eta}| > 1}",sCent[i].Data()));
		hK0s_10bins[i]->SetTitle(Form("v_{2} {2,|#Delta#it{#eta}| > 1} %s; #it{p}_{T} (GeV/#it{c}); v_{2} {2,|#Delta#it{#eta}| > 1}",sCent[i].Data()));
		hQMK0s[i]->SetStats(0);
		hQMK0s[i]->SetMinimum(0);
		hQMK0s[i]->SetMaximum(0.3);

		lCompK0s[i] = new TList();
		lCompK0s[i]->Add(hQMK0s[i]);
		lCompK0s[i]->Add(hK0s_10bins[i]);
		lCompK0s[i]->Add(hK0s_5bins[i]);
		lCompK0s[i]->Add(hK0s_Nobins[i]);
		lCompK0s[i]->Add(hK0s_Nobins2[i]);
		
		hLambda_10bins[i]->SetTitle(Form("v_{2} {2,|#Delta#it{#eta}| > 1} %s; #it{p}_{T} (GeV/#it{c}); v_{2} {2,|#Delta#it{#eta}| > 1}",sCent[i].Data()));
		hQMLambda[i]->SetTitle(Form("v_{2} {2,|#Delta#it{#eta}| > 1} %s; #it{p}_{T} (GeV/#it{c}); v_{2} {2,|#Delta#it{#eta}| > 1}",sCent[i].Data()));
		hLambda_10bins[i]->SetStats(0);
		hLambda_10bins[i]->SetMinimum(0);
		hLambda_10bins[i]->SetMaximum(0.3);
		
		lCompLambda[i] = new TList();
		lCompLambda[i]->Add(hLambda_10bins[i]);
		lCompLambda[i]->Add(hLambda_5bins[i]);
		lCompLambda[i]->Add(hLambda_Nobins[i]);
		lCompLambda[i]->Add(hQMLambda[i]);
	}

	// 

	TCanvas* cCompRatio[2];
	TCanvas* cCompRatioLambda[2];
		//cCompRatio[1] = CompareRatio(hK0s[1],hQMK0s[1]);

	for(Int_t i(0); i < 2; i++)
	{
		cCompRatio[i] = CompareHistos(lCompK0s[i],sLabel,0.,0.3,kFALSE);
		cCompRatio[i] = CompareRatios(4,lCompK0s[i],sLabel,0.,0.3);
		cCompRatio[i]->Print(Form("%s/comp_K0s_%d.pdf",sPath.Data(),i),"pdf");
		cCompRatioLambda[i] = CompareHistos(lCompLambda[i],sLabelLambda,-0.1,0.4,kFALSE);
		cCompRatioLambda[i]->Print(Form("%s/comp_Lambda_%d.pdf",sPath.Data(),i),"pdf");
	}
}