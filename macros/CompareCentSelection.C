TH1D* ConvertHist(TH1D* hist, TH1D* stat, TH1D* sys);

void CompareCentSelection()
{
	TString sOutputPath = "~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/centComparison/";


	//TFile* fInputNew = new TFile("~/NBI/Flow/results/7-GFK-PbPb-OldCent/Flow.root","READ");
	
	//TFile* fInputOld = new TFile("~/NBI/Flow/results/7-GFK-PbPb-OldCent/Flow_OldCent.root","READ");

	TFile* fInputOld_PileON_PeriodON = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_Old_PileON_PeriodON.root","READ");
	TFile* fInputOld_PileOFF_PeriodON = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_Old_PileOFF_PeriodON.root","READ");
	TFile* fInputOld_PileON_PeriodOFF = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_Old_PileON_PeriodOFF.root","READ");
	TFile* fInputOld_PileOFF_PeriodOFF = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_Old_PileOFF_PeriodOFF.root","READ");
	
	TFile* fInputNew_PileON_PeriodON = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_New_PileON_PeriodON.root","READ");
	TFile* fInputNew_PileOFF_PeriodON = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_New_PileOFF_PeriodON.root","READ");
	TFile* fInputNew_PileON_PeriodOFF = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_New_PileON_PeriodOFF.root","READ");
	TFile* fInputNew_PileOFF_PeriodOFF = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_New_PileOFF_PeriodOFF.root","READ");
	
	TFile* fQA_Old = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/merge/AnalysisResults.root","READ");
	TFile* fQA_New = new TFile("~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/merge/AnalysisResults.root","READ");
	
	Color_t colNew = kRed;
	Color_t colOld = kGreen+2;
	Color_t colPRL = kBlack;

	Color_t col_ON_ON = kBlue;
	Color_t col_ON_OFF = kRed;
	Color_t col_OFF_ON = kGreen+2;
	Color_t col_OFF_OFF = kMagenta;

	TLatex latex = TLatex();
	//latex.UseNDC();

	/*
	printf("=================== fInputNew =====================\n");
	fInputNew->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_New = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	hRef_v22_Gap10_New->SetLineColor(colNew);
	hRef_v22_Gap10_New->SetMarkerColor(colNew);

	printf("=================== fInpuOld =====================\n");
	fInputOld->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_Old = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	hRef_v22_Gap10_Old->SetLineColor(colOld);
	hRef_v22_Gap10_Old->SetMarkerColor(colOld);
	*/

	printf("================ fInputOld_variations =================\n");
	fInputOld_PileON_PeriodON->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_Old_ON_ON= (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	//hRef_v22_Gap10_Old_ON_ON->SetLineColor(colOld_OFF_OFF);
	//hRef_v22_Gap10_Old_ON_ON->SetMarkerColor(colOld_OFF_OFF);

	fInputOld_PileOFF_PeriodOFF->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_Old_OFF_OFF = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	//hRef_v22_Gap10_Old_OFF_OFF->SetLineColor(colOld_OFF_OFF);
	//hRef_v22_Gap10_Old_OFF_OFF->SetMarkerColor(colOld_OFF_OFF);

	fInputOld_PileON_PeriodOFF->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_Old_ON_OFF = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	//hRef_v22_Gap10_Old_ON_OFF->SetLineColor(colOld_OFF_OFF);
	//hRef_v22_Gap10_Old_ON_OFF->SetMarkerColor(colOld_OFF_OFF);

	fInputOld_PileOFF_PeriodON->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_Old_OFF_ON = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	//hRef_v22_Gap10_Old_OFF_ON->SetLineColor(colOld_OFF_ON);
	//hRef_v22_Gap10_Old_OFF_ON->SetMarkerColor(colOld_OFF_ON);

	printf("================ fInputNew_variations =================\n");
	fInputNew_PileON_PeriodON->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_New_ON_ON= (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	
	fInputNew_PileOFF_PeriodOFF->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_New_OFF_OFF = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	
	fInputNew_PileON_PeriodOFF->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_New_ON_OFF = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	
	fInputNew_PileOFF_PeriodON->cd();
	gDirectory->ls();
	TH1D* hRef_v22_Gap10_New_OFF_ON = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	
	// loading PRL 116
	printf("=================== Loading PRL116 ===================\n");
	TFile* fInputPRL = new TFile("~/NBI/Flow/hepData/PRL116/HEPData-ins1419244-1-root.root","READ");
	fInputPRL->cd("Table 3"); // 2.76 PbPb
	gDirectory->ls();
	TH1D* hPRL116_v22_Gap10_orig = (TH1D*) gDirectory->Get("Hist1D_y1");
	hPRL116_v22_Gap10_orig->SetLineColor(colPRL);
	hPRL116_v22_Gap10_orig->SetMarkerColor(colPRL);
	TH1D* hPRL116_v22_Gap10_stat = (TH1D*) gDirectory->Get("Hist1D_y1_e1");
	TH1D* hPRL116_v22_Gap10_sys = (TH1D*) gDirectory->Get("Hist1D_y1_e2");
	TGraphAsymmErrors* gPRL116_v22_Gap10 = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1");
	gPRL116_v22_Gap10->SetLineColor(colPRL);
	gPRL116_v22_Gap10->SetMarkerColor(colPRL);
	TGraphAsymmErrors* gPRL116_v24 = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y2");
	gPRL116_v24->SetLineColor(colPRL);
	gPRL116_v24->SetMarkerColor(colPRL);
	TH1D* hPRL116_v24 = (TH1D*) gDirectory->Get("Hist1D_y2");
	hPRL116_v24->SetLineColor(colPRL);
	hPRL116_v24->SetMarkerColor(colPRL);
	TH1D* hPRL116_v22_Gap10 = ConvertHist(hPRL116_v22_Gap10_orig,hPRL116_v22_Gap10_stat,hPRL116_v22_Gap10_sys);
	printf("=================== Loading done =====================\n");
	

	hRef_v22_Gap10_Old_ON_ON->SetLineColor(colOld);
	hRef_v22_Gap10_Old_ON_ON->SetMarkerColor(colOld);
	hRef_v22_Gap10_New_ON_ON->SetLineColor(colNew);
	hRef_v22_Gap10_New_ON_ON->SetMarkerColor(colNew);

	TCanvas* cOldNew = new TCanvas("cOldNew","cOldNew");
	cOldNew->cd();
	hRef_v22_Gap10_Old_ON_ON->Draw();
	hRef_v22_Gap10_New_ON_ON->Draw("same");

	cOldNew->SaveAs(Form("%s/Cent_OLD_NEW.pdf",sOutputPath.Data()));


	printf("=================== Comparing OLD cent variations =====================\n");

	TLegend* legComp = new TLegend(0.5,0.13,0.89,0.4);
	legComp->SetBorderSize(0);
	legComp->AddEntry(hRef_v22_Gap10_Old_ON_ON,"PileUp ON && Period ON","pel");
	legComp->AddEntry(hRef_v22_Gap10_Old_ON_OFF,"PileUp ON && Period OFF","pel");
	legComp->AddEntry(hRef_v22_Gap10_Old_OFF_ON,"PileUp OFF && Period ON","pel");
	legComp->AddEntry(hRef_v22_Gap10_Old_OFF_OFF,"PileUp OFF && Period OFF","pel");
	legComp->AddEntry(gPRL116_v22_Gap10,"PRL 116","pel");

	// ratios
	TH1D* hRatio_Old_ON_ON = (TH1D*) hRef_v22_Gap10_Old_ON_ON->Clone(); // unity by def.
	hRatio_Old_ON_ON->Divide(hRef_v22_Gap10_Old_ON_ON);
	hRatio_Old_ON_ON->SetLineColor(col_ON_ON);
	hRatio_Old_ON_ON->SetMarkerColor(col_ON_ON);

	TH1D* hRatio_Old_ON_OFF = (TH1D*) hRef_v22_Gap10_Old_ON_ON->Clone();
	hRatio_Old_ON_OFF->Divide(hRef_v22_Gap10_Old_ON_OFF);
	hRatio_Old_ON_OFF->SetLineColor(col_ON_OFF);
	hRatio_Old_ON_OFF->SetMarkerColor(col_ON_OFF);

	TH1D* hRatio_Old_OFF_ON = (TH1D*) hRef_v22_Gap10_Old_ON_ON->Clone();
	hRatio_Old_OFF_ON->Divide(hRef_v22_Gap10_Old_OFF_ON);
	hRatio_Old_OFF_ON->SetLineColor(col_OFF_ON);
	hRatio_Old_OFF_ON->SetMarkerColor(col_OFF_ON);

	TH1D* hRatio_Old_OFF_OFF = (TH1D*) hRef_v22_Gap10_Old_ON_ON->Clone();
	hRatio_Old_OFF_OFF->Divide(hRef_v22_Gap10_Old_OFF_OFF);
	hRatio_Old_OFF_OFF->SetLineColor(col_OFF_OFF);
	hRatio_Old_OFF_OFF->SetMarkerColor(col_OFF_OFF);

	TH1D* hRatio_Old_PRL116 = (TH1D*) hRef_v22_Gap10_Old_ON_ON->Clone();
	hRatio_Old_PRL116->Divide(hPRL116_v22_Gap10);
	hRatio_Old_PRL116->SetLineColor(colPRL);
	hRatio_Old_PRL116->SetMarkerColor(colPRL);
	hRatio_Old_PRL116->SetMarkerStyle(4);

	//loading centrality dist
	//fQA_Old->ls();
	fQA_Old->cd("flowPID_FB768_Old_PileON_PeriodON");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_Old_PileON_PeriodON");
	TH1D* hCentDist_ON_ON = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_ON_ON->SetLineColor(col_ON_ON);
	hCentDist_ON_ON->SetMarkerColor(col_ON_ON);

	fQA_Old->cd("flowPID_FB768_Old_PileON_PeriodOFF");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_Old_PileON_PeriodOFF");
	TH1D* hCentDist_ON_OFF = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_ON_OFF->SetLineColor(col_ON_OFF);
	hCentDist_ON_OFF->SetMarkerColor(col_ON_OFF);

	fQA_Old->cd("flowPID_FB768_Old_PileOFF_PeriodOFF");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_Old_PileOFF_PeriodOFF");
	TH1D* hCentDist_OFF_OFF = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_OFF_OFF->SetLineColor(col_OFF_OFF);
	hCentDist_OFF_OFF->SetMarkerColor(col_OFF_OFF);

	fQA_Old->cd("flowPID_FB768_Old_PileOFF_PeriodON");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_Old_PileOFF_PeriodON");
	TH1D* hCentDist_OFF_ON = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_OFF_ON->SetLineColor(col_OFF_ON);
	hCentDist_OFF_ON->SetMarkerColor(col_OFF_ON);

	// plotting
	TCanvas* cOld = new TCanvas("cOld","Old",1000,1000);
	cOld->Divide(2,2);
	cOld->cd(1);
	hRef_v22_Gap10_Old_ON_ON->SetLineColor(col_ON_ON);
	hRef_v22_Gap10_Old_ON_ON->Draw();
	hRef_v22_Gap10_Old_ON_OFF->SetLineColor(col_ON_OFF);
	hRef_v22_Gap10_Old_ON_OFF->Draw("same");
	hRef_v22_Gap10_Old_OFF_ON->SetLineColor(col_OFF_ON);
	hRef_v22_Gap10_Old_OFF_ON->Draw("same");
	hRef_v22_Gap10_Old_OFF_OFF->SetLineColor(col_OFF_OFF);
	hRef_v22_Gap10_Old_OFF_OFF->Draw("same");
	hPRL116_v22_Gap10->Draw("same");
	legComp->SetHeader("Old centrality framework (AliCentrality)");
	legComp->Draw();

	cOld->cd(2);
	hRatio_Old_ON_ON->SetMinimum(0.92);
	hRatio_Old_ON_ON->SetMaximum(1.08);
	hRatio_Old_ON_ON->Draw();
	hRatio_Old_ON_OFF->Draw("same");
	hRatio_Old_OFF_ON->Draw("same");
	hRatio_Old_OFF_OFF->Draw("same");
	hRatio_Old_PRL116->Draw("same");


	cOld->cd(3);
	hCentDist_ON_ON->Draw();
	hCentDist_ON_OFF->Draw("same");
	hCentDist_OFF_ON->Draw("same");
	hCentDist_OFF_OFF->Draw("same");

	cOld->SaveAs(Form("%s/CentOLD.pdf",sOutputPath.Data()));
	
	cOld->cd(3);
	//hCentDist_ON_ON->Draw();
	//hCentDist_ON_OFF->Draw("same");
	//hCentDist_OFF_ON->Draw("same");
	hCentDist_OFF_OFF->Draw();
	cOld->SaveAs(Form("%s/CentOLD_one.pdf",sOutputPath.Data()));

	printf("=================== Comparing NEW cent variations =====================\n");
	
	// ratios
	TH1D* hRatio_New_ON_ON = (TH1D*) hRef_v22_Gap10_New_ON_ON->Clone(); // unity by def.
	hRatio_New_ON_ON->Divide(hRef_v22_Gap10_New_ON_ON);
	hRatio_New_ON_ON->SetLineColor(col_ON_ON);
	hRatio_New_ON_ON->SetMarkerColor(col_ON_ON);

	TH1D* hRatio_New_ON_OFF = (TH1D*) hRef_v22_Gap10_New_ON_ON->Clone();
	hRatio_New_ON_OFF->Divide(hRef_v22_Gap10_New_ON_OFF);
	hRatio_New_ON_OFF->SetLineColor(col_ON_OFF);
	hRatio_New_ON_OFF->SetMarkerColor(col_ON_OFF);

	TH1D* hRatio_New_OFF_ON = (TH1D*) hRef_v22_Gap10_New_ON_ON->Clone();
	hRatio_New_OFF_ON->Divide(hRef_v22_Gap10_New_OFF_ON);
	hRatio_New_OFF_ON->SetLineColor(col_OFF_ON);
	hRatio_New_OFF_ON->SetMarkerColor(col_OFF_ON);

	TH1D* hRatio_New_OFF_OFF = (TH1D*) hRef_v22_Gap10_New_ON_ON->Clone();
	hRatio_New_OFF_OFF->Divide(hRef_v22_Gap10_New_OFF_OFF);
	hRatio_New_OFF_OFF->SetLineColor(col_OFF_OFF);
	hRatio_New_OFF_OFF->SetMarkerColor(col_OFF_OFF);

	TH1D* hRatio_New_PRL116 = (TH1D*) hRef_v22_Gap10_New_ON_ON->Clone();
	hRatio_New_PRL116->Divide(hPRL116_v22_Gap10);
	hRatio_New_PRL116->SetLineColor(colPRL);
	hRatio_New_PRL116->SetMarkerColor(colPRL);
	hRatio_New_PRL116->SetMarkerStyle(4);

	//loading centrality dist
	//fQA_New->ls();
	fQA_New->cd("flowPID_FB768_New_PileON_PeriodON");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_New_PileON_PeriodON");
	TH1D* hCentDist_ON_ON = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_ON_ON->SetLineColor(col_ON_ON);
	hCentDist_ON_ON->SetMarkerColor(col_ON_ON);

	fQA_New->cd("flowPID_FB768_New_PileON_PeriodOFF");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_New_PileON_PeriodOFF");
	TH1D* hCentDist_ON_OFF = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_ON_OFF->SetLineColor(col_ON_OFF);
	hCentDist_ON_OFF->SetMarkerColor(col_ON_OFF);

	fQA_New->cd("flowPID_FB768_New_PileOFF_PeriodOFF");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_New_PileOFF_PeriodOFF");
	TH1D* hCentDist_OFF_OFF = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_OFF_OFF->SetLineColor(col_OFF_OFF);
	hCentDist_OFF_OFF->SetMarkerColor(col_OFF_OFF);

	fQA_New->cd("flowPID_FB768_New_PileOFF_PeriodON");
	TList* listTemp = (TList*) gDirectory->Get("Events_flowPID_FB768_New_PileOFF_PeriodON");
	TH1D* hCentDist_OFF_ON = (TH1D*) listTemp->FindObject("fCentDistUnitBin");
	hCentDist_OFF_ON->SetLineColor(col_OFF_ON);
	hCentDist_OFF_ON->SetMarkerColor(col_OFF_ON);

	// plotting
	TCanvas* cNew = new TCanvas("cNew","New",1000,1000);
	cNew->Divide(2,2);
	cNew->cd(1);
	hRef_v22_Gap10_New_ON_ON->SetLineColor(col_ON_ON);
	hRef_v22_Gap10_New_ON_ON->Draw();
	hRef_v22_Gap10_New_ON_OFF->SetLineColor(col_ON_OFF);
	hRef_v22_Gap10_New_ON_OFF->Draw("same");
	hRef_v22_Gap10_New_OFF_ON->SetLineColor(col_OFF_ON);
	hRef_v22_Gap10_New_OFF_ON->Draw("same");
	hRef_v22_Gap10_New_OFF_OFF->SetLineColor(col_OFF_OFF);
	hRef_v22_Gap10_New_OFF_OFF->Draw("same");
	hPRL116_v22_Gap10->Draw("same");
	legComp->SetHeader("New centrality framework (AliMultSelection)");
	legComp->Draw();

	cNew->cd(2);
	hRatio_New_ON_ON->SetMinimum(0.92);
	hRatio_New_ON_ON->SetMaximum(1.08);
	hRatio_New_ON_ON->Draw();
	hRatio_New_ON_OFF->Draw("same");
	hRatio_New_OFF_ON->Draw("same");
	hRatio_New_OFF_OFF->Draw("same");
	hRatio_New_PRL116->Draw("same");

	cNew->cd(3);
	hCentDist_ON_ON->Draw();
	hCentDist_ON_OFF->Draw("same");
	hCentDist_OFF_ON->Draw("same");
	hCentDist_OFF_OFF->Draw("same");
	
	cNew->SaveAs(Form("%s/CentNEW.pdf",sOutputPath.Data()));



	//hRef_v22_Gap10_Old_OFF_OFF->Draw("same");
}
// ======================================================================
TH1D* ConvertHist(TH1D* hist, TH1D* stat, TH1D* sys)
{
	TH1D* hOut = (TH1D*) hist->Clone();
	const Short_t iNumBins = hist->GetNbinsX();
	
	Double_t dCont = 0, dStat = 0, dSys = 0;
	for(Short_t i(1); i < iNumBins+1; i++)
	{
		dCont = hist->GetBinContent(i);
		dStat = stat->GetBinContent(i);
		dSys = sys->GetBinContent(i);

		hOut->SetBinContent(i,dCont);
		hOut->SetBinError(i,TMath::Sqrt(dStat*dStat + dSys*dSys));
	}


	return hOut;

}