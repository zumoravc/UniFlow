
#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"





ValidateFlow(
		TString sOutPath = "~/NBI/Flow/temp/comp/",
	 	TString sInputFile = "~/NBI/Flow/temp/Flow-10samples.root",
		Bool_t bYou = kFALSE,
		Bool_t bYouRatio = kFALSE,
		Bool_t bPID = kFALSE
	)
{

	// GENERAL STUFF

	Color_t colMine = kBlue;
	Color_t colPRL = kGreen+2;
	Color_t colYou = kRed;

	TLine lineUnityRef = TLine(0,1.,80,1.);
	lineUnityRef.SetLineStyle(kDashed);
	lineUnityRef.SetLineColor(kBlack);

	TLine lineUnityDif = TLine(0.2,1.,5.,1.);
	lineUnityDif.SetLineStyle(kDashed);
	lineUnityDif.SetLineColor(kBlack);


	TLatex latex;
	latex.SetNDC();

	// loading MY input
	printf("=================== Loading Mine =====================\n");
	TFile* fInput = new TFile(sInputFile.Data(),"READ");
	fInput->cd();
	fInput->ls();
	
	TH1D* hRef_v22_Gap10 = (TH1D*) gDirectory->Get("hRef_n22_Gap10");
	hRef_v22_Gap10->SetLineColor(colMine);
	hRef_v22_Gap10->SetMarkerColor(colMine);
	TH1D* hRef_v32_Gap10 = (TH1D*) gDirectory->Get("hRef_n32_Gap10");
	hRef_v32_Gap10->SetLineColor(colMine);
	hRef_v32_Gap10->SetMarkerColor(colMine);
	TH1D* hRef_v24 = (TH1D*) gDirectory->Get("hRef_n24_Gap-10");
	hRef_v24->SetLineColor(colMine);
	hRef_v24->SetMarkerColor(colMine);

	TList* listTrack = (TList*) gDirectory->Get("Tracks_n2_Gap10");
	TH1D* hTracks_v22_Gap10_Cent05 = (TH1D*) listTrack->FindObject("fTracks_n22_gap10_cent0_number0_0_px_desampled");
	TH1D* hTracks_v22_Gap10_Cent3040 = (TH1D*) listTrack->FindObject("fTracks_n22_gap10_cent4_number0_0_px_desampled");
	printf("======================================================\n");

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
	TH1D* hPRL116_v24_orig = (TH1D*) gDirectory->Get("Hist1D_y2");
	hPRL116_v24_orig->SetLineColor(colPRL);
	hPRL116_v24_orig->SetMarkerColor(colPRL);
	TH1D* hPRL116_v24_stat = (TH1D*) gDirectory->Get("Hist1D_y2_e1");
	TH1D* hPRL116_v24_sys = (TH1D*) gDirectory->Get("Hist1D_y2_e2");


	
	// loading PRL 107
	printf("=================== Loading PRL107 ===================\n");
	TFile* fInputPRL107 = new TFile("~/NBI/Flow/hepData/PRL107/HEPData-ins900651-1-root.root","READ");
	fInputPRL107->cd("Table 22"); // v22 Gap 10
	gDirectory->ls();
	TGraphAsymmErrors* gPRL107_v22_Gap10_centO5 = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1");
	gPRL107_v22_Gap10_centO5->SetLineColor(colPRL);
	gPRL107_v22_Gap10_centO5->SetMarkerColor(colPRL);
	
	fInputPRL107->cd("Table 14"); // v22 Gap 10
	//gDirectory->ls();
	TGraphAsymmErrors* gPRL107_v22_Gap10_cent3040 = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y1");
	gPRL107_v22_Gap10_cent3040->SetLineColor(colPRL);
	gPRL107_v22_Gap10_cent3040->SetMarkerColor(colPRL);
	

	printf("======================================================\n");

	// loading You 
	printf("=================== Loading You ======================\n");
	TFile* fInputYouRef= new TFile("~/NBI/Flow/results/you/FB768/Analysis_vn_2468.root","READ");
	fInputYouRef->cd();
	fInputYouRef->ls();
	TH1D* hYouRef_v22_Gap10 = (TH1D*)	gDirectory->Get("fv22Gap10");
	hYouRef_v22_Gap10->SetLineColor(colYou);
	hYouRef_v22_Gap10->SetMarkerColor(colYou);
	TH1D* hYouRef_v32_Gap10 = (TH1D*)	gDirectory->Get("fv32Gap10");
	hYouRef_v32_Gap10->SetLineColor(colYou);
	hYouRef_v32_Gap10->SetMarkerColor(colYou);
	TH1D* hYouRef_v24 = (TH1D*) gDirectory->Get("fv24");
	hYouRef_v24->SetLineColor(colYou);
	hYouRef_v24->SetMarkerColor(colYou);

	TFile* fInputYouPt2 = new TFile("/Users/vpacik/NBI/Flow/results/you/FB768/Analysis_v2pt.root","READ");
	fInputYouPt2->cd();
	TH1D* hYouDif_v22_Gap10_cent05 = (TH1D*) gDirectory->Get("hisv22Gap10_cent0");
	hYouDif_v22_Gap10_cent05->SetLineColor(colYou);
	hYouDif_v22_Gap10_cent05->SetMarkerColor(colYou);
	TH1D* hYouDif_v22_Gap10_cent3040 = (TH1D*) gDirectory->Get("hisv22Gap10_cent4");
	hYouDif_v22_Gap10_cent3040->SetLineColor(colYou);
	hYouDif_v22_Gap10_cent3040->SetMarkerColor(colYou);
	printf("======================================================\n");

	// ratio

	TH1D* hRatio_YouRef_v22_Gap10 = (TH1D*) hRef_v22_Gap10->Clone("hRatio_YouRef_v22_Gap10");
	hRatio_YouRef_v22_Gap10->Divide(hYouRef_v22_Gap10);
	hRatio_YouRef_v22_Gap10->SetLineColor(colYou);
	hRatio_YouRef_v22_Gap10->SetMarkerColor(colYou);
	
	TH1D* hRatio_YouRef_v32_Gap10 = (TH1D*) hRef_v32_Gap10->Clone("hRatio_YouRef_v32_Gap10");
	hRatio_YouRef_v32_Gap10->Divide(hYouRef_v32_Gap10);
	hRatio_YouRef_v32_Gap10->SetLineColor(colYou);
	hRatio_YouRef_v32_Gap10->SetMarkerColor(colYou);

	TH1D* hPRL116_v22_Gap10 = ConvertHist(hPRL116_v22_Gap10_orig,hPRL116_v22_Gap10_stat,hPRL116_v22_Gap10_sys);
	TH1D* hRatio_PRL116_v22_Gap10 = (TH1D*) hRef_v22_Gap10->Clone("hRatio_PRL116_v22_Gap10");
	hRatio_PRL116_v22_Gap10->Divide(hPRL116_v22_Gap10);
	hRatio_PRL116_v22_Gap10->SetLineColor(colPRL);
	hRatio_PRL116_v22_Gap10->SetMarkerColor(colPRL);

	TH1D* hPRL116_v24 = ConvertHist(hPRL116_v24_orig,hPRL116_v24_stat,hPRL116_v24_sys);
	TH1D* hRatio_PRL116_v24 = (TH1D*) hRef_v24->Clone("hRatio_PRL116_v24");
	hRatio_PRL116_v24->Divide(hPRL116_v24);
	hRatio_PRL116_v24->SetLineColor(colPRL);
	hRatio_PRL116_v24->SetMarkerColor(colPRL);

	TH1D* hRatio_YouRef_v24 = (TH1D*) hRef_v24->Clone("hRatio_YouRef_v24");
	hRatio_YouRef_v24->Divide(hYouRef_v24);
	hRatio_YouRef_v24->SetLineColor(colYou);
	hRatio_YouRef_v24->SetMarkerColor(colYou);
	
	TH1D* hRatio_YouDif_v22_Gap10_cent05 = (TH1D*) hTracks_v22_Gap10_Cent05->Clone("hRatio_YouDif_v22_Gap10_cent05");
	hRatio_YouDif_v22_Gap10_cent05->Divide(hYouDif_v22_Gap10_cent05);
	hRatio_YouDif_v22_Gap10_cent05->SetLineColor(colYou);
	hRatio_YouDif_v22_Gap10_cent05->SetMarkerColor(colYou);
	
	TH1D* hRatio_YouDif_v22_Gap10_cent3040 = (TH1D*) hTracks_v22_Gap10_Cent3040->Clone("hRatio_YouDif_v22_Gap10_cent3040");
	hRatio_YouDif_v22_Gap10_cent3040->Divide(hYouDif_v22_Gap10_cent3040);
	hRatio_YouDif_v22_Gap10_cent3040->SetLineColor(colYou);
	hRatio_YouDif_v22_Gap10_cent3040->SetMarkerColor(colYou);

	TLegend* legAbs = new TLegend(0.32,0.12,0.65,0.35);
	legAbs->SetBorderSize(0);
	legAbs->AddEntry(hRef_v22_Gap10, "Mine","pel");
	legAbs->AddEntry(hPRL116_v22_Gap10, "PRL116","pel");
	legAbs->AddEntry(hYouRef_v22_Gap10, "You","pel");

	TLegend* legDiff = new TLegend(0.32,0.12,0.75,0.35);
	legDiff->SetBorderSize(0);
	legDiff->AddEntry(hRef_v22_Gap10, "Mine","pel");
	legDiff->AddEntry(hPRL116_v22_Gap10, "PRL107","pel");
	legDiff->AddEntry(hYouDif_v22_Gap10_cent05, "You","pel");

	TLegend* legRatio = new TLegend(0.12,0.7,0.5,0.89);
	legRatio->SetHeader("Ratio Mine / X");
	legRatio->SetBorderSize(0);
	//legRatio->AddEntry(hRef_v22_Gap10, "Mine","pel");
	legRatio->AddEntry(hPRL116_v22_Gap10, "PRL116","pel");
	legRatio->AddEntry(hYouRef_v22_Gap10, "You","pel");

	TLegend* legRatioPRL107 = new TLegend(0.12,0.7,0.5,0.89);
	legRatioPRL107->SetHeader("Ratio Mine / X");
	legRatioPRL107->SetBorderSize(0);
	//legRatio->AddEntry(hRef_v22_Gap10, "Mine","pel");
	legRatioPRL107->AddEntry(hPRL116_v22_Gap10, "PRL107","pel");
	legRatioPRL107->AddEntry(hYouRef_v22_Gap10, "You","pel");

	// Plotting
	TCanvas* cPRL = new TCanvas();
	cPRL->Divide(2,1);
	
	// v22
	cPRL->cd(1);
	hRef_v22_Gap10->Draw();
	hYouRef_v22_Gap10->Draw("same");
	//hPRL116_v22_Gap10->SetLineColor(kRed);
	//hPRL116_v22_Gap10->Draw("same");
	gPRL116_v22_Gap10->Draw("P");
	legAbs->Draw();
	latex.DrawLatex(0.3,0.4,"v_{2}{2,Gap 1.0}");

	cPRL->cd(2);
	hRatio_PRL116_v22_Gap10->SetMinimum(0.9);
	hRatio_PRL116_v22_Gap10->SetMaximum(1.1);
	hRatio_PRL116_v22_Gap10->Draw();
	hRatio_YouRef_v22_Gap10->Draw("same");
	legRatio->Draw();
	lineUnityRef.Draw();

	cPRL->SaveAs(Form("%s/Ref/v22_Gap10.pdf",sOutPath.Data()));
	
	// v32
	 // to be revisited - some error
	cPRL->cd(1);
	hRef_v32_Gap10->Draw();
	hYouRef_v32_Gap10->Draw("same");
	//hPRL116_v22_Gap10->SetLineColor(kRed);
	//hPRL116_v22_Gap10->Draw("same");
	//gPRL116_v22_Gap10->Draw("P");
	legAbs->Draw();
	latex.DrawLatex(0.5,0.4,"v_{3}{2,Gap 1.0}");

	cPRL->cd(2);
	//hRatio_PRL116_v22_Gap10->SetMinimum(0.9);
	//hRatio_PRL116_v22_Gap10->SetMaximum(1.1);
	//hRatio_PRL116_v22_Gap10->Draw();
	hRatio_YouRef_v32_Gap10->SetMinimum(0.9);
	hRatio_YouRef_v32_Gap10->SetMaximum(1.1);
	hRatio_YouRef_v32_Gap10->Draw();
	legRatio->Draw();
	lineUnityRef.Draw();

	cPRL->SaveAs(Form("%s/Ref/v32_Gap10.pdf",sOutPath.Data()));


	// v24
	cPRL->cd(1);
	hRef_v24->SetMinimum(0.);
	hRef_v24->SetMaximum(0.13);
	hRef_v24->Draw();
	hYouRef_v24->Draw("same");
	gPRL116_v24->Draw("P");
	legAbs->Draw();
	latex.DrawLatex(0.3,0.4,"v_{2}{4}");

	cPRL->cd(2);
	hRatio_PRL116_v24->SetMinimum(0.9);
	hRatio_PRL116_v24->SetMaximum(1.1);
	hRatio_PRL116_v24->Draw();
	hRatio_YouRef_v24->Draw("same");
	lineUnityRef.Draw();
	legRatio->Draw();

	cPRL->SaveAs(Form("%s/Ref/v24.pdf",sOutPath.Data()));

	// pT diff v22 cent 0-5
	cPRL->cd(1);
	hTracks_v22_Gap10_Cent05->Draw();
	hYouDif_v22_Gap10_cent05->Draw("same");
	gPRL107_v22_Gap10_centO5->Draw("P");
	//legDiff->SetHeader("v_{2}{2, Gap10} 0-5%");
	latex.DrawLatex(0.3,0.4,"v_{2}{2,Gap 1.0} 0-5%");
	legDiff->Draw();

	cPRL->cd(2);
	hRatio_YouDif_v22_Gap10_cent05->SetMinimum(0.9);
	hRatio_YouDif_v22_Gap10_cent05->SetMaximum(1.1);
	hRatio_YouDif_v22_Gap10_cent05->Draw();
	lineUnityDif.Draw();
	legRatioPRL107->Draw();
	cPRL->SaveAs(Form("%s/Tracks/v22_Gap10_cent05.pdf",sOutPath.Data()));
	
	// pT diff v22 cent 30-40
	cPRL->cd(1);
	hTracks_v22_Gap10_Cent3040->Draw();
	hYouDif_v22_Gap10_cent3040->Draw("same");
	gPRL107_v22_Gap10_cent3040->Draw("P");
	//legDiff->SetHeader("v_{2}{2, Gap10} 30-40%");
	legDiff->Draw();
	latex.DrawLatex(0.3,0.4,"v_{2}{2,Gap 1.0} 30-40%");

	cPRL->cd(2);
	hRatio_YouDif_v22_Gap10_cent3040->SetMinimum(0.9);
	hRatio_YouDif_v22_Gap10_cent3040->SetMaximum(1.1);
	hRatio_YouDif_v22_Gap10_cent3040->Draw();
	lineUnityDif.Draw();
	legRatio->Draw();
	cPRL->SaveAs(Form("%s/Tracks/v22_Gap10_cent3040.pdf",sOutPath.Data()));

}

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



