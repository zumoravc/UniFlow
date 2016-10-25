//TCanvas* CompareRatio(TH1D* fNom, TH1D* fDenom);

TCanvas* CompareRatios(
		const Int_t iNEntries,
		TList* lList,
		TString* sLabel,
		Double_t dMinY,
		Double_t dMaxY
	)
{
	Color_t colors[] = {kBlack, kRed, kBlue,kYellow+2,kViolet-4,kGreen-1};
	
	//TCanvas* cCan = new TCanvas("temp","temp",800,400);
	TCanvas* cCan = new TCanvas("cCan","Canvas",800,400);
	cCan->Divide(2,1);
	
	TLegend* tLeg = new TLegend(0.6,0.13,0.9,0.4);
	tLeg->SetBorderSize(0);

	TLatex* latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);
	latex->SetTextSize(0.05);
	
	cCan->cd(1);
	TH1D* hNom[iNEntries] = {0x0};
	TH1D* hRatio[iNEntries] = {0x0};
	for(Int_t i(0); i < iNEntries; i++)
	{
		//hNom[i] = (TH1D*) (lList->At(i))->Clone(Form("hNom%d",i));
		hNom[i] = (TH1D*) (lList->At(i))->Clone(Form("hNom%d",i));
		hNom[i]->SetStats(0);
		hNom[i]->SetLineColor(colors[i]);
		hNom[i]->SetMarkerColor(colors[i]);
		hNom[i]->SetMarkerStyle(20);
		tLeg->AddEntry(hNom[i],sLabel[i].Data(),"epl");
		hNom[i]->Draw("same");
		
		hRatio[i] = (TH1D*) (lList->At(i))->Clone(Form("hRatio%d",i));
		hRatio[i]->SetLineColor(colors[i]);
		hRatio[i]->SetMarkerColor(colors[i]);
		hRatio[i]->SetMarkerStyle(20);
	}
 
 	tLeg->Draw("same");

 	cCan->cd(2);


	for(Int_t i(1); i < iNEntries; i++)
	{
		hRatio[i]->SetStats(0);
		hRatio[i]->SetMinimum(0.5);
		hRatio[i]->SetMaximum(1.5);
		hRatio[i]->Divide(hNom[0]);
		hRatio[i]->Draw("same");
	}

	latex->DrawLatex(0.2,0.8,"Me/QM");

 	TLine* lUnity = new TLine((hNom[0]->GetXaxis())->GetXmin(),1.,(hNom[0]->GetXaxis())->GetXmax(),1.);
	lUnity->SetLineColor(kBlack);
	lUnity->SetLineStyle(2);
	lUnity->Draw("same");

	return cCan;
}


/*


TCanvas* CompareRatio(TH1D* fNom, TH1D* fDenom)
{
	TCanvas* cCan = new TCanvas();
	cCan->Divide(2,1);

	TString sLabel[] = {"Me","Other"};

	TList* list = new TList();
	list->Add(fNom);
	list->Add(fDenom);

	gPad = CompareHistos(list,sLabel,0,0.5,0);

	Double_t dVal,dErr;
	Double_t dNomVal, dDenomVal, dNomErr, dDenomErr;
	TH1D* hRatio = (TH1D*) fNom->Clone("hRatio");
	for(Int_t i(1); i < hRatio->GetNbinsX()+1; i++)
	{
		dNomVal = fNom->GetBinContent(i);
		dNomErr = fNom->GetBinError(i);
		dDenomVal = fDenom->GetBinContent(i);
		dDenomErr = fDenom->GetBinError(i);

		dVal = dNomVal / dDenomVal;
		dErr = TMath::Power(dNomErr/dDenomVal,2) + TMath::Power(dNomVal*dDenomErr/(dDenomVal*dDenomVal),2);
		hRatio->SetBinContent(i,dVal);
		hRatio->SetBinError(i,TMath::Sqrt(dErr));
	}
	
	hRatio->SetMinimum(0.5);
	hRatio->SetMaximum(1.5);
	//hRatio->SetTitle("Me/Other");
	hRatio->SetLineColor(kRed);
	hRatio->SetMarkerColor(kRed);


	TLatex* latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(42);
	latex->SetTextSize(0.05);

	TLine* lUnity = new TLine((fNom->GetXaxis())->GetXmin(),1.,(fNom->GetXaxis())->GetXmax(),1.);
	lUnity->SetLineColor(kBlack);
	lUnity->SetLineStyle(2);

	TCanvas* cCan = new TCanvas();
	cCan->Divide(2,1);
	cCan->cd(1);
	fNom->Draw("ep");
	fDenom->Draw("same");

	cCan->cd(2);
	hRatio->Draw();
	lUnity->Draw("same");

	latex->DrawLatex(0.15,0.85,"Me/Other");


	return cCan;
}
*/