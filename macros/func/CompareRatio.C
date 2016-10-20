TCanvas* CompareRatio(TH1D* fNom, TH1D* fDenom)
{
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