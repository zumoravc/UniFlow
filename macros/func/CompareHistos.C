#include <iostream>

#include "TObject.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"


TCanvas* CompareHistos(
	TList* lList,
	TString* sLabel,
	Double_t dMinY,
	Double_t dMaxY,
	Bool_t bLog
)
{
	Color_t colors[] = {kBlack, kRed, kBlue, kMagenta, kYellow+2,kGreen,kViolet-4};

	const Int_t iNEntries = lList->GetEntries();
	
	//lList->ls();

	TLegend* tLeg = new TLegend(0.5,0.6,0.7,0.8);
	tLeg->SetBorderSize(0);

	TCanvas* cCan = new TCanvas();
	cCan->cd();
	
	if(bLog)
		gPad->SetLogy();
		

	TH1* hTemp;

	for(Int_t i = 0; i < iNEntries; i++)
	{	
		hTemp = static_cast<TH1*>(lList->At(i));
		hTemp->SetMinimum(dMinY);
		hTemp->SetMaximum(dMaxY);
		hTemp->SetLineColor(colors[i]);
		hTemp->SetMarkerColor(colors[i]);
		hTemp->SetMarkerStyle(20);
		hTemp->Draw("same");
		tLeg->AddEntry(hTemp,sLabel[i].Data(),"epl");
	}

	tLeg->Draw("same");
	return cCan;
}