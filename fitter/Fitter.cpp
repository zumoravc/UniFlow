/* Fitter
 *
 * Class implemented for PID flow fitting purposes.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"

#include "RooPolynomial.h"
#include "RooGaussian.h"


class Fitter
{
public:
		Fitter();
    ~Fitter();

    void	Run(); // start the fitting procedure
    void	TestRooFit(); // RooFit testing tutorial
    void 	RooFitInvMass(); // invariant mass fitting by RooFit
    void  FitInvMass(); // classical invariant mass fitting via ROOT fit
    void	AttachInvMass(TH1* hInvMass) { fhInvMass = (TH1*) hInvMass->Clone(); } // attach input inv. mass histo for fitting
    void	AttachFlowMass(TH1* hFlowMass) { fhFlowMass = (TH1*) hFlowMass->Clone(); } // attach input flow mass histo for fitting
    void	SetRangeMassMin(Double_t mass) { fdRangeMassMin = mass; } // set minimum range of inv. mass (X axis)
    void	SetRangeMassMax(Double_t mass) { fdRangeMassMax = mass; } // set maximum range of inv. mass (X axis)

    Double_t GetChi2() { return fChi2; }

	private:
		TH1*	fhInvMass; //! invariant mass histogram 
		TH1*	fhFlowMass; //! flow mass histogram 

		TF1* 	ffInvMassFirst; // InvMass: first shot

		Double_t fdRangeMassMin; // Minimum histogram range X (inv mass)
		Double_t fdRangeMassMax; // Maximum histogram range X (inv mass)

		Double_t fChi2; // chi2 of final fit

};

//_____________________________________________________________________________
Fitter::Fitter() 
{
	fhInvMass = 0x0;
	fhFlowMass = 0x0;
	
	fdRangeMassMax = 1;
	fdRangeMassMin = 0;
	fChi2 = -1;
}
//_____________________________________________________________________________
Fitter::~Fitter() 
{
	delete fhInvMass;
	delete fhFlowMass;
}
//_____________________________________________________________________________
void Fitter::Run()
{
	return;
}
//_____________________________________________________________________________
/*
*	RooFit tutorial as available at https://root.cern.ch/download/doc/RooFit_Users_Manual_2.91-33.pdf
*/
void Fitter::TestRooFit()
{
	using namespace RooFit;
	RooRealVar x("x","x",-10,10);
	RooRealVar mean("mean","Mean of Gauss",0,-10,10);
	RooRealVar sigma("sigma","Width of Gauss",3,-10,10);
	RooGaussian gauss("gauss","gauss(x,mean,sigma)",x,mean,sigma);

	RooPlot* frame = x.frame();
	gauss.plotOn(frame);
	sigma = 2;
	gauss.plotOn(frame,LineColor(kRed));
	frame->Draw();

	RooDataHist data("data","dataset with x",x,fhInvMass);

	RooPlot* histFrame = x.frame();
	data.plotOn(histFrame,DataError(RooAbsData::SumW2));
	histFrame->Draw();

	gauss.fitTo(data);



	return;
}
//_____________________________________________________________________________
void Fitter::RooFitInvMass()
{
	using namespace RooFit;

	RooRealVar mass("mass","M_{inv}",fdRangeMassMin,fdRangeMassMax);
	RooDataHist invMass("invMass","InvMass",mass,fhInvMass);

	RooRealVar mean("mean","mean",0.5,0.4,0.6);
	RooRealVar sigma("sigma","sigma",0.01,0.,0.02);
	RooGaussian sig("sig","signal",mass,mean,sigma);

	RooRealVar cBg0("cBg0","bg slope", 1, -100, 100);
	RooRealVar cBg1("cBg1","bg kvadr", 1, -100, 100);
	RooRealVar cBg2("cBg2","bg cube", 1, -100, 100);
	RooPolynomial bg("bg","background",mass,RooArgList(cBg0,cBg1,cBg2));

	RooRealVar frac("frac","signal fraction", 0.5, 0., 1.);
	

	RooAddPdf model("model","model",RooArgList(sig,bg),frac);

	printf("--------------------- Model fitting\n");

	model.fitTo(invMass);
	
	
	TCanvas* cCan = new TCanvas();
	RooPlot* frame = mass.frame();
	invMass.plotOn(frame,DataError(RooAbsData::SumW2));
	model.plotOn(frame);
	model.plotOn(frame,Components(bg),LineStyle(kDashed));	
	frame->Draw();

	RooRealVar numSig("numSig","Number of signal",500, 0, 200000);
	RooRealVar numBg("numBg","Number of background",500, 0, 200000);

	printf("----------------------- Extended Model fitting\n");
	RooAddPdf modelExt("modelExt","Extended model",RooArgList(sig,bg),RooArgList(numSig,numBg));
	modelExt.fitTo(invMass);

	TCanvas* cCan2 = new TCanvas();
	RooPlot* frame2 = mass.frame();
	invMass.plotOn(frame2,DataError(RooAbsData::SumW2));
	modelExt.plotOn(frame2,LineColor(kRed));
	modelExt.plotOn(frame2,Components(bg),LineColor(kRed),LineStyle(kDashed));	
	frame2->Draw();
}
//_____________________________________________________________________________
void Fitter::FitInvMass()
{
	// from V0sExtractFlow.C

	TCanvas* cCan = new TCanvas();
	

	TH1* hInvMass = (TH1*) fhInvMass;

	const Short_t dNumSigma = 5;
	Double_t dMassPeak = 0., dMassPeakLimit = 0.;
	Double_t dFitRange[2] = {0};
	TString sPartName;

	sPartName = TString("K0s");
	dMassPeak = 0.49;
	dMassPeakLimit = 0.03;
	dFitRange[0] = 0.4;
	dFitRange[1] = 0.6;

	printf("=== Fitting %s ===\n", sPartName.Data());

	TH1D* hInvMass_temp = (TH1D*) hInvMass->Clone("hInvMass_temp"); // cloning inv mass dist for peak window fitting
	TH1D* hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side"); // cloning inv mass dist for sidebands fitting
	TH1D* hInvMass_subs = (TH1D*) hInvMass->Clone("hInvMass_subs"); // cloning inv mass dist for BG subtracttion
	TH1D* hInvMassRebin = (TH1D*) hInvMass->Clone("hInvMassRebin"); // cloning inv mass dist for rebining

	// ==== Fitting inv mass dist =====
	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass first shot\n");
	printf("-----------------------------------------\n");
	
	TF1* fFitInvMass = new TF1("fFitInvMass","gaus(0)+pol3(3)",dFitRange[0],dFitRange[1]); 
	fFitInvMass->SetParNames("Amp","Mean","Sigma");
	fFitInvMass->SetNpx(100000);
	fFitInvMass->SetParameter(0,hInvMass->GetMaximum());
	fFitInvMass->SetParameter(1,dMassPeak);
	fFitInvMass->SetParameter(2,0.005);
	fFitInvMass->SetParLimits(2,0.,0.01);
	//fFitInvMass->SetLineColor(kGreen);
	fFitInvMass->SetLineWidth(2);
	//hInvMass->SetMinimum(0);
	//hInvMass->SetMaximum(4000);	
	hInvMass->Fit("fFitInvMass","R0");

	// * TODO - fit verification procedure * //

	// extracting mean & sigma
	Double_t dMean = fFitInvMass->GetParameter(1);
	Double_t dSigma = fFitInvMass->GetParameter(2);
	
	Double_t dMassWindow[2] = {dMean - dNumSigma*dSigma, dMean + dNumSigma*dSigma}; // setting inv. mass peak
	printf("dMassWindow: %g - %g\n", dMassWindow[0],dMassWindow[1]);


	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hInvMass->GetBinCenter(i) > dMassWindow[0] && hInvMass->GetBinCenter(i) < dMassWindow[1])
		{
			hInvMass_side->SetBinError(i,9999999999999); 
		}
	}

	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass in sidebins\n");
	printf("-----------------------------------------\n");
	
	TF1* fFitInvMassBg = new TF1("fFitInvMassBg","pol3",dFitRange[0],dFitRange[1]);
	fFitInvMassBg->SetLineWidth(2);
	fFitInvMassBg->SetLineStyle(2);
	fFitInvMassBg->SetLineColor(kGreen+2);
	hInvMass_side->SetMinimum(0);
	hInvMass_side->SetMaximum(100000);
  hInvMass_side->Fit("fFitInvMassBg","R");
  
  /*	
	for(Int_t i(1); i < hInvMass->GetNbinsX()+1; i++)
	{
		// subtracting the BG
		hInvMass_subs->SetBinContent(i,hInvMass->GetBinContent(i) - fFitInvMassBg->Eval(hInvMass_subs->GetBinCenter(i)));
	}

	TH1D* hInvMass_RatioSigTot = (TH1D*) hInvMass_subs->Clone("hInvMass_RatioSigTot");
	hInvMass_RatioSigTot->Divide(hInvMass);
	
	printf("-----------------------------------------\n");
	printf("Fitting Inv Mass BG subtracted\n");
	printf("-----------------------------------------\n");
	TF1* fFitInvMassSubs = new TF1("fFitInvMassSubs","gaus(0)+pol3(3)",dFitRange[0],dFitRange[1]); 
	fFitInvMassSubs->SetParameter(1,dMassPeak);
	fFitInvMassSubs->SetParameter(2,0.005);
	fFitInvMassSubs->SetParLimits(2,0.,0.02);

	fFitInvMassSubs->SetLineWidth(2);
	fFitInvMassSubs->SetLineStyle(2);
	fFitInvMassSubs->SetLineColor(kRed+2);
	fFitInvMassSubs->SetNpx(1000000);
	hInvMass_subs->Fit("fFitInvMassSubs","R");
	hInvMass_subs->SetMaximum(2000);
	*/

	hInvMass->Draw();
	fFitInvMass->Draw("same");
	fFitInvMassBg->Draw("same");
}
