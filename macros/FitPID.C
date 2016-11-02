/* Fit PID 
 *
 * Macro for extraction of identified hadrons yields and flow.
 * Using Fitter class.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "Fitter.cpp"
#include "TFile.h"
#include "TH1D.h"

void FitPID(
	
		const TString sInput = "~/NBI/Flow/results/V0s/15-sampling/sampling_test/MassDist_V0s_Gap10.root",
		const TString sOutput = "~/NBI/Flow/results/V0s/15-sampling/sampling_test/plots_Gap10",
		const TString sTag = "flowPID_JHEP",
		const TString sEtaGap = "Gap10",
		const TString sOutputFormat = "png"
		
	)
{
	/*
	const Int_t iNumPtBins = 22; // pT bins
	const Int_t iNumCentBins = 9; // centrality bins
	// bins edges
	Double_t fPtBinEdges[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4}; 
	Double_t fCentBinEdges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};

	// =======================================
	
	//gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareRatio.C");
	//gROOT->LoadMacro("~/NBI/Flow/macros/func/CompareHistos.C");

	TFile* fInput = new TFile(sInput.Data(),"READ");
	TFile* fOutput = new TFile(Form("%s/../PtFlow_V0s_%s.root",sOutput.Data(),sEtaGap.Data()),"RECREATE");
	
	if(!fInput->IsOpen())
		return;

	fInput->cd();

	// loading input
	printf("==== Input Loading");
	TList* lInputInvMassK0s = (TList*) gDirectory->Get("lInvMass_K0s");
	TList* lInputInvMassLambda = (TList*) gDirectory->Get("lInvMass_Lambda");
	TList* lInputFlowMassK0s = (TList*) gDirectory->Get("lFlowMass_K0s");
	TList* lInputFlowMassLambda = (TList*) gDirectory->Get("lFlowMass_Lambda");
	

	TH1D* hInvMassK0s[iNumCentBins][iNumPtBins];
	TH1D* hInvMassLambda[iNumCentBins][iNumPtBins];
	TH1D* hFlowMassK0s[iNumCentBins][iNumPtBins];
	TH1D* hFlowMassLambda[iNumCentBins][iNumPtBins];
	for(Int_t i(0); i < iNumCentBins; i++)
	{
		for(Int_t j(0); j < iNumPtBins; j++)
		{
			hInvMassK0s[i][j] = (TH1D*) lInputInvMassK0s->FindObject(Form("hInvMass_K0s_Cent%d_pt%d",i,j))->Clone(Form("hInvMassK0s_Cent%d_pt%d",i,j));
			hInvMassLambda[i][j] = (TH1D*) lInputInvMassLambda->FindObject(Form("hInvMass_Lambda_Cent%d_pt%d",i,j))->Clone(Form("hInvMassLambda_Cent%d_pt%d",i,j));
			hFlowMassK0s[i][j] = (TH1D*) lInputFlowMassK0s->FindObject(Form("hFlowMass_K0s_Cent%d_pt%d",i,j))->Clone(Form("hFlowMassK0s_Cent%d_pt%d",i,j));
			hFlowMassLambda[i][j] = (TH1D*) lInputFlowMassLambda->FindObject(Form("hFlowMass_Lambda_Cent%d_pt%d",i,j))->Clone(Form("hFlowMassLambda_Cent%d_pt%d",i,j));
		}
	}

	printf(": DONE ====\n");

	*/

	// fitting 


	Fitter* fit = new Fitter();
	//fit->AttachInvMass(hInvMassK0s[3][3]);
	

	printf("Chi2 %g \n",fit->GetChi2());

}



int main()
{
	FitPID();
	return 0;
}