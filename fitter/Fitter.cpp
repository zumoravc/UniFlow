/* Fitter
 *
 * Class implemented for PID flow fitting purposes.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "TH1.h"

class Fitter
{
public:
		Fitter();
    ~Fitter();

    void	Run(); // main macro which 
    void	AttachInvMass(TH1* hInvMass) { fhInvMass = (TH1*) hInvMass->Clone(); } // attach input inv. mass histo for fitting
    void	AttachFlowMass(TH1* hFlowMass) { fhFlowMass = (TH1*) hFlowMass->Clone(); } // attach input flow mass histo for fitting

    Double_t GetChi2() { return fChi2; }

	private:
		TH1*	fhInvMass; //! invariant mass histogram 
		TH1*	fhFlowMass; //! flow mass histogram 

		Double_t fChi2; // chi2 of final fit
};

//_____________________________________________________________________________
Fitter::Fitter() 
{
	fhInvMass = 0x0;
	fhFlowMass = 0x0;
	
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


