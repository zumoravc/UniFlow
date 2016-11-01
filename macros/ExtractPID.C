/* Extract PID 
 *
 * Macro for extraction of identified hadrons yields and flow.
 * Using FPFitter class.
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */


#include "TROOT.h"

void ExtractPID()
{
	gROOT->LoadMacro("~/NBI/Flow/cfitter/FPFitter.cxx++g"); // compiling CFitter locally

	FPFitter* fit = new FPFitter();
	printf("Chi2 %g \n",fit->GetChi2());

}



int main()
{
	ExtractPID();
	return 0;
}