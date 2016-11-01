#include "TROOT.h"

void ExtractPID()
{
	gROOT->LoadMacro("~/NBI/Flow/cfitter/CFitter.cxx++g"); // compiling CFitter locally

	CFitter* fit = new CFitter();
	printf("Chi2 %g \n",fit->GetChi2());

}



int main()
{
	ExtractPID();
	return 0;
}