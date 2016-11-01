#include "TROOT.h"

#include "FPFitter.h"
#include "FPFitter.cxx"

void testFitter()
{
	gROOT->LoadMacro("FPFitter.cxx++g"); // compiling CFitter locally


	FPFitter* fit = new FPFitter();
	printf("Chi2 %g \n",fit->GetChi2());

}

int main()
{
	testFitter();
	return 0;
}
