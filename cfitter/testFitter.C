#include "TROOT.h"

#include "CFitter.h"
#include "CFitter.cxx"

void testFitter()
{
	gROOT->LoadMacro("CFitter.cxx++g"); // compiling CFitter locally


	CFitter* fit = new CFitter();
	printf("Chi2 %g \n",fit->GetChi2());

}

int main()
{
	testFitter();
	return 0;
}
