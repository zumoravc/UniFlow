#include "TROOT.h"

#include "CFitter.h"

void testFitter()
{
	gROOT->LoadMacro("CFitter.cxx+");


	CFitter* fit = new CFitter();
	printf("Chi2 %f",fit->GetChi2());

}