/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */

#ifndef FPFitter_H
#define FPFitter_H

class FPFitter
{
	public:
		FPFitter();
    ~FPFitter();

    Double_t GetChi2() { return fChi2; }

	private:
		Double_t fChi2; // chi2 of final fit

	ClassDef(FPFitter, 0);
};

#endif

