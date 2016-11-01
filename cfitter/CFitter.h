/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */

#ifndef CFitter_H
#define CFitter_H

class CFitter
{
	public:
		CFitter();
    ~CFitter();

    Double_t GetChi2() { return fChi2; }

	private:
		Double_t fChi2; // chi2 of final fit

	ClassDef(CFitter, 0);
};

#endif

