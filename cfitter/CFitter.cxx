/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* CFitter
 *
 * Class implemented (primarily for PID flow) fitting purposes of identified hardon extraction
 * 
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
 */

#include "CFitter.h"

class CFitter;    

//ClassImp(CFitter) // classimp: necessary for root

//_____________________________________________________________________________
CFitter::CFitter() 
{
	fChi2 = -1;
}
//_____________________________________________________________________________
CFitter::~CFitter()
{
}
//_____________________________________________________________________________