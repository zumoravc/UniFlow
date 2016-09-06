/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlowPID_H
#define AliAnalysisTaskFlowPID_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowPID : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskFlowPID();
                                AliAnalysisTaskFlowPID(const char *name);
        virtual                 ~AliAnalysisTaskFlowPID();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void										SetAODAnalysis(Bool_t aodAnalysis) { fAODAnalysis = aodAnalysis; }
        void										SetPbPbAnalysis(Bool_t pbpb) { fPbPb = pbpb; }
        void										SetPVtxZ(Double_t dPVtxCutZ) { fPVtxCutZ = dPVtxCutZ; }
        void										SetCentEdgeLow(Double_t dCentLow) { fCentEdgeLow = dCentLow; }
        void										SetCentEdgeUp(Double_t dCentHigh) { fCentEdgeUp = dCentHigh; }
    private:
        Bool_t IsEventSelected(const AliAODEvent* event);
        void EventQA(const AliAODEvent* event);
        
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! main output list
        TList*									fOutputListQA;	//! QA output list

        Bool_t									fAODAnalysis;		//! is AOD analysis?
        Bool_t									fPbPb;					//! is PbPb analysis?
        Double_t 								fCentEdgeLow;		//! centrality low edge
        Double_t 								fCentEdgeUp;		//! centrality upper edge
        Double_t 								fPVtxCutZ; 			//! PV z cut
        
        //std histos
        TH1F*                   fHistPt;        //! dummy histogram

        // QA histos
        TH1D* 									fEventCounter;  //! event rejection tracker
        TH1D* 									fQAPVz;					//! PV z distance distribution
        TH1D*										fQANumTracks;		//! number of AOD tracks distribution


        AliAnalysisTaskFlowPID(const AliAnalysisTaskFlowPID&); // not implemented
        AliAnalysisTaskFlowPID& operator=(const AliAnalysisTaskFlowPID&); // not implemented

        ClassDef(AliAnalysisTaskFlowPID, 1);
};

#endif
