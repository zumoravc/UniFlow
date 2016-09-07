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

        void										SetAODAnalysis(Bool_t aod) { fAODAnalysis = aod; }
        void										SetPbPbAnalysis(Bool_t pbpb) { fPbPb = pbpb; }
        void										SetCentEdgeLow(Double_t cent) { fCentEdgeLow = cent; }
        void										SetCentEdgeUp(Double_t cent) { fCentEdgeUp = cent; }
        void										SetPVtxZ(Double_t z) { fPVtxCutZ = z; }
        void										SetTrackEtaMax(Double_t eta) { fTrackEtaMax = eta; }
        void										SetTrackPtMax(Double_t pt) { fTrackPtMax = pt; }
        void										SetTrackPtMin(Double_t pt) { fTrackPtMin = pt; }
        void										SetNumTPCclsMin(UShort_t tpcCls) { fNumTPCclsMin = tpcCls; }
        void										SetTrackFilterBit(UInt_t filter) { fTrackFilterBit = filter; }
    private:
        Bool_t IsEventSelected(const AliAODEvent* event);
				Bool_t IsTrackSelected(const AliAODTrack* track);
       	void EventQA(const AliAODEvent* event);
        
        TList*                  fOutputList;    //! main output list
        TList*									fOutputListQA;	//! QA output list
        AliAODEvent*            fAOD;           //! input event
        AliAODTrack*						fTrack;					//! AOD track

        Bool_t									fAODAnalysis;		//! is AOD analysis?
        Bool_t									fPbPb;					//! is PbPb analysis?
        Double_t 								fCentEdgeLow;		//! centrality low edge
        Double_t 								fCentEdgeUp;		//! centrality upper edge
        Double_t 								fPVtxCutZ; 			//! PV z cut
        Double_t 								fTrackEtaMax; 	//! Maximum pseudorapidity range
        Double_t								fTrackPtMax;		//! Maximal track pT
        Double_t								fTrackPtMin;		//! Manimal track pT
        UShort_t								fNumTPCclsMin;	//! Minimal number of TPC clusters used for track reconstruction
        UInt_t									fTrackFilterBit;//! Required track filter bit 
        
        //std histos
        TH1D*										fEventMult;			 //! selected events multiplicity distribution
        TH1D*                   fTracksPt;       //! selected tracks pT distribution
        TH1D*                   fTracksEta;      //! selected tracks eta distribution

        // QA histos
        TH1D* 									fEventCounter;  //! event rejection tracker
        TH1D* 									fQAPVz;					//! PV z distance distribution
        TH1D*										fQANumTracks;		//! number of AOD tracks distribution


        AliAnalysisTaskFlowPID(const AliAnalysisTaskFlowPID&); // not implemented
        AliAnalysisTaskFlowPID& operator=(const AliAnalysisTaskFlowPID&); // not implemented

        ClassDef(AliAnalysisTaskFlowPID, 1);
};

#endif
