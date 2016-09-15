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
        void										SetPVtxZ(Double_t z) { fPVtxCutZ = z; }
        void										SetTrackEtaMax(Double_t eta) { fTrackEtaMax = eta; }
        void										SetTrackPtMax(Double_t pt) { fTrackPtMax = pt; }
        void										SetTrackPtMin(Double_t pt) { fTrackPtMin = pt; }
        void										SetNumTPCclsMin(UShort_t tpcCls) { fNumTPCclsMin = tpcCls; }
        void										SetTrackFilterBit(UInt_t filter) { fTrackFilterBit = filter; }
        void										SetDiffFlow(Bool_t diff) { fDiffFlow = diff; }
       
        const static Int_t 			fNumPtBins = 12;			//! number of pT bins used for pT-differential flow
        static Double_t					fPtBinEdges[fNumPtBins+1];				//! pointer for array of pT bin edges
        const static Int_t 			fNumCentBins = 9;			//! number of centrality bins used for pT-differential flow (so far independently of reference flow)
        static Double_t					fCentBinEdges[fNumCentBins+1];				//! pointer for array of pT bin edges

    private:
        Bool_t IsEventSelected(const AliAODEvent* event);
				Bool_t IsTrackSelected(const AliAODTrack* track);
				Bool_t IsV0aK0s(const AliAODv0* v0);
				Bool_t IsV0aLambda(const AliAODv0* v0);
				Bool_t IsV0aALambda(const AliAODv0* v0);
				Bool_t IsV0Selected(const AliAODv0* v0);
				void EventQA(const AliAODEvent* event);
       	void EstimateCentrality(AliVEvent* ev);
    
			  Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1); 
        Bool_t plpMV(const AliVEvent *event);

        Short_t GetPtBinIndex(const Double_t dPt);

        //cuts & selection
        Bool_t									fAODAnalysis;		//! is AOD analysis?
        Bool_t									fPbPb;					//! is PbPb analysis?
				Bool_t       						fLHC10h;             // flag to LHC10h data?
				Short_t									fCentFlag;			//! centrality flag
        Double_t 								fPVtxCutZ; 			//! PV z cut
        Double_t 								fTrackEtaMax; 	//! Maximum pseudorapidity range
        Double_t								fTrackPtMax;		//! Maximal track pT
        Double_t								fTrackPtMin;		//! Manimal track pT
        UShort_t								fNumTPCclsMin;	//! Minimal number of TPC clusters used for track reconstruction
        UInt_t									fTrackFilterBit;//! Required track filter bit 
        Bool_t 									fDiffFlow;			//! Do differential flow ? (or reference only)
        Bool_t 									fPID;						//! Do PID (so far V0s) ? 
				//cuts & selection: V0 reconstruction
				Bool_t 									fCutV0onFly;		//! V0 reconstruction method: is On-the-fly? (or offline)
				Bool_t									fCutV0refitTPC; //! Check TPC refit of V0 daughters ?
				Bool_t									fCutV0rejectKinks; //! Reject Kink V0 daughter tracks ?
				Double_t								fCutV0MinCPAK0s;	//! min cosine of pointing angle of K0s candidate to PV
				Double_t								fCutV0MinCPALambda;	//! min cosine of pointing angle of K0s candidate to PV
				Double_t								fCutV0MaxDCAtoPV;	//! max DCA of V0 daughter to PV
				Double_t								fCutV0MaxDCADaughters;	//! max DCA of V0 daughters among themselves
				Double_t								fCutV0MaxDecayRadius; //! max distance between PV and secondary vertex in transverse plane

        // members
        AliAODEvent*            fAOD;           //! input event
        AliAODTrack*						fTrack;					//! AOD track
        Double_t 								fTrackPt;				//! track pT
        Double_t 								fTrackPhi;				//! track phi
        Double_t 								fTrackEta;				//! track eta
        Int_t 									fNumV0s;					//! number of V0s in given event
        AliAODv0*			 					fV0;					//! V0 candidate
        Bool_t 									fV0candK0s;		//! Is V0 a K0s candidate ?
        Bool_t 									fV0candLambda;		//! Is V0 a Lambda or Anti-Lambda (ALambda) candidate ?
        Bool_t 									fV0candALambda;		//! Is V0 a Lambda or Anti-Lambda (ALambda) candidate ?
        Double_t 								fV0MaxMassK0s;		//! Upper limit of K0s inv. mass window
        Double_t 								fV0MinMassK0s;		//! Lower limit of K0s inv. mass window
        Double_t 								fV0MaxMassLambda;		//! Upper limit of Lambda inv. mass window
        Double_t 								fV0MinMassLambda;		//! Upper limit of Lambda inv. mass window
        Short_t									fCentBinIndex;					//! event centrality bin index indicator
        Float_t									fCentPercentile;					//! event centrality bin index indicator
        Short_t									fPtBinIndex;		//! track pT bin index indicator
        TComplex								fQvec2;					//! complex flow vector Q (n = 2)
        TComplex								fQvec3;					//! complex flow vector Q (n = 3)
        TComplex								fQvec4;					//! complex flow vector Q (n = 4)
        TComplex								fQvec5;					//! complex flow vector Q (n = 5)
        TComplex								fQvec2Gap00P;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap00N;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap04P;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap04N;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap08P;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap08N;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap10P;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fQvec2Gap10N;					//! complex flow vector Q (n = 5) with eta gap
        TComplex								fPvec2[fNumPtBins];	//! complex vector p (n = 2) for pT-differential flow 
        TComplex								fPvec2Gap00P[fNumPtBins];	//! complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec2Gap04P[fNumPtBins];	//! complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec2Gap08P[fNumPtBins];	//! complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec2Gap10P[fNumPtBins];	//! complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec3[fNumPtBins];	//! complex vector p (n = 2) for pT-differential flow 
        
        
        TList*                  fOutputList;    //! main output list
        TList*									fOutputListQA;	//! QA output list

        //std histos
        TH1D*										fEventMult;			 //! selected events multiplicity distribution
        TH1D*   					      fCentralityDis;     //! event centrality distribution
        TH2D*							      fCentSPDvsV0M;      //! V0M vs SPD
        TH1D*										fMultTracksSelected; //! multiplicity of selected tracks in a given event
        TH2D*										fTracksPtCent;		//! selected tracks pT vs event centrality
        TH1D*                   fTracksPt;       //! selected tracks pT distribution
        TH1D*                   fTracksEta;      //! selected tracks eta distribution
        TH1D* 									fTracksPhi;			 //! selected tracks phi distribution
        TH1D*										fV0sMult;				//! multiplicity of V0s in selected events
        TProfile*								fRefCorTwo2;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v2
        TProfile*								fRefCorTwo3;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v3
        TProfile*								fRefCorTwo4;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v4
        TProfile*								fRefCorTwo5;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap00;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap04;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap08;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap10;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
				TProfile*								fDiffCorTwo2[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap00[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap04[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap08[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap10[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo3[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>


        // QA histos
        TH1D* 									fEventCounter;  //! event rejection tracker
        TH1D*										fV0sCounter;		//! V0s counter
        TH1D* 									fQAPVz;					//! PV z distance distribution
        TH1D*										fQANumTracks;		//! number of AOD tracks distribution
        TH1D*										fQATrackPt;			//! pT dist of all tracks in all events
        TH1D*										fQATrackEta;		//! eta dist of all tracks in all events
        TH1D*										fQATrackPhi;		//! phi dist of all tracks in all events
        TH1D*										fQATrackFilterMap;//! filter bit of all tracks


        AliAnalysisTaskFlowPID(const AliAnalysisTaskFlowPID&); // not implemented
        AliAnalysisTaskFlowPID& operator=(const AliAnalysisTaskFlowPID&); // not implemented

        ClassDef(AliAnalysisTaskFlowPID, 1);
};

#endif
