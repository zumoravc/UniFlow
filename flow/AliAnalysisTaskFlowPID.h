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
        // event & tracks setters
        void										SetAODAnalysis(Bool_t aod) { fAODAnalysis = aod; }
        void										SetPbPbAnalysis(Bool_t pbpb) { fPbPb = pbpb; }
        void										SetPeriod10h(Bool_t period) { fLHC10h = period; }
        void										SetCentFlag(Short_t flag) { fCentFlag = flag; }
        void										SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        void										SetTrackEtaMax(Double_t eta) { fTrackEtaMax = eta; }
        void										SetTrackPtMax(Double_t pt) { fTrackPtMax = pt; }
        void										SetTrackPtMin(Double_t pt) { fTrackPtMin = pt; }
        void										SetNumTPCclsMin(UShort_t tpcCls) { fNumTPCclsMin = tpcCls; }
        void										SetTrackFilterBit(UInt_t filter) { fTrackFilterBit = filter; }
        void										SetDiffFlow(Bool_t diff) { fDiffFlow = diff; }
        void										SetPID(Bool_t pid) { fPID = pid; }
				// V0s setters
        void										SetV0sOnFly(Bool_t onFly) { fCutV0onFly = onFly; }
        void										SetV0sTPCRefit(Bool_t refit) { fCutV0refitTPC = refit; }
        void										SetV0sRejectKinks(Bool_t reject) { fCutV0rejectKinks = reject; }
        void										SetV0sDCAPVMin(Double_t dca) { fCutV0MinDCAtoPV = dca; }
        void										SetV0sDCAPVMax(Double_t dca) { fCutV0MaxDCAtoPV = dca; }
        void										SetV0sDCADaughtersMax(Double_t dca) { fCutV0MaxDCADaughters = dca; }
        void										SetV0sDecayRadiusMin(Double_t radius) { fCutV0MinDecayRadius = radius; }
        void										SetV0sDecayRadiusMax(Double_t radius) { fCutV0MaxDecayRadius = radius; }
        void										SetV0sDaughterPtMin(Double_t pt) { fCutV0DaughterPtMin = pt; }
        void										SetV0sDaughterEtaMax(Double_t eta) { fCutV0DaughterEtaMax = eta; }
        void										SetV0sMotherEtaMax(Double_t eta) { fCutV0MotherEtaMax = eta; }
        void                                        SetV0sMotherRapMax(Double_t rap) { fCutV0MotherRapMax = rap; }
        void                                        SetV0sMotherPtMin(Double_t pt) { fCutV0MotherPtMin = pt; }
        void										SetV0sMotherPtMax(Double_t pt) { fCutV0MotherPtMax = pt; }
        void										SetV0sK0sCPAMin(Double_t cpa) { fCutV0MinCPAK0s = cpa; }
        void										SetV0sLambdaCPAMin(Double_t cpa) { fCutV0MinCPALambda = cpa; }
        void										SetV0sK0sNumTauMax(Double_t nTau) { fCutV0NumTauK0sMax = nTau; }
        void										SetV0sLambdaNumTauMax(Double_t nTau) { fCutV0NumTauLambdaMax = nTau; }
        void										SetV0sK0sArmenterosAlphaMin(Double_t alpha) { fCutV0K0sArmenterosAlphaMin = alpha; }
        void                                        SetV0sProtonNumSigmaMax(Double_t nSigma) { fCutV0ProtonNumSigmaMax = nSigma; }
        void										SetV0sProtonPIDPtMax(Double_t pt) { fCutV0ProtonPIDPtMax = pt; }


        const static Int_t 			fNumPtBins = 10;			// number of pT bins used for pT-differential flow
        static Double_t					fPtBinEdges[fNumPtBins+1];				// pointer for array of pT bin edges
        const static Int_t      fNumMinvFlowBinsK0s = 12;  // number of inv. mass bin for differential flow plots (K0s)
        static Double_t         fMinvFlowBinEdgesK0s[fNumMinvFlowBinsK0s+1]; // pointer to array of Minv bin edges (K0s)
        const static Int_t      fNumMinvFlowBinsLambda = 10;  // number of inv. mass bin for differential flow plots ((A)Lambda)
        static Double_t         fMinvFlowBinEdgesLambda[fNumMinvFlowBinsLambda+1]; // pointer to array of Minv bin edges ((A)Lambda)
        const static Int_t 			fNumCentBins = 9;			// number of centrality bins used for pT-differential flow (so far independently of reference flow)
        static Double_t					fCentBinEdges[fNumCentBins+1];				// pointer for array of pT bin edges

    private:
        Bool_t                  IsEventSelected(const AliAODEvent* event);
				Bool_t                  IsTrackSelected(const AliAODTrack* track);
                void            FilterTracks();
                void            FilterV0s();

				void                  IsV0aK0s(const AliAODv0* v0);
				void                 IsV0aLambda(const AliAODv0* v0);
				Bool_t                  IsV0Selected(const AliAODv0* v0);
                Bool_t          AreRefFlowVectorsFilled(const Float_t dEtaGap = -1, const Short_t iHarm = -1);
                void            FillRefFlowVectors(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2);
                void            EstimateRefCumulant(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2, TProfile* profile = 0x0);
                void            EstimateV0Cumulant(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2, TProfile2D* profilePos = 0x0, TProfile2D* profileNeg = 0x0);
                void            EstimateV0CumulantWithMass(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2, TProfile2D* profilePos = 0x0, TProfile2D* profileNeg = 0x0);

				void                    EventQA(const AliAODEvent* event);
				void            FillV0sQA(const AliAODv0* v0, const Short_t iQAindex);
       	void                    EstimateCentrality(AliVEvent* ev);
    
	    	Double_t                GetWDist(const AliVVertex* v0, const AliVVertex* v1); 
        Bool_t                  plpMV(const AliVEvent *event);

        Short_t                 GetPtBinIndex(const Double_t dPt);
        Short_t                 GetMinvFlowBinIndexK0s(const Double_t dMass);
        Short_t                 GetMinvFlowBinIndexLambda(const Double_t dMass);

        //cuts & selection: event & track
        Bool_t									fAODAnalysis;		// is AOD analysis?
        Bool_t									fPbPb;					// is PbPb analysis?
				Bool_t       						fLHC10h;        // flag to LHC10h data?
				Short_t									fCentFlag;			// centrality flag
        Double_t 								fPVtxCutZ; 			// (cm) PV z cut
        Double_t 								fTrackEtaMax; 	// () Maximum pseudorapidity range
        Double_t								fTrackPtMax;		// (GeV/c) Maximal track pT
        Double_t								fTrackPtMin;		// (GeV/c) Minimal track pT
        UShort_t								fNumTPCclsMin;	// () Minimal number of TPC clusters used for track reconstruction
        UInt_t									fTrackFilterBit;// Required track filter bit 
        Bool_t 									fDiffFlow;			// Do differential flow ? (or reference only)
        Bool_t 									fPID;						// Do PID (so far V0s) ? 
				//cuts & selection: V0 reconstruction
				Bool_t 									fCutV0onFly;		// V0 reconstruction method: is On-the-fly? (or offline)
				Bool_t									fCutV0refitTPC; // Check TPC refit of V0 daughters ?
				Bool_t									fCutV0rejectKinks; // Reject Kink V0 daughter tracks ?
				Double_t                fCutV0MinDCAtoPV;   // (cm) min DCA of V0 daughter to PV
        Double_t								fCutV0MaxDCAtoPV;	// (cm) max DCA of V0 daughter to PV
				Double_t								fCutV0MaxDCADaughters;	// (cm) max DCA of V0 daughters among themselves
        Double_t                fCutV0MinDecayRadius; // (cm) min distance of secondary vertex from z-axis in transverse plane
				Double_t								fCutV0MaxDecayRadius; // (cm) max distance of secondary vertex from z-axis in transverse plane
        Double_t                fCutV0DaughterPtMin; // (GeV/c) minimum pT of V0 daughters
        Double_t                fCutV0DaughterEtaMax; // () max value of Eta of V0 daughters
        Double_t                fCutV0MotherEtaMax; // () max eta value of V0 mother
        Double_t                fCutV0MotherRapMax; // () max rapidity value of V0 mother
        Double_t                fCutV0MotherPtMin; // () min transverse momentum value of V0 mother
        Double_t                fCutV0MotherPtMax; // () max transverse momentum value of V0 mother
        Double_t                fCutV0MinCPAK0s;    // () min cosine of pointing angle of K0s candidate to PV
        Double_t                fCutV0MinCPALambda; // () min cosine of pointing angle of K0s candidate to PV
        Double_t                fCutV0NumTauK0sMax; // (c*tau) max number of c*tau (K0s)
        Double_t                fCutV0NumTauLambdaMax; // (c*tau) max number of c*tau ((A)Lambda)
        Double_t								fCutV0K0sArmenterosAlphaMin; // (alpha) max Armenteros alpha for K0s
        Double_t                                fCutV0ProtonNumSigmaMax;    // (sigmaTPC) --- both MUST be on --- max number of TPC sigma for proton PID (Lambda candidates)
        Double_t								fCutV0ProtonPIDPtMax;	// (GeV/c) --- both MUST be on --- max pT of proton for PID (Lambda candidates) - only protons with smaller will be checked for num sigma TPC
        // members
        AliAODEvent*            fAOD;           //! input event
        AliPIDResponse*					fPIDResponse;		//! PID response
        AliTPCPIDResponse               fTPCPIDResponse;       // TPC PID response
        Float_t                         fEtaCutFlag;    // value of eta cut during flow vectors filling  (actual value)
        Short_t                         fHarmFlag;    // value of harmonics during flow vectors filling  (actual value)
        AliAODTrack*						fTrack;					//! AOD track
        Double_t 								fTrackPt;				// track pT
        Double_t 								fTrackPhi;				// track phi
        Double_t 								fTrackEta;				// track eta
        Int_t 									fNumV0s;					// number of V0s in given event
        AliAODv0*			 					fV0;					//! V0 candidate
        Bool_t 									fV0candK0s;		// Is V0 a K0s candidate ? flag 
        Bool_t                                  fV0candLambda;      // Is V0 a Lambda  candidate ? flag
        Bool_t 									fV0candALambda;		// Is V0 a Anti-Lambda (ALambda) candidate ? flag
        Double_t 								fV0MaxMassK0s;		// Upper limit of K0s inv. mass window
        Double_t 								fV0MinMassK0s;		// Lower limit of K0s inv. mass window
        Double_t 								fV0MaxMassLambda;		// Upper limit of Lambda inv. mass window
        Double_t 								fV0MinMassLambda;		// Upper limit of Lambda inv. mass window
        Short_t									fCentBinIndex;					// event centrality bin index indicator
        Float_t									fCentPercentile;					// event centrality bin index indicator
        Short_t                 fPtBinIndex;        // track pT bin index indicator
        Short_t									fMinvFlowBinIndex;		// track pT bin index indicator
        
        TComplex                                fVecRefPos;                 // complex flow vector Q for RFPs in positive eta
        TComplex                                fVecRefNeg;                 // complex flow vector Q for RFPs in negative eta
        Int_t                                   fCountRefPos;           // counter for number of RFPs in positive eta
        Int_t                                   fCountRefNeg;           // counter for number of RFPs in positive eta

        TComplex								fQvec2;					// complex flow vector Q (n = 2)
        TComplex								fQvec3;					// complex flow vector Q (n = 3)
        TComplex								fQvec4;					// complex flow vector Q (n = 4)
        TComplex								fQvec5;					// complex flow vector Q (n = 5)
        TComplex								fQvec2Gap00P;				// complex flow vector Q (n = 2) with eta gap
        TComplex								fQvec2Gap00N;					// complex flow vector Q (n = 2) with eta gap
        TComplex                fQvec2Gap04P;                   // complex flow vector Q (n = 2) with eta gap
        TComplex                fQvec2Gap04N;                   // complex flow vector Q (n = 2) with eta gap
        TComplex                fQvec2Gap08P;                   // complex flow vector Q (n = 2) with eta gap
        TComplex                fQvec2Gap08N;                   // complex flow vector Q (n = 2) with eta gap
        TComplex								fQvec2Gap09P;					// complex flow vector Q (n = 2) with eta gap
        TComplex								fQvec2Gap09N;					// complex flow vector Q (n = 2) with eta gap
        TComplex                fQvec2Gap10P;                   // complex flow vector Q (n = 2) with eta gap
        TComplex								fQvec2Gap10N;					// complex flow vector Q (n = 2) with eta gap
        TComplex								fPvec2[fNumPtBins];	// complex vector p (n = 2) for pT-differential flow 
        TComplex								fPvec2Gap00P[fNumPtBins];	// complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec2Gap04P[fNumPtBins];	// complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec2Gap08P[fNumPtBins];	// complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec2Gap10P[fNumPtBins];	// complex vector p (n = 2) for pT-differential flow with eta gap
        TComplex								fPvec3[fNumPtBins];	// complex vector p (n = 2) for pT-differential flow 

        TComplex                fVvec2Gap00P_K0s[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap00N_K0s[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap09P_K0s[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap09N_K0s[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap00P_Lambda[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap00N_Lambda[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap09P_Lambda[fNumPtBins][fNumMinvFlowBinsK0s]; //
        TComplex                fVvec2Gap09N_Lambda[fNumPtBins][fNumMinvFlowBinsK0s]; //
        
        TClonesArray            fArrTracksFiltered; // container for filtered tracks
        TClonesArray            fArrV0sK0sFiltered; // container for filtered V0 K0s candidates
        TClonesArray            fArrV0sLambdaFiltered; // container for filtered V0 Lambda candidates
        TClonesArray            fArrV0sALambdaFiltered; // container for filtered V0 ALambda candidates

        TList*                  fOutList;    //! main output list
        TList*                  fOutListV0s;    //! main output list
        TList*                  fOutListQA;	//! QA output list

        // testing / temporary histos
        TH1D*                   fTestTracksPt;      //! test histo
        TH1D*                   fTestTracksMult;      //! test histo

        //std histos
        TH1D*										fEventMult;			 //! selected events multiplicity distribution
        TH1D*   					      fCentralityDis;     //! event centrality distribution
        TH2D*							      fCentSPDvsV0M;      //! V0M vs SPD
        TH1D*										fMultTracksSelected; //! multiplicity of selected tracks in a given event
        TH2D*										fTracksPtCent;		//! selected tracks pT vs event centrality
        TH1D*                   fTracksPt;       //! selected tracks pT distribution
        TH1D*                   fTracksEta;      //! selected tracks eta distribution
        TH1D* 									fTracksPhi;			 //! selected tracks phi distribution
        TH1D* 									fTracksCharge;			 //! selected tracks charge distribution
        TProfile*								fRefCorTwo2;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v2
        TProfile*								fRefCorTwo3;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v3
        TProfile*								fRefCorTwo4;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v4
        TProfile*								fRefCorTwo5;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap00;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap04;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*               fRefCorTwo2Gap08;                //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*                               fRefCorTwo2Gap09;                //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap09_test;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
        TProfile*								fRefCorTwo2Gap10;			 	 //! event averaged 2-particle correlation for reference flow <<2>> v5
				TProfile*								fDiffCorTwo2[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap00[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap04[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap08[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo2Gap10[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>
				TProfile*								fDiffCorTwo3[fNumCentBins];			 //! event averaged 2-particle correlation for differential flow <<2'>>

				TProfile2D* 			fV0sDiffTwo2Gap00P_K0s[fNumCentBins];      //! selected K0s candidates Minv, pT v2 profile
				TProfile2D*             fV0sDiffTwo2Gap00N_K0s[fNumCentBins];      //! selected K0s candidates Minv, pT v2 profile
                TProfile2D*             fV0sDiffTwo2Gap09P_K0s_test[fNumCentBins];      //! selected K0s candidates Minv, pT v2 profile
                TProfile2D*             fV0sDiffTwo2Gap09N_K0s_test[fNumCentBins];
  				TProfile2D*             fV0sDiffTwo2Gap09P_K0s[fNumCentBins];      //! selected K0s candidates Minv, pT v2 profile
				TProfile2D*             fV0sDiffTwo2Gap09N_K0s[fNumCentBins];      //! selected K0s candidates Minv, pT v2 profile
				TProfile2D*             fV0sDiffTwo2Gap00P_Lambda[fNumCentBins];      //! selected (Anti)Lambda candidates Minv, pT v2 profile
				TProfile2D*             fV0sDiffTwo2Gap00N_Lambda[fNumCentBins];      //! selected (Anti)Lambda candidates Minv, pT v2 profile
				TProfile2D*             fV0sDiffTwo2Gap09P_Lambda[fNumCentBins];      //! selected (Anti)Lambda candidates Minv, pT v2 profile
				TProfile2D*             fV0sDiffTwo2Gap09N_Lambda[fNumCentBins];      //! selected (Anti)Lambda candidates Minv, pT v2 profile

				// V0s histos
        TH1D*										fV0sMult;				//! multiplicity of V0s in selected events
        TH1D*                                       fV0sInvMassK0sGap00;        //! selected K0s inv. mass distribution (pT & cent integrated)
        TH1D*                                       fV0sInvMassK0sGap09;        //! selected K0s inv. mass distribution (pT & cent integrated)
        TH1D*										fV0sInvMassK0sGap09_test;		//! selected K0s inv. mass distribution (pT & cent integrated)
        TH1D*                                       fV0sInvMassLambdaGap00;     //! selected Lambda candidates inv. mass distribution (pT & cent integrated)
        TH1D*										fV0sInvMassLambdaGap09;		//! selected Lambda candidates inv. mass distribution (pT & cent integrated)
        TH2D*										fV0sK0sGap00[fNumCentBins];							//! selected K0s distribution (InvMass, pT)
        TH2D*                                       fV0sK0sGap09[fNumCentBins];                         //! selected K0s distribution (InvMass, pT)
        TH2D*										fV0sK0sGap09_test[fNumCentBins];							//! selected K0s distribution (InvMass, pT)
        TH2D*										fV0sLambdaGap00[fNumCentBins];							//! selected K0s distribution (InvMass, pT)
        TH2D*										fV0sLambdaGap09[fNumCentBins];							//! selected K0s distribution (InvMass, pT)
        // QA histos // index 0: before / 1: after cuts
        TH1D* 									fEventCounter;  //! event rejection tracker
        TH1D* 									fQAPVz;					//! PV z distance distribution
        TH1D*										fQANumTracks;		//! number of AOD tracks distribution
        TH1D*										fQATrackPt;			//! pT dist of all tracks in all events
        TH1D*										fQATrackEta;		//! eta dist of all tracks in all events
        TH1D*										fQATrackPhi;		//! phi dist of all tracks in all events
        TH1D*										fQATrackFilterMap;//! filter bit of all tracks
        	// QA V0s
        const static Short_t                        fQAV0sNumSteps = 2;        // number of various steps (0 before cuts / 1 after cuts / 2 testing) 
        TH1D*										fQAV0sCounter;		//! V0s counter
        TH1D*										fQAV0sCounterK0s;		//! K0s counter
        TH1D*										fQAV0sCounterLambda;		//! Lambda counter
        TH1D*										fQAV0sRecoMethod[fQAV0sNumSteps];	//! offline/online V0 reconstruction method
        TH1D*										fQAV0sTPCRefit[fQAV0sNumSteps];	//! TPC refit true/false
        TH1D*										fQAV0sKinks[fQAV0sNumSteps];	//! V0 kinks true/false
        TH1D*										fQAV0sDCAtoPV[fQAV0sNumSteps];	//! V0 DCA to PV
        TH1D*										fQAV0sDCADaughters[fQAV0sNumSteps];	//! DCA between V0 daughters
        TH1D*										fQAV0sDecayRadius[fQAV0sNumSteps];	//! Distance between PV and Secondary vertex in transverse plane
        TH1D*                                       fQAV0sDaughterPt[fQAV0sNumSteps];    //! pT dist of V0 daughters
        TH1D*										fQAV0sDaughterPhi[fQAV0sNumSteps];	//! pT dist of V0 daughters
        TH1D*                                       fQAV0sDaughterEta[fQAV0sNumSteps];   //! pseudorapidity dist of V0 daughters
        TH2D*										fQAV0sDaughterTPCdEdxPt[fQAV0sNumSteps];	//! TPC dEdx vs pT of V0 daughters
        TH1D*                                       fQAV0sMotherPt[fQAV0sNumSteps];  //! pT dist of V0s
        TH1D*										fQAV0sMotherPhi[fQAV0sNumSteps];	//! azimuthal dist of V0s
        TH1D*                                       fQAV0sMotherEta[fQAV0sNumSteps]; //! pseudorapidity dist of V0s
        TH1D*                                       fQAV0sMotherRapK0s[fQAV0sNumSteps];  //! rapidity dist of V0s (K0s mass hypothesis)
        TH1D*                                       fQAV0sMotherRapLambda[fQAV0sNumSteps];   //! rapidity dist of V0s (Lambda mass hypothesis)
        TH1D*                                       fQAV0sInvMassK0s[fQAV0sNumSteps];    //! inv. mass dist of V0s (K0s mass hypothesis)
        TH1D*										fQAV0sInvMassLambda[fQAV0sNumSteps];	//! inv. mass dist of V0s ((A)Lambda mass hypothesis)
        TH1D*										fQAV0sCPAK0s[fQAV0sNumSteps];	//! cosine of pointing angle of K0s candidates
        TH1D*										fQAV0sCPALambda[fQAV0sNumSteps];	//! cosine of pointing angle of Lambda candidates
        TH1D*										fQAV0sNumTauK0s[fQAV0sNumSteps];	//! number of c*tau of K0s candidates
        TH1D*										fQAV0sNumTauLambda[fQAV0sNumSteps];	//! number of c*tau of Lambda candidates
        TH2D*										fQAV0sArmenterosK0s[fQAV0sNumSteps];	//! Armenteros-Podolanski plot for K0s candidates
        TH2D*										fQAV0sArmenterosLambda[fQAV0sNumSteps];	//! Armenteros-Podolanski plot for K0s candidates
        TH2D*                                       fQAV0sProtonNumSigmaPtLambda[fQAV0sNumSteps];  //! number of TPC sigmas and pT of proton (Lambda candidates)
        

        AliAnalysisTaskFlowPID(const AliAnalysisTaskFlowPID&); // not implemented
        AliAnalysisTaskFlowPID& operator=(const AliAnalysisTaskFlowPID&); // not implemented

        ClassDef(AliAnalysisTaskFlowPID, 3);
};

#endif
