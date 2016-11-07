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
        // event setters
        void					SetAODAnalysis(Bool_t aod) { fAODAnalysis = aod; }
        void					SetPbPbAnalysis(Bool_t pbpb) { fPbPb = pbpb; }
        void					SetPeriod10h(Bool_t period) { fLHC10h = period; }
        void					SetCentFlag(Short_t flag) { fCentFlag = flag; }
        void					SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        void                    SetSampling(Bool_t sample) { fSampling = sample; }
        void					SetDiffFlow(Bool_t diff) { fDiffFlow = diff; }
        void					SetPID(Bool_t pid) { fPID = pid; }
        // track setters
        void                    SetTrackEtaMax(Double_t eta) { fTrackEtaMax = eta; }
        void                    SetTrackPtMax(Double_t pt) { fTrackPtMax = pt; }
        void                    SetTrackPtMin(Double_t pt) { fTrackPtMin = pt; }
        void                    SetNumTPCclsMin(UShort_t tpcCls) { fNumTPCclsMin = tpcCls; }
        void                    SetTrackFilterBit(UInt_t filter) { fTrackFilterBit = filter; }
        void                    SetPionNumSigmasMax(Double_t numSigmas) { fCutPionNumSigmaMax = numSigmas; }
        void                    SetKaonNumSigmasMax(Double_t numSigmas) { fCutKaonNumSigmaMax = numSigmas; }
        void                    SetProtonNumSigmasMax(Double_t numSigmas) { fCutProtonNumSigmaMax = numSigmas; }
        // V0s setters
        void					SetV0sOnFly(Bool_t onFly) { fCutV0onFly = onFly; }
        void					SetV0sTPCRefit(Bool_t refit) { fCutV0refitTPC = refit; }
        void					SetV0sRejectKinks(Bool_t reject) { fCutV0rejectKinks = reject; }
        void					SetV0sDCAPVMin(Double_t dca) { fCutV0MinDCAtoPV = dca; }
        void					SetV0sDCAPVMax(Double_t dca) { fCutV0MaxDCAtoPV = dca; }
        void					SetV0sDCADaughtersMax(Double_t dca) { fCutV0MaxDCADaughters = dca; }
        void					SetV0sDecayRadiusMin(Double_t radius) { fCutV0MinDecayRadius = radius; }
        void					SetV0sDecayRadiusMax(Double_t radius) { fCutV0MaxDecayRadius = radius; }
        void					SetV0sDaughterPtMin(Double_t pt) { fCutV0DaughterPtMin = pt; }
        void					SetV0sDaughterEtaMax(Double_t eta) { fCutV0DaughterEtaMax = eta; }
        void					SetV0sMotherEtaMax(Double_t eta) { fCutV0MotherEtaMax = eta; }
        void                    SetV0sMotherRapMax(Double_t rap) { fCutV0MotherRapMax = rap; }
        void                    SetV0sMotherPtMin(Double_t pt) { fCutV0MotherPtMin = pt; }
        void					SetV0sMotherPtMax(Double_t pt) { fCutV0MotherPtMax = pt; }
        void					SetV0sK0sCPAMin(Double_t cpa) { fCutV0MinCPAK0s = cpa; }
        void					SetV0sLambdaCPAMin(Double_t cpa) { fCutV0MinCPALambda = cpa; }
        void					SetV0sK0sNumTauMax(Double_t nTau) { fCutV0NumTauK0sMax = nTau; }
        void					SetV0sLambdaNumTauMax(Double_t nTau) { fCutV0NumTauLambdaMax = nTau; }
        void					SetV0sK0sArmenterosAlphaMin(Double_t alpha) { fCutV0K0sArmenterosAlphaMin = alpha; }
        void                    SetV0sProtonNumSigmaMax(Double_t nSigma) { fCutV0ProtonNumSigmaMax = nSigma; }
        void					SetV0sProtonPIDPtMax(Double_t pt) { fCutV0ProtonPIDPtMax = pt; }


        const static Int_t 		fNumPtBins = 22;			// number of pT bins used for pT-differential flow
        static Double_t			fPtBinEdges[fNumPtBins+1];				// pointer for array of pT bin edges
        const static Int_t      fNumMinvFlowBinsK0s = 12;  // number of inv. mass bin for differential flow plots (K0s)
        static Double_t         fMinvFlowBinEdgesK0s[fNumMinvFlowBinsK0s+1]; // pointer to array of Minv bin edges (K0s)
        const static Int_t      fNumMinvFlowBinsLambda = 11;  // number of inv. mass bin for differential flow plots ((A)Lambda)
        static Double_t         fMinvFlowBinEdgesLambda[fNumMinvFlowBinsLambda+1]; // pointer to array of Minv bin edges ((A)Lambda)
        const static Int_t 		fNumCentBins = 9;			// number of centrality bins used for pT-differential flow (so far independently of reference flow)
        static Double_t			fCentBinEdges[fNumCentBins+1];				// pointer for array of pT bin edges
        const static Int_t      fNumHarmonics = 1; // number of harmonics
        static Int_t            fHarmonics[fNumHarmonics]; // values of used harmonics
        const static Int_t      fNumEtaGap = 8; // number of harmonics
        static Double_t         fEtaGap[fNumEtaGap]; // values of used harmonics
    private:
        Bool_t                  IsEventSelected(const AliAODEvent* event);
		Bool_t                  IsTrackSelected(const AliAODTrack* track);
        Bool_t                  IsTrackPion(const AliAODTrack* track); 
        Bool_t                  IsTrackKaon(const AliAODTrack* track); 
        Bool_t                  IsTrackProton(const AliAODTrack* track); 
		void                    IsV0aK0s(const AliAODv0* v0);
		void                    IsV0aLambda(const AliAODv0* v0);
		Bool_t                  IsV0Selected(const AliAODv0* v0);
        
        void                    FilterTracks(); // filter all input tracks (including PID tracks for pi,K,p)
        void                    FilterPIDTracks(); // obsolete / not used
        void                    FilterV0s();

        Bool_t                  AreRefFlowVectorsFilled(const Float_t dEtaGap = -1, const Short_t iHarm = -1);
        void                    FillRefFlowVectors(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2);
        void                    EstimateRefCumulant(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2, TProfile* profile = 0x0);
        void                    EstimateRefPtDiffCumulant(const Float_t dEtaGap = 0.9, const Short_t iHarm = 2, TProfile* profilePos = 0x0, TProfile* profileNeg = 0x0);
        void                    EstimateV0Cumulant(const Short_t iEtaGapIndex = 0, const Short_t iHarmonicsIndex = 0, const Short_t iSampleIndex = 0);
                
		void                    FillEventQA(const AliAODEvent* event, const Short_t iQAindex);
        void                    FillTrackQA(const AliAODTrack* track, const Short_t iQAindex);
        void                    FillV0sQA(const AliAODv0* v0, const Short_t iQAindex);
       	void                    EstimateCentrality(AliVEvent* ev);
    
    	Double_t                GetWDist(const AliVVertex* v0, const AliVVertex* v1); 
        Bool_t                  plpMV(const AliVEvent *event);

        Short_t                 GetPtBinIndex(const Double_t dPt);
        Short_t                 GetMinvFlowBinIndexK0s(const Double_t dMass);
        Short_t                 GetMinvFlowBinIndexLambda(const Double_t dMass);

        //cuts & selection: analysis
        Bool_t                  fAODAnalysis;       // is AOD analysis?
        Bool_t                  fPbPb;                  // is PbPb analysis?
        Bool_t                  fLHC10h;        // flag to LHC10h data?
        Short_t                 fCentFlag;          // centrality flag
        Bool_t                  fSampling;      // Do random sampling ? (estimation of vn stat. uncertanity)
        Bool_t                  fDiffFlow;          // Do differential flow ? (or reference only)
        Bool_t                  fPID;                       // Do PID (so far V0s) ? 
        //cuts & selection: events & tracks
        Double_t                fPVtxCutZ;          // (cm) PV z cut
        UInt_t                  fTrackFilterBit; // tracks filter bit 
        UShort_t                fNumTPCclsMin;  // () Minimal number of TPC clusters used for track reconstruction
        Double_t                fTrackEtaMax;   // () Maximum pseudorapidity range
        Double_t                fTrackPtMax;        // (GeV/c) Maximal track pT
        Double_t                fTrackPtMin;        // (GeV/c) Minimal track pT
        Double_t                fCutPionNumSigmaMax;
        Double_t                fCutKaonNumSigmaMax;
        Double_t                fCutProtonNumSigmaMax;
        //cuts & selection: V0 reconstruction
		Bool_t 					fCutV0onFly;		// V0 reconstruction method: is On-the-fly? (or offline)
		Bool_t					fCutV0refitTPC; // Check TPC refit of V0 daughters ?
		Bool_t					fCutV0rejectKinks; // Reject Kink V0 daughter tracks ?
		Double_t                fCutV0MinDCAtoPV;   // (cm) min DCA of V0 daughter to PV
        Double_t				fCutV0MaxDCAtoPV;	// (cm) max DCA of V0 daughter to PV
		Double_t				fCutV0MaxDCADaughters;	// (cm) max DCA of V0 daughters among themselves
        Double_t                fCutV0MinDecayRadius; // (cm) min distance of secondary vertex from z-axis in transverse plane
		Double_t				fCutV0MaxDecayRadius; // (cm) max distance of secondary vertex from z-axis in transverse plane
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
        Double_t				fCutV0K0sArmenterosAlphaMin; // (alpha) max Armenteros alpha for K0s
        Double_t                fCutV0ProtonNumSigmaMax;    // (sigmaTPC) --- both MUST be on --- max number of TPC sigma for proton PID (Lambda candidates)
        Double_t				fCutV0ProtonPIDPtMax;	// (GeV/c) --- both MUST be on --- max pT of proton for PID (Lambda candidates) - only protons with smaller will be checked for num sigma TPC
        
        // members
        AliAODEvent*            fAOD;           //! input event
        AliPIDResponse*			fPIDResponse;		//! PID response
        AliTPCPIDResponse       fTPCPIDResponse;       // TPC PID response
        Float_t                 fEtaCutFlag;    // value of eta cut during flow vectors filling  (actual value)
        Short_t                 fHarmFlag;    // value of harmonics during flow vectors filling  (actual value)

        TClonesArray            fArrTracksFiltered; // container for filtered tracks
        TClonesArray            fArrPionFiltered; // container for filtered pions
        TClonesArray            fArrKaonFiltered; // container for filtered kaons
        TClonesArray            fArrProtonFiltered; // container for filtered protons
        TClonesArray            fArrV0sK0sFiltered; // container for filtered V0 K0s candidates
        TClonesArray            fArrV0sLambdaFiltered; // container for filtered V0 Lambda candidates
        TClonesArray            fArrV0sALambdaFiltered; // container for filtered V0 ALambda candidates

        Short_t                 fCentBinIndex;                  // event centrality bin index indicator
        Float_t                 fCentPercentile;                    // event centrality bin index indicator
        Short_t                 fPtBinIndex;        // track pT bin index indicator
        Short_t                 fMinvFlowBinIndex;      // track pT bin index indicator
        Short_t                 fSampleBinIndex;         // sampling bin index indicator
        
        TComplex                fVecRefPos;                 // complex flow vector Q for RFPs in positive eta (or all if no eta gap)
        TComplex                fVecRefNeg;                 // complex flow vector Q for RFPs in negative eta
        Int_t                   fCountRefPos;           // counter for number of RFPs in positive eta (or all if no eta gap)
        Int_t                   fCountRefNeg;           // counter for number of RFPs in positive eta

        Bool_t 					fV0candK0s;		// Is V0 a K0s candidate ? flag 
        Bool_t                  fV0candLambda;      // Is V0 a Lambda  candidate ? flag
        Bool_t 					fV0candALambda;		// Is V0 a Anti-Lambda (ALambda) candidate ? flag
        Double_t 				fV0MaxMassK0s;		// Upper limit of K0s inv. mass window
        Double_t 				fV0MinMassK0s;		// Lower limit of K0s inv. mass window
        Double_t 				fV0MaxMassLambda;		// Upper limit of Lambda inv. mass window
        Double_t 				fV0MinMassLambda;		// Upper limit of Lambda inv. mass window
    
        TList*                  fOutList;    //! main output list
        TList*                  fOutListV0s;    //! main output list
        TList*                  fOutListQA;	//! QA output list

        //std histos
        TH1D*					fEventMult;			 //! selected events multiplicity distribution
        TH1D*   				fCentralityDis;     //! event centrality distribution
        TH2D*					fCentSPDvsV0M;      //! V0M vs SPD
        TH1D*					fMultTracksSelected; //! multiplicity of selected tracks in a given event
        TH2D*					fTracksPtCent;		//! selected tracks pT vs event centrality
        TH1D*                   fTracksPt;       //! selected tracks pT distribution
        TH1D*                   fTracksEta;      //! selected tracks eta distribution
        TH1D* 				    fTracksPhi;			 //! selected tracks phi distribution
        TH1D* 					fTracksCharge;			 //! selected tracks charge distribution
        

        const static Int_t      fNumSampleBins = 5; // number of sampling bins 
        TProfile*				fTracksRefTwo[fNumHarmonics][fNumEtaGap][fNumSampleBins];			 	 //! event averaged 2-particle correlation for reference flow <<2>> v2
        TProfile*               fTracksDiffTwoPos[fNumCentBins][fNumHarmonics][fNumEtaGap][fNumSampleBins];          //! event averaged 2-particle correlation for differential flow <<2'>>
        TProfile*               fTracksDiffTwoNeg[fNumCentBins][fNumHarmonics][fNumEtaGap][fNumSampleBins];          //! event averaged 2-particle correlation for differential flow <<2'>>
        
        TProfile2D*             fV0sDiffTwoPos_K0s[fNumCentBins][fNumHarmonics][fNumEtaGap][fNumSampleBins];      //! selected K0s candidates Minv, pT v2 profile
        TProfile2D*             fV0sDiffTwoNeg_K0s[fNumCentBins][fNumHarmonics][fNumEtaGap][fNumSampleBins];      //! selected K0s candidates Minv, pT v2 profile
        TProfile2D*             fV0sDiffTwoPos_Lambda[fNumCentBins][fNumHarmonics][fNumEtaGap][fNumSampleBins];      //! selected (Anti)Lambda candidates Minv, pT v2 profile
        TProfile2D*             fV0sDiffTwoNeg_Lambda[fNumCentBins][fNumHarmonics][fNumEtaGap][fNumSampleBins];      //! selected (Anti)Lambda candidates Minv, pT v2 profile

		// V0s histos
        TH1D*                   fV0sInvMassK0s[fNumEtaGap]; //!
        TH1D*                   fV0sInvMassLambda[fNumEtaGap]; //!
        TH2D*                   fV0sPtInvMassK0s[fNumCentBins][fNumHarmonics][fNumEtaGap];                         //! selected K0s distribution (InvMass, pT)
        TH2D*                   fV0sPtInvMassLambda[fNumCentBins][fNumHarmonics][fNumEtaGap];                         //! selected K0s distribution (InvMass, pT)
        
        // QA histos // index 0: before / 1: after cuts
        const static Short_t    fQANumSteps = 2;        // number of various steps (0 before cuts / 1 after cuts / 2 testing) 
        TH1D*                   fEventCounter;  //! event rejection tracker
        TH1D*                   fSampleCounter; //! distribution of events in sampling bins
        // QA events
        TH1D*                   fQAEventsPVz[fQANumSteps];                 //! PV z distance distribution
        TH1D*                   fQAEventsNumContrPV[fQANumSteps];                 //! number of contributors to PV 
        TH1D*                   fQAEventsNumSPDContrPV[fQANumSteps];                 //! number of SPD contributors to PV 
        TH1D*                   fQAEventsDistPVSPD[fQANumSteps];                 //! z-distance between AOD & SPD PV
        // QA tracks
        TH1D*                   fQATracksMult[fQANumSteps];       //! number of AOD tracks distribution
        TH1D*                   fQATracksPt[fQANumSteps];         //! pT dist of all tracks in all events
        TH1D*                   fQATracksEta[fQANumSteps];        //! eta dist of all tracks in all events
        TH1D*                   fQATracksPhi[fQANumSteps];        //! phi dist of all tracks in all events
        TH1D*                   fQATracksFilterMap[fQANumSteps];//! filter bit of all tracks
        TH1D*                   fQATracksNumTPCcls[fQANumSteps];//! dist of track number of TPC clusters
        TH1D*                   fQATracksDCAxy[fQANumSteps]; //! dist of tracks dca in transverse plane
        TH1D*                   fQATracksDCAz[fQANumSteps]; //! dist of tracks dca in z
        TH1D*                   fQATracksTPCstatus[fQANumSteps]; //! based on AliPIDResponse::CheckPIDStatus();
        TH1D*                   fQATracksTOFstatus[fQANumSteps]; //! based on AliPIDResponse::CheckPIDStatus();
        TH2D*                   fQATracksTPCdEdx[fQANumSteps]; //! TPC PID information
        TH2D*                   fQATracksTOF[fQANumSteps]; //! TOF PID information
        TH2D*                   fQATracksTOFbeta[fQANumSteps]; //! TOF PID information
        // QA PID tracks
        TH2D*                   fQAPionPID[fQANumSteps]; //! number of sigmas TPC and TOF PID
        TH2D*                   fQAKaonPID[fQANumSteps]; //! number of sigmas TPC and TOF PID
        TH2D*                   fQAProtonPID[fQANumSteps]; //! number of sigmas TPC and TOF PID
        // QA V0s
        TH1D*					fQAV0sCounter;		//! V0s counter
        TH1D*					fQAV0sCounterK0s;		//! K0s counter
        TH1D*					fQAV0sCounterLambda;		//! Lambda counter
        TH1D*					fQAV0sRecoMethod[fQANumSteps];	//! offline/online V0 reconstruction method
        TH1D*					fQAV0sTPCRefit[fQANumSteps];	//! TPC refit true/false
        TH1D*					fQAV0sKinks[fQANumSteps];	//! V0 kinks true/false
        TH1D*					fQAV0sDCAtoPV[fQANumSteps];	//! V0 DCA to PV
        TH1D*					fQAV0sDCADaughters[fQANumSteps];	//! DCA between V0 daughters
        TH1D*					fQAV0sDecayRadius[fQANumSteps];	//! Distance between PV and Secondary vertex in transverse plane
        TH1D*                   fQAV0sDaughterPt[fQANumSteps];    //! pT dist of V0 daughters
        TH1D*					fQAV0sDaughterPhi[fQANumSteps];	//! pT dist of V0 daughters
        TH1D*                   fQAV0sDaughterEta[fQANumSteps];   //! pseudorapidity dist of V0 daughters
        TH2D*					fQAV0sDaughterTPCdEdxPt[fQANumSteps];	//! TPC dEdx vs pT of V0 daughters
        TH1D*                   fQAV0sMotherPt[fQANumSteps];  //! pT dist of V0s
        TH1D*					fQAV0sMotherPhi[fQANumSteps];	//! azimuthal dist of V0s
        TH1D*                   fQAV0sMotherEta[fQANumSteps]; //! pseudorapidity dist of V0s
        TH1D*                   fQAV0sMotherRapK0s[fQANumSteps];  //! rapidity dist of V0s (K0s mass hypothesis)
        TH1D*                   fQAV0sMotherRapLambda[fQANumSteps];   //! rapidity dist of V0s (Lambda mass hypothesis)
        TH1D*                   fQAV0sInvMassK0s[fQANumSteps];    //! inv. mass dist of V0s (K0s mass hypothesis)
        TH1D*					fQAV0sInvMassLambda[fQANumSteps];	//! inv. mass dist of V0s ((A)Lambda mass hypothesis)
        TH1D*					fQAV0sCPAK0s[fQANumSteps];	//! cosine of pointing angle of K0s candidates
        TH1D*					fQAV0sCPALambda[fQANumSteps];	//! cosine of pointing angle of Lambda candidates
        TH1D*					fQAV0sNumTauK0s[fQANumSteps];	//! number of c*tau of K0s candidates
        TH1D*					fQAV0sNumTauLambda[fQANumSteps];	//! number of c*tau of Lambda candidates
        TH2D*					fQAV0sArmenterosK0s[fQANumSteps];	//! Armenteros-Podolanski plot for K0s candidates
        TH2D*					fQAV0sArmenterosLambda[fQANumSteps];	//! Armenteros-Podolanski plot for K0s candidates
        TH2D*                   fQAV0sProtonNumSigmaPtLambda[fQANumSteps];  //! number of TPC sigmas and pT of proton (Lambda candidates)
        
        AliAnalysisTaskFlowPID(const AliAnalysisTaskFlowPID&); // not implemented
        AliAnalysisTaskFlowPID& operator=(const AliAnalysisTaskFlowPID&); // not implemented

        ClassDef(AliAnalysisTaskFlowPID, 4);
};

#endif
