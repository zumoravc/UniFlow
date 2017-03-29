/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUniFlow_H
#define AliAnalysisTaskUniFlow_H

#include "AliAnalysisTaskSE.h"

enum RunMode {kTest, kFull}; // task running mode (NOT GRID MODE)
enum ColSystem {kPP, kPPb, kPbPb}; // tag for collisional system
enum AnalType {kAOD, kESD}; // tag for analysis type
enum DataPeriod {kNon, k16k, k16l, k16q, k16r, k16s, k16t}; // tag for data period

class AliAnalysisTaskUniFlow : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskUniFlow(); // constructor
                                AliAnalysisTaskUniFlow(const char *name); // named (primary) constructor
        virtual                 ~AliAnalysisTaskUniFlow(); // destructor

        virtual void            UserCreateOutputObjects(); //
        virtual void            UserExec(Option_t* option); // main methond - called for each event
        virtual void            Terminate(Option_t* option); // called after all events are processed
        // event setters
        void                    SetRunMode(RunMode mode = kFull) { fRunMode = mode; }
        void                    SetNumEventsAnalyse(Short_t num) { fNumEventsAnalyse = num; }
        void                    SetColisionSystem(ColSystem colSystem = kPP) {fColSystem = colSystem; }
        void					          SetAnalysisType(AnalType type = kAOD) { fAnalType = type; }
        void                    SetPeriod(DataPeriod period = kNon) { fPeriod = period; }
        void                    SetTrigger(Short_t trigger = 0) { fTrigger = trigger; }
        void                    SetSampling(Bool_t sample = kTRUE) { fSampling = sample; }
        void                    SetNumberOfSamples(Short_t numSamples = 10) { fNumSamples = numSamples; }
        void                    SetProcessCharged(Bool_t filter = kTRUE) { fProcessCharged = filter; }
        void                    SetProcessPID(Bool_t filter = kTRUE) { fProcessPID = filter; }
        void                    SetProcessV0s(Bool_t filter = kTRUE) { fProcessV0s = filter; }

        void					SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        // track setters
        void                    SetChargedEtaMax(Double_t eta) { fCutChargedEtaMax = eta; }
        void                    SetChargedPtMax(Double_t pt) { fCutChargedPtMax = pt; }
        void                    SetChargedPtMin(Double_t pt) { fCutChargedPtMin = pt; }
        void                    SetChargedDCAzMax(Double_t dcaz) {  fCutChargedDCAzMax = dcaz; }
        void                    SetChargedDCAxyMax(Double_t dcaxy) {  fCutChargedDCAxyMax = dcaxy; }
        void                    SetChargedNumTPCclsMin(UShort_t tpcCls) { fCutChargedNumTPCclsMin = tpcCls; }
        void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; }
        // PID (pi,K,p) setters
        void                    SetPionNumSigmasMax(Double_t numSigmas) { fCutPionNumSigmaMax = numSigmas; }
        void                    SetKaonNumSigmasMax(Double_t numSigmas) { fCutKaonNumSigmaMax = numSigmas; }
        void                    SetProtonNumSigmasMax(Double_t numSigmas) { fCutProtonNumSigmaMax = numSigmas; }
        void                    SetPIDBayesProbPionMin(Double_t probPi) { fCutBayesPIDPionMin = probPi; }
        void                    SetPIDBayesProbKaonMin(Double_t probK) { fCutBayesPIDKaonMin = probK; }
        void                    SetPIDBayesProbProtonMin(Double_t probP) { fCutBayesPIDProtonMin = probP; }
        // V0s setters
        void					          SetV0sOnFly(Bool_t onFly) { fCutV0onFly = onFly; }
        void					          SetV0sTPCRefit(Bool_t refit) { fCutV0refitTPC = refit; }
        void					          SetV0sRejectKinks(Bool_t reject) { fCutV0rejectKinks = reject; }
        void					          SetV0sDCAPVMin(Double_t dca) { fCutV0MinDCAtoPV = dca; }
        void					          SetV0sDCAPVMax(Double_t dca) { fCutV0MaxDCAtoPV = dca; }
        void					          SetV0sDCADaughtersMin(Double_t dca) { fCutV0MinDCADaughters = dca; }
        void					          SetV0sDCADaughtersMax(Double_t dca) { fCutV0MaxDCADaughters = dca; }
        void					          SetV0sDecayRadiusMin(Double_t radius) { fCutV0MinDecayRadius = radius; }
        void					          SetV0sDecayRadiusMax(Double_t radius) { fCutV0MaxDecayRadius = radius; }
        void					          SetV0sDaughterPtMin(Double_t pt) { fCutV0DaughterPtMin = pt; }
        void					          SetV0sDaughterPtMax(Double_t pt) { fCutV0DaughterPtMax = pt; }
        void					          SetV0sDaughterEtaMax(Double_t eta) { fCutV0DaughterEtaMax = eta; }
        void					          SetV0sMotherEtaMax(Double_t eta) { fCutV0MotherEtaMax = eta; }
        void                    SetV0sMotherRapMax(Double_t rap) { fCutV0MotherRapMax = rap; }
        void                    SetV0sMotherPtMin(Double_t pt) { fCutV0MotherPtMin = pt; }
        void					          SetV0sK0sInvMassMin(Double_t mass) { fCutV0K0sInvMassMin = mass; }
        void					          SetV0sK0sInvMassMax(Double_t mass) { fCutV0K0sInvMassMax = mass; }
        void					          SetV0sLambdaInvMassMin(Double_t mass) { fCutV0LambdaInvMassMin = mass; }
        void					          SetV0sLambdaInvMassMax(Double_t mass) { fCutV0LambdaInvMassMax = mass; }
        void					          SetV0sMotherPtMax(Double_t pt) { fCutV0MotherPtMax = pt; }
        void					          SetV0sK0sCPAMin(Double_t cpa) { fCutV0MinCPAK0s = cpa; }
        void					          SetV0sLambdaCPAMin(Double_t cpa) { fCutV0MinCPALambda = cpa; }
        void					          SetV0sK0sNumTauMax(Double_t nTau) { fCutV0NumTauK0sMax = nTau; }
        void					          SetV0sLambdaNumTauMax(Double_t nTau) { fCutV0NumTauLambdaMax = nTau; }
        void					          SetV0sK0sArmenterosAlphaMin(Double_t alpha) { fCutV0K0sArmenterosAlphaMin = alpha; }
        void					          SetV0sLambdaArmenterosAlphaMax(Double_t alpha) { fCutV0LambdaArmenterosAlphaMax = alpha; }
        void                    SetV0sProtonNumSigmaMax(Double_t nSigma) { fCutV0ProtonNumSigmaMax = nSigma; }
        void					          SetV0sProtonPIDPtMin(Double_t pt) { fCutV0ProtonPIDPtMin = pt; }
        void					          SetV0sProtonPIDPtMax(Double_t pt) { fCutV0ProtonPIDPtMax = pt; }

    private:
        Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
        Bool_t                  EventSelection(); // main method for event selection (specific event selection is applied within)
        Bool_t                  IsEventSelected_2016(); // event selection for LHC2016 pp & pPb data
        void                    FillEventsQA(const Short_t iQAindex); // filling QA plots related to event selection

        Bool_t                  Filtering(); // main (envelope) method for filtering all POIs in event
        Bool_t                  FilterCharged(); // charged tracks filtering
        Bool_t                  IsChargedSelected(const AliAODTrack* track = 0x0); // charged track selection
        Bool_t                  FilterPID(); // pi,K,p filtering
        Bool_t                  FilterV0s(); // K0s, Lambda, ALambda filtering
        Bool_t                  IsV0Selected(const AliAODv0* v0 = 0x0); // general (common) V0 selection
        Bool_t                  IsV0aK0s(const AliAODv0* v0 = 0x0); // V0 selection: K0s specific
        Bool_t                  IsV0aLambda(const AliAODv0* v0 = 0x0); // V0 selection: (A)Lambda specific
        void                    FillQACharged(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for charged track selection
        void                    FillQAV0s(const Short_t iQAindex, const AliAODv0* v0 = 0x0, const Bool_t bIsK0s = kFALSE, const Bool_t bIsLambda = kFALSE); // filling QA plots for V0s candidates

        Bool_t                  ProcessEvent(); // main (envelope) method for processing events passing selection
        void                    ListParameters();

        Short_t                 GetSamplingIndex(); // returns sampling index based on sampling selection (number of samples)
        Short_t                 GetCentralityIndex(); // returns centrality index based centrality estimator or number of selected tracks

        // properties
        AliAODEvent*            fEventAOD; //! AOD event countainer
        AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
        Bool_t                  fInit; // initialization check
        Short_t                 fIndexSampling; // sampling index (randomly generated)
        Short_t                 fIndexCentrality; // centrality bin index (based on centrality est. or number of selected tracks)
        Short_t                 fEventCounter; // event counter (used for local test runmode purpose)
        Short_t                 fNumEventsAnalyse; // [50] number of events to be analysed / after passing selection (only in test mode)

        // selected POIs containers
        TClonesArray*           fArrCharged; //! container for filtered (all) charged tracks
        TClonesArray*           fArrChargedRPF; //! container for filtered RFPs tracks
        TClonesArray*           fArrChargedPOI; //! container for filtered POIs tracks
        TClonesArray*           fArrPion; //! container for filtered pions
        TClonesArray*           fArrKaon; //! container for filtered kaons
        TClonesArray*           fArrProton; //! container for filtered protons
        TClonesArray*           fArrK0s; //! container for filtered K0s candidates
        TClonesArray*           fArrLambda; //! container for filtered Lambda candidates
        TClonesArray*           fArrALambda; //! container for filtered ALambda candidates

        //cuts & selection: analysis
        RunMode                 fRunMode; // running mode (not grid related)
        AnalType                fAnalType; // analysis type: AOD / ESD
        ColSystem               fColSystem; // collisional system
        DataPeriod              fPeriod; // period of analysed data sample (e.g. LHC16k, ...)
        Bool_t                  fSampling;      // Do random sampling ? (estimation of vn stat. uncertanity)
        Short_t                 fTrigger; // physics selection trigger
        Short_t                 fNumSamples; // overall number of samples (from random sampling) used
        Bool_t                  fProcessCharged; // flag for processing charged tracks (both RPF and POIs)
        Bool_t                  fProcessPID; // flag for processing PID tracks (pi,K,p)
        Bool_t                  fProcessV0s; // flag for processing V0 candidates (K0s, Lambda/ALambda)

        //cuts & selection: events & tracks
        Float_t                 fPVtxCutZ;          // (cm) PV z cut

        UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
        UShort_t                fCutChargedNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
        Float_t                 fCutChargedEtaMax; // (-) Maximum pseudorapidity range
        Float_t                 fCutChargedPtMax; // (GeV/c) Maximal track pT
        Float_t                 fCutChargedPtMin; // (GeV/c) Minimal track pT
        Float_t                 fCutChargedDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
        Float_t                 fCutChargedDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
        // cuts & selection: PID selection
        Double_t                fCutPionNumSigmaMax;
        Double_t                fCutKaonNumSigmaMax;
        Double_t                fCutProtonNumSigmaMax;
        Double_t                fCutBayesPIDPionMin; // minimal value of Bayes PID probability for pion
        Double_t                fCutBayesPIDKaonMin; // minimal value of Bayes PID probability for Kaon
        Double_t                fCutBayesPIDProtonMin; // minimal value of Bayes PID probability for proton
        //cuts & selection: V0 reconstruction
		    Bool_t 					        fCutV0onFly;		// V0 reconstruction method: is On-the-fly? (or offline)
    		Bool_t					        fCutV0refitTPC; // Check TPC refit of V0 daughters ?
    		Bool_t					        fCutV0rejectKinks; // Reject Kink V0 daughter tracks ?
    		Double_t                fCutV0MinDCAtoPV;   // (cm) min DCA of V0 daughter to PV
        Double_t				        fCutV0MaxDCAtoPV;	// (cm) max DCA of V0 daughter to PV
  		  Double_t				        fCutV0MinDCADaughters;	// (cm) min DCA of V0 daughters among themselves
  		  Double_t				        fCutV0MaxDCADaughters;	// (cm) max DCA of V0 daughters among themselves
        Double_t                fCutV0MinDecayRadius; // (cm) min distance of secondary vertex from z-axis in transverse plane
  		  Double_t				        fCutV0MaxDecayRadius; // (cm) max distance of secondary vertex from z-axis in transverse plane
        Double_t                fCutV0DaughterPtMin; // (GeV/c) min pT of V0 daughters
        Double_t                fCutV0DaughterPtMax; // (GeV/c) max pT of V0 daughters
        Double_t                fCutV0DaughterEtaMax; // (-) max value of Eta of V0 daughters
        Double_t                fCutV0MotherEtaMax; // (-) max eta value of V0 mother
        Double_t                fCutV0MotherRapMax; // (-) max rapidity value of V0 mother
        Double_t                fCutV0MotherPtMin; // (GeV/c) min transverse momentum value of V0 mother
        Double_t                fCutV0MotherPtMax; // (GeV/c) max transverse momentum value of V0 mother
        Double_t                fCutV0MinCPAK0s;    // (-) min cosine of pointing angle of K0s candidate to PV
        Double_t                fCutV0MinCPALambda; // (-) min cosine of pointing angle of K0s candidate to PV
        Double_t                fCutV0NumTauK0sMax; // (c*tau) max number of c*tau (K0s)
        Double_t                fCutV0NumTauLambdaMax; // (c*tau) max number of c*tau ((A)Lambda)
        Double_t				        fCutV0K0sArmenterosAlphaMin; // (alpha) min Armenteros alpha for K0s
        Double_t                fCutV0K0sInvMassMin; // [0.4] (GeV/c2) min inv. mass window for selected K0s candidates
        Double_t                fCutV0K0sInvMassMax; // [0.6] (GeV/c2) max inv. mass window for selected K0s candidates
        Double_t                fCutV0LambdaInvMassMin; // [1.08] (GeV/c2) min inv. mass window for selected (Anti)Lambda candidates
        Double_t                fCutV0LambdaInvMassMax; // [1.16] (GeV/c2) max inv. mass window for selected (Anti)Lambda candidates
        Double_t                fCutV0LambdaArmenterosAlphaMax; // (alpha) max Armenteros alpha for (Anti)Lambda
        Double_t                fCutV0ProtonNumSigmaMax;    // (sigmaTPC) max number of TPC sigma for proton PID (Lambda candidates)
        Double_t				        fCutV0ProtonPIDPtMin;	// (GeV/c) min pT of proton for PID (Lambda candidates) - only protons within pT range will be checked for num sigma TPC
        Double_t				        fCutV0ProtonPIDPtMax;	// (GeV/c) max pT of proton for PID (Lambda candidates) - only protons within pT range will be checked for num sigma TPC

        // output lists
        TList*      fOutListEvents; //! events list
        TList*      fOutListCharged; //! charged tracks list
        TList*      fOutListV0s; //! V0s candidates list

        // histograms
        TH1D*           fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
        TH1D*           fhEventCounter; //! counter following event selection
        TH1D*           fhChargedCounter; //! counter following charged track selection
        TH1D*           fhV0sCounter; //! counter following V0s selection
        TH1D*           fhV0sCounterK0s; //! counter following K0s selection
        TH1D*           fhV0sCounterLambda; //! counter following (Anti-)Lambda selection

        static const Short_t   fiNumIndexQA = 2; // 0: before cuts // 1: after cuts
        // QA: events
        TH1D*           fhQAEventsPVz[fiNumIndexQA]; //!
        TH1D*           fhQAEventsNumContrPV[fiNumIndexQA]; //!
        TH1D*           fhQAEventsNumSPDContrPV[fiNumIndexQA]; //!
        TH1D*           fhQAEventsDistPVSPD[fiNumIndexQA]; //!
        TH1D*           fhQAEventsSPDresol[fiNumIndexQA]; //!
        // QA: charged tracks
        TH1D*           fhQAChargedMult[fiNumIndexQA];       //! number of AOD charged tracks distribution
        TH1D*           fhQAChargedPt[fiNumIndexQA];         //! pT dist of charged tracks
        TH1D*           fhQAChargedEta[fiNumIndexQA];        //! eta dist of charged tracks
        TH1D*           fhQAChargedPhi[fiNumIndexQA];        //! phi dist of charged tracks
        TH1D*           fhQAChargedCharge[fiNumIndexQA];     //! charge dist of charged tracks
        TH1D*           fhQAChargedFilterBit[fiNumIndexQA];  //! filter bit distribution of charged tracks
        TH1D*           fhQAChargedNumTPCcls[fiNumIndexQA];  //! dist of track number of TPC clusters
        TH1D*           fhQAChargedDCAxy[fiNumIndexQA];      //! dist of Charged DCA in transverse plane
        TH1D*           fhQAChargedDCAz[fiNumIndexQA];       //! dist of charged DCA in z coordinate
        TH1D*           fhQAChargedTPCstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
        TH1D*           fhQAChargedTOFstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
        TH2D*           fhQAChargedTPCdEdx[fiNumIndexQA];    //! TPC PID information
        TH2D*           fhQAChargedTOFbeta[fiNumIndexQA];    //! TOF PID information
        // QA: V0s candidates
        TH1D*					fhQAV0sRecoMethod[fiNumIndexQA];	//! offline/online V0 reconstruction method
        TH1D*					fhQAV0sTPCRefit[fiNumIndexQA];	//! TPC refit true/false
        TH1D*					fhQAV0sKinks[fiNumIndexQA];	//! V0 kinks true/false
        TH1D*					fhQAV0sDCAtoPV[fiNumIndexQA];	//! V0 DCA to PV
        TH1D*					fhQAV0sDCADaughters[fiNumIndexQA];	//! DCA between V0 daughters
        TH1D*					fhQAV0sDecayRadius[fiNumIndexQA];	//! Distance between PV and Secondary vertex in transverse plane
        TH1D*         fhQAV0sInvMassK0s[fiNumIndexQA];    //! inv. mass dist of V0s (K0s mass hypothesis)
        TH1D*					fhQAV0sInvMassLambda[fiNumIndexQA];	//! inv. mass dist of V0s ((A)Lambda mass hypothesis)
        TH1D*         fhQAV0sMotherPt[fiNumIndexQA];  //! pT dist of V0s
        TH1D*					fhQAV0sMotherPhi[fiNumIndexQA];	//! azimuthal dist of V0s
        TH1D*         fhQAV0sMotherEta[fiNumIndexQA]; //! pseudorapidity dist of V0s
        TH1D*         fhQAV0sMotherCharge[fiNumIndexQA]; //! charge distribution of mothers
        TH1D*         fhQAV0sMotherRapK0s[fiNumIndexQA];  //! rapidity dist of V0s (K0s mass hypothesis)
        TH1D*         fhQAV0sMotherRapLambda[fiNumIndexQA]; //! rapidity dist of V0s (Lambda mass hypothesis)
        TH1D*         fhQAV0sDaughterPt[fiNumIndexQA];    //! pT dist of V0 daughters
        TH1D*					fhQAV0sDaughterPhi[fiNumIndexQA];	//! pT dist of V0 daughters
        TH1D*         fhQAV0sDaughterEta[fiNumIndexQA];   //! pseudorapidity dist of V0 daughters
        TH1D*         fhQAV0sDaughterCharge[fiNumIndexQA]; //! charge distribution of daughters
        TH2D*					fhQAV0sDaughterTPCdEdxP[fiNumIndexQA];	//! TPC dEdx vs p of V0 daughters
        TH1D*					fhQAV0sCPAK0s[fiNumIndexQA];	//! cosine of pointing angle of K0s candidates
        TH1D*					fhQAV0sCPALambda[fiNumIndexQA];	//! cosine of pointing angle of Lambda candidates
        TH1D*					fhQAV0sNumTauK0s[fiNumIndexQA];	//! number of c*tau of K0s candidates
        TH1D*					fhQAV0sNumTauLambda[fiNumIndexQA];	//! number of c*tau of Lambda candidates
        TH2D*					fhQAV0sArmenterosK0s[fiNumIndexQA];	//! Armenteros-Podolanski plot for K0s candidates
        TH2D*					fhQAV0sArmenterosLambda[fiNumIndexQA];	//! Armenteros-Podolanski plot for K0s candidates
        TH2D*         fhQAV0sProtonNumSigmaPtLambda[fiNumIndexQA];  //! number of TPC sigmas and pT of proton (Lambda candidates)


        AliAnalysisTaskUniFlow(const AliAnalysisTaskUniFlow&); // not implemented
        AliAnalysisTaskUniFlow& operator=(const AliAnalysisTaskUniFlow&); // not implemented

        ClassDef(AliAnalysisTaskUniFlow, 0);
};

#endif
