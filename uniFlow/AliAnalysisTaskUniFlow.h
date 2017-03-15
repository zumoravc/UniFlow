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
        void                    SetColisionSystem(ColSystem colSystem = kPP) {fColSystem = colSystem; }
        void					          SetAnalysisType(AnalType type = kAOD) { fAnalType = type; }
        void                    SetPeriod(DataPeriod period = kNon) { fPeriod = period; }
        void                    SetTrigger(Short_t trigger = 0) { fTrigger = trigger; }
        void                    SetSampling(Bool_t sample = kTRUE) { fSampling = sample; }
        void                    SetNumberOfSamples(Short_t numSamples = 10) { fNumSamples = numSamples; }
        void                    SetFilterCharged(Bool_t filter = kTRUE) { fProcessCharged = filter; }
        void                    SetProcessPID(Bool_t filter = kTRUE) { fProcessPID = filter; }
        void                    SetProcessV0s(Bool_t filter = kTRUE) { fProcessV0s = filter; }

        void					SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        // track setters
        void                    SetTrackEtaMax(Double_t eta) { fTrackEtaMax = eta; }
        void                    SetTrackPtMax(Double_t pt) { fTrackPtMax = pt; }
        void                    SetTrackPtMin(Double_t pt) { fTrackPtMin = pt; }
        void                    SetTrackDCAzMax(Double_t dcaz) {  fTracksDCAzMax = dcaz; }
        void                    SetNumTPCclsMin(UShort_t tpcCls) { fNumTPCclsMin = tpcCls; }
        void                    SetTrackFilterBit(UInt_t filter) { fTrackFilterBit = filter; }
        void                    SetPionNumSigmasMax(Double_t numSigmas) { fCutPionNumSigmaMax = numSigmas; }
        void                    SetKaonNumSigmasMax(Double_t numSigmas) { fCutKaonNumSigmaMax = numSigmas; }
        void                    SetProtonNumSigmasMax(Double_t numSigmas) { fCutProtonNumSigmaMax = numSigmas; }
        void                    SetPIDBayesProbPionMin(Double_t probPi) { fCutBayesPIDPionMin = probPi; }
        void                    SetPIDBayesProbKaonMin(Double_t probK) { fCutBayesPIDKaonMin = probK; }
        void                    SetPIDBayesProbProtonMin(Double_t probP) { fCutBayesPIDProtonMin = probP; }
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

    private:
        Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
        Bool_t                  EventSelection(); // main method for event selection (specific event selection is applied within)
        Bool_t                  IsEventSelected_2016(); // event selection for LHC2016 pp & pPb data
        Bool_t                  Filtering(); // main (envelope) method for filtering all POIs in event
        Bool_t                  FilterPID(); // pi,K,p filtering
        Bool_t                  FilterV0s(); // K0s, Lambda, ALambda filtering
        Bool_t                  ProcessEvent(); // main (envelope) method for processing events passing selection
        void                    ListParameters();

        Short_t                 GetSamplingIndex(); // returns sampling index based on sampling selection (number of samples)
        Short_t                 GetCentralityIndex(); // returns centrality index based centrality estimator or number of selected tracks

        // properties
        AliAODEvent*            fEventAOD; // AOD event countainer
        Bool_t                  fInit; // initialization check
        Short_t                 fIndexSampling; // sampling index (randomly generated)
        Short_t                 fIndexCentrality; // centrality bin index (based on centrality est. or number of selected tracks)
        Short_t                 fEventCounter; // event counter (used for local test runmode purpose)

        // selected POIs containers
        TClonesArray*           fArrChargedRPF; // container for filtered RFPs tracks
        TClonesArray*           fArrChargedPOI; // container for filtered POIs tracks
        TClonesArray*           fArrPion; // container for filtered pions
        TClonesArray*           fArrKaon; // container for filtered kaons
        TClonesArray*           fArrProton; // container for filtered protons
        TClonesArray*           fArrK0s; // container for filtered K0s candidates
        TClonesArray*           fArrLambda; // container for filtered Lambda candidates
        TClonesArray*           fArrALambda; // container for filtered ALambda candidates

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
        Double_t                fPVtxCutZ;          // (cm) PV z cut
        UInt_t                  fTrackFilterBit; // tracks filter bit
        UShort_t                fNumTPCclsMin;  // () Minimal number of TPC clusters used for track reconstruction
        Double_t                fTrackEtaMax;   // () Maximum pseudorapidity range
        Double_t                fTrackPtMax;        // (GeV/c) Maximal track pT
        Double_t                fTrackPtMin;        // (GeV/c) Minimal track pT
        Double_t                fTracksDCAzMax; //  (cm) Maximal DCA cuts for tracks (pile-up rejection suggested for LHC16)
        Double_t                fCutPionNumSigmaMax;
        Double_t                fCutKaonNumSigmaMax;
        Double_t                fCutProtonNumSigmaMax;
        Double_t                fCutBayesPIDPionMin; // minimal value of Bayes PID probability for pion
        Double_t                fCutBayesPIDKaonMin; // minimal value of Bayes PID probability for Kaon
        Double_t                fCutBayesPIDProtonMin; // minimal value of Bayes PID probability for proton
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

        // output lists
        TList*      fOutListEvents; //! events list

        // histograms
        TH1D*           fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
        TH1D*           fhEventCounter; //! counter following event selection

        AliAnalysisTaskUniFlow(const AliAnalysisTaskUniFlow&); // not implemented
        AliAnalysisTaskUniFlow& operator=(const AliAnalysisTaskUniFlow&); // not implemented

        ClassDef(AliAnalysisTaskUniFlow, 0);
};

#endif
