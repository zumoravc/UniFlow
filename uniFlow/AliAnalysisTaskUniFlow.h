/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUniFlow_H
#define AliAnalysisTaskUniFlow_H

#include "AliAnalysisTaskSE.h"

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
        void                    SetColisionSystem(ColSystem colSystem) {fColSystem = colSystem; }
        void					          SetAnalysisType(AnalType type) { fAnalType = type; }
        void                    SetPeriod(DataPeriod period) { fPeriod = period; }

        void                    SetOldPeriod10h(Bool_t period) { fLHC10h = period; }

        void                    SetTrigger(Short_t trigger) { fTrigger = trigger; }
        void                    SetRejectPileUpSPD(Bool_t pileSPD) { fRejectPileFromSPD = pileSPD; }
        void                    SetUsePlpMV(Bool_t plpMV) { fUsePlpMV = plpMV; }
        void                    SetRejectOutOfBunchPU(Bool_t oobpu) { fRejectOutOfBunchPU = oobpu; }
        void					          SetUseIsPileUpFromSPD(Bool_t IsPileUp) { fUseIsPileUpFromSPD = IsPileUp; }
        void					SetCentFlag(Short_t flag) { fCentFlag = flag; }
        void					SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
        void                    SetSampling(Bool_t sample) { fSampling = sample; }
        void                    SetDoFlow(Bool_t doFlow) { fDoFlow = doFlow; }
        void                    SetDiffFlow(Bool_t diff) { fDiffFlow = diff; }
        void                    SetPID(Bool_t pid) { fPID = pid; }
        void                    SetUseBayesPID(Bool_t pidBay) { fUseBayesPID = pidBay; }
        void                    SetTracksScan(Bool_t qa) { fTracksScan = qa; }
        void					SetDoV0s(Bool_t pidV0s) { fDoV0s = pidV0s; }
        void                    SetDoFlowGenFramKatarina(Bool_t genFlow) { fDoGenFramKat = genFlow; } // do Gen Frame with Katarina's code
        void                    SetDoOldFlow(Bool_t oldFlow) { fOldFlow = oldFlow; } // do old flow (before Katarinas GFK)
        void                    SetUseOldCent(Bool_t oldCent) { fUseOldCent = oldCent; } // use old cent framework
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



        //const static Int_t      fNumPtBins = 29;            // number of pT bins used for pT-differential flow // mine
        const static Int_t 		  fNumPtBins = 33;			// number of pT bins used for pT-differential flow // you, katarina
        static Double_t			    fPtBinEdges[fNumPtBins+1];				// pointer for array of pT bin edges
        const static Int_t      fNumMinvFlowBinsK0s = 12;  // number of inv. mass bin for differential flow plots (K0s)
        static Double_t         fMinvFlowBinEdgesK0s[fNumMinvFlowBinsK0s+1]; // pointer to array of Minv bin edges (K0s)
        const static Int_t      fNumMinvFlowBinsLambda = 11;  // number of inv. mass bin for differential flow plots ((A)Lambda)
        static Double_t         fMinvFlowBinEdgesLambda[fNumMinvFlowBinsLambda+1]; // pointer to array of Minv bin edges ((A)Lambda)
        const static Int_t 		  fNumCentBins = 9;			// number of centrality bins used for pT-differential flow (so far independently of reference flow)
        static Double_t			    fCentBinEdges[fNumCentBins+1];				// pointer for array of pT bin edges
        const static Int_t      fNumHarmonics = 1; // number of harmonics
        static Int_t            fHarmonics[fNumHarmonics]; // values of used harmonics
        const static Int_t      fNumEtaGap = 3; // number of harmonics
        static Double_t         fEtaGap[fNumEtaGap]; // values of used harmonics
        const static Int_t      fNumScanFB = 14; // number of scanned FB
        static Short_t          fTracksScanFB[fNumScanFB]; // values of scanned FBs
        const static Int_t      fMaxNumHarmonics = 8; // maximal number of harmonics for Q,p,q vector arrays
        const static Int_t      fMaxNumWeights = 8; // maximal number of weights for Q,p,q vector arrays

    private:
        Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
        Bool_t                  EventSelection(const AliAODEvent* event); // main method for event selection (specific event selection is applied within)
        Bool_t                  ProcessEvent(const AliAODEvent* event); // main (envelope) method for processing events passing selection
        void                    ListParameters();

        // properties
        Bool_t                  fInit; // initialization check


        //cuts & selection: analysis
        AnalType                fAnalType; // analysis type: AOD / ESD
        ColSystem               fColSystem; // collisional system
        DataPeriod              fPeriod; // period of analysed data sample (e.g. LHC16k, ...)

        Bool_t                  fLHC10h;        // flag to LHC10h data?
        Short_t                 fTrigger;   // switch for pp trigger selection / 0: INT7 / 1: kHighMultV0 / 2: kHighMultSPD
        Bool_t                  fRejectPileFromSPD;   // switch for rejection based on is PileFromSPD
        Bool_t                  fUseIsPileUpFromSPD;   // use IsPileupFromSPD method in (old) event selection
        Bool_t                  fUsePlpMV;   // use plpMV method in (old) event selection
        Bool_t                  fRejectOutOfBunchPU; // use AliUtils out-of-bunch pile-up rejection (provided by Christian)
        Short_t                 fCentFlag;          // centrality flag
        Bool_t                  fSampling;      // Do random sampling ? (estimation of vn stat. uncertanity)
        Bool_t                  fDoFlow;        // Do flow analysis (if kFALSE: only selection, filtering and QA)
        Bool_t                  fDoV0s;                       // Do V0s analysis?
        Bool_t                  fDiffFlow;          // Do differential flow ? (or reference only)
        Bool_t                  fTracksScan;      // do scan Pt, Phi, eta of various track FBs
        Bool_t                  fPID;                       // Do PID ?
        Bool_t                  fUseBayesPID;                       // Do Bayes PID ?
        Bool_t                  fDoGenFramKat; // switch gen frame. by Katarina
        Bool_t                  fOldFlow; // switch to old flow analysis
        Bool_t                  fUseOldCent; // switch for old (Run1) centrality selection
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

        AliAnalysisTaskUniFlow(const AliAnalysisTaskUniFlow&); // not implemented
        AliAnalysisTaskUniFlow& operator=(const AliAnalysisTaskUniFlow&); // not implemented

        ClassDef(AliAnalysisTaskUniFlow, 0);
};

#endif
