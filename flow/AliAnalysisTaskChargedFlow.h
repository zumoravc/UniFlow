#ifndef ALIANALYSISTASKCHARGEDFLOW_H
#define ALIANALYSISTASKCHARGEDFLOW_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TComplex.h>

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskChargedFlow : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskChargedFlow();
  AliAnalysisTaskChargedFlow(const char *name);

  virtual ~AliAnalysisTaskChargedFlow();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t* option);
  virtual void   Terminate(Option_t* ); 

  Double_t GetVtxCut() {return fVtxCut;}
  Double_t GetEtaCut() {return fEtaCut;}     
  Double_t GetMinPt()  {return fMinPt;}
  //Int_t    GetNPtBins(){return fNPtBins;}
 
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetFilterbit(UInt_t filterbit){fFilterbit = filterbit;}	
  virtual void  SetNoClus(Int_t noclus){fNoClus = noclus;}
  virtual void  SetEff(Int_t effsys){fEffsys = effsys;}
  virtual void  SetCharge(Int_t charge){fULS = charge;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
  virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
  virtual void  SetNHarmonic(Short_t nHarm){fNHarm = nHarm;}
  virtual void  SetIsSample(Bool_t IsSample){fSample = IsSample;}
  virtual void  SetFlagLHC10h(Bool_t IsLHC10h){fLHC10h = IsLHC10h;}	
  virtual void  SetFlagPileUp(Bool_t IsPileUP){fPileUp = IsPileUP;}
  virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
  virtual void  SetNPtBins(Int_t nPtB){fNPtBins = nPtB;}
  virtual void  SetPtBins(Double_t ptBins[15]) {for (Int_t i = 0; i < 15; i++) fPtBins[i] = ptBins[i];}
  virtual void  SetCentFlag(Short_t nCent){fCent = nCent;}
	virtual void	SetDifferentialFlowFlag(Bool_t differentialFlow){fDifferential = differentialFlow;}
		

 private:
  AliAnalysisTaskChargedFlow(const AliAnalysisTaskChargedFlow&);
  AliAnalysisTaskChargedFlow& operator=(const AliAnalysisTaskChargedFlow&);

  virtual void AnalyzeAOD(AliAODEvent* aod, Short_t centrV0);  
  virtual Float_t GetVertex(AliVEvent* ev) const;
  virtual Int_t GetTPCMult(AliVEvent* ev) const;
  virtual Int_t GetGlobalMult(AliVEvent* ev) const;
  Short_t  GetCentrCode(AliVEvent* ev);
  Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  Bool_t plpMV(const AliVEvent *event);
  void ProcessMCTruth(AliAODEvent* aod, Short_t cV0);

	double GetWeight(double phi, double eta);

	TComplex Qvector[20][20];
	TComplex QvectorM[20][20];
	TComplex QvectorP[20][20];
	TComplex Qvector0M[20][20];
	TComplex Qvector0P[20][20];
	TComplex Qvector4M[20][20];
	TComplex Qvector4P[20][20];
	TComplex Qvector8M[20][20];
	TComplex Qvector8P[20][20];
	TComplex pvector[20][20];
	TComplex pvectorM[20][20];
	TComplex pvectorP[20][20];
	TComplex qvector[20][20];
	TComplex pvector0M[20][20];
	TComplex pvector0P[20][20];
	TComplex pvector4M[20][20];
	TComplex pvector4P[20][20];
	TComplex pvector8M[20][20];
	TComplex pvector8P[20][20];
	
	TComplex Q(int n, int p);
	TComplex QGapM(int n, int p);
	TComplex QGapP(int n, int p);
	TComplex QGap0M(int n, int p);
	TComplex QGap0P(int n, int p);
	TComplex QGap4M(int n, int p);
	TComplex QGap4P(int n, int p);
	TComplex QGap8M(int n, int p);
	TComplex QGap8P(int n, int p);
	TComplex p(int n, int p);
	TComplex pGapM(int n, int p);
	TComplex pGapP(int n, int p);
	TComplex q(int n, int p);
	TComplex pGap0M(int n, int p);
	TComplex pGap0P(int n, int p);
	TComplex pGap4M(int n, int p);
	TComplex pGap4P(int n, int p);
	TComplex pGap8M(int n, int p);
	TComplex pGap8P(int n, int p);
	void ResetQ(const int nMaxHarm, const int nMaxPow);

	TComplex Two(int n1, int n2);
	TComplex TwoTest(int n1, int n2);
	TComplex TwoGap(int n1, int n2);
	TComplex TwoGap0(int n1, int n2);
	TComplex TwoGap4(int n1, int n2);
	TComplex TwoGap8(int n1, int n2);
	TComplex TwoDiff(int n1, int n2);
	TComplex TwoDiffGapQM(int n1, int n2);
	TComplex TwoDiffGapQP(int n1, int n2);
	TComplex TwoDiffGap0QP(int n1, int n2);
	TComplex TwoDiffGap0QM(int n1, int n2);
	TComplex TwoDiffGap4QP(int n1, int n2);
	TComplex TwoDiffGap4QM(int n1, int n2);
	TComplex TwoDiffGap8QP(int n1, int n2);
	TComplex TwoDiffGap8QM(int n1, int n2);
	TComplex TwoDiffPtGap0(int n1, int n2);
	TComplex TwoDiffPtGap4(int n1, int n2);
	TComplex TwoDiffPtGap8(int n1, int n2);
	TComplex TwoDiffPtGap10(int n1, int n2);
	/*
	TComplex Three(int n1, int n2, int n3);
	TComplex ThreeDiff(int n1, int n2, int n3);
	TComplex Four(int n1, int n2, int n3, int n4);
	TComplex FourDiff(int n1, int n2, int n3, int n4);
	TComplex FourPtDiff(int n1, int n2, int n3, int n4);
	TComplex Five(int n1, int n2, int n3, int n4, int n5);
	TComplex Six(int n1, int n2, int n3, int n4, int n5, int n6);
	TComplex Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
	TComplex Eight_1(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
	TComplex Eight_2(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
	TComplex Eight_3(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
	TComplex Eight_4(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8);
  */
	AliAODEvent* fAOD;                //! AOD object
  TString      fAnalysisType;       //  "ESD" or "AOD
  Bool_t       fAnalysisMC;         // Real(kFALSE) or MC(kTRUE) flag

  // Cuts and options
  Double_t     fVtxCut;             // Vtx cut on z position in cm
  UInt_t       fFilterbit;          // filter bit
  Double_t     fEtaCut;             // Eta cut used to select particles
  Int_t        fNoClus;				// No of TPC clusters
  Int_t        fEffsys;				// efficiency correction
  Int_t        fULS;                            // charge, 1:all, 2:pp,  3: mm
  Double_t     fMinPt;              // Min pt - for histogram limits
  Double_t     fMaxPt;              // Max pt - for histogram limits
  Short_t      fNHarm;              // harmonic number
  Short_t      fSample;              // number of sample
  Bool_t       fLHC10h;             // flag to LHC10h data
  Bool_t       fPileUp;             // flag for pileup
  Int_t        fNPtBins;            // number of pt bins
  Short_t      fCent;               // centrality flag
  Double_t     fPtBins[15];         // pt bins edges
	Bool_t			 fDifferential; 				//flag for differential flow
 
  // Output objects
  TList*        fListOfObjects;     //! Output list of objects
  TH1I*         fVtx;               //! Event vertex info
  TH1F*         fVtxBeforeCuts;     //! Vertex z dist before cuts
  TH1F*         fVtxAfterCuts;      //! Vertex z dist after cuts
  TH1F*         fCentralityDis;     //! centrality dist
  TH1F*         fPhiDis;     //! phi dist
  TH1F*         fPhiDisGapM;     //! phi dist
  TH1F*         fPhiDisGapP;     //! phi dist
  TH1F*         fPhiDisGap0M;     //! phi dist
  TH1F*         fPhiDisGap0P;     //! phi dist
  TH1F*         fPhiDisGap4M;     //! phi dist
  TH1F*         fPhiDisGap4P;     //! phi dist
  TH1F*         fPhiDisGap8M;     //! phi dist
  TH1F*         fPhiDisGap8P;     //! phi dist
  TH1F*         fPhiDispT2;     //! phi dist
  TH1F*         fPhiDispT3;     //! phi dist
  TH1F*         fPhiDispT5;     //! phi dist
	TH1F*					fPhiDisVtx10;		//!
	TH1F*					fPhiDisVtx9;		//!
	TH1F*					fPhiDisVtx8;		//!
	TH1F*					fPhiDisVtx7;		//!
  TH1F*         fPhiDisCent[10];     //! phi dist
	TH1F*					fPhiPt[15];		//! phi in different pT bins
  TH1F*         fEtaDis;     //! eta dist
  TH2F*         fMultCorBeforeCuts; //! correlation between TPC and global track multiplicity before cuts from flow package(centrality)
  TH2F*         fMultCorAfterCuts;  //! correlation between TPC and global track multiplicity after cuts from flow package(centrality)
  TH2F*         fCentSPDvsV0M;      //! V0M vs SPD
  TH2F*         fMultvsCentr;       //! multiplicity vs centrality

	TH1F*					fWeight;	//!
  
    TProfile*    fcn2Re[6][12];  //!
    TProfile*    fcn2GapRe[6][12];  //!
    TProfile*    fcn2Gap0Re[6][12];  //!
    TProfile*    fcn2Gap4Re[6][12];  //!
    TProfile*    fcn2Gap8Re[6][12];  //!
    TProfile*    fcn4Re[6][12];  //!
    TProfile*    fcn6Re[6][12];  //!
    TProfile*    fcn8Re[6][12];  //!

	//..Standard Candles
	TProfile*			fsc4242Re[12];	//!
	TProfile*			fsc3232Re[12];	//!
	TProfile*			fsc4343Re[12];	//!
//	TProfile*			fsc5252Re[12];	//!	
//	TProfile*			fsc5353Re[12];	//!

	//..pT differential flow

	TProfile*			fdn2Re[6][10][12];	//! 
	TProfile*			fdn2GapQMRe[6][10][12];	//! 
	TProfile*			fdn2GapQPRe[6][10][12];	//! 
	TProfile*			fdn2Gap0QMRe[6][10][12];	//! 
	TProfile*			fdn2Gap0QPRe[6][10][12];	//! 
	TProfile*			fdn2Gap4QMRe[6][10][12];	//! 
	TProfile*			fdn2Gap4QPRe[6][10][12];	//! 
	TProfile*			fdn2Gap8QMRe[6][10][12];	//! 
	TProfile*			fdn2Gap8QPRe[6][10][12];	//! 
	TProfile*			fdn4Re[6][10][12];	//! 

	TProfile*			fcn2PtGap0Re[6][10][12]; //!
	TProfile*			fcn2PtGap4Re[6][10][12]; //!
	TProfile*			fcn2PtGap8Re[6][10][12]; //!
	TProfile*			fcn2PtGap10Re[6][10][12]; //!

	TProfile*			fcn4PtRe[6][10][12]; //! 

	ClassDef(AliAnalysisTaskChargedFlow, 1);    //Analysis task
};

#endif
