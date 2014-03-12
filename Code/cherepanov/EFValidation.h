#ifndef EFValidation_h
#define EFValidation_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"

class EFValidation : public Selection {

 public:
  EFValidation(TString Name_, TString id_);
  virtual ~EFValidation();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasMuon,
	     hasTau,
	     NCuts};



 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  virtual void  Finish();
 private:
  // Selection Variables


  TRandom *ran;

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;


  std::vector<TH1D> ZMass0;
  
  std::vector<TH1D> TruthA1Pt;
  std::vector<TH1D> TruthA1Phi;
  std::vector<TH1D> TruthA1Eta;
  
  std::vector<TH1D> TruthMuPt;
  std::vector<TH1D> TruthMuPhi;
  std::vector<TH1D> TruthMuEta;
  
  std::vector<TH1D> TruthTauMuPt;
  std::vector<TH1D> TruthTauMuPhi;
  std::vector<TH1D> TruhTauMuEta;
  
  std::vector<TH1D> TruthTauA1Pt;
  std::vector<TH1D> TruthTauA1Phi;
  std::vector<TH1D> TruthTauA1Eta;
 
  std::vector<TH1D> EFitTauA1Pt;
  std::vector<TH1D> EFitTauA1Phi;
  std::vector<TH1D> EFitTauA1Eta;

  std::vector<TH1D> EFitTauMuPt;
  std::vector<TH1D> EFitTauMuPhi;
  std::vector<TH1D> EFitTauMuEta;



  std::vector<TH1D> EFitTauA1PtAmbPoint;
  std::vector<TH1D> EFitTauA1PhiAmbPoint;
  std::vector<TH1D> EFitTauA1EtaAmbPoint;

  std::vector<TH1D> EFitTauMuPtAmbPoint;
  std::vector<TH1D> EFitTauMuPhiAmbPoint;
  std::vector<TH1D> EFitTauMuEtaAmbPoint;



  
  std::vector<TH1D> TauA1PtResolKFM;
  std::vector<TH1D> TauA1PhiResolKFM;
  std::vector<TH1D> TauA1EtaResolKFM;
  
 
  std::vector<TH1D> TauA1PtResolKFP;
  std::vector<TH1D> TauA1PhiResolKFP;
  std::vector<TH1D> TauA1EtaResolKFP;
  
 
  std::vector<TH1D> RecoA1Pt;
  std::vector<TH1D> RecoA1Phi;
  std::vector<TH1D> RecoA1Eta;
  
  std::vector<TH1D> RecoMuPt;
  std::vector<TH1D> RecoMuPhi;
  std::vector<TH1D> RecoMuEta;
 
  std::vector<TH1D> A1Mass;
  
  std::vector<TH1D> TruthA1PtAfterFit;
  std::vector<TH1D> TruthA1EAfterFit;
  std::vector<TH1D> TruthA1PhiAfterFit;
  std::vector<TH1D> TruthA1EtaAfterFit;
 
  std::vector<TH1D> TruthMuPtAfterFit;
  std::vector<TH1D> TruthMuEAfterFit;
  std::vector<TH1D> TruthMuPhiAfterFit;
  std::vector<TH1D> TruthMuEtaAfterFit;
 
 
 
 
  std::vector<TH1D> TauA1PtResolution;
  std::vector<TH1D> TauA1EResolution;
  std::vector<TH1D> TauA1PhiResolution;
  std::vector<TH1D> TauA1EtaResolution;
 
  std::vector<TH1D> TauMuPtResolution;
  std::vector<TH1D> TauMuEResolution;
  std::vector<TH1D> TauMuPhiResolution;
  std::vector<TH1D> TauMuEtaResolution;
 
 
  std::vector<TH2D> TauA1PhiResolVsProb;
  std::vector<TH2D> TauMuPhiResolVsProb;
  std::vector<TH2D> TauA1EtaResolVsProb;
  std::vector<TH2D> TauMuEtaResolVsProb;
  std::vector<TH2D> TauA1PtResolVsProb;
  std::vector<TH2D> TauMuPtResolVsProb;
  std::vector<TH2D> TauA1EResolVsProb;
  std::vector<TH2D> TauMuEResolVsProb;
      
  std::vector<TH2D> TauA1PtResolVsPt;
  std::vector<TH2D> TauMuPtResolVsPt;
  std::vector<TH2D> TauA1PtResolVsEta;
  std::vector<TH2D> TauMuPtResolVsEta;
  std::vector<TH2D> TauA1EtaResolVsEta;
  std::vector<TH2D> TauMuEtaResolVsEta;
 
 
  std::vector<TH2D> TauA1PtResolVsZPt;
  std::vector<TH2D> TauMuPtResolVsZPt;
 
  std::vector<TH2D> TauA1EtaResolVsZPt;
  std::vector<TH2D> TauMuEtaResolVsZPt;
  
 

  
 
  std::vector<TH2D> TauA1PhiResolVsPVSVSignificance;
  std::vector<TH2D> TauA1EtaResolVsPVSVSignificance;
  std::vector<TH2D> TauA1PtResolVsPVSVSignificance;
  std::vector<TH2D> TauA1EResolVsPVSVSignificance;
  
  std::vector<TH2D> TauA1PtResolVsPtAmbPoint;
  std::vector<TH2D> TauMuPtResolVsPtAmbPoint;
  std::vector<TH2D> TauA1PtResolVsEtaAmbPoint;
  std::vector<TH2D> TauMuPtResolVsEtaAmbPoint;
  std::vector<TH2D> TauA1EtaResolVsEtaAmbPoint;
  std::vector<TH2D> TauMuEtaResolVsEtaAmbPoint;
  
  std::vector<TH2D> TauA1PtResolVsZPtAmbPoint;
  std::vector<TH2D> TauMuPtResolVsZPtAmbPoint;
  
  std::vector<TH2D> TauA1EtaResolVsZPtAmbPoint;
  std::vector<TH2D> TauMuEtaResolVsZPtAmbPoint;
  
  std::vector<TH1D> TauA1PtResolutionAmbPoint;
  std::vector<TH1D> TauA1PhiResolutionAmbPoint;
  std::vector<TH1D> TauA1EtaResolutionAmbPoint;
  
  std::vector<TH1D> TauMuPtResolutionAmbPoint;
  std::vector<TH1D> TauMuPhiResolutionAmbPoint;
  std::vector<TH1D> TauMuEtaResolutionAmbPoint;
  
  std::vector<TH2D> TauA1PhiResolVsPhiRotSignificance;
  std::vector<TH2D> TauA1EtaResolVsPhiRotSignificance;
  std::vector<TH2D> TauA1PtResolVsPhiRotSignificance;
  std::vector<TH2D> TauA1EResolVsPhiRotSignificance;
  
  
  std::vector<TH2D> TauA1PhiResolVsPVSVSignificanceAmbPoint;
  std::vector<TH2D> TauA1EtaResolVsPVSVSignificanceAmbPoint;
  std::vector<TH2D> TauA1PtResolVsPVSVSignificanceAmbPoint;
  std::vector<TH2D> TauA1EResolVsPVSVSignificanceAmbPoint;

  std::vector<TH2D> Chi2Dim;
  std::vector<TH2D> csum2Dim;

  std::vector<TH2D> csumProb2Dim;
  std::vector<TH1D> ProbabilityOfCorrect;
  std::vector<TH1D> ProbabilityOfCorrectndf2;
  std::vector<TH1D> ProbabilityOfCorrectndf3;
  std::vector<TH1D> Chi2OfCorrect;
  std::vector<TH2D> ProbabilityOfCorrectVsZPt;

  std::vector<TH1D> idAllTaus;
  std::vector<TH1D> idIdentifiedTau;
  std::vector<TH1D> idPassedEF;
  std::vector<TH1D> CorrectAmbiguityTruth;

  std::vector<TH1D> MeasuredPhi;
  std::vector<TH1D> pvsvsignificance;
  std::vector<TH1D> PhiRotationSignificanceUnPhysicalTaus;
  std::vector<TH1D> ProbabilityOfAmbiguityPoint;


  std::vector<TH1D> EfficiencyOverA1Pt;
  std::vector<TH1D> EfficiencyOverA1Phi;
  std::vector<TH1D> EfficiencyOverA1Eta;


  std::vector<TH1D> EfficiencyOverMuPt;
  std::vector<TH1D> EfficiencyOverMuPhi;
  std::vector<TH1D> EfficiencyOverMuEta;

  //--- dependency plots
  std::vector<TH1D> TauA1PtResolVsPVSVSignificanceRMS;
  std::vector<TH1D> TauMuPtResolVsPtRMS;
  std::vector<TH1D> TauA1PtResolVsEtaRMS;
  std::vector<TH1D> TauMuPtResolVsZPtRMS;
  std::vector<TH1D> TauMuPtResolVsZPtMean;

  std::vector<TH1D> TauA1PtResolVsPtRMS;
  std::vector<TH2D> ProbabilityOfCorrectVsChi2;
  std::vector<TH1D> PhiRotSignificanceCutEfficiency;
  std::vector<TH1D> PVSVSignificanceCutEfficiency;
  std::vector<TH1D> ChannelEfficiency;


  //-- pull plots

  std::vector<TH1D> TauMuPxPull;
  std::vector<TH1D> TauMuPyPull;
  std::vector<TH1D> TauMuPzPull;

  std::vector<TH1D> TauA1PxPull;
  std::vector<TH1D> TauA1PyPull;
  std::vector<TH1D> TauA1PzPull;

  std::vector<TH1D> TauA1SFPxPull;
  std::vector<TH1D> TauA1SFPyPull;
  std::vector<TH1D> TauA1SFPzPull;

  std::vector<TH1D> SVxPull;
  std::vector<TH1D> SVyPull;
  std::vector<TH1D> SVzPull;



  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
