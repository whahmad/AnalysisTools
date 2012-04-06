#ifndef Ztotautau_ControlSample_h
#define Ztotautau_ControlSample_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Ztotautau_ControlSample : public Selection {

 public:
  Ztotautau_ControlSample(TString Name_, TString id_);
  virtual ~Ztotautau_ControlSample();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasTag,
	     TagPtmin,
             TagPtmax,
	     TagIso,
	     NJets,
	     MaxTracksinJet,
	     MinTracksinJet,
	     charge,
	     deltaPhi,
	     MT,
	     MET,
	     JetPt,
	     etaq,
	     sumcosdeltaphi,
	     ZMassmax,
             ZMassmin,
	     HT,
	     NCuts};

  enum Channel{muontag,electontag,rhotag,threepiontag,NChannels};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH2D> TagEtaPT;

  std::vector<TH1D> TauCandFound;
  std::vector<TH2D> TauCandEtaPhi;

  std::vector<TH1D> TauCandPhi;
  std::vector<TH1D> TauCandPhiRes;
  std::vector<TH1D> TauCandEta;
  std::vector<TH1D> TauCandEtaRes;
  std::vector<TH1D> TauCandE;
  std::vector<TH1D> TauCandERes;

  std::vector<TH1D> TauCandMass;
  std::vector<TH1D> TauCandPiMass;
  std::vector<TH1D> TauCandNuMass;

  std::vector<TH1D> TauSolutionResult;
  std::vector<TH1D> EstimatedTauE;
  std::vector<TH1D> EstimatedTauPhi;
  std::vector<TH1D> EstimatedTauEta;

  std::vector<TH1D> EstimatedTauERes;
  std::vector<TH1D> EstimatedTauPhiRes;
  std::vector<TH1D> EstimatedTauEtaRes;

  std::vector<TH1D> KFTau_Fit_chiprob;
  std::vector<TH1D> KFTau_Fit_a1mass;
  std::vector<TH1D> KFTau_Fit_chi2;
  std::vector<TH1D> KFTau_Fit_ndf;
  std::vector<TH1D> KFTau_Fit_ambiguity;
  std::vector<TH1D> KFTau_Fit_csum;
  std::vector<TH1D> KFTau_Fit_iterations;
  std::vector<TH1D> KFTau_Fit_TauEnergyFraction;
  std::vector<TH1D> KFTau_Fit_PV_PV_significance;
  std::vector<TH1D> KFTau_Fit_SV_PV_significance;

  int channel;
  double jeteta;
};
#endif
