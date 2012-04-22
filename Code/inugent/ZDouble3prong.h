#ifndef ZDouble3prong_h
#define ZDouble3prong_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ZDouble3prong : public Selection {

 public:
  ZDouble3prong(TString Name_, TString id_);
  virtual ~ZDouble3prong();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     NTauKinFit,
	     NTauPt,
	     NTauEta,
	     TauIso,
	     TauTauVertex,
	     TauChiProb,
             MET,
	     deltaPhi,
	     charge,
             ZMassmax,
             ZMassmin,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH2D> TagEtaPT;
  std::vector<TH1D> ZDouble3prongMT;

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

  double tau_pt,tau_eta,jet_pt,jet_eta,Iso_dr;

};
#endif
