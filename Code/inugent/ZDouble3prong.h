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
	     NJets,
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
  double tau_pt,tau_eta,jet_pt,jet_eta;

};
#endif
