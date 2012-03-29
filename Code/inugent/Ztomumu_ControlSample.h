#ifndef Ztomumu_ControlSample_h
#define Ztomumu_ControlSample_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Ztomumu_ControlSample : public Selection {

 public:
  Ztomumu_ControlSample(TString Name_, TString id_);
  virtual ~Ztomumu_ControlSample();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     NMu,
	     NMuPt,
	     NMuEta,
	     MuIso,
	     MuMuVertex,
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
    
  double mu_pt,mu_eta,jet_pt,jet_eta;

};
#endif
