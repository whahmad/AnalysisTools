#ifndef Ztotautau_hadmu_ControlSample_h
#define Ztotautau_hadmu_ControlSample_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Ztotautau_hadmu_ControlSample : public Selection {

 public:
  Ztotautau_hadmu_ControlSample(TString Name_, TString id_);
  virtual ~Ztotautau_hadmu_ControlSample();

  virtual void  Configure();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
	     MuonisGlob,
	     MuonPt,
	     TauPt,
	     TauIsRef,
	     MuonIso,
	     TauIsIso,
	     MET,
	     deltaPhi,
	     ZMassV,
	     ZMassHPS,
	     tauPhi,
	     charge,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;



};
#endif
