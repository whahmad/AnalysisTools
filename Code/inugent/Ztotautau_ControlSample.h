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
	     PrimeVtx,
	     hasTag,
	     TagPtmin,
             TagPtmax,
	     TagIso,
	     NJets,
	     JetPt,
	     deltaPhi,
	     MET,
	     MT,
	     TauAvgMETPhi,
	     PInBalance,
	     ZMassmax,
             ZMassmin,
	     charge,
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
  int channel;

};
#endif
