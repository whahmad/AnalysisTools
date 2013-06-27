#ifndef TauLifeTime_h
#define TauLifeTime_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TauLifeTime : public Selection {

 public:
  TauLifeTime(TString Name_, TString id_);
  virtual ~TauLifeTime();

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
	     JetTrackPtMax,
	     ZMassmax,
             ZMassmin,
	     HT,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;

  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
