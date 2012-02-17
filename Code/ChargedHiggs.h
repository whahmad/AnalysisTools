#ifndef ChargedHiggs_h
#define ChargedHiggs_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ChargedHiggs : public Selection {

 public:
  ChargedHiggs(TString Name_, TString id_);
  virtual ~ChargedHiggs();

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
	     ZMassmax,
             ZMassmin,
	     charge,
	     NCuts};

  enum tagtype{muontag,electrontag};

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
