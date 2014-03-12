// Code written by Vladimir Cherepanov
// RWTH Aachen
#ifndef EFExample_h
#define EFExample_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TRandom.h"

class EFExample : public Selection {

 public:
  EFExample(TString Name_, TString id_);
  virtual ~EFExample();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasMuon,
	     hasTau,
	     NCuts};



 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
 private:
  // Selection Variables



  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;


  
  std::vector<TH1D> TauA1PtPhysical;
  std::vector<TH1D> TauA1PtAmbiguityPoint;

  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;


};
#endif
