#ifndef Tau_momentum_calculation_h
#define Tau_momentum_calculation_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Tau_momentum_calculation : public Selection {

 public:
  Tau_momentum_calculation(TString Name_, TString id_);
  virtual ~Tau_momentum_calculation();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     MuonisGlob,
	     TauIsQuality,
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
