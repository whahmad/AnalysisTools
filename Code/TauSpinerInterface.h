#ifndef TauSpinerInterface_h
#define TauSpinerInterface_h

#include "SimpleParticle.h"
#include <vector>

using namespace TauSpinner;

class TauSpinerInterface {

 public:
  TauSpinerInterface();
  ~TauSpinerInterface();

  enum TauSpinerType {Spin=0,UnSpin,FlipSpin,hplus,hminus};
  enum TauSpinerSignalCharge {tauminus=-1,tauplus=1};

  double Get(TauSpinerType type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2);
  void SetTauSignalCharge(int tsc){signalcharge=tsc;}

 private:
  void Initialize();
  static int signalcharge;
  static bool initialized;

};
#endif
