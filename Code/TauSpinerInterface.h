#ifndef TauSpinerInterface_h
#define TauSpinerInterface_h

#include "SimpleParticle.h"
#include <vector>

class TauSpinerInterface {

 public:
  TauSpinerInterface();
  ~TauSpinerInterface();

  enum TauSpinerType {Spin=0,UnSpin,FlipSpin,LPolarization};
  double Get(TauSpinerType type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2);

};
#endif
