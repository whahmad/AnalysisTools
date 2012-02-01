#ifndef TauSpinerInterface_h
#define TauSpinerInterface_h

#include "SimpleParticle.h"
#include <vector>

class TauSpinerInterface {

 public:
  TauSpinerInterface();
  ~TauSpinerInterface();

  enum TauSpinerType {Spin=0,UnSpin,FlipSpin,LPolarization,hplus,hminus};
  double Get(TauSpinerType type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2);

 private:
  double tautauHelicityState(SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2,std::vector<SimpleParticle> &sp_tau1_daughters, std::vector<SimpleParticle> &sp_tau2_daughters);
};
#endif
