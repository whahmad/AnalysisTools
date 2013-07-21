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
	     TagPtMin,
	     TagIso,
//	     numIsoTags,
	     TauPt,
	     TauEta,
	     TauIsIsolated,
         TauFit,
	     deltaPhi,
	     Charge,
	     //ZMassmax,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> TauFlightLength;
  std::vector<TH1D> TauFlightLengthTransverse;
  std::vector<TH1D> TauMomentum;
  std::vector<TH1D> TauMomentumTransverse;
  std::vector<TH1D> TauLife;
  std::vector<TH1D> TauLifeTransverse;
  std::vector<TH1D> ResTauFlightLength;
  std::vector<TH1D> ResTauFlightLengthTransverse;
  std::vector<TH1D> ResTrueTauMomentum;
  std::vector<TH1D> ResTauMomentumTransverse;
  std::vector<TH1D>  ResTauMomentum;
  std::vector<TH1D> ResTrueTauMomentumTransverse;
  std::vector<TH1D> ResTrueTauLife;
  std::vector<TH1D> ResTrueTauLifeTransverse;
  std::vector<TH1D> ResTauLife;
  std::vector<TH1D> ResTauLifeTransverse;

  std::vector<TH1D> TauMomentumAvg;
  std::vector<TH1D> TauMomentumTransverseAvg;
  std::vector<TH1D> ResTauMomentumAvg;
  std::vector<TH1D> ResTauMomentumTransverseAvg;
  std::vector<TH1D> TauLifeAvg;
  std::vector<TH1D> TauLifeTransverseAvg;
  std::vector<TH1D> ResTauLifeAvg;
  std::vector<TH1D> ResTauLifeTransverseAvg;

  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
