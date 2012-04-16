#ifndef A_QCD_h
#define A_QCD_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class A_QCD : public Selection {

 public:
  A_QCD(TString Name_, TString id_); 
  virtual ~A_QCD();

  virtual void  Configure();
  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     MuonisGlob,
	     TauIsQuality,
	     MuonPt,
	     TauPtCut,
	     MET,
	     MuonIso,
	     charge,
	     TauIsIso,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
 
 private:
  // Selection Variables








 


};
#endif
