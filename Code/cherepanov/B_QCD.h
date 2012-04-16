#ifndef B_QCD_h
#define B_QCD_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class B_QCD : public Selection {

 public:
  B_QCD(TString Name_, TString id_); 
  virtual ~B_QCD();

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
