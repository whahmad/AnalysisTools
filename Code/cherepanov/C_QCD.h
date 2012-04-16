#ifndef C_QCD_h
#define C_QCD_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class C_QCD : public Selection {

 public:
  C_QCD(TString Name_, TString id_); 
  virtual ~C_QCD();

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
