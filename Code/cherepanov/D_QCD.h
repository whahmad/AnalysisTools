#ifndef D_QCD_h
#define D_QCD_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class D_QCD : public Selection {

 public:
  D_QCD(TString Name_, TString id_); 
  virtual ~D_QCD();

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
