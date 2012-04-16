#ifndef TemplateQCD_Using_KinematicFit_h
#define TemplateQCD_Using_KinematicFit_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TauSpinerInterface.h"

#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"


class TemplateQCD_Using_KinematicFit : public Selection {

 public:
  TemplateQCD_Using_KinematicFit(TString Name_, TString id_); 
  virtual ~TemplateQCD_Using_KinematicFit();

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

  std::vector<TH1D> NGoodVtx;

  std::vector<TH1D> xmuon;
  std::vector<TH1D> xmuonT;

 


};
#endif
