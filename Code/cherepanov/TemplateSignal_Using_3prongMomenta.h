#ifndef TemplateSignal_Using_3prongMomenta_h
#define TemplateSignal_Using_3prongMomenta_h

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


class TemplateSignal_Using_3prongMomenta : public Selection {

 public:
  TemplateSignal_Using_3prongMomenta(TString Name_, TString id_); 
  virtual ~TemplateSignal_Using_3prongMomenta();

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
  std::vector<TH1D> TauPt;
  std::vector<TH1D> ZMass; 

  std::vector<TH1D> xmuon_plus;
  std::vector<TH1D> xmuon_minus;
  std::vector<TH1D> xmuon_plus_truth;
  std::vector<TH1D> xmuon_minus_truth;
  std::vector<TH1D> xmuon_plus_MC;
  std::vector<TH1D> xmuon_minus_MC;    



  std::vector<TH1D> xmuonT_plus;
  std::vector<TH1D> xmuonT_minus;
  std::vector<TH1D> xmuonT_plus_truth;
  std::vector<TH1D> xmuonT_minus_truth;
  std::vector<TH1D> xmuonT_plus_MC;
  std::vector<TH1D> xmuonT_minus_MC;    

  std::vector<TH1D> ResolMu;
  std::vector<TH1D> ResolMuTruth;
  std::vector<TH1D> ResolMuT;
  std::vector<TH1D> ResolMuTruthT;



 


};
#endif
