#ifndef TemplateSignal_Using_KinematicFit_h
#define TemplateSignal_Using_KinematicFit_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"




class TemplateSignal_Using_KinematicFit : public Selection {

 public:
  TemplateSignal_Using_KinematicFit(TString Name_, TString id_);
  virtual ~TemplateSignal_Using_KinematicFit();

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

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> Tau3PiPt;
  std::vector<TH1D> Tau1theta;
  std::vector<TH1D> Tau2theta;  
  std::vector<TH1D> Tau1Deltatheta;
  std::vector<TH1D> Tau2Deltatheta;  
  std::vector<TH1D> DeltaTheta_1st;
  std::vector<TH1D> DeltaTheta_2nd;
  std::vector<TH1D> Tau1Pt;
  std::vector<TH1D> Tau2Pt;
  std::vector<TH1D> Energy_resolution1;
  std::vector<TH1D> Energy_resolution2;
  std::vector<TH1D> TruthResolution;
  std::vector<TH1D> TransverseEnergy_resolution;

  std::vector<TH1D> xmuon_plus;
  std::vector<TH1D> xmuonT_plus;
  std::vector<TH1D> xmuonTruthT_plus;
  std::vector<TH1D> xmuonTruth_plus;
  std::vector<TH1D> xmuonMC_plus;
  std::vector<TH1D> xmuonTMC_plus;



  std::vector<TH1D> xmuon_minus;
  std::vector<TH1D> xmuonT_minus;
  std::vector<TH1D> xmuonTruthT_minus;
  std::vector<TH1D> xmuonTruth_minus;
  std::vector<TH1D> xmuonMC_minus;
  std::vector<TH1D> xmuonTMC_minus;   

  
  std::vector<TH1D> TruthDeltaTheta;




};
#endif
