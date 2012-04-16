#ifndef StudyResolution_h
#define StudyResolution_h

#include "Selection.h"
#include <vector>
#include "TString.h"

#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"




class StudyResolution : public Selection {

 public:
  StudyResolution(TString Name_, TString id_);
  virtual ~StudyResolution();

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

  std::vector<TH1D> KinTauRes;
  std::vector<TH1D> KinTauPhiRes;
  std::vector<TH1D> KinTauEtaRes;
  std::vector<TH2D> KinTauResVsE;

  std::vector<TH2D> KinTauResVsAmbig;
  std::vector<TH2D> KinResVsProb;
  std::vector<TH2D> KinResVsMuonTauRes;
  std::vector<TH1D> KinTaudeltaR;


  std::vector<TH2D> ResVsDeltaTheta;
  std::vector<TH2D> ResVsPt;
  std::vector<TH2D> ResVsSign;
  std::vector<TH2D> ResVsProb;
  std::vector<TH2D> ResVsDeltaTheta1;
  std::vector<TH2D> ResVsDeltaTheta2;
  std::vector<TH2D> ResVsDeltaThetas;
  std::vector<TH2D> ResVsAmbig;
  std::vector<TH2D> ResVsTTAmbig;

 



};
#endif
