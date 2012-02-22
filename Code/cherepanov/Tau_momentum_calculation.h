#ifndef Tau_momentum_calculation_h
#define Tau_momentum_calculation_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class Tau_momentum_calculation : public Selection {

 public:
  Tau_momentum_calculation(TString Name_, TString id_);
  virtual ~Tau_momentum_calculation();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     MuonisGlob,
	     TauIsQuality,
	     MuonPt,
	     TauPt,
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
  std::vector<TH1D> ZTruthPt;
  std::vector<TH1D> Tau3PiPt;
  std::vector<TH1D> Tau1Pt;
  std::vector<TH1D> Tau2Pt;
  std::vector<TH1D> Tau1E;
  std::vector<TH1D> Tau2E;  
  std::vector<TH1D> Tau1theta;
  std::vector<TH1D> Tau2theta;  
  std::vector<TH1D> Tau1Deltatheta;
  std::vector<TH1D> Tau2Deltatheta;  
  std::vector<TH1D> DeltaTheta_1st;
  std::vector<TH1D> DeltaTheta_2nd;
  std::vector<TH1D> Energy_resolution1;
  std::vector<TH1D> Energy_resolution2;
  std::vector<TH1D> TransverseEnergy_resolution;
  std::vector<TH1D> TauMuEnergyRatio;
  std::vector<TH1D> TauMuTransverseEnergyRatio;

  std::vector<TH1D> TruthRatioTrTauVsTrMuon;
  std::vector<TH1D> TruthRatio;
  std::vector<TH1D> TruthResolution;
  std::vector<TH1D> TruthDeltaTheta;


  std::vector<TH1D> Energy_resolution_cuts;

  std::vector<TH2D> ResVsDeltaMu;

};
#endif
