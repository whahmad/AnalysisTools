#ifndef Momentum_calculation_using_muon_angle_h
#define Momentum_calculation_using_muon_angle_h

#include "Selection.h"
#include <vector>
#include "TString.h"


class Momentum_calculation_using_muon_angle : public Selection {

 public:
  Momentum_calculation_using_muon_angle(TString Name_, TString id_); 
  virtual ~Momentum_calculation_using_muon_angle();

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
  std::vector<TH1D> ZTruthPt;
  std::vector<TH1D> TauPt;
  std::vector<TH1D> Tautheta;
  std::vector<TH1D> Energy_resolution;
  std::vector<TH1D> TransverseEnergy_resolution;
  std::vector<TH1D> TauMuEnergyRatio;
  std::vector<TH1D> TauMuTransverseEnergyRatio;
  std::vector<TH1D> Energy_resolution_cuts;
  std::vector<TH1D> ZMass; 
  std::vector<TH1D> resol;
  std::vector<TH1D> TruthRecoiledZ;  
  std::vector<TH1D> TruthRatio;  
  std::vector<TH1D> TruthTransverseRatio;  



  std::vector<TH2D> ResVsZPt;

 


};
#endif
