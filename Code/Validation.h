#ifndef Validation_h
#define Validation_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Validation : public Selection {

 public:
  Validation(TString Name_, TString id_);
  virtual ~Validation();

  virtual void  Configure();

  enum cuts {TriggerOk=0,PrimeVtx,GoodMuon,TauPt,QC,LooseIso,ET,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:


  // Selection Variables
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> MuonM;
  std::vector<TH1D> SecondMuonPt;

  std::vector<TH1D> MET_et;
  std::vector<TH1D> MET_phi;
  std::vector<TH1D> KFTau_TauFitPt;
  std::vector<TH1D> PFTau_DecMode;

  std::vector<TH1D> discrQC;          
  std::vector<TH1D> discrFT;          




  std::vector<TH1D> Muon_emEt03_hist;
  std::vector<TH1D> Muon_hadEt03_hist;
  std::vector<TH1D> Muon_sumPt03_hist;
  std::vector<TH1D> Muon_emEt05_hist;
  std::vector<TH1D> Muon_hadEt05_hist;
  std::vector<TH1D> Muon_sumPt05_hist;
  std::vector<TH1D> Muon_emVetoEt05_hist;
  std::vector<TH1D> Muon_hadVetoEt05_hist;
  std::vector<TH1D> Muon_sumVetoPt05_hist;

  std::vector<TH1D> Muon_isGlobalMuon_hist;
  std::vector<TH1D> Muon_isStandAloneMuon_hist;

  std::vector<TH1D> Muon_nJets05_hist;	
  std::vector<TH1D> Muon_nTracks05_hist;



  std::vector<TH1D> MuonIsol03;
  std::vector<TH1D> MuonIsol05;          

  std::vector<TH1D> TransvereMass;          


  std::vector<TH1D> NFitTaus;          
  std::vector<TH1D> NQCTaus;          
  std::vector<TH1D> NGlobalMuons;           

  std::vector<TH1D> PFTau_isTightIsolation_hist;
  std::vector<TH1D> PFTau_isMediumIsolation_hist;
  std::vector<TH1D> PFTau_isLooseIsolation_hist;
  std::vector<TH1D> PFTau_hpsDecayMode_hist;
  std::vector<TH1D> PFTau_Charge_hist;

 std::vector<TH1D> VisibleMass;     
 std::vector<TH1D> FullMass;
 std::vector<TH1D> DeltaPhi;  
 std::vector<TH1D> DeltaRTauMu;
 std::vector<TH1D> DeltaPhiTauMuNeutrino;
 std::vector<TH1D> DeltaVtxZ;
 std::vector<TH1D> GlobMuonEta;
 std::vector<TH1D> TauEta;


 std::vector<TH1D> KFTau_Fit_iterations_hist;
 std::vector<TH1D> KFTau_Fit_chi2_hist;
 std::vector<TH1D> KFTau_TauFit_pt_hist;
 std::vector<TH1D> KFTau_TauFit_phi_hist;
 std::vector<TH1D> KFTau_TauVis_phi_hist;
 std::vector<TH1D> KFTau_Neutrino_p4_hist;
 std::vector<TH1D> KFTau_Neutrino_phi_hist;
 std::vector<TH1D> KFTau_TauVis_a1_hist;
 std::vector<TH1D> KFTau_Fit_csum_hist;
 std::vector<TH1D> KFTau_TauFit_Highestpt_hist;
 std::vector<TH1D> KFTau_TauFit_Highestpt_phi_hist;

 std::vector<TH1D> PFTau_Tau_phi_hist;
 std::vector<TH1D> PFTau_Tau_pt_hist;
 std::vector<TH1D> PFTau_Tau_Highestpt_hist;
 std::vector<TH1D> PFTau_Tau_Highestpt_phi_hist;

 
 std::vector<TH1D> KinFitTau_HPSMode_hist;
 std::vector<TH1D> KinFitTau_HPSLooseIso_hist;
 std::vector<TH1D> KinFitTau_HPSMediumIso_hist;


};
#endif
