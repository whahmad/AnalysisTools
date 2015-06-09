#ifndef TauSpinStudy_h
#define TauSpinStudy_h
      
#include "Selection.h"
#include <vector>
#include "TString.h"

class TauSpinStudy : public Selection {

 public:
  TauSpinStudy(TString Name_, TString id_);
  virtual ~TauSpinStudy();

  virtual void  Configure();
  virtual void  Finish();
  enum cuts {isZtautautopimu=0,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> LongitudinalPolarization;
  std::vector<TH1D> LongitudinalPolarization_Spin;
  std::vector<TH1D> LongitudinalPolarization_UnSpin;
  std::vector<TH1D> LongitudinalPolarization_FlipSpin;

  std::vector<TH1D> mu_ExoverEtau;
  std::vector<TH1D> mu_ExoverEtau_hplus;
  std::vector<TH1D> mu_ExoverEtau_hminus;
  std::vector<TH1D> mu_ExoverEtau_Spin;
  std::vector<TH1D> mu_ExoverEtau_UnSpin;
  std::vector<TH1D> mu_ExoverEtau_FlipSpin;
  std::vector<TH1D> mu_WT_Spin;
  std::vector<TH1D> mu_WT_UnSpin;
  std::vector<TH1D> mu_WT_FlipSpin;
  std::vector<TH1D> mu_PtRatio_hplus;
  std::vector<TH1D> mu_PtRatio_hminus;

  std::vector<TH1D> murec_ExoverEtaurec;
  std::vector<TH1D> murec_ExoverEtaurec_hplus;
  std::vector<TH1D> murec_ExoverEtaurec_hminus;
  std::vector<TH1D> murec_ExoverEtaurec_Spin;
  std::vector<TH1D> murec_ExoverEtaurec_UnSpin;
  std::vector<TH1D> murec_ExoverEtaurec_FlipSpin;
  std::vector<TH1D> murec_WT_Spin;
  std::vector<TH1D> murec_WT_UnSpin;
  std::vector<TH1D> murec_WT_FlipSpin;
  std::vector<TH1D> murec_PtRatio_hplus;
  std::vector<TH1D> murec_PtRatio_hminus;
 
  std::vector<TH1D> pi_ExoverEtau;
  std::vector<TH1D> pi_ExoverEtau_hplus;
  std::vector<TH1D> pi_ExoverEtau_hminus;
  std::vector<TH1D> pi_ExoverEtau_Spin;
  std::vector<TH1D> pi_ExoverEtau_UnSpin;
  std::vector<TH1D> pi_ExoverEtau_FlipSpin;
  std::vector<TH1D> pi_WT_Spin;
  std::vector<TH1D> pi_WT_UnSpin;
  std::vector<TH1D> pi_WT_FlipSpin;
  std::vector<TH1D> pi_PtRatio_hplus;
  std::vector<TH1D> pi_PtRatio_hminus;

  std::vector<TH1D> pirec_ExoverEtaurec;
  std::vector<TH1D> pirec_ExoverEtaurec_hplus;
  std::vector<TH1D> pirec_ExoverEtaurec_hminus;
  std::vector<TH1D> pirec_ExoverEtaurec_Spin;
  std::vector<TH1D> pirec_ExoverEtaurec_UnSpin;
  std::vector<TH1D> pirec_ExoverEtaurec_FlipSpin;
  std::vector<TH1D> pirec_WT_Spin;
  std::vector<TH1D> pirec_WT_UnSpin;
  std::vector<TH1D> pirec_WT_FlipSpin;
  std::vector<TH1D> pirec_PtRatio_hplus;
  std::vector<TH1D> pirec_PtRatio_hminus;

  std::vector<TH1D> pi_zs;
  std::vector<TH1D> pi_zs_hplus;
  std::vector<TH1D> pi_zs_hminus;
  std::vector<TH1D> pi_zs_Spin;
  std::vector<TH1D> pi_zs_UnSpin;
  std::vector<TH1D> pi_zs_FlipSpin;

  std::vector<TH1D> pi_Mvis;
  std::vector<TH1D> pi_Mvis_hplus;
  std::vector<TH1D> pi_Mvis_hminus;
  std::vector<TH1D> pi_Mvis_Spin;
  std::vector<TH1D> pi_Mvis_UnSpin;
  std::vector<TH1D> pi_Mvis_FlipSpin;

  // Selection Variables
  //Reco objects
  std::vector<TH1D> Mu_pt;
  std::vector<TH1D> Mu_eta;
  std::vector<TH1D> Mu_phi;
  std::vector<TH1D> Mu_theta;
  
  std::vector<TH1D> Pion_pt;
  std::vector<TH1D> Pion_eta; 
  std::vector<TH1D> Pion_phi; 
  std::vector<TH1D> Pion_theta; 
  
  std::vector<TH1D> MET_eta;
  std::vector<TH1D> MET_phi;

  std::vector<TH1D> PionRec_pt;
  std::vector<TH1D> PionRec_eta;
  std::vector<TH1D> PionRec_phi;
  std::vector<TH1D> Pion_mdR_pt;
  std::vector<TH1D> Pion_mdR_eta;
  std::vector<TH1D> Pion_mdR_phi;

  std::vector<TH1D> MuonRec_pt;
  std::vector<TH1D> MuonRec_eta;
  std::vector<TH1D> MuonRec_phi;
  std::vector<TH1D> Muon_mdR_pt;
  std::vector<TH1D> Muon_mdR_eta;
  std::vector<TH1D> Muon_mdR_phi;
  
  std::vector<TH1D> M_TauPi_rec;
  std::vector<TH1D> Pt_TauPi_rec;
  std::vector<TH1D> Eta_TauPi_rec;
  std::vector<TH1D> Phi_TauPi_rec;
  std::vector<TH1D> Theta_TauPi_rec;
  
  std::vector<TH1D> M_TauMu_rec;
  std::vector<TH1D> Pt_TauMu_rec;
  std::vector<TH1D> Eta_TauMu_rec;
  std::vector<TH1D> Phi_TauMu_rec;
  std::vector<TH1D> Theta_TauMu_rec;
  
  //Gen objects
  //Z particle
  std::vector<TH1D> M_Z_gen;
  std::vector<TH1D> Pt_Z_gen;
  std::vector<TH1D> Eta_Z_gen;
  std::vector<TH1D> Phi_Z_gen;
  std::vector<TH1D> Theta_Z_gen;
  //Z particle
  std::vector<TH1D> M_Z1_gen;
  std::vector<TH1D> Pt_Z1_gen;
  std::vector<TH1D> Eta_Z1_gen;
  std::vector<TH1D> Phi_Z1_gen;
  std::vector<TH1D> Theta_Z1_gen;
  //Z-> TauPi
  std::vector<TH1D> M_TauPi_gen;
  std::vector<TH1D> Pt_TauPi_gen;
  std::vector<TH1D> Eta_TauPi_gen;
  std::vector<TH1D> Phi_TauPi_gen;
  std::vector<TH1D> Theta_TauPi_gen;
  //Z-> TauMu
  std::vector<TH1D> M_TauMu_gen;
  std::vector<TH1D> Pt_TauMu_gen;
  std::vector<TH1D> Eta_TauMu_gen;
  std::vector<TH1D> Phi_TauMu_gen;
  std::vector<TH1D> Theta_TauMu_gen;
  // TauPi->Pi
  std::vector<TH1D> M_Pi_gen;
  std::vector<TH1D> Pt_Pi_gen;
  std::vector<TH1D> Eta_Pi_gen;
  std::vector<TH1D> Phi_Pi_gen;
  std::vector<TH1D> Theta_Pi_gen;
  //TauMu->Mu
  std::vector<TH1D> M_Mu_gen;
  std::vector<TH1D> Pt_Mu_gen; 
  std::vector<TH1D> Eta_Mu_gen;
  std::vector<TH1D> Phi_Mu_gen;
  std::vector<TH1D> Theta_Mu_gen;
  //TauPi-> Nutp
  std::vector<TH1D> M_Nutp_gen;
  std::vector<TH1D> Pt_Nutp_gen;
  std::vector<TH1D> Eta_Nutp_gen;
  std::vector<TH1D> Phi_Nutp_gen;
  std::vector<TH1D> Theta_Nutp_gen;
  //TauMu-> Nutm
  std::vector<TH1D> M_Nutm_gen;
  std::vector<TH1D> Pt_Nutm_gen; 
  std::vector<TH1D> Eta_Nutm_gen;
  std::vector<TH1D> Phi_Nutm_gen;
  std::vector<TH1D> Theta_Nutm_gen;
  //TauMu-> Num
  std::vector<TH1D> M_Num_gen;
  std::vector<TH1D> Pt_Num_gen;
  std::vector<TH1D> Eta_Num_gen;
  std::vector<TH1D> Phi_Num_gen; 
  std::vector<TH1D> Theta_Num_gen;
    
  // DiTau (Z-> TauPi+TauMu) 
  std::vector<TH1D> M_DiTau_gen;
  std::vector<TH1D> Pt_DiTau_gen;
  std::vector<TH1D> dEta_DiTau_gen;
  std::vector<TH1D> dPhi_DiTau_gen;
  std::vector<TH1D> dTheta_DiTau_gen;
  std::vector<TH1D> dM_DiTau_gen;
  std::vector<TH1D> dR_DiTau_gen;
  // Physics objects(Pi; std::vector<TH1D>Mu) 
  std::vector<TH1D> M_MuPi_gen; 
  std::vector<TH1D> Pt_MuPi_gen;
  std::vector<TH1D> dEta_MuPi_gen;
  std::vector<TH1D> dPhi_MuPi_gen;
  std::vector<TH1D> dTheta_MuPi_gen;
  std::vector<TH1D> dR_MuPi_gen;
  
  //Delta_R(TauPi_rec , TauPi_gen)
  //std::vector<TH1D> dR_recTauPi_genTauPi;
  //std::vector<TH1D> dR_recTauMu_genTauMu;

  //Delta_R(Mu_rec , Mu_gen) , Delta_R(Pi_rec , Pi_gen)
  std::vector<TH1D> dR_recMu_genMu;
  std::vector<TH1D> mindR_recMu_genMu;
  std::vector<TH1D> dR_recPi_genPi;
  std::vector<TH1D> mindR_recPi_genPi;

  bool verbose; 

  double Zstoa(double zs);

  int zsbins;
  float zsmin,zsmax;


};
#endif
   
