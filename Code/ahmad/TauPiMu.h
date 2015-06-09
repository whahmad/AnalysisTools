#ifndef TauPiMu_h
#define TauPiMu_h

#include "Selection.h"     
#include <vector>
#include "TString.h"

class TauPiMu : public Selection {
  
 public:
  TauPiMu(TString Name_, TString id_);
  virtual ~TauPiMu();
  
  virtual void  Configure();
  virtual void  Finish();     
  
  enum cuts {ZTTtoPiMu=0,NCuts};
  
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  
 private:
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
  
  /* std::vector<TH1D> Tau_pt; */
  /* std::vector<TH1D> Tau_eta; */
  /* std::vector<TH1D> Tau_phi; */
  /* std::vector<TH1D> Tau_theta; */
  
  std::vector<TH1D> MET_eta;
  std::vector<TH1D> MET_phi;
  
  //Gen objects
  //Z particle
  std::vector<TH1D> M_Z_gen;
  std::vector<TH1D> Pt_Z_gen;
  std::vector<TH1D> Eta_Z_gen;
  std::vector<TH1D> Phi_Z_gen;
  std::vector<TH1D> Theta_Z_gen;
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
  std::vector<TH1D> dR_recPi_genPi;
  

};
#endif
 
