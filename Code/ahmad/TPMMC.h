#ifndef TPMMC_h
#define TPMMC_h
         
#include "Selection.h"
#include <vector>
#include "TString.h"

class TPMMC : public Selection {

 public:
  TPMMC(TString Name_, TString id_);
  virtual ~TPMMC();

  virtual void  Configure();
  virtual void  Finish();
  TF1 DeltaR_had(double pt);
  TF1 DeltaR_lep(double pt);
  enum cuts {isZtautautopimu=0,NCuts};
  
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  
 private:
  // Selection Variables
  std::vector<TH1D> LongitudinalPolarization;
  std::vector<TH1D> LongitudinalPolarization_Spin;

  std::vector<TH1D> spin_WT;

  //// Gen Taus /////
  std::vector<TH1D> M_TauPi_gen;
  std::vector<TH1D> Pt_TauPi_gen;
  std::vector<TH1D> Eta_TauPi_gen;
  std::vector<TH1D> Phi_TauPi_gen;
  std::vector<TH1D> Theta_TauPi_gen;
  
  std::vector<TH1D> M_TauMu_gen;   
  std::vector<TH1D> Pt_TauMu_gen;
  std::vector<TH1D> Eta_TauMu_gen;
  std::vector<TH1D> Phi_TauMu_gen;
  std::vector<TH1D> Theta_TauMu_gen;

  std::vector<TH1D> Mz_gen;
  std::vector<TH1D> Mz_gen1;

  //// Rec Taus /////
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
 
  /////////////////////  dR  ////////////////////////////////
  //Delta_R(Mu_rec , Mu_gen) , Delta_R(Pi_rec , Pi_gen)
  std::vector<TH1D> dR_recMu_genMu;
  std::vector<TH1D> mindR_recMu_genMu;
  std::vector<TH1D> dR_recPi_genPi;
  std::vector<TH1D> mindR_recPi_genPi;

  // MET infos
  std::vector<TH1D> MetX;
  std::vector<TH1D> MetY;
  std::vector<TH1D> MetX_gen;
  std::vector<TH1D> MetY_gen;
  
  // Di-Tau mass
  std::vector<TH1D> mtautau11_gen;
  std::vector<TH1D> mtautau12_gen;
  std::vector<TH1D> mtautau21_gen;
  std::vector<TH1D> mtautau22_gen;
  std::vector<TH1D> mtautau3_gen; 
  std::vector<TH1D> mtautau4_gen; 
  std::vector<TH1D> mtautau5_gen; 
  std::vector<TH1D> Lw_gen; 
  std::vector<TH1D> mtautau_gen; 

  // Di-Tau mass
  std::vector<TH1D> mtautau11;
  std::vector<TH1D> mtautau12;
  std::vector<TH1D> mtautau21;
  std::vector<TH1D> mtautau22;
  std::vector<TH1D> mtautau3; 
  std::vector<TH1D> mtautau4; 
  std::vector<TH1D> mtautau5; 
  std::vector<TH1D> Lw; 
  std::vector<TH1D> mtautau; 


  std::vector<TH2D> Mzgenmmc_gen;
  std::vector<TH2D> Mzrecmmc_gen;
  std::vector<TH2D> Mzrecmmc_genmmc;
  
     
  bool verbose; 
    
};
#endif
   
