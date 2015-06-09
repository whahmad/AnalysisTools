#ifndef MissingMassCalc_h
#define MissingMassCalc_h
      
#include "Selection.h"
#include <vector>
#include "TString.h"

class MissingMassCalc : public Selection {

 public:
  MissingMassCalc(TString Name_, TString id_);
  virtual ~MissingMassCalc();

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
  ////// Gen-Infos /////// 
 /////// Lab_Frame ///////
  // X1 = Pi_Pt / Tau_Pt
  std::vector<TH1D> pi_PtRatio;  
  std::vector<TH1D> pi_PtRatio_hplus;
  std::vector<TH1D> pi_PtRatio_hminus;
  std::vector<TH1D> pi_PtRatio_Spin;
  std::vector<TH1D> pi_PtRatio_hplus_Spin;
  std::vector<TH1D> pi_PtRatio_hminus_Spin;
  // X2 = Pi_Pt / Z_Mt
  std::vector<TH1D> pi_PtoverZmt;
  std::vector<TH1D> pi_PtoverZmt_hplus;
  std::vector<TH1D> pi_PtoverZmt_hminus;
  std::vector<TH1D> pi_PtoverZmt_Spin;
  std::vector<TH1D> pi_PtoverZmt_hplus_Spin;
  std::vector<TH1D> pi_PtoverZmt_hminus_Spin;
  
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  std::vector<TH1D> pi_EoverEtau;
  std::vector<TH1D> pi_EoverEtau_hplus;
  std::vector<TH1D> pi_EoverEtau_hminus;
  std::vector<TH1D> pi_EoverEtau_Spin;
  std::vector<TH1D> pi_EoverEtau_hplus_Spin;
  std::vector<TH1D> pi_EoverEtau_hminus_Spin;
  // X4 = Pi_E / Z_M
  std::vector<TH1D> pi_EoverZm;
  std::vector<TH1D> pi_EoverZm_hplus;
  std::vector<TH1D> pi_EoverZm_hminus;
  std::vector<TH1D> pi_EoverZm_Spin;
  std::vector<TH1D> pi_EoverZm_hplus_Spin;
  std::vector<TH1D> pi_EoverZm_hminus_Spin;  

  ////// Rec-Infos /////// 
  /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  std::vector<TH1D> pirec_PtRatio;  
  std::vector<TH1D> pirec_PtRatio_hplus;
  std::vector<TH1D> pirec_PtRatio_hminus;
  std::vector<TH1D> pirec_PtRatio_Spin;
  std::vector<TH1D> pirec_PtRatio_hplus_Spin;
  std::vector<TH1D> pirec_PtRatio_hminus_Spin;
  // X2 = Pirec_Pt / Zrec_Mt
  std::vector<TH1D> pirec_PtoverZmt;
  std::vector<TH1D> pirec_PtoverZmt_hplus;
  std::vector<TH1D> pirec_PtoverZmt_hminus;
  std::vector<TH1D> pirec_PtoverZmt_Spin;
  std::vector<TH1D> pirec_PtoverZmt_hplus_Spin;
  std::vector<TH1D> pirec_PtoverZmt_hminus_Spin;
  
  ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  std::vector<TH1D> pirec_EoverEtau;
  std::vector<TH1D> pirec_EoverEtau_hplus;
  std::vector<TH1D> pirec_EoverEtau_hminus;
  std::vector<TH1D> pirec_EoverEtau_Spin;
  std::vector<TH1D> pirec_EoverEtau_hplus_Spin;
  std::vector<TH1D> pirec_EoverEtau_hminus_Spin;
  // X4 = Pirec_E / Zrec_M
  std::vector<TH1D> pirec_EoverZm;
  std::vector<TH1D> pirec_EoverZm_hplus;
  std::vector<TH1D> pirec_EoverZm_hminus;
  std::vector<TH1D> pirec_EoverZm_Spin;
  std::vector<TH1D> pirec_EoverZm_hplus_Spin;
  std::vector<TH1D> pirec_EoverZm_hminus_Spin;  

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

  // DR(vis,mis)
 std::vector<TH1D> dR_Mu_MisNum_10_15;
 std::vector<TH1D> dR_Mu_MisNum_15_20;
 std::vector<TH1D> dR_Mu_MisNum_20_25;
 std::vector<TH1D> dR_Mu_MisNum_25_30;
 std::vector<TH1D> dR_Mu_MisNum_30_35;
 std::vector<TH1D> dR_Mu_MisNum_35_40;
 std::vector<TH1D> dR_Mu_MisNum_40_45;
 std::vector<TH1D> dR_Mu_MisNum_45_50;
 std::vector<TH1D> dR_Mu_MisNum_50_55;
 std::vector<TH1D> dR_Mu_MisNum_55_60;
 std::vector<TH1D> dR_Mu_MisNum_60_65;
 std::vector<TH1D> dR_Mu_MisNum_65_70;
 std::vector<TH1D> dR_Mu_MisNum_70_75;
 std::vector<TH1D> dR_Mu_MisNum_75_80;
 std::vector<TH1D> dR_Mu_MisNum_80_85;
 std::vector<TH1D> dR_Mu_MisNum_85_90;
 std::vector<TH1D> dR_Mu_MisNum_90_95;
 std::vector<TH1D> dR_Mu_MisNum_95_100;

 std::vector<TH1D> dR_Pi_MisNup_10_15;
 std::vector<TH1D> dR_Pi_MisNup_15_20;
 std::vector<TH1D> dR_Pi_MisNup_20_25;
 std::vector<TH1D> dR_Pi_MisNup_25_30;
 std::vector<TH1D> dR_Pi_MisNup_30_35;
 std::vector<TH1D> dR_Pi_MisNup_35_40;
 std::vector<TH1D> dR_Pi_MisNup_40_45;
 std::vector<TH1D> dR_Pi_MisNup_45_50;
 std::vector<TH1D> dR_Pi_MisNup_50_55;
 std::vector<TH1D> dR_Pi_MisNup_55_60;
 std::vector<TH1D> dR_Pi_MisNup_60_65;
 std::vector<TH1D> dR_Pi_MisNup_65_70;
 std::vector<TH1D> dR_Pi_MisNup_70_75;
 std::vector<TH1D> dR_Pi_MisNup_75_80;
 std::vector<TH1D> dR_Pi_MisNup_80_85;
 std::vector<TH1D> dR_Pi_MisNup_85_90;
 std::vector<TH1D> dR_Pi_MisNup_90_95;
 std::vector<TH1D> dR_Pi_MisNup_95_100; 


  // MET infos
  std::vector<TH1D> MetX;
  std::vector<TH1D> MetY;
  std::vector<TH1D> MetX_gen;
  std::vector<TH1D> MetY_gen;
  
  // Di-Tau mass
  std::vector<TH1D> mtautau11;
  std::vector<TH1D> mtautau12;
  std::vector<TH1D> mtautau21;
  std::vector<TH1D> mtautau22;
  std::vector<TH1D> mtautau3; 
  std::vector<TH1D> mtautau4; 
  std::vector<TH1D> mtautau5; 

  std::vector<TH1D> Lw; 
  //  std::vector<TH1D> mtt_est; 
  std::vector<TH1D> mtautau; 
     
  bool verbose; 
  
};
#endif
   
