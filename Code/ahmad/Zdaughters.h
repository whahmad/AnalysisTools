#ifndef Zdaughters_h
#define Zdaughters_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "ReferenceScaleFactors.h"

class Zdaughters : public Selection {

 public:
  Zdaughters(TString Name_, TString id_);
  virtual ~Zdaughters();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
	     NMuId,
	     NMuKin,
	     NMuIso,
	     NTauId,
	     NTauKin,
	     NTauIso,
	     ChargeSum,
	     MT_MuMET,
	     TauDecayMode,
	     TauFLSigma,
	     NCuts};
  
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  /* ReferenceScaleFactors *RSF; */

  /* // cut values */
  /* double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta, cMu_dRHltMatch; */
  /* double cTau_pt, cTau_eta, cMuTau_dR, cTau_IsoRaw, cTau_dRHltMatch; */
  /* std::vector<TString> cTriggerNames; */

  /* double OneProngNoPiWeight; */
  /* //dummy values */
  /* int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy; */
  /* double MTDummy, MvisDummy, TauFLSigmaDummy; */
   
  /* int Charge; */
  /* double SB_lowerLimit, SB_upperLimit; */
  /* bool Scaleby_Counting; */
  /* TString tau_corr; */
  /* bool verbose, Use_Embedded; */
  /* bool selectMuon_Id(unsigned i, unsigned vertex); */
  /* bool selectMuon_Kinematics(unsigned i); */
  /* bool selectMuon_Isolation(unsigned i); */
  /* bool selectMuon_AntiIsolation(unsigned i); */

  /* bool selectPFTau_Id(unsigned i); */
  /* bool selectPFTau_Id(unsigned i, std::vector<int> muonCollection); */
  /* bool selectPFTau_Isolation(unsigned i); */
  /* bool selectPFTau_Kinematics(unsigned i); */
  /* double Reconstruct_hadronicTauEnergy(unsigned i); */


  // Selection Variables
  
  // Histograms
   
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  
  /* std::vector<TH1D> NMtTauMET; */
  /* std::vector<TH1D> NMvis, Mvis_SignalOnly, Mvis_SignalOnly_genMu, Mvis_SignalOnly_genA1, Mvis_SignalOnly_genTaumu, Mvis_SignalOnly_genTauh; */
  /* std::vector<TH1D> NSignal_SB_WJets; */
  /* std::vector<TH1D> NSB_Data; */
  /* std::vector<TH1D> NSB; */
 
  /* std::vector<std::vector<TH1D>* > Extradist1d_OS, Extradist1d_SS; */

  //Reco objects
  std::vector<TH1D> Mu_pt, Mu_eta, Mu_phi, Mu_theta, Tau_pt, Tau_eta, Tau_phi,Tau_theta, MET_eta,  MET_phi;

  /* std::vector<TH1D> TauFL, TauFLSigned, TauFLSigmaSigned, TauFLSigmaUnsigned; */
  /* std::vector<TH1D> A1mass, A1mass10GeV; */

  /* std::vector<TH1D> Mvis3Prong, Mvis1Prong, MvisIncl; */
  /* std::vector<TH1D> MTMuMET3Prong, MTMuMET1Prong, MTMuMETIncl; */

  /* //Gen objects */
  /* std::vector<TH1D> dR_selTauh_genTauh, dR_selMu_genMu; */
  /* std::vector<TH1D> POCAPV_Mag; */
  /* std::vector<TH1D> Phi_SVPV, Phi_genTauh, Theta_SVPV, Theta_genTauh, dPhi_SVPV_genTauh, dTheta_SVPV_genTauh, Angle_SVPV_genTauh; */
  /* std::vector<TH1D> Phi_POCAPV, Phi_genTaumu, Theta_POCAPV, Theta_genTaumu, dPhi_POCAPV_genTaumu, dTheta_POCAPV_genTaumu; */
  /* std::vector<TH1D> dPhi_MinusSVPV_genTaumu, dTheta_MinusSVPV_genTaumu, Angle_MinusSVPV_genTaumu; */
  /* std::vector<TH1D> GJ_Tauh, GJ_Taumu; */
  /* std::vector<TH1D> dPhi_DiTauGen, Pt_DiTauGen, Pt_ZGen, M_ZGen, M_DiTauPtBalance, dM_DiTau, dPt_GenTaumuPtBalance, dP_GenTaumuPtBalance, dP_GenTauh; */
  /* std::vector<TH2D> dP_GenTauMuPtBalance_vs_dPTauh, Pt_vs_dPhi_DiTauGen; */
  /* std::vector<TH2D> TauFLSigmaCut_vs_Res, TauFLSigma_vs_Res; */

  /* std::vector<TH1D> NQCD; */
  /* std::vector<TH1D> QCD_MT_MuMET_A, QCD_MT_MuMET_B, QCD_MT_MuMET_C, QCD_MT_MuMET_D; */
  /* std::vector<TH1D> QCD_MET_A, QCD_MET_B, QCD_MET_C, QCD_MET_D; */

   
  //Gen objects
  //Z particle
  std::vector<TH1D> M_Z_gen,     Pt_Z_gen,      Eta_Z_gen,     Phi_Z_gen,      Theta_Z_gen;
  //Z-> TauPi
  std::vector<TH1D> M_TauPi_gen, Pt_TauPi_gen,  Eta_TauPi_gen, Phi_TauPi_gen,  Theta_TauPi_gen;
  //Z-> TauMu
  std::vector<TH1D> M_TauMu_gen, Pt_TauMu_gen,  Eta_TauMu_gen, Phi_TauMu_gen,  Theta_TauMu_gen;
  // TauPi->Pi
  std::vector<TH1D> M_Pi_gen,    Pt_Pi_gen,     Eta_Pi_gen,    Phi_Pi_gen,     Theta_Pi_gen;
  //TauMu->Mu
  std::vector<TH1D> M_Mu_gen,    Pt_Mu_gen,     Eta_Mu_gen,    Phi_Mu_gen,     Theta_Mu_gen;
  //TauPi-> Nutp
  std::vector<TH1D> M_Nutp_gen,  Pt_Nutp_gen,   Eta_Nutp_gen,  Phi_Nutp_gen,   Theta_Nutp_gen;
  //TauMu-> Nutm
  std::vector<TH1D> M_Nutm_gen,  Pt_Nutm_gen,   Eta_Nutm_gen,  Phi_Nutm_gen,   Theta_Nutm_gen;
  //TauMu-> Num
  std::vector<TH1D> M_Num_gen,   Pt_Num_gen,    Eta_Num_gen,   Phi_Num_gen,    Theta_Num_gen;
  
  
  // DiTau (Z-> TauPi+TauMu) 
  std::vector<TH1D> M_DiTau_gen, Pt_DiTau_gen,  dEta_DiTau_gen, dPhi_DiTau_gen,  dTheta_DiTau_gen,  dM_DiTau_gen, dR_DiTau_gen;
  // Physics objects(Pi,Mu) 
  std::vector<TH1D> M_MuPi_gen,  Pt_MuPi_gen,   dEta_MuPi_gen,  dPhi_MuPi_gen,   dTheta_MuPi_gen,   dR_MuPi_gen;
  
  //Delta_R(TauPi_rec,TauPi_gen)
  std::vector<TH1D> dR_recTauPi_genTauPi,   dR_recTauMu_genTauMu;
  //Delta_R(Mu_rec,Mu_gen)
  std::vector<TH1D> dR_recMu_genMu,         dR_recPi_genPi;

 private:

};
#endif
