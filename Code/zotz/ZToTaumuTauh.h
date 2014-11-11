#ifndef ZToTaumuTauh_h
#define ZToTaumuTauh_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "ReferenceScaleFactors.h"

class ZToTaumuTauh : public Selection {

 public:
  ZToTaumuTauh(TString Name_, TString id_);
  virtual ~ZToTaumuTauh();

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
	  NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  ReferenceScaleFactors *RSF;

  bool verbose;

  // cut values
  double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta, cMu_dRHltMatch;
  double cTau_pt, cTau_eta, cMuTau_dR, cTau_dRHltMatch, cTau_IsoRaw;
  std::vector<TString> cTriggerNames;

  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> NMtTauMET;
  std::vector<TH1D> NMvis;
  std::vector<TH1D> NSignal_SB_WJets;
  std::vector<TH1D> NSB_Data;
  std::vector<TH1D> NSB;
  std::vector<TH1D> Mu_pt, Mu_eta, Mu_phi, Tau_pt, Tau_eta, Tau_phi;
  std::vector<TH1D> NQCD_MT_MuMET_A, NQCD_MT_MuMET_B, NQCD_MT_MuMET_C, NQCD_MT_MuMET_D;
  std::vector<TH1D> NQCD_MET_A, NQCD_MET_B, NQCD_MET_C, NQCD_MET_D;

  int Charge;

  double Signal_upperLimit;
  double SB_lowerLimit, SB_upperLimit;

  bool Scaleby_Counting;

  bool selectMuon_Id(unsigned i, unsigned vertex);
  bool selectMuon_Kinematics(unsigned i);
  bool selectMuon_Isolation(unsigned i);
  bool selectMuon_AntiIsolation(unsigned i);

  bool selectPFTau_Id(unsigned i);
  bool selectPFTau_Id(unsigned i, std::vector<int> muonCollection);
  bool selectPFTau_Isolation(unsigned i);
  bool selectPFTau_Kinematics(unsigned i);

 private:

};
#endif
