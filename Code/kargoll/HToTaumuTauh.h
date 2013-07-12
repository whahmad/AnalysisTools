#ifndef HToTaumuTauh_h
#define HToTaumuTauh_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class HToTaumuTauh : public Selection {

 public:
  HToTaumuTauh(TString Name_, TString id_);
  virtual ~HToTaumuTauh();

  virtual void  Configure();

  enum cuts {
	  TriggerOk=0,
	  PrimeVtx,
	  NMuId,
	  NMuKin,
	  DiMuonVeto,
	  NTauId,
	  NTauIso,
	  NTauKin,
	  OppCharge,
	  TriLeptonVeto,
	  NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  virtual void  Finish();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NVtxFullSelection;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> VtxZ;
  std::vector<TH1D> VtxRho;
  std::vector<TH1D> VtxPhi;
  std::vector<TH1D> VtxNdof;
  std::vector<TH1D> VtxIsfake;

  std::vector<TH1D> MuDxy;
  std::vector<TH1D> MuDz;
  std::vector<TH1D> MuRelIso;
  std::vector<TH1D> MuPt;
  std::vector<TH1D> MuEta;

  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauEta;

  std::vector<TH1D> MuVetoDPtSelMuon;
  std::vector<TH1D> MuVetoInvM;
  std::vector<TH1D> MuVetoPtPositive;
  std::vector<TH1D> MuVetoPtNegative;
  std::vector<TH1D> MuVetoDRTau;

  std::vector<TH1D> NMuonTriLepVeto;
  std::vector<TH1D> NElecTriLepVeto;

  std::vector<TH1D> MuCharge;
  std::vector<TH1D> TauCharge;

  std::vector<TH1D> MuTauDR;
  std::vector<TH1D> MuTauDPhi;
  std::vector<TH1D> MuTauDEta;
  std::vector<TH1D> MuTauDPt;
  std::vector<TH1D> MuTauRelDPt;


  // cut values
  double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta;
  double cTau_pt, cTau_eta;

  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);

  bool selectVertex(unsigned i);

  double Muon_AbsIso(unsigned int i);
  double Muon_RelIso(unsigned int i);
  bool isTightMuon(unsigned i);
  bool isTightMuon(unsigned i, unsigned i_vtx);
  bool selectMuon_Id(unsigned i, unsigned vertex);
  bool selectMuon_Kinematics(unsigned i);

  bool selectMuon_diMuonVeto(unsigned i, unsigned i_vtx);
  bool selectMuon_triLeptonVeto(unsigned i, int selectedMuon, unsigned i_vtx);

  bool isLooseMVAElectron(unsigned i);
  double Electron_RelIso(unsigned i);
  bool selectElectron_triLeptonVeto(unsigned i, unsigned i_vtx);

  bool selectPFTau_Id(unsigned i, unsigned i_muon);
  bool selectPFTau_Iso(unsigned i);
  bool selectPFTau_Kinematics(unsigned i);

};
#endif
