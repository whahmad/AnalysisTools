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
	  MT,
	  CatCut1,
	  CatCut2,
	  CatCut3,
	  CatCut4,
	  NCuts};

  	// cuts in categories
	enum cuts_VBF {
		VbfNJet	= CatCut1,
		VbfDeltaEta,
		VbfNJetRapGap,
		VbfJetInvM,
		VbfNCuts
	};
	enum cuts_OneJetHigh {
		OneJetNJet = CatCut1,
		OneJetNoVBF,
		OneJetNBtagJets,
		OneJetTauPt,
		OneJetNCuts
	};
	enum cuts_ZeroJetHigh {
		ZeroJetNJet = CatCut1,
		ZeroJetNBtagJets,
		ZeroJetTauPt,
		ZeroJetNCuts
	};
	enum cuts_NoCategory {
		NoCategoryNCuts = CatCut1
	};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  virtual void  Finish();

 private:
  // Selection Variables
  std::vector<TH1D> NCatFired;

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
  std::vector<TH1D> MuPhi;

  std::vector<TH1D> MuSelPt;
  std::vector<TH1D> MuSelEta;
  std::vector<TH1D> MuSelPhi;
  std::vector<TH1D> MuSelFakesTauID;

  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPhi;

  std::vector<TH1D> TauSelPt;
  std::vector<TH1D> TauSelEta;
  std::vector<TH1D> TauSelPhi;

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

  std::vector<TH1D> MetPt;
  std::vector<TH1D> MetPhi;


  // cut values
  double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta;
  double cTau_pt, cTau_eta, cMuTau_dR;
  double cMuTriLep_pt, cMuTriLep_eta, cEleTriLep_pt, cEleTriLep_eta;
  std::vector<TString> cTriggerNames;
  double cVBFJet_eta, cVBFJet_pt;
  double cCat_jetPt, cCat_jetEta, cCat_btagPt, cCat_splitTauPt;

  // flag for category to run
  TString categoryFlag;

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
  bool selectElectron_triLeptonVeto(unsigned i, unsigned i_vtx, std::vector<int>);

  bool selectPFTau_Id(unsigned i);
  bool selectPFTau_Id(unsigned i, std::vector<int>);
  bool selectPFTau_Iso(unsigned i);
  bool selectPFTau_Kinematics(unsigned i);

  bool selectPFJet_VBF(unsigned i);
  bool selectPFJet_Categories(unsigned i);
  bool selectBJet_Categories(unsigned i);

  inline double transverseMass(double pt1, double phi1, double pt2, double phi2){
	  return sqrt(2 * pt1 * pt2 * (1 - cos(phi1 - phi2)));
  }

  // categories
  std::vector<float> cut_VBF, cut_OneJet, cut_ZeroJet, cut_NoCategory;

  void configure_VBF();
  bool category_VBF();

  void configure_OneJetHigh();
  bool category_OneJetHigh(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection, bool passedVBF);

  void configure_OneJetLow();
  bool category_OneJetLow(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection, bool passedVBF);

  void configure_ZeroJetHigh();
  bool category_ZeroJetHigh(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection);

  void configure_ZeroJetLow();
  bool category_ZeroJetLow(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection);

  void configure_NoCategory();
  bool category_NoCategory();

};
#endif
