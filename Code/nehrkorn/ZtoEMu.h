#ifndef ZtoEMu_h
#define ZtoEMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ZtoEMu : public Selection {

 public:
  ZtoEMu(TString Name_, TString id_);
  virtual ~ZtoEMu();

  virtual void  Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
	     qualitycuts,
	     SameVtx,
	     NMuPt,
	     NMuEta,
		 NMu,
	     NEPt,
	     NEEta,
		 NE,
		 drMuE,
		 diMuonVeto,
		 triLeptonVeto,
		 looseMuonVeto,
		 charge,
		 jetVeto,
		 MtMu,
	     ptBalance,
	     ZMassmax,
	     ZMassmin,
	     Phimin,
	     Phimax,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;

  std::vector<TH1D> RelIsoE;
  std::vector<TH1D> RelIsoMu;
  std::vector<TH1D> EPt;
  std::vector<TH1D> MuPt;
  std::vector<TH1D> mtMu;
  std::vector<TH1D> mtE;
  std::vector<TH1D> pzeta;
  std::vector<TH1D> pzetaDQM;
  std::vector<TH1D> etaMu;
  std::vector<TH1D> etaE;
  std::vector<TH1D> jetsum;
  std::vector<TH1D> NJets;
  std::vector<TH1D> chargesum;
  std::vector<TH1D> chargesumID;
  std::vector<TH1D> drmue;
  std::vector<TH1D> drmueID;
  std::vector<TH1D> deltaphi;
  std::vector<TH1D> deltaphiID;
  std::vector<TH1D> ptbal;
  std::vector<TH1D> chargesumIDsigned;
  std::vector<TH1D> chargesumsigned;
  std::vector<TH1D> FirstJetPt;
  std::vector<TH1D> SecondJetPt;
  std::vector<TH1D> ThirdJetPt;
  std::vector<TH1D> FourthJetPt;
  
  std::vector<TH1D> invmass;
  std::vector<TH1D> invmass_zmass;
  std::vector<TH1D> invmass_ptbalance;
  std::vector<TH1D> invmass_mtmu;
  std::vector<TH1D> invmass_jetveto;
  std::vector<TH1D> invmass_charge;
  std::vector<TH1D> invmass_loosemuonveto;
  std::vector<TH1D> invmass_dremu;
  std::vector<TH1D> invmass_only_object_id;
  
  std::vector<TH1D> nm2_charge;
  std::vector<TH1D> nm2_jetveto;
  std::vector<TH1D> nm2_mtmu;
  std::vector<TH1D> nm2_ptbalance;
  std::vector<TH1D> nm2_drmue;

  std::vector<TH1D> no_jetveto_mtmu;
  std::vector<TH1D> no_charge_ptbalance;
  std::vector<TH1D> no_charge;
  
  std::vector<TH1D> drmue_plus_charge;
  std::vector<TH1D> drmue_plus_jetveto;
  std::vector<TH1D> drmue_plus_mtmu;
  std::vector<TH1D> drmue_plus_ptbalance;
  
  std::vector<TH1D> mproj1;
  std::vector<TH1D> mproj2;
  std::vector<TH1D> mproj3;
  std::vector<TH1D> mproj4;
  std::vector<TH1D> mproj5;
  std::vector<TH1D> mproj6;
  std::vector<TH1D> mproj7;
  std::vector<TH1D> mproj8;
  std::vector<TH1D> mproj9;
  std::vector<TH1D> phiproj1;
  std::vector<TH1D> phiproj2;
  std::vector<TH1D> phiproj3;
  std::vector<TH1D> phiproj4;
  std::vector<TH1D> phiproj5;
  std::vector<TH1D> phiproj6;
  std::vector<TH1D> phiproj7;
  std::vector<TH1D> phiproj8;
  std::vector<TH1D> phiproj9;
  
  std::vector<TH1D> m1;
  std::vector<TH1D> m2;
  std::vector<TH1D> m3;
  std::vector<TH1D> phi1;
  std::vector<TH1D> phi2;
  std::vector<TH1D> phi3;
  
  std::vector<TH1D> phi1_nopt;
  std::vector<TH1D> phi2_nopt;
  std::vector<TH1D> phi3_nopt;
  
  std::vector<TH1D> phi1_ptdiff;
  std::vector<TH1D> phi2_ptdiff;
  std::vector<TH1D> phi3_ptdiff;
  
  std::vector<TH1D> ptbal2;
  
  std::vector<TH1D> frMu;
  std::vector<TH1D> frE;
  
  std::vector<TH1D> pzetaCut;
  std::vector<TH1D> pzetaStatus;
  std::vector<TH1D> metCut;
  std::vector<TH1D> metStatus;
  
  std::vector<TH2D> InvmassVsDeltaPhi;
  std::vector<TH2D> PtDiffVsDeltaPhi;

  double mu_pt,mu_eta,e_pt,e_eta,jet_pt,jet_eta,jet_sum,zmin,zmax,phimin,phimax,mtmu,ptbalance,dRmue;
  int n_mu,n_e;
  double pex,pey,pmux,pmuy,phie,phimu;
  double combpt;
  double aemu; //angle between electron and muon
  double beta; //angle between combined pt and bisector of electron and muon
  double gamma; //angle between MET and bisector of electron and muon
  double phismall; //smaller angle (electron or muon)
  double pvis,pmiss;
  
  double calculatePzeta(int muiterator, int eiterator,std::vector<unsigned int> vec1, std::vector<unsigned int> vec2);
  double calculatePzetaDQM(int muiterator, int eiterator);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  
  bool isTightMuon(unsigned int i);
  bool isTightMuon(unsigned int i, unsigned int j);
  bool isLooseMuon(unsigned int i);
  bool isFakeMuon(unsigned int i);
  bool isFakeMuon(unsigned int i, unsigned int j);
  double Muon_RelIso(unsigned int i);
  double Muon_AbsIso(unsigned int i);
  
  bool isMVAElectron(unsigned int i);
  bool isTightElectron(unsigned int i);
  bool isTightElectron(unsigned int i, unsigned int j);
  bool isLooseElectron(unsigned int i);
  bool isFakeElectron(unsigned int i);
  bool isFakeElectron(unsigned int i, unsigned int j);
  double Electron_RelIso(unsigned int i);
  double Electron_Aeff(double Eta);
  
  double MuonSF(unsigned int i);
  double MuonDataSF(unsigned int i);
  double ElectronSF(unsigned int i);
  double ElectronDataSF(unsigned int i);
  double ElectronEffRecHit(unsigned int i);
  
  double Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type);
  
  bool MVA_ID;
  TFile* FRFile;
  TH2D* ElectronFakeRate;
  TH2D* MuonFakeRate;
  double fakeRate;
  double fakeRateMu;
  double fakeRateE;
  
  bool twod;

};
#endif

