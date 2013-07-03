#ifndef ZtoEMu_Skim_h
#define ZtoEMu_Skim_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ZtoEMu_Skim : public Selection {

 public:
  ZtoEMu_Skim(TString Name_, TString id_);
  virtual ~ZtoEMu_Skim();

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
		 charge,
		 MtMu,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  double mu_pt,mu_eta,e_pt,e_eta,jet_pt,jet_eta,jet_sum,zmin,zmax;
  int n_mu,n_e;
  
  double cosphi2d(double px1, double py1, double px2, double py2);
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

};
#endif

