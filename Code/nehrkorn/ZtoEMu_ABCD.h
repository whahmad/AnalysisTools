#ifndef ZtoEMu_ABCD_h
#define ZtoEMu_ABCD_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TGraphAsymmErrors.h"

class ZtoEMu_ABCD : public Selection {

 public:
  ZtoEMu_ABCD(TString Name_, TString id_);
  virtual ~ZtoEMu_ABCD();

  virtual void  Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
		 NMu,
		 NE,
		 ptthreshold,
		 diMuonVeto,
 		 triLeptonVeto,
 		 charge,
 		 mass,
 		 jetVeto,
		 mtmu,
		 ptbal,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  // wjets
  std::vector<TH1D> wjets_a;
  std::vector<TH1D> wjets_b;
  std::vector<TH1D> wjets_c;
  std::vector<TH1D> wjets_d;
  std::vector<TH1D> wa;
  std::vector<TH1D> wb;
  std::vector<TH1D> wc;
  std::vector<TH1D> wd;

  std::vector<TH1D> wjets_control_a;
  std::vector<TH1D> wjets_control_b;
  std::vector<TH1D> wjets_control_c;
  std::vector<TH1D> wjets_control_d;

  std::vector<TH2D> wjets;

  double mu_ptlow,mu_pthigh,mu_eta,e_ptlow,e_pthigh,e_eta,jet_pt,jet_eta,jet_sum,zmin,zmax;
  int n_mu,n_e;

  double cosphi2d(double px1, double py1, double px2, double py2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  bool matchTrigger(unsigned int i, double dr, std::string trigger, std::string object);
  int matchTruth(TLorentzVector tvector);
  bool matchTruth(TLorentzVector tvector, int pid, double dr);

  bool isTightMuon(unsigned int i);
  bool isTightMuon(unsigned int i, unsigned int j);
  bool isHiggsMuon(unsigned int i, unsigned int j);
  bool isLooseMuon(unsigned int i);
  bool isFakeMuon(unsigned int i);
  bool isFakeMuon(unsigned int i, unsigned int j);
  double Muon_RelIso(unsigned int i);
  double Muon_AbsIso(unsigned int i);

  bool isTrigPreselElectron(unsigned int i);
  bool isTrigNoIPPreselElectron(unsigned int i);
  bool isMVATrigElectron(unsigned int i);
  bool isMVATrigNoIPElectron(unsigned int i);
  bool isMVANonTrigElectron(unsigned int i, unsigned int j);
  bool isHiggsElectron(unsigned int i, unsigned int j);
  bool isTightElectron(unsigned int i);
  bool isTightElectron(unsigned int i, unsigned int j);
  bool isFakeElectron(unsigned int i);
  bool isFakeElectron(unsigned int i, unsigned int j);
  double Electron_RelIso(unsigned int i);
  double Electron_Aeff_R04(double Eta);
  double Electron_Aeff_R03(double Eta);

  double MuonIDeff(unsigned int i);
  double MuonHiggsIDeff(unsigned int i);
  double ElectronIDeff(unsigned int i, std::string id);
  double ElectronTrigIDeff(unsigned int i);
  double ElectronNonTrigIDeff(unsigned int i);
  double ElectronHiggsIDeff(unsigned int i);

  double Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type);
  double rundependentJetPtCorrection(double jeteta, int runnumber);

  TFile* FRFile;
  TFile* EmbEffFile;
  TH2D* ElectronFakeRate;
  TH2D* MuonFakeRate;
  TH2D* EmbEff;

  TFile* MuIdEffFile;
  TFile* MuIsoEffFile;
  TFile* ETrigIdEffFile;
  TFile* ENonTrigIdEffFile;

  TH2D* ElectronTrigEff;
  TH2D* ElectronNonTrigEff;
  TGraphAsymmErrors* MuIdEff09;
  TGraphAsymmErrors* MuIdEff12;
  TGraphAsymmErrors* MuIdEff21;
  TGraphAsymmErrors* MuIdEff24;
  TGraphAsymmErrors* MuIsoEff09;
  TGraphAsymmErrors* MuIsoEff12;
  TGraphAsymmErrors* MuIsoEff21;
  TGraphAsymmErrors* MuIsoEff24;

};
#endif

