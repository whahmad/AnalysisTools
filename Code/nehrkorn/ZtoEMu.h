#ifndef ZtoEMu_h
#define ZtoEMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

class ZtoEMu : public Selection {

 public:
  ZtoEMu(TString Name_, TString id_);
  virtual ~ZtoEMu();

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
		 jetVeto,
		 MtMu,
	     ptBalance,
	     ZMassmax,
	     ZMassmin,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> RelIsoE;
  std::vector<TH1D> RelIsoMu;
  std::vector<TH1D> EPt;
  std::vector<TH1D> MuPt;
  std::vector<TH1D> mtMu;
  std::vector<TH1D> mtE;
  std::vector<TH1D> etaMu;
  std::vector<TH1D> etaE;
  std::vector<TH1D> jetsum;
  std::vector<TH1D> NJets;
  std::vector<TH1D> NJetsLoose;
  std::vector<TH1D> NJetsMedium;
  std::vector<TH1D> NJetsTight;
  std::vector<TH1D> chargesum;
  std::vector<TH1D> drmue;
  std::vector<TH1D> deltaphi;
  std::vector<TH1D> ptbal;
  std::vector<TH1D> chargesumsigned;
  std::vector<TH1D> FirstJetPt;
  std::vector<TH1D> SecondJetPt;
  std::vector<TH1D> jeteta;
  
  std::vector<TH1D> invmass_zmass;
  std::vector<TH1D> invmass_ptbalance;
  std::vector<TH1D> invmass_mtmu;
  std::vector<TH1D> invmass_jetveto;
  std::vector<TH1D> invmass_vetos;
  std::vector<TH1D> invmass_only_object_id;
  
  std::vector<TH1D> nm0_met;
  std::vector<TH1D> nm0_jetsum;
  std::vector<TH1D> nm0_onejet;
  std::vector<TH1D> nm0_mtmu;
  std::vector<TH1D> nm0_ptbalance;
  
  std::vector<TH1D> NPV;
  
  std::vector<TH1D> met;
  std::vector<TH1D> met_xycorr;
  std::vector<TH1D> met_uncorr;
  std::vector<TH1D> onejet;
  std::vector<TH1D> mte_mtmu;
  std::vector<TH1D> leadingjet;
  std::vector<TH1D> subleadingjet;
  std::vector<TH1D> sumjets;

  // comparison of generators

  std::vector<TH1D> zpt;
  std::vector<TH1D> zeta;
  std::vector<TH1D> zmass;
  std::vector<TH1D> leadingjet_pt;
  std::vector<TH1D> subleadingjet_pt;
  std::vector<TH1D> leadingjet_eta;
  std::vector<TH1D> subleadingjet_eta;
  std::vector<TH1D> jetsumcustom;

  double mu_ptlow,mu_pthigh,mu_eta,e_ptlow,e_pthigh,e_eta,jet_pt,jet_eta,jet_sum,zmin,zmax,mtmu,ptbalance;
  int n_mu,n_e;
  bool doHiggsObjects;
  bool runOverSkim;
  
  double calculatePzeta(int muiterator, int eiterator);
  double calculatePzetaDQM(int muiterator, int eiterator);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  bool matchTrigger(unsigned int i, double dr, std::string trigger, std::string object);
  int matchTruth(TLorentzVector tvector);
  bool matchTruth(TLorentzVector tvector, int pid, double dr);
  int findBin(TGraphAsymmErrors* graph, double xval);
  
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
  bool isLooseElectron(unsigned int i);
  bool isFakeElectron(unsigned int i);
  bool isFakeElectron(unsigned int i, unsigned int j);
  double Electron_RelIso(unsigned int i);
  double Electron_Aeff_R04(double Eta);
  double Electron_Aeff_R03(double Eta);
  
  double MuonIDeff(unsigned int i);
  double MuonIDerrUp(unsigned int i);
  double MuonIDerrDown(unsigned int i);
  double MuonHiggsIDeff(unsigned int i);
  double MuonTriggerEff(unsigned int i);
  double MuonTriggerErr(unsigned int i);
  double ElectronIDeff(unsigned int i, std::string id);
  double ElectronIDerr(unsigned int i, std::string id);
  double ElectronTrigIDeff(unsigned int i);
  double ElectronTrigIDerr(unsigned int i);
  double ElectronNonTrigIDeff(unsigned int i);
  double ElectronNonTrigIDerr(unsigned int i);
  double ElectronHiggsIDeff(unsigned int i);
  double ElectronTriggerEff(unsigned int i);
  double ElectronTriggerErr(unsigned int i);
  
  double ElectronMassScale(unsigned int i);
  double ZPtReweight(double zpt);
  double rundependentJetPtCorrection(double jeteta, int runnumber);

  //double JECuncertainty(unsigned int i, TString datamc);

  double Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type);
  
  TFile* FRFile;
  TFile* EmbEffFile;
  TH2D* ElectronFakeRate;
  TH2D* MuonFakeRate;
  TH2D* EmbEff;
  double fakeRate;
  double fakeRateMu;
  double fakeRateE;
  
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

  TF1* gause;
  TF1* gausmu;
  TH1D* eres;
  TH1D* mures;
  double eleres;
  double muonres;

};
#endif

