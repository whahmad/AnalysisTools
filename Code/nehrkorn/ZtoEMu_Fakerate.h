#ifndef ZtoEMu_Fakerate_h
#define ZtoEMu_Fakerate_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ZtoEMu_Fakerate : public Selection {

 public:
  ZtoEMu_Fakerate(TString Name_, TString id_);
  virtual ~ZtoEMu_Fakerate();

  virtual void  Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	  	 PrimeVtx,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  
  // fakerate
  std::vector<TH2D> tightmu;
  std::vector<TH2D> tightmu_rebin;
  std::vector<TH2D> tighte;
  std::vector<TH2D> tighte_rebin;
  std::vector<TH2D> fakemu;
  std::vector<TH2D> fakemu_rebin;
  std::vector<TH2D> fakee;
  std::vector<TH2D> fakee_rebin;
  std::vector<TH2D> mueff;
  std::vector<TH2D> eeff;

  // trigger eff
  std::vector<TH1D> mudr;
  std::vector<TH1D> mupt;
  std::vector<TH1D> edr;
  std::vector<TH1D> ept;
  
  std::vector<TH2D> muleg_numerator;
  std::vector<TH2D> muleg_denominator;
  std::vector<TH2D> eleg_numerator;
  std::vector<TH2D> eleg_denominator;

  double mu_pt,mu_eta,e_pt,e_eta;
  int n_mu,n_e;
  
  double cosphi2d(double px1, double py1, double px2, double py2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  
  bool isTightMuon(unsigned int i);
  bool isTightMuon(unsigned int i, unsigned int j);
  bool isLooseMuon(unsigned int i);
  bool isFakeMuon(unsigned int i);
  bool isFakeMuon(unsigned int i, unsigned int j);
  double Muon_RelIso(unsigned int i);
  double Muon_AbsIso(unsigned int i);
  
  bool isMVATrigElectron(unsigned int i);
  bool isMVATrigNoIPElectron(unsigned int i);
  bool isMVANonTrigElectron(unsigned int i, unsigned int j);
  bool isTightElectron(unsigned int i);
  bool isTightElectron(unsigned int i, unsigned int j);
  bool isFakeElectron(unsigned int i);
  bool isFakeElectron(unsigned int i, unsigned int j);
  double Electron_RelIso(unsigned int i);
  double Electron_Aeff_R04(double Eta);
  double Electron_Aeff_R03(double Eta);
  
  double Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type);
  int getxbin(double pt);
  int getybin(double eta, std::string object);
  
  TFile* outfile;

};
#endif

