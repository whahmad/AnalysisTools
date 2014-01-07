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

  enum cuts {PrimeVtx,
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
  std::vector<TH2D> muleg_numerator;
  std::vector<TH2D> muleg_numerator_rebin;
  std::vector<TH2D> muleg_denominator;
  std::vector<TH2D> muleg_denominator_rebin;
  std::vector<TH2D> eleg_numerator;
  std::vector<TH2D> eleg_numerator_rebin;
  std::vector<TH2D> eleg_denominator;
  std::vector<TH2D> eleg_denominator_rebin;
  std::vector<TH2D> muleg_eff;
  std::vector<TH2D> eleg_eff;

  //control plots
  std::vector<TH1D> tagmupt;
  std::vector<TH1D> tagmueta;
  std::vector<TH1D> tagept;
  std::vector<TH1D> tageeta;
  std::vector<TH1D> probemupt;
  std::vector<TH1D> probemueta;
  std::vector<TH1D> probeept;
  std::vector<TH1D> probeeeta;
  std::vector<TH1D> drtagmuprobee;
  std::vector<TH1D> drtageprobemu;
  std::vector<TH1D> drmumu;
  std::vector<TH1D> dree;
  std::vector<TH1D> ptbaltagmuprobee;
  std::vector<TH1D> ptbaltageprobemu;
  std::vector<TH1D> ptbalmumu;
  std::vector<TH1D> ptbalee;
  std::vector<TH1D> mttagmu;
  std::vector<TH1D> mttage;
  std::vector<TH1D> mtprobemu;
  std::vector<TH1D> mtprobee;
  std::vector<TH1D> mtagmuprobee;
  std::vector<TH1D> mtageprobemu;
  std::vector<TH1D> mmumu;
  std::vector<TH1D> mee;

  std::vector<TH1D> sip;
  std::vector<TH1D> nmu;
  std::vector<TH1D> ne;

  double mu_pt,mu_eta,e_pt,e_eta;
  int n_mu,n_e;
  
  double cosphi2d(double px1, double py1, double px2, double py2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  
  void triggerMatch(std::vector<unsigned int> tags, std::vector<unsigned int> probes, std::string tagtrigger, std::string probetrigger, std::string tagname, TH2D &denominator, TH2D &numerator,
		  TH1D &tagmupt, TH1D &tagmueta, TH1D &tagept, TH1D & tageeta, TH1D &probemupt, TH1D &probemueta, TH1D &probeept, TH1D &probeeeta,
		  TH1D &drtagmuprobee, TH1D &drtageprobemu, TH1D &ptbaltagmuprobee, TH1D &ptbaltageprobemu, TH1D &mttagmu, TH1D &mttage, TH1D &mtprobemu, TH1D &mtprobee, TH1D &mtagmuprobee, TH1D &mtageprobemu, double w);
  void triggerMatch(std::vector<unsigned int> tags, std::vector<unsigned int> probes, std::string tagtrigger, std::string probetrigger, std::string tagname, TH2D &denominator, TH2D &numerator, double w);

  void doubleMuE(std::vector<unsigned int> objects, std::string particle,std::string tagtrigger, std::string probetrigger, TH2D &denominator, TH2D &numerator,
		  TH1D &tagpt, TH1D &tageta, TH1D &probept, TH1D &probeeta, TH1D &dr, TH1D &ptbal, TH1D &mttag, TH1D &mtprobe, TH1D &m, double w);

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

