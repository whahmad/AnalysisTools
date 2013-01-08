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

  enum cuts {//TriggerOk=0, 
	     PrimeVtx,
	     NMuPt,
	     NMuEta,
		 NMu,
	     MuIso,
	     NEPt,
	     NEEta,
		 NE,
	     EIso,
		 jetVeto,
         charge,
         deltaPhi,
	     cosdeltaPhi,
	     ptBalance,
         MET,
         ZMassmax,
         ZMassmin,
		 Mt,
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
  std::vector<TH1D> phiMtMu;
  std::vector<TH1D> mtE;
  std::vector<TH1D> phiMtE;
  std::vector<TH1D> NJets;
  std::vector<TH1D> pzeta;
  std::vector<TH1D> pzetaDQM;
  std::vector<TH2D> METvsMt;
    
  double mu_pt,mu_eta,e_pt,e_eta,jet_pt,jet_eta;
  int n_mu,n_e;
  double pex,pey,pmux,pmuy,phie,phimu;
  double combpt;
  double aemu; //angle between electron and muon
  double beta; //angle between combined pt and bisector of electron and muon
  double gamma; //angle between MET and bisector of electron and muon
  double phi1; //smaller angle (electron or muon)
  double pvis,pmiss;
  
  double calculatePzeta(int iterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2);
  double calculatePzetaDQM(int iterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2);

};
#endif

