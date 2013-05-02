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
  //virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
	     qualitycuts,
	     NMuPt,
	     NMuEta,
		 NMu,
	     NEPt,
	     NEEta,
		 NE,
		 drMuE,
		 charge,
		 MuRelIso,
		 looseMuonVeto,
		 SameVtx,
		 jetVeto,
		 MtMu,
	     ptBalance,
	     ZMassmax,
	     ZMassmin,
	     MET,
	     MtE,
         deltaPhi,
         cosdeltaPhi,
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
  std::vector<TH1D> JetPt;
  std::vector<TH1D> JetFromVtxPt;
  std::vector<TH1D> drejet;
  std::vector<TH1D> invmass;
  std::vector<TH1D> etaMu;
  std::vector<TH1D> etaE;
  std::vector<TH1D> jetsum;
  std::vector<TH1D> chargesum;
  std::vector<TH1D> chargesumID;
  std::vector<TH1D> chargesumMuID;
  std::vector<TH1D> chargesumEID;
  std::vector<TH1D> chargeE;
  std::vector<TH1D> chargeEID;
  std::vector<TH1D> chargeMu;
  std::vector<TH1D> chargeMuID;
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
  std::vector<TH1D> MuPoca;
  std::vector<TH1D> EPoca;
  std::vector<TH1D> EMuPoca;
  std::vector<TH1D> EMuPocaDiff;
  std::vector<TH1D> NMuonsInRange;
  std::vector<TH1D> NElectronsInRange;
  std::vector<TH1D> phiMu;
  std::vector<TH1D> MuDXY;
  std::vector<TH1D> MuDZ;
  std::vector<TH1D> EDXY;
  std::vector<TH1D> EDZ;
  std::vector<TH1D> mudxypass;
  std::vector<TH1D> mudzpass;
  std::vector<TH1D> edxypass;
  std::vector<TH1D> edzpass;
  std::vector<TH1D> etachargemu;
  std::vector<TH1D> etachargee;
  std::vector<TH1D> jetsumpass;
  
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
  
  std::vector<TH2D> ABCD;
  std::vector<TH2D> ABCDEF;
  std::vector<TH2D> chargemue;
  std::vector<TH2D> InvmassVsDeltaPhi;
  std::vector<TH2D> EPtVsCosTheta;
  std::vector<TH2D> MuPtVsCosTheta;
  std::vector<TH2D> EPtVsEta;
  std::vector<TH2D> MuPtVsEta;

  double mu_pt,mu_eta,mu_reliso,e_pt,e_eta,jet_pt,jet_eta,jet_sum;
  int n_mu,n_e;
  double pex,pey,pmux,pmuy,phie,phimu;
  double combpt;
  double aemu; //angle between electron and muon
  double beta; //angle between combined pt and bisector of electron and muon
  double gamma; //angle between MET and bisector of electron and muon
  double phi1; //smaller angle (electron or muon)
  double pvis,pmiss;
  
  double calculatePzeta(int muiterator, int eiterator,std::vector<unsigned int> vec1, std::vector<unsigned int> vec2);
  double calculatePzetaDQM(int muiterator, int eiterator,std::vector<unsigned int> vec1, std::vector<unsigned int> vec2);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  double jetdxy(int vtx_idx, int leadingtrack_idx, int jet_idx);
  double jetdz(int vtx_idx, int leadingtrack_idx, int jet_idx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);

};
#endif

