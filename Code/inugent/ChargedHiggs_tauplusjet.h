#ifndef ChargedHiggs_h
#define ChargedHiggs_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class ChargedHiggs : public Selection {

 public:
  ChargedHiggs(TString Name_, TString id_);
  virtual ~ChargedHiggs();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     PrimeVtx,
	     NTauKinFit,
	     NTauPt,
	     NTauEta,
             N1Jets,
	     N2Jets,
	     N3Jets,
	     NJets,
             NBJets,
             MET,
             HT,
	     etaq,
	     HadWMass,
	     HadTopMass,
	     TauMETTopMT,
	     TauMETdphi,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH2D> TagEtaPT;
  std::vector<TH1D> ChargedHiggsMT;
  double tau_pt,tau_eta,jet_pt,jet_eta;

};
#endif
