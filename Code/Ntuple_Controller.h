//Ntuple_Controller.h HEADER FILE

#ifndef Ntuple_Controller_h
#define Ntuple_Controller_h


// Root include files
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixT.h"

// Include files (C & C++ libraries)
#include<iostream>
#include <vector>
#include <string.h>

#include "NtupleReader.h"

#include "HistoConfig.h"
///////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//*
//*   Class: Ntuple_Controller
//*   
//*   Purpose: The purpose of this class is to provide a interface to the
//*            Ntuple
//*
//*   Designed by: Ian Nugent
//*
//*
//*****************************************************************************
///////////////////////////////////////////////////////////////////////////////


class Ntuple_Controller{
 private:
  NtupleReader *Ntp;
  TFile *newfile;
  TTree *SkimmedTree;
  int nbytes;
  int jentry;
  int nb;
  bool copyTree;

  // Ntuple Access Functions
  virtual void Branch_Setup(TString B_Name, int type);
  virtual void Branch_Setup(){}

  // Functions to configure objects
  virtual void ConfigureObjects(); 
  void doElectrons();
  void doPhotons();
  void doJets();
  void doMuons();
  void doTaus();
  void doMET();
  unsigned int ObjEvent;

  // Object Variables
  std::vector<TLorentzVector> electrons_default;
  std::vector<TLorentzVector> photons_default;
  std::vector<TLorentzVector> jets_default;
  std::vector<TLorentzVector> muons_default;
  std::vector<TLorentzVector> taus_default;
  TLorentzVector              met_default;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> muons;
  std::vector<TLorentzVector> taus;
  TLorentzVector              met;

  // Systematic controls variables
  int theSys;
  HistoConfig HConfig;
 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles);

  // Destructor
  ~Ntuple_Controller() ;

  enum TrackQuality {
    undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
    confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
    qualitySize = 7
  };
  enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};

  // Ntuple Access Functions 
  virtual Int_t Get_Entries();
  virtual void Get_Event(int _jentry);
  virtual Int_t Get_EventIndex();
  virtual TString Get_File_Name();

  //Ntuple Cloning Functions
  virtual void CloneTree(TString n);
  virtual void SaveCloneTree();
  inline void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}

  // Systematic controls
  enum    Systematic {Default=0,NSystematics};
  enum    MCDataID {Data=0,MC_DY=1,MC_WJ=2,MC_QCD=3,DY_Signal=10230530,Signal=998,Unknown=999};
  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}
  double  GetSystematicWeight();


  // Data/MC switch and thin
  bool isData(){return Ntp->Event_isRealData;}
  void ThinTree();

  // Physics Varible Get Functions
  // Event Varibles
  int GetMCID();
  unsigned int RunNumber(){return Ntp->Event_RunNumber;;}
  unsigned int EventNumber(){ return Ntp->Event_EventNumber;}
  int BunchCrossing(){ return Ntp->Event_bunchCrossing;}
  int OrbitNumber(){ return Ntp->Event_orbitNumber;}
  unsigned int LuminosityBlock(){return Ntp->Event_luminosityBlock;}
  int           PileupInfo_NumInteractions_nm1(){Ntp->PileupInfo_NumInteractions_nm1;}
  int           PileupInfo_NumInteractions_n0(){Ntp->PileupInfo_NumInteractions_n0;}
  int           PileupInfo_NumInteractions_np1(){Ntp->PileupInfo_NumInteractions_np1;}
  double        EvtWeight3D(){return Ntp->EvtWeight3D;}



  // Vertex Information
  unsigned int NVtx(){return Ntp->Vtx_ndof->size();}
  TVector3     Vtx(unsigned int i){return TVector3(Ntp->Vtx_x->at(i),Ntp->Vtx_y->at(i),Ntp->Vtx_z->at(i));}
  float        Vtx_chi2(unsigned int i){return Ntp->Vtx_chi2->at(i);}
  float        Vtx_nTrk(unsigned int i){return Ntp->Vtx_nTrk->at(i);}
  float        Vtx_ndof(unsigned int i){return Ntp->Vtx_ndof->at(i);}
  TMatrixF     Vtx_Cov(unsigned int i);
  std::vector<int>  Vtx_Track_idx(unsigned int i){return Ntp->Vtx_Track_idx->at(i);}
  float Vtx_isFake(unsigned int i){return Ntp->Vtx_isFake->at(i);}
  bool isVtxGood(unsigned int i); 


  // Muon information
  unsigned int   NMuons(){return Ntp->Muon_p4->size();}
  TLorentzVector Muons_p4(unsigned int i){return TLorentzVector(Ntp->Muon_p4->at(i).at(1),Ntp->Muon_p4->at(i).at(2),Ntp->Muon_p4->at(i).at(3),Ntp->Muon_p4->at(i).at(0));}
  TVector3       Muon_Poca(unsigned int i){return TVector3(Ntp->Muon_Poca->at(i).at(0),Ntp->Muon_Poca->at(i).at(1),Ntp->Muon_Poca->at(i).at(2));}
  bool           Muon_isGlobalMuon(unsigned int i){return Ntp->Muon_isGlobalMuon->at(i);}
  bool           Muon_isStandAloneMuon(unsigned int i){return Ntp->Muon_isStandAloneMuon->at(i);}
  bool           Muon_isTrackerMuon(unsigned int i){return Ntp->Muon_isTrackerMuon->at(i);}
  bool           Muon_isCaloMuon(unsigned int i){return Ntp->Muon_isCaloMuon->at(i);}
  bool           Muon_isIsolationValid(unsigned int i){return Ntp->Muon_isIsolationValid->at(i);}
  bool           Muon_isQualityValid(unsigned int i){return Ntp->Muon_isQualityValid->at(i);}
  bool           Muon_isTimeValid(unsigned int i){return Ntp->Muon_isTimeValid->at(i);}
  float          Muon_emEt03(unsigned int i){return Ntp->Muon_emEt03->at(i);}
  float          Muon_emVetoEt03(unsigned int i){return Ntp->Muon_emVetoEt03->at(i);}
  float          Muon_hadEt03(unsigned int i){return Ntp->Muon_hadEt03->at(i);}
  float          Muon_hadVetoEt03(unsigned int i){return Ntp->Muon_hadVetoEt03->at(i);}
  float          Muon_nJets03(unsigned int i){return Ntp->Muon_nJets03->at(i);}
  float          Muon_nTracks03(unsigned int i){return Ntp->Muon_nTracks03->at(i);}
  float          Muon_sumPt03(unsigned int i){return Ntp->Muon_sumPt03->at(i);}
  float          Muon_trackerVetoPt03(unsigned int i){return Ntp->Muon_trackerVetoPt03->at(i);}
  float          Muon_emEt05(unsigned int i){return Ntp->Muon_emEt05->at(i);}
  float          Muon_emVetoEt05(unsigned int i){return Ntp->Muon_emVetoEt05->at(i);}
  float          Muon_hadEt05(unsigned int i){return Ntp->Muon_hadEt05->at(i);}
  float          Muon_hadVetoEt05(unsigned int i){return Ntp->Muon_hadVetoEt05->at(i);}
  float          Muon_nJets05(unsigned int i){return Ntp->Muon_nJets05->at(i);}
  float          Muon_nTracks05(unsigned int i){return Ntp->Muon_nTracks05->at(i);}
  float          Muon_sumPt05(unsigned int i){return Ntp->Muon_sumPt05->at(i);}
  float          Muon_trackerVetoPt05(unsigned int i){return Ntp->Muon_trackerVetoPt05->at(i);}
  unsigned int   Muon_Track_idx(unsigned int i){return Ntp->Muon_Track_idx->at(i);}
  float          Muon_hitPattern_pixelLayerwithMeas(unsigned int i){return Ntp->Muon_hitPattern_pixelLayerwithMeas->at(i);}
  float          Muon_numberOfMatchedStations(unsigned int i){return Ntp->Muon_numberOfMatchedStations->at(i);}
  float          Muon_normChi2(unsigned int i){return Ntp->Muon_normChi2->at(i);}
  float          Muon_hitPattern_numberOfValidMuonHits(unsigned int i){return Ntp->Muon_hitPattern_numberOfValidMuonHits->at(i);}
  float          Muon_innerTrack_numberofValidHits(unsigned int i){return Ntp->Muon_innerTrack_numberofValidHits->at(i);}
  float          Muon_numberOfMatches(unsigned int i){return Ntp->Muon_numberOfMatches->at(i);}
  int            Muon_numberOfChambers(unsigned int i){return Ntp->Muon_numberOfChambers->at(i);}
  int            Muon_Charge(unsigned int i){return Ntp->Muon_Charge->at(i);}
  bool           isGoodMuon(unsigned int i);
  bool           isGoodMuon_nooverlapremoval(unsigned int i);
  
  //Base Tau Information (PF)
   unsigned int      NPFTaus(){return Ntp->PFTau_p4->size();}
   TLorentzVector    PFTau_p4(unsigned int i){return TLorentzVector(Ntp->PFTau_p4->at(i).at(1),Ntp->PFTau_p4->at(i).at(2),Ntp->PFTau_p4->at(i).at(3),Ntp->PFTau_p4->at(i).at(0));}
   bool              PFTau_isTightIsolation(unsigned int i){return Ntp->PFTau_isTightIsolation->at(i);}
   bool              PFTau_isMediumIsolation(unsigned int i){return Ntp->PFTau_isMediumIsolation->at(i);}
   bool              PFTau_isLooseIsolation(unsigned int i){return Ntp->PFTau_isLooseIsolation->at(i);}
   int               PFTau_hpsDecayMode(unsigned int i){return Ntp->PFTau_hpsDecayMode->at(i);}
   int               PFTau_Charge(unsigned int i){return Ntp->PFTau_Charge->at(i);}
   std::vector<int>  PFTau_Track_idx(unsigned int i){return Ntp->PFTau_Track_idx->at(i);}

   // Kinematic Fit Tau Information
   unsigned int     NKFTau(){return Ntp->KFTau_TauVis_p4->size();}
   TLorentzVector   KFTau_TauVis_p4(unsigned int i){return TLorentzVector(Ntp->KFTau_TauVis_p4->at(i).at(1),Ntp->KFTau_TauVis_p4->at(i).at(2),Ntp->KFTau_TauVis_p4->at(i).at(3),Ntp->KFTau_TauVis_p4->at(i).at(0));}
   TLorentzVector   KFTau_TauFit_p4(unsigned int i){return TLorentzVector(Ntp->KFTau_TauFit_p4->at(i).at(1),Ntp->KFTau_TauFit_p4->at(i).at(2),Ntp->KFTau_TauFit_p4->at(i).at(3),Ntp->KFTau_TauFit_p4->at(i).at(0));}
   TLorentzVector   KFTau_Neutrino_p4(unsigned int i){return TLorentzVector(Ntp->KFTau_Neutrino_p4->at(i).at(1),Ntp->KFTau_Neutrino_p4->at(i).at(2),Ntp->KFTau_Neutrino_p4->at(i).at(3),Ntp->KFTau_Neutrino_p4->at(i).at(0));}
   bool             KFTau_discriminatorByKFit(unsigned int i){return Ntp->KFTau_discriminatorByKFit->at(i);}
   bool             KFTau_discriminatorByQC(unsigned int i){return Ntp->KFTau_discriminatorByQC->at(i);}
   int              KFTau_nKinTaus(){return Ntp->KFTau_nKinTaus;}
   unsigned int     KFTau_MatchedHPS_idx(unsigned int i){return Ntp->KFTau_MatchedHPS_idx->at(i);}
   std::vector<int> KFTau_Track_idx(unsigned int i){return Ntp->KFTau_Track_idx->at(i);}
   int      KFTau_indexOfFitInfo(unsigned int i){return Ntp->KFTau_indexOfFitInfo->at(i);}
   TVector3 KFTau_Fit_TauPrimVtx(unsigned int i){return TVector3(Ntp->KFTau_Fit_TauPrimVtx->at(i).at(0),
								 Ntp->KFTau_Fit_TauPrimVtx->at(i).at(1),
								 Ntp->KFTau_Fit_TauPrimVtx->at(i).at(2));}
   int      KFTau_Fit_IndexToPrimVertexVector(unsigned int i){return Ntp->KFTau_Fit_IndexToPrimVertexVector->at(i);}
   float    KFTau_Fit_chi2(unsigned int i){return Ntp->KFTau_Fit_chi2->at(i);}
   float    KFTau_Fit_ndf(unsigned int i){return Ntp->KFTau_Fit_ndf->at(i);}
   int      KFTau_Fit_ambiguity(unsigned int i){return Ntp->KFTau_Fit_ambiguity->at(i);}
   int      KFTau_Fit_charge(unsigned int i){return Ntp->KFTau_Fit_charge->at(i);}
   int      KFTau_Fit_csum(unsigned int i){return Ntp->KFTau_Fit_csum->at(i);}
   int      KFTau_Fit_iterations(unsigned int i){return Ntp->KFTau_Fit_iterations->at(i);}


   // Jet Information
   unsigned int       NPFJets(){return Ntp->PFJet_p4->size();}
   TLorentzVector     PFJet_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_p4->at(i).at(1),Ntp->PFJet_p4->at(i).at(2),Ntp->PFJet_p4->at(i).at(3),Ntp->PFJet_p4->at(i).at(0));}
   float              PFJet_chargedEmEnergy(unsigned int i){return Ntp->PFJet_chargedEmEnergy->at(i);}
   float              PFJet_chargedHadronEnergy(unsigned int i){return Ntp->PFJet_chargedHadronEnergy->at(i);}
   float              PFJet_chargedHadronMultiplicity(unsigned int i){return Ntp->PFJet_chargedHadronMultiplicity->at(i);}
   float              PFJet_chargedMuEnergy(unsigned int i){return Ntp->PFJet_chargedMuEnergy->at(i);}
   float              PFJet_chargedMultiplicity(unsigned int i){return Ntp->PFJet_chargedMultiplicity->at(i);}
   float              PFJet_electronEnergy(unsigned int i){return Ntp->PFJet_electronEnergy->at(i);}
   float              PFJet_electronMultiplicity(unsigned int i){return Ntp->PFJet_electronMultiplicity->at(i);}
   float              PFJet_HFEMEnergy(unsigned int i){return Ntp->PFJet_HFEMEnergy->at(i);}
   float              PFJet_HFEMMultiplicity(unsigned int i){return Ntp->PFJet_HFEMMultiplicity->at(i);}
   float              PFJet_HFHadronEnergy(unsigned int i){return Ntp->PFJet_HFHadronEnergy->at(i);}
   float              PFJet_HFHadronMultiplicity(unsigned int i){return Ntp->PFJet_HFHadronMultiplicity->at(i);}
   float              PFJet_muonEnergy(unsigned int i){return Ntp->PFJet_muonEnergy->at(i);}
   float              PFJet_muonMultiplicity(unsigned int i){return Ntp->PFJet_muonMultiplicity->at(i);}
   float              PFJet_neutralEmEnergy(unsigned int i){return Ntp->PFJet_neutralEmEnergy->at(i);}
   float              PFJet_neutralHadronEnergy(unsigned int i){return Ntp->PFJet_neutralHadronEnergy->at(i);}
   float              PFJet_neutralHadronMultiplicity(unsigned int i){return Ntp->PFJet_neutralHadronMultiplicity->at(i);}
   float              PFJet_photonEnergy(unsigned int i){return Ntp->PFJet_photonEnergy->at(i);}
   float              PFJet_photonMultiplicity(unsigned int i){return Ntp->PFJet_photonMultiplicity->at(i);}
   float              PFJet_jetArea(unsigned int i){return Ntp->PFJet_jetArea->at(i);}
   float              PFJet_maxDistance(unsigned int i){return Ntp->PFJet_maxDistance->at(i);}
   int                PFJet_nConstituents(unsigned int i){return Ntp->PFJet_nConstituents->at(i);}
   float              PFJet_pileup(unsigned int i){return Ntp->PFJet_pileup->at(i);}
   float              PFJet_etaetaMoment(unsigned int i){return Ntp->PFJet_etaetaMoment->at(i);}
   float              PFJet_etaphiMoment(unsigned int i){return Ntp->PFJet_etaphiMoment->at(i);}
   std::vector<int>   PFJet_Track_idx(unsigned int i){return Ntp->PFJet_Track_idx->at(i);}
   unsigned int       PFJet_MatchedHPS_idx(unsigned int i){return Ntp->PFJet_MatchedHPS_idx->at(i);}
   int                PFJet_numberOfDaughters(unsigned int i){return Ntp->PFJet_numberOfDaughters->at(i);}
   float              PFJet_chargedEmEnergyFraction(unsigned int i){return Ntp->PFJet_chargedEmEnergyFraction->at(i);}
   float              PFJet_chargedHadronEnergyFraction(unsigned int i){return Ntp->PFJet_chargedHadronEnergyFraction->at(i);}
   float              PFJet_neutralHadronEnergyFraction(unsigned int i){return Ntp->PFJet_neutralHadronEnergyFraction->at(i);}
   float              PFJet_PFJet_neutralEmEnergyFraction(unsigned int i){return Ntp->PFJet_PFJet_neutralEmEnergyFraction->at(i);}
   bool               isGoodJet(unsigned int i);
   bool               isGoodJet_nooverlapremoval(unsigned int i);
   bool               isJetID(unsigned int i);


   //MET information
   double             MET_et(){return Ntp->MET_et;}
   double             MET_phi(){return Ntp->MET_phi;}
   double             MET_sumET(){return Ntp->MET_sumET;}

   //Track Information
   unsigned int      NTracks(){return Ntp->Track_p4->size();}
   TLorentzVector    Track_p4(unsigned int i){return TLorentzVector(Ntp->Track_p4->at(i).at(1),Ntp->Track_p4->at(i).at(2),Ntp->Track_p4->at(i).at(3),Ntp->Track_p4->at(i).at(0));}
   TVector3          Track_Poca(unsigned int i){return TVector3(Ntp->Track_Poca->at(i).at(0),Ntp->Track_Poca->at(i).at(1),Ntp->Track_Poca->at(i).at(2));}
   int               Track_charge(unsigned int i){return Ntp->Track_charge->at(i);}
   float             Track_chi2(unsigned int i){return Ntp->Track_chi2->at(i);}
   float             Track_ndof(unsigned int i){return Ntp->Track_ndof->at(i);}
   unsigned short    Track_numberOfLostHits(unsigned int i){return Ntp->Track_numberOfLostHits->at(i);}
   unsigned short    Track_numberOfValidHits(unsigned int i){return Ntp->Track_numberOfValidHits->at(i);}
   unsigned int      Track_qualityMask(unsigned int i){return Ntp->Track_qualityMask->at(i);}
   float             Track_par(unsigned int i, TrackPar par){return Ntp->Track_par->at(i).at(par);}
   TMatrixF          Track_parCov(unsigned int i);
   float             Track_parCov(unsigned int i, TrackPar par1, TrackPar par2);

   // MC Informaation
   // Singal particles (Z0,W+/-,H0,H+/-)
   unsigned int               NMCSignalParticles(){return Ntp->MCSignalParticle_p4->size();}
   TLorentzVector             MCSignalParticle_p4(unsigned int i){return TLorentzVector(Ntp->MCSignalParticle_p4->at(i).at(1),Ntp->MCSignalParticle_p4->at(i).at(2),Ntp->MCSignalParticle_p4->at(i).at(3),Ntp->MCSignalParticle_p4->at(i).at(0));}
   int                        MCSignalParticle_pdgid(unsigned int i){return Ntp->MCSignalParticle_pdgid->at(i);}
   int                        MCSignalParticle_charge(unsigned int i){return Ntp->MCSignalParticle_charge->at(i);}
   TVector3                   MCSignalParticle_Poca(unsigned int i){return TVector3(Ntp->MCSignalParticle_Poca->at(i).at(0),Ntp->MCSignalParticle_Poca->at(i).at(1),Ntp->MCSignalParticle_Poca->at(i).at(2));}
   std::vector<unsigned int>  MCSignalParticle_Tauidx(unsigned int i){return Ntp->MCSignalParticle_Tauidx->at(i);}

   // Tau decays (Tau is first element of vector)
   int NMCTaus(){return Ntp->MCTauandProd_p4->size();}
   TLorentzVector MCTau_p4(unsigned int i){return MCTauandProd_p4(i,0);}
   int MCTauandProd_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
   int MCTauandProd_charge(unsigned int i){return MCTauandProd_charge(i,0);}
   unsigned int MCTau_JAK(unsigned int i){return Ntp->MCTau_JAK->at(i);}
   unsigned int MCTau_DecayBitMask(unsigned int i){return Ntp->MCTau_DecayBitMask->at(i);}
   //Tau and decay products
   int NMCTauDecayProducts(unsigned int i){if(0<=i && i<NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
   TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
   int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
   unsigned int MCTauandProd_midx(unsigned int i, unsigned int j){return Ntp->MCTauandProd_midx->at(i).at(j);}
   int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}

};
#endif

