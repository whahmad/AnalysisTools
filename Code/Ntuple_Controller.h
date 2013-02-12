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
#include "TauSpinerInterface.h"
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"
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

  bool verbose;

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

  // Interfaces
  TauSpinerInterface TauSpinerInt;
  HistoConfig HistoC;

 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles);

  // Destructor
  ~Ntuple_Controller() ;

  //TauSpiner function
    double TauSpinerGet(TauSpinerInterface::TauSpinerType SpinType);
   void TauSpinerSetSignal(int signalcharge){TauSpinerInt.SetTauSignalCharge(signalcharge);}
  enum TrackQuality {
    undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
    confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
    qualitySize = 7
  };
  enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};

  enum KFTauPar{KFTau_vx = 0, KFTau_vy, KFTau_vz, KFTau_px, KFTau_py, KFTau_pz, KFTau_m, NKFTau_par};

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

  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}


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
  float          Muon_RelIso(unsigned int i);

  bool           Muon_isPFMuon(unsigned int i){return Ntp->Muon_isPFMuon->at(i);}                                                     
  float          Muon_sumChargedHadronPt03(unsigned int i){return Ntp->Muon_sumChargedHadronPt03->at(i);}                             // sum-pt of charged Hadron					                               
  float          Muon_sumChargedParticlePt03(unsigned int i){return Ntp->Muon_sumChargedParticlePt03->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			    
  float          Muon_sumNeutralHadronEt03(unsigned int i){return Ntp->Muon_sumNeutralHadronEt03->at(i);}			      // sum pt of neutral hadrons					    
  float          Muon_sumNeutralHadronEtHighThreshold03(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold03->at(i);}   // sum pt of neutral hadrons with a higher threshold		    
  float          Muon_sumPhotonEt03(unsigned int i){return Ntp->Muon_sumPhotonEt03->at(i);}					      // sum pt of PF photons					    
  float          Muon_sumPhotonEtHighThreshold03(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold03->at(i);}		      // sum pt of PF photons with a higher threshold		    
  float          Muon_sumPUPt03(unsigned int i){return Ntp->Muon_sumPUPt03->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)  
																      								    
  float          Muon_sumChargedHadronPt04(unsigned int i){return Ntp->Muon_sumChargedHadronPt04->at(i);}			      // sum-pt of charged Hadron					    
  float          Muon_sumChargedParticlePt04(unsigned int i){return Ntp->Muon_sumChargedParticlePt04->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			    
  float          Muon_sumNeutralHadronEt04(unsigned int i){return Ntp->Muon_sumNeutralHadronEt04->at(i);}			      // sum pt of neutral hadrons					    
  float          Muon_sumNeutralHadronEtHighThreshold04(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold04->at(i);}   // sum pt of neutral hadrons with a higher threshold		    
  float          Muon_sumPhotonEt04(unsigned int i){return Ntp->Muon_sumPhotonEt04->at(i);}					      // sum pt of PF photons					    
  float          Muon_sumPhotonEtHighThreshold04(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold04->at(i);}		      // sum pt of PF photons with a higher threshold		    
  float          Muon_sumPUPt04(unsigned int i){return Ntp->Muon_sumPUPt04->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)  

  int            Muon_numberofValidPixelHits(unsigned int i){return Ntp->Muon_numberofValidPixelHits->at(i);}
  int            Muon_trackerLayersWithMeasurement(unsigned int i){return Ntp->Muon_trackerLayersWithMeasurement->at(i);}



  //Base Tau Information (PF)
   unsigned int      NPFTaus(){return Ntp->PFTau_p4->size();}
   TLorentzVector    PFTau_p4(unsigned int i){return TLorentzVector(Ntp->PFTau_p4->at(i).at(1),Ntp->PFTau_p4->at(i).at(2),Ntp->PFTau_p4->at(i).at(3),Ntp->PFTau_p4->at(i).at(0));}
   TVector3          PFTau_Poca(unsigned int i){return TVector3(Ntp->PFTau_Poca->at(i).at(0),Ntp->PFTau_Poca->at(i).at(1),Ntp->PFTau_Poca->at(i).at(2));}
   bool              PFTau_isTightIsolation(unsigned int i){return Ntp->PFTau_isTightIsolation->at(i);}
   bool              PFTau_isMediumIsolation(unsigned int i){return Ntp->PFTau_isMediumIsolation->at(i);}
   bool              PFTau_isLooseIsolation(unsigned int i){return Ntp->PFTau_isLooseIsolation->at(i);}
   int               PFTau_hpsDecayMode(unsigned int i){return Ntp->PFTau_hpsDecayMode->at(i);}
   int               PFTau_Charge(unsigned int i){return Ntp->PFTau_Charge->at(i);}

   bool              PFTau_isTightIsolationDBSumPtCorr(unsigned int i){return Ntp->PFTau_isTightIsolationDBSumPtCorr->at(i);}
   bool              PFTau_isMediumIsolationDBSumPtCorr(unsigned int i){return Ntp->PFTau_isMediumIsolationDBSumPtCorr->at(i);}
   bool              PFTau_isLooseIsolationDBSumPtCorr(unsigned int i){return Ntp->PFTau_isLooseIsolationDBSumPtCorr->at(i);}
   bool              PFTau_isVLooseIsolationDBSumPtCorr(unsigned int i){return Ntp->PFTau_isVLooseIsolationDBSumPtCorr->at(i);}
   
   bool              PFTau_isHPSAgainstElectronsLoose(unsigned int i){return Ntp->PFTau_isHPSAgainstElectronsLoose->at(i);}
   bool              PFTau_isHPSAgainstElectronsMedium(unsigned int i){return Ntp->PFTau_isHPSAgainstElectronsMedium->at(i);}
   bool              PFTau_isHPSAgainstElectronsTight(unsigned int i){return Ntp->PFTau_isHPSAgainstElectronsTight->at(i);}
   bool              PFTau_isHPSAgainstMuonLoose(unsigned int i){return Ntp->PFTau_isHPSAgainstMuonLoose->at(i);}
   bool              PFTau_isHPSAgainstMuonTight(unsigned int i){return Ntp->PFTau_isHPSAgainstMuonTight->at(i);}
   bool              PFTau_isHPSByDecayModeFinding(unsigned int i){return Ntp->PFTau_isHPSByDecayModeFinding->at(i);}
   


   bool              PFTau_isHPSAgainstMuonLoose2(unsigned int i){return Ntp->PFTau_isHPSAgainstMuonLoose2->at(i);}
   bool	             PFTau_isHPSAgainstMuonMedium2(unsigned int i){return Ntp->PFTau_isHPSAgainstMuonMedium2->at(i);}
   bool	             PFTau_isHPSAgainstMuonTight2(unsigned int i){return Ntp->PFTau_isHPSAgainstMuonTight2->at(i);}
													
   bool	             PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection->at(i);}
													
   bool	             PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByTightIsolationMVA(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA->at(i);}
													
   bool	             PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2->at(i);}
   bool	             PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2->at(i);}

   std::vector<int>  PFTau_Track_idx(unsigned int i){return Ntp->PFTau_Track_idx->at(i);}




   // Kinematic Fit Tau Information
   unsigned int     NKFTau(){return Ntp->KFTau_TauVis_p4->size();}
   bool             isGoodKFTau(unsigned int i, unsigned int j);
   int KFTauNSolutions1(unsigned int i) {return  Ntp->KFTau_TauVis_p4->at(i).size();}
   int KFTauNSolutions11(unsigned int i) {return  Ntp->Muon_Charge->size();}
   int KFTauNSolutions2(unsigned int i) {return  Ntp->KFTau_Fit_charge->size();}

  
   TLorentzVector   KFTau_TauVis_p4(unsigned int i,unsigned int j){return TLorentzVector(Ntp->KFTau_TauVis_p4->at(i).at(j).at(1),Ntp->KFTau_TauVis_p4->at(i).at(j).at(2),Ntp->KFTau_TauVis_p4->at(i).at(j).at(3),Ntp->KFTau_TauVis_p4->at(i).at(j).at(0));}
   TLorentzVector   KFTau_TauFit_p4(unsigned int i,unsigned int j){return TLorentzVector(Ntp->KFTau_TauFit_p4->at(i).at(j).at(1),Ntp->KFTau_TauFit_p4->at(i).at(j).at(2),Ntp->KFTau_TauFit_p4->at(i).at(j).at(3),Ntp->KFTau_TauFit_p4->at(i).at(j).at(0));}
   TLorentzVector   KFTau_TauFitInitial_p4(unsigned int i,unsigned int j){return TLorentzVector(Ntp->KFTau_TauFitInitial_p4->at(i).at(j).at(1),Ntp->KFTau_TauFitInitial_p4->at(i).at(j).at(2),Ntp->KFTau_TauFitInitial_p4->at(i).at(j).at(3),Ntp->KFTau_TauFitInitial_p4->at(i).at(j).at(0));}
   TLorentzVector   KFTau_Neutrino_p4(unsigned int i,unsigned int j){return TLorentzVector(Ntp->KFTau_Neutrino_p4->at(i).at(j).at(1),Ntp->KFTau_Neutrino_p4->at(i).at(j).at(2),Ntp->KFTau_Neutrino_p4->at(i).at(j).at(3),Ntp->KFTau_Neutrino_p4->at(i).at(j).at(0));}
   TLorentzVector   KFTau_NeutrinoInitial_p4(unsigned int i,unsigned int j){return TLorentzVector(Ntp->KFTau_NeutrinoInitial_p4->at(i).at(j).at(1),Ntp->KFTau_NeutrinoInitial_p4->at(i).at(j).at(2),Ntp->KFTau_NeutrinoInitial_p4->at(i).at(j).at(3),Ntp->KFTau_NeutrinoInitial_p4->at(i).at(j).at(0));}
   TLorentzVector   KFTau_a1Initial_p4(unsigned int i){return TLorentzVector(Ntp->KFTau_a1Initial_p4->at(i).at(1),Ntp->KFTau_a1Initial_p4->at(i).at(2),Ntp->KFTau_a1Initial_p4->at(i).at(3),Ntp->KFTau_a1Initial_p4->at(i).at(0));}
   int              KFTau_discriminatorByKFit(unsigned int i,unsigned int j){return Ntp->KFTau_discriminatorByKFit->at(i).at(j);}
   int              KFTau_discriminatorByQC(unsigned int i,unsigned int j){return Ntp->KFTau_discriminatorByQC->at(i).at(j);}
   int              KFTau_nKinTaus(){return Ntp->KFTau_nKinTaus;}
   unsigned int     KFTau_MatchedHPS_idx(unsigned int i){return Ntp->KFTau_MatchedHPS_idx->at(i);}
   std::vector<int> KFTau_Track_idx(unsigned int i){return Ntp->KFTau_Track_idx->at(i);}


   int      KFTau_indexOfFitInfo(unsigned int i){return Ntp->KFTau_indexOfFitInfo->at(i);}

   TVector3 KFTau_Fit_TauPrimVtx(unsigned int i){return TVector3(Ntp->KFTau_Fit_TauPrimVtx->at(i).at(0),Ntp->KFTau_Fit_TauPrimVtx->at(i).at(1),Ntp->KFTau_Fit_TauPrimVtx->at(i).at(2));}

   TVector3 KFTau_Fit_PrimaryVertex(unsigned int i){return TVector3(Ntp->KFTau_Fit_PrimaryVertex->at(i).at(0),Ntp->KFTau_Fit_PrimaryVertex->at(i).at(1),Ntp->KFTau_Fit_PrimaryVertex->at(i).at(2));}
   TVector3 KFTau_Fit_InitialPrimaryVertex(unsigned int i){return TVector3(Ntp->KFTau_Fit_InitialPrimaryVertex->at(i).at(0),Ntp->KFTau_Fit_InitialPrimaryVertex->at(i).at(1),Ntp->KFTau_Fit_InitialPrimaryVertex->at(i).at(2));}

   TVector3 KFTau_Fit_InitialSecondaryVertex(unsigned int i){return TVector3(Ntp->KFTau_Fit_InitialSecondaryVertex->at(i).at(0),Ntp->KFTau_Fit_InitialSecondaryVertex->at(i).at(1),Ntp->KFTau_Fit_InitialSecondaryVertex->at(i).at(2));}
   TVector3 KFTau_Fit_SecondaryVertex(unsigned int i,unsigned int j){return TVector3(Ntp->KFTau_Fit_SecondaryVertex->at(i).at(j).at(0),Ntp->KFTau_Fit_SecondaryVertex->at(i).at(j).at(1),Ntp->KFTau_Fit_SecondaryVertex->at(i).at(j).at(2));}

   //   float    KFTau_Fit_chi2(unsigned int i){return Ntp->KFTau_Fit_chi2->at(i);}
   float    KFTau_Fit_ndf(unsigned int i,unsigned int j){return Ntp->KFTau_Fit_ndf->at(i).at(j);}
   //  int      KFTau_Fit_ambiguity(unsigned int i){return Ntp->KFTau_Fit_ambiguity->at(i);}
   int      KFTau_Fit_charge(unsigned int i){return Ntp->KFTau_Fit_charge->at(i);}

   int      KFTau_Fit_csum(unsigned int i,unsigned int j){return Ntp->KFTau_Fit_csum->at(i).at(j);}
   int      KFTau_Fit_iterations(unsigned int i,unsigned int j){return Ntp->KFTau_Fit_iterations->at(i).at(j);}

   double   KFTau_Fit_TauEnergyFraction(unsigned int i,unsigned int j){return Ntp->KFTau_Fit_TauEnergyFraction->at(i).at(j);}
   double   KFTau_Fit_RefitVisibleMass(unsigned int i, unsigned int j){return Ntp->KFTau_Fit_RefitVisibleMass->at(i).at(j);}
   double   KFTau_Fit_Chi2Prob(unsigned int i, unsigned int j){return Ntp->KFTau_Fit_Chi2Prob->at(i).at(j);}
   float    KFTau_Fit_chi2(unsigned int i, unsigned int j){return Ntp->KFTau_Fit_chi2->at(i).at(j);}

   float    KFTau_Fit_BDTVal(unsigned int i, unsigned int j){return Ntp->KFTau_Fit_BDTVal->at(i).at(j);}


   double   KFTau_Fit_PV_PV_significance(unsigned int i, unsigned int j){return Ntp->KFTau_Fit_PV_PV_significance->at(i).at(j);}
   double   KFTau_Fit_SV_PV_significance(unsigned int i, unsigned int j){return Ntp->KFTau_Fit_SV_PV_significance->at(i).at(j);}

   unsigned int KFTau_NDaughter(unsigned int i){return Ntp->KFTau_Daughter_pdgid->at(i).size();}
   int KFTau_Daughter_pdgid(unsigned int i, unsigned int j){return Ntp->KFTau_Daughter_pdgid->at(i).at(j);}
   int KFTau_Daughter_charge(unsigned int i, unsigned int j){return Ntp->KFTau_Daughter_charge->at(i).at(j);}
   float KFTau_Daughter_ambiguity(unsigned int i, unsigned int j){return Ntp->KFTau_Daughter_ambiguity->at(i).at(j);}
   float KFTau_Daughter_par(unsigned int i, unsigned int j,int par){return Ntp->KFTau_Daughter_par->at(i).at(j).at(par);}
   float KFTau_Daughter_parCov(unsigned int i, unsigned int j,int par1,int par2);
   float KFTau_Daughter_inputpar(unsigned int i, unsigned int j,int par){return Ntp->KFTau_Daughter_par->at(i).at(j).at(par);}
   float KFTau_Daughter_inputparCov(unsigned int i, unsigned int j,int par1,int par2);

   TVector3  KFTau_RotatedVtx(unsigned int i);
   TMatrixF  KFTau_RotatedVtx_Cov(unsigned int i);

   TVector3     KFTau_ReducedVtx(){unsigned int i=0;// only 1 reduced vertex for now
     return TVector3(Ntp->ReducedVtx_x->at(i),Ntp->ReducedVtx_y->at(i),Ntp->ReducedVtx_z->at(i));
   }
   float   KFTau_ReducedVtx_chi2(){unsigned int i=0;return Ntp->ReducedVtx_chi2->at(i);}
   float   KFTau_ReducedVtx_nTrk(){unsigned int i=0;return Ntp->ReducedVtx_nTrk->at(i);}
   float   KFTau_ReducedVtx_ndof(){unsigned int i=0;return Ntp->ReducedVtx_ndof->at(i);}
   TMatrixF  KFTau_ReducedVtx_Cov();
   std::vector<int> KFTau_ReducedVtx_Track_idx(){unsigned int i=0;return Ntp->ReducedVtx_Track_idx->at(i);}
   float   KFTau_ReducedVtx_isFake(){unsigned int i=0;return Ntp->ReducedVtx_isFake->at(i);}

   TVector3  KFTau_SecondayVtx(unsigned int i);
   TMatrixF  KFTau_SecondaryVtx_Cov(unsigned int i);

   TVector3  KFTau_InitialSecondaryVtx(unsigned int i);
   TMatrixF  KFTau_InitialSecondaryVtx_Cov(unsigned int i);

/*   std::vector<std::vector<std::vector<float> > >  KFTau_Fit_SecondaryVertex; */
/*   std::vector<std::vector<float> >  KFTau_Fit_InitialSecondaryVertex; */


/*   std::vector<std::vector<float> > KFTau_Fit_TauEnergyFraction; */
/*   std::vector<std::vector<float> > KFTau_Fit_RefitVisibleMass; */
/*   std::vector<std::vector<float> > KFTau_Fit_Chi2Prob; */
/*   std::vector<std::vector<float> > KFTau_Fit_chi2; */
/*   std::vector<std::vector<float> > KFTau_Fit_BDTVal; */
/*   std::vector<std::vector<float> > KFTau_Fit_PV_PV_significance; */
/*   std::vector<std::vector<float> > KFTau_Fit_SV_PV_significance; */

/*   std::vector<std::vector<int> > KFTau_Daughter_pdgid; */
/*   std::vector<std::vector<int> > KFTau_Daughter_charge; */
/*   std::vector<std::vector<float> > KFTau_Daughter_ambiguity; */

/*   std::vector<std::vector<std::vector<float> > > KFTau_Daughter_par; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_Daughter_parCov; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_Daughter_inputpar; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_Daughter_inputparCov; */

/*   std::vector<float> ReducedVtx_chi2; */
/*   std::vector<float> ReducedVtx_nTrk; */
/*   std::vector<float> ReducedVtx_ndof; */
/*   std::vector<float> ReducedVtx_y; */
/*   std::vector<float> ReducedVtx_x; */
/*   std::vector<float> ReducedVtx_z; */
/*   std::vector<std::vector<std::vector<float> > >  ReducedVtx_Cov; */
/*   std::vector<std::vector<int> >    ReducedVtx_Track_idx; */
/*   std::vector<float> ReducedVtx_isFake; */









/*     std::vector<std::vector<int> > KFTau_discriminatorByKFit; */
/*   std::vector<std::vector<int> > KFTau_discriminatorByQC; */


/*   std::vector<std::vector<std::vector<float> > > KFTau_TauVis_p4; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_TauFit_p4; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_TauFitInitial_p4; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_Neutrino_p4; */
/*   std::vector<std::vector<std::vector<float> > > KFTau_NeutrinoInitial_p4; */
/*   std::vector<std::vector<float> >  KFTau_a1Initial_p4; */
/*   std::vector<std::vector<int> > KFTau_Track_idx; */
/*   std::vector<int> KFTau_indexOfFitInfo; */
/*   std::vector<std::vector<float> >	KFTau_Fit_ndf; */
/*   // std::vector<int>	KFTau_Fit_ambiguity; */
/*   std::vector<int>	KFTau_Fit_charge; */
/*   std::vector<std::vector<int> >	KFTau_Fit_csum;      */
/*   std::vector<std::vector<int> >     KFTau_Fit_iterations; */
/*   std::vector<std::vector<float> > KFTau_Fit_TauPrimVtx; */


/*   std::vector<std::vector<float> > KFTau_Fit_PrimaryVertex; */
/*   std::vector<std::vector<float> > KFTau_Fit_InitialPrimaryVertex; */


   // Jet Information
   unsigned int       NPFJets(){return Ntp->PFJet_p4->size();}
   TLorentzVector     PFJet_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_p4->at(i).at(1),Ntp->PFJet_p4->at(i).at(2),Ntp->PFJet_p4->at(i).at(3),Ntp->PFJet_p4->at(i).at(0));}
   TVector3           PFJet_Poca(unsigned int i){return TVector3(Ntp->PFJet_Poca->at(i).at(0),Ntp->PFJet_Poca->at(i).at(1),Ntp->PFJet_Poca->at(i).at(2));}
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
   float              PFJet_neutralEmEnergyFraction(unsigned int i){return Ntp->PFJet_neutralEmEnergyFraction->at(i);}
   bool               isGoodJet(unsigned int i);
   bool               isGoodJet_nooverlapremoval(unsigned int i);
   bool               isJetID(unsigned int i);
 

   //MET information
   double             MET_et(){return Ntp->MET_et;}
   double             MET_phi(){return Ntp->MET_phi;}
   double             MET_sumET(){return Ntp->MET_sumET;}
   double             MET_ex(){return Ntp->MET_et*cos(Ntp->MET_phi);}
   double             MET_ey(){return Ntp->MET_et*sin(Ntp->MET_phi);}

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
   int MCTau_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
   int MCTau_charge(unsigned int i){return MCTauandProd_charge(i,0);}
   unsigned int MCTau_JAK(unsigned int i){return Ntp->MCTau_JAK->at(i);}
   unsigned int MCTau_DecayBitMask(unsigned int i){return Ntp->MCTau_DecayBitMask->at(i);}
   //Tau and decay products
   int NMCTauDecayProducts(unsigned int i){if(0<=i && i<NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
   TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
   int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
   unsigned int MCTauandProd_midx(unsigned int i, unsigned int j){return Ntp->MCTauandProd_midx->at(i).at(j);}
   int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}
   bool hasSignalTauDecay(PdtPdgMini::PdgPDTMini parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &idx);

   bool jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx);
   bool muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx);
   bool muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx);


   // Electrons
   unsigned int       NElectrons(){return Ntp->Electron_p4->size();}
   TLorentzVector     Electron_p4(unsigned int i){return TLorentzVector(Ntp->Electron_p4->at(i).at(1),Ntp->Electron_p4->at(i).at(2),Ntp->Electron_p4->at(i).at(3),Ntp->Electron_p4->at(i).at(0));}
   TVector3           Electron_Poca(unsigned int i){return TVector3(Ntp->Electron_Poca->at(i).at(0),Ntp->Electron_Poca->at(i).at(1),Ntp->Electron_Poca->at(i).at(2));}
   float   Electron_Charge(unsigned int i){return Ntp->Electron_Charge->at(i);}
   float   Electron_Gsf_deltaEtaEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaEleClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaEtaSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSeedClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx->at(i);}
   float   Electron_Gsf_deltaPhiEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiEleClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaPhiSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSeedClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx->at(i);}
   float   Electron_Gsf_dr03EcalRecHitSumE(unsigned int i){return Ntp->Electron_Gsf_dr03EcalRecHitSumE->at(i);}
   float   Electron_Gsf_dr03HcalDepth1TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalDepth1TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03HcalDepth2TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalDepth2TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03HcalTowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalTowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03TkSumPt(unsigned int i){return Ntp->Electron_Gsf_dr03TkSumPt->at(i);}
   bool    Electron_Gsf_passingCutBasedPreselection(unsigned int i){return Ntp->Electron_Gsf_passingCutBasedPreselection->at(i);}
   bool    Electron_Gsf_passingMvaPreselection(unsigned int i){return Ntp->Electron_Gsf_passingMvaPreselection->at(i);}
   int     Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits(unsigned int i){return Ntp->Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits->at(i);}
   float   Electron_supercluster_e(unsigned int i){return Ntp->Electron_supercluster_e->at(i);}
   float   Electron_supercluster_phi(unsigned int i){return Ntp->Electron_supercluster_phi->at(i);}
   float   Electron_supercluster_eta(unsigned int i){return Ntp->Electron_supercluster_eta->at(i);}
   float   Electron_supercluster_centroid_x(unsigned int i){return Ntp->Electron_supercluster_centroid_x->at(i);}
   float   Electron_supercluster_centroid_y(unsigned int i){return Ntp->Electron_supercluster_centroid_y->at(i);}
   float   Electron_supercluster_centroid_z(unsigned int i){return Ntp->Electron_supercluster_centroid_z->at(i);}
   unsigned int Electron_Track_idx(unsigned int i){return Ntp->Electron_Track_idx->at(i);}


   float    Electron_ecalRecHitSumEt03(unsigned int i){return Ntp->Electron_ecalRecHitSumEt03->at(i);}
   float    Electron_hcalDepth1TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt03->at(i);}
   float    Electron_hcalDepth1TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc03->at(i);}
   float    Electron_hcalDepth2TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt03->at(i);}
   float    Electron_hcalDepth2TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc03->at(i);}
   float    Electron_tkSumPt03(unsigned int i){return Ntp->Electron_tkSumPt03->at(i);}
   float    Electron_ecalRecHitSumEt04(unsigned int i){return Ntp->Electron_ecalRecHitSumEt04->at(i);}
   float    Electron_hcalDepth1TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt04->at(i);}
   float    Electron_hcalDepth1TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc04->at(i);}
   float    Electron_hcalDepth2TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt04->at(i);}
   float    Electron_hcalDepth2TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc04->at(i);}
   float    Electron_tkSumPt04(unsigned int i){return Ntp->Electron_tkSumPt04->at(i);}
   
   
   float    Electron_chargedHadronIso(unsigned int i){return Ntp->Electron_chargedHadronIso->at(i);}
   float    Electron_neutralHadronIso(unsigned int i){return Ntp->Electron_neutralHadronIso->at(i);}
   float    Electron_photonIso(unsigned int i){return Ntp->Electron_photonIso->at(i);}
     	 						    
   
   float    Electron_sigmaIetaIeta(unsigned int i){return Ntp->Electron_sigmaIetaIeta->at(i);}
   float    Electron_hadronicOverEm(unsigned int i){return Ntp->Electron_hadronicOverEm->at(i);}
   float    Electron_fbrem(unsigned int i){return Ntp->Electron_fbrem->at(i);}
   float    Electron_eSuperClusterOverP(unsigned int i){return Ntp->Electron_eSuperClusterOverP->at(i);}
   float    Electron_ecalEnergy(unsigned int i){return Ntp->Electron_ecalEnergy->at(i);}
   float    Electron_trackMomentumAtVtx(unsigned int i){return Ntp->Electron_trackMomentumAtVtx->at(i);}
   float    Electron_numberOfMissedHits(unsigned int i){return Ntp->Electron_numberOfMissedHits->at(i);}
   bool     Electron_HasMatchedConversions(unsigned int i){return Ntp->Electron_HasMatchedConversions->at(i);}





   // Trigger
   bool         TriggerAccept(TString n);
   unsigned int HLTPrescale(TString n);
   unsigned int L1SEEDPrescale(TString n);
   bool         GetTriggerIndex(TString n, unsigned int &i);
   unsigned int NHLTTriggers(){return Ntp->HTLTriggerName->size();}
   std::string  HTLTriggerName(unsigned int i){return Ntp->HTLTriggerName->at(i);}
   bool         TriggerAccept(unsigned int i){return Ntp->TriggerAccept->at(i);}
   bool         TriggerError(unsigned int i){return Ntp->TriggerError->at(i);}
   bool         TriggerWasRun(unsigned int i){return Ntp->TriggerWasRun->at(i);}
   unsigned int HLTPrescale(unsigned int i){return Ntp->HLTPrescale->at(i);}
   unsigned int NHLTL1GTSeeds(unsigned int i){return Ntp->NHLTL1GTSeeds->at(i);}
   unsigned int L1SEEDPrescale(unsigned int i){return Ntp->L1SEEDPrescale->at(i);}
   bool         L1SEEDInvalidPrescale(unsigned int i){return Ntp->L1SEEDInvalidPrescale->at(i);}
   float        MuonTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->MuonTriggerMatch->at(i).size()) return Ntp->MuonTriggerMatch->at(i).at(j); return 999;}
   float        ElectronTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->ElectronTriggerMatch->at(i).size()) return Ntp->ElectronTriggerMatch->at(i).at(j);return 999;}
   float        JetTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->JetTriggerMatch->at(i).size()) return Ntp->JetTriggerMatch->at(i).at(j);return 999;}
   float        TauTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->TauTriggerMatch->at(i).size()) return Ntp->TauTriggerMatch->at(i).at(j);return 999;}
   unsigned int NHLTTriggerObject(unsigned int i){return Ntp->HLTTrigger_objs_Eta->at(i).size();}
   TLorentzVector HLTTriggerObject_p4(unsigned int i, unsigned int j){
     TLorentzVector L(0,0,0,0); 
     if(j<Ntp->HLTTrigger_objs_Eta->at(i).size())L.SetPtEtaPhiM(Ntp->HLTTrigger_objs_Pt->at(i).at(j),Ntp->HLTTrigger_objs_Eta->at(i).at(j), Ntp->HLTTrigger_objs_Phi->at(i).at(j),0.0);
     return L;
   }


};

#endif

