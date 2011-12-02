//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  2 09:37:16 2011 by ROOT version 5.26/00c
// from TChain t/
//////////////////////////////////////////////////////////

#ifndef NtupleReader_h
#define NtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TFile.h>
#include <TTreeCache.h>
#include <vector>
#include <string>
#include<iostream>


class NtupleReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   std::vector<float>   *Vtx_chi2;
   std::vector<float>   *Vtx_nTrk;
   std::vector<float>   *Vtx_ndof;
   std::vector<float>   *Vtx_x;
   std::vector<float>   *Vtx_y;
   std::vector<float>   *Vtx_z;
   std::vector<std::vector<float> > *Vtx_Cov;
   std::vector<std::vector<int> > *Vtx_Track_idx;
   std::vector<std::vector<float> > *Muon_p4;
   std::vector<std::vector<double> > *Muon_Poca;
   std::vector<bool>    *Muon_isGlobalMuon;
   std::vector<bool>    *Muon_isStandAloneMuon;
   std::vector<bool>    *Muon_isTrackerMuon;
   std::vector<bool>    *Muon_isCaloMuon;
   std::vector<bool>    *Muon_isIsolationValid;
   std::vector<bool>    *Muon_isQualityValid;
   std::vector<bool>    *Muon_isTimeValid;
   std::vector<float>   *Muon_emEt03;
   std::vector<float>   *Muon_emVetoEt03;
   std::vector<float>   *Muon_hadEt03;
   std::vector<float>   *Muon_hadVetoEt03;
   std::vector<float>   *Muon_nJets03;
   std::vector<float>   *Muon_nTracks03;
   std::vector<float>   *Muon_sumPt03;
   std::vector<float>   *Muon_trackerVetoPt03;
   std::vector<float>   *Muon_emEt05;
   std::vector<float>   *Muon_emVetoEt05;
   std::vector<float>   *Muon_hadEt05;
   std::vector<float>   *Muon_hadVetoEt05;
   std::vector<float>   *Muon_nJets05;
   std::vector<float>   *Muon_nTracks05;
   std::vector<float>   *Muon_sumPt05;
   std::vector<float>   *Muon_trackerVetoPt05;
   std::vector<unsigned int> *Muon_Track_idx;
   std::vector<std::vector<float> > *PFTau_p4;
   std::vector<bool>    *PFTau_isTightIsolation;
   std::vector<bool>    *PFTau_isMediumIsolation;
   std::vector<bool>    *PFTau_isLooseIsolation;
   std::vector<int>     *PFTau_hpsDecayMode;
   std::vector<int>     *PFTau_Charge;
   std::vector<std::vector<int> > *PFTau_Track_idx;
   std::vector<bool>    *KFTau_discriminatorByKFit;
   std::vector<bool>    *KFTau_discriminatorByQC;
   Int_t           KFTau_nKinTaus;
   std::vector<std::vector<float> > *KFTau_TauVis_p4;
   std::vector<std::vector<float> > *KFTau_TauFit_p4;
   std::vector<std::vector<float> > *KFTau_Neutrino_p4;
   std::vector<unsigned int> *KFTau_MatchedHPS_idx;
   std::vector<std::vector<int> > *KFTau_Track_idx;
   std::vector<int>     *KFTau_indexOfFitInfo;
   std::vector<std::vector<float> > *KFTau_Fit_TauPrimVtx;
   std::vector<int>     *KFTau_Fit_IndexToPrimVertexVector;
   std::vector<float>   *KFTau_Fit_chi2;
   std::vector<float>   *KFTau_Fit_ndf;
   std::vector<int>     *KFTau_Fit_ambiguity;
   std::vector<int>     *KFTau_Fit_charge;
   std::vector<int>     *KFTau_Fit_csum;
   std::vector<int>     *KFTau_Fit_iterations;
   std::vector<std::vector<float> > *PFJet_p4;
   std::vector<float>   *PFJet_chargedEmEnergy;
   std::vector<float>   *PFJet_chargedHadronEnergy;
   std::vector<float>   *PFJet_chargedHadronMultiplicity;
   std::vector<float>   *PFJet_chargedMuEnergy;
   std::vector<float>   *PFJet_chargedMultiplicity;
   std::vector<float>   *PFJet_electronEnergy;
   std::vector<float>   *PFJet_electronMultiplicity;
   std::vector<float>   *PFJet_HFEMEnergy;
   std::vector<float>   *PFJet_HFEMMultiplicity;
   std::vector<float>   *PFJet_HFHadronEnergy;
   std::vector<float>   *PFJet_HFHadronMultiplicity;
   std::vector<float>   *PFJet_muonEnergy;
   std::vector<float>   *PFJet_muonMultiplicity;
   std::vector<float>   *PFJet_neutralEmEnergy;
   std::vector<float>   *PFJet_neutralHadronEnergy;
   std::vector<float>   *PFJet_neutralHadronMultiplicity;
   std::vector<float>   *PFJet_photonEnergy;
   std::vector<float>   *PFJet_photonMultiplicity;
   std::vector<float>   *PFJet_jetArea;
   std::vector<float>   *PFJet_maxDistance;
   std::vector<int>     *PFJet_nConstituents;
   std::vector<float>   *PFJet_pileup;
   std::vector<float>   *PFJet_etaetaMoment;
   std::vector<float>   *PFJet_etaphiMoment;
   std::vector<std::vector<int> > *PFJet_Track_idx;
   std::vector<unsigned int> *PFJet_MatchedHPS_idx;
   Double_t        MET_et;
   Double_t        MET_phi;
   Double_t        MET_sumET;
   UInt_t          Event_EventNumber;
   UInt_t          Event_RunNumber;
   Int_t           Event_bunchCrossing;
   Int_t           Event_orbitNumber;
   UInt_t          Event_luminosityBlock;
   Bool_t          Event_isRealData;
   std::vector<std::vector<float> > *Track_p4;
   std::vector<std::vector<double> > *Track_Poca;
   std::vector<int>     *Track_charge;
   std::vector<float>   *Track_chi2;
   std::vector<float>   *Track_ndof;
   std::vector<unsigned short> *Track_numberOfLostHits;
   std::vector<unsigned short> *Track_numberOfValidHits;
   std::vector<unsigned int> *Track_qualityMask;
   std::vector<std::vector<float> > *MCSignalParticle_p4;
   std::vector<int>     *MCSignalParticle_pdgid;
   std::vector<int>     *MCSignalParticle_charge;
   std::vector<std::vector<double> > *MCSignalParticle_Poca;
   std::vector<std::vector<unsigned int> > *MCSignalParticle_Tauidx;
   std::vector<std::vector<std::vector<float> > > *MCTauandProd_p4;
   std::vector<std::vector<int> > *MCTauandProd_pdgid;
   std::vector<std::vector<unsigned int> > *MCTauandProd_midx;
   std::vector<std::vector<int> > *MCTauandProd_charge;
   std::vector<unsigned int> *MCTau_JAK;
   std::vector<unsigned int> *MCTau_DecayBitMask;

   // List of branches
   TBranch        *b_Vtx_chi2;   //!
   TBranch        *b_Vtx_nTrk;   //!
   TBranch        *b_Vtx_ndof;   //!
   TBranch        *b_Vtx_x;   //!
   TBranch        *b_Vtx_y;   //!
   TBranch        *b_Vtx_z;   //!
   TBranch        *b_Vtx_Cov;   //!
   TBranch        *b_Vtx_Track_idx;   //!
   TBranch        *b_Muon_p4;   //!
   TBranch        *b_Muon_Poca;   //!
   TBranch        *b_Muon_isGlobalMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isCaloMuon;   //!
   TBranch        *b_Muon_isIsolationValid;   //!
   TBranch        *b_Muon_isQualityValid;   //!
   TBranch        *b_Muon_isTimeValid;   //!
   TBranch        *b_Muon_emEt03;   //!
   TBranch        *b_Muon_emVetoEt03;   //!
   TBranch        *b_Muon_hadEt03;   //!
   TBranch        *b_Muon_hadVetoEt03;   //!
   TBranch        *b_Muon_nJets03;   //!
   TBranch        *b_Muon_nTracks03;   //!
   TBranch        *b_Muon_sumPt03;   //!
   TBranch        *b_Muon_trackerVetoPt03;   //!
   TBranch        *b_Muon_emEt05;   //!
   TBranch        *b_Muon_emVetoEt05;   //!
   TBranch        *b_Muon_hadEt05;   //!
   TBranch        *b_Muon_hadVetoEt05;   //!
   TBranch        *b_Muon_nJets05;   //!
   TBranch        *b_Muon_nTracks05;   //!
   TBranch        *b_Muon_sumPt05;   //!
   TBranch        *b_Muon_trackerVetoPt05;   //!
   TBranch        *b_Muon_Track_idx;   //!
   TBranch        *b_PFTau_p4;   //!
   TBranch        *b_PFTau_isTightIsolation;   //!
   TBranch        *b_PFTau_isMediumIsolation;   //!
   TBranch        *b_PFTau_isLooseIsolation;   //!
   TBranch        *b_PFTau_hpsDecayMode;   //!
   TBranch        *b_PFTau_Charge;   //!
   TBranch        *b_PFTau_Track_idx;   //!
   TBranch        *b_KFTau_discriminatorByKFit;   //!
   TBranch        *b_KFTau_discriminatorByQC;   //!
   TBranch        *b_KFTau_nKinTaus;   //!
   TBranch        *b_KFTau_TauVis_p4;   //!
   TBranch        *b_KFTau_TauFit_p4;   //!
   TBranch        *b_KFTau_Neutrino_p4;   //!
   TBranch        *b_KFTau_MatchedHPS_idx;   //!
   TBranch        *b_KFTau_Track_idx;   //!
   TBranch        *b_KFTau_indexOfFitInfo;   //!
   TBranch        *b_KFTau_Fit_TauPrimVtx;   //!
   TBranch        *b_KFTau_Fit_IndexToPrimVertexVector;   //!
   TBranch        *b_KFTau_Fit_chi2;   //!
   TBranch        *b_KFTau_Fit_ndf;   //!
   TBranch        *b_KFTau_Fit_ambiguity;   //!
   TBranch        *b_KFTau_Fit_charge;   //!
   TBranch        *b_KFTau_Fit_csum;   //!
   TBranch        *b_KFTau_Fit_iterations;   //!
   TBranch        *b_PFJet_p4;   //!
   TBranch        *b_PFJet_chargedEmEnergy;   //!
   TBranch        *b_PFJet_chargedHadronEnergy;   //!
   TBranch        *b_PFJet_chargedHadronMultiplicity;   //!
   TBranch        *b_PFJet_chargedMuEnergy;   //!
   TBranch        *b_PFJet_chargedMultiplicity;   //!
   TBranch        *b_PFJet_electronEnergy;   //!
   TBranch        *b_PFJet_electronMultiplicity;   //!
   TBranch        *b_PFJet_HFEMEnergy;   //!
   TBranch        *b_PFJet_HFEMMultiplicity;   //!
   TBranch        *b_PFJet_HFHadronEnergy;   //!
   TBranch        *b_PFJet_HFHadronMultiplicity;   //!
   TBranch        *b_PFJet_muonEnergy;   //!
   TBranch        *b_PFJet_muonMultiplicity;   //!
   TBranch        *b_PFJet_neutralEmEnergy;   //!
   TBranch        *b_PFJet_neutralHadronEnergy;   //!
   TBranch        *b_PFJet_neutralHadronMultiplicity;   //!
   TBranch        *b_PFJet_photonEnergy;   //!
   TBranch        *b_PFJet_photonMultiplicity;   //!
   TBranch        *b_PFJet_jetArea;   //!
   TBranch        *b_PFJet_maxDistance;   //!
   TBranch        *b_PFJet_nConstituents;   //!
   TBranch        *b_PFJet_pileup;   //!
   TBranch        *b_PFJet_etaetaMoment;   //!
   TBranch        *b_PFJet_etaphiMoment;   //!
   TBranch        *b_PFJet_Track_idx;   //!
   TBranch        *b_PFJet_MatchedHPS_idx;   //!
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_sumET;   //!
   TBranch        *b_Event_EventNumber;   //!
   TBranch        *b_Event_RunNumber;   //!
   TBranch        *b_Event_bunchCrossing;   //!
   TBranch        *b_Event_orbitNumber;   //!
   TBranch        *b_Event_luminosityBlock;   //!
   TBranch        *b_Event_isRealData;   //!
   TBranch        *b_Track_p4;   //!
   TBranch        *b_Track_Poca;   //!
   TBranch        *b_Track_charge;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_Track_ndof;   //!
   TBranch        *b_Track_numberOfLostHits;   //!
   TBranch        *b_Track_numberOfValidHits;   //!
   TBranch        *b_Track_qualityMask;   //!
   TBranch        *b_MCSignalParticle_p4;   //!
   TBranch        *b_MCSignalParticle_pdgid;   //!
   TBranch        *b_MCSignalParticle_charge;   //!
   TBranch        *b_MCSignalParticle_Poca;   //!
   TBranch        *b_MCSignalParticle_Tauidx;   //!
   TBranch        *b_MCTauandProd_p4;   //!
   TBranch        *b_MCTauandProd_pdgid;   //!
   TBranch        *b_MCTauandProd_midx;   //!
   TBranch        *b_MCTauandProd_charge;   //!
   TBranch        *b_MCTau_JAK;   //!
   TBranch        *b_MCTau_DecayBitMask;   //!

   NtupleReader(TTree *tree=0);
   virtual ~NtupleReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleReader_cxx
NtupleReader::NtupleReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("t");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("t","");
      chain->Add("/net/scratch_cms/institut_3b/cherepanov/MC_DY_SkimmedTauNtuple.root/t");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

NtupleReader::~NtupleReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupleReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupleReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtupleReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Vtx_chi2 = 0;
   Vtx_nTrk = 0;
   Vtx_ndof = 0;
   Vtx_x = 0;
   Vtx_y = 0;
   Vtx_z = 0;
   Vtx_Cov = 0;
   Vtx_Track_idx = 0;
   Muon_p4 = 0;
   Muon_Poca = 0;
   Muon_isGlobalMuon = 0;
   Muon_isStandAloneMuon = 0;
   Muon_isTrackerMuon = 0;
   Muon_isCaloMuon = 0;
   Muon_isIsolationValid = 0;
   Muon_isQualityValid = 0;
   Muon_isTimeValid = 0;
   Muon_emEt03 = 0;
   Muon_emVetoEt03 = 0;
   Muon_hadEt03 = 0;
   Muon_hadVetoEt03 = 0;
   Muon_nJets03 = 0;
   Muon_nTracks03 = 0;
   Muon_sumPt03 = 0;
   Muon_trackerVetoPt03 = 0;
   Muon_emEt05 = 0;
   Muon_emVetoEt05 = 0;
   Muon_hadEt05 = 0;
   Muon_hadVetoEt05 = 0;
   Muon_nJets05 = 0;
   Muon_nTracks05 = 0;
   Muon_sumPt05 = 0;
   Muon_trackerVetoPt05 = 0;
   Muon_Track_idx = 0;
   PFTau_p4 = 0;
   PFTau_isTightIsolation = 0;
   PFTau_isMediumIsolation = 0;
   PFTau_isLooseIsolation = 0;
   PFTau_hpsDecayMode = 0;
   PFTau_Charge = 0;
   PFTau_Track_idx = 0;
   KFTau_discriminatorByKFit = 0;
   KFTau_discriminatorByQC = 0;
   KFTau_TauVis_p4 = 0;
   KFTau_TauFit_p4 = 0;
   KFTau_Neutrino_p4 = 0;
   KFTau_MatchedHPS_idx = 0;
   KFTau_Track_idx = 0;
   KFTau_indexOfFitInfo = 0;
   KFTau_Fit_TauPrimVtx = 0;
   KFTau_Fit_IndexToPrimVertexVector = 0;
   KFTau_Fit_chi2 = 0;
   KFTau_Fit_ndf = 0;
   KFTau_Fit_ambiguity = 0;
   KFTau_Fit_charge = 0;
   KFTau_Fit_csum = 0;
   KFTau_Fit_iterations = 0;
   PFJet_p4 = 0;
   PFJet_chargedEmEnergy = 0;
   PFJet_chargedHadronEnergy = 0;
   PFJet_chargedHadronMultiplicity = 0;
   PFJet_chargedMuEnergy = 0;
   PFJet_chargedMultiplicity = 0;
   PFJet_electronEnergy = 0;
   PFJet_electronMultiplicity = 0;
   PFJet_HFEMEnergy = 0;
   PFJet_HFEMMultiplicity = 0;
   PFJet_HFHadronEnergy = 0;
   PFJet_HFHadronMultiplicity = 0;
   PFJet_muonEnergy = 0;
   PFJet_muonMultiplicity = 0;
   PFJet_neutralEmEnergy = 0;
   PFJet_neutralHadronEnergy = 0;
   PFJet_neutralHadronMultiplicity = 0;
   PFJet_photonEnergy = 0;
   PFJet_photonMultiplicity = 0;
   PFJet_jetArea = 0;
   PFJet_maxDistance = 0;
   PFJet_nConstituents = 0;
   PFJet_pileup = 0;
   PFJet_etaetaMoment = 0;
   PFJet_etaphiMoment = 0;
   PFJet_Track_idx = 0;
   PFJet_MatchedHPS_idx = 0;
   Track_p4 = 0;
   Track_Poca = 0;
   Track_charge = 0;
   Track_chi2 = 0;
   Track_ndof = 0;
   Track_numberOfLostHits = 0;
   Track_numberOfValidHits = 0;
   Track_qualityMask = 0;
   MCSignalParticle_p4 = 0;
   MCSignalParticle_pdgid = 0;
   MCSignalParticle_charge = 0;
   MCSignalParticle_Poca = 0;
   MCSignalParticle_Tauidx = 0;
   MCTauandProd_p4 = 0;
   MCTauandProd_pdgid = 0;
   MCTauandProd_midx = 0;
   MCTauandProd_charge = 0;
   MCTau_JAK = 0;
   MCTau_DecayBitMask = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Vtx_chi2", &Vtx_chi2, &b_Vtx_chi2);
   fChain->SetBranchAddress("Vtx_nTrk", &Vtx_nTrk, &b_Vtx_nTrk);
   fChain->SetBranchAddress("Vtx_ndof", &Vtx_ndof, &b_Vtx_ndof);
   fChain->SetBranchAddress("Vtx_x", &Vtx_x, &b_Vtx_x);
   fChain->SetBranchAddress("Vtx_y", &Vtx_y, &b_Vtx_y);
   fChain->SetBranchAddress("Vtx_z", &Vtx_z, &b_Vtx_z);
   fChain->SetBranchAddress("Vtx_Cov", &Vtx_Cov, &b_Vtx_Cov);
   fChain->SetBranchAddress("Vtx_Track_idx", &Vtx_Track_idx, &b_Vtx_Track_idx);
   fChain->SetBranchAddress("Muon_p4", &Muon_p4, &b_Muon_p4);
   fChain->SetBranchAddress("Muon_Poca", &Muon_Poca, &b_Muon_Poca);
   fChain->SetBranchAddress("Muon_isGlobalMuon", &Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", &Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isCaloMuon", &Muon_isCaloMuon, &b_Muon_isCaloMuon);
   fChain->SetBranchAddress("Muon_isIsolationValid", &Muon_isIsolationValid, &b_Muon_isIsolationValid);
   fChain->SetBranchAddress("Muon_isQualityValid", &Muon_isQualityValid, &b_Muon_isQualityValid);
   fChain->SetBranchAddress("Muon_isTimeValid", &Muon_isTimeValid, &b_Muon_isTimeValid);
   fChain->SetBranchAddress("Muon_emEt03", &Muon_emEt03, &b_Muon_emEt03);
   fChain->SetBranchAddress("Muon_emVetoEt03", &Muon_emVetoEt03, &b_Muon_emVetoEt03);
   fChain->SetBranchAddress("Muon_hadEt03", &Muon_hadEt03, &b_Muon_hadEt03);
   fChain->SetBranchAddress("Muon_hadVetoEt03", &Muon_hadVetoEt03, &b_Muon_hadVetoEt03);
   fChain->SetBranchAddress("Muon_nJets03", &Muon_nJets03, &b_Muon_nJets03);
   fChain->SetBranchAddress("Muon_nTracks03", &Muon_nTracks03, &b_Muon_nTracks03);
   fChain->SetBranchAddress("Muon_sumPt03", &Muon_sumPt03, &b_Muon_sumPt03);
   fChain->SetBranchAddress("Muon_trackerVetoPt03", &Muon_trackerVetoPt03, &b_Muon_trackerVetoPt03);
   fChain->SetBranchAddress("Muon_emEt05", &Muon_emEt05, &b_Muon_emEt05);
   fChain->SetBranchAddress("Muon_emVetoEt05", &Muon_emVetoEt05, &b_Muon_emVetoEt05);
   fChain->SetBranchAddress("Muon_hadEt05", &Muon_hadEt05, &b_Muon_hadEt05);
   fChain->SetBranchAddress("Muon_hadVetoEt05", &Muon_hadVetoEt05, &b_Muon_hadVetoEt05);
   fChain->SetBranchAddress("Muon_nJets05", &Muon_nJets05, &b_Muon_nJets05);
   fChain->SetBranchAddress("Muon_nTracks05", &Muon_nTracks05, &b_Muon_nTracks05);
   fChain->SetBranchAddress("Muon_sumPt05", &Muon_sumPt05, &b_Muon_sumPt05);
   fChain->SetBranchAddress("Muon_trackerVetoPt05", &Muon_trackerVetoPt05, &b_Muon_trackerVetoPt05);
   fChain->SetBranchAddress("Muon_Track_idx", &Muon_Track_idx, &b_Muon_Track_idx);
   fChain->SetBranchAddress("PFTau_p4", &PFTau_p4, &b_PFTau_p4);
   fChain->SetBranchAddress("PFTau_isTightIsolation", &PFTau_isTightIsolation, &b_PFTau_isTightIsolation);
   fChain->SetBranchAddress("PFTau_isMediumIsolation", &PFTau_isMediumIsolation, &b_PFTau_isMediumIsolation);
   fChain->SetBranchAddress("PFTau_isLooseIsolation", &PFTau_isLooseIsolation, &b_PFTau_isLooseIsolation);
   fChain->SetBranchAddress("PFTau_hpsDecayMode", &PFTau_hpsDecayMode, &b_PFTau_hpsDecayMode);
   fChain->SetBranchAddress("PFTau_Charge", &PFTau_Charge, &b_PFTau_Charge);
   fChain->SetBranchAddress("PFTau_Track_idx", &PFTau_Track_idx, &b_PFTau_Track_idx);
   fChain->SetBranchAddress("KFTau_discriminatorByKFit", &KFTau_discriminatorByKFit, &b_KFTau_discriminatorByKFit);
   fChain->SetBranchAddress("KFTau_discriminatorByQC", &KFTau_discriminatorByQC, &b_KFTau_discriminatorByQC);
   fChain->SetBranchAddress("KFTau_nKinTaus", &KFTau_nKinTaus, &b_KFTau_nKinTaus);
   fChain->SetBranchAddress("KFTau_TauVis_p4", &KFTau_TauVis_p4, &b_KFTau_TauVis_p4);
   fChain->SetBranchAddress("KFTau_TauFit_p4", &KFTau_TauFit_p4, &b_KFTau_TauFit_p4);
   fChain->SetBranchAddress("KFTau_Neutrino_p4", &KFTau_Neutrino_p4, &b_KFTau_Neutrino_p4);
   fChain->SetBranchAddress("KFTau_MatchedHPS_idx", &KFTau_MatchedHPS_idx, &b_KFTau_MatchedHPS_idx);
   fChain->SetBranchAddress("KFTau_Track_idx", &KFTau_Track_idx, &b_KFTau_Track_idx);
   fChain->SetBranchAddress("KFTau_indexOfFitInfo", &KFTau_indexOfFitInfo, &b_KFTau_indexOfFitInfo);
   fChain->SetBranchAddress("KFTau_Fit_TauPrimVtx", &KFTau_Fit_TauPrimVtx, &b_KFTau_Fit_TauPrimVtx);
   fChain->SetBranchAddress("KFTau_Fit_IndexToPrimVertexVector", &KFTau_Fit_IndexToPrimVertexVector, &b_KFTau_Fit_IndexToPrimVertexVector);
   fChain->SetBranchAddress("KFTau_Fit_chi2", &KFTau_Fit_chi2, &b_KFTau_Fit_chi2);
   fChain->SetBranchAddress("KFTau_Fit_ndf", &KFTau_Fit_ndf, &b_KFTau_Fit_ndf);
   fChain->SetBranchAddress("KFTau_Fit_ambiguity", &KFTau_Fit_ambiguity, &b_KFTau_Fit_ambiguity);
   fChain->SetBranchAddress("KFTau_Fit_charge", &KFTau_Fit_charge, &b_KFTau_Fit_charge);
   fChain->SetBranchAddress("KFTau_Fit_csum", &KFTau_Fit_csum, &b_KFTau_Fit_csum);
   fChain->SetBranchAddress("KFTau_Fit_iterations", &KFTau_Fit_iterations, &b_KFTau_Fit_iterations);
   fChain->SetBranchAddress("PFJet_p4", &PFJet_p4, &b_PFJet_p4);
   fChain->SetBranchAddress("PFJet_chargedEmEnergy", &PFJet_chargedEmEnergy, &b_PFJet_chargedEmEnergy);
   fChain->SetBranchAddress("PFJet_chargedHadronEnergy", &PFJet_chargedHadronEnergy, &b_PFJet_chargedHadronEnergy);
   fChain->SetBranchAddress("PFJet_chargedHadronMultiplicity", &PFJet_chargedHadronMultiplicity, &b_PFJet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_chargedMuEnergy", &PFJet_chargedMuEnergy, &b_PFJet_chargedMuEnergy);
   fChain->SetBranchAddress("PFJet_chargedMultiplicity", &PFJet_chargedMultiplicity, &b_PFJet_chargedMultiplicity);
   fChain->SetBranchAddress("PFJet_electronEnergy", &PFJet_electronEnergy, &b_PFJet_electronEnergy);
   fChain->SetBranchAddress("PFJet_electronMultiplicity", &PFJet_electronMultiplicity, &b_PFJet_electronMultiplicity);
   fChain->SetBranchAddress("PFJet_HFEMEnergy", &PFJet_HFEMEnergy, &b_PFJet_HFEMEnergy);
   fChain->SetBranchAddress("PFJet_HFEMMultiplicity", &PFJet_HFEMMultiplicity, &b_PFJet_HFEMMultiplicity);
   fChain->SetBranchAddress("PFJet_HFHadronEnergy", &PFJet_HFHadronEnergy, &b_PFJet_HFHadronEnergy);
   fChain->SetBranchAddress("PFJet_HFHadronMultiplicity", &PFJet_HFHadronMultiplicity, &b_PFJet_HFHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_muonEnergy", &PFJet_muonEnergy, &b_PFJet_muonEnergy);
   fChain->SetBranchAddress("PFJet_muonMultiplicity", &PFJet_muonMultiplicity, &b_PFJet_muonMultiplicity);
   fChain->SetBranchAddress("PFJet_neutralEmEnergy", &PFJet_neutralEmEnergy, &b_PFJet_neutralEmEnergy);
   fChain->SetBranchAddress("PFJet_neutralHadronEnergy", &PFJet_neutralHadronEnergy, &b_PFJet_neutralHadronEnergy);
   fChain->SetBranchAddress("PFJet_neutralHadronMultiplicity", &PFJet_neutralHadronMultiplicity, &b_PFJet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("PFJet_photonEnergy", &PFJet_photonEnergy, &b_PFJet_photonEnergy);
   fChain->SetBranchAddress("PFJet_photonMultiplicity", &PFJet_photonMultiplicity, &b_PFJet_photonMultiplicity);
   fChain->SetBranchAddress("PFJet_jetArea", &PFJet_jetArea, &b_PFJet_jetArea);
   fChain->SetBranchAddress("PFJet_maxDistance", &PFJet_maxDistance, &b_PFJet_maxDistance);
   fChain->SetBranchAddress("PFJet_nConstituents", &PFJet_nConstituents, &b_PFJet_nConstituents);
   fChain->SetBranchAddress("PFJet_pileup", &PFJet_pileup, &b_PFJet_pileup);
   fChain->SetBranchAddress("PFJet_etaetaMoment", &PFJet_etaetaMoment, &b_PFJet_etaetaMoment);
   fChain->SetBranchAddress("PFJet_etaphiMoment", &PFJet_etaphiMoment, &b_PFJet_etaphiMoment);
   fChain->SetBranchAddress("PFJet_Track_idx", &PFJet_Track_idx, &b_PFJet_Track_idx);
   fChain->SetBranchAddress("PFJet_MatchedHPS_idx", &PFJet_MatchedHPS_idx, &b_PFJet_MatchedHPS_idx);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumET", &MET_sumET, &b_MET_sumET);
   fChain->SetBranchAddress("Event_EventNumber", &Event_EventNumber, &b_Event_EventNumber);
   fChain->SetBranchAddress("Event_RunNumber", &Event_RunNumber, &b_Event_RunNumber);
   fChain->SetBranchAddress("Event_bunchCrossing", &Event_bunchCrossing, &b_Event_bunchCrossing);
   fChain->SetBranchAddress("Event_orbitNumber", &Event_orbitNumber, &b_Event_orbitNumber);
   fChain->SetBranchAddress("Event_luminosityBlock", &Event_luminosityBlock, &b_Event_luminosityBlock);
   fChain->SetBranchAddress("Event_isRealData", &Event_isRealData, &b_Event_isRealData);
   fChain->SetBranchAddress("Track_p4", &Track_p4, &b_Track_p4);
   fChain->SetBranchAddress("Track_Poca", &Track_Poca, &b_Track_Poca);
   fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("Track_ndof", &Track_ndof, &b_Track_ndof);
   fChain->SetBranchAddress("Track_numberOfLostHits", &Track_numberOfLostHits, &b_Track_numberOfLostHits);
   fChain->SetBranchAddress("Track_numberOfValidHits", &Track_numberOfValidHits, &b_Track_numberOfValidHits);
   fChain->SetBranchAddress("Track_qualityMask", &Track_qualityMask, &b_Track_qualityMask);
   fChain->SetBranchAddress("MCSignalParticle_p4", &MCSignalParticle_p4, &b_MCSignalParticle_p4);
   fChain->SetBranchAddress("MCSignalParticle_pdgid", &MCSignalParticle_pdgid, &b_MCSignalParticle_pdgid);
   fChain->SetBranchAddress("MCSignalParticle_charge", &MCSignalParticle_charge, &b_MCSignalParticle_charge);
   fChain->SetBranchAddress("MCSignalParticle_Poca", &MCSignalParticle_Poca, &b_MCSignalParticle_Poca);
   fChain->SetBranchAddress("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx, &b_MCSignalParticle_Tauidx);
   fChain->SetBranchAddress("MCTauandProd_p4", &MCTauandProd_p4, &b_MCTauandProd_p4);
   fChain->SetBranchAddress("MCTauandProd_pdgid", &MCTauandProd_pdgid, &b_MCTauandProd_pdgid);
   fChain->SetBranchAddress("MCTauandProd_midx", &MCTauandProd_midx, &b_MCTauandProd_midx);
   fChain->SetBranchAddress("MCTauandProd_charge", &MCTauandProd_charge, &b_MCTauandProd_charge);
   fChain->SetBranchAddress("MCTau_JAK", &MCTau_JAK, &b_MCTau_JAK);
   fChain->SetBranchAddress("MCTau_DecayBitMask", &MCTau_DecayBitMask, &b_MCTau_DecayBitMask);
   Notify();
}

Bool_t NtupleReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupleReader_cxx
