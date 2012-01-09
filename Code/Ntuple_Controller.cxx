//Ntuple_Controller.cxx IMPLEMENTATION FILE

#ifndef Ntuple_Controller_cxx
#define Ntuple_Controller_cxx

#include "Ntuple_Controller.h"
#include "Tools.h"
#include "PDG_Var.h"

// External code
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,ObjEvent(-1)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("t");
  std::cout << "Ntuple_Controller" << RootFiles.size() << std::endl;
  for(int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;  
  if(chain==0){
    std::cout << "failed" << std::endl;
  }
  std::cout << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  std::cout << "Ntuple Configured" << std::endl;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  std::cout << Ntp->fChain->GetEntries() << std::endl;
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Int_t ientry = Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  std::cout << "Ntuple_Controller::~Ntuple_Controller()" << std::endl;
  delete Ntp;
  std::cout << "Ntuple_Controller::~Ntuple_Controller() complete" << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
    std::cout << "Starting D3PD cloning " << std::endl;  
    newfile = new TFile(n+".root","recreate");
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  std::cout << "Ntuple_Controller::SaveCloneTree()"<< std::endl;
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  std::cout << "Ntuple_Controller::SaveCloneTree() Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  std::cout << "Ntuple_Controller::ThinTree" << std::endl;

  std::cout << "Ntuple_Controller::ThinTree complete" << std::endl;  
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}


//Physics get Functions
int Ntuple_Controller::GetMCID(){
  if((Ntp->DataMC_Type)==10230530){
    for(int i=0;i<NMCSignalParticles();i++){
      if(abs(MCSignalParticle_pdgid(i))==PdtPdgMini::Z0){
	if(fabs(MCSignalParticle_p4(i).M()-PDG_Var::Z_mass())<3*PDG_Var::Z_width()){
	  return Signal;
	}
      }
    }
    return Ntp->DataMC_Type;
  }
  if(HConfig.hasID(Ntp->DataMC_Type))return Ntp->DataMC_Type;
  return (Ntp->DataMC_Type%100); 
}

TMatrixF     Ntuple_Controller::Vtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      if(j*dim+k<Ntp->Vtx_Cov->at(i).size()){
	M[j][k]=Ntp->Vtx_Cov->at(i).at(j).at(k);
	M[k][j]=Ntp->Vtx_Cov->at(i).at(j).at(k);
      }
    }
  }
}

bool Ntuple_Controller::isVtxGood(unsigned int i){
  if(0<=i && i<NVtx()){
    if(Vtx_Track_idx(i).size()>4)return true;
  }
  return false;
}

bool Ntuple_Controller::isGoodMuon(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ΔR(μ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isGoodMuon_nooverlapremoval(i)){
    for(unsigned int j=0;j<NPFJets();j++){
      if(isGoodJet_nooverlapremoval(j)){
	if(Tools::dr(Muons_p4(i),PFJet_p4(j))<0.3) return false;
      }
    }
    return true;
  }
  return false;
}


bool Ntuple_Controller::isGoodMuon_nooverlapremoval(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger and no overlpar removal applied
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
  //  GlobalMuon && TrackerMuon applied applied                                         
  //  pT > 20 GeV                                                                       
  //  fasb(eta) < 2.4                                                                   
  //  normChi2 < 10                                                                     
  //  n Tracker Hits > 10                                                               
  //  n Muon Hit > 0                                                                    
  //  Transverse IP of the muon wrt the beam spot (cm) < 0.02                           
  //  relIso < 0.20                                                                      
  //  fabs( muon.vertex().z() - PV.z()) not applied                           
  //  muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() not applied
  //  numberOfMatchedStations() not applied                              
  if(Muon_isGlobalMuon(i) && Muon_isStandAloneMuon(i)){
    if(Muons_p4(i).Pt()>17.0){
      if(fabs(Muons_p4(i).Eta())<2.4){
	if(Muon_normChi2(i)<10.0){
	  if(Muon_innerTrack_numberofValidHits(i)>10){
	    if(Muon_hitPattern_numberOfValidMuonHits(i)>0){
	      if((Muon_emEt03(i)+Muon_hadEt03(i)+Muon_sumPt03(i))/Muons_p4(i).Pt()<0.2){
		return true;
	      }
	    }
	  }
	}
      }
    }
  }
  return false;
}





bool Ntuple_Controller::isGoodJet(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel 
  //  isGoodJet_nooverlapremoval(i) with:
  //  deltaR jet-electron cleaning < 0.4 (2 selected lepton only) < 0.3 (e+jets only) 0.3
  //  deltaR jet-muon cleaning < 0.4(2 selected lepton only) < 0.3 (mu+jets only), 0.1 for PF and JET 0.3
  if(isGoodJet_nooverlapremoval(i)){
    if(isGoodMuon_nooverlapremoval(i)){
      for(unsigned int j=0;j<NMuons();j++){
	if(isGoodMuon_nooverlapremoval(j)){
	  if(Tools::dr(Muons_p4(j),PFJet_p4(i))<0.4) return false;
	}
      }
      return true;
    }

  }
}

bool Ntuple_Controller::isGoodJet_nooverlapremoval(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel
  //  Jet ID defined in isJetID(i)               
  //  cut         dilepton    l+jet
  //  corrected pT 30 GeV 30 GeV    
  //  residual correction (data) applied applied
  //  abs(eta) < 2.5 < 2.4                      
  //  jet ID applied applied                    
  if(isJetID(i)){
    if(PFJet_p4(i).Pt()>15.0){
      if(fabs(PFJet_p4(i).Eta())<2.4){
	return true;
      }
    }
  }
  return false;
}

bool Ntuple_Controller::isJetID(unsigned int i){
  //  Top Dilepton Jet selection with pt and iso matching the muon and tau.
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
  //  Jet ID :
  //  number of constituents>1 (patJet->numberOfDaughters())
  //  CEF<0.99 (patJet->chargedEmEnergyFraction())
  //  NHF<0.99 (patJet->neutralHadronEnergyFraction())
  //  NEF<0.99 (patJet->neutralEmEnergyFraction())
  //  if |η|<2.4, CHF>0 (patJet->chargedHadronEnergyFraction())
  //  if |η|<2.4, NCH>0 (patJet->chargedMultiplicity()) 
  /////////////////////////////////////////////////////////////////////////
  // apply jet ID
  bool JetID_ok=false;
  if(PFJet_numberOfDaughters(i)>1){
    if(PFJet_chargedEmEnergyFraction(i)<0.99){
      if(PFJet_neutralHadronEnergyFraction(i)<0.99){
	if(PFJet_PFJet_neutralEmEnergyFraction(i)<0.99){
	  if(fabs(PFJet_p4(i).Eta())<2.4){
	    if(PFJet_chargedHadronEnergyFraction(i)>0){
	      if(PFJet_chargedMultiplicity(i)>0){
		return true;
	      }
	    }
	  }
	  else{
	    return true;
	  }
	}
      }
    }
  }
  return false;
}


TMatrixF     Ntuple_Controller::Track_parCov(unsigned int i){
  unsigned int dim=5;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      if(j*dim+k<Ntp->Track_parCov->at(i).size()){
        M[j][k]=Ntp->Track_parCov->at(i).at(j).at(k);
        M[k][j]=Ntp->Track_parCov->at(i).at(j).at(k);
      }
    }
  }
}

float     Ntuple_Controller::Track_parCov(unsigned int i, TrackPar par1, TrackPar par2){
  if(par1>par2)return Ntp->Track_parCov->at(i).at(par1).at(par2);
  return Ntp->Track_parCov->at(i).at(par2).at(par1);
}

#endif



