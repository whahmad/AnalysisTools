//Ntuple_Controller.cxx IMPLEMENTATION FILE


#include "Ntuple_Controller.h"
#include "Tools.h"
#include "PDG_Var.h"

// External code
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,ObjEvent(-1)
  ,verbose(false)
  ,TauSpinerInt()
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

  if((Ntp->DataMC_Type)==DataMCType::DY_ll_Signal && HistoC.hasID(DataMCType::DY_ll_Signal)){
    for(int i=0;i<NMCSignalParticles();i++){
      if(abs(MCSignalParticle_pdgid(i))==PdtPdgMini::Z0){
	if(fabs(MCSignalParticle_p4(i).M()-PDG_Var::Z_mass())<3*PDG_Var::Z_width()){
	  return DataMCType::Signal;
	}
      }
    }
    return Ntp->DataMC_Type;
  }
  if(Ntp->DataMC_Type>100){
    if(HistoC.hasID(Ntp->DataMC_Type%100)){
      //std::cout << "MODULO OPERATION RETURNS: " << Ntp->DataMC_Type%100 << std::endl;
      return Ntp->DataMC_Type%100;
    }
  }
  if(HConfig.hasID(Ntp->DataMC_Type)){
	  //std::cout << "Ntp->DataMC_Type: " << Ntp->DataMC_Type << std::endl;
	  return Ntp->DataMC_Type;
  }
}

TMatrixF     Ntuple_Controller::Vtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      M[j][k]=Ntp->Vtx_Cov->at(i).at(j).at(k);
      M[k][j]=Ntp->Vtx_Cov->at(i).at(j).at(k);
    }
  }
  return M;
}

bool Ntuple_Controller::isVtxGood(unsigned int i){
  if(0<=i && i<NVtx()){
    if(Vtx_Track_idx(i).size()>4)return true;
  }
  return false;
}

/*int Ntuple_Controller::Vtx_Corr_Track_idx(unsigned int i){
	for(unsigned j=0;j<Vtx_nTrk(i);j++){
		for(unsigned k=0;k<NTracks();k++){
			if(Vtx_TracksP4(i,j).DeltaR(Track_p4(k))<0.00001){
				return k;
			}
		}
	}
	return -1;
}*/

std::vector<int> Ntuple_Controller::Vtx_Corr_Track_idx(unsigned int i){
	std::vector<int> track_idx;
	track_idx.clear();
	for(unsigned j=0;j<Vtx_nTrk(i);j++){
		for(unsigned k=0;k<NTracks();k++){
			if(Vtx_TracksP4(i,j).DeltaR(Track_p4(k))<0.00001){
				track_idx.push_back(k);
			}
		}
	}
	return track_idx;
}

bool Ntuple_Controller::isGoodMuon(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ?R(µ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isTightMuon(i)){
    unsigned int jet_idx=0;
    return !muonhasJetOverlap(i,jet_idx);
  }
  return false;
}

bool Ntuple_Controller::isGoodMuon(unsigned int i, unsigned int j){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ?R(µ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isTightMuon(i,j)){
    unsigned int jet_idx=0;
    //return !muonhasJetOverlap(i,jet_idx);
    if(!muonhasJetOverlap(i,jet_idx) && !muonhasJetMatch(i,jet_idx)){
    	//std::cout << "MUON PASSES ID" << std::endl;
    	return true;
    }else{
    	return false;
    }
  }
  return false;
}

bool Ntuple_Controller::muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Muons_p4(muon_idx),PFJet_p4(j))>0.2 && Tools::dr(Muons_p4(muon_idx),PFJet_p4(j))<0.4){ jet_idx=j;return true;}
    }
  }
  return false;
}

bool Ntuple_Controller::muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Muons_p4(muon_idx),PFJet_p4(j))<0.2){ jet_idx=j;return true;}
    }
  }
  return false;
}



bool Ntuple_Controller::isTightMuon(unsigned int i, unsigned int j){
  // Tight muon ID from Muon POG
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
  // Criteria:
  // Muon is GlobalMuon
  // Muon is PFMuon
  // normChi2 < 10
  // number of Muon chamber hits > 0
  // number of matched Muon stations > 1
  // dxy < 0.2 cm, dxy = [-(vx-vtx_x)*py + (vy-vtx_y)*px]/pt
  // dz < 0.5 cm, dz = (vz-vtx_z) - [(vx-vtx_x)*px+(vy-vtx_y)*py]/pt * pz/pt
  // number of pixel hits > 0
  // number of tracker layers with hits > 5
  if(Muon_isGlobalMuon(i) &&
		  Muon_isPFMuon(i) &&
		  Muon_normChi2(i)<10.0 &&
		  Muon_hitPattern_numberOfValidMuonHits(i)>0 &&
		  Muon_numberOfMatchedStations(i)>1 &&
		  fabs((-(Muon_Poca(i).X()-Vtx(j).X())*Muons_p4(i).Py()+(Muon_Poca(i).Y()-Vtx(j).Y())*Muons_p4(i).Px())/Muons_p4(i).Pt())<0.2 &&
		  fabs(Muon_Poca(i).Z()-Vtx(j).Z()-((Muon_Poca(i).X()-Vtx(j).X())*Muons_p4(i).Px()+(Muon_Poca(i).Y()-Vtx(j).Y())*Muons_p4(i).Py())*Muons_p4(i).Pz()/pow(Muons_p4(i).Pt(),2))<0.5 &&
		  Muon_numberofValidPixelHits(i)>0 &&
		  Muon_trackerLayersWithMeasurement(i)>5
		  ){
	  return true;
  }
  return false;
}


bool Ntuple_Controller::isTightMuon(unsigned int i){
  // Tight muon ID from Muon POG
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
  // Criteria:
  // Muon is GlobalMuon
  // Muon is PFMuon
  // normChi2 < 10
  // number of Muon chamber hits > 0
  // number of matched Muon stations > 1
  // dxy < 0.2 cm, dxy = [-(vx-vtx_x)*py + (vy-vtx_y)*px]/pt
  // dz < 0.5 cm, dz = (vz-vtx_z) - [(vx-vtx_x)*px+(vy-vtx_y)*py]/pt * pz/pt
  // number of pixel hits > 0
  // number of tracker layers with hits > 5
  if(Muon_isGlobalMuon(i) &&
		  Muon_isPFMuon(i) &&
		  Muon_normChi2(i)<10.0 &&
		  Muon_hitPattern_numberOfValidMuonHits(i)>0 &&
		  Muon_numberOfMatchedStations(i)>1 &&
		  fabs((-(Muon_Poca(i).X()-Vtx(0).X())*Muons_p4(i).Py()+(Muon_Poca(i).Y()-Vtx(0).Y())*Muons_p4(i).Px())/Muons_p4(i).Pt())<0.2 &&
		  fabs(Muon_Poca(i).Z()-Vtx(0).Z()-((Muon_Poca(i).X()-Vtx(0).X())*Muons_p4(i).Px()+(Muon_Poca(i).Y()-Vtx(0).Y())*Muons_p4(i).Py())*Muons_p4(i).Pz()/pow(Muons_p4(i).Pt(),2))<0.5 &&
		  Muon_numberofValidPixelHits(i)>0 &&
		  Muon_trackerLayersWithMeasurement(i)>5
		  ){
	  return true;
  }
  return false;
}


bool Ntuple_Controller::isLooseMuon(unsigned int i){
	// Loose muon ID from Muon POG
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
	// Criteria:
	// Muon is GlobalMuon || Muon is TrackerMuon
	// Muon is PFMuon
	if(Muon_isPFMuon(i) &&
			(Muon_isGlobalMuon(i) || Muon_isTrackerMuon(i))
			){
		return true;
	}
	return false;
}


bool Ntuple_Controller::isSoftMuon(unsigned int i){
	// Loose muon ID from Muon POG
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
	// Criteria:
	// number of matched Muon stations > 1
	// Number of tracker layers with hits > 5
	// Number of pixel layers > 1
	// normChi2 < 1.8
	// dxy < 3.0 cm
	// dz < 30.0 cm
	if(Muon_numberOfMatchedStations(i)>1 &&
			Muon_trackerLayersWithMeasurement(i)>5 &&
			Muon_hitPattern_pixelLayerwithMeas(i)>1 &&
			Muon_normChi2(i)<1.8 &&
			fabs((-(Muon_Poca(i).X()-Vtx(0).X())*Muons_p4(i).Py()+(Muon_Poca(i).Y()-Vtx(0).Y())*Muons_p4(i).Px())/Muons_p4(i).Pt())<3.0 &&
			fabs(Muon_Poca(i).Z()-Vtx(0).Z()-((Muon_Poca(i).X()-Vtx(0).X())*Muons_p4(i).Px()+(Muon_Poca(i).Y()-Vtx(0).Y())*Muons_p4(i).Py())*Muons_p4(i).Pz()/pow(Muons_p4(i).Pt(),2))<30.0
			){
		return true;
	}
	return false;
}


bool Ntuple_Controller::isFakeMuon(unsigned int i, unsigned int j){
	if(Muon_isGlobalMuon(i) &&
			Muons_p4(i).Pt()>10 &&
			fabs(Muons_p4(i).Eta())<2.1 &&
			fabs((-(Muon_Poca(i).X()-Vtx(j).X())*Muons_p4(i).Py()+(Muon_Poca(i).Y()-Vtx(j).Y())*Muons_p4(i).Px())/Muons_p4(i).Pt())<0.2
			){
		if(Muons_p4(i).Pt()>20 &&
				Muon_RelIso(i)<0.4
				){
			return true;
		}else if(Muons_p4(i).Pt()<=20 &&
				Muon_sumChargedHadronPt04(i)+std::max(0.,Muon_sumNeutralHadronEt04(i)+Muon_sumPhotonEt04(i)-0.5*Muon_sumPUPt04(i))<8.
				){
			return true;
		}
	}
	return false;
}


float Ntuple_Controller::Muon_RelIso(unsigned int i){
	// Loose relative isolation from Muon POG
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
	// Criteria:
	// sum of Et from charged hadrons, neutral hadrons and photons divided by pt of muon < 0.2
	// cone size < 0.4
	return (Muon_sumChargedHadronPt04(i)+std::max(0.,Muon_sumNeutralHadronEt04(i)+Muon_sumPhotonEt04(i)-0.5*Muon_sumPUPt04(i)))/Muons_p4(i).Pt();
}



bool Ntuple_Controller::isGoodElectron(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ?R(µ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isTightElectron(i)){
    unsigned int jet_idx=0;
    return !electronhasJetOverlap(i,jet_idx);
  }
  return false;
}
bool Ntuple_Controller::isGoodElectron(unsigned int i, unsigned int j){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ?R(µ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isTightElectron(i,j)){
    unsigned int jet_idx=0;
    //return !electronhasJetOverlap(i,jet_idx);
    if(!electronhasJetOverlap(i,jet_idx) && !electronhasJetMatch(i,jet_idx)){
    	//std::cout << "ELECTRON PASSES ID" << std::endl;
    	return true;
    }else{
    	return false;
    }
  }
  return false;
}

bool Ntuple_Controller::electronhasJetOverlap(unsigned int electron_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Electron_p4(electron_idx),PFJet_p4(j))>0.2 && Tools::dr(Electron_p4(electron_idx),PFJet_p4(j))<0.4){ jet_idx=j;return true;}
    }
  }
  return false;
}

bool Ntuple_Controller::electronhasJetMatch(unsigned int electron_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Electron_p4(electron_idx),PFJet_p4(j))<0.2){ jet_idx=j;return true;}
    }
  }
  return false;
}


bool Ntuple_Controller::isTightElectron(unsigned int i){
	// Tight electron ID from EGamma POG for barrel (endcaps)
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Conversion_Rejection
	// Criteria:
	// |eta_supercluster|<1.442 (1.566<|eta_supercluster|<2.5)
	// dEtaIn < 0.004 (0.005)
	// dPhiIn < 0.03 (0.02)
	// sigmaIEtaIEta < 0.01 (0.03) ???????
	// H/E < 0.12 (0.10)
	// d0 < 0.02 (0.02), d0 = -dxy
	// dz < 0.1 (0.1)
	// fabs(1/E-1/p) < 0.05 (0.05)
	// PFIsolation/pT < 0.10 (0.10 pT>20 Gev, 0.07 pT<20GeV)
	// Conversion rejection: vertex fit probability < 1e-6 (1e-6)
	// Conversion rejection: missing hits = 0 (0)
	if(fabs(Electron_supercluster_eta(i))<1.442){ //barrel
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.004 &&
				Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.03 &&
				Electron_sigmaIetaIeta(i)<0.01 &&
				Electron_hadronicOverEm(i)<0.12 &&
				fabs((-(Electron_Poca(i).X()-Vtx(0).X())*Electron_p4(i).Py()+(Electron_Poca(i).Y()-Vtx(0).Y())*Electron_p4(i).Px())/Electron_p4(i).Pt())<0.02 &&
				fabs(Electron_Poca(i).Z()-Vtx(0).Z()-((Electron_Poca(i).X()-Vtx(0).X())*Electron_p4(i).Px()+(Electron_Poca(i).Y()-Vtx(0).Y())*Electron_p4(i).Py())*Electron_p4(i).Pz()/pow(Electron_p4(i).Pt(),2))<0.1 &&
				fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))<0.05 &&
				Electron_RelIso(i)<0.10 &&
				!Electron_HasMatchedConversions(i) &&
				Electron_numberOfMissedHits(i)==0
				){
			return true;
		}
	}else if(fabs(Electron_supercluster_eta(i))>1.566 && fabs(Electron_supercluster_eta(i))<2.5){ //endcaps
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.005 &&
				Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.02 &&
				Electron_sigmaIetaIeta(i)<0.03 &&
				Electron_hadronicOverEm(i)<0.10 &&
				fabs((-(Electron_Poca(i).X()-Vtx(0).X())*Electron_p4(i).Py()+(Electron_Poca(i).Y()-Vtx(0).Y())*Electron_p4(i).Px())/Electron_p4(i).Pt())<0.02 &&
				fabs(Electron_Poca(i).Z()-Vtx(0).Z()-((Electron_Poca(i).X()-Vtx(0).X())*Electron_p4(i).Px()+(Electron_Poca(i).Y()-Vtx(0).Y())*Electron_p4(i).Py())*Electron_p4(i).Pz()/pow(Electron_p4(i).Pt(),2))<0.1 &&
				fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))<0.05 &&
				!Electron_HasMatchedConversions(i) &&
				Electron_numberOfMissedHits(i)==0
				){
			if(Electron_p4(i).Pt()>=20.0 && Electron_RelIso(i)<0.10){
				return true;
			}else if(Electron_p4(i).Pt()<20.0 && Electron_RelIso(i)<0.07){
				return true;
			}
		}
	}
	return false;
}



bool Ntuple_Controller::isTightElectron(unsigned int i, unsigned int j){
	// Tight electron ID from EGamma POG for barrel (endcaps)
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Conversion_Rejection
	// Criteria:
	// |eta_supercluster|<1.442 (1.566<|eta_supercluster|<2.5)
	// dEtaIn < 0.004 (0.005)
	// dPhiIn < 0.03 (0.02)
	// sigmaIEtaIEta < 0.01 (0.03) ???????
	// H/E < 0.12 (0.10)
	// d0 < 0.02 (0.02), d0 = -dxy
	// dz < 0.1 (0.1)
	// fabs(1/E-1/p) < 0.05 (0.05)
	// PFIsolation/pT < 0.10 (0.10 pT>20 Gev, 0.07 pT<20GeV)
	// Conversion rejection: vertex fit probability < 1e-6 (1e-6)
	// Conversion rejection: missing hits = 0 (0)
	if(fabs(Electron_supercluster_eta(i))<1.442){ //barrel
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.004 &&
				Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.03 &&
				Electron_sigmaIetaIeta(i)<0.01 &&
				Electron_hadronicOverEm(i)<0.12 &&
				fabs((-(Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Py()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Px())/Electron_p4(i).Pt())<0.02 &&
				fabs(Electron_Poca(i).Z()-Vtx(j).Z()-((Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Px()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Py())*Electron_p4(i).Pz()/pow(Electron_p4(i).Pt(),2))<0.1 &&
				fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))<0.05 &&
				Electron_RelIso(i)<0.10 &&
				!Electron_HasMatchedConversions(i) &&
				Electron_numberOfMissedHits(i)==0
				){
			return true;
		}
	}else if(fabs(Electron_supercluster_eta(i))>1.566 && fabs(Electron_supercluster_eta(i))<2.5){ //endcaps
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.005 &&
				Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.02 &&
				Electron_sigmaIetaIeta(i)<0.03 &&
				Electron_hadronicOverEm(i)<0.10 &&
				fabs((-(Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Py()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Px())/Electron_p4(i).Pt())<0.02 &&
				fabs(Electron_Poca(i).Z()-Vtx(j).Z()-((Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Px()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Py())*Electron_p4(i).Pz()/pow(Electron_p4(i).Pt(),2))<0.1 &&
				fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))<0.05 &&
				!Electron_HasMatchedConversions(i) &&
				Electron_numberOfMissedHits(i)==0
				){
			if(Electron_p4(i).Pt()>=20.0 && Electron_RelIso(i)<0.10){
				return true;
			}else if(Electron_p4(i).Pt()<20.0 && Electron_RelIso(i)<0.07){
				return true;
			}
		}
	}
	return false;
}


bool Ntuple_Controller::isFakeElectron(unsigned int i, unsigned int j){
	// Fakeable electrons partially from H->tautau analysis (AN2013_011_v6)
	// Criteria:
	// |eta_supercluster|<1.442 (1.566<|eta_supercluster|<2.5)
	// dEtaIn < 0.007 (0.009)
	// dPhiIn < 0.15 (0.10)
	// sigmaIEtaIEta < 0.01 (0.03)
	// d0 < 0.03 (0.03), d0 = -dxy
	// dz < 0.1 (0.1)
	// PFIsolation/pT < 0.2
	// Conversion rejection: vertex fit probability < 1e-6 (1e-6)
	// Conversion rejection: missing hits = 0 (0)
	if(fabs(Electron_supercluster_eta(i))<1.442){
		if(Electron_p4(i).Pt()>20 &&
				fabs(Electron_Poca(i).Z()-Vtx(j).Z()-((Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Px()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Py())*Electron_p4(i).Pz()/pow(Electron_p4(i).Pt(),2))<0.1 &&
				fabs((-(Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Py()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Px())/Electron_p4(i).Pt())<0.03 &&
				!Electron_HasMatchedConversions(i) &&
				Electron_sigmaIetaIeta(i)<0.01 &&
				Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.15 &&
				Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.007 &&
				Electron_RelIso(i)<0.2
				){
			return true;
		}
	}else if(fabs(Electron_supercluster_eta(i))>=1.442 && fabs(Electron_supercluster_eta(i))<2.5){
		if(Electron_p4(i).Pt()>20 &&
				fabs(Electron_Poca(i).Z()-Vtx(j).Z()-((Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Px()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Py())*Electron_p4(i).Pz()/pow(Electron_p4(i).Pt(),2))<0.1 &&
				fabs((-(Electron_Poca(i).X()-Vtx(j).X())*Electron_p4(i).Py()+(Electron_Poca(i).Y()-Vtx(j).Y())*Electron_p4(i).Px())/Electron_p4(i).Pt())<0.03 &&
				!Electron_HasMatchedConversions(i) &&
				Electron_sigmaIetaIeta(i)<0.03 &&
				Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.10 &&
				Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.009 &&
				Electron_RelIso(i)<0.2
				){
			return true;
		}
	}
	return false;
}


//Effective area estimated on Z->ee in 2012 data for dr=0.3
float Ntuple_Controller::Electron_Aeff(float Eta){
	double eta=fabs(Eta);
	if(eta<1.0)return 0.13;
	else if(eta>1.0 && eta<1.479)return 0.14;
	else if(eta>1.479 && eta<2.0)return 0.07;
	else if(eta>2.0 && eta<2.2)return 0.09;
	else if(eta>2.2 && eta<2.3)return 0.11;
	else if(eta>2.3 && eta<2.4)return 0.11;
	else if(eta>2.4)return 0.14;
}

float Ntuple_Controller::Electron_RelIso(unsigned int i){
  return (Electron_chargedHadronIso(i)+std::max((float)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff(Electron_p4(i).Eta())))/Electron_p4(i).Pt();
}



bool Ntuple_Controller::isGoodJet(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel 
  //  isGoodJet_nooverlapremoval(i) with:
  //  deltaR jet-electron cleaning < 0.4 (2 selected lepton only) < 0.3 (e+jets only) 0.3
  //  deltaR jet-muon cleaning < 0.4(2 selected lepton only) < 0.3 (mu+jets only), 0.1 for PF and JET 0.3
  if(isGoodJet_nooverlapremoval(i)){
    unsigned int muon_idx;
    return !jethasMuonOverlap(i,muon_idx);
  }
  return false;
}

bool Ntuple_Controller::jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx){
  for(unsigned int j=0;j<NMuons();j++){
    if(isTightMuon(j) && Muon_RelIso(j)<0.2){
      if(Tools::dr(Muons_p4(j),PFJet_p4(jet_idx))<0.4){ muon_idx=j;return true;}
    }
  }
  return false;
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
  //  if |?|<2.4, CHF>0 (patJet->chargedHadronEnergyFraction())
  //  if |?|<2.4, NCH>0 (patJet->chargedMultiplicity()) 
  /////////////////////////////////////////////////////////////////////////
  // apply jet ID
  bool JetID_ok=false;
  if(PFJet_numberOfDaughters(i)>1){
    if(PFJet_chargedEmEnergyFraction(i)<0.99){
      if(PFJet_neutralHadronEnergyFraction(i)<0.99){
	if(PFJet_neutralEmEnergyFraction(i)<0.99){
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

/*int Ntuple_Controller::PFJet_Corr_Track_idx(unsigned int i){
	for(unsigned j=0;j<PFJet_nTrk(i);j++){
		for(unsigned k=0;k<NTracks();k++){
			if(PFJet_TracksP4(i,j).DeltaR(Track_p4(k))<0.00001){
				return k;
			}
		}
	}
	return -1;
}*/

std::vector<int> Ntuple_Controller::PFJet_Corr_Track_idx(unsigned int i){
	std::vector<int> track_idx;
	track_idx.clear();
	for(unsigned j=0;j<PFJet_nTrk(i);j++){
		for(unsigned k=0;k<NTracks();k++){
			if(PFJet_TracksP4(i,j).DeltaR(Track_p4(k))<0.00001){
				track_idx.push_back(k);
			}
		}
	}
	return track_idx;
}

float     Ntuple_Controller::Track_parCov(unsigned int i, TrackPar par1, TrackPar par2){
  if(par1>par2)return Ntp->Track_parCov->at(i).at(par1).at(par2);
  return Ntp->Track_parCov->at(i).at(par2).at(par1);
}




double Ntuple_Controller::TauSpinerGet(TauSpinerInterface::TauSpinerType SpinType){
  if(!isData()){
    std::vector<SimpleParticle> tau_daughters, tau_daughters2;
    SimpleParticle tau, tau2;
    for(int i=0; i<NMCSignalParticles();i++){
      // check for signal Boson
      if(MCSignalParticle_pdgid(i)==25 || 
	 MCSignalParticle_pdgid(i)==36 || 
	 MCSignalParticle_pdgid(i)==22 || 
	 MCSignalParticle_pdgid(i)==23){
	if(MCSignalParticle_Tauidx(i).size()==2){
	  SimpleParticle X(MCSignalParticle_p4(i).Px(),
			   MCSignalParticle_p4(i).Py(),
			   MCSignalParticle_p4(i).Pz(),
			   MCSignalParticle_p4(i).E(),
			   MCSignalParticle_pdgid(i));
	  bool tau1good(false),tau2good(false);
	  //first tau
	  unsigned int tauidx=MCSignalParticle_Tauidx(i).at(0);
	  if(verbose)std::cout  << "tau 1 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      if(verbose)std::cout <<  tauidx << " " << t 
				   << "pdgid " << MCTauandProd_pdgid(tauidx,t)
				   << " px " << MCTauandProd_p4(tauidx,t).Px()
				   << " py " << MCTauandProd_p4(tauidx,t).Py()
				   << " pz " << MCTauandProd_p4(tauidx,t).Pz()
				   << " E " << MCTauandProd_p4(tauidx,t).E() << std::endl;
	      if(t==0){
		if(verbose) std::cout << "isTau" << std::endl;
		tau=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
				   MCTauandProd_p4(tauidx,t).Py(),
				   MCTauandProd_p4(tauidx,t).Pz(),
				   MCTauandProd_p4(tauidx,t).E(),
				   MCTauandProd_pdgid(tauidx,t));
	      }
	      else{
		if(verbose) std::cout << "isDaughter" << std::endl;
		if(tau_daughters.size()>0)tau1good=true;
		tau_daughters.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
						       MCTauandProd_p4(tauidx,t).Py(),
						       MCTauandProd_p4(tauidx,t).Pz(),
						       MCTauandProd_p4(tauidx,t).E(),
						       MCTauandProd_pdgid(tauidx,t)));
	      }
	    }
	    // second tau
	    tauidx=MCSignalParticle_Tauidx(i).at(1);
	    if(verbose)std::cout  << "tau 2 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	    if(tauidx<NMCTaus()){
	      for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
		if(verbose)std::cout <<  tauidx << " " << t
				     << "pdgid " << MCTauandProd_pdgid(tauidx,t)
				     << " px " << MCTauandProd_p4(tauidx,t).Px()
				     << " py " << MCTauandProd_p4(tauidx,t).Py()
				     << " pz " << MCTauandProd_p4(tauidx,t).Pz()
				     << " E " << MCTauandProd_p4(tauidx,t).E() << std::endl;
		
		if(t==0){
		  if(verbose) std::cout << "isTau" << std::endl;
		  tau2=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
				      MCTauandProd_p4(tauidx,t).Py(),
				      MCTauandProd_p4(tauidx,t).Pz(),
				      MCTauandProd_p4(tauidx,t).E(),
				      MCTauandProd_pdgid(tauidx,t));
		}
		else{
		  if(verbose) std::cout << "isDaughter" << std::endl;
		  if(tau_daughters.size()>0)tau2good=true;
		  tau_daughters2.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
							  MCTauandProd_p4(tauidx,t).Py(),
							  MCTauandProd_p4(tauidx,t).Pz(),
							  MCTauandProd_p4(tauidx,t).E(),
							  MCTauandProd_pdgid(tauidx,t)));
		}
	      }
	    }
	    if(tau1good && tau2good){
	      if(verbose)std::cout  << "Two Taus found: " << tau_daughters.size() << " " << tau_daughters2.size() << std::endl;
	      return TauSpinerInt.Get(SpinType,X,tau,tau_daughters,tau2,tau_daughters2);
	    }
	  }
	}
      }
    }
  }
  return 1.0;
}





bool Ntuple_Controller::hasSignalTauDecay(PdtPdgMini::PdgPDTMini parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &tau_idx){
  for(int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==parent_pdgid){
      for(int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	if(MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
	  std::cout << "Warning INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
	  return false;
	}
      }
      for(int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	unsigned int tauidx=MCSignalParticle_Tauidx(i).at(j);
	if(verbose)std::cout << "MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << " " << Ntp->MCTau_JAK->size() << " " << Ntp->MCTauandProd_pdgid->size() << std::endl;
	if(MCTau_JAK(tauidx)==tau_jak){ tau_idx=tauidx;Boson_idx=i;return true;}
      }
    }
  }
  return false;
}



bool Ntuple_Controller::isGoodKFTau(unsigned int i, unsigned int j){
  if(KFTau_discriminatorByKFit(i,j)){
    if(KFTau_discriminatorByQC(i,j)){
      if(PFTau_hpsDecayMode(KFTau_MatchedHPS_idx(i)) == 10){
	if(PFTau_isVLooseIsolationDBSumPtCorr(KFTau_MatchedHPS_idx(i))){
	  //if(PFTau_isMediumIsolationDBSumPtCorr(KFTau_MatchedHPS_idx(i))){
	  return true;
	}
      }
    }
  }
  return false;
}


 
bool Ntuple_Controller::TriggerAccept(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return TriggerAccept(i);
  return false;
}
unsigned int Ntuple_Controller::HLTPrescale(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return HLTPrescale(i);
  return 1;
}
unsigned int Ntuple_Controller::L1SEEDPrescale(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return L1SEEDPrescale(i);
  return 1;
}

bool Ntuple_Controller::GetTriggerIndex(TString n, unsigned int &i){
  for(i=0; i<Ntp->HTLTriggerName->size();i++){
    TString name=HTLTriggerName(i);
    if(name.Contains(n))return true;
  }
  return false;
}



float Ntuple_Controller::KFTau_Daughter_parCov(unsigned int i, unsigned int j,int par1,int par2){
  TMatrixF M(NKFTau_par,NKFTau_par);
  unsigned int idx(0);
  for(unsigned int k=0;k<NKFTau_par;k++){
    for(unsigned int l=0;l<=k;l++){
      if(idx<Ntp->KFTau_Daughter_parCov->at(i).at(j).size()){
        M[l][k]=Ntp->KFTau_Daughter_parCov->at(i).at(j).at(idx);
	M[k][l]=Ntp->KFTau_Daughter_parCov->at(i).at(j).at(idx);
      }
      idx++;
    }
  }
  return  M[par1][par2];
}

float Ntuple_Controller::KFTau_Daughter_inputparCov(unsigned int i, unsigned int j,int par1,int par2){
  TMatrixF M(NKFTau_par,NKFTau_par);
  unsigned int idx(0);
  for(unsigned int k=0;k<NKFTau_par;k++){
    for(unsigned int l=0;l<=k;l++){
      if(idx<Ntp->KFTau_Daughter_inputparCov->at(i).at(j).size()){
        M[l][k]=Ntp->KFTau_Daughter_inputparCov->at(i).at(j).at(idx);
        M[k][l]=Ntp->KFTau_Daughter_inputparCov->at(i).at(j).at(idx);
      }
      idx++;
    }
  }
  return  M[par1][par2];
}

TVector3  Ntuple_Controller::KFTau_RotatedVtx(unsigned int i){
  for(unsigned int j=0; j<KFTau_NDaughter(i);j++){
    if(abs(KFTau_Daughter_pdgid(i,j))==abs(PdtPdgMini::tau_plus)){
      return TVector3(Ntp->KFTau_Daughter_inputpar->at(i).at(j).at(KFTau_vx),
		     Ntp->KFTau_Daughter_inputpar->at(i).at(j).at(KFTau_vy),
		     Ntp->KFTau_Daughter_inputpar->at(i).at(j).at(KFTau_vz));
    }
  }
  return TVector3(0,0,0);
}

TMatrixF  Ntuple_Controller::KFTau_RotatedVtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0; j<KFTau_NDaughter(i);j++){
    if(abs(KFTau_Daughter_pdgid(i,j))==abs(PdtPdgMini::tau_plus)){
      for(unsigned int k=0;k<dim;k++){
	for(unsigned int l=0;l<=k;l++){
	  if(k*NKFTau_par+l<Ntp->KFTau_Daughter_inputparCov->at(i).at(j).size()){
	    M[k][l]=KFTau_Daughter_inputparCov(i,j,k,l);
	    M[l][k]=KFTau_Daughter_inputparCov(i,j,k,l);
	  }
	}
      }
      return M;
    }
  }
  return M;
}

TMatrixF  Ntuple_Controller::KFTau_ReducedVtx_Cov(){
  unsigned int i=0;
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      M[j][k]=Ntp->ReducedVtx_Cov->at(i).at(j).at(k);
      M[k][j]=Ntp->ReducedVtx_Cov->at(i).at(j).at(k);
    }
  }
  return M;
}


TVector3  Ntuple_Controller::KFTau_SecondayVtx(unsigned int i){
  for(unsigned int j=0; j<KFTau_NDaughter(i);j++){
    if(abs(KFTau_Daughter_pdgid(i,j))==abs(PdtPdgMini::tau_plus)){
      return TVector3(Ntp->KFTau_Daughter_par->at(i).at(j).at(KFTau_vx),
		      Ntp->KFTau_Daughter_par->at(i).at(j).at(KFTau_vy),
		      Ntp->KFTau_Daughter_par->at(i).at(j).at(KFTau_vz));
    }
  }
  return TVector3(0,0,0);
}


TMatrixF  Ntuple_Controller::KFTau_SecondaryVtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0; j<KFTau_NDaughter(i);j++){
    if(abs(KFTau_Daughter_pdgid(i,j))==abs(PdtPdgMini::tau_plus)){
      for(unsigned int k=0;k<dim;k++){
        for(unsigned int l=0;l<=k;l++){
          if(k*NKFTau_par+l<Ntp->KFTau_Daughter_parCov->at(i).at(j).size()){
            M[k][l]=KFTau_Daughter_parCov(i,j,k,l);
            M[l][k]=KFTau_Daughter_parCov(i,j,k,l);
          }
        }
      }
      return M;
    }
  }
  return M;
}



TVector3  Ntuple_Controller::KFTau_InitialSecondaryVtx(unsigned int i){
  for(unsigned int j=0; j<KFTau_NDaughter(i);j++){
    if(abs(KFTau_Daughter_pdgid(i,j))==abs(PdtPdgMini::pi_plus)){
      return TVector3(Ntp->KFTau_Daughter_inputpar->at(i).at(j).at(KFTau_vx),
		      Ntp->KFTau_Daughter_inputpar->at(i).at(j).at(KFTau_vy),
		      Ntp->KFTau_Daughter_inputpar->at(i).at(j).at(KFTau_vz));
    }
  }
  return TVector3(0,0,0);
}

TMatrixF  Ntuple_Controller::KFTau_InitialSecondaryVtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0; j<KFTau_NDaughter(i);j++){
    if(abs(KFTau_Daughter_pdgid(i,j))==abs(PdtPdgMini::pi_plus)){
      for(unsigned int k=0;k<dim;k++){
        for(unsigned int l=0;l<=k;l++){
          if(k*NKFTau_par+l<Ntp->KFTau_Daughter_inputparCov->at(i).at(j).size()){
            M[k][l]=KFTau_Daughter_inputparCov(i,j,k,l);
            M[l][k]=KFTau_Daughter_inputparCov(i,j,k,l);
          }
        }
      }
      return M;
    }
  }
  return M;
}

