#include "TauPiMu.h"        
#include "TLorentzVector.h"
#include <cstdlib>
#include <iostream>    
#include "TVector3.h"  
#include "TF1.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "HistoConfig.h"
#include "SkimConfig.h"
#include "PDG_Var.h"
#include "Tools.h"
#include <cmath>
#include <TMath.h> 
#include <TROOT.h>
 
TauPiMu::TauPiMu(TString Name_, TString id_):
  Selection(Name_,id_)
{
  //  verbose=true;
}

TauPiMu::~TauPiMu(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    std::cout << "TauPiMu::~TauPiMu Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TauPiMu::~TauPiMu()" << std::endl;
}

void  TauPiMu::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==ZTTtoPiMu)         cut.at(ZTTtoPiMu)=1;
  }
  
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back(""); 
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    
    if(i==ZTTtoPiMu) {
      title.at(i)="Is $Z\\rightarrow\\tau\\tau\\rightarrow\\mu\\pi$ MC (bool)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Is Z#rightarrow#tau#tau#rightarrow#mu#pi MC (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isZtautautoPiMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isZtautautoPiMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }      
    //-----------
  }
  
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  
  // Setup Extra Histograms
  // Reco Objects
  // Tau_pt=HConfig.GetTH1D(Name+"_Tau_pt","Tau_pt",50,0,100,"Tau p_{t}","Events");
  // Tau_phi=HConfig.GetTH1D(Name+"_Tau_phi","Tau_phi",30,-TMath::Pi(),TMath::Pi(),"Tau #phi","Events");
  // Tau_eta=HConfig.GetTH1D(Name+"_Tau_eta","Tau_eta",50,-2.5,2.5,"Tau #eta","Events");
  // Tau_theta=HConfig.GetTH1D(Name+"_Tau_theta","Tau_theta",50,-2.5,2.5,"Tau #theta","Events");

  Pion_pt=HConfig.GetTH1D(Name+"_Pion_pt","Pion_pt",50,0,100,"Pion p_{t}","Events");
  Pion_phi=HConfig.GetTH1D(Name+"_Pion_phi","Pion_phi",30,-TMath::Pi(),TMath::Pi(),"Pion #phi","Events");
  Pion_eta=HConfig.GetTH1D(Name+"_Pion_eta","Pion_eta",50,-2.5,2.5,"Pion #eta","Events");
  Pion_theta=HConfig.GetTH1D(Name+"_Pion_theta","Pion_theta",50,-2.5,2.5,"Pion #theta","Events");

  Mu_pt=HConfig.GetTH1D(Name+"_Mu_pt","Mu_pt",50,0,100,"Muon p_{t}","Events");
  Mu_phi=HConfig.GetTH1D(Name+"_Mu_phi","Mu_phi",30,-TMath::Pi(),TMath::Pi(),"Muon #phi","Events");
  Mu_eta=HConfig.GetTH1D(Name+"_Mu_eta","Mu_eta",50,-2.5,2.5,"Muon #eta","Events");
  Mu_theta=HConfig.GetTH1D(Name+"_Mu_theta","Mu_theta",50,-2.5,2.5,"Muon #theta","Events");
  
  MET_phi=HConfig.GetTH1D(Name+"_MET_phi","MET_phi",30,-TMath::Pi(),TMath::Pi(),"MET #phi","Events");
  MET_eta=HConfig.GetTH1D(Name+"_MET_eta","MET_eta",50,-2.5,2.5,"MET #eta","Events");
  
  // Generator Objects
  //Z particle
  M_Z_gen=HConfig.GetTH1D(Name+"_M_Z_gen","M_Z_gen",180,50.,160.,"M_Z_gen","Events");
  Pt_Z_gen=HConfig.GetTH1D(Name+"_Pt_Z_gen","Pt_Z_gen",30,0.,80.,"Pt_Z_gen","Events");
  Eta_Z_gen=HConfig.GetTH1D(Name+"_Eta_Z_gen","Eta_Z_gen",50,-2.5,2.5,"Eta_Z_gen","Events");
  Phi_Z_gen=HConfig.GetTH1D(Name+"_Phi_Z_gen","Phi_Z_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Z_gen","Events");
  Theta_Z_gen=HConfig.GetTH1D(Name+"_Theta_Z_gen","Theta_Z_gen",32,0.,TMath::Pi(),"Theta_Z_gen","Events");
  //Z-> TauPi
  M_TauPi_gen=HConfig.GetTH1D(Name+"_M_TauPi_gen","M_TauPi_gen",180,0.,4.,"M_TauPi_gen","Events");
  Pt_TauPi_gen=HConfig.GetTH1D(Name+"_Pt_TauPi_gen","Pt_TauPi_gen",30,0.,70.,"Pt_TauPi_gen","Events");
  Eta_TauPi_gen=HConfig.GetTH1D(Name+"_Eta_TauPi_gen","Eta_TauPi_gen",50,-2.5,2.5,"Eta_TauPi_gen","Events");
  Phi_TauPi_gen=HConfig.GetTH1D(Name+"_Phi_TauPi_gen","Phi_TauPi_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_gen","Events");
  Theta_TauPi_gen=HConfig.GetTH1D(Name+"_Theta_TauPi_gen","Theta_TauPi_gen",32,0.,TMath::Pi(),"Theta_TauPi_gen","Events");
  //Z-> TauMu
  M_TauMu_gen=HConfig.GetTH1D(Name+"_M_TauMu_gen","M_TauMu_gen",180,0.,4.,"M_TauMu_gen","Events");
  Pt_TauMu_gen=HConfig.GetTH1D(Name+"_Pt_TauMu_gen","Pt_TauMu_gen",30,0.,70.,"Pt_TauMu_gen","Events");
  Eta_TauMu_gen=HConfig.GetTH1D(Name+"_Eta_TauMu_gen","Eta_TauMu_gen",50,-2.5,2.5,"Eta_TauMu_gen","Events");
  Phi_TauMu_gen=HConfig.GetTH1D(Name+"_Phi_TauMu_gen","Phi_TauMu_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_gen","Events");
  Theta_TauMu_gen=HConfig.GetTH1D(Name+"_Theta_TauMu_gen","Theta_TauMu_gen",32,0.,TMath::Pi(),"Theta_TauMu_gen","Events");
  // TauPi->Pi
  M_Pi_gen=HConfig.GetTH1D(Name+"_M_Pi_gen","M_Pi_gen",180,0.,0.3,"M_Pi_gen","Events");
  Pt_Pi_gen=HConfig.GetTH1D(Name+"_Pt_Pi_gen","Pt_Pi_gen",30,0.,70.,"Pt_Pi_gen","Events");
  Eta_Pi_gen=HConfig.GetTH1D(Name+"_Eta_Pi_gen","Eta_Pi_gen",50,-2.5,2.5,"Eta_Pi_gen","Events");
  Phi_Pi_gen=HConfig.GetTH1D(Name+"_Phi_Pi_gen","Phi_Pi_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Pi_gen","Events");
  Theta_Pi_gen=HConfig.GetTH1D(Name+"_Theta_Pi_gen","Theta_Pi_gen",32,0.,TMath::Pi(),"Theta_Pi_gen","Events");
  //TauMu->Mu
  M_Mu_gen=HConfig.GetTH1D(Name+"_M_Mu_gen","M_Mu_gen",180,0.,0.21,"M_Mu_gen","Events");
  Pt_Mu_gen=HConfig.GetTH1D(Name+"_Pt_Mu_gen","Pt_Mu_gen",30,0.,70.,"Pt_Mu_gen","Events");
  Eta_Mu_gen=HConfig.GetTH1D(Name+"_Eta_Mu_gen","Eta_Mu_gen",50,-2.5,2.5,"Eta_Mu_gen","Events");
  Phi_Mu_gen=HConfig.GetTH1D(Name+"_Phi_Mu_gen","Phi_Mu_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Mu_gen","Events");
  Theta_Mu_gen=HConfig.GetTH1D(Name+"_Theta_Mu_gen","Theta_Mu_gen",32,0.,TMath::Pi(),"Theta_Mu_gen","Events");
  //TauPi-> Nutp
  M_Nutp_gen=HConfig.GetTH1D(Name+"_M_Nutp_gen","M_Nutp_gen",180,0.,0.000000005,"M_Nutp_gen","Events");
  Pt_Nutp_gen=HConfig.GetTH1D(Name+"_Pt_Nutp_gen","Pt_Nutp_gen",30,0.,70.,"Pt_Nutp_gen","Events");
  Eta_Nutp_gen=HConfig.GetTH1D(Name+"_Eta_Nutp_gen","Eta_Nutp_gen",50,-2.5,2.5,"Eta_Nutp_gen","Events");
  Phi_Nutp_gen=HConfig.GetTH1D(Name+"_Phi_Nutp_gen","Phi_Nutp_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Nutp_gen","Events");
  Theta_Nutp_gen=HConfig.GetTH1D(Name+"_Theta_Nutp_gen","Theta_Nutp_gen",32,0.,TMath::Pi(),"Theta_Nutp_gen","Events");
  //TauMu-> Nutm
  M_Nutm_gen=HConfig.GetTH1D(Name+"_M_Nutm_gen","M_Nutm_gen",180,0.,0.000000005,"M_Nutm_gen","Events");
  Pt_Nutm_gen=HConfig.GetTH1D(Name+"_Pt_Nutm_gen","Pt_Nutm_gen",30,0.,70.,"Pt_Nutm_gen","Events");
  Eta_Nutm_gen=HConfig.GetTH1D(Name+"_Eta_Nutm_gen","Eta_Nutm_gen",50,-2.5,2.5,"Eta_Nutm_gen","Events");
  Phi_Nutm_gen=HConfig.GetTH1D(Name+"_Phi_Nutm_gen","Phi_Nutm_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Nutm_gen","Events");
  Theta_Nutm_gen=HConfig.GetTH1D(Name+"_Theta_Nutm_gen","Theta_Nutm_gen",32,0.,TMath::Pi(),"Theta_Nutm_gen","Events");
  //TauMu-> Num
  M_Num_gen=HConfig.GetTH1D(Name+"_M_Num_gen","M_Num_gen",180,0.,0.000000005,"M_Num_gen","Events");
  Pt_Num_gen=HConfig.GetTH1D(Name+"_Pt_Num_gen","Pt_Num_gen",30,0.,70.,"Pt_Num_gen","Events");
  Eta_Num_gen=HConfig.GetTH1D(Name+"_Eta_Num_gen","Eta_Num_gen",50,-2.5,2.5,"Eta_Num_gen","Events");
  Phi_Num_gen=HConfig.GetTH1D(Name+"_Phi_Num_gen","Phi_Num_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Num_gen","Events");
  Theta_Num_gen=HConfig.GetTH1D(Name+"_Theta_Num_gen","Theta_Num_gen",32,0.,TMath::Pi(),"Theta_Num_gen","Events");
  // DiTau (Z-> TauPi+TauMu) 
  M_DiTau_gen=HConfig.GetTH1D(Name+"_M_DiTau_gen","M_DiTau_gen",180,50.,150.,"M_DiTau_gen","Events");
  Pt_DiTau_gen=HConfig.GetTH1D(Name+"_Pt_DiTau_gen","Pt_DiTau_gen",30,0.,80.,"Pt_DiTau_gen","Events");
  dEta_DiTau_gen=HConfig.GetTH1D(Name+"_dEta_DiTau_gen","dEta_DiTau_gen",50,-2.5,2.5,"dEta_DiTau_gen","Events");
  dPhi_DiTau_gen=HConfig.GetTH1D(Name+"_dPhi_DiTau_gen","dPhi_DiTau_gen",32,-TMath::Pi(),TMath::Pi(),"dPhi_DiTau_gen","Events");
  dTheta_DiTau_gen=HConfig.GetTH1D(Name+"_dTheta_DiTau_gen","dTheta_DiTau_gen",32,0.,TMath::Pi(),"dTheta_DiTau_gen","Events");
  dM_DiTau_gen=HConfig.GetTH1D(Name+"_dM_DiTau_gen","dM_DiTau_gen",180,0.,50.,"dM_DiTau_Z_gen","Events");
  dR_DiTau_gen=HConfig.GetTH1D(Name+"_dR_DiTau_gen","dR_DiTau_gen",100,0.,8.,"dR_DiTau_gen","Events");
  // Physics objects(Pi,Mu)  
  M_MuPi_gen=HConfig.GetTH1D(Name+"_M_MuPi_gen","M_MuPi_gen",180,0.,0.300,"M_MuPi_gen","Events");
  Pt_MuPi_gen=HConfig.GetTH1D(Name+"_Pt_MuPi_gen","Pt_MuPi_gen",30,0.,70.,"Pt_MuPi_gen","Events");
  dEta_MuPi_gen=HConfig.GetTH1D(Name+"_dEta_MuPi_gen","dEta_MuPi_gen",50,-2.5,2.5,"dEta_MuPi_gen","Events");
  dPhi_MuPi_gen=HConfig.GetTH1D(Name+"_dPhi_MuPi_gen","dPhi_MuPi_gen",32,-TMath::Pi(),TMath::Pi(),"dPhi_MuPi_gen","Events");
  dTheta_MuPi_gen=HConfig.GetTH1D(Name+"_dTheta_MuPi_gen","dTheta_MuPi_gen",32,0.,TMath::Pi(),"dTheta_MuPi_gen","Events");
  dR_MuPi_gen=HConfig.GetTH1D(Name+"_dR_MuPi_gen","dR_MuPi_gen",100,0.,8.,"dR_MuPi_gen","Events");
  // //Delta_R(TauPi_rec,TauPi_gen
  // dR_recTauPi_genTauPi=HConfig.GetTH1D(Name+"_dR_recTauPi_genTauPi","dR_recTauPi_genTauPi",100,0.,8.,"dR_recTauPi_genTauPi","Events");
  // dR_recTauMu_genTauMu=HConfig.GetTH1D(Name+"_dR_recTauMu_genTauMu","dR_recTauMu_genTauMu",100,0.,8.,"dR_recTauMu_genTauMu","Events");
  
  //Delta_R(Mu_rec,Mu_gen), Delta_R(Pi_rec,Pi_gen)
  dR_recMu_genMu=HConfig.GetTH1D(Name+"_dR_recMu_genMu","dR_recMu_genMu",100,0.,8.,"dR_recMu_genMu","Events");
  dR_recPi_genPi=HConfig.GetTH1D(Name+"_dR_recPi_genPi","dR_recPi_genPi",100,0.,8.,"dR_recPi_genPi","Events");
  
  
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}


void  TauPiMu::Store_ExtraDist(){
  
  //Reco objets
  // Extradist1d.push_back(&Tau_pt);               Extradist1d.push_back(&Tau_eta);          
  // Extradist1d.push_back(&Tau_phi);              Extradist1d.push_back(&Tau_theta);          
  
  Extradist1d.push_back(&Pion_pt);              Extradist1d.push_back(&Pion_eta);          
  Extradist1d.push_back(&Pion_phi);             Extradist1d.push_back(&Pion_theta);
  
  Extradist1d.push_back(&Mu_pt);                Extradist1d.push_back(&Mu_eta);          
  Extradist1d.push_back(&Mu_phi);               Extradist1d.push_back(&Mu_theta);          
  
  Extradist1d.push_back(&MET_phi);              Extradist1d.push_back(&MET_eta);          
  
  //Generator objets   
  Extradist1d.push_back(&M_Z_gen);              Extradist1d.push_back(&Pt_Z_gen);           Extradist1d.push_back(&Eta_Z_gen);
  Extradist1d.push_back(&Phi_Z_gen);            Extradist1d.push_back(&Theta_Z_gen);
  
  Extradist1d.push_back(&M_TauPi_gen);          Extradist1d.push_back(&Pt_TauPi_gen);       Extradist1d.push_back(&Eta_TauPi_gen);
  Extradist1d.push_back(&Phi_TauPi_gen);        Extradist1d.push_back(&Theta_TauPi_gen);
  
  Extradist1d.push_back(&M_TauMu_gen);          Extradist1d.push_back(&Pt_TauMu_gen);       Extradist1d.push_back(&Eta_TauMu_gen);
  Extradist1d.push_back(&Phi_TauMu_gen);        Extradist1d.push_back(&Theta_TauMu_gen);
  
  Extradist1d.push_back(&M_Pi_gen);             Extradist1d.push_back(&Pt_Pi_gen);          Extradist1d.push_back(&Eta_Pi_gen);
  Extradist1d.push_back(&Phi_Pi_gen);           Extradist1d.push_back(&Theta_Pi_gen);
  
  Extradist1d.push_back(&M_Mu_gen);             Extradist1d.push_back(&Pt_Mu_gen);          Extradist1d.push_back(&Eta_Mu_gen);
  Extradist1d.push_back(&Phi_Mu_gen);           Extradist1d.push_back(&Theta_Mu_gen);
  
  Extradist1d.push_back(&M_Nutp_gen);           Extradist1d.push_back(&Pt_Nutp_gen);        Extradist1d.push_back(&Eta_Nutp_gen);
  Extradist1d.push_back(&Phi_Nutp_gen);         Extradist1d.push_back(&Theta_Nutp_gen);

  Extradist1d.push_back(&M_Nutm_gen);           Extradist1d.push_back(&Pt_Nutm_gen);        Extradist1d.push_back(&Eta_Nutm_gen);
  Extradist1d.push_back(&Phi_Nutm_gen);         Extradist1d.push_back(&Theta_Nutm_gen);
  
  Extradist1d.push_back(&M_Num_gen);            Extradist1d.push_back(&Pt_Num_gen);         Extradist1d.push_back(&Eta_Num_gen);
  Extradist1d.push_back(&Phi_Num_gen);          Extradist1d.push_back(&Theta_Num_gen);
  
  Extradist1d.push_back(&M_DiTau_gen);          Extradist1d.push_back(&Pt_DiTau_gen);       Extradist1d.push_back(&dEta_DiTau_gen);
  Extradist1d.push_back(&dPhi_DiTau_gen);       Extradist1d.push_back(&dTheta_DiTau_gen);   Extradist1d.push_back(&dM_DiTau_gen);
  Extradist1d.push_back(&dR_DiTau_gen);
  
  Extradist1d.push_back(&M_MuPi_gen);           Extradist1d.push_back(&Pt_MuPi_gen);        Extradist1d.push_back(&dEta_MuPi_gen);
  Extradist1d.push_back(&dPhi_MuPi_gen);        Extradist1d.push_back(&dTheta_MuPi_gen);    Extradist1d.push_back(&dR_MuPi_gen);
  
  //Extradist1d.push_back(&dR_recTauPi_genTauPi); Extradist1d.push_back(&dR_recTauMu_genTauMu);
  
  Extradist1d.push_back(&dR_recMu_genMu);       Extradist1d.push_back(&dR_recPi_genPi);
  
} // End of Store_ExtraDist


void  TauPiMu::doEvent() {
  unsigned int t;                    
  int id(Ntp->GetMCID()); //MCid
  if(!HConfig.GetHisto(Ntp->isData(),id,t)) {t=2;}//{ std::cout << "failed to find id" <<std::endl; return;}
  
  value.at(ZTTtoPiMu)=0;
  pass.at(ZTTtoPiMu) = true;
  if(pass.at(ZTTtoPiMu))value.at(ZTTtoPiMu)=1;
  
  double wobs=1;
  double w=1;
  // if(!Ntp->isData()) { w = Ntp->PUWeight(); }
  // else{w=1;}
  //if(verbose)  std::cout << Ntp->GetMCID() << " " << Npassed.size() << " " << t  << Ntp->NMCTaus() << std::endl;
  
  
  bool status=AnalysisCuts(t,w,wobs); 
  // Add plots
  if(status) {
    
    TLorentzVector Pion; 
    TLorentzVector Muon; 
    bool hasPionrec = false; 
    bool hasMurec   = false; 
    double dR_Pions=-1;
    double dR_Muons=-1;     
    // unsigned int HighestPtPionIndex=0;
    // unsigned int HighestPtMuonIndex=0;
    // float PionPt=0;
    // float muonPt=0;
    
    
    TLorentzVector GenTauMu, GenTauPi, GenZ, GenPi, GenMu, GenNutp, GenNutm, GenNum;
    TLorentzVector DiTau, MuPi;
    
    double Delta_M_DiTau, dEta_genTaus, dPhi_genTaus, dTheta_genTaus, dR_genTaus;
    double dEta_genMuPi, dPhi_genMuPi, dTheta_genMuPi, dR_genMu_genPi;
    
    bool hasZ    = false;  bool hasTaumu = false;  bool hasTaupi = false; 
    bool hasMu   = false;  bool hasPi    = false; 
    bool hasNutp = false;  bool hasNutm  = false;  bool hasNum   = false;  
    
    if(id != DataMCType::Data && (id == DataMCType::DY_tautau )) { //DY_tautau
      if(Ntp->NMCTaus() == 2) {  //  NMCTaus() == 2       
	
	////////////////////////////
	////     Reco Study     ////    
	
	/// Pion and Muon not classified 
	for(unsigned int iPion=0; iPion < Ntp->NPFTaus(); iPion++) {
	  //std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
	  Pion = Ntp->PFTau_p4(iPion);
	  hasPionrec  = true;
	}
	for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++) {
	  //std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
	  Muon = Ntp->Muon_p4(iMuon);
	  hasMurec = true;
	}
	
	// // Pion and Muon classified with highest_Pt ones
	//    if(Ntp->NPFTaus()!=0) {  
	//      for(unsigned int iPion=0;iPion < Ntp->NPFTaus();iPion++) {
	
	//    	std::cout<< "NPions==" << Ntp->NPFTaus() <<std::endl;	
	
	//    	if(Ntp->PFTau_p4(iPion).Pt() > PionPt) {
	//    	  PionPt = Ntp->PFTau_p4(iPion).Pt();
	//    	  HighestPtPionIndex = iPion;
	// 	}
	// 	std::cout<< "Pion_Pt==" << HighestPtPionIndex <<std::endl;
	// 	Pion = Ntp->PFTau_p4(HighestPtPionIndex);
	// 	hasPionrec  = true;
	//      }
	//    }
	//    if(Ntp->NMuons()!=0) {
	//      for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++) {
	
	//    	std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
	
	//    	if(Ntp->Muon_p4(iMuon).Pt() > muonPt) {
	//    	  muonPt = Ntp->Muon_p4(iMuon).Pt();
	//    	  HighestPtMuonIndex = iMuon;
	// 	}
	// 	std::cout<< "Muon_Pt==" << HighestPtMuonIndex <<std::endl;
	// 	Muon = Ntp->Muon_p4(HighestPtMuonIndex);
	  // 	hasMurec  = true;
	  //      }
	  //    }
	
	
	if (hasPionrec && hasMurec) {
	  
	  Pion_pt.at(t).Fill(Pion.Pt(),w);  
	  Pion_phi.at(t).Fill(Pion.Phi(),w);
	  Pion_eta.at(t).Fill(Pion.Eta(),w);
	  Pion_theta.at(t).Fill(Pion.Theta(),w);
	  
	  
	  Mu_pt.at(t).Fill(Muon.Pt(),w);
	  Mu_phi.at(t).Fill(Muon.Phi(),w);
	  Mu_eta.at(t).Fill(Muon.Eta(),w);
	  Mu_theta.at(t).Fill(Muon.Theta(),w);
	}
	
	////////////////////////////
	//// Generator Study my ////
	// if(id != DataMCType::Data && (id == DataMCType::DY_tautau )) { //DY_tautau
	//if(Ntp->NMCTaus() == 2) {  //  NMCTaus() == 2  
	
	for(unsigned int i=0; i<Ntp->NMCParticles(); i++) { // loop over MC Particles 
	  if(Ntp->MCParticle_pdgid(i) == PDGInfo::Z0) {     // Z bososn 
	    GenZ = Ntp->MCParticle_p4(i);
	    hasZ = true;
	  }	
	}   
	for(int i=0; i<Ntp->NMCTaus(); i++) {  //loop over NMCTaus
	  if(Ntp->MCTau_JAK(i) == TauDecay::JAK_MUON) {         //Tau->Muon   //JAK_MUON=2
	    GenTauMu = Ntp->MCTau_p4(i);
	    hasTaumu=true;
	    
	    for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++) {
	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::mu_minus) {
		GenMu = Ntp->MCTauandProd_p4(i,j); 
		hasMu = true; 
	      } 
	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_mu) {
		GenNum = Ntp->MCTauandProd_p4(i,j);
		hasNum = true; 
	      } 
	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_tau) {
		GenNutm = Ntp->MCTauandProd_p4(i,j);
		hasNutm = true;
	      }
	    }
	  } //end Tau->Muon   //JAK_MUON=2
	  
	  if(Ntp->MCTau_JAK(i) == TauDecay::JAK_PION) {//Tau->Pion  //JAK_PION=3
	    GenTauPi = Ntp->MCTau_p4(i);
	    hasTaupi = true;
	    
	    for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++) {
	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::pi_plus) {
		GenPi = Ntp->MCTauandProd_p4(i,j); 
		hasPi = true;
	      }
	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_tau) {
		GenNutp = Ntp->MCTauandProd_p4(i,j);
		hasNutp = true;
	      }
	    }
	  }  //end Tau->Pion  //JAK_PION=3
	  
	} //end loop over NMCTaus
	
      }//end Taus && NMCTaus() == 2
      
      //if(hasZ && hasTaumu && hasTaupi && hasMu && hasPi && hasNum && hasNutm && hasNutp) { // fill variables
      if(hasZ && hasMu && hasPi ) { // fill variables
	
	M_Z_gen.at(t).Fill(GenZ.M(),w);              Pt_Z_gen.at(t).Fill(GenZ.Pt(),w);                  Eta_Z_gen.at(t).Fill(GenZ.Eta(),w);
	Phi_Z_gen.at(t).Fill(GenZ.Phi(),w);          Theta_Z_gen.at(t).Fill(GenZ.Theta(),w);
	
	M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);      Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);          Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);
	Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w);  Theta_TauMu_gen.at(t).Fill(GenTauMu.Theta(),w);
	M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);      Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);          Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w); 
	Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w);  Theta_TauPi_gen.at(t).Fill(GenTauPi.Theta(),w); 
	
	M_Mu_gen.at(t).Fill(GenMu.M(),w);            Pt_Mu_gen.at(t).Fill(GenMu.Pt(),w);                Eta_Mu_gen.at(t).Fill(GenMu.Eta(),w); 
	Phi_Mu_gen.at(t).Fill(GenMu.Phi(),w);        Theta_Mu_gen.at(t).Fill(GenMu.Theta(),w); 
	M_Pi_gen.at(t).Fill(GenPi.M(),w);   	     Pt_Pi_gen.at(t).Fill(GenPi.Pt(),w);                Eta_Pi_gen.at(t).Fill(GenPi.Eta(),w); 
	Phi_Pi_gen.at(t).Fill(GenPi.Phi(),w);        Theta_Pi_gen.at(t).Fill(GenPi.Theta(),w);    	
	
	M_Nutp_gen.at(t).Fill(GenNutp.M(),w);        Pt_Nutp_gen.at(t).Fill(GenNutp.Pt(),w);            Eta_Nutp_gen.at(t).Fill(GenNutp.Eta(),w); 
	Phi_Nutp_gen.at(t).Fill(GenNutp.Phi(),w);    Theta_Nutp_gen.at(t).Fill(GenNutp.Theta(),w); 
	M_Nutm_gen.at(t).Fill(GenNutm.M(),w);        Pt_Nutm_gen.at(t).Fill(GenNutm.Pt(),w);            Eta_Nutm_gen.at(t).Fill(GenNutm.Eta(),w); 
	Phi_Nutm_gen.at(t).Fill(GenNutm.Phi(),w);    Theta_Nutm_gen.at(t).Fill(GenNutm.Theta(),w); 
	M_Num_gen.at(t).Fill(GenNum.M(),w);          Pt_Num_gen.at(t).Fill(GenNum.Pt(),w);             	Eta_Num_gen.at(t).Fill(GenNum.Eta(),w); 
	Phi_Num_gen.at(t).Fill(GenNum.Phi(),w);      Theta_Num_gen.at(t).Fill(GenNum.Theta(),w); 
	
	//DiTau_gen
	DiTau = GenTauPi + GenTauMu;
	M_DiTau_gen.at(t).Fill(DiTau.M(),w);                   Pt_DiTau_gen.at(t).Fill(DiTau.Pt(),w);
	Delta_M_DiTau = DiTau.M() - GenZ.M();                  dM_DiTau_gen.at(t).Fill(Delta_M_DiTau,w);	
	dEta_genTaus = GenTauPi.Eta() - GenTauMu.Eta();        dEta_DiTau_gen.at(t).Fill(dEta_genTaus,w);
	dPhi_genTaus = GenTauPi.Phi() - GenTauMu.Phi();        dPhi_DiTau_gen.at(t).Fill(dPhi_genTaus,w);
	dTheta_genTaus = GenTauPi.Theta() - GenTauMu.Theta();  dTheta_DiTau_gen.at(t).Fill(dTheta_genTaus,w);
	dR_genTaus = sqrt(pow(dEta_genTaus,2) + pow(dPhi_genTaus,2));      
	dR_DiTau_gen.at(t).Fill(dR_genTaus,w);
	//MuPi_gen
	MuPi = GenMu + GenPi;
	M_MuPi_gen.at(t).Fill(GenMu.M()+GenPi.M(),w);          Pt_MuPi_gen.at(t).Fill(GenMu.Pt()+GenPi.Pt(),w);
	dEta_genMuPi = GenMu.Eta() - GenPi.Eta();              dEta_MuPi_gen.at(t).Fill(dEta_genMuPi,w);
	dPhi_genMuPi = GenMu.Phi() - GenPi.Phi();              dPhi_MuPi_gen.at(t).Fill(dPhi_genMuPi,w); 
	dTheta_genMuPi = GenMu.Theta() - GenPi.Theta();        dTheta_MuPi_gen.at(t).Fill(dTheta_genMuPi,w);
	
	dR_genMu_genPi = sqrt(pow(GenMu.Eta() - GenPi.Eta(),2) + pow( GenMu.Phi() - GenPi.Phi(),2));      
	//dR_genMu_genPi = sqrt(pow(dEta_genMuPi,2) + pow(dPhi_genMuPi,2));      
	dR_MuPi_gen.at(t).Fill(dR_genMu_genPi,w);
        
	// DeltaR(Tau_gen,Tau_rec)  , DeltaR(Mu_gen,Mu_rec), DeltaR(Pi_gen,Pi_rec)
	//dR_recTauPi_genTauPi.at(t).Fill(,w);            
	//dR_recTauMu_genTauMu.at(t).Fill(,w); 
	//dR_recMu_genMu.at(t).Fill(,w); 
	//dR_recPi_genPi.at(t).Fill(,w);
	
      } // end fill variables
      
      /////////////////////////////////////////
      // DR(Pirec , Pigen)   , DR(Murec ,Mugen)
      ////////////////////////////////////////
      
      //if (hasMu && hasPi && hasPionrec && hasMurec) {
      
      // dR_Pions = sqrt(pow(Pion.Eta() - GenPi.Eta(),2) + pow(Pion.Phi() - GenPi.Phi(),2));
      //dR_Muons = sqrt(pow(Muon.Eta() - GenMu.Eta(),2) + pow(Muon.Phi() - GenMu.Phi(),2));
	
      // std::cout<< "DR_Pion==" << dR_Pions <<std::endl;
      // std::cout<< "DR_Muon==" << dR_Muons <<std::endl;
      
      // dR_recPi_genPi.at(t).Fill(dR_Pions,w);   
      //dR_recMu_genMu.at(t).Fill(dR_Muons,w);
      //}
      
   
      //if (hasMu && hasPi && hasPionrec && hasMurec) {
      if (hasMu && hasPi) {
	//dR_Pions=Tools::dr(Ntp->PFTau_p4(i),GenPi);
	//dR_Muons=Tools::dr(Ntp->Muon_p4(i),GenMu);
	
	for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
	  std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
	  dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi) ;
	  std::cout<< "DR_Pion==" << dR_Pions <<std::endl;
	  dR_recPi_genPi.at(t).Fill(dR_Pions,w); 
	}
	for(unsigned int i=0; i < Ntp->NMuons(); i++) {
	  std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
	  dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu) ;
	  std::cout<< "DR_Muon==" << dR_Muons <<std::endl;
	  dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
	}
      }
      
    } //end DY_tautau
  } // End if(status)
}// End doEvent() 


void  TauPiMu::Finish(){
  
  std::cout << "Starting Finish" << std::endl;
  
  std::cout << "Starting Selection::Finish" << std::endl;
  
  Selection::Finish();
}












///////////////////////////
//Backup
///////////////////////////


    // Tau and Muon classified with highest_Pt ones
    // if(Ntp->NPFTaus()!=0) {  
    //   for(unsigned int iPion=0;iPion < Ntp->NPFTaus();iPion++) {
	
    // 	std::cout<< "NPFTaus==" << Ntp->NPFTaus() <<std::endl;	
	
    // 	if(Ntp->PFTau_p4(iPion).Pt() > tauPt) {
    // 	  tauPt = Ntp->PFTau_p4(iPion).Pt();
    // 	  HighestPtPionIndex = iPion;
	  
    // 	  std::cout<< "Pion_Pt==" << HighestPtPionIndex <<std::endl;
	  
    // 	  Pion = Ntp->PFTau_p4(HighestPtPionIndex);
    // 	  hasTaurec  = true;
    // 	}
    //   }
    // }
    // if(Ntp->NMuons()!=0) {
    //   for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++) {
	
    // 	std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
	
    // 	if(Ntp->Muon_p4(iMuon).Pt() > muonPt) {
    // 	  muonPt = Ntp->Muon_p4(iMuon).Pt();
    // 	  HighestPtMuonIndex = iMuon;
	  
    // 	  std::cout<< "Muon_Pt==" << HighestPtMuonIndex <<std::endl;
	  
    // 	  Muon = Ntp->Muon_p4(HighestPtMuonIndex);
    // 	  hasMurec  = true;
    // 	}
    //   }
    // }


    // for(unsigned int iPi=0; iPi<Ntp->PFTau_PionsP4.size(); iPi++) {
    // 	std::cout<< "NPFTau_Pions==" << Ntp->PFTau_PionsP4.size() <<std::endl;
    // 	Poin = Ntp->PFTau_PionsP4(iPi);
    // 	hasPionrec  = true;
    //   }

	  // Pion_pt.at(t).Fill(Pion.Pt(),w);  
	  // Pion_phi.at(t).Fill(Pion.Phi(),w);
	  // Pion_eta.at(t).Fill(Pion.Eta(),w);
	  // Pion_theta.at(t).Fill(Pion.Theta(),w); 


      // for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
      // 	std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
      
      // 	if(Tools::dr(Ntp->PFTau_p4(i),GenPi)>=0){
      // 	  dR_Pions=Tools::dr(Ntp->PFTau_p4(i),GenPi);
      // 	  std::cout<< "DR_Pion==" << dR_Pions <<std::endl;
      // 	  // dR_recPi_genPi.at(t).Fill(dR_Pions,w); 
      // 	}
      // }
      // for(unsigned int i=0; i < Ntp->NMuons(); i++) {
      // 	std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
      
      // 	if(Tools::dr(Ntp->Muon_p4(i),GenMu)>=0){
      // 	  dR_Muons=Tools::dr(Ntp->Muon_p4(i),GenMu);
      // 	  std::cout<< "DR_Muon==" << dR_Muons <<std::endl;
      // 	  //dR_recMu_genMu.at(t).Fill(dR_Muons,w);  
      // 	}	
      // }
