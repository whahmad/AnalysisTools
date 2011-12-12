#include "Validation.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

Validation::Validation(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

Validation::~Validation(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Validation::~Validation Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Validation::~Validation()" << std::endl;
}

void  Validation::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    cut.at(TriggerOk)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=1;
    if(i==TauPt)        cut.at(TauPt) = 15;
    if(i==QC)           cut.at(QC) = 1;
    if(i==LooseIso)     cut.at(LooseIso) = 1;
    if(i==ET)           cut.at(ET) = 20;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==TauPt){
      title.at(i)="TauPt";
      hlabel="TauPt ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i==QC){
      title.at(i)="QC";
      hlabel="QC ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_QC_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_QC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==LooseIso){
      title.at(i)="LooseIso";
      hlabel="LooseIso ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LooseIso_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LooseIso_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==ET){
      title.at(i)="ET";
      hlabel="ET ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ET_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ET_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }

  } 



  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Accumulative Cuts Passed","Events");
  MuonM=HConfig.GetTH1D(Name+"_muonPt","muonPt",56,-0.5,55.5,"Muon Pt","Events");
  MET_et=HConfig.GetTH1D(Name+"_MET","MET",56,-0.5,55.5,"MET","Events");
  MET_phi=HConfig.GetTH1D(Name+"_METPHI","METPHI",70,-3.5,3.5,"METPHI","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  KFTau_TauFitPt=HConfig.GetTH1D(Name+"_TauPt","TauPt",56,-0.5,55.5,"Tau Pt","Events");
  PFTau_DecMode=HConfig.GetTH1D(Name+"_DecMode","DecMode",11,-0.5,10.5,"DecMode","Events");
  Muon_emEt03_hist=HConfig.GetTH1D(Name+"_Muon_emEt03","Muon_emEt03",41,-0.5,20.5,"Muon_emEt03","Events");	
  Muon_hadEt03_hist=HConfig.GetTH1D(Name+"_Muon_hadEt03","Muon_hadEt03",41,-0.5,20.5,"Muon_hadEt03","Events");	
  Muon_sumPt03_hist=HConfig.GetTH1D(Name+"_Muon_sumPt03","Muon_sumPt03",41,-0.5,20.5,"Muon_sumPt03","Events");	
  Muon_emEt05_hist=HConfig.GetTH1D(Name+"_Muon_emEt05","Muon_emEt05",41,-0.5,20.5,"Muon_emEt05","Events");	
  Muon_hadEt05_hist=HConfig.GetTH1D(Name+"_Muon_hadEt05","Muon_hadEt05",41,-0.5,20.5,"Muon_hadEt05","Events");	
  Muon_sumPt05_hist=HConfig.GetTH1D(Name+"_Muon_sumPt05","Muon_sumPt05",41,-0.5,20.5,"Muon_sumPt05","Events");	

  Muon_emVetoEt05_hist=HConfig.GetTH1D(Name+"_Muon_emVetoEt05","Muon_emVetoEt05",41,-0.5,20.5,"Muon_emVetoEt05","Events");	
  Muon_hadVetoEt05_hist=HConfig.GetTH1D(Name+"_Muon_hadVetoEt05","Muon_hadVetoEt05",41,-0.5,20.5,"Muon_hadVetoEt05","Events");	
  Muon_sumVetoPt05_hist=HConfig.GetTH1D(Name+"_Muon_sumVetoPt05","Muon_sumVetoPt05",41,-0.5,20.5,"Muon_sumVetoPt05","Events");	

  SecondMuonPt=HConfig.GetTH1D(Name+"_SecondMuonPt","SecondMuonPt",51,0.5,50.5,"SecondMuonPt","Events");	
  MuonIsol03=HConfig.GetTH1D(Name+"_MuonIsol03","MuonIsol03",40,-0.05,1.55,"MuonIsol03","Events");
  MuonIsol05=HConfig.GetTH1D(Name+"_MuonIsol05","MuonIsol05",40,-0.05,1.55,"MuonIsol05","Events");
 
  NFitTaus =HConfig.GetTH1D(Name+"_NFitTaus","NFitTaus",11,-0.05,10.05,"NFitTaus","Events");
  NQCTaus  =HConfig.GetTH1D(Name+"_NQCTaus","NQCTaus",11,-0.05,10.05,"NQCTaus" ,"Events");  
  NGlobalMuons  =HConfig.GetTH1D(Name+"_NGlobalMuons","NGlobalMuons",11,-0.05,10.05,"NGlobalMuons","Events");

  Muon_isGlobalMuon_hist =HConfig.GetTH1D(Name+"_Muon_isGlobalMuon","Muon_isGlobalMuon",2,0,2,"Muon_isGlobalMuon","Events");
  Muon_isStandAloneMuon_hist =HConfig.GetTH1D(Name+"_Muon_isStandAloneMuon","Muon_isStandAloneMuon",2,0,2,"Muon_isStandAloneMuon","Events");
  discrQC =HConfig.GetTH1D(Name+"_discrQC","discrQC",2,0,2,"discrQC","Events");
  discrFT =HConfig.GetTH1D(Name+"_discrFT","discrFT",2,0,2,"discrFT","Events");

  VisibleMass=HConfig.GetTH1D(Name+"_VisibleMass","VisibleMass",51,19.5,120.5,"VisibleMass","Events");
  FullMass=HConfig.GetTH1D(Name+"_FullMass","FullMass",51,19.5,120.5,"FullMass","Events");
  DeltaPhi=HConfig.GetTH1D(Name+"_DeltaPhi","DeltaPhi",40,0,6.28,"DeltaPhi","Events");
  DeltaRTauMu=HConfig.GetTH1D(Name+"_DeltaRTauMu","DeltaRTauMu",40,0,6,"DeltaRTauMu","Events");
  DeltaPhiTauMuNeutrino=HConfig.GetTH1D(Name+"_DeltaPhiTauMuNeutrino","DeltaPhiTauMuNeutrino",40,0,2,"DeltaPhiTauMuNeutrino","Events");
  Muon_nJets05_hist=HConfig.GetTH1D(Name+"_Muon_nJets05_hist","Muon_nJets05_hist",11,-0.05,10.05,"Muon_nJets05_hist","Events");
  Muon_nTracks05_hist=HConfig.GetTH1D(Name+"_Muon_nTracks05_hist","Muon_nTracks05_hist",11,-0.05,10.05,"Muon_nTracks05_hist","Events");
  DeltaVtxZ=HConfig.GetTH1D(Name+"_DeltaVtxZ","DeltaVtxZ",100,-2,2,"DeltaVtxZ","Events");
  GlobMuonEta=HConfig.GetTH1D(Name+"_GlobMuonEta","GlobMuonEta",51,-4,4,"GlobMuonEta","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",51,-4,4,"TauEta","Events");



  PFTau_isTightIsolation_hist=HConfig.GetTH1D(Name+"_PFTau_isTightIsolation_hist","PFTau_isTightIsolation_hist",2,0,2,"PFTau_isTightIsolation_hist","Events");
  PFTau_isMediumIsolation_hist=HConfig.GetTH1D(Name+"PFTau_isMediumIsolation_hist","PFTau_isMediumIsolation_hist",2,0,2 ,"PFTau_isMediumIsolation_hist","Events");
  PFTau_isLooseIsolation_hist=HConfig.GetTH1D(Name+"_PFTau_isLooseIsolation_hist","PFTau_isLooseIsolation_hist",2,0,2 ,"PFTau_isLooseIsolation_hist","Events");
  PFTau_hpsDecayMode_hist=HConfig.GetTH1D(Name+"_PFTau_hpsDecayMode_hist","PFTau_hpsDecayMode_hist",11,-0.05,10.05,"PFTau_hpsDecayMode_hist","Events");	     
  PFTau_Charge_hist=HConfig.GetTH1D(Name+"_PFTau_Charge_hist","PFTau_Charge_hist",4,-1,3,"PFTau_Charge_hist","Events");
  
  PFTau_Tau_Highestpt_hist=HConfig.GetTH1D(Name+"_PFTau_Tau_Highestpt_hist","PFTau_Tau_Highestpt_hist",50,0,50,"PFTau_Tau_Highestpt_hist","Events");      
  PFTau_Tau_Highestpt_phi_hist=HConfig.GetTH1D(Name+"_PFTau_Tau_Highestpt_phi_hist","PFTau_Tau_Highestpt_phi_hist",50,-3.14,3.14,"PFTau_Tau_Highestpt_phi_hist","Events");      



  KFTau_Fit_iterations_hist=HConfig.GetTH1D(Name+"_KFTau_Fit_iterations_hist","KFTau_Fit_iterations_hist",10,0,10,"KFTau_Fit_iterations_hist","Events");
  KFTau_Fit_chi2_hist=HConfig.GetTH1D(Name+"_KFTau_Fit_chi2_hist","KFTau_Fit_chi2_hist",20,0,5,"KFTau_Fit_chi2_hist","Events");	       
  KFTau_TauFit_pt_hist=HConfig.GetTH1D(Name+"_KFTau_TauFit_pt_hist","KFTau_TauFit_pt_hist",50,0,50,"KFTau_TauFit_pt_hist","Events");      
  KFTau_TauFit_phi_hist=HConfig.GetTH1D(Name+"_KFTau_TauFit_phi_hist","KFTau_TauFit_phi_hist",50,-3.14,3.14,"KFTau_TauFit_phi_hist","Events");      
  KFTau_TauVis_phi_hist=HConfig.GetTH1D(Name+"_KFTau_TauVis_phi_hist","KFTau_TauVis_phi_hist",50,-3.14,3.14,"KFTau_TauVis_phi_hist","Events");      
  KFTau_Neutrino_p4_hist=HConfig.GetTH1D(Name+"_KFTau_Neutrino_p4_hist","KFTau_Neutrino_p4_hist",50,0,50 ,"KFTau_Neutrino_p4_hist" ,"Events");   
  KFTau_Neutrino_phi_hist=HConfig.GetTH1D(Name+"_KFTau_Neutrino_phi_hist","KFTau_Neutrino_phi_hist",50,-3.14,3.14 ,"KFTau_Neutrino_phi_hist" ,"Events");   
  KFTau_TauVis_a1_hist=HConfig.GetTH1D(Name+"_KFTau_TauVis_a1_hist","KFTau_TauVis_a1_hist",100,0.3,1.7,"KFTau_TauVis_a1_hist","Events");      
  KFTau_Fit_csum_hist=HConfig.GetTH1D(Name+"_KFTau_Fit_csum_hist","KFTau_Fit_csum_hist",20,0,0.05,"KFTau_Fit_csum_hist","Events");	       
					   

  KFTau_TauFit_Highestpt_hist=HConfig.GetTH1D(Name+"_KFTau_TauFit_Highestpt_hist","KFTau_TauFit_Highestpt_hist",50,0,50,"KFTau_TauFit_Highestpt_hist","Events");      
  KFTau_TauFit_Highestpt_phi_hist=HConfig.GetTH1D(Name+"_KFTau_TauFit_Highestpt_phi_hist","KFTau_TauFit_Highestpt_phi_hist",50,-3.14,3.14,"KFTau_TauFit_Highestpt_phi_hist","Events");      
  PFTau_Tau_phi_hist=HConfig.GetTH1D(Name+"_PFTau_Tau_phi_hist","PFTau_Tau_phi_hist",50,-3.14,3.14,"PFTau_Tau_phi_hist","Events");      
  PFTau_Tau_pt_hist =HConfig.GetTH1D(Name+"_PFTau_Tau_pt_hist","PFTau_Tau_pt_hist",50,0,50,"PFTau_Tau_pt_hist","Events");   
  PFTau_Tau_Highestpt_hist =HConfig.GetTH1D(Name+"_PFTau_Tau_Highestpt_hist","PFTau_Tau_Highestpt_hist",50,0,50,"PFTau_Tau_Highestpt_hist","Events");   

 
					       
  KinFitTau_HPSMode_hist=HConfig.GetTH1D(Name+"_KinFitTau_HPSMode_hist","KinFitTau_HPSMode_hist",11,-0.05,10.05,"KinFitTau_HPSMode_hist","Events");     
  KinFitTau_HPSLooseIso_hist=HConfig.GetTH1D(Name+"_KinFitTau_HPSLooseIso_hist","KinFitTau_HPSLooseIso_hist",2,0,2,"KinFitTau_HPSLooseIso_hist","Events"); 
  KinFitTau_HPSMediumIso_hist=HConfig.GetTH1D(Name+"_KinFitTau_HPSMediumIso_hist","KinFitTau_HPSMediumIso_hist",2,0,2,"KinFitTau_HPSMediumIso_hist","Events");











  TransvereMass=HConfig.GetTH1D(Name+"_TransvereMass","TransvereMass",121,-0.5,120.5,"TransvereMass","Events");	
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  Validation::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&MuonM);
 Extradist1d.push_back(&MET_et);
 Extradist1d.push_back(&MET_phi);
 Extradist1d.push_back(&KFTau_TauFitPt);
 Extradist1d.push_back(&PFTau_DecMode);
 Extradist1d.push_back(&Muon_emEt03_hist);	 
 Extradist1d.push_back(&Muon_hadEt03_hist);	 
 Extradist1d.push_back(&Muon_sumPt03_hist);	 
 Extradist1d.push_back(&Muon_emEt05_hist);	 
 Extradist1d.push_back(&Muon_hadEt05_hist);	 
 Extradist1d.push_back(&Muon_sumPt05_hist);
 Extradist1d.push_back(&Muon_emVetoEt05_hist);	 
 Extradist1d.push_back(&Muon_hadVetoEt05_hist);	 
 Extradist1d.push_back(&Muon_sumVetoPt05_hist);	 

 Extradist1d.push_back(&MuonIsol03);		 
 Extradist1d.push_back(&MuonIsol05);          
 Extradist1d.push_back(&TransvereMass);  
 Extradist1d.push_back(&NFitTaus);          
 Extradist1d.push_back(&NQCTaus);          
 Extradist1d.push_back(&NGlobalMuons);          
 Extradist1d.push_back(&VisibleMass);
 Extradist1d.push_back(&FullMass);
 Extradist1d.push_back(&DeltaPhi);
 Extradist1d.push_back(&DeltaRTauMu);
 Extradist1d.push_back(&DeltaPhiTauMuNeutrino);
 Extradist1d.push_back(&DeltaVtxZ);


 Extradist1d.push_back(&Muon_isGlobalMuon_hist);
 Extradist1d.push_back(&Muon_isStandAloneMuon_hist);
 Extradist1d.push_back(&Muon_hadVetoEt05_hist);
 Extradist1d.push_back(&Muon_nTracks05_hist);
 Extradist1d.push_back(&PFTau_isTightIsolation_hist);
 Extradist1d.push_back(&PFTau_isMediumIsolation_hist);
 Extradist1d.push_back(&PFTau_isLooseIsolation_hist);
 Extradist1d.push_back(&PFTau_hpsDecayMode_hist);
 Extradist1d.push_back(&PFTau_Charge_hist);
 
 Extradist1d.push_back(&PFTau_Tau_Highestpt_hist);

 Extradist1d.push_back(&GlobMuonEta);
 Extradist1d.push_back(&TauEta);

 Extradist1d.push_back(&KFTau_Fit_iterations_hist);  
 Extradist1d.push_back(&KFTau_Fit_chi2_hist);	       
 Extradist1d.push_back(&KFTau_TauFit_pt_hist);       
 Extradist1d.push_back(&KFTau_TauFit_phi_hist); 
 Extradist1d.push_back(&KFTau_TauVis_phi_hist); 
 Extradist1d.push_back(&KFTau_Neutrino_p4_hist); 
 Extradist1d.push_back(&KFTau_Neutrino_phi_hist); 
 Extradist1d.push_back(&KFTau_TauVis_a1_hist);       
 Extradist1d.push_back(&KFTau_Fit_csum_hist);	       
					       
 Extradist1d.push_back(&KFTau_TauFit_Highestpt_hist);

 Extradist1d.push_back(&PFTau_Tau_phi_hist);					     
 Extradist1d.push_back(&PFTau_Tau_pt_hist);
 Extradist1d.push_back(&PFTau_Tau_Highestpt_phi_hist);

  
 Extradist1d.push_back(&KinFitTau_HPSMode_hist);     
 Extradist1d.push_back(&KinFitTau_HPSLooseIso_hist); 
 Extradist1d.push_back(&KinFitTau_HPSMediumIso_hist);
 Extradist1d.push_back(&KFTau_TauFit_Highestpt_phi_hist);


 Extradist1d.push_back(&SecondMuonPt);

 Extradist1d.push_back(&discrQC);
 Extradist1d.push_back(&discrFT);
}

void  Validation::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;


    for(unsigned int iTau=0;iTau < Ntp->NKFTau();iTau++){
      value.at(TauPt)= Ntp->KFTau_TauFit_p4(iTau).Pt();
      pass.at(TauPt)=true;//(value.at(TauPt)<=cut.at(TauPt));

      value.at(QC) = Ntp->KFTau_discriminatorByQC(iTau);
      pass.at(QC)=(value.at(QC)==cut.at(QC));


        if( Ntp->KFTau_discriminatorByKFit(iTau) ){
	  value.at(LooseIso) = Ntp->PFTau_isLooseIsolation(Ntp->KFTau_MatchedHPS_idx(iTau));
	  pass.at(LooseIso)=true;//(value.at(LooseIso)==cut.at(LooseIso));
	}
    }
    
      value.at(ET) = Ntp->MET_et();
      pass.at(ET)=(value.at(ET)<=cut.at(ET));

 


  
  double wobs=1;
  double w;

  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs);
 
  ///////////////////////////////////////////////////////////
  // Add plots
  //  if(status){
  if(true){
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
  



    unsigned int HighestPtMuonIndex=0;
    unsigned int SecondPtMuonIndex=0;
    unsigned int iMuon=0;
    float muonPt=0;
    int NGlMuons=0;
    int NQTaus =0;





    for(iMuon=0;iMuon<Ntp->NMuons(); iMuon++){

      Muon_isGlobalMuon_hist.at(t).Fill(Ntp->Muon_isGlobalMuon(iMuon));
      Muon_isStandAloneMuon_hist.at(t).Fill(Ntp->Muon_isStandAloneMuon(iMuon));

      if(Ntp->Muon_isGlobalMuon(iMuon)){  	NGlMuons++;       }

      if(Ntp->Muons_p4(iMuon).Pt() > muonPt){
	muonPt = Ntp->Muons_p4(iMuon).Pt();
	SecondPtMuonIndex = HighestPtMuonIndex;
	HighestPtMuonIndex = iMuon;
      }
    }


    unsigned int HighestPtTauIndex=0;
    float tauPt=0;
    unsigned int iTau=0;

    unsigned int HighestPtJetIndex=0;
    float jetPtBuffer=0;


    for(iTau=0;iTau < Ntp->NKFTau();iTau++){
      if(Ntp->KFTau_discriminatorByQC(iTau)){NQTaus++;}
      discrQC.at(t).Fill(Ntp->KFTau_discriminatorByQC(iTau));       
      discrFT.at(t).Fill(Ntp->KFTau_discriminatorByKFit(iTau));             

      if( Ntp->KFTau_indexOfFitInfo(iTau)!=-1){
	KFTau_Fit_iterations_hist.at(t).Fill(Ntp->KFTau_Fit_iterations(Ntp->KFTau_indexOfFitInfo(iTau)),w);
	KFTau_Fit_chi2_hist.at(t).Fill(Ntp->KFTau_Fit_chi2(Ntp->KFTau_indexOfFitInfo(iTau)),w);
	KFTau_Fit_csum_hist.at(t).Fill(Ntp->KFTau_Fit_csum(Ntp->KFTau_indexOfFitInfo(iTau)),w);
       }

      //  if(Ntp->KFTau_discriminatorByKFit(iTau) and Ntp->KFTau_discriminatorByQC(iTau) )
{
	KFTau_TauFit_pt_hist.at(t).Fill(Ntp->KFTau_TauFit_p4(iTau).Pt(),w);
	KFTau_TauFit_phi_hist.at(t).Fill(Ntp->KFTau_TauFit_p4(iTau).Phi(),w);
	KFTau_TauVis_phi_hist.at(t).Fill(Ntp->KFTau_TauVis_p4(iTau).Phi(),w);
	KFTau_Neutrino_p4_hist.at(t).Fill(Ntp->KFTau_Neutrino_p4(iTau).Pt(),w);
	KFTau_Neutrino_phi_hist.at(t).Fill(Ntp->KFTau_Neutrino_p4(iTau).Phi(),w);
	KFTau_TauVis_a1_hist.at(t).Fill(Ntp->KFTau_TauVis_p4(iTau).M(),w);
      }
      
        if(Ntp->NKFTau()>0 && Ntp->NPFTaus()>0 &&  Ntp->NPFTaus()>=Ntp->KFTau_MatchedHPS_idx(iTau) && Ntp->KFTau_discriminatorByKFit(iTau) ){
	  KinFitTau_HPSMode_hist.at(t).Fill(Ntp->PFTau_hpsDecayMode(Ntp->KFTau_MatchedHPS_idx(iTau)));
	  KinFitTau_HPSLooseIso_hist.at(t).Fill(Ntp->PFTau_isLooseIsolation(Ntp->KFTau_MatchedHPS_idx(iTau)));
	  KinFitTau_HPSMediumIso_hist.at(t).Fill(Ntp->PFTau_isMediumIsolation(Ntp->KFTau_MatchedHPS_idx(iTau)));
	  
	}
      
    
      if(Ntp->KFTau_TauFit_p4(iTau).Pt() > tauPt){
	tauPt = Ntp->KFTau_TauFit_p4(iTau).Pt();
	HighestPtTauIndex = iTau;
      }

    }
    KFTau_TauFit_Highestpt_hist.at(t).Fill(Ntp->KFTau_TauFit_p4(HighestPtTauIndex).Pt(),w);
    KFTau_TauFit_Highestpt_phi_hist.at(t).Fill(Ntp->KFTau_TauFit_p4(HighestPtTauIndex).Phi(),w);

    unsigned int iPfTau=0;
    for(iPfTau=0; iPfTau < Ntp->NPFTaus(); ++iPfTau  ){

      PFTau_Tau_phi_hist.at(t).Fill(Ntp->PFTau_p4(iPfTau).Phi(),w);
      PFTau_Tau_pt_hist.at(t).Fill(Ntp->PFTau_p4(iPfTau).Pt(),w);
      PFTau_isTightIsolation_hist.at(t).Fill(Ntp->PFTau_isTightIsolation(iPfTau),w);
      PFTau_isMediumIsolation_hist.at(t).Fill(Ntp->PFTau_isMediumIsolation(iPfTau),w);
      PFTau_isLooseIsolation_hist.at(t).Fill(Ntp->PFTau_isLooseIsolation(iPfTau),w);
      PFTau_hpsDecayMode_hist.at(t).Fill(Ntp->PFTau_hpsDecayMode(iPfTau),w);
      PFTau_Charge_hist.at(t).Fill(Ntp->PFTau_Charge(iPfTau),w);
      

      if(Ntp->PFTau_p4(iPfTau).Pt() > jetPtBuffer){
	jetPtBuffer = Ntp->PFTau_p4(iPfTau).Pt();
	HighestPtJetIndex= iPfTau;
      }
    }
    
    
    PFTau_Tau_Highestpt_hist.at(t).Fill(Ntp->PFTau_p4(HighestPtJetIndex).Pt(),w);
    PFTau_Tau_Highestpt_phi_hist.at(t).Fill(Ntp->PFTau_p4(HighestPtJetIndex).Phi(),w);

    NFitTaus.at(t).Fill(Ntp->KFTau_nKinTaus(),w);
    NQCTaus.at(t).Fill(NQTaus,w);  
    NGlobalMuons.at(t).Fill(NGlMuons,w);
    
    if(Ntp->NMuons()!=0){
      if(Ntp->Muon_isGlobalMuon(HighestPtMuonIndex)){
	
      MuonM.at(t).Fill(Ntp->Muons_p4(HighestPtMuonIndex).Pt(),w);
      Muon_emEt03_hist.at(t).Fill(Ntp->Muon_emEt03(HighestPtMuonIndex),w);  
      Muon_hadEt03_hist.at(t).Fill(Ntp->Muon_hadEt03(HighestPtMuonIndex),w); 
      Muon_sumPt03_hist.at(t).Fill(Ntp->Muon_sumPt03(HighestPtMuonIndex),w); 
      Muon_emEt05_hist.at(t).Fill(Ntp->Muon_emEt05(HighestPtMuonIndex),w);  
      Muon_hadEt05_hist.at(t).Fill(Ntp->Muon_hadEt05(HighestPtMuonIndex),w); 
      Muon_sumPt05_hist.at(t).Fill(Ntp->Muon_sumPt05(HighestPtMuonIndex),w); 
      Muon_emVetoEt05_hist.at(t).Fill(Ntp->Muon_emVetoEt05(HighestPtMuonIndex),w);  
      Muon_hadVetoEt05_hist.at(t).Fill(Ntp->Muon_hadVetoEt05(HighestPtMuonIndex),w); 
      Muon_sumVetoPt05_hist.at(t).Fill(Ntp->Muon_trackerVetoPt05(HighestPtMuonIndex),w); 

      Muon_nJets05_hist.at(t).Fill(Ntp->Muon_nJets05(HighestPtMuonIndex),w); 
      Muon_nTracks05_hist.at(t).Fill(Ntp->Muon_nTracks05(HighestPtMuonIndex),w); 

      MuonIsol03.at(t).Fill((Ntp->Muon_emEt03(HighestPtMuonIndex) + Ntp->Muon_hadEt03(HighestPtMuonIndex) + Ntp->Muon_sumPt03(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt(),w);	       
      MuonIsol05.at(t).Fill((Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt(),w);	   
      TransvereMass.at(t).Fill(sqrt(2*Ntp->Muons_p4(HighestPtMuonIndex).Pt()*Ntp->MET_et()*(1  - cos(Ntp->Muons_p4(HighestPtMuonIndex).Phi() - Ntp->MET_phi())) ),w);
      }
      
    }

 
    
    if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){

      if(Ntp->Muon_isGlobalMuon(HighestPtMuonIndex) && Ntp->KFTau_discriminatorByQC(HighestPtTauIndex)){
	TLorentzVector TauVis = Ntp->KFTau_TauVis_p4(HighestPtTauIndex);
	TLorentzVector TauFit = Ntp->KFTau_TauFit_p4(HighestPtTauIndex);
	TLorentzVector Neutrino= Ntp->KFTau_Neutrino_p4(HighestPtTauIndex);
	TLorentzVector MuoGlo = Ntp->Muons_p4(HighestPtMuonIndex);
	TLorentzVector ZVis = TauVis + MuoGlo;
	TLorentzVector ZFit = TauFit + MuoGlo;
   
 	GlobMuonEta.at(t).Fill(MuoGlo.Eta(),w);
 	TauEta.at(t).Fill(TauFit.Eta(),w);
      if( Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex)!=-1){
	DeltaVtxZ.at(t).Fill(Ntp->KFTau_Fit_TauPrimVtx(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex)).Z() - Ntp->Muon_Poca(HighestPtMuonIndex).Z(),w);
      }
	
  

    
  
	VisibleMass.at(t).Fill(ZVis.M(),w);
	FullMass.at(t).Fill(ZFit.M(),w);
	DeltaPhi.at(t).Fill(fabs(MuoGlo.Phi() - TauFit.Phi()),w);
	DeltaRTauMu.at(t).Fill(sqrt(pow(MuoGlo.Phi() - TauFit.Phi(),2) + pow(MuoGlo.Eta() - TauFit.Eta(),2)),w);
	DeltaPhiTauMuNeutrino.at(t).Fill(fabs(1 - fabs(TauVis.Phi() - (Ntp->MET_phi() - Neutrino.Phi()))/3.1415),w);


      }
    }
 
    if( Ntp->NMuons() > 1){
      TLorentzVector SecondMuon = Ntp->Muons_p4(SecondPtMuonIndex);
      SecondMuonPt.at(t).Fill(Ntp->Muons_p4(SecondPtMuonIndex).Pt(),w);

    }

//      for(unsigned int k=0;k<Ntp->NPFTaus();k++){
//        PFTau_DecMode.at(t).Fill(Ntp->PFTau_hpsDecayMode(k),w);
//      }
    std::cout<<"debug 6  "<<std::endl;
    MET_et.at(t).Fill(Ntp->MET_et(),w);
    MET_phi.at(t).Fill(Ntp->MET_phi(),w);

    for(unsigned int ii=0;ii<Ntp->NKFTau();ii++){
      if(Ntp->KFTau_discriminatorByQC(ii))  KFTau_TauFitPt.at(t).Fill(Ntp->KFTau_TauVis_p4(ii).Pt(),w);
    }


  }
}




