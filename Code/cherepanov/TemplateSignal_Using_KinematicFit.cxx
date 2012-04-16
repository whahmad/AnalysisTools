#include "TemplateSignal_Using_KinematicFit.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"





TemplateSignal_Using_KinematicFit::TemplateSignal_Using_KinematicFit(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

TemplateSignal_Using_KinematicFit::~TemplateSignal_Using_KinematicFit(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TemplateSignal_Using_KinematicFit::~TemplateSignal_Using_KinematicFit Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TemplateSignal_Using_KinematicFit::~TemplateSignal_Using_KinematicFit()" << std::endl;
}

void  TemplateSignal_Using_KinematicFit::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==MuonisGlob)         cut.at(MuonisGlob)=1;
    if(i==TauIsQuality)       cut.at(TauIsQuality)=1;
    if(i==MuonPt)             cut.at(MuonPt)=24; //18
    if(i==TauPtCut)           cut.at(TauPtCut)=24;
    if(i==MET)                cut.at(MET)=30;
    if(i==MuonIso)            cut.at(MuonIso)=0.2;
    if(i==TauIsIso)           cut.at(TauIsIso)=1;
    if(i==charge)             cut.at(charge)=-1;


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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==MuonisGlob){
      title.at(i)="Muon is global ";
      hlabel="MuonisGlob";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==TauIsQuality){
      title.at(i)="TauIsQuality ";
      hlabel="TauIsQuality ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsQuality_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsQuality_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==MuonPt){
      title.at(i)="$MuonPt > $";
      title.at(i)+=cut.at(MuonPt);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#mu_{pT}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonPt_",htitle,20,0,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonPt_",htitle,20,0,50,hlabel,"Events"));
    }
   else if(i==TauPtCut){
      title.at(i)="$TauPtCut > $";
      title.at(i)+=cut.at(TauPtCut);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#tau_{pT}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPtCut_",htitle,20,0,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPtCut_",htitle,20,0,50,hlabel,"Events"));
    }
   else if(i==MET){
      title.at(i)="$MET < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MET";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,10,0,60,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,10,0,60,hlabel,"Events"));
    }  
   else if(i==MuonIso){
      title.at(i)="$MuonIso < $";
      title.at(i)+=cut.at(MuonIso);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#MuonIso";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIso_",htitle,30,0,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIso_",htitle,30,0,1,hlabel,"Events"));
    }
   else if(i==charge){
      title.at(i)="$opposite charge$";
      title.at(i)+=cut.at(charge);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
    } 

   else if(i==TauIsIso){
      title.at(i)="$TauIsIso == $";
      title.at(i)+=cut.at(TauIsIso);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="TauIsIso";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIso_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIso_",htitle,2,-0.5,1.5,hlabel,"Events"));
    } 


  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Accumulative Cuts Passed","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  Tau3PiPt=HConfig.GetTH1D(Name+"_Tau3PiPt","Tau3PiPt",30,0,100,"Tau3PiPt","Events");
  Tau1Pt =HConfig.GetTH1D(Name+"_Tau1Pt","Tau1Pt",20,0,80,"Tau1Pt","Events");
  Tau2Pt =HConfig.GetTH1D(Name+"_Tau2Pt","Tau2Pt",20,0,80,"Tau2Pt","Events");
  Tau1theta=HConfig.GetTH1D(Name+"_Tau1theta","Tau1theta",20,0,3.14,"Tau1theta","Events");
  Tau2theta=HConfig.GetTH1D(Name+"_Tau2theta","Tau2theta",20,0,3.14,"Tau2theta","Events");  
  Tau1Deltatheta=HConfig.GetTH1D(Name+"_Tau1Deltatheta","Tau1Deltatheta",20,0,3.14,"Tau1Deltatheta","Events");
  Tau2Deltatheta=HConfig.GetTH1D(Name+"_Tau2Deltatheta","Tau2Deltatheta",20,0,3.14,"Tau2Deltatheta","Events");  
  DeltaTheta_1st=HConfig.GetTH1D(Name+"_DeltaTheta_1st","DeltaTheta_1st",20,0,3.14,"DeltaTheta_1st","Events");
  DeltaTheta_2nd=HConfig.GetTH1D(Name+"_DeltaTheta_2nd","DeltaTheta_2nd",20,0,3.14,"DeltaTheta_2nd","Events");
  Energy_resolution1=HConfig.GetTH1D(Name+"_Energy_resolution1","Energy_resolution1",30,-100,100,"Energy_resolution1","Events");
  Energy_resolution2=HConfig.GetTH1D(Name+"_Energy_resolution2","Energy_resolution2",30,-100,100,"Energy_resolution2","Events");
  TransverseEnergy_resolution=HConfig.GetTH1D(Name+"_TransverseEnergy_resolution","TransverseEnergy_resolution",30,100,100,"TransverseEnergy_resolution","Events");


  xmuon_plus=HConfig.GetTH1D(Name+"_xmuon_plus","xmuon_plus",15,0.2,1.4,"xmuon_plus","Events");
  xmuon_minus=HConfig.GetTH1D(Name+"_xmuon_minus","xmuon_minus",15,0.2,1.4,"xmuon_minus","Events");

  xmuonT_plus=HConfig.GetTH1D(Name+"_xmuonT_plus","xmuonT_plus",15,0,2,"xmuonT_plus","Events");
  xmuonT_minus=HConfig.GetTH1D(Name+"_xmuonT_minus","xmuonT_minus",15,0,2,"xmuonT_minus","Events");

  xmuonTruthT_plus=HConfig.GetTH1D(Name+"_xmuonTruthT_plus","xmuonTruthT_plus",15,0,2,"xmuonTruthT_plus","Events");
  xmuonTruthT_minus=HConfig.GetTH1D(Name+"_xmuonTruthT_minus","xmuonTruthT_minus",15,0,2,"xmuonTruthT_minus","Events");

  xmuonTruth_plus=HConfig.GetTH1D(Name+"_xmuonTruth_plus","xmuonTruth_plus",15,0.2,1.4,"xmuonTruth_plus","Events");
  xmuonTruth_minus=HConfig.GetTH1D(Name+"_xmuonTruth_minus","xmuonTruth_minus",15,0.2,1.4,"xmuonTruth_minus","Events");

  xmuonMC_plus=HConfig.GetTH1D(Name+"_xmuonMC_plus","xmuonMC_plus",15,0.2,1.4,"xmuonMC_plus","Events");
  xmuonMC_minus=HConfig.GetTH1D(Name+"_xmuonMC_minus","xmuonMC_minus",15,0.2,1.4,"xmuonMC_minus","Events");

  xmuonTMC_plus=HConfig.GetTH1D(Name+"_xmuonTMC_plus","xmuonTMC_plus",15,0,2,"xmuonTMC_plus","Events");
  xmuonTMC_minus=HConfig.GetTH1D(Name+"_xmuonTMC_minus","xmuonTMC_minus",15,0,2,"xmuonTMC_minus","Events");



  TruthResolution=HConfig.GetTH1D(Name+"_TruthResolution","TruthResolution",30,-100,100,"TruthResolution","Events");
  TruthDeltaTheta=HConfig.GetTH1D(Name+"_TruthDeltaTheta","TruthDeltaTheta",20,0,3.14,"TruthDeltaTheta","Events");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}


// TString name,TString title,int nbinsx, double minx, double maxx, 
// 				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis)

void  TemplateSignal_Using_KinematicFit::Store_ExtraDist(){

 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&Tau3PiPt);
 Extradist1d.push_back(&Tau1Pt);
 Extradist1d.push_back(&Tau2Pt);
 Extradist1d.push_back(&Tau1theta);	    
 Extradist1d.push_back(&Tau2theta);      
 Extradist1d.push_back(&Tau1Deltatheta); 
 Extradist1d.push_back(&Tau2Deltatheta); 
 Extradist1d.push_back(&DeltaTheta_1st);
 Extradist1d.push_back(&DeltaTheta_2nd);
 Extradist1d.push_back(&Energy_resolution1);
 Extradist1d.push_back(&Energy_resolution2);
 Extradist1d.push_back(&TransverseEnergy_resolution);

 Extradist1d.push_back(&xmuon_plus);	      
 Extradist1d.push_back(&xmuonT_plus);      
 Extradist1d.push_back(&xmuonTruthT_plus); 
 Extradist1d.push_back(&xmuonTruth_plus);  
 Extradist1d.push_back(&xmuonMC_plus);     
 Extradist1d.push_back(&xmuonTMC_plus);    
				      
 Extradist1d.push_back(&xmuon_minus);      
 Extradist1d.push_back(&xmuonT_minus);     
 Extradist1d.push_back(&xmuonTruthT_minus);
 Extradist1d.push_back(&xmuonTruth_minus); 
 Extradist1d.push_back(&xmuonMC_minus);    
 Extradist1d.push_back(&xmuonTMC_minus);   

 Extradist1d.push_back(&TruthResolution);
 Extradist1d.push_back(&TruthDeltaTheta);

}

void  TemplateSignal_Using_KinematicFit::doEvent(){
 
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
 

  unsigned int HighestPtMuonIndex=0;
  unsigned int HighestPtTauIndex=0;
  unsigned int SecondPtMuonIndex =0;
  float muonPt=0;
  float tauPt=0;
  unsigned int iTau=0;
  unsigned int iMuon=0;
  for(iTau=0;iTau < Ntp->NKFTau();iTau++){
    if(Ntp->KFTau_TauFit_p4(iTau).Pt() > tauPt){
      tauPt = Ntp->KFTau_TauFit_p4(iTau).Pt();
      HighestPtTauIndex = iTau; 
    }
  }
  
  if(Ntp->NMuons()!=0){
    for(iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
      if(Ntp->Muons_p4(iMuon).Pt() > muonPt){
	muonPt = Ntp->Muons_p4(iMuon).Pt();
	SecondPtMuonIndex = HighestPtMuonIndex;
	HighestPtMuonIndex = iMuon;
      }
    }
  }
  
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;

  if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){

  value.at(MuonisGlob) = Ntp->Muon_isGlobalMuon(HighestPtMuonIndex);
  pass.at(MuonisGlob)  = (value.at(MuonisGlob)==cut.at(MuonisGlob));

  value.at(TauIsQuality) =  Ntp->KFTau_discriminatorByKFit(HighestPtTauIndex);//Ntp->KFTau_discriminatorByQC(HighestPtTauIndex);
  pass.at(TauIsQuality)  =  (value.at(TauIsQuality)==cut.at(TauIsQuality));

  value.at(TauPtCut) = Ntp->KFTau_TauFit_p4(HighestPtTauIndex).Pt();
  pass.at(TauPtCut)=(value.at(TauPtCut)>=cut.at(TauPtCut));

  value.at(MuonPt) = Ntp->Muons_p4(HighestPtMuonIndex).Pt();
  pass.at(MuonPt)=(value.at(MuonPt) >=cut.at(MuonPt));

  value.at(MET) =   sqrt(2*Ntp->Muons_p4(HighestPtMuonIndex).Pt()*Ntp->MET_et()*(1  - cos(Ntp->Muons_p4(HighestPtMuonIndex).Phi() - Ntp->MET_phi())) );
  pass.at(MET)=(value.at(MET)<=cut.at(MET));

  value.at(TauIsIso) =   Ntp->PFTau_isMediumIsolation(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex));
  pass.at(TauIsIso)=(value.at(TauIsIso) ==cut.at(TauIsIso));


  value.at(MuonIso) = (Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt();
  pass.at(MuonIso)=(value.at(MuonIso)<=cut.at(MuonIso));


  value.at(charge) =Ntp->PFTau_Charge(HighestPtTauIndex)*Ntp->Muon_Charge(HighestPtMuonIndex);
  pass.at(charge)=(Ntp->PFTau_Charge(HighestPtTauIndex)*Ntp->Muon_Charge(HighestPtMuonIndex) == -1);


  }

  double wobs=1;
  double w;

  double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
  double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
  double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);


  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}
  bool status=AnalysisCuts(t,w,wobs); 
 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){

    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
    if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){

    TLorentzVector tau_3pi = Ntp->KFTau_TauFit_p4(HighestPtTauIndex);
    TLorentzVector muon    = Ntp->Muons_p4(HighestPtMuonIndex);

    TLorentzVector tau1,tau2,taumu1,taumu2,Z;
    double e1,e2,pz1,pz2;

    TLorentzVector z1,z2;

    double zmass = 91.1876,mtau = 1.777;

    double A =  0.5*(zmass*zmass - 2*mtau*mtau) -tau_3pi.Pt()* tau_3pi.Pt() ;
    double B =  mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt() + A*A/tau_3pi.Pz()/tau_3pi.Pz();
    double root = sqrt(A*A*tau_3pi.E()*tau_3pi.E() - B*tau_3pi.Pz()*tau_3pi.Pz()*(mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()));



    e1 = (A*tau_3pi.E() + root)/(mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt());
    e2 = (A*tau_3pi.E() - root)/(mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt());

    pz1 = sqrt(e1*e1 -  (mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()));
    pz2 = sqrt(e2*e2 -  (mtau*mtau + tau_3pi.Pt()* tau_3pi.Pt()));
    
    tau1.SetXYZM(-tau_3pi.Px(),-tau_3pi.Py(),pz1,mtau); 
    tau2.SetXYZM(-tau_3pi.Px(),-tau_3pi.Py(),pz2,mtau); 

    z1 = tau1+tau_3pi;
    z2 = tau2+tau_3pi;
    Tau3PiPt.at(t).Fill(tau_3pi.Pt(),w);
    Tau1Pt.at(t).Fill(tau1.Pt(),w);
    Tau2Pt.at(t).Fill(tau2.Pt(),w);
    Tau1theta.at(t).Fill(tau1.Theta(),w);	    
    Tau2theta.at(t).Fill(tau2.Theta(),w);      
    Tau1Deltatheta.at(t).Fill(tau1.Theta() - muon.Theta(),w); 
    Tau2Deltatheta.at(t).Fill(tau2.Theta() - muon.Theta(),w); 
    
    if(fabs(tau1.Theta() - muon.Theta()) < fabs(tau2.Theta() - muon.Theta())){taumu1 = tau1; taumu2 = tau2;}
    else{taumu1 = tau2; taumu2 = tau1;}
    DeltaTheta_1st.at(t).Fill(taumu1.Theta() - muon.Theta(),w); 
    DeltaTheta_2nd.at(t).Fill(taumu2.Theta() - muon.Theta(),w); 
    //    if(taumu1.E()!=0 && taumu1.Et()!=0){
    xmuon_plus.at(t).Fill(muon.E()/taumu1.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
    xmuonT_plus.at(t).Fill(muon.Et()/taumu1.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 


    xmuon_minus.at(t).Fill(muon.E()/taumu1.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
    xmuonT_minus.at(t).Fill(muon.Et()/taumu1.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 


    
    //------------------------------- comment out for work with data

    for(int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
      TLorentzVector TruthTauMu;
      TLorentzVector TruthTauPi;
      
      TLorentzVector TruthMu;
      TLorentzVector RecoiledTauMu;
      TLorentzVector RecoiledZ;
      int TauMuIndex =0;
      
      double e1sim,e2sim,pz1sim,pz2sim;
      TLorentzVector z1sim,z2sim;
      TLorentzVector tau1sim,tau2sim,taumu1sim,taumu2sim,Zsim;
      
      
      bool signal =false;
      if(Ntp->MCTau_JAK(0) == 2 and Ntp->MCTau_JAK(1) ==5 ){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); signal =true;TauMuIndex = 0;}
      else if(Ntp->MCTau_JAK(0) == 5 and Ntp->MCTau_JAK(1) ==2){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); signal =true;TauMuIndex = 1;}
      
      
      if(signal){
	for(int iProd1 =0; iProd1 < Ntp->NMCTauDecayProducts(TauMuIndex); iProd1++ ){
	 	  if(abs( Ntp->MCTauandProd_pdgid(TauMuIndex,iProd1))==13){
		    TruthMu = Ntp->MCTauandProd_p4(TauMuIndex,iProd1);
		  }
	}
	    
	    double A1 =  0.5*(zmass*zmass - 2*mtau*mtau) -TruthTauPi.Pt()*TruthTauPi.Pt() ;
	    double B1 =  mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt() + A1*A1/TruthTauPi.Pz()/TruthTauPi.Pz();
	    double root1 = sqrt(A1*A1*TruthTauPi.E()*TruthTauPi.E() - B1*TruthTauPi.Pz()*TruthTauPi.Pz()*(mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt()));
	    
	    e1sim = (A1*TruthTauPi.E() + root1)/(mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt());
	    e2sim = (A1*TruthTauPi.E() - root1)/(mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt());
	    
	    pz1sim = sqrt(e1sim*e1sim -  (mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt()));
	    pz2sim = sqrt(e2sim*e2sim -  (mtau*mtau + TruthTauPi.Pt()* TruthTauPi.Pt()));
	 
	 
	    tau1sim.SetXYZM(-TruthTauPi.Px(),-TruthTauPi.Py(),pz1sim,mtau); 
	    tau2sim.SetXYZM(-TruthTauPi.Px(),-TruthTauPi.Py(),pz2sim,mtau); 

	    z1sim = tau1sim+TruthTauPi;
	    z2sim = tau2sim+TruthTauPi;
	    if(fabs(tau1sim.Theta() - TruthMu.Theta()) < fabs(tau2sim.Theta() - TruthMu.Theta())){taumu1sim = tau1sim; taumu2sim = tau2sim;}
	    else{taumu1sim = tau2sim; taumu2sim = tau1sim;}


	    xmuonTruthT_plus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	    xmuonTruth_plus.at(t).Fill(TruthMu.E()/taumu1sim.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
  
	    xmuonTruthT_minus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	    xmuonTruth_minus.at(t).Fill(TruthMu.E()/taumu1sim.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 



	    xmuonMC_plus.at(t).Fill(TruthMu.E()/TruthTauMu.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
	    xmuonTMC_plus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT); 
  
	    xmuonMC_minus.at(t).Fill(TruthMu.E()/TruthTauMu.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 
	    xmuonTMC_minus.at(t).Fill(TruthMu.Et()/TruthTauMu.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT); 

	    TruthResolution.at(t).Fill(taumu1sim.E() - TruthTauMu.E(),w); 
	    TruthDeltaTheta.at(t).Fill(taumu1sim.Theta() - TruthMu.Theta(),w); 


      }
      
      if(signal){
	Energy_resolution1.at(t).Fill(taumu1.E() -TruthTauMu.E(),w); 
	Energy_resolution2.at(t).Fill(taumu2.E() -TruthTauMu.E(),w); 
	TransverseEnergy_resolution.at(t).Fill(taumu1.Et() -TruthTauMu.Et(),w); 
      }
    
      
    }

    }
    
    //      //------------------------------- comment out for work with data
    
  }
  
}
 




 

