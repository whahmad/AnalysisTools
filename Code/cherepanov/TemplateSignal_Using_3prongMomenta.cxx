#include "TemplateSignal_Using_3prongMomenta.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"

TemplateSignal_Using_3prongMomenta::TemplateSignal_Using_3prongMomenta(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

TemplateSignal_Using_3prongMomenta::~TemplateSignal_Using_3prongMomenta(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TemplateSignal_Using_3prongMomenta::~TemplateSignal_Using_3prongMomenta Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TemplateSignal_Using_3prongMomenta::~TemplateSignal_Using_3prongMomenta()" << std::endl;
}

void  TemplateSignal_Using_3prongMomenta::Configure(){
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
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",20,0,80,"TauPt","Events");

  ResolMu=HConfig.GetTH1D(Name+"_ResolMu","ResolMu",30,-50,50,"ResolMu","Events");
  ResolMuTruth=HConfig.GetTH1D(Name+"_ResolMuTruth","ResolMuTruth",30,-50,50,"ResolMuTruth","Events");
  ResolMuT=HConfig.GetTH1D(Name+"_ResolMuT","ResolMuT",30,-50,50,"ResolMuT","Events");
  ResolMuTruthT=HConfig.GetTH1D(Name+"_ResolMuTruthT","ResolMuTruthT",30,-50,50,"ResolMuTruthT","Events");
  


  xmuon_plus=HConfig.GetTH1D(Name+"_xmuon_plus","xmuon_plus",15,0.2,1.4,"xmuon_plus","Events");
  xmuon_minus=HConfig.GetTH1D(Name+"_xmuon_minus","xmuon_minus",15,0.2,1.4,"xmuon_minus","Events");
  xmuonT_plus=HConfig.GetTH1D(Name+"_xmuonT_plus","xmuonT_plus",15,0,2,"xmuonT_plus","Events");
  xmuonT_minus=HConfig.GetTH1D(Name+"_xmuonT_minus","xmuonT_minus",15,0,2,"xmuonT_minus","Events");
 

  xmuon_plus_truth=HConfig.GetTH1D(Name+"_xmuon_plus_truth","xmuon_plus_truth",15,0.2,1.4,"xmuon_plus_truth","Events");
  xmuon_minus_truth=HConfig.GetTH1D(Name+"_xmuon_minus_truth","xmuon_minus_truth",15,0.2,1.4,"xmuon_minus_truth","Events");
  xmuonT_plus_truth=HConfig.GetTH1D(Name+"_xmuonT_plus_truth","xmuonT_plus_truth",15,0,2,"xmuonT_plus_truth","Events");
  xmuonT_minus_truth=HConfig.GetTH1D(Name+"_xmuonT_minus_truth","xmuonT_minus_truth",15,0,2,"xmuonT_minus_truth","Events");


  xmuon_plus_MC=HConfig.GetTH1D(Name+"_xmuon_plus_MC","xmuon_plus_MC",15,0,1,"xmuon_plus_MC","Events");
  xmuon_minus_MC=HConfig.GetTH1D(Name+"_xmuon_minus_MC","xmuon_minus_MC",15,0,1,"xmuon_minus_MC","Events");
  xmuonT_plus_MC=HConfig.GetTH1D(Name+"_xmuonT_plus_MC","xmuonT_plus_MC",15,0,1,"xmuonT_plus_MC","Events");
  xmuonT_minus_MC=HConfig.GetTH1D(Name+"_xmuonT_minus_MC","xmuonT_minus_MC",15,0,1,"xmuonT_minus_MC","Events");



  ZMass =HConfig.GetTH1D(Name+"_ZMass","ZMass",20,40,160,"ZMass","Events");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
} 


// TString name,TString title,int nbinsx, double minx, double maxx, 
// 				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis)

void  TemplateSignal_Using_3prongMomenta::Store_ExtraDist(){


 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&TauPt);
 Extradist1d.push_back(&ZMass);
 
 Extradist1d.push_back(&xmuon_plus);		
 Extradist1d.push_back(&xmuon_minus);	
 Extradist1d.push_back(&xmuon_plus_truth);	
 Extradist1d.push_back(&xmuon_minus_truth);	
 Extradist1d.push_back(&xmuon_plus_MC);	
 Extradist1d.push_back(&xmuon_minus_MC);    	
 
 Extradist1d.push_back(&xmuonT_plus);	
 Extradist1d.push_back(&xmuonT_minus);	
 Extradist1d.push_back(&xmuonT_plus_truth);	
 Extradist1d.push_back(&xmuonT_minus_truth);	
 Extradist1d.push_back(&xmuonT_plus_MC);	
 Extradist1d.push_back(&xmuonT_minus_MC);    
 Extradist1d.push_back(&ResolMu ); 
 Extradist1d.push_back(&ResolMuTruth); 
 Extradist1d.push_back(&ResolMuT); 	  
 Extradist1d.push_back(&ResolMuTruthT); 



}

void  TemplateSignal_Using_3prongMomenta::doEvent(){
 
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
  
  
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
         if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
    if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){

    TLorentzVector tau_3pi = Ntp->KFTau_TauFit_p4(HighestPtTauIndex);

    TLorentzVector muon    = Ntp->Muons_p4(HighestPtMuonIndex);

    TLorentzVector tau,Z;
    double zmass = 91.1876,mtau = 1.777;
    
    double energy;
    double theta = muon.Theta();

    double xmu;



    xmu = muon.Et()/tau_3pi.Et();




    ZMass.at(t).Fill(Z.M(),w); 
    xmuon_plus.at(t).Fill(xmu,w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
    xmuon_minus.at(t).Fill(xmu,w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);

    //------------------------------- comment out for work with data

  for(int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
    TLorentzVector TruthTauMu;
    TLorentzVector TruthTauPi;
    TLorentzVector TruthMu;

    TLorentzVector RecoiledTauMu;
    TLorentzVector RecoiledZ;
    TLorentzVector a1;
    double xmutruth;
    int TauMuIndex =0;
    int TauPiIndex =0;

    bool signal =false;
    if(Ntp->MCTau_JAK(0) == 2 and Ntp->MCTau_JAK(1) ==5 ){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); signal =true;TauMuIndex = 0;TauPiIndex=1;}
    else if(Ntp->MCTau_JAK(0) == 5 and Ntp->MCTau_JAK(1) ==2){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); signal =true;TauMuIndex = 1;TauPiIndex=0;}



    a1.SetXYZM(0,0,0,0);
    if(signal){
      for(int iProd1 =0; iProd1 < Ntp->NMCTauDecayProducts(TauMuIndex); iProd1++ ){
        if(abs( Ntp->MCTauandProd_pdgid(TauMuIndex,iProd1))==13){
	  TruthMu = Ntp->MCTauandProd_p4(TauMuIndex,iProd1);
	}
        if(abs( Ntp->MCTauandProd_pdgid(TauPiIndex,iProd1))==211){
	  a1+= Ntp->MCTauandProd_p4(TauPiIndex,iProd1);
	}

      }
      xmutruth = TruthMu.Et()/a1.Et();

      ResolMu.at(t).Fill(muon.E()/xmu - TruthTauMu.E(),w);
      ResolMuTruth.at(t).Fill(TruthMu.E()/xmutruth - TruthTauMu.E(),w);
      ResolMuT.at(t).Fill(muon.Et()/xmu - TruthTauMu.Et(),w);
      ResolMuTruthT.at(t).Fill(TruthMu.Et()/xmutruth - TruthTauMu.Et(),w);

      xmuon_plus_MC.at(t).Fill(TruthMu.Et()/TruthTauPi.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
      xmuon_minus_MC.at(t).Fill(TruthMu.Et()/TruthTauPi.Et(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);
    }

    
    
  }

//------------------------------- comment out for work with data

    }

  }
}




 
