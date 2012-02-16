#include "Momentum_calculation_using_muon_angle.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"

Momentum_calculation_using_muon_angle::Momentum_calculation_using_muon_angle(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

Momentum_calculation_using_muon_angle::~Momentum_calculation_using_muon_angle(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Momentum_calculation_using_muon_angle::~Momentum_calculation_using_muon_angle Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Momentum_calculation_using_muon_angle::~Momentum_calculation_using_muon_angle()" << std::endl;
}

void  Momentum_calculation_using_muon_angle::Configure(){
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
  ZTruthPt=HConfig.GetTH1D(Name+"_ZTruthPt","ZTruthPt",30,0,100,"ZtruthPt","Events");
  TauPt =HConfig.GetTH1D(Name+"_TauPt","TauPt",20,0,80,"TauPt","Events");
  Tautheta=HConfig.GetTH1D(Name+"_Tautheta","Tautheta",20,0,3.14,"Tautheta","Events");
  Energy_resolution=HConfig.GetTH1D(Name+"_Energy_resolution","Energy_resolution",30,-100,100,"Energy_resolution","Events");
  TransverseEnergy_resolution=HConfig.GetTH1D(Name+"_TransverseEnergy_resolution","TransverseEnergy_resolution",30,-100,100,"TransverseEnergy_resolution","Events");
  Energy_resolution_cuts=HConfig.GetTH1D(Name+"_Energy_resolution_cuts","Energy_resolution_cuts",30,-100,100,"Energy_resolution_cuts","Events");
  ResVsZPt=HConfig.GetTH2D(Name+"_ResVsZPt","ResVsZPt",20,0,40,10,0,0.5,"ResVsZPt","Events");
  TauMuEnergyRatio=HConfig.GetTH1D(Name+"_TauMuEnergyRatio","TauMuEnergyRatio",10,0,1,"TauMuEnergyRatio","Events");
  TauMuTransverseEnergyRatio=HConfig.GetTH1D(Name+"_TauMuTransverseEnergyRatio","TauMuTransverseEnergyRatio",10,0,2,"TauMuTransverseEnergyRatio","Events");
  ZMass =HConfig.GetTH1D(Name+"_ZMass","ZMass",20,40,160,"ZMass","Events");
  resol=HConfig.GetTH1D(Name+"_resol","resol",15,40,160,"resol","Events");
  TruthRecoiledZ=HConfig.GetTH1D(Name+"_TruthRecoiledZ","TruthRecoiledZ",20,0,100,"TruthRecoiledZ","Events");

  TruthRatio=HConfig.GetTH1D(Name+"_TruthRatio","TruthRatio",10,0,1,"TruthRatio","Events");
  TruthTransverseRatio=HConfig.GetTH1D(Name+"_TruthTransverseRatio","TruthTransverseRatio",10,0,1,"TruthTransverseRatio","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
} 


// TString name,TString title,int nbinsx, double minx, double maxx, 
// 				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis)

void  Momentum_calculation_using_muon_angle::Store_ExtraDist(){


 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&ZTruthPt);
 Extradist1d.push_back(&TauPt);
 Extradist1d.push_back(&Tautheta);	    
 Extradist1d.push_back(&Energy_resolution);
 Extradist1d.push_back(&Energy_resolution_cuts);
 Extradist1d.push_back(&TransverseEnergy_resolution);
 Extradist2d.push_back(&ResVsZPt);
 Extradist1d.push_back(&TauMuEnergyRatio);
 Extradist1d.push_back(&TauMuTransverseEnergyRatio);
 Extradist1d.push_back(&ZMass);
  Extradist1d.push_back(&TruthRecoiledZ);
  Extradist1d.push_back(&resol);


Extradist1d.push_back(&TruthRatio);
Extradist1d.push_back(&TruthTransverseRatio);

}

void  Momentum_calculation_using_muon_angle::doEvent(){
 
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
   std::cout<<" Momentum Calculation Using Muon angle "<<std::endl;
  std::cout<<"DataMCType   ==================++>" << Ntp->GetMCID()<<std::endl;
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

    if( Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex)!=-1 && Ntp->NTracks() >= Ntp->Muon_Track_idx(HighestPtMuonIndex)){
      value.at(charge) =Ntp->KFTau_Fit_charge(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex))*Ntp->Track_charge(Ntp->Muon_Track_idx(HighestPtMuonIndex));

      //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
      pass.at(charge)=(Ntp->KFTau_Fit_charge(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex))*Ntp->Track_charge(Ntp->Muon_Track_idx(HighestPtMuonIndex)) == -1);
    }
    

  }

  double wobs=1;
  double w;







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
    if(sin(theta)!=0){
      energy = sqrt(mtau*mtau*sin(theta)*sin(theta) + tau_3pi.Pt()*tau_3pi.Pt())/sin(theta);
      std::cout<<"sin(30) " << sin(30) << "sin(30*3.14/180) "<<sin(30*3.14/180)<<std::endl;
    }else energy = -1;
    tau.SetXYZM(-tau_3pi.Px(), -tau_3pi.Py(),sqrt(energy*energy - mtau*mtau)*cos(theta),mtau);
    Z  = tau_3pi + tau;
    ZMass.at(t).Fill(Z.M(),w); 


    TauMuEnergyRatio.at(t).Fill(muon.E()/tau.E(),w); 
    TauMuTransverseEnergyRatio.at(t).Fill(muon.Et()/tau.Et(),w); 

    //------------------------------- comment out for work with data
//      std::cout<<"---------------------------------- " <<std::endl;
//      printf("energy1 %f   first tau Pt %f\n",tau.E(),tau_3pi.E());
//      std::cout<<"Z pt " <<Z.Pt() << "tau PT  "<<tau.Pt() <<std::endl;




//  for(int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
//    TLorentzVector TruthTauMu;
//    TLorentzVector TruthTauPi;
//    TLorentzVector TruthMu;
//    TLorentzVector RecoiledTauMu;
//    TLorentzVector RecoiledZ;
//    int TauMuIndex =0;



//    bool signal =false;
//    if(Ntp->MCTau_JAK(0) == 2 and Ntp->MCTau_JAK(1) ==5 ){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); signal =true;TauMuIndex = 0;}
//    else if(Ntp->MCTau_JAK(0) == 5 and Ntp->MCTau_JAK(1) ==2){TruthTauMu = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(1)); TruthTauPi = Ntp->MCTau_p4(Ntp->MCSignalParticle_Tauidx(iz).at(0)); signal =true;TauMuIndex = 1;}




//    if(signal){
//      Energy_resolution.at(t).Fill((tau.E() -TruthTauMu.E()),w); 
//      TransverseEnergy_resolution.at(t).Fill((tau.Et() -TruthTauMu.Et()),w); 

//    }


//    if(signal and Ntp->MCSignalParticle_p4(iz).Pt() < 10 ){Energy_resolution_cuts.at(t).Fill(tau.E() -TruthTauMu.E(),w);  }


//    if(signal){
//      for(int iProd1 =0; iProd1 < Ntp->NMCTauDecayProducts(TauMuIndex); iProd1++ ){
//        std::cout<<"--pdgId 1  " << Ntp->MCTauandProd_pdgid(TauMuIndex,iProd1) <<std::endl;
//        if(abs( Ntp->MCTauandProd_pdgid(TauMuIndex,iProd1))==14){
// 	 TruthMu = Ntp->MCTauandProd_p4(TauMuIndex,iProd1);
// 	 	 if(TruthMu.Theta()!=0){
// 		   RecoiledTauMu.SetXYZM(-TruthTauPi.Px(), -TruthTauPi.Py(), sqrt(   mtau*mtau*sin(TruthMu.Theta())*sin(TruthMu.Theta()) + TruthTauPi.Pt()*TruthTauPi.Pt())/sin(TruthMu.Theta())  ,mtau);
// 		   RecoiledZ = RecoiledTauMu + TruthTauPi;
// 		   std::cout<<"-----RecoiledZ  " << RecoiledZ.M() <<std::endl;

// 		 }
		 
//  		 TruthRecoiledZ.at(t).Fill(RecoiledTauMu.E() -TruthTauMu.E(),w); 
//  		 resol.at(t).Fill(RecoiledZ.M(),w); 
// 		 TruthRatio.at(t).Fill(TruthMu.E()/RecoiledTauMu.E(),w); 
// 		 TruthTransverseRatio.at(t).Fill(TruthMu.Et()/RecoiledTauMu.Et(),w); 

//        }


//      }
//    }

//    ZTruthPt.at(t).Fill(Ntp->MCSignalParticle_p4(iz).Pt(),w);


//  }

//      //------------------------------- comment out for work with data

    }

  }
} 




 
