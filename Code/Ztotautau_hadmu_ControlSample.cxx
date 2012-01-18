#include "Ztotautau_hadmu_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

Ztotautau_hadmu_ControlSample::Ztotautau_hadmu_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

Ztotautau_hadmu_ControlSample::~Ztotautau_hadmu_ControlSample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Ztotautau_hadmu_ControlSample::~Ztotautau_hadmu_ControlSample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Ztotautau_hadmu_ControlSample::~Ztotautau_hadmu_ControlSample()" << std::endl;
}

void  Ztotautau_hadmu_ControlSample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==MuonisGlob)         cut.at(MuonisGlob)=1;
    if(i==MuonPt)             cut.at(MuonPt)=25; //18
    if(i==TauPt)              cut.at(TauPt)=25;
    if(i==TauIsRef)           cut.at(TauIsRef)=1;
    if(i==MuonIso)            cut.at(MuonIso)=0.2;
    if(i==TauIsIso)           cut.at(TauIsIso)=1;
    if(i==deltaPhi)           cut.at(deltaPhi)=1;
    if(i==ZMassV)             cut.at(ZMassV)=1;
    if(i==ZPt)                cut.at(ZPt)=1;
    if(i==ZMassHPS)           cut.at(ZMassHPS)=1;
    if(i==MET)                cut.at(MET)=25;
    if(i==tauPhi)             cut.at(tauPhi)=1;
    if(i==charge)             cut.at(charge)=-1;
    if(i==decayMode)          cut.at(decayMode)=-1;



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
    //-----------
   else if(i==MuonisGlob){
      title.at(i)="$MuonIsGlob == $";
      title.at(i)+=cut.at(MuonisGlob);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MuonIsGlobal (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonisGlob_",htitle,2,-0.5,1.5,hlabel,"Events"));
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
   else if(i==MuonIso){
      title.at(i)="$MuonIso < $";
      title.at(i)+=cut.at(MuonIso);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#MuonIso";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIso_",htitle,40,0,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIso_",htitle,40,0,1,hlabel,"Events"));
    }

   else  if(i==TauIsRef){
      title.at(i)="$TauIsRef = $";
      title.at(i)+=cut.at(TauIsRef);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#TauIsRef (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsRef_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsRef_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==TauPt){
      title.at(i)="$TauPt > $";
      title.at(i)+=cut.at(TauPt);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#tau_{pT}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,30,0,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,30,0,50,hlabel,"Events"));
    }
   else if(i==ZPt){
      title.at(i)="$ZPt > $";
      title.at(i)+=cut.at(ZPt);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Z_{pT}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZPt_",htitle,30,0,40,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZPt_",htitle,30,0,40,hlabel,"Events"));
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
   else if(i==MET){
      title.at(i)="$MET < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="MET";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,20,0,60,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,20,0,60,hlabel,"Events"));
    } 
   else if(i==deltaPhi){
      title.at(i)="$DeltaPhi < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="deltaPhi";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,20,0,6,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,20,0,6,hlabel,"Events"));
    } 
   else if(i==ZMassV){
      title.at(i)="$ZMassV  $";
      title.at(i)+=cut.at(ZMassV);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="ZMassV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassV_",htitle,20,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassV_",htitle,20,0,200,hlabel,"Events"));
    } 
   else if(i==ZMassHPS){
      title.at(i)="$ZMassHPS  $";
      title.at(i)+=cut.at(ZMassHPS);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="ZMassHPS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassHPS_",htitle,20,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassHPS_",htitle,20,0,200,hlabel,"Events"));
    } 



   else if(i==tauPhi){
      title.at(i)="$tauPhis  $";
      title.at(i)+=cut.at(tauPhi);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="tauPhi";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_tauPhi_",htitle,20,-3,3,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_tauPhi_",htitle,20,-3,3,hlabel,"Events"));
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

   else if(i==decayMode){
      title.at(i)="$opposite charge$";
      title.at(i)+=cut.at(decayMode);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="decayMode";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_decayMode_",htitle,2000000,10130500,12130500,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_decayMode_",htitle,2000000,10130500,12130500,hlabel,"Events"));
    } 


    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Accumulative Cuts Passed","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  Ztotautau_hadmu_ControlSample::Store_ExtraDist(){

 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
}

void  Ztotautau_hadmu_ControlSample::doEvent(){
 
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
  
 

    unsigned int HighestPtMuonIndex=0;
    unsigned int SecondPtMuonIndex=0;
    float muonPt=0;
    int NGlMuons=0;
    int NQTaus =0;
    unsigned int HighestPtTauIndex=0;
    float tauPt=0;
    unsigned int iTau=0;

    for(iTau=0;iTau < Ntp->NKFTau();iTau++){
	if(Ntp->KFTau_TauFit_p4(iTau).Pt() > tauPt){
	  tauPt = Ntp->KFTau_TauFit_p4(iTau).Pt();
	  HighestPtTauIndex = iTau; 
	}
    }

  if(Ntp->NMuons()!=0){
    unsigned int iMuon=0;
   
    std::cout<<Ntp->NMuons() <<std::endl;

    for(iMuon=0; iMuon<Ntp->NMuons(); iMuon++){
    
      if(Ntp->Muon_isGlobalMuon(iMuon)){  	NGlMuons++;       }
      
      if(Ntp->Muons_p4(iMuon).Pt() > muonPt){
	muonPt = Ntp->Muons_p4(iMuon).Pt();
	SecondPtMuonIndex = HighestPtMuonIndex;
	HighestPtMuonIndex = iMuon;
      }
    }
    

  }


 if(Ntp->NMuons()!=0){
  value.at(MuonisGlob) = Ntp->Muon_isGlobalMuon(HighestPtMuonIndex);
  pass.at(MuonisGlob)=(value.at(MuonisGlob)==1);

  value.at(MuonPt) = Ntp->Muons_p4(HighestPtMuonIndex).Pt();
  pass.at(MuonPt)=(value.at(MuonPt) >=cut.at(MuonPt));

  value.at(MuonIso) = (Ntp->Muon_emEt05(HighestPtMuonIndex) + Ntp->Muon_hadEt05(HighestPtMuonIndex) + Ntp->Muon_sumPt05(HighestPtMuonIndex))/Ntp->Muons_p4(HighestPtMuonIndex).Pt();
  pass.at(MuonIso)=(value.at(MuonIso)<=cut.at(MuonIso));

 }
 


  value.at(TauIsRef) = Ntp->KFTau_discriminatorByQC(HighestPtTauIndex);
  pass.at(TauIsRef)=(value.at(TauIsRef) == 1);

  value.at(TauIsIso) =   Ntp->PFTau_isLooseIsolation(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex));
  pass.at(TauIsIso)=(value.at(TauIsIso) ==cut.at(TauIsIso));


 
  value.at(TauPt) = Ntp->KFTau_TauFit_p4(HighestPtTauIndex).Pt();
  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));

  //  value.at(TauPt) = Ntp->KFTau_TauFit_p4(HighestPtTauIndex).Pt();
  //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
  //pass.at(TauPt)=true;



  if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){
    TLorentzVector TauVis = Ntp->KFTau_TauVis_p4(HighestPtTauIndex);
    TLorentzVector TauFit = Ntp->KFTau_TauFit_p4(HighestPtTauIndex);
    TLorentzVector HPSMom = Ntp->PFTau_p4(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex));
    TLorentzVector Neutrino= Ntp->KFTau_Neutrino_p4(HighestPtTauIndex);
    TLorentzVector MuoGlo = Ntp->Muons_p4(HighestPtMuonIndex);
    TLorentzVector ZVis = TauVis + MuoGlo;
    TLorentzVector ZFit = TauFit + MuoGlo;
    TLorentzVector ZHPS = HPSMom + MuoGlo;


    value.at(MET) =   sqrt(2*Ntp->Muons_p4(HighestPtMuonIndex).Pt()*Ntp->MET_et()*(1  - cos(Ntp->Muons_p4(HighestPtMuonIndex).Phi() - Ntp->MET_phi())) );
    pass.at(MET)=(value.at(MET)<=cut.at(MET));

  
    value.at(deltaPhi) = fabs(MuoGlo.Phi() - HPSMom.Phi());
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(deltaPhi)=true;
    
    value.at(ZMassV) = ZVis.M();
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(ZMassV)=true;

    value.at(ZMassHPS) = ZHPS.M();
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(ZMassHPS)=true;
 
    value.at(ZPt) = fabs(MuoGlo.Pt() - TauFit.Pt());
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(ZPt)=true;
   
    value.at(tauPhi) = TauFit.Phi();
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(tauPhi)=true;

    value.at(decayMode) = Ntp->GetMCID();
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(decayMode)=true;
    

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
    std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
  }
}




