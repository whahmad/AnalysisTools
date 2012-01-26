#include "Ztotautau_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

Ztotautau_ControlSample::Ztotautau_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

Ztotautau_ControlSample::~Ztotautau_ControlSample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Ztotautau_ControlSample::~Ztotautau_ControlSample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Ztotautau_ControlSample::~Ztotautau_ControlSample()" << std::endl;
}

void  Ztotautau_ControlSample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==hasTag)             cut.at(hasTag)=1;
    if(i==TagPt)              cut.at(TagPt)=18;
    if(i==TagIso)             cut.at(TagIso)=1;
    if(i==JetPt)              cut.at(JetPt)==17;
    if(i==deltaPhi)           cut.at(deltaPhi)==TMath::Pi()*3.0/4.0;
    if(i==MET)                cut.at(MET)=20;
    if(i==MT)                 cut.at(MT)=30;
    if(i==ZMassV)             cut.at(ZMassV)=1;
    if(i==charge)             cut.at(charge)=1;
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
   else if(i==hasTag){
      title.at(i)="has Tag";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="has Tag (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
   else if(i==TagPt){
      title.at(i)="$P_{T}^{Tag}>$";
      title.at(i)+=cut.at(TagPt);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="P_{T}^{Tag} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPt_",htitle,20,0,50,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPt_",htitle,20,0,50,hlabel,"Events"));
    }
   else if(i==TagIso){
      title.at(i)="Relative Isolation (Tag) $<$";
      title.at(i)+=cut.at(TagIso);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Relative Isolation (Tag)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagIso_",htitle,40,0,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagIso_",htitle,40,0,1,hlabel,"Events"));
    }
   else if(i==JetPt){
     title.at(i)="$P_{T}^{Jet}>$";
     title.at(i)+=cut.at(JetPt);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="P_{T}^{Tag} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPt_",htitle,20,0,50,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPt_",htitle,20,0,50,hlabel,"Events"));
   }
   else if(i==MET){
      title.at(i)="$E_{T}^{Miss} < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="E_{T}^{Miss} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,20,0,60,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,20,0,60,hlabel,"Events"));
    } 
   else if(i==MT){
     title.at(i)="$M_{T} < $";
     title.at(i)+=cut.at(MT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_",htitle,20,0,60,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_",htitle,20,0,60,hlabel,"Events"));
   }
   else if(i==deltaPhi){
      title.at(i)="$\\Delta\\phi(Tag,Jet) < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta#phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,20,0,6,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,20,0,6,hlabel,"Events"));
    } 
   else if(i==ZMassV){
      title.at(i)="$M_{Z} $";
      title.at(i)+=cut.at(ZMassV);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{Z} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassV_",htitle,20,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassV_",htitle,20,0,200,hlabel,"Events"));
    } 
   else if(i==charge){
      title.at(i)="Charge";
      title.at(i)+=cut.at(charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,3,-1.5,1.5,hlabel,"Events"));
    } 

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Accumulative Cuts Passed","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  PmuoverEtau=HConfig.GetTH1D(Name+"_PmuoverEtau","PmuoverEtau",100,0.0,1.0,"P_{#mu}/P_{#tau}","Events");
  PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_PmuoverEtau_hplus","PmuoverEtau_hplus",100,0.0,1.0,"P_{#mu}/P_{#tau}","Events");
  PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_PmuoverEtau_hminus","PmuoverEtau_hminus",100,0.0,1.0,"P_{#mu}/P_{#tau}","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  Ztotautau_ControlSample::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
}

void  Ztotautau_ControlSample::doEvent(){
 
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



  if(Ntp->NKFTau()!=0 && Ntp->NMuons()!=0){
    TLorentzVector TauVis = Ntp->KFTau_TauVis_p4(HighestPtTauIndex);
    TLorentzVector TauFit = Ntp->KFTau_TauFit_p4(HighestPtTauIndex);
    TLorentzVector HPSMom = Ntp->PFTau_p4(Ntp->KFTau_MatchedHPS_idx(HighestPtTauIndex));
    TLorentzVector Neutrino= Ntp->KFTau_Neutrino_p4(HighestPtTauIndex);
    TLorentzVector MuoGlo = Ntp->Muons_p4(HighestPtMuonIndex);
    TLorentzVector ZVis = TauVis + MuoGlo;
    TLorentzVector ZFit = TauFit + MuoGlo;
    TLorentzVector ZHPS = HPSMom + MuoGlo;


    

    if( Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex)!=-1 && Ntp->NTracks() >= Ntp->Muon_Track_idx(HighestPtMuonIndex)){
      value.at(charge) =Ntp->KFTau_Fit_charge(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex))*Ntp->Track_charge(Ntp->Muon_Track_idx(HighestPtMuonIndex));

      //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
      pass.at(charge)=(Ntp->KFTau_Fit_charge(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex))*Ntp->Track_charge(Ntp->Muon_Track_idx(HighestPtMuonIndex)) == -1);
    }
    
    
  }



  PmuoverEtau.at(t).Fill(1);
  PmuoverEtau_hplus.at(t).Fill(1);
  PmuoverEtau_hminus.at(t).Fill(1);

  

  double wobs=1;
  double w;

  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
    w*=Ntp->TauSpinerWeight(TauSpinerInterface::FlipSpin);
  }
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




