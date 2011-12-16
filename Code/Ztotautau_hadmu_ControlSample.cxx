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
    if(i==NMuons)             cut.at(NMuons)=1;
    if(i==NMuonswithOverLap)  cut.at(NMuonswithOverLap)=1;
    if(i==NTau)               cut.at(NTau)=1;
    if(i==TriggerMatch)       cut.at(TriggerMatch)=1;
    if(i==dphiMuonTau)        cut.at(dphiMuonTau)=1;
    if(i==dpocaMuonTau)       cut.at(GoodVertexMatch)=1;
    if(i==ZMass)              cut.at(ZMass)=1;
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
    else if(i==NMuons){
      title.at(i)="Number of Muon $(N>$";
      title.at(i)+=cut.at(NMuons);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==NMuonswithOverLap){
      title.at(i)="Number of Muon with Overlap Removal$(N>$";
      title.at(i)+=cut.at(NMuonswithOverLap);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Muon with Overlap Removal";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==NTau){
      title.at(i)="Number of Taus $(N>$";
      title.at(i)+=cut.at(NTau);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==TriggerMatch){
      title.at(i)="Trigger Matched";
      hlabel="Trigger Matched";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatched_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatched_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==dphiMuonTau){
      title.at(i)="$d\\phi(\\tau,\\mu)>$";
      title.at(i)+=cut.at(dphiMuonTau);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="d#phi(#tau,#mu) (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,32,-TMath::Pi(),TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,32,-TMath::Pi(),TMath::Pi(),hlabel,"Events"));
    }
    if(i==dpocaMuonTau){
      title.at(i)="$dL(\\tau,\\mu)>$";
      title.at(i)+=cut.at(dphiMuonTau);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="dL(#tau,#mu) (mm)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,100,-25,25,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,100,-25,25,hlabel,"Events"));
    }
    else if(i==GoodVertexMatch){
      title.at(i)="Good Vertex Matched";
      hlabel="Good Vertex Matched";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatched_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatched_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    if(i==ZMass){
      title.at(i)="$M_{Z}>$";
      title.at(i)+=cut.at(ZMass);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{Z} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,40,0,200,hlabel,"Events"));
    }
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
  
  double wobs=1;
  double w=Ntp->EvtWeight3D();
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
  }
}




