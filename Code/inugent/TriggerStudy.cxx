#include "TriggerStudy.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "TLorentzVector.h"

TriggerStudy::TriggerStudy(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

TriggerStudy::~TriggerStudy(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TriggerStudy::~TriggerStudy Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TriggerStudy::~TriggerStudy()" << std::endl;
}

void  TriggerStudy::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    cut.at(TriggerOk)=1;
    if(i==NTau)         cut.at(NTau)=2;
    if(i==TriggerMatch) cut.at(TriggerMatch)=1;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    Accumdist.push_back(std::vector<TH1D>());
    Nminus1dist.push_back(std::vector<TH1D>());

    TString c="_Cut_";
    if(i<10)c+="0";
    c+=i;
    if(i==TriggerOk){
      title.at(i)="TriggerOk";
      htitle=title.at(i);
      hlabel="TriggerOk";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    if(i==NTau){
      title.at(i)="NTau";
      htitle=title.at(i);
      hlabel="NTau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTau_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTau_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    if(i==TriggerMatch){
      title.at(i)="TriggerMatch";
      htitle=title.at(i);
      hlabel="TriggerMatch";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerMatch_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerMatch_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  for(unsigned int h=0;h<2;h++){
    TString NP;
    if(h==0)NP="DataParking";
    if(h==1)NP="1ProngTrigger";
    N_pt.push_back(HConfig.GetTH1D(Name+NP+"_N_pt","N_pt",20,0,200,"P_{T} (GeV)","Number of Trigger Obj. (#taus)"));
    N_eta.push_back(HConfig.GetTH1D(Name+NP+"_N_eta","N_eta",20,-5.0,5.0,"#eta","Number of Trigger Obj. (#taus)"));
    N_2d.push_back(HConfig.GetTH2D(Name+NP+"_N_2d","N_2d",20,0,200,20,-5,5,"P_{T} (GeV)","#eta"));
    HpsMode.push_back(HConfig.GetTH1D(Name+NP+"_HpsMode","HpsMode",21,-0.5,20.5,"HPS Mode","Number of #taus"));
  }
  Mtautau=HConfig.GetTH1D(Name+"_Mtautau","Mtautau",20,0.0,200,"M_{#tau#tau} (GeV)","Number of Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  TriggerStudy::Store_ExtraDist(){
  for(unsigned int h=0;h<HpsMode.size();h++){
    Extradist1d.push_back(&N_pt.at(h));
    Extradist1d.push_back(&N_eta.at(h));
    Extradist1d.push_back(&HpsMode.at(h));
    Extradist2d.push_back(&N_2d.at(h));
  }
  Extradist1d.push_back(&Mtautau);
}

void  TriggerStudy::doEvent(){
  unsigned int t;
  int id(2);
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << id <<std::endl; return;}

  if(verbose)std::cout << "void  TriggerStudy::doEvent() A" << std::endl;

  TString DataParking="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_v";
  TString OneProng="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1_v"; 

  bool FlagDataParking=Ntp->TriggerAccept(DataParking);
  bool Flag1Prong=Ntp->TriggerAccept(OneProng);

  /*  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
    std::cout <<  Ntp->HTLTriggerName(i) << " " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << " "  << std::endl;
    }*/

  //////////////////////////////////////////////////////////////////////////
  //
  // Make sure only events with the 1prong or DataParking trigger are used
  //
  pass.at(TriggerOk)=(Flag1Prong || FlagDataParking);
  value.at(TriggerOk)=0;
  if(FlagDataParking) value.at(TriggerOk)+=1;
  if(Flag1Prong) value.at(TriggerOk)+=2;

  //////////////////////////////////////////////////////////////////////////                                                                                                     
  //
  // Find HPS Taus
  //
  std::vector<unsigned int> tau_idx;
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_isTightIsolationDBSumPtCorr(i) && Ntp->PFTau_isHPSAgainstElectronsTight(i) && Ntp->PFTau_isHPSAgainstMuonTight(i) && Ntp->PFTau_p4(i).Pt()>20){
      tau_idx.push_back(i);
    }
  }

  for(unsigned int i=0; i<tau_idx.size();i++){
    for(unsigned int j=i+1; j<tau_idx.size();j++){
      if(fabs(Ntp->PFTau_p4(tau_idx.at(i)).DeltaPhi(Ntp->PFTau_p4(tau_idx.at(j))))<TMath::Pi()*7.0/8.0) tau_idx.erase(tau_idx.begin()+j);
    }
  }
  value.at(NTau)=tau_idx.size();
  pass.at(NTau)=value.at(NTau)>=cut.at(NTau);
 
  //////////////////////////////////////////////////////////////////////////
  //
  // Check that HPS taus were the triggered taus
  //
  unsigned int tautrig_DataParking=0;
  unsigned int tautrig_OneProng=0;
  float drmatch=0.2;
  unsigned int A(0),B(0);
  for(unsigned int i=0; i<tau_idx.size();i++){
    if(Ntp->GetTriggerIndex(DataParking,tautrig_DataParking)){
      if(Ntp->MuonTriggerMatch(tautrig_DataParking,tau_idx.at(i))<drmatch)A++;
    }
    if(Ntp->GetTriggerIndex(OneProng,tautrig_OneProng)){
      if(Ntp->MuonTriggerMatch(tautrig_OneProng,tau_idx.at(i))<drmatch)B++;
    }
  }
  if(A==2 && B==2 && value.at(TriggerOk)==3) value.at(TriggerMatch)=3;
  else if(A==2 && (value.at(TriggerOk)==1 || value.at(TriggerOk)==3)) value.at(TriggerMatch)=1;
  else if(B==2 && (value.at(TriggerOk)==2 || value.at(TriggerOk)==3)) value.at(TriggerMatch)=2;
  pass.at(TriggerMatch)=true;
  //////////////////////////////////////////////////////////
  double wobs(1),w(1);
  if(!Ntp->isData()){

  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  if(verbose)std::cout << "void  TriggerStud::doEvent() L" << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status && tau_idx.size()>=2){
    TLorentzVector LV=Ntp->PFTau_p4(tau_idx.at(0));LV+=Ntp->PFTau_p4(tau_idx.at(1));
    Mtautau.at(t).Fill(LV.M(),w);
    for(unsigned int h=0; h<HpsMode.size();h++){
      for(unsigned int i=0; i<tau_idx.size();i++){
	if((h==0 && FlagDataParking) ||(h==1 && Flag1Prong)){
	  HpsMode.at(h).at(t).Fill(Ntp->PFTau_hpsDecayMode(tau_idx.at(i)),w);
	  unsigned int trig=0;
	  float tau_pt(0),tau_eta(0);
	  if(i==0){
	    unsigned int trig=tautrig_OneProng;
	    if(h==0)trig=tautrig_DataParking;
	    if(h==1)trig=tautrig_OneProng;
	    for(unsigned int jidx=0; jidx<Ntp->NHLTTriggerObject(trig);jidx++){
	      tau_pt=Ntp->HLTTriggerObject_p4(trig,jidx).Pt();
	      tau_eta=Ntp->HLTTriggerObject_p4(trig,jidx).Eta();
	      N_pt.at(h).at(t).Fill(tau_pt,w);
	      N_eta.at(h).at(t).Fill(tau_eta,w);
	      N_2d.at(h).at(t).Fill(tau_pt,tau_eta,w);  
	    }
	  }
	}
      }
    }
  }
}

void TriggerStudy::Finish(){
  Selection::Finish();
}

