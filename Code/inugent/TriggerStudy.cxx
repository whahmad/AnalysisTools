#include "TriggerStudy.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"

TriggerStudy::TriggerStudy(TString Name_, TString id_):
  Selection(Name_,id_)
{
  //  verbose=true;
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
    if(i==isSignal)          cut.at(isSignal)=1;
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
  
    if(i==isSignal){
      title.at(i)="isSignal";
      htitle=title.at(i);
      hlabel="isSignal";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isSignal_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isSignal_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  for(unsigned int h=0;h<4;h++){
    TString NP;
    if(h==0)NP="All";
    if(h==1)NP="1prong1prong";
    if(h==2)NP="1prong3prong";
    if(h==3)NP="3prong3prong";
    N_pt.push_back(HConfig.GetTH1D(Name+NP+"_N_pt","N_pt",40,0,200,"P_{T} (GeV)","Events"));
    n_pt.push_back(HConfig.GetTH1D(Name+NP+"_n_pt","n_pt",40,0,200,"P_{T} (GeV)","Events"));
    Eff_pt.push_back(HConfig.GetTH1D(Name+NP+"_Eff_pt","Eff_pt",40,0,200,"P_{T} (GeV)","Efficiency"));
    N_eta.push_back(HConfig.GetTH1D(Name+NP+"_N_eta","N_eta",40,-2.0,2.0,"#eta","Events"));
    n_eta.push_back(HConfig.GetTH1D(Name+NP+"_n_eta","n_eta",40,-2.0,2.0,"#eta","Events"));
    Eff_eta.push_back(HConfig.GetTH1D(Name+NP+"_Eff_eta","Eff_eta",40,-2.0,2.0,"#eta","Efficiency"));
  }

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  TriggerStudy::Store_ExtraDist(){
  for(unsigned int h=0;h<Eff_pt.size();h++){
    Extradist1d.push_back(&N_pt.at(h));
    Extradist1d.push_back(&N_eta.at(h));
    Extradist1d.push_back(&n_pt.at(h));
    Extradist1d.push_back(&n_eta.at(h));
    Extradist1d.push_back(&Eff_pt.at(h));
    Extradist1d.push_back(&Eff_eta.at(h));
  }
}

void  TriggerStudy::doEvent(){
  unsigned int t;
  int id(2);
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << id <<std::endl; return;}

  if(verbose)std::cout << "void  TriggerStudy::doEvent() A" << std::endl;

  // cout << "ntaus" <<Ntp->NMCTaus() << endl;
  value.at(isSignal)=0;
  if(Ntp->NMCTaus()==2){
    //cout<< Ntp->MCTau_JAK(0) <<  " " << Ntp->MCTau_JAK(1) << endl;
    if(Ntp->MCTau_JAK(0)!=TauDecay::JAK_ELECTRON && Ntp->MCTau_JAK(1)!=TauDecay::JAK_ELECTRON &&
       Ntp->MCTau_JAK(0)!=TauDecay::JAK_MUON && Ntp->MCTau_JAK(1)!=TauDecay::JAK_MUON){
      pass.at(isSignal)=true;
      value.at(isSignal)=1;
    }
  }
  for(unsigned p=0;p<Ntp->NHLTTriggers();p++){
    cout << Ntp->HTLTriggerName(p) <<  " Prescale: " << Ntp->HLTPrescale(p) << " " << Ntp->L1SEEDPrescale(p) <<endl;
  }

  //////////////////////////////////////////////////////////
  double wobs(1),w(1);
  if(!Ntp->isData()){
    if(verbose)std::cout << "void  TriggerStudy::doEvent() J" << std::endl;
    //w*=Ntp->EvtWeight3D();
    if(verbose)std::cout << "void  TriggerStudy::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  if(verbose)std::cout << "void  TriggerStud::doEvent() L" << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    TauDecay TD;
    int idx=0;
    unsigned int taubit1(Ntp->MCTau_DecayBitMask(0)), taubit2(Ntp->MCTau_DecayBitMask(1));
    if( TD.nProng(taubit1)==1 && TD.nProng(taubit2)==1){
      idx=1;
    }
    else if(TD.nProng(taubit1)==1 && TD.nProng(taubit2)==3 ||
	    TD.nProng(taubit1)==3 && TD.nProng(taubit2)==1){
      idx=2;
    }
    else if( TD.nProng(taubit1)==3 && TD.nProng(taubit2)==3){
      idx=3;
    }
    bool flag=Ntp->TriggerAccept("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v2");
    for(unsigned int tau=0;tau<Ntp->NMCTaus();tau++){
      double tau_pt=Ntp->MCTau_p4(tau).Pt();
      double tau_eta=Ntp->MCTau_p4(tau).Eta(); 
      N_pt.at(0).at(t).Fill(tau_pt,w);
      if(flag)n_pt.at(0).at(t).Fill(tau_pt,w);
      N_eta.at(0).at(t).Fill(tau_eta,w);
      if(flag)n_eta.at(0).at(t).Fill(tau_eta,w);
      if(idx!=0){
	N_pt.at(idx).at(t).Fill(tau_pt,w);
	if(flag)n_pt.at(idx).at(t).Fill(tau_pt,w);
	N_eta.at(idx).at(t).Fill(tau_eta,w);
	if(flag)n_eta.at(idx).at(t).Fill(tau_eta,w);
      }
    }
  }
}


void TriggerStudy::Finish(){
  for(unsigned int i=0;i<N_pt.size();i++){
    for(unsigned int j=0;j<N_pt.at(i).size();j++){
      Eff_pt.at(i).at(j).Reset();
      Eff_pt.at(i).at(j).Divide(&n_pt.at(i).at(j),&N_pt.at(i).at(j),1.0,1.0,"B");
      Eff_eta.at(i).at(j).Reset();
      Eff_eta.at(i).at(j).Divide(&n_eta.at(i).at(j),&N_eta.at(i).at(j),1.0,1.0,"B");
    }
  }
  Selection::Finish();
}

