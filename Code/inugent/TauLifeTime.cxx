#include "TauLifeTime.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSolver.h"
//#include "TCanvas.h"
//#include "NTuple_Controller.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "math.h"

TauLifeTime::TauLifeTime(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
{
  verbose=true;//false;
}

TauLifeTime::~TauLifeTime(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "TauLifeTime::~TauLifeTime Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "TauLifeTime::~TauLifeTime()" << std::endl;
}

void  TauLifeTime::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasTag)	      cut.at(hasTag)=1;
    if(i==TagPtMin)           cut.at(TagPtMin)=20;
    if(i==TagIso)             cut.at(TagIso)=0.2;
//    if(i==numIsoTags)         cut.at(numIsoTags)=1;
    if(i==TauPt)              cut.at(TauPt)=20;
    if(i==TauEta)             cut.at(TauEta)=2.0;
    if(i==TauIsIsolated)      cut.at(TauIsIsolated)=true;
    if(i==TauFit)             cut.at(TauFit)=1;
    if(i==deltaPhi)           cut.at(deltaPhi)=TMath::Pi()*3.0/4.0;
    if(i==Charge)             cut.at(Charge)=0.0;
    //if(i==ZMassmax)           cut.at(ZMassmax)=80;
  //std::cout << "Setting cut no. i=" << i << std::endl;
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
      title.at(i)="$N_{\\mu} passing the Trigger$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu} passing the Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }        
    else if(i==hasTag){
      title.at(i)="$Number of good muons$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TagPtMin){
      title.at(i)="$N_{\\mu} with P_{T}>$";
      title.at(i)+=cut.at(TagPtMin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtMin_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtMin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TagIso){
      title.at(i)="$N_{\\mu} with rel. isolation <=$";
      title.at(i)+=cut.at(TagIso);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauPt){
      title.at(i)="$N_{\\tau} with P_{T} >$";
      title.at(i)+=cut.at(TauPt);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauEta){
      title.at(i)="$N_{\\tau} | \\eta <$";
      title.at(i)+=cut.at(TauEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauIsIsolated){
      title.at(i)="$N_{\\tau} | isolated$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsolated_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsolated_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
else if(i==TauFit){
      title.at(i)="$N_{\\tau} | fitted good$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==deltaPhi){
      title.at(i)="$\\Delta \\phi distribution >$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta #phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
    }
    else if(i==Charge){
      title.at(i)="$|Q_{\\tau}+Q_{\\mu}|-0.5<$";
      title.at(i)+=cut.at(Charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="|Q_{#tau}+Q_{#mu}|";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Charge_",htitle,44,0,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Charge_",htitle,44,0,2.2,hlabel,"Events"));
    }
    //-----------
   /*else if(i==ZMassmax){
     title.at(i)="$M_{Z}<$";
     title.at(i)+=cut.at(ZMassmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{Z} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
   }*/
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  TauFlightLength=HConfig.GetTH1D(Name+"_TauFlightLength","TauFlightLength",62,-0.55,2.55,"L (cm)","Events");
  TauFlightLengthTransverse=HConfig.GetTH1D(Name+"_TauFlightLengthTransverse","TauFlightLengthTransverse",62,-0.55,2.55,"L_{T} (cm)","Events");
  TauMomentum=HConfig.GetTH1D(Name+"_TauMomentum","TauMomentum",50,0,400,"|P| (GeV)","Events");
  TauMomentumTransverse=HConfig.GetTH1D(Name+"_TauMomentumTransverse","TauMomentum",50,0,400,"|P_{T}| (GeV)","Events");
  TauLife=HConfig.GetTH1D(Name+"_TauLife","TauLife",1000,-0.5,2,"#tau_{#tau} (ps)","Events");
  TauLifeTransverse=HConfig.GetTH1D(Name+"_TauLifeTransverse","TauLifeTransverse",1000,-0.5,2,"transverse #tau_{#tau} (ps)","Events");
  ResTauFlightLength=HConfig.GetTH1D(Name+"_ResTauFlightLength","ResTauFlightLength",100,-2,2,"L (cm)","Events");
  ResTauFlightLengthTransverse=HConfig.GetTH1D(Name+"_ResTauFlightLengthTransverse","ResTauFlightLengthTransverse",100,-2,2,"L_{T} (cm)","Events");
  ResTauMomentum=HConfig.GetTH1D(Name+"_ResTauMomentum","ResTauMomentum",80,-200,200,"|P| (GeV)","Events");
  ResTauMomentumTransverse=HConfig.GetTH1D(Name+"_ResTauMomentumTransverse","ResTauMomentum",80,-200,200,"|P_{T}| (GeV)","Events");
  ResTauLife=HConfig.GetTH1D(Name+"_ResTauLife","ResTauLife",1000,-1,1,"#tau_{#tau} (ps)","Events");
  ResTauLifeTransverse=HConfig.GetTH1D(Name+"_ResTauLifeTransverse","ResTauLifeTransverse",1000,-1,1,"transverse #tau_{#tau} (ps)","Events");

  //
  TauMomentumAvg=HConfig.GetTH1D(Name+"_TauMomentumAvg","TauMomentumAvg",50,0,400,"|P| (GeV)","Events");
  TauMomentumTransverseAvg=HConfig.GetTH1D(Name+"_TauMomentumTransverseAvg","TauMomentumAvg",50,0,400,"|P_{T}| (GeV)","Events");
  TauLifeAvg=HConfig.GetTH1D(Name+"_TauLifeAvg","TauLifeAvg",1000,-0.5,2,"#tau_{#tau} (ps)","Events");
  TauLifeTransverseAvg=HConfig.GetTH1D(Name+"_TauLifeTransverseAvg","TauLifeTransverseAvg",1000,-0.5,2,"transverse #tau_{#tau} (ps)","Events");
  ResTauMomentumAvg=HConfig.GetTH1D(Name+"_ResTauMomentumAvg","ResTauMomentumAvg",80,-200,200,"|P| (GeV)","Events");
  ResTauMomentumTransverseAvg=HConfig.GetTH1D(Name+"_ResTauMomentumTransverseAvg","ResTauMomentumAvg",80,-200,200,"|P_{T}| (GeV)","Events");
  ResTauLifeAvg=HConfig.GetTH1D(Name+"_ResTauLifeAvg","ResTauLifeAvg",1000,-1,1,"#tau_{#tau} (ps)","Events");
  ResTauLifeTransverseAvg=HConfig.GetTH1D(Name+"_ResTauLifeTransverseAvg","ResTauLifeTransverseAvg",1000,-1,1,"transverse #tau_{#tau} (ps)","Events");



  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}



void  TauLifeTime::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&TauFlightLength);
 Extradist1d.push_back(&TauFlightLengthTransverse);
 Extradist1d.push_back(&TauMomentum);
 Extradist1d.push_back(&TauMomentumTransverse);
 Extradist1d.push_back(&TauLife);
 Extradist1d.push_back(&TauLifeTransverse);
 Extradist1d.push_back(&ResTauFlightLength);
 Extradist1d.push_back(&ResTauFlightLengthTransverse);
 Extradist1d.push_back(&ResTauMomentum);
 Extradist1d.push_back(&ResTauMomentumTransverse);
 Extradist1d.push_back(&ResTauLife);
 Extradist1d.push_back(&ResTauLifeTransverse);

 Extradist1d.push_back(&TauMomentumAvg);
 Extradist1d.push_back(&TauMomentumTransverseAvg);
 Extradist1d.push_back(&ResTauMomentumAvg);
 Extradist1d.push_back(&ResTauMomentumTransverseAvg);
 Extradist1d.push_back(&TauLifeAvg);
 Extradist1d.push_back(&TauLifeTransverseAvg);
 Extradist1d.push_back(&ResTauLifeAvg);
 Extradist1d.push_back(&ResTauLifeTransverseAvg);
}

void  TauLifeTime::doEvent(){
  unsigned int t;
  int id(0);
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i) && Ntp->PFTau_TIP_hassecondaryVertex(i) && Ntp->PFTau_TIP_hasA1Momentum(i)){
      std::vector<bool> Tau_FitOk;
      std::vector<TLorentzVector> Tau_sol;
      for(unsigned int j=0;j< MultiProngTauSolver::NAmbiguity;j++){
	LorentzVectorParticle theTau;
	std::vector<LorentzVectorParticle> daughter;
	double LC_chi2(0);
	std::cout << "void  TauLifeTime::doEvent() " << j << std::endl;
	Ntp->ThreeProngTauFit(i,j,theTau,daughter,LC_chi2);
      }
    }
  }
}

