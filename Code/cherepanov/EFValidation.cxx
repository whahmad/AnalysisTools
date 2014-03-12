
#include "EFValidation.h"
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
#include "TMath.h"

EFValidation::EFValidation(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
		//  ,ran()
{
  verbose=false;
}

EFValidation::~EFValidation(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "EFValidation::~EFValidation Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "EFValidation::~EFValidation()" << std::endl;
}

void  EFValidation::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasMuon)	      cut.at(hasMuon)=0;
    if(i==hasTau)	      cut.at(hasTau)=0;

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

    else if(i==hasMuon){
      title.at(i)="$Number of good muons$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }


    else if(i==hasTau){
      title.at(i)="$Number of good taus$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTau_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTau_",htitle,6,-0.5,5.5,hlabel,"Events"));
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

  idAllTaus=HConfig.GetTH1D(Name+"_idAllTaus","idAllTaus",4,0.5,4.5,"idAllTaus","Events");
  idIdentifiedTau=HConfig.GetTH1D(Name+"_idIdentifiedTau","idIdentifiedTau",4,0.5,4.5," idIdentifiedTau","Events");
  idPassedEF=HConfig.GetTH1D(Name+"_idPassedEF","idPassedEF",4,0.5,4.5,"idPassedEF ","Events");


  TruthA1Pt=HConfig.GetTH1D(Name+"_TruthA1Pt","Truth A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  TruthA1Phi=HConfig.GetTH1D(Name+"_TruthA1Phi","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  TruthA1Eta=HConfig.GetTH1D(Name+"_TruthA1Eta","Truth A1 Eta",50,-2,2,"A1 Eta","");
  
  TruthMuPt=HConfig.GetTH1D(Name+"_TruthMuPt","Truth Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  TruthMuPhi=HConfig.GetTH1D(Name+"_TruthMuPhi","Truth Mu Phi",50,-3.14,3.14,"#mu Phi, (rad)","");
  TruthMuEta=HConfig.GetTH1D(Name+"_TruthMuEta","Truth Mu Eta",50,-2,2,"#mu Eta","");
  
  TruthTauMuPt=HConfig.GetTH1D(Name+"_TruthTauMuPt","Truth TauMu Pt",50,10,60,"#tau_{#mu} Pt, (GeV)","");
  TruthTauMuPhi=HConfig.GetTH1D(Name+"_TruthTauMuPhi","Truth TauMu Phi",50,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  TruhTauMuEta=HConfig.GetTH1D(Name+"_TruthTauMuEta","Truth TauMu Eta",50,-2,2,"#tau_{#mu} Eta","");
  
  TruthTauA1Pt=HConfig.GetTH1D(Name+"_TruthTauA1Pt","Truth TauA1 Pt",50,10,60,"#tau_{a1} Pt, (GeV)","");
  TruthTauA1Phi=HConfig.GetTH1D(Name+"_TruthTauA1Phi","Truth TauA1 Phi",50,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  TruthTauA1Eta=HConfig.GetTH1D(Name+"_TruthTauA1Eta","Truth TauA1 Eta",50,-2,2,"#tau_{a1} Eta","");
  
  
  TauA1PtResolKFM=HConfig.GetTH1D(Name+"_TauA1PtResolKFM","Pt resolution of #tau_{a1} (TPF ambiguity minus)",50,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1PhiResolKFM=HConfig.GetTH1D(Name+"_TauA1PhiResolKFM","#phi resolution of #tau_{a1} (TPF ambiguity minus)",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolKFM=HConfig.GetTH1D(Name+"_TauA1EtaResolKFM","#eta resolution of #tau_{a1} (TPF ambiguity minus)",50,-1,1,"#delta #eta (mc - rec) ","");
  

  TauA1PtResolKFP=HConfig.GetTH1D(Name+"_TauA1PtResolKFP","Pt resolution of #tau_{a1} (TPF ambiguity plus)",50,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1PhiResolKFP=HConfig.GetTH1D(Name+"_TauA1PhiResolKFP","#phi resolution of #tau_{a1} (TPF ambiguity plus)",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolKFP=HConfig.GetTH1D(Name+"_TauA1EtaResolKFP","#eta resolution of #tau_{a1} (TPF ambiguity plus)",50,-1,1,"#delta #eta (mc - rec) ","");


  RecoA1Pt=HConfig.GetTH1D(Name+"_RecoA1Pt","Reco A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  RecoA1Phi=HConfig.GetTH1D(Name+"_RecoA1Phi","Reco A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  RecoA1Eta=HConfig.GetTH1D(Name+"_RecoA1Eta","Reco A1 Eta",50,-2,2,"A1 Eta","");

  RecoMuPt=HConfig.GetTH1D(Name+"_RecoMuPt","Reco Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  RecoMuPhi=HConfig.GetTH1D(Name+"_RecoMuPhi","Reco Mu Phi",50,-3.14,3.14,"#mu Phi, (rad)","");
  RecoMuEta=HConfig.GetTH1D(Name+"_RecoMuEta","Reco Mu Eta",50,-2,2,"#mu Eta","");


  EFitTauA1Pt=HConfig.GetTH1D(Name+"_FitTauA1Pt"," TauA1 Pt",50,10,60,"#tau_{a1} Pt, (GeV)","");
  EFitTauA1Phi=HConfig.GetTH1D(Name+"_FitTauA1Phi"," TauA1 Phi",50,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  EFitTauA1Eta=HConfig.GetTH1D(Name+"_FitTauA1Eta"," TauA1 Eta",50,-2,2,"#tau_{a1} Eta","");

  EFitTauMuPt=HConfig.GetTH1D(Name+"_FitTauMuPt"," TauMu Pt",50,10,60,"#tau_{#mu} Pt, (GeV)","");
  EFitTauMuPhi=HConfig.GetTH1D(Name+"_FitTauMuPhi"," TauMu Phi",50,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  EFitTauMuEta=HConfig.GetTH1D(Name+"_FitTauMuEta"," TauMu Eta",50,-2,2,"#tau_{#mu} Eta","");



  EFitTauA1PtAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1PtAmbPoint"," TauA1 Pt",50,10,60,"#tau_{a1} Pt, (GeV)","");
  EFitTauA1PhiAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1PhiAmbPoint"," TauA1 Phi",50,-3.14,3.14,"#tau_{a1} Phi, (rad)","");
  EFitTauA1EtaAmbPoint=HConfig.GetTH1D(Name+"_EFitTauA1EtaAmbPoint"," TauA1 Eta",50,-2,2,"#tau_{a1} Eta","");

  EFitTauMuPtAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuPtAmbPoint"," TauMu Pt",50,10,60,"#tau_{#mu} Pt, (GeV)","");
  EFitTauMuPhiAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuPhiAmbPoint"," TauMu Phi",50,-3.14,3.14,"#tau_{#mu} Phi, (rad)","");
  EFitTauMuEtaAmbPoint=HConfig.GetTH1D(Name+"_EFitTauMuEtaAmbPoint"," TauMu Eta",50,-2,2,"#tau_{#mu} Eta","");



  A1Mass=HConfig.GetTH1D(Name+"_A1Mass","a1 mass",50,0.8,1.777,"M_{a1}, GeV","");
  
  TruthA1PtAfterFit=HConfig.GetTH1D(Name+"_TruthA1PtAfterFit","Truth A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  TruthA1EAfterFit=HConfig.GetTH1D(Name+"_TruthA1EAfterFit","Truth A1 E",50,10,100,"A1 E, (GeV)","");
  TruthA1PhiAfterFit=HConfig.GetTH1D(Name+"_TruthA1PhiAfterFit","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  TruthA1EtaAfterFit=HConfig.GetTH1D(Name+"_TruthA1EtaAfterFit","Truth A1 Eta",50,-2,2,"A1 Eta","");
  
  TruthMuPtAfterFit=HConfig.GetTH1D(Name+"_TruthMuPtAfterFit","Truth Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  TruthMuEAfterFit=HConfig.GetTH1D(Name+"_TruthMuEAfterFit","Truth Mu E",50,10,100,"#mu E, (GeV)","");
  TruthMuPhiAfterFit=HConfig.GetTH1D(Name+"_TruthMuPhiAfterFit","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  TruthMuEtaAfterFit=HConfig.GetTH1D(Name+"_TruthMuEtaAfterFit","Truth A1 Eta",50,-2,2,"A1 Eta","");
  
  EfficiencyOverA1Pt=HConfig.GetTH1D(Name+"_EfficiencyOverA1Pt","Truth A1 Pt",50,10,60,"A1 Pt, (GeV)","");
  EfficiencyOverA1Phi=HConfig.GetTH1D(Name+"_EfficiencyOverA1Phi","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  EfficiencyOverA1Eta=HConfig.GetTH1D(Name+"_EfficiencyOverA1Eta","Truth A1 Eta",50,-2,2,"A1 Eta","");


  EfficiencyOverMuPt=HConfig.GetTH1D(Name+"_EfficiencyOverMuPt","Truth Mu Pt",50,10,60,"#mu Pt, (GeV)","");
  EfficiencyOverMuPhi=HConfig.GetTH1D(Name+"_EfficiencyOverMuPhi","Truth A1 Phi",50,-3.14,3.14,"A1 Phi, (rad)","");
  EfficiencyOverMuEta=HConfig.GetTH1D(Name+"_EfficiencyOverMuEta","Truth A1 Eta",50,-2,2,"A1 Eta","");


  
  TauA1PtResolution=HConfig.GetTH1D(Name+"_TauA1PtResolution","Pt resolution of #tau_{a1}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1EResolution=HConfig.GetTH1D(Name+"_TauA1EResolution","E resolution of #tau_{a1}",80,-100,100,"#delta E (mc - rec), GeV ","");
  TauA1PhiResolution=HConfig.GetTH1D(Name+"_TauA1PhiResolution","#phi resolution of #tau_{a1}",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolution=HConfig.GetTH1D(Name+"_TauA1EtaResolution","#eta resolution of #tau_{a1}",50,-1,1,"#delta #eta (mc - rec) ","");
  
  TauMuPtResolution=HConfig.GetTH1D(Name+"_TauMuPtResolution","Pt resolution of #tau_{#mu}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauMuEResolution=HConfig.GetTH1D(Name+"_TauMuEResolution","E resolution of #tau_{#mu}",80,-100,100,"#delta E (mc - rec), GeV ","");
  TauMuPhiResolution=HConfig.GetTH1D(Name+"_TauMuPhiResolution","#phi resolution of #tau_{#mu}",50,-2,2,"#delta #phi (mc - rec) ","");
  TauMuEtaResolution=HConfig.GetTH1D(Name+"_TauMuEtaResolution","#eta resolution of #tau_{#mu}",50,-2,2,"#delta #eta (mc - rec) ","");

  TauA1PtResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauA1PtResolutionAmbPoint","Pt resolution of #tau_{a1}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauA1PhiResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauA1PhiResolutionAmbPoint","#phi resolution of #tau_{a1}",50,-1,1,"#delta #phi (mc - rec) ","");
  TauA1EtaResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauA1EtaResolutionAmbPoint","#eta resolution of #tau_{a1}",50,-1,1,"#delta #eta (mc - rec) ","");

  TauMuPtResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauMuPtResolutionAmbPoint","Pt resolution of #tau_{#mu}",80,-100,100,"#delta Pt (mc - rec), GeV ","");
  TauMuPhiResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauMuPhiResolutionAmbPoint","#phi resolution of #tau_{#mu}",50,-2,2,"#delta #phi (mc - rec) ","");
  TauMuEtaResolutionAmbPoint=HConfig.GetTH1D(Name+"_TauMuEtaResolutionAmbPoint","#eta resolution of #tau_{#mu}",50,-2,2,"#delta #eta (mc - rec) ","");


  TauA1PhiResolVsProb=HConfig.GetTH2D(Name+"_TauA1PhiResolVsProb","Phi resolution of #tau_{a1} vs Probability  ",200,0,1,60,-1,1,"Probability ","#delta #phi (mc - rec)");
  TauMuPhiResolVsProb=HConfig.GetTH2D(Name+"_TauMuPhiResolVsProb","Phi resolution of #tau_{#mu} vs Probability ",200,0,1,60,-1,1,"Probability ","#delta #phi (mc - rec)");

  TauA1EtaResolVsProb=HConfig.GetTH2D(Name+"_TauA1EtaResolVsProb","Eta resolution of #tau_{a1} vs Probability  ",200,0,1,60,-1,1,"Probability ","#delta #eta (mc - rec)");
  TauMuEtaResolVsProb=HConfig.GetTH2D(Name+"_TauMuEtaResolVsProb","Eta resolution of #tau_{#mu} vs Probability  ",200,0,1,60,-1,1,"Probability ","#delta #eta (mc - rec)");

  TauA1PtResolVsProb=HConfig.GetTH2D(Name+"_TauA1PtResolVsProb","Pt resolution of #tau_{a1} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta pt (mc - rec)");
  TauMuPtResolVsProb=HConfig.GetTH2D(Name+"_TauMuPtResolVsProb","Pt resolution of #tau_{#mu} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta pt (mc - rec)");

  TauA1EResolVsProb=HConfig.GetTH2D(Name+"_TauA1EResolVsProb","E resolution of #tau_{a1} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta E (mc - rec)");
  TauMuEResolVsProb=HConfig.GetTH2D(Name+"_TauMuEResolVsProb","E resolution of #tau_{#mu} vs Probability  ",200,0,1,80,-100,100,"Probability ","#delta E (mc - rec)");


  TauA1PtResolVsPt=HConfig.GetTH2D(Name+"_TauA1PtResolVsPt","Pt resolution of #tau_{a1} vs a1 Pt",50,10,60,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsPt=HConfig.GetTH2D(Name+"_TauMuPtResolVsPt","Pt resolution of #tau_{#mu} vs #mu Pt",50,10,60,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");

  TauA1PtResolVsEta=HConfig.GetTH2D(Name+"_TauA1PtResolVsEta","Pt resolution of #tau_{a1} vs a1 Eta",50,-2,2,80,-100,100,"a1 Eta ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsEta=HConfig.GetTH2D(Name+"_TauMuPtResolVsEta","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,80,-100,100,"#mu Eta ","#delta Pt (mc - rec), GeV");
 
  TauA1EtaResolVsEta=HConfig.GetTH2D(Name+"_TauA1EtaResolVsEta","#eta resolution of #tau_{a1} vs a1 Eta",50,-2,2,50,-1,1,"a1 Eta ","#delta #eta (mc - rec)");
  TauMuEtaResolVsEta=HConfig.GetTH2D(Name+"_TauMuEtaResolVsEta","#eta resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,50,-1,1,"#mu Eta ","#delta #eta (mc - rec)");

  
  TauA1PtResolVsZPt=HConfig.GetTH2D(Name+"_TauA1PtResolVsZPt","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsZPt=HConfig.GetTH2D(Name+"_TauMuPtResolVsZPt","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
 
  TauA1EtaResolVsZPt=HConfig.GetTH2D(Name+"_TauA1EtaResolVsZPt","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"a1 Eta ","#delta #eta (mc - rec)");
  TauMuEtaResolVsZPt=HConfig.GetTH2D(Name+"_TauMuEtaResolVsZPt","Pt resolution of #tau_{#mu} vs  Z Pt",50,0,20,80,-100,100,"#mu Eta ","#delta #eta (mc - rec)");

  TauMuPtResolVsZPtMean=HConfig.GetTH1D(Name+"_TauMuPtResolVsZPtMean","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsZPtRMS=HConfig.GetTH1D(Name+"_TauMuPtResolVsZPtRMS","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauA1PtResolVsEtaRMS=HConfig.GetTH1D(Name+"_TauA1PtResolVsEtaRMS","Pt resolution of #tau_{a1} vs a1 Eta",50,-2,2,"a1 Eta ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsPtRMS=HConfig.GetTH1D(Name+"_TauMuPtResolVsPtRMS","Pt resolution of #tau_{#mu} vs #mu Pt",50,10,60,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");


  TauA1PtResolVsPtRMS=HConfig.GetTH1D(Name+"_TauA1PtResolVsPtRMS","Pt resolution of #tau_{a1} vs a1 Pt",50,10,60,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");


  TauA1PhiResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1PhiResolVsPVSVSignificance","Phi resolution of #tau_{a1} vs PVSV significance ",50,0,20,60,-1,1,"#sigma (PV-SV) ","#delta #phi (mc - rec)");
  TauA1EtaResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1EtaResolVsPVSVSignificance","Eta resolution of #tau_{a1} vs PVSV significance ",50,0,20,60,-1,1,"#sigma (PV-SV) ","#delta #eta (mc - rec)");
  TauA1PtResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1PtResolVsPVSVSignificance","Pt resolution of #tau_{a1} vs PVSV significance ",50,0,20,80,-100,100,"#sigma (PV-SV) ","#delta #pt (mc - rec), GeV");
  TauA1EResolVsPVSVSignificance=HConfig.GetTH2D(Name+"_TauA1EResolVsPVSVSignificance","E resolution of #tau_{a1} vs PVSV significance ",50,0,20,80,-100,100,"#sigma (PV-SV) ","#delta E (mc - rec), GeV");
  
  PVSVSignificanceCutEfficiency=HConfig.GetTH1D(Name+"_PVSVSignificanceCutEfficiency","PVSVSignificanceCutEfficiency",50,0,20,"#sigma (PV-SV)","#epsilon of a cut");
  TauA1PtResolVsPVSVSignificanceRMS=HConfig.GetTH1D(Name+"_TauA1PtResolVsPVSVSignificanceRMS","Pt resolution of #tau_{a1} vs PVSV significance ",50,0,20,"#sigma (PV-SV) ","#delta #pt (mc - rec), GeV");


  TauA1PtResolVsPtAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsPtAmbPoint","Pt resolution of #tau_{a1} vs a1 Pt",50,10,60,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsPtAmbPoint=HConfig.GetTH2D(Name+"_TauMuPtResolVsPtAmbPoint","Pt resolution of #tau_{#mu} vs #mu Pt",50,10,60,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauA1PtResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsEtaAmbPoint","Pt resolution of #tau_{a1} vs a1 Eta",50,-2,2,80,-100,100,"a1 Eta ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauMuPtResolVsEtaAmbPoint","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,80,-100,100,"#mu Eta ","#delta Pt (mc - rec), GeV");
  TauA1EtaResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauA1EtaResolVsEtaAmbPoint","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,50,-2,2,"#mu Eta ","#delta Pt (mc - rec), GeV");
  TauMuEtaResolVsEtaAmbPoint=HConfig.GetTH2D(Name+"_TauMuEtaResolVsEtaAmbPoint","Pt resolution of #tau_{#mu} vs  #mu Eta",50,-2,2,50,-2,2,"#mu Eta ","#delta Pt (mc - rec), GeV");

  TauA1PtResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsZPtAmbPoint","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"#tau_{a1} Pt, GeV ","#delta Pt (mc - rec), GeV");
  TauMuPtResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauMuPtResolVsZPtAmbPoint","Pt resolution of #tau_{#mu} vs Z  Pt",50,0,20,80,-100,100,"#tau_{#mu} Pt, GeV ","#delta Pt (mc - rec), GeV");

  TauA1EtaResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauA1EtaResolVsZPtAmbPoint","Pt resolution of #tau_{a1} vs Z Pt",50,0,20,80,-100,100,"a1 Eta ","#delta #eta (mc - rec)");
  TauMuEtaResolVsZPtAmbPoint=HConfig.GetTH2D(Name+"_TauMuEtaResolVsZPtAmbPoint","Pt resolution of #tau_{#mu} vs  Z Pt",50,0,20,80,-100,100,"#mu Eta ","#delta #eta (mc - rec)");


  TauA1PhiResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1PhiResolVsPhiRotSignificance","Phi resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (phi rotation) ","#delta #phi (mc - rec)");
  TauA1EtaResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1EtaResolVsPhiRotSignificance","Eta resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (phi rotation) ","#delta #eta (mc - rec)");
  TauA1PtResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1PtResolVsPhiRotSignificance","Pt resolution of #tau_{a1} vs PVSV significance ",50,0,10,80,-100,100,"#sigma (phi rotation) ","#delta #pt (mc - rec), GeV");
  TauA1EResolVsPhiRotSignificance=HConfig.GetTH2D(Name+"_TauA1EResolVsPhiRotSignificance","E resolution of #tau_{a1} vs PVSV significance ",50,0,10,80,-100,100,"#sigma (phi rotation) ","#delta E (mc - rec), GeV");
  
  PhiRotSignificanceCutEfficiency=HConfig.GetTH1D(Name+"_PhiRotSignificanceCutEfficiency","PhiRotSignificanceCutEfficiency",50,0,10,"#sigma (phi rotation) ","#epsilon of a cut");


  TauA1PhiResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1PhiResolVsPVSVSignificanceAmbPoint","Phi resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (PV-SV) ","#delta #phi (mc - rec)");
  TauA1EtaResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1EtaResolVsPVSVSignificanceAmbPoint","Eta resolution of #tau_{a1} vs PVSV significance ",50,0,10,60,-1,1,"#sigma (PV-SV) ","#delta #phi (mc - rec)");
  TauA1PtResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1PtResolVsPVSVSignificanceAmbPoint","Pt resolution of #tau_{a1} vs PVSV significanc ",50,0,10,80,-100,100,"#sigma (PV-SV)","#delta #pt (mc - rec), GeV");
  TauA1EResolVsPVSVSignificanceAmbPoint=HConfig.GetTH2D(Name+"_TauA1EResolVsPVSVSignificanceAmbPoint","E resolution of #tau_{a1} vs PVSV significance ",50,0,10,80,-100,100,"#sigma (PV-SV) ","#delta E (mc - rec), GeV");


  
  ProbabilityOfCorrect=HConfig.GetTH1D(Name+"_ProbabilityOfCorrect","ProbabilityOfCorrect",300,0,1,"Probability","Events");

  ProbabilityOfCorrectndf2=HConfig.GetTH1D(Name+"_ProbabilityOfCorrectndf2","ProbabilityOfCorrectndf2",300,0,1,"Probability ndf 2","Events");
  ProbabilityOfCorrectndf3=HConfig.GetTH1D(Name+"_ProbabilityOfCorrectndf3","ProbabilityOfCorrectndf3",300,0,1,"Probability ndf 3","Events");
  Chi2OfCorrect=HConfig.GetTH1D(Name+"_Chi2OfCorrect","Chi2OfCorrect",300,0,50,"Probability","Events");


  ProbabilityOfAmbiguityPoint=HConfig.GetTH1D(Name+"_ProbabilityOfAmbiguityPoint","ProbabilityOfAmbiguityPoint",300,0,1,"Probability","Events");
  Chi2Dim=HConfig.GetTH2D(Name+"_Chi2Dim","Chi2Dim",100,0,1,100,0,1,"Chi2Dim","Events");
  csum2Dim=HConfig.GetTH2D(Name+"_csum2Dim","csum2Dim",200,0,1,200,0,1,"csum, GeV","Events");
  CorrectAmbiguityTruth=HConfig.GetTH1D(Name+"_CorrectAmbiguityTruth","CorrectAmbiguityTruth",2,0.5,2.5,"Correct Ambiguity Truth ","Events");

  MeasuredPhi=HConfig.GetTH1D(Name+"_MeasuredPhi","Measured TauA1 phi",40,-3.14,3.14,"  #phi, (rad) ","");

  pvsvsignificance=HConfig.GetTH1D(Name+"_pvsvsignificance","pvsv significance",40,0,10,"  #sigma (PV-SV) ","");
  PhiRotationSignificanceUnPhysicalTaus=HConfig.GetTH1D(Name+"_PhiRotationSignificanceUnPhysicalTaus","Significance of #phi angle rotation",20,0,10,"#sigma(#phi)","");


  csumProb2Dim=HConfig.GetTH2D(Name+"_csumProb2Dim","csumProb2Dim",200,0,0.5,200,0,0.2,"csum, GeV","probability");
  ProbabilityOfCorrectVsZPt=HConfig.GetTH2D(Name+"_ProbabilityOfCorrectVsZPt","ProbabilityOfCorrectVsZPt",400,0,0.5,100,0,50,"probability","Z pt, GeV");
  ProbabilityOfCorrectVsChi2=HConfig.GetTH2D(Name+"_ProbabilityOfCorrectVsChi2","ProbabilityOfCorrectVsChi2",400,0,1,100,0,50,"probability","#chi^{2}");


  ChannelEfficiency=HConfig.GetTH1D(Name+"_ChannelEfficiency","ChannelEfficiency",4,0.5,4.5,"  0 - signal, 1 - 4#pi, 2 - K#pi#pi, 3 - KK#pi ","");


  TauMuPxPull=HConfig.GetTH1D(Name+"_TauMuPxPull","TauMuPxPull",50,-100,100,"TauMu Px, pull ","");
  TauMuPyPull=HConfig.GetTH1D(Name+"_TauMuPyPull","TauMuPyPull",50,-100,100,"TauMu Py, pull ","");
  TauMuPzPull=HConfig.GetTH1D(Name+"_TauMuPzPull","TauMuPzPull",50,-100,100,"TauMu Pz, pull ","");

  TauA1PxPull=HConfig.GetTH1D(Name+"_TauA1PxPull","TauA1PxPull",50,-100,100,"TauA1 Px, pull ","");
  TauA1PyPull=HConfig.GetTH1D(Name+"_TauA1PyPull","TauA1PyPull",50,-100,100,"TauA1 Py, pull ","");
  TauA1PzPull=HConfig.GetTH1D(Name+"_TauA1PzPull","TauA1PzPull",50,-100,100,"TauA1 Pz, pull ","");

  TauA1SFPxPull=HConfig.GetTH1D(Name+"_TauA1SFPxPull","TauA1SFPxPull",50,-100,100,"TauA1 Px, pull ","");
  TauA1SFPyPull=HConfig.GetTH1D(Name+"_TauA1SFPyPull","TauA1SFPyPull",50,-100,100,"TauA1 Py, pull ","");
  TauA1SFPzPull=HConfig.GetTH1D(Name+"_TauA1SFPzPull","TauA1SFPzPull",50,-100,100,"TauA1 Pz, pull ","");

  SVxPull=HConfig.GetTH1D(Name+"_SVxPull","SVxPull",50,-5,5,"SVx, pull ","");
  SVyPull=HConfig.GetTH1D(Name+"_SVyPull","SVyPull",50,-5,5,"SVy, pull ","");
  SVzPull=HConfig.GetTH1D(Name+"_SVzPull","SVzPull",50,-5,5,"SVz, pull ","");



  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
  //  ran = new TRandom();
}



void  EFValidation::Store_ExtraDist(){


 Extradist1d.push_back(&TruthA1Pt);
 Extradist1d.push_back(&TruthA1Phi);
 Extradist1d.push_back(&TruthA1Eta);
 
 Extradist1d.push_back(&TruthMuPt);
 Extradist1d.push_back(&TruthMuPhi);
 Extradist1d.push_back(&TruthMuEta);

 Extradist1d.push_back(&TruthTauMuPt);
 Extradist1d.push_back(&TruthTauMuPhi);
 Extradist1d.push_back(&TruhTauMuEta);
 
 Extradist1d.push_back(&TruthTauA1Pt);
 Extradist1d.push_back(&TruthTauA1Phi);
 Extradist1d.push_back(&TruthTauA1Eta);
 
 
 Extradist1d.push_back(&TauA1PtResolKFM);
 Extradist1d.push_back(&TauA1PhiResolKFM);
 Extradist1d.push_back(&TauA1EtaResolKFM);
 
 
 Extradist1d.push_back(&TauA1PtResolKFP);
 Extradist1d.push_back(&TauA1PhiResolKFP);
 Extradist1d.push_back(&TauA1EtaResolKFP);
 
 
 Extradist1d.push_back(&RecoA1Pt);
 Extradist1d.push_back(&RecoA1Phi);
 Extradist1d.push_back(&RecoA1Eta);
 
 Extradist1d.push_back(&RecoMuPt);
 Extradist1d.push_back(&RecoMuPhi);
 Extradist1d.push_back(&RecoMuEta);
 
 Extradist1d.push_back(&A1Mass);
 
 Extradist1d.push_back(&TruthA1PtAfterFit);
 Extradist1d.push_back(&TruthA1EAfterFit);
 Extradist1d.push_back(&TruthA1PhiAfterFit);
 Extradist1d.push_back(&TruthA1EtaAfterFit);
 
 Extradist1d.push_back(&TruthMuPtAfterFit);
 Extradist1d.push_back(&TruthMuEAfterFit);
 Extradist1d.push_back(&TruthMuPhiAfterFit);
 Extradist1d.push_back(&TruthMuEtaAfterFit);
 
 
 Extradist1d.push_back(&EFitTauA1Pt);
 Extradist1d.push_back(&EFitTauA1Phi);
 Extradist1d.push_back(&EFitTauA1Eta);

 Extradist1d.push_back(&EFitTauMuPt);
 Extradist1d.push_back(&EFitTauMuPhi);
 Extradist1d.push_back(&EFitTauMuEta);



 Extradist1d.push_back(&EFitTauA1PtAmbPoint);
 Extradist1d.push_back(&EFitTauA1PhiAmbPoint);
 Extradist1d.push_back(&EFitTauA1EtaAmbPoint);

 Extradist1d.push_back(&EFitTauMuPtAmbPoint);
 Extradist1d.push_back(&EFitTauMuPhiAmbPoint);
 Extradist1d.push_back(&EFitTauMuEtaAmbPoint);

 
 Extradist1d.push_back(&TauA1PtResolution);
 Extradist1d.push_back(&TauA1EResolution);
 Extradist1d.push_back(&TauA1PhiResolution);
 Extradist1d.push_back(&TauA1EtaResolution);
 
 Extradist1d.push_back(&TauMuPtResolution);
 Extradist1d.push_back(&TauMuEResolution);
 Extradist1d.push_back(&TauMuPhiResolution);
 Extradist1d.push_back(&TauMuEtaResolution);
 
 
 Extradist2d.push_back(&TauA1PhiResolVsProb);
 Extradist2d.push_back(&TauMuPhiResolVsProb);
 Extradist2d.push_back(&TauA1EtaResolVsProb);
 Extradist2d.push_back(&TauMuEtaResolVsProb);
 Extradist2d.push_back(&TauA1PtResolVsProb);
 Extradist2d.push_back(&TauMuPtResolVsProb);
 Extradist2d.push_back(&TauA1EResolVsProb);
 Extradist2d.push_back(&TauMuEResolVsProb);
      
 Extradist2d.push_back(&TauA1PtResolVsPt);
 Extradist2d.push_back(&TauMuPtResolVsPt);
 Extradist2d.push_back(&TauA1PtResolVsEta);
 Extradist2d.push_back(&TauMuPtResolVsEta);
 Extradist2d.push_back(&TauA1EtaResolVsEta);
 Extradist2d.push_back(&TauMuEtaResolVsEta);
 
 
 Extradist2d.push_back(&TauA1PtResolVsZPt);
 Extradist2d.push_back(&TauMuPtResolVsZPt);
 
 Extradist2d.push_back(&TauA1EtaResolVsZPt);
 Extradist2d.push_back(&TauMuEtaResolVsZPt);
 
 

 
 Extradist2d.push_back(&TauA1PhiResolVsPVSVSignificance);
 Extradist2d.push_back(&TauA1EtaResolVsPVSVSignificance);
 Extradist2d.push_back(&TauA1PtResolVsPVSVSignificance);
 Extradist2d.push_back(&TauA1EResolVsPVSVSignificance);
 
 Extradist2d.push_back(&TauA1PtResolVsPtAmbPoint);
 Extradist2d.push_back(&TauMuPtResolVsPtAmbPoint);
 Extradist2d.push_back(&TauA1PtResolVsEtaAmbPoint);
 Extradist2d.push_back(&TauMuPtResolVsEtaAmbPoint);
 Extradist2d.push_back(&TauA1EtaResolVsEtaAmbPoint);
 Extradist2d.push_back(&TauMuEtaResolVsEtaAmbPoint);
 
 Extradist2d.push_back(&TauA1PtResolVsZPtAmbPoint);
 Extradist2d.push_back(&TauMuPtResolVsZPtAmbPoint);
 
 Extradist2d.push_back(&TauA1EtaResolVsZPtAmbPoint);
 Extradist2d.push_back(&TauMuEtaResolVsZPtAmbPoint);
 
 Extradist1d.push_back(&TauA1PtResolutionAmbPoint);
 Extradist1d.push_back(&TauA1PhiResolutionAmbPoint);
 Extradist1d.push_back(&TauA1EtaResolutionAmbPoint);
 
 Extradist1d.push_back(&TauMuPtResolutionAmbPoint);
 Extradist1d.push_back(&TauMuPhiResolutionAmbPoint);
 Extradist1d.push_back(&TauMuEtaResolutionAmbPoint);
 
 Extradist2d.push_back(&TauA1PhiResolVsPhiRotSignificance);
 Extradist2d.push_back(&TauA1EtaResolVsPhiRotSignificance);
 Extradist2d.push_back(&TauA1PtResolVsPhiRotSignificance);
 Extradist2d.push_back(&TauA1EResolVsPhiRotSignificance);
 
 
 Extradist2d.push_back(&TauA1PhiResolVsPVSVSignificanceAmbPoint);
 Extradist2d.push_back(&TauA1EtaResolVsPVSVSignificanceAmbPoint);
 Extradist2d.push_back(&TauA1PtResolVsPVSVSignificanceAmbPoint);
 Extradist2d.push_back(&TauA1EResolVsPVSVSignificanceAmbPoint);
 
 Extradist2d.push_back(&Chi2Dim);
 Extradist2d.push_back(&csum2Dim);
 Extradist1d.push_back(&ProbabilityOfCorrect);


 Extradist1d.push_back(&idAllTaus);
 Extradist1d.push_back(&idIdentifiedTau);
 Extradist1d.push_back(&idPassedEF);
 Extradist1d.push_back(&CorrectAmbiguityTruth);
 Extradist1d.push_back(&MeasuredPhi);
 Extradist1d.push_back(&pvsvsignificance);

 Extradist1d.push_back(&PhiRotationSignificanceUnPhysicalTaus);

 Extradist2d.push_back(&csumProb2Dim);
 Extradist2d.push_back(&ProbabilityOfCorrectVsZPt);
 Extradist1d.push_back(&ProbabilityOfAmbiguityPoint);



 Extradist1d.push_back(&TauA1PtResolVsPVSVSignificanceRMS);
 Extradist1d.push_back(&TauA1PtResolVsPtRMS);
 Extradist1d.push_back(&TauMuPtResolVsZPtRMS);
 Extradist1d.push_back(&TauA1PtResolVsEtaRMS);
 Extradist1d.push_back(&TauMuPtResolVsPtRMS);
 Extradist1d.push_back(&TauMuPtResolVsZPtMean);



 Extradist1d.push_back(&EfficiencyOverA1Pt);
 Extradist1d.push_back(&EfficiencyOverA1Phi);
 Extradist1d.push_back(&EfficiencyOverA1Eta);


 Extradist1d.push_back(&EfficiencyOverMuPt);
 Extradist1d.push_back(&EfficiencyOverMuPhi);
 Extradist1d.push_back(&EfficiencyOverMuEta);

 Extradist1d.push_back(&ProbabilityOfCorrectndf2);
 Extradist1d.push_back(&ProbabilityOfCorrectndf3);
 Extradist1d.push_back(&Chi2OfCorrect);


 Extradist2d.push_back(&ProbabilityOfCorrectVsChi2);
 Extradist1d.push_back(&PhiRotSignificanceCutEfficiency);
 Extradist1d.push_back(&PVSVSignificanceCutEfficiency);


 Extradist1d.push_back(&ChannelEfficiency);


 Extradist1d.push_back(&TauMuPxPull);
 Extradist1d.push_back(&TauMuPyPull);
 Extradist1d.push_back(&TauMuPzPull);

 Extradist1d.push_back(&TauA1PxPull);
 Extradist1d.push_back(&TauA1PyPull);
 Extradist1d.push_back(&TauA1PzPull);


 Extradist1d.push_back(&TauA1SFPxPull);
 Extradist1d.push_back(&TauA1SFPyPull);
 Extradist1d.push_back(&TauA1SFPzPull);

 Extradist1d.push_back(&SVxPull);
 Extradist1d.push_back(&SVyPull);
 Extradist1d.push_back(&SVzPull);

}

void  EFValidation::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
    std::cout << "EFValidation Ntp->GetMCID(): " << Ntp->GetMCID() << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){/* std::cout << "failed to find id" <<std::endl;*/ return;}

  if(verbose)std::cout << "void  EFValidation::doEvent() Mu" << std::endl;


//   // check out triggers
//   for(unsigned int itr = 0; itr <  Ntp->NHLTTriggers(); itr ++){
//     std::cout << "Trigger names  " <<Ntp->HTLTriggerName(itr) << std::endl;
//   }

  int idToFill(0);
  if(id ==998)idToFill = 1;
  if(id ==10230833)idToFill = 2;
  if(id ==10231433)idToFill = 3;
  if(id ==10231833)idToFill = 4;
  idAllTaus.at(t).Fill(idToFill,1);


  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("eta2p1_LooseIsoPFTau"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);

//   for(unsigned int iTrig =0; iTrig <Ntp-> NHLTTriggers(); iTrig ++){

//     std::cout<<" Triggers  "<< Ntp->HTLTriggerName(iTrig)<<std::endl;
//     }
  // Apply Selection
  std::vector<unsigned int> mu_idx_good, mu_idx_pt, mu_idx_iso;
  unsigned mu_idx(999);
  //double mu_pt(0);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
    //std::cout << "Checking if muon number = " << i << " is good..." << std::endl;
    if(Ntp->isGoodMuon(i) && fabs(Ntp->Muons_p4(i).Eta())<2.4 && Ntp->Muons_p4(i).Pt() > 17){
      mu_idx_good.push_back(i);
    }  
  }



  value.at(hasMuon)=mu_idx_good.size();
  pass.at(hasMuon)=(value.at(hasMuon)>cut.at(hasMuon));
  
  std::vector<unsigned int> tau_idx_pt, tau_idx_eta, tau_idx_iso, tau_idx_good ;
  unsigned int tau_idx(999);
  for(unsigned int i=0;i<Ntp->NPFTaus();i++){
    if(Ntp->PFTau_p4(i).Pt()> 17){
      if(fabs(Ntp->PFTau_p4(i).Eta())<2.7){
	if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i)  ){ 
	//	if(Ntp->PFTau_hpsDecayMode(i)==10 && Ntp->PFTau_isHPSByDecayModeFinding(i)  ){ 
	  tau_idx_good.push_back(i);
	}
      }
    }
  }



  value.at(hasTau)=tau_idx_good.size();
  pass.at(hasTau)=(value.at(hasTau)>cut.at(hasTau));

  std::cout<<"ntaus  "<<tau_idx_good.size()<<std::endl;
  std::cout<<"nmuons  "<<mu_idx_good.size()<<std::endl;

   tau_idx =   tau_idx_good.at(0);
   mu_idx=   mu_idx_good.at(0);



  double wobs(1),w(1),w1(1);
   if(!Ntp->isData()){
     w1*=Ntp->EvtWeight3D();
   }
   else{w1=1;}

  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;
    //std::cout << "tau_idx: " << tau_idx << "; mu_idx: " << mu_idx << std::endl;
    //std::cout << "tau_idx_iso.size: " << tau_idx_iso.size() << " mu_idx_iso.size: " << mu_idx_iso.size() << std::endl;
    int isInPhysicalRegion(0);
    double significan(0);
    if (tau_idx!=999 && mu_idx!=999){
      //>>>>>>>>> search for a matched jet 
      unsigned int matchedJet(0);
      for(unsigned int ijet =0; ijet < Ntp->NPFJets(); ijet++){
	
	double delta = 999.;
	if(Ntp->PFTau_p4(tau_idx).DeltaR(Ntp->PFJet_p4(ijet)) < delta ){
	  delta = Ntp->PFTau_p4(tau_idx).DeltaR(Ntp->PFJet_p4(ijet));
	  matchedJet=ijet;
	}
      }
      //<<<<<<<<< search for a matched jet


      //>>>>>>>>>>>>>>>>>>>>>>>Fill histos for efficiency 
  
      if(Ntp->CheckDecayID(2,5)){
	std::cout<<" 1 "<<std::endl;
	TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	
	TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
    
	TruthA1Pt.at(t).Fill(TruthA1.Pt(),1);
	TruthA1Phi.at(t).Fill(TruthA1.Phi(),1);
	TruthA1Eta.at(t).Fill(TruthA1.Eta(),1);

	TruthMuPt.at(t).Fill(TruthMu.Pt(),1);
	TruthMuPhi.at(t).Fill(TruthMu.Phi(),1);
	TruthMuEta.at(t).Fill(TruthMu.Eta(),1);

	TruthTauMuPt.at(t).Fill(TruthMu.Pt(),1);
	TruthTauMuPhi.at(t).Fill(TruthMu.Phi(),1);
	TruhTauMuEta.at(t).Fill(TruthMu.Eta(),1);
	
	TruthTauA1Pt.at(t).Fill(TruthA1.Pt(),1);
	TruthTauA1Phi.at(t).Fill(TruthA1.Phi(),1);
	TruthTauA1Eta.at(t).Fill(TruthA1.Eta(),1);

	double angle = TruthTauA1.Angle(TruthA1.Vect());
	double rootAmb=sqrt((TruthA1.M()*TruthA1.M()+TruthA1.P()*TruthA1.P())*(pow(TruthA1.M()*TruthA1.M()-TruthTauA1.M()*TruthTauA1.M(),2)-4*TruthTauA1.M()*TruthTauA1.M()*TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle)));
	double PtruthMinus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  -rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	double PtruthPlus  = ((TruthA1.M()*TruthA1.M()+TruthTauA1.M()*TruthTauA1.M() )*TruthA1.P()*cos(angle)  +rootAmb )/(TruthA1.M()*TruthA1.M() + TruthA1.P()*TruthA1.P()*sin(angle)*sin(angle))/2;
	
	int TruthAmibiga =0;
	
	if( fabs(TruthTauA1.P()  -PtruthMinus ) < fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =1;
	if( fabs(TruthTauA1.P()  -PtruthMinus ) > fabs(TruthTauA1.P()  - PtruthPlus))TruthAmibiga =2;
	CorrectAmbiguityTruth.at(t).Fill(TruthAmibiga,1.);
      }
  

      //<<<<<<<<<<<<<<<<<<<<<<<Fill histos for efficiency 
      
      int idToFill(0);
      if(id ==998)idToFill = 1;
      if(id ==10230833)idToFill = 2;
      if(id ==10231433)idToFill = 3;
      if(id ==10231833)idToFill = 4;
      if(Ntp->PFTau_hpsDecayMode(tau_idx) == 10)idIdentifiedTau.at(t).Fill(idToFill,1);

      std::vector<bool> Tau_FitOk; 
      std::vector<double> Tau_FitChi2;
      std::vector<TLorentzVector> TauA1_TPF;


      std::vector<bool> EventFit_Ok;
      std::vector<double> SignificanceOfTauRotation;

      std::vector<double> DeltaPhiSign;

      std::vector<double> deltaphiRotation;

      std::vector<TLorentzVector> TauA1_EF;
      std::vector<TLorentzVector> TauMu_EF;
      std::vector<TLorentzVector> TauMuNetrino_EF;
      std::vector<TLorentzVector> TauA1Neutrino_TPF;

      std::vector<TLorentzVector> Z_EF;

      std::vector<double> Chi2EventFit;
      std::vector<double> Chi2ProbabilityEventFit;
      std::vector<double> Chi2ProbabilityEventFitndf2;
      std::vector<double> Chi2ProbabilityEventFitndf3;


      std::vector<double> NiterationsEventFit;
      std::vector<double> csum_GF;

      double flightsignificanceA1Side;
      double flightsignificanceMuSide;
      std::vector<LorentzVectorParticle> theTauAmbiguityVector;

      std::vector<LorentzVectorParticle> theTauA1TLV;
      std::vector<LorentzVectorParticle> theTauMuTLV;

      TVector3 pv;
      TMatrixTSym<double> PVcov;
      TVector3 sv;
      TMatrixTSym<double> SVcov;
    
      TLorentzVector A1Reco = Ntp->PFTau_a1_lvp(tau_idx).LV();
      TLorentzVector MuReco = Ntp->Muons_p4(mu_idx);

      RecoA1Pt.at(t).Fill(A1Reco.Pt(),1);
      RecoA1Phi.at(t).Fill(A1Reco.Phi(),1);
      RecoA1Eta.at(t).Fill(A1Reco.Eta(),1);

      RecoMuPt.at(t).Fill(MuReco.Pt(),1);
      RecoMuPhi.at(t).Fill(MuReco.Phi(),1);
      RecoMuEta.at(t).Fill(MuReco.Eta(),1);

      A1Mass.at(t).Fill(A1Reco.M(),1);

      //////////////////////////////
      //loop over ambiguity points
      MeasuredPhi.at(t).Fill(Ntp->MeasuredTauDirection(tau_idx).Phi(),1);

      for(unsigned int j=0;j< MultiProngTauSolver::NAmbiguity;j++){

	LorentzVectorParticle theTau;
	LorentzVectorParticle theZ;
	std::vector<LorentzVectorParticle> daughter;
	std::vector<LorentzVectorParticle> theZdaughter;

	LorentzVectorParticle  LVPEF_TauMu;
	LorentzVectorParticle  LVPEF_TauA1;
	       
	TLorentzVector NeutrinoA1(0,0,0,0);
	TLorentzVector NeutrinoMu(0,0,0,0);
	       
	TLorentzVector TauA1EventFit(0,0,0,0);
	TLorentzVector TauA1ThreeProngFit(0,0,0,0);
	TLorentzVector TauMuEventFit(0,0,0,0);
	TLorentzVector Z_sol(0,0,0,0);
	double deltaphi(0);
	       
	double LC_Eventchi2(0);
	double LC_Eventchi2Probability(0);
	double LC_Eventchi2Probabilityndf2(0);
	double LC_Eventchi2Probabilityndf3(0);


	double LC_chi2(0);
	double phisign(0);
	int NiterationsEF(0);
	double csum(0);
	double SEA1(0);
	double SEMU(0);

	bool  ThreeProngFitSuccess =false;
	bool  EventFitSuccess =false;
	ThreeProngFitSuccess=Ntp->ThreeProngTauFit(tau_idx,j,theTau,daughter,LC_chi2,phisign);
	MeasuredPhi.at(t).Fill(Ntp->MeasuredTauDirection(tau_idx).Phi(),1);
	if(ThreeProngFitSuccess){
	  //>>>>>>>>>>> Fill Resol plots dor TPF
	  if(Ntp->CheckDecayID(2,5)){

	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    if(j==1){
	      TauA1PtResolKFM.at(t).Fill(TruthTauA1.Pt() - theTau.LV().Pt(),1);
	      TauA1PhiResolKFM.at(t).Fill(TruthTauA1.Phi() - theTau.LV().Phi(),1);
	      TauA1EtaResolKFM.at(t).Fill(TruthTauA1.Eta() - theTau.LV().Eta(),1);
	    }
	    if(j==2){
	      TauA1PtResolKFP.at(t).Fill(TruthTauA1.Pt() -  theTau.LV().Pt(),1);
	      TauA1PhiResolKFP.at(t).Fill(TruthTauA1.Phi() - theTau.LV().Phi(),1);
	      TauA1EtaResolKFP.at(t).Fill(TruthTauA1.Eta() - theTau.LV().Eta(),1);
	    }
	  }
	  
	  
	  NeutrinoA1=daughter.at(1).LV();
	  TauA1ThreeProngFit=theTau.LV();
	  
	  EventFitSuccess = Ntp->EventFit(tau_idx,mu_idx,theTau,theZ,theZdaughter,LC_Eventchi2,NiterationsEF,csum);
	  if(EventFitSuccess){
	    LC_Eventchi2Probability = TMath::Prob(LC_Eventchi2,1);
	    LC_Eventchi2Probabilityndf2 = TMath::Prob(LC_Eventchi2,2);
	    LC_Eventchi2Probabilityndf3 = TMath::Prob(LC_Eventchi2,3);


	    TauA1EventFit =theZdaughter.at(0).LV();
	    TauMuEventFit =theZdaughter.at(1).LV();

	    LVPEF_TauMu = theZdaughter.at(0);
	    LVPEF_TauA1 = theZdaughter.at(1);
	      

	    NeutrinoMu=TauMuEventFit-Ntp->Muons_p4(mu_idx);
	    Z_sol=theZ.LV();
	  }
	}

	SignificanceOfTauRotation.push_back(phisign);
	TauA1Neutrino_TPF.push_back(NeutrinoA1);
	TauMuNetrino_EF.push_back(NeutrinoMu);

	Tau_FitOk.push_back(ThreeProngFitSuccess);
	Tau_FitChi2.push_back(LC_chi2);
  
	Chi2ProbabilityEventFit.push_back(LC_Eventchi2Probability);
	Chi2ProbabilityEventFitndf2.push_back(LC_Eventchi2Probabilityndf2);
	Chi2ProbabilityEventFitndf3.push_back(LC_Eventchi2Probabilityndf3);

	EventFit_Ok.push_back(EventFitSuccess);
  	Z_EF.push_back(Z_sol);
 	Chi2EventFit.push_back(LC_Eventchi2);
	theTauAmbiguityVector.push_back(theTau);

 	theTauA1TLV.push_back(LVPEF_TauA1);
 	theTauMuTLV.push_back(LVPEF_TauMu);


 	TauA1_EF.push_back(TauA1EventFit);
 	TauMu_EF.push_back(TauMuEventFit);
	NiterationsEventFit.push_back(NiterationsEF);
	csum_GF.push_back(csum);
	DeltaPhiSign.push_back(phisign);   
	deltaphiRotation.push_back(deltaphi);
      }
      
      PhiRotationSignificanceUnPhysicalTaus.at(t).Fill(SignificanceOfTauRotation.at(0));
    


      int AmbiguitySolution =0;
      bool AmbPoint(false);
      if(Ntp->AmbiguitySolver(Tau_FitOk,EventFit_Ok,Chi2ProbabilityEventFit,AmbiguitySolution, AmbPoint)){

 	double flightsignificance=Ntp->PFTau_FlightLenght_significance2(Ntp->PFTau_TIP_primaryVertex_pos(tau_idx),Ntp->PFTau_TIP_primaryVertex_cov(tau_idx),theTauAmbiguityVector.at(AmbiguitySolution).Vertex(),theTauAmbiguityVector.at(AmbiguitySolution).VertexCov());
 	pvsvsignificance.at(t).Fill(flightsignificance,1);


	int idToFill(0);
	if(id ==998)idToFill = 1;
	if(id ==10230833)idToFill = 2;
	if(id ==10231433)idToFill = 3;
	if(id ==10231833)idToFill = 4;
	idPassedEF.at(t).Fill(idToFill,1);


	TLorentzVector EventFitTauA1 = TauA1_EF.at(AmbiguitySolution);
	TLorentzVector EventFitTauMu = TauMu_EF.at(AmbiguitySolution);

	TLorentzVector EventFitZ = EventFitTauA1+EventFitTauMu;

	if(!AmbPoint){
	  EFitTauA1Pt.at(t).Fill(EventFitTauA1.Pt(),1);
	  EFitTauA1Phi.at(t).Fill(EventFitTauA1.Phi(),1);
	  EFitTauA1Eta.at(t).Fill(EventFitTauA1.Eta(),1);

	  EFitTauMuPt.at(t).Fill(EventFitTauMu.Pt(),1);
	  EFitTauMuPhi.at(t).Fill(EventFitTauMu.Phi(),1);
	  EFitTauMuEta.at(t).Fill(EventFitTauMu.Eta(),1);
	}


	if(AmbPoint){
	  EFitTauA1PtAmbPoint.at(t).Fill(EventFitTauA1.Pt(),1);
	  EFitTauA1PhiAmbPoint.at(t).Fill(EventFitTauA1.Phi(),1);
	  EFitTauA1EtaAmbPoint.at(t).Fill(EventFitTauA1.Eta(),1);

	  EFitTauMuPtAmbPoint.at(t).Fill(EventFitTauMu.Pt(),1);
	  EFitTauMuPhiAmbPoint.at(t).Fill(EventFitTauMu.Phi(),1);
	  EFitTauMuEtaAmbPoint.at(t).Fill(EventFitTauMu.Eta(),1);
	}

	  if(Ntp->CheckDecayID(2,5)){
	    TLorentzVector TruthA1 = Ntp->GetTruthTauProductLV(5, 20213);
	    TLorentzVector TruthMu = Ntp->GetTruthTauProductLV(2, 13);
	
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    

	    TVector3 TruthTauA1SV = Ntp->GetTruthTauProductVertex(5,211);
	    std::cout<<" TruthTauA1SV  "<<TruthTauA1SV.X()<<" TruthTauA1SV  "<<TruthTauA1SV.Y()<<" TruthTauA1SV  "<<TruthTauA1SV.Z()<<std::endl;

	    std::cout<<" reco sv  "<<theTauMuTLV.at(AmbiguitySolution).Vertex().X()<<" TruthTauA1SV  "<<theTauMuTLV.at(AmbiguitySolution).Vertex().Y()<<" TruthTauA1SV  "<<theTauMuTLV.at(AmbiguitySolution).Vertex().Z()<<std::endl;
	    SVxPull.at(t).Fill((TruthTauA1SV.X() - theTauMuTLV.at(AmbiguitySolution).Vertex().X() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(0,0)),1);
	    SVyPull.at(t).Fill((TruthTauA1SV.Y() - theTauMuTLV.at(AmbiguitySolution).Vertex().Y() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(1,1)),1);
	    SVzPull.at(t).Fill((TruthTauA1SV.Z() - theTauMuTLV.at(AmbiguitySolution).Vertex().Z() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(2,2)),1);

	    TruthA1PtAfterFit.at(t).Fill(TruthA1.Pt(),1);
	    TruthA1EAfterFit.at(t).Fill(TruthA1.E(),1);
	    TruthA1PhiAfterFit.at(t).Fill(TruthA1.Phi(),1);
	    TruthA1EtaAfterFit.at(t).Fill(TruthA1.Eta(),1);

	    TruthMuPtAfterFit.at(t).Fill(TruthMu.Pt(),1);
	    TruthMuEAfterFit.at(t).Fill(TruthMu.E(),1);
	    TruthMuPhiAfterFit.at(t).Fill(TruthMu.Phi(),1);
	    TruthMuEtaAfterFit.at(t).Fill(TruthMu.Eta(),1);



	    if(!AmbPoint){
	      TauA1PtResolution.at(t).Fill(TruthTauA1.Pt()-EventFitTauA1.Pt(),1);
	      TauA1EResolution.at(t).Fill(TruthTauA1.E()-EventFitTauA1.E(),1);
	      TauA1PhiResolution.at(t).Fill(TruthTauA1.Phi()-EventFitTauA1.Phi(),1);
	      TauA1EtaResolution.at(t).Fill(TruthTauA1.Eta()-EventFitTauA1.Eta(),1);
	      
	      TauMuPtResolution.at(t).Fill(TruthTauMu.Pt()-EventFitTauMu.Pt(),1);
	      TauMuEResolution.at(t).Fill(TruthTauMu.E()-EventFitTauMu.E(),1);
	      TauMuPhiResolution.at(t).Fill(TruthTauMu.Phi()-EventFitTauMu.Phi(),1);
	      TauMuEtaResolution.at(t).Fill(TruthTauMu.Eta()-EventFitTauMu.Eta(),1);


	      TauA1PhiResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauMuPhiResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.Phi() -EventFitTauMu.Phi(),1);
	      TauA1EtaResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.Eta() -EventFitTauMu.Eta(),1);
	      TauA1PtResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1EResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauA1.E() - EventFitTauA1.E(),1);
	      TauMuEResolVsProb.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),TruthTauMu.E() - EventFitTauMu.E(),1);
	    

	      TauA1PtResolVsPt.at(t).Fill(A1Reco.Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsPt.at(t).Fill(MuReco.Pt(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1PtResolVsEta.at(t).Fill(A1Reco.Eta(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsEta.at(t).Fill(MuReco.Eta(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1EtaResolVsEta.at(t).Fill(A1Reco.Eta(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsEta.at(t).Fill(MuReco.Eta(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);


	      TauA1PtResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Pt()-EventFitTauMu.Pt(),1);

	      TauA1EtaResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsZPt.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);


 	      std::cout<<"theTauMuTLV  "<<theTauMuTLV.size()<<std::endl;
 	      std::cout<<"theTauA1TLV  "<<theTauA1TLV.size()<<std::endl;


	      std::cout<<" (TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px()  "<< TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px() <<"   " <<sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(3,3))<<std::endl;
	      std::cout<<" (TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px()  "<< TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px() <<"   " <<theTauMuTLV.at(AmbiguitySolution).Covariance(3,3)<<std::endl;

 	      std::cout<<" EventFitTauA1.Px()   "<< EventFitTauA1.Px() <<"   " <<" EventFitTauA1.Py()   "<<EventFitTauA1 .Py() <<"   " <<" EventFitTauA1.Pz()   "<< EventFitTauA1.Pz() <<"   " <<std::endl;
 	      std::cout<<" EventFitTauMu.Px()   "<< EventFitTauMu.Px() <<"   " <<" EventFitTauMu.Py()   "<< EventFitTauMu.Py() <<"   " <<" EventFitTauMu.Pz()   "<< EventFitTauMu.Pz() <<"   " <<std::endl;
 	      std::cout<<" EventFit Zmass    "<< (EventFitTauMu +  EventFitTauA1).M()<< " pt1 " << (EventFitTauMu +  EventFitTauA1).Pt() <<" pt2 "  <<EventFitTauMu.Pt()  -  EventFitTauA1.Pt() <<std::endl;


 	      TauA1SFPxPull.at(t).Fill((TruthTauA1.Px() - theTauAmbiguityVector.at(AmbiguitySolution).LV().Px() )/sqrt(theTauAmbiguityVector.at(AmbiguitySolution).Covariance(3,3)),1);
 	      TauA1SFPyPull.at(t).Fill((TruthTauA1.Py() - theTauAmbiguityVector.at(AmbiguitySolution).LV().Py() )/sqrt(theTauAmbiguityVector.at(AmbiguitySolution).Covariance(4,4)),1);
 	      TauA1SFPzPull.at(t).Fill((TruthTauA1.Pz() - theTauAmbiguityVector.at(AmbiguitySolution).LV().Pz() )/sqrt(theTauAmbiguityVector.at(AmbiguitySolution).Covariance(5,5)),1);


 	      TauMuPxPull.at(t).Fill((TruthTauMu.Px() - theTauMuTLV.at(AmbiguitySolution).LV().Px() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(3,3)),1);
 	      TauMuPyPull.at(t).Fill((TruthTauMu.Py() - theTauMuTLV.at(AmbiguitySolution).LV().Py() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(4,4)),1);
 	      TauMuPzPull.at(t).Fill((TruthTauMu.Pz() - theTauMuTLV.at(AmbiguitySolution).LV().Pz() )/sqrt(theTauMuTLV.at(AmbiguitySolution).Covariance(5,5)),1);

 	      TauA1PxPull.at(t).Fill((TruthTauA1.Px() - theTauA1TLV.at(AmbiguitySolution).LV().Px() )/sqrt(theTauA1TLV.at(AmbiguitySolution).Covariance(3,3)),1);
 	      TauA1PyPull.at(t).Fill((TruthTauA1.Py() - theTauA1TLV.at(AmbiguitySolution).LV().Py() )/sqrt(theTauA1TLV.at(AmbiguitySolution).Covariance(4,4)),1);
 	      TauA1PzPull.at(t).Fill((TruthTauA1.Pz() - theTauA1TLV.at(AmbiguitySolution).LV().Pz() )/sqrt(theTauA1TLV.at(AmbiguitySolution).Covariance(5,5)),1);



	      TauA1PhiResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauA1EtaResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauA1PtResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauA1EResolVsPVSVSignificance.at(t).Fill(flightsignificance,TruthTauA1.E() - EventFitTauA1.E(),1);
	    }
	    if(AmbPoint){

	      ProbabilityOfAmbiguityPoint.at(t).Fill(Chi2ProbabilityEventFit.at(AmbiguitySolution),1);

	      TauA1PtResolVsPtAmbPoint.at(t).Fill(A1Reco.Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsPtAmbPoint.at(t).Fill(MuReco.Pt(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1PtResolVsEtaAmbPoint.at(t).Fill(A1Reco.Eta(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsEtaAmbPoint.at(t).Fill(MuReco.Eta(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);
	      TauA1EtaResolVsEtaAmbPoint.at(t).Fill(A1Reco.Eta(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsEtaAmbPoint.at(t).Fill(MuReco.Eta(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);

	      TauA1PtResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauMuPtResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Pt() -EventFitTauMu.Pt(),1);

	      TauA1EtaResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauMuEtaResolVsZPtAmbPoint.at(t).Fill((TruthTauA1+TruthTauMu).Pt(),TruthTauMu.Eta() - EventFitTauMu.Eta(),1);

	      TauA1PtResolutionAmbPoint.at(t).Fill(TruthTauA1.Pt()-EventFitTauA1.Pt() ,1);
	      TauA1PhiResolutionAmbPoint.at(t).Fill(TruthTauA1.Phi()-EventFitTauA1.Phi(),1);
	      TauA1EtaResolutionAmbPoint.at(t).Fill(TruthTauA1.Eta()-EventFitTauA1.Eta(),1);
	      
	      TauMuPtResolutionAmbPoint.at(t).Fill(TruthTauMu.Pt()-EventFitTauMu.Pt(),1);
	      TauMuPhiResolutionAmbPoint.at(t).Fill(TruthTauMu.Phi()-EventFitTauMu.Phi(),1);
	      TauMuEtaResolutionAmbPoint.at(t).Fill(TruthTauMu.Eta()-EventFitTauMu.Eta(),1);

	      TauA1PhiResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauA1EtaResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauA1PtResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauA1EResolVsPhiRotSignificance.at(t).Fill(SignificanceOfTauRotation.at(AmbiguitySolution),TruthTauA1.E() - EventFitTauA1.E(),1);


	      TauA1PhiResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.Phi() - EventFitTauA1.Phi(),1);
	      TauA1EtaResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.Eta() - EventFitTauA1.Eta(),1);
	      TauA1PtResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.Pt() - EventFitTauA1.Pt(),1);
	      TauA1EResolVsPVSVSignificanceAmbPoint.at(t).Fill(flightsignificance,TruthTauA1.E() - EventFitTauA1.E(),1);

	    }
	  }
      }
      




      if(Tau_FitOk.at(1) && Tau_FitOk.at(2)){
	if(EventFit_Ok.at(1) && EventFit_Ok.at(2)){
	  if(Ntp->CheckDecayID(2,5)){
	    TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	    TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	    

	    int CorrectAmbiguity; int IncorrectAmbiguity;
	    if(fabs(TauA1_EF.at(1).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(2).E() - TruthTauA1.E())){  CorrectAmbiguity=1; IncorrectAmbiguity =2;}
	    if(fabs(TauA1_EF.at(2).E() - TruthTauA1.E()) < fabs(TauA1_EF.at(1).E() - TruthTauA1.E())){  CorrectAmbiguity=2;IncorrectAmbiguity =1; }
	    //----------------
	    Chi2Dim.at(t).Fill(TMath::Prob(Chi2EventFit.at(CorrectAmbiguity),1),TMath::Prob(Chi2EventFit.at(IncorrectAmbiguity),1),1);
	    csum2Dim.at(t).Fill(csum_GF.at(CorrectAmbiguity),csum_GF.at(IncorrectAmbiguity),1);
	    csumProb2Dim.at(t).Fill(csum_GF.at(CorrectAmbiguity),Chi2ProbabilityEventFit.at(CorrectAmbiguity),1);
	    if((TruthTauA1 + TruthTauMu).Pt() < 3)	    ProbabilityOfCorrect.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),1);



	    ProbabilityOfCorrectVsChi2.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),Chi2EventFit.at(CorrectAmbiguity),1);
	    ProbabilityOfCorrectndf2.at(t).Fill(Chi2ProbabilityEventFitndf2.at(CorrectAmbiguity),1);
	    ProbabilityOfCorrectndf3.at(t).Fill(Chi2ProbabilityEventFitndf3.at(CorrectAmbiguity),1);
	    Chi2OfCorrect.at(t).Fill(Chi2EventFit.at(CorrectAmbiguity),1);
	    std::cout<<" chiw EFVaildation code --->   "<< Chi2EventFit.at(CorrectAmbiguity) << std::endl;
	    ProbabilityOfCorrectVsZPt.at(t).Fill(Chi2ProbabilityEventFit.at(CorrectAmbiguity),(TruthTauA1 + TruthTauMu).Pt(),1);
	  }
	}
      }


    }
  }
}




void  EFValidation::Finish(){
  unsigned int t=1;


 
  for(int iBin = 1; iBin < TauA1PtResolVsPtRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsPt.at(t).GetNbinsX();
    TauA1PtResolVsPtRMS.at(t).SetBinContent(iBin,TauA1PtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMS());
    TauA1PtResolVsPtRMS.at(t).SetBinError(iBin,TauA1PtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMSError());
  }


  for(int iBin = 1; iBin < TauMuPtResolVsPtRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauMuPtResolVsPt.at(t).GetNbinsX();
    TauMuPtResolVsPtRMS.at(t).SetBinContent(iBin,TauMuPtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMS());
    TauMuPtResolVsPtRMS.at(t).SetBinError(iBin,TauMuPtResolVsPt.at(t).ProjectionY(" ",iBin, NMax)->GetRMSError());
  }

  for(int iBin = 1; iBin < TauA1PtResolVsEtaRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsEta.at(t).GetNbinsX();
    TauA1PtResolVsEtaRMS.at(t).SetBinContent(iBin,TauA1PtResolVsEta.at(t).ProjectionY(" ",iBin, iBin )->GetRMS());
    TauA1PtResolVsEtaRMS.at(t).SetBinError(iBin,TauA1PtResolVsEta.at(t).ProjectionY(" ",iBin, iBin)->GetRMSError());
  }


  for(int iBin = 1; iBin < TauMuPtResolVsZPtRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauMuPtResolVsZPt.at(t).GetNbinsX();
    TauMuPtResolVsZPtRMS.at(t).SetBinContent(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetRMS());
    TauMuPtResolVsZPtRMS.at(t).SetBinError(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetRMSError());
  }

  for(int iBin = 1; iBin < TauMuPtResolVsZPtMean.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauMuPtResolVsZPt.at(t).GetNbinsX();
    TauMuPtResolVsZPtMean.at(t).SetBinContent(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetMean());
    TauMuPtResolVsZPtMean.at(t).SetBinError(iBin,TauMuPtResolVsZPt.at(t).ProjectionY(" ",1, iBin)->GetMeanError());
  }



  for(int iBin = 1; iBin < TauA1PtResolVsPVSVSignificanceRMS.at(t).GetNbinsX() + 1; iBin++){
    double NMax = TauA1PtResolVsPVSVSignificance.at(t).GetNbinsX();
    TauA1PtResolVsPVSVSignificanceRMS.at(t).SetBinContent(iBin,TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",iBin, NMax)->GetRMS());
    TauA1PtResolVsPVSVSignificanceRMS.at(t).SetBinError(iBin,TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",iBin, NMax)->GetRMSError());
  }


   for(int iBin = 1; iBin < TauA1PtResolVsPhiRotSignificance.at(t).GetNbinsX() + 1; iBin++){
     double NMax = TauA1PtResolVsPhiRotSignificance.at(t).GetNbinsX();
     
     double  nom = TauA1PtResolVsPhiRotSignificance.at(t).ProjectionY(" ",0, iBin)->GetEntries();
     double  denom = TauA1PtResolVsPhiRotSignificance.at(t).ProjectionY(" ",0, NMax)->GetEntries();


     if(denom!=0){     PhiRotSignificanceCutEfficiency.at(t).SetBinContent(iBin,nom/denom);}
     else{PhiRotSignificanceCutEfficiency.at(t).SetBinContent(iBin,0);}

     
   }

    for(int iBin = 1; iBin < TauA1PtResolVsPVSVSignificance.at(t).GetNbinsX() + 1; iBin++){

      double NMax = TauA1PtResolVsPVSVSignificance.at(t).GetNbinsX();


      double  nom = TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",iBin,NMax)->GetEntries();
      double  denom = TauA1PtResolVsPVSVSignificance.at(t).ProjectionY(" ",0, NMax)->GetEntries();

      if(denom!=0){     PVSVSignificanceCutEfficiency.at(t).SetBinContent(iBin,nom/denom);}
      else{PVSVSignificanceCutEfficiency.at(t).SetBinContent(iBin,0);}
      
    }



    ChannelEfficiency.at(1).Divide(&idPassedEF.at(1), &idAllTaus.at(1));
    ChannelEfficiency.at(2).Divide(&idPassedEF.at(2), &idAllTaus.at(2));
    ChannelEfficiency.at(3).Divide(&idPassedEF.at(3), &idAllTaus.at(3));
    ChannelEfficiency.at(4).Divide(&idPassedEF.at(4), &idAllTaus.at(4));
    EfficiencyOverA1Pt.at(t).Divide(&TruthA1PtAfterFit.at(t),&TruthA1Pt.at(t));
    EfficiencyOverA1Phi.at(t).Divide(&TruthA1PhiAfterFit.at(t),&TruthA1Phi.at(t));
    EfficiencyOverA1Eta.at(t).Divide(&TruthA1EtaAfterFit.at(t),&TruthA1Eta.at(t));


    EfficiencyOverMuPt.at(t).Divide(&TruthMuPtAfterFit.at(t),&TruthMuPt.at(t));
    EfficiencyOverMuPhi.at(t).Divide(&TruthMuPhiAfterFit.at(t),&TruthMuPhi.at(t));
    EfficiencyOverMuEta.at(t).Divide(&TruthMuEtaAfterFit.at(t),&TruthMuEta.at(t));






  Selection::Finish();
}



