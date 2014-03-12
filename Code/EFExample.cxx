// Code written by Vladimir Cherepanov
// RWTH Aachen 
#include "EFExample.h"
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

EFExample::EFExample(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
		//  ,ran()
{
  verbose=false;
}

EFExample::~EFExample(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "EFExample::~EFExample Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "EFExample::~EFExample()" << std::endl;
}

void  EFExample::Configure(){
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

  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  // Setup Extra Histograms

  TauA1PtPhysical=HConfig.GetTH1D(Name+"_TauA1PtPhysical","",50,10,60,"pT_{#tau}, (GeV)","");
  TauA1PtAmbiguityPoint=HConfig.GetTH1D(Name+"_TauA1PtAmbiguityPoint","",50,10,60,"pT_{#tau}, (GeV)","");



  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);

}



void  EFExample::Store_ExtraDist(){


 Extradist1d.push_back(&TauA1PtPhysical);
 Extradist1d.push_back(&TauA1PtAmbiguityPoint);

}

void  EFExample::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){/* std::cout << "failed to find id" <<std::endl;*/ return;}

  if(verbose)std::cout << "void  EFExample::doEvent() Mu" << std::endl;



  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("eta2p1_LooseIsoPFTau"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);


  // Apply  basic Selection
  std::vector<unsigned int> mu_idx_good, mu_idx_pt, mu_idx_iso;
  unsigned mu_idx(999);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
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
	  tau_idx_good.push_back(i);
	}
      }
    }
  }



  value.at(hasTau)=tau_idx_good.size();
  pass.at(hasTau)=(value.at(hasTau)>cut.at(hasTau));


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
    if (tau_idx!=999 && mu_idx!=999){

      std::vector<bool> Tau_FitOk; 
      std::vector<double> Tau_FitChi2;
      std::vector<bool> EventFit_Ok;
      std::vector<TLorentzVector> TauA1_EF;
      std::vector<TLorentzVector> TauMu_EF;
      
      std::vector<double> Chi2EventFit;
      std::vector<double> Chi2ProbabilityEventFit;

      TVector3 pv;
      TMatrixTSym<double> PVcov;
      TVector3 sv;
      TMatrixTSym<double> SVcov;
      
      
      //////////////////////////////
      //loop over ambiguity points
      for(unsigned int j=0;j< MultiProngTauSolver::NAmbiguity;j++){
	
	LorentzVectorParticle theTau;
	LorentzVectorParticle theZ;
	std::vector<LorentzVectorParticle> daughter;
	std::vector<LorentzVectorParticle> theZdaughter;
	       
	TLorentzVector TauA1EventFit(0,0,0,0);
	TLorentzVector TauA1ThreeProngFit(0,0,0,0);
	TLorentzVector TauMuEventFit(0,0,0,0);
	double LC_Eventchi2(0);
	double LC_Eventchi2Probability(0);
	int NiterationsEF(0);
	double csum(0);
	double LC_chi2(0);

	bool  ThreeProngFitSuccess =false;
	bool  EventFitSuccess =false;
	ThreeProngFitSuccess=Ntp->ThreeProngTauFit(tau_idx,j,theTau,daughter,LC_chi2);
	if(ThreeProngFitSuccess){
   	  TauA1ThreeProngFit=theTau.LV();
  
	  EventFitSuccess = Ntp->EventFit(tau_idx,mu_idx,theTau,theZ,theZdaughter,LC_Eventchi2,NiterationsEF,csum);
	  if(EventFitSuccess){
	    LC_Eventchi2Probability = TMath::Prob(LC_Eventchi2,1);
	    TauA1EventFit =theZdaughter.at(0).LV();
	    TauMuEventFit =theZdaughter.at(1).LV();

	  }
	}


	Tau_FitOk.push_back(ThreeProngFitSuccess);
	Tau_FitChi2.push_back(LC_chi2);
  
	Chi2ProbabilityEventFit.push_back(LC_Eventchi2Probability);

	EventFit_Ok.push_back(EventFitSuccess);
 	TauA1_EF.push_back(TauA1EventFit);
 	TauMu_EF.push_back(TauMuEventFit);
      }
      
      int AmbiguitySolution =0;
      bool AmbPoint(false);
      if(Ntp->AmbiguitySolver(Tau_FitOk,EventFit_Ok,Chi2ProbabilityEventFit,AmbiguitySolution, AmbPoint)){
	
	TLorentzVector EventFitTauA1 = TauA1_EF.at(AmbiguitySolution);
	TLorentzVector EventFitTauMu = TauMu_EF.at(AmbiguitySolution);


	if(!AmbPoint){
	  TauA1PtPhysical.at(t).Fill(EventFitTauA1.Pt(),1);
	}
	
	if(AmbPoint){
	  TauA1PtAmbiguityPoint.at(t).Fill(EventFitTauA1.Pt(),1);
	}
      }

    }
  }
}



