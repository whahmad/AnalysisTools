#include "ZtoEMu_Skim.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

#include "Parameters.h"
#include "TMath.h"

#include <TFile.h>

ZtoEMu_Skim::ZtoEMu_Skim(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(10)
  ,e_pt(10)
  ,mu_eta(2.1)
  ,e_eta(2.3)
  ,jet_pt(30)
  ,jet_eta(2.4)
  ,jet_sum(70)
  ,zmin(88)
  ,zmax(94)
{
    //verbose=true;
	MVA_ID = true;
	FRFile = new TFile("/net/scratch_cms/institut_3b/nehrkorn/FakeRates_2012_19ifb_rereco.root");
	EmbEffFile = new TFile("/net/scratch_cms/institut_3b/nehrkorn/RecHitElectronEfficiencies.root");
	ElectronFakeRate = (TH2D*)(FRFile->Get("ElectronFakeRateHist"));
	MuonFakeRate = (TH2D*)(FRFile->Get("MuonFakeRateHist"));
	EmbEff = (TH2D*)(EmbEffFile->Get("hPtEtaSFL"));

}

ZtoEMu_Skim::~ZtoEMu_Skim(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu_Skim::~ZtoEMu_Skim Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu_Skim::~ZtoEMu_Skim()" << std::endl;
}

void  ZtoEMu_Skim::Configure(){

  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==NMu)                cut.at(NMu)=1;
    if(i==NMuPt)              cut.at(NMuPt)=1;
    if(i==NMuEta)             cut.at(NMuEta)=1;
    if(i==NE)                 cut.at(NE)=1;
    if(i==NEPt)               cut.at(NEPt)=1;
    if(i==NEEta)              cut.at(NEEta)=1;
    if(i==SameVtx)            cut.at(SameVtx)=true;
    if(i==charge)             cut.at(charge)=0;
	if(i==MtMu)               cut.at(MtMu)=40;
	if(i==drMuE)              cut.at(drMuE)=0.2;
	if(i==qualitycuts)        cut.at(qualitycuts)=true;
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
  
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>=$";
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
    }

    else if(i==NMu){
      title.at(i)="Number $\\mu =$";
      title.at(i)+=cut.at(NMu);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NE){
      title.at(i)="Number $e =$";
      title.at(i)+=cut.at(NE);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NMuPt){
      title.at(i)="Number of $\\mu$ [$P_{T}^{\\mu}>$";
      title.at(i)+=mu_pt;
      title.at(i)+="(GeV)] $>=$";
      title.at(i)+=cut.at(NMuPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NMuPt_","P_{T,#mu} (N-1 Distribution)",100,0,200,"P_{T,#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NMuPt_","P_{T,#mu} (Accumulative Distribution)",100,0,200,"P_{T,#mu} (GeV)","Events");
    }
    else if(i==NEPt){
      title.at(i)="Number of $e$ [$P_{T}^{e}>$";
      title.at(i)+=e_pt;
      title.at(i)+="(GeV)] $>=$";
      title.at(i)+=cut.at(NEPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NEPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NEPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NEPt_","P_{T,e} (N-1 Distribution)",100,0,200,"P_{T,e} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NEPt_","P_{T,e} (Accumulative Distribution)",100,0,200,"P_{T,e} (GeV)","Events");
    }
    else if(i==NMuEta){
      title.at(i)="Number of $\\mu$ [$|\\eta^{\\mu}|<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",mu_eta);
      title.at(i)+=buffer;
      title.at(i)+="] $>=$";
      title.at(i)+=cut.at(NMuEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NMuEta_","#eta{#mu} (N-1 Distribution)",100,-7,7,"#eta_{#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NMuEta_","#eta_{#mu} (Accumulative Distribution)",100,-7,7,"#eta_{#mu} (GeV)","Events");
    }
    else if(i==NEEta){
      title.at(i)="Number of $e$ [$|\\eta^{e}|<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",e_eta);
      title.at(i)+=buffer;
      title.at(i)+="] $>=$";
      title.at(i)+=cut.at(NEEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NEEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NEEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NEEta_","#eta{e} (N-1 Distribution)",100,-7,7,"#eta_{e} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NEEta_","#eta_{e} (Accumulative Distribution)",100,-7,7,"#eta_{e} (GeV)","Events");
    }
    else if(i==charge){
      title.at(i)="$e-\\mu$ Charge = ";
      title.at(i)+=cut.at(charge);
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="e-#mu Charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,21,-5.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,21,-5.5,5.5,hlabel,"Events"));
    }
	else if(i==MtMu){
      title.at(i)="$m_{T}^{\\mu,Miss} < $";
      title.at(i)+=cut.at(MtMu);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="m_{T}^{#mu,Miss} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MtMu_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MtMu_",htitle,40,0,200,hlabel,"Events"));
    }
	else if(i==drMuE){
	  title.at(i)="$dR(e,\\mu) > $";
	  char buffer[50];
	  sprintf(buffer,"%5.2f",cut.at(drMuE));
	  title.at(i)+=buffer;
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="dR(e,#mu)";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_drMuE_",htitle,400,0.,10.,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_drMuE_",htitle,400,0.,10.,hlabel,"Events"));
	}
	else if(i==SameVtx){
	  title.at(i)="$e$ and $\\mu$ from same vtx";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SameVtx_",htitle,20,0,0.1,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SameVtx_",htitle,20,0,0.1,hlabel,"Events"));
	}
    else if(i==qualitycuts){
  	  title.at(i)="quality cuts";
  	  htitle=title.at(i);
  	  htitle.ReplaceAll("$","");
  	  htitle.ReplaceAll("\\","#");
  	  hlabel="";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_qualitycuts_",htitle,20,0,0.1,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_qualitycuts_",htitle,20,0,0.1,hlabel,"Events"));
  	}
    
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_Skim::Store_ExtraDist(){
}

void  ZtoEMu_Skim::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
  
	double embedd_weight = 1.;
	if(id==34){
		embedd_weight = Ntp->EmbeddedWeight();
		if(Ntp->EmbeddedWeight()!=Ntp->TauSpinnerWeight()*Ntp->SelEffWeight()*Ntp->RadiationCorrWeight()*Ntp->MinVisPtFilter()*Ntp->KinWeightPt()*Ntp->KinWeightEta()*Ntp->KinWeightMassPt()){
			std::cout << "Product = " << Ntp->TauSpinnerWeight()*Ntp->SelEffWeight()*Ntp->RadiationCorrWeight()*Ntp->MinVisPtFilter()*Ntp->KinWeightPt()*Ntp->KinWeightEta()*Ntp->KinWeightMassPt() << std::endl;
			std::cout << "Embedded weight = " << embedd_weight << std::endl;
		}
	}

  value.at(TriggerOk)=0;
  /*std::cout << "Number of HLT Triggers:" << Ntp->NHLTTriggers() << std::endl;
  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
	  std::cout << "######################################" << std::endl;
	  std::cout << "Trigger: " << Ntp->HTLTriggerName(i) << std::endl;
	  std::cout << "Acceptet (1) or not (0): " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << std::endl;
	  std::cout << "######################################" << std::endl;
  }*/
  if(Ntp->TriggerAccept("HLT_IsoMu24_eta2p1"))value.at(TriggerOk)=1;
  if(Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloId")
		  && !Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloId")
		  ){
	  value.at(TriggerOk)=1;
	  mu_pt = 20;
	  e_pt = 10;
  }
  if(Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloId")
		  && !Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloId")
		  ){
	  value.at(TriggerOk)=1;
	  mu_pt = 10;
	  e_pt = 20;
  }
  if(Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloId")
		  && Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloId")
		  ){
	  value.at(TriggerOk)=1;
	  mu_pt = 20;
	  e_pt = 20;
  }
  if(id==34){
	  value.at(TriggerOk)=1;
	  mu_pt = 10;
	  e_pt = 10;
  }
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));

  // Apply Selection
  if(verbose) std::cout << "find primary vertex" << std::endl;
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));

  unsigned int pVtx(999);
  if(nGoodVtx>=cut.at(PrimeVtx))pVtx=nGoodVtx;
  
  ///////////////////////////////////////////////
  //
  // Quality cuts
  //
  if(verbose)std::cout << "Quality cuts" << std::endl;
  value.at(qualitycuts)=false;
  std::vector<int> qualitymuons;
  std::vector<int> qualityelectrons;
  qualitymuons.clear();
  qualityelectrons.clear();
  
  // loop over muons & electrons and save the ones with 30<=pt<=70 passing medium object ID
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muons_p4(i).Pt()>=mu_pt
			  && Ntp->Muons_p4(i).Pt()<=70.
			  && isFakeMuon(i)
			  ){
		  qualitymuons.push_back(i);
	  }
  }
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>=e_pt
			  && Ntp->Electron_p4(i).Et()<=70.
			  ){
		  if(fabs(Ntp->Electron_supercluster_eta(i))<0.8 && Ntp->Electron_MVA_discriminator(i)>0.7){
			  qualityelectrons.push_back(i);
		  }else if(fabs(Ntp->Electron_supercluster_eta(i))>=0.8 && fabs(Ntp->Electron_supercluster_eta(i))<1.479 && Ntp->Electron_MVA_discriminator(i)>0.9){
			  qualityelectrons.push_back(i);
		  }else if(fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && Ntp->Electron_MVA_discriminator(i)>0.8375){
			  qualityelectrons.push_back(i);
		  }
	  }
  }
  // require exactly one such muon & electron
  if(qualitymuons.size()==1 && qualityelectrons.size()==1)value.at(qualitycuts)=true;
  pass.at(qualitycuts)=(value.at(qualitycuts));
  
  ///////////////////////////////////////////////
  //
  // Vertex constraint
  //
  if(verbose)std::cout << "Vertex constraint" << std::endl;
  int pos = 0;
  int posmu = 0;
  int pose = 0;
  int mutrack=0;
  int etrack=0;
  value.at(SameVtx)=false;
  
  // loop over all tracks and find the one with lowest dR(e/mu,track). this track number will be saved
  if(qualitymuons.size()>0){
	  for(unsigned i=0;i<Ntp->NTracks();i++){
		  if(Ntp->Track_p4(i).DeltaR(Ntp->Muons_p4(qualitymuons.at(0)))<Ntp->Track_p4(mutrack).DeltaR(Ntp->Muons_p4(qualitymuons.at(0))))mutrack=i;
	  }
  }
  if(qualityelectrons.size()>0){
	  for(unsigned i=0;i<Ntp->NTracks();i++){
		  if(Ntp->Track_p4(i).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))<Ntp->Track_p4(etrack).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))){
			  if(Ntp->Track_p4(i).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))<0.1 || Ntp->Track_p4(etrack).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0)))<0.1){
				  etrack=i;
			  }
		  }
	  }
  }
  
  // fall-back solution in case electron & muon are assigned to same track
  if(qualitymuons.size()>0 && qualityelectrons.size()>0 && mutrack==etrack){
	  if(verbose){
		  std::cout << "muon and electron from same track " << mutrack << std::endl;
		  std::cout << "dR(e,mu) = " << Ntp->Muons_p4(qualitymuons.at(0)).DeltaR(Ntp->Electron_p4(qualityelectrons.at(0))) << std::endl;
	  }
	  mutrack=-1;
	  etrack=-1;
  }
  
  if(verbose)std::cout << "looping over vertices" << std::endl;
  int muvtx=-1;
  int evtx=-1;
  // loop over vertices & check which vertex has tracks from muon/electron assigned to them
  for(unsigned i=0;i<Ntp->NVtx();i++){
	  for(unsigned j=0;j<Ntp->Vtx_Track_idx(i).size();j++){
		  //if(Ntp->Vtx_Track_idx(i).at(j)!=-1)std::cout << "Vtx track id = " << Ntp->Vtx_Track_idx(i).at(j) << std::endl;
		  if(mutrack>=0 && mutrack==Ntp->Vtx_Track_idx(i).at(j)){
			  muvtx=i;
		  }
		  if(etrack>=0 && etrack==Ntp->Vtx_Track_idx(i).at(j)){
			  evtx=i;
		  }
	  }
  }
  
  // check if muon & electron are assigned to the same vertex
  if(muvtx==evtx && muvtx!=-1 && evtx!=-1){
	  value.at(SameVtx)=true;
  }
  // if one track was not assigned to a vertex, use dxy/dz as fall-back solution
  else if(muvtx==-1 && evtx!=-1){
	  if(mutrack!=-1 && qualitymuons.size()>0){
		  if(dxy(Ntp->Muons_p4(qualitymuons.at(0)),Ntp->Muon_Poca(qualitymuons.at(0)),Ntp->Vtx(evtx))<0.2){
			  if(dz(Ntp->Muons_p4(qualitymuons.at(0)),Ntp->Muon_Poca(qualitymuons.at(0)),Ntp->Vtx(evtx))<0.5){
				  muvtx=evtx;
				  value.at(SameVtx)=true;
			  }
		  }
	  }
  }
  else if(evtx==-1 && muvtx!=-1){
	  if(etrack!=-1 && qualityelectrons.size()>0){
		  if(dxy(Ntp->Electron_p4(qualityelectrons.at(0)),Ntp->Electron_Poca(qualityelectrons.at(0)),Ntp->Vtx(muvtx))<0.02){
			  if(dz(Ntp->Electron_p4(qualityelectrons.at(0)),Ntp->Electron_Poca(qualityelectrons.at(0)),Ntp->Vtx(muvtx))<0.1){
				  evtx=muvtx;
				  value.at(SameVtx)=true;
			  }
		  }
	  }
  }
  pass.at(SameVtx)=(value.at(SameVtx));
  if(value.at(SameVtx)){
	  pos=muvtx;
	  posmu=muvtx;
	  pose=evtx;
  }else{
	  pos=0;
	  if(muvtx!=-1){
		  posmu=muvtx;
	  }else{
		  posmu=0;
	  }
	  if(evtx!=-1){
		  pose=evtx;
	  }else{
		  pose=0;
	  }
  }

  ///////////////////////////////////////////////
  //
  // Mu Cuts
  //
  if(verbose) std::cout << "Muon cuts" << std::endl;
  std::vector<unsigned int> GoodMuons;
  bool fakeMuon = false;

  // muon ID cuts (including eta-dependent isolation)
  for(unsigned i=0;i<qualitymuons.size();i++){
	  if(isTightMuon(qualitymuons.at(i),posmu) &&
			  dxy(Ntp->Muons_p4(qualitymuons.at(i)),Ntp->Muon_Poca(qualitymuons.at(i)),Ntp->Vtx(posmu))<0.02 &&
			  dz(Ntp->Muons_p4(qualitymuons.at(i)),Ntp->Muon_Poca(qualitymuons.at(i)),Ntp->Vtx(posmu))<0.2
			  ){
		  if(fabs(Ntp->Muons_p4(qualitymuons.at(i)).Eta())<1.479 && Muon_RelIso(qualitymuons.at(i))<0.15){
			  GoodMuons.push_back(qualitymuons.at(i));
		  }else if(fabs(Ntp->Muons_p4(qualitymuons.at(i)).Eta())>=1.479 && Muon_RelIso(qualitymuons.at(i))<0.1){
			  GoodMuons.push_back(qualitymuons.at(i));
		  }
	  }else if(isFakeMuon(qualitymuons.at(i),posmu) && Ntp->isData() && MVA_ID){
		  fakeMuon = true;
		  GoodMuons.push_back(qualitymuons.at(i));
		  fakeRateMu = Fakerate(Ntp->Muons_p4(qualitymuons.at(i)),MuonFakeRate,"muon");
	  }
  }

  // muon pt cut
  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuPt).push_back(Ntp->Muons_p4(GoodMuons.at(i)).Pt());
    if(Ntp->Muons_p4(GoodMuons.at(i)).Pt()<mu_pt){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  
  value.at(NMuPt)=GoodMuons.size();
  pass.at(NMuPt)=(value.at(NMuPt)>=cut.at(NMuPt));
  
  // muon eta cut
  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuEta).push_back(Ntp->Muons_p4(GoodMuons.at(i)).Eta());
    if(fabs(Ntp->Muons_p4(GoodMuons.at(i)).Eta())>mu_eta){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  
  value.at(NMuEta)=GoodMuons.size();
  pass.at(NMuEta)=(value.at(NMuEta)>=cut.at(NMuEta));
  
  value.at(NMu)=GoodMuons.size();
  pass.at(NMu)=(value.at(NMu)==cut.at(NMu));
  
  unsigned int muidx1(999),muInfoidx1(999);
  if(GoodMuons.size()>=1){muidx1=GoodMuons.at(0);muInfoidx1=muidx1;}
  if(verbose)std::cout << "void  ZtoEMu_Skim::doEvent() E " << muidx1 <<" "<< muInfoidx1 << std::endl;

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  bool fakeElectron = false;
  
  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<qualityelectrons.size();i++){
	  if(!Ntp->Electron_HasMatchedConversions(qualityelectrons.at(i)) &&
			  Ntp->Electron_numberOfMissedHits(qualityelectrons.at(i))==0 &&
			  dz(Ntp->Electron_p4(qualityelectrons.at(i)),Ntp->Electron_Poca(qualityelectrons.at(i)),Ntp->Vtx(pose))<0.2 &&
			  dxy(Ntp->Electron_p4(qualityelectrons.at(i)),Ntp->Electron_Poca(qualityelectrons.at(i)),Ntp->Vtx(pose))<0.02
			  ){
		  if(MVA_ID){
			  if(isMVAElectron(qualityelectrons.at(i))){
				  GoodElectrons.push_back(qualityelectrons.at(i));
			  }else if(isFakeElectron(qualityelectrons.at(i),pose) && Ntp->isData()){
				  fakeElectron = true;
				  GoodElectrons.push_back(qualityelectrons.at(i));
				  fakeRateE = Fakerate(Ntp->Electron_p4(qualityelectrons.at(i)),ElectronFakeRate,"electron");
			  }
		  }else{
		  	if(isTightElectron(qualityelectrons.at(i),pose))GoodElectrons.push_back(qualityelectrons.at(i));
		  }
	  }
  }
  
  // electron pt cut
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEPt).push_back(Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
    if(Ntp->Electron_p4(GoodElectrons.at(i)).Pt()<e_pt){
      GoodElectrons.erase(GoodElectrons.begin()+i);
      i--;
    }
  }
  value.at(NEPt)=GoodElectrons.size();
  pass.at(NEPt)=(value.at(NEPt)>=cut.at(NEPt));
  
  // electron eta cut
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEEta).push_back(Ntp->Electron_supercluster_eta(GoodElectrons.at(i)));
    if(fabs(Ntp->Electron_supercluster_eta(GoodElectrons.at(i)))>e_eta){
      GoodElectrons.erase(GoodElectrons.begin()+i);
      i--;
    }
  }
  value.at(NEEta)=GoodElectrons.size();
  pass.at(NEEta)=(value.at(NEEta)>=cut.at(NEEta));
  
  value.at(NE)=GoodElectrons.size();
  pass.at(NE)=(value.at(NE)==cut.at(NE));
  
  unsigned int eidx1(999),eInfoidx1(999);
  if(GoodElectrons.size()>=1){eidx1=GoodElectrons.at(0);eInfoidx1=eidx1;}
  if(verbose)std::cout << "void  ZtoEMu_Skim::doEvent() E " << eidx1 <<" "<< eInfoidx1 << std::endl;
  
  ///////////////////////////////////////////////
  //
  // dR cleaning of e & mu
  //
  if(verbose)std::cout << "dR cleaning" << std::endl;
  
  value.at(drMuE)=10.;
  if(muidx1!=999 && eidx1!=999){
	  if(Ntp->Muon_Track_idx(GoodMuons.at(0))==Ntp->Electron_Track_idx(GoodElectrons.at(0))){
		  value.at(drMuE)=10.;
	  }
	  value.at(drMuE)=Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0)));
  }
  pass.at(drMuE)=(value.at(drMuE)>cut.at(drMuE));

  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=-5;
  if(eidx1!=999 && muidx1!=999){
	  value.at(charge)=Ntp->Electron_Charge(GoodElectrons.at(0))+Ntp->Muon_Charge(GoodMuons.at(0));
  }
  pass.at(charge)=(value.at(charge)==cut.at(charge));
  
  ///////////////////////////////////////////////
  //
  // Mt Mu cut
  //
  if(verbose) std::cout << "Mt Mu cut" << std::endl;
  value.at(MtMu)=999;
  if(muidx1!=999){
	  value.at(MtMu)=sqrt(2*Ntp->Muons_p4(GoodMuons.at(0)).Pt()*Ntp->MET_Corr_et()*(1-cosphi2d(Ntp->Muons_p4(GoodMuons.at(0)).Px(),Ntp->Muons_p4(GoodMuons.at(0)).Py(),Ntp->MET_Corr_ex(),Ntp->MET_Corr_ey())));
  }
  pass.at(MtMu)=(value.at(MtMu)<cut.at(MtMu));
  
  ////////////////////////////////////////////////
  //
  // QCD
  //
  if(MVA_ID){
	  fakeRate = 0.;
	  if(!pass.at(charge)
		  && Ntp->isData()
		  ){
		  if(fakeMuon && !fakeElectron){
			  fakeRate = fakeRateMu;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
		  }else if(fakeElectron && !fakeMuon){
			  fakeRate = fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
		  }else if(fakeMuon && fakeElectron){
			  fakeRate = fakeRateMu*fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
		  }else{
			  fakeRate = 1.;
		  }
	  }else if(pass.at(charge)
		  && Ntp->isData()
		  && !fakeMuon
		  && !fakeElectron
		  ){
		  fakeRate = 1.;
	  }
  }else{
	fakeRate = 1.;
  }
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=34){
    w*=Ntp->EvtWeight3D();
    if(verbose)std::cout << "void  ZtoEMu_Skim::doEvent() k" << w << " " << wobs << std::endl;
  }else if(Ntp->GetMCID()==34){
		w*=embedd_weight;
		if(pass.at(NE))w*=ElectronEffRecHit(GoodElectrons.at(0));
		if(pass.at(NMu))w*=MuonDataSF(GoodMuons.at(0));
  }
  else{w=1*fakeRate;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_Skim::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu_Skim::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu_Skim::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu_Skim::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_Skim::isTightMuon(unsigned int i){
	if(verbose)std::cout << "isTightMuon(unsigned int i)" << std::endl;
	if(Ntp->Muon_isGlobalMuon(i) &&
			Ntp->Muon_isPFMuon(i) &&
			Ntp->Muon_normChi2(i)<10.0 &&
			Ntp->Muon_hitPattern_numberOfValidMuonHits(i)>0 &&
			Ntp->Muon_numberOfMatchedStations(i)>1 &&
			Ntp->Muon_numberofValidPixelHits(i)>0 &&
			Ntp->Muon_trackerLayersWithMeasurement(i)>5
			  ){
		return true;
	}
	return false;
}

bool ZtoEMu_Skim::isTightMuon(unsigned int i, unsigned int j){
	if(isTightMuon(i) &&
			dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))<0.2 &&
			dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))<0.5
			  ){
		return true;
	}
	return false;
}

bool ZtoEMu_Skim::isLooseMuon(unsigned int i){
	if(Ntp->Muon_isPFMuon(i) &&
			(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))
			){
		return true;
	}
	return false;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int i){
	if(Ntp->Muon_isGlobalMuon(i) &&
			Ntp->Muons_p4(i).Pt()>10 &&
			fabs(Ntp->Muons_p4(i).Eta())<2.1
			){
		if(Ntp->Muons_p4(i).Pt()>20 &&
				Muon_RelIso(i)<0.4
				){
			return true;
		}else if(Ntp->Muons_p4(i).Pt()<=20 &&
				Muon_AbsIso(i)<8.
				){
			return true;
		}
	}
	return false;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int i, unsigned int j){
	if(isFakeMuon(i)
			&& dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))<0.2
			){
		return true;
	}
	return false;
}

double ZtoEMu_Skim::Muon_RelIso(unsigned int i){
	return (Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i)))/Ntp->Muons_p4(i).Pt();
}

double ZtoEMu_Skim::Muon_AbsIso(unsigned int i){
	return Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i));
}

//////////////////////////////
//
// Electron related functions
//

bool ZtoEMu_Skim::isMVAElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(mvapt<20){
		if(mvaeta<0.8 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_discriminator(i)>0.925){
			return true;
		}else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_discriminator(i)>0.915){
			return true;
		}else if(mvaeta>=1.479 && Electron_RelIso(i)<0.10 && Ntp->Electron_MVA_discriminator(i)>0.965){
			return true;
		}
	}else if(mvapt>=20){
		if(mvaeta<0.8 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_discriminator(i)>0.905){
			return true;
		}else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_discriminator(i)>0.955){
			return true;
		}else if(mvaeta>=1.479 && Electron_RelIso(i)<0.10 && Ntp->Electron_MVA_discriminator(i)>0.975){
			return true;
		}
	}
	return false;
}

bool ZtoEMu_Skim::isTightElectron(unsigned int i){
	if(verbose)std::cout << "isTightElectron(unsigned int i)" << std::endl;
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.442){ //barrel
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.004 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.03 &&
				Ntp->Electron_sigmaIetaIeta(i)<0.01 &&
				Ntp->Electron_hadronicOverEm(i)<0.12 &&
				fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))<0.05 &&
				Ntp->Electron_RelIso(i)<0.10 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_numberOfMissedHits(i)==0
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(i))>1.566 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){ //endcaps
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.005 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.02 &&
				Ntp->Electron_sigmaIetaIeta(i)<0.03 &&
				Ntp->Electron_hadronicOverEm(i)<0.10 &&
				fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))<0.05 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_numberOfMissedHits(i)==0
				){
			if(Ntp->Electron_p4(i).Pt()>=20.0 && Electron_RelIso(i)<0.10){
				return true;
			}else if(Ntp->Electron_p4(i).Pt()<20.0 && Electron_RelIso(i)<0.07){
				return true;
			}
		}
	}
	return false;
}

bool ZtoEMu_Skim::isTightElectron(unsigned int i, unsigned int j){
	if(isTightElectron(i)
			&& dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.02
			&& dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.1
			){
		return true;
	}
	return false;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.442){
		if(Ntp->Electron_p4(i).Pt()>20 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_sigmaIetaIeta(i)<0.01 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.15 &&
				Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.007 &&
				Electron_RelIso(i)<0.2
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(i))>=1.442 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_p4(i).Pt()>20 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_sigmaIetaIeta(i)<0.03 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.10 &&
				Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.009 &&
				Electron_RelIso(i)<0.2
				){
			return true;
		}
	}
	return false;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int i, unsigned int j){
	if(isFakeElectron(i)
			&& dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.1
			&& dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.03
			){
		return true;
	}
	return false;
}

double ZtoEMu_Skim::Electron_RelIso(unsigned int i){
	return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();
}

double ZtoEMu_Skim::Electron_Aeff(double Eta){
	double eta=fabs(Eta);
	if(eta<1.0)return 0.13;
	else if(eta>1.0 && eta<1.479)return 0.14;
	else if(eta>1.479 && eta<2.0)return 0.07;
	else if(eta>2.0 && eta<2.2)return 0.09;
	else if(eta>2.2 && eta<2.3)return 0.11;
	else if(eta>2.3 && eta<2.4)return 0.11;
	else if(eta>2.4)return 0.14;
}

//////////////////////////////
//
// Trigger & ID efficiencies
//

double ZtoEMu_Skim::MuonSF(unsigned int i){
	double pt = Ntp->Muons_p4(i).Pt();
	double eta = fabs(Ntp->Muons_p4(i).Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9829*0.9771;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9745*0.9746;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9943*0.9644;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9158*0.9891;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9850*0.9548;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9852*0.9701;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9743*0.9766;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9333*0.9892;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9951*0.9648;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9610*0.9836;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9716*0.9820;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9459*0.9909;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9869*0.9676;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9779*0.9817;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9665*0.9886;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9501*0.9883;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9959*0.9883;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9881*0.9833;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9932*0.9910;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9391*0.9900;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9986*0.9826;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9540*0.9841;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9549*0.9900;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9386*0.9886;
		}
	}
}

double ZtoEMu_Skim::MuonDataSF(unsigned int i){
	double pt = Ntp->Muons_p4(i).Pt();
	double eta = fabs(Ntp->Muons_p4(i).Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9701*0.5981;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9419*0.6578;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9303*0.6738;
		}else if(eta>=1.6 && eta<2.1){
			return 0.8623*0.6246;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9720*0.6740;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9305*0.7309;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9267*0.7416;
		}else if(eta>=1.6 && eta<2.1){
			return 0.8995*0.6954;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9764*0.7533;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9439*0.7915;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9366*0.7997;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9134*0.7567;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9725*0.8141;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9405*0.8364;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9218*0.8462;
		}else if(eta>=1.6 && eta<2.1){
			return 0.8824*0.8051;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9785*0.8606;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9342*0.8680;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9184*0.8745;
		}else if(eta>=1.6 && eta<2.1){
			return 0.8990*0.8399;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9679*0.9255;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9310*0.9249;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9092*0.9291;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9016*0.9025;
		}
	}
}

double ZtoEMu_Skim::ElectronSF(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9548*0.7654;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9015*0.7693;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9017*0.5719;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9830*0.8394;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9672*0.8457;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9463*0.7024;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9707*0.8772;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9731*0.8530;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9691*0.7631;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9768*0.9006;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9870*0.8874;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9727*0.8092;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 1.0047*0.9261;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9891*0.9199;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9858*0.8469;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 1.0063*0.9514;
		}else if(eta>=0.8 && eta<1.5){
			return 1.0047*0.9445;
		}else if(eta>=1.5 && eta<2.3){
			return 1.0015*0.9078;
		}
	}
}

double ZtoEMu_Skim::ElectronDataSF(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.7270*0.3436;
		}else if(eta>=0.8 && eta<1.5){
			return 0.7380*0.3481;
		}else if(eta>=1.5 && eta<2.3){
			return 0.6899*0.1104;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.8752*0.5196;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9059*0.5235;
		}else if(eta>=1.5 && eta<2.3){
			return 0.8635*0.2431;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9142*0.6442;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9484*0.5535;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9356*0.2888;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9368*0.7191;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9630*0.6472;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9466*0.3746;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9499*0.7819;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9642*0.7224;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9735*0.4527;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9689*0.8650;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9809*0.8201;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9802*0.6015;
		}
	}
}

double ZtoEMu_Skim::ElectronEffRecHit(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = Ntp->Electron_supercluster_eta(i);
	
	if(pt>199.99)pt=199.9;
	eta=fabs(eta);
	if(eta>2.49)eta=2.49;
	if(pt<10)return 0;
	
	Float_t eff=0;
	Int_t bin = EmbEff->FindFixBin(pt,eta);
	eff = EmbEff->GetBinContent(bin);
	
	return eff;
}

//////////////////////////////
//
// Calculate fakerate
//

double ZtoEMu_Skim::Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type){
	
	double eta1, eta2;
	int ptbin = 0;
	int etabin = 0;
	double FakePt = vec.Pt();
	double FakeEta = vec.Eta();
	double fakerate = 0.;
	
	if(type=="muon"){
		eta1 = 2.1;
		eta2 = 1.2;
	}else if(type=="electron"){
		eta1 = 2.;
		eta2 = 1.479;
	}
	
	if(FakePt<15.)ptbin = 1;
	else if(FakePt>=15. && FakePt<20.)ptbin = 2;
	else if(FakePt>=20. && FakePt<25.)ptbin = 3;
	else if(FakePt>=25. && FakePt<30.)ptbin = 4;
	else if(FakePt>=30.)ptbin = 5;
	
	if(FakeEta<-eta1)etabin = 1;
	else if(FakeEta>=-eta1 && FakeEta<-eta2)etabin = 2;
	else if(FakeEta>=-eta2 && FakeEta<-0.8)etabin = 3;
	else if(FakeEta>=-0.8 && FakeEta<0.)etabin = 4;
	else if(FakeEta>=0. && FakeEta<0.8)etabin = 5;
	else if(FakeEta>=0.8 && FakeEta<eta2)etabin = 6;
	else if(FakeEta>=eta2 && FakeEta<eta1)etabin = 7;
	else if(FakeEta>=eta1)etabin = 8;
	
	if(ptbin==0 || etabin==0){
		fakerate = 0;
	}else{
		fakerate = fakeRateHist->GetBinContent(ptbin,etabin);
	}
	
	return fakerate;
}

//////////////////////////////
//
// Finish function
//

void ZtoEMu_Skim::Finish(){
	unsigned int dytautau = 6;
	//if(Nminus0.at(0).at(dytautau).Integral()!=0)if(HConfig.GetHisto(false,DataMCType::DY_embedded,dytautau))ScaleAllHistOfType(dytautau,19789.302*1121.0/869489./Nminus0.at(0).at(dytautau).Integral());
	Selection::Finish();
}
