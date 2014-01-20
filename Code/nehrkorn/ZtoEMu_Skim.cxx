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
  ,mu_ptlow(10)
  ,mu_pthigh(20)
  ,e_ptlow(10)
  ,e_pthigh(20)
  ,mu_eta(2.4)
  ,e_eta(2.5)
  ,jet_pt(18)
  ,jet_eta(4.)
  ,jet_sum(70)
  ,zmin(88)
  ,zmax(94)
{
    //verbose=true;
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
    if(i==NE)                 cut.at(NE)=1;
    if(i==ptthreshold)        cut.at(ptthreshold)=1;
    if(i==charge)             cut.at(charge)=0;
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
      title.at(i)="Number $\\mu >=$";
      title.at(i)+=cut.at(NMu);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NE){
      title.at(i)="Number $e >=$";
      title.at(i)+=cut.at(NE);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==ptthreshold){
        title.at(i)="ptthreshold ";
        hlabel="ptthreshold ";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ptthreshold_",htitle,17,-0.5,16.5,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ptthreshold_",htitle,17,-0.5,16.5,hlabel,"Events"));
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
    
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NPV=HConfig.GetTH1D(Name+"_NPV","NPV",60,0.,60.,"Number of vertices");
  mupt=HConfig.GetTH1D(Name+"_mupt","mupt",40,0.,100.,"p_{T}^{#mu}");
  mueta=HConfig.GetTH1D(Name+"_mueta","mueta",20,-2.5,2.5,"#eta_{#mu}");
  ept=HConfig.GetTH1D(Name+"_ept","ept",40,0.,100.,"p_{T}^{e}");
  eeta=HConfig.GetTH1D(Name+"_eeta","eeta",20,-2.5,2.5,"#eta_{e}");
  mupt_n0=HConfig.GetTH1D(Name+"_mupt_n0","mupt_n0",40,0.,100.,"p_{T}^{#mu}");
  mueta_n0=HConfig.GetTH1D(Name+"_mueta_n0","mueta_n0",20,-2.5,2.5,"#eta_{#mu}");
  ept_n0=HConfig.GetTH1D(Name+"_ept_n0","ept_n0",40,0.,100.,"p_{T}^{e}");
  eeta_n0=HConfig.GetTH1D(Name+"_eeta_n0","eeta_n0",20,-2.5,2.5,"#eta_{e}");
  muptw=HConfig.GetTH1D(Name+"_muptw","muptw",40,0.,100.,"p_{T}^{#mu}");
  muetaw=HConfig.GetTH1D(Name+"_muetaw","muetaw",20,-2.5,2.5,"#eta_{#mu}");
  eptw=HConfig.GetTH1D(Name+"_eptw","eptw",40,0.,100.,"p_{T}^{e}");
  eetaw=HConfig.GetTH1D(Name+"_eetaw","eetaw",20,-2.5,2.5,"#eta_{e}");
  goodmuons=HConfig.GetTH1D(Name+"_goodmuons","goodmuons",20,0.,20.,"Number of tight muons");
  fakemuons=HConfig.GetTH1D(Name+"_fakemuons","fakemuons",20,0.,20.,"Number of fake muons");
  goodelectrons=HConfig.GetTH1D(Name+"_goodelectrons","goodelectrons",20,0.,20.,"Number of tight electrons");
  fakeelectrons=HConfig.GetTH1D(Name+"_fakeelectrons","fakeelectrons",20,0.,20.,"Number of fake electrons");
  nontriggr20=HConfig.GetTH1D(Name+"_nontriggr20","nontriggr20",100,-1.,1.,"nontrig MVA discriminator (p_{t}#geq20 GeV)");
  nontrigsm20=HConfig.GetTH1D(Name+"_nontrigsm20","nontrigsm20",100,-1.,1.,"nontrig MVA discriminator (p_{t}<20 GeV)");
  triggr20=HConfig.GetTH1D(Name+"_triggr20","triggr20",100,-1.,1.,"trig MVA discriminator (p_{t}#geq20 GeV)");
  trigsm20=HConfig.GetTH1D(Name+"_trigsm20","trigsm20",100,-1.,1.,"trig MVA discriminator (p_{t}<20 GeV)");
  trignoipgr20=HConfig.GetTH1D(Name+"_trignoipgr20","trignoipgr20",100,-1.,1.,"trignoip MVA discriminator (p_{t}#geq20 GeV)");
  trignoipsm20=HConfig.GetTH1D(Name+"_trignoipsm20","trignoipsm20",100,-1.,1.,"trignoip MVA discriminator (p_{t}<20 GeV)");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_Skim::Store_ExtraDist(){
	Extradist1d.push_back(&NPV);
	Extradist1d.push_back(&mupt);
	Extradist1d.push_back(&mueta);
	Extradist1d.push_back(&ept);
	Extradist1d.push_back(&eeta);
	Extradist1d.push_back(&muptw);
	Extradist1d.push_back(&muetaw);
	Extradist1d.push_back(&eptw);
	Extradist1d.push_back(&eetaw);
	Extradist1d.push_back(&mupt_n0);
	Extradist1d.push_back(&mueta_n0);
	Extradist1d.push_back(&ept_n0);
	Extradist1d.push_back(&eeta_n0);
	Extradist1d.push_back(&goodmuons);
	Extradist1d.push_back(&fakemuons);
	Extradist1d.push_back(&goodelectrons);
	Extradist1d.push_back(&fakeelectrons);
	Extradist1d.push_back(&nontriggr20);
	Extradist1d.push_back(&nontrigsm20);
	Extradist1d.push_back(&triggr20);
	Extradist1d.push_back(&trigsm20);
	Extradist1d.push_back(&trignoipgr20);
	Extradist1d.push_back(&trignoipsm20);
}

void  ZtoEMu_Skim::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}

  value.at(TriggerOk)=0;
  /*std::cout << "Number of HLT Triggers:" << Ntp->NHLTTriggers() << std::endl;
  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
	  std::cout << "######################################" << std::endl;
	  std::cout << "Trigger: " << Ntp->HTLTriggerName(i) << std::endl;
	  std::cout << "Acceptet (1) or not (0): " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << std::endl;
	  std::cout << "######################################" << std::endl;
  }*/
  if(Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") || Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")){
	  value.at(TriggerOk)=1;
  }
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));

  // Apply Selection
    
  ///////////////////////////////////////////////
  //
  // Vertex selection
  //
  if(verbose)std::cout << "Vertex selection" << std::endl;
  unsigned int nGoodVtx=0;
  int vertex = -1;
  for(unsigned i=0;i<Ntp->NVtx();i++){
	  if(isGoodVtx(i)){
		  if(vertex==-1)vertex=i;
		  nGoodVtx++;
	  }
  }
  if(verbose)std::cout << "selected vertex: " << vertex << std::endl;
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));

  ///////////////////////////////////////////////
  //
  // Mu Cuts
  //
  if(verbose) std::cout << "Muon cuts" << std::endl;
  std::vector<unsigned int> GoodMuons;
  std::vector<unsigned int> Fakemuons;
  std::vector<double> FakeRateMuVec;
  
  // muon ID cuts
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muons_p4(i).Pt()>mu_ptlow
			  && fabs(Ntp->Muons_p4(i).Eta())<mu_eta
			  && (matchTrigger(i,0.5,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon") || matchTrigger(i,0.5,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon"))
			  && vertex>=0
			  ){
		  if(isTightMuon(i,vertex)
				  && Muon_RelIso(i)<0.12
				  ){
			  if(Ntp->GetMCID()==32 && !fabs(matchTruth(Ntp->Muons_p4(i))==13)) continue;
			  GoodMuons.push_back(i);
		  }else if(isFakeMuon(i,vertex)
				  && Ntp->isData()
				  ){
			  Fakemuons.push_back(i);
			  //GoodMuons.push_back(i);
			  FakeRateMuVec.push_back(Fakerate(Ntp->Muons_p4(i),MuonFakeRate,"muon"));
		  }
	  }
  }
  
  value.at(NMu)=GoodMuons.size();
  if(GoodMuons.size()==0)value.at(NMu)=Fakemuons.size();
  pass.at(NMu)=(value.at(NMu)>=cut.at(NMu));
  
  unsigned int muidx(999);
  double hardestmu(0);
  if(GoodMuons.size()>1){
	  for(unsigned i=0;i<GoodMuons.size();i++){
		  if(Ntp->Muons_p4(GoodMuons.at(i)).Pt()>hardestmu){
			  hardestmu = Ntp->Muons_p4(GoodMuons.at(i)).Pt();
			  muidx = GoodMuons.at(i);
		  }
	  }
  }
  if(GoodMuons.size()==1){muidx=GoodMuons.at(0);}
  if(GoodMuons.size()==0){
	  if(Fakemuons.size()>1){
		  for(unsigned i=0;i<Fakemuons.size();i++){
			  if(Ntp->Muons_p4(Fakemuons.at(i)).Pt()>hardestmu){
				  hardestmu = Ntp->Muons_p4(Fakemuons.at(i)).Pt();
				  muidx = Fakemuons.at(i);
			  }
		  }
	  }
	  if(Fakemuons.size()==1){muidx=Fakemuons.at(0);}
  }

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  std::vector<unsigned int> Fakeelectrons;
  std::vector<double> FakeRateEVec;
  bool hasMuonTrack = false;
  bool notThisOne = false;
  bool matchRecoMuon = false;
  
  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>e_ptlow
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && (matchTrigger(i,0.5,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron") || matchTrigger(i,0.5,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron"))
			  && vertex>=0
			  ){
		  for(unsigned j=0;j<GoodMuons.size();j++){
			  if(Ntp->Electron_Track_idx(i)==Ntp->Muon_Track_idx(GoodMuons.at(j))) hasMuonTrack = true;
		  }
		  if(GoodMuons.size()==0){
			  for(unsigned j=0;j<Fakemuons.size();j++){
				  if(Ntp->Electron_Track_idx(i)==Ntp->Muon_Track_idx(Fakemuons.at(j))) hasMuonTrack = true;
			  }
		  }
		  if(hasMuonTrack) continue;
		  for(unsigned j=0;j<Ntp->NMuons();j++){
			  for(unsigned k=0;k<GoodMuons.size();k++){
				  if(j==k){
					  notThisOne = true;
				  }
				  if(notThisOne)break;
			  }
			  if(!notThisOne && Ntp->Electron_p4(i).DeltaR(Ntp->Muons_p4(j))<0.3) matchRecoMuon = true;
		  }
		  if(matchRecoMuon) continue;
		  if(isMVANonTrigElectron(i,vertex)){
			  if(Ntp->GetMCID()==31 && !(fabs(matchTruth(Ntp->Electron_p4(i)))==11 || fabs(matchTruth(Ntp->Electron_p4(i)))==13 || fabs(matchTruth(Ntp->Electron_p4(i)))==22)) continue;
			  GoodElectrons.push_back(i);
		  }else if(isFakeElectron(i,vertex)
				  && Ntp->isData()
				  ){
			  Fakeelectrons.push_back(i);
			  //GoodElectrons.push_back(i);
			  FakeRateEVec.push_back(Fakerate(Ntp->Electron_p4(i),ElectronFakeRate,"electron"));
		  }
	  }
  }
  
  value.at(NE)=GoodElectrons.size();
  if(GoodElectrons.size()==0)value.at(NE)=Fakeelectrons.size();
  pass.at(NE)=(value.at(NE)>=cut.at(NE));
  
  unsigned int eidx(999);
  double hardeste(0);
  if(GoodElectrons.size()>1){
	  for(unsigned i=0;i<GoodElectrons.size();i++){
		  if(Ntp->Electron_p4(GoodElectrons.at(i)).Et()>hardeste){
			  hardeste = Ntp->Electron_p4(GoodElectrons.at(i)).Et();
			  eidx = GoodElectrons.at(i);
		  }
	  }
  }
  if(GoodElectrons.size()==1){eidx=GoodElectrons.at(0);}
  if(GoodElectrons.size()==0){
	  if(Fakeelectrons.size()>1){
		  for(unsigned i=0;i<Fakeelectrons.size();i++){
			  if(Ntp->Electron_p4(Fakeelectrons.at(i)).Et()>hardeste){
				  hardeste = Ntp->Electron_p4(Fakeelectrons.at(i)).Et();
				  eidx = Fakeelectrons.at(i);
			  }
		  }
	  }
	  if(Fakeelectrons.size()==1){eidx=Fakeelectrons.at(0);}
  }
  
  ///////////////////////////////////////////////
  //
  // pt thresholds
  //
  if(verbose) std::cout << "Setting pt thresholds" << std::endl;
  value.at(ptthreshold)=0;
  if(muidx!=999 && eidx!=999){
	  value.at(ptthreshold)=1;
	  if(Ntp->Muons_p4(muidx).Pt()<=mu_ptlow || Ntp->Electron_p4(eidx).Et()<=e_ptlow) value.at(ptthreshold)=0;
	  if(Ntp->Muons_p4(muidx).Pt()<mu_pthigh && Ntp->Electron_p4(eidx).Et()<e_pthigh) value.at(ptthreshold)=0;
	  if(Ntp->Muons_p4(muidx).Pt()<mu_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
	  }else if(Ntp->Electron_p4(eidx).Et()<e_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
	  }
  }
  pass.at(ptthreshold)=(value.at(ptthreshold)==cut.at(ptthreshold));

  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=-5;
  if(eidx!=999 && muidx!=999){
	  value.at(charge)=Ntp->Electron_Charge(eidx)+Ntp->Muon_Charge(muidx);
  }
  pass.at(charge)=(value.at(charge)==cut.at(charge));
  
  ////////////////////////////////////////////////
  //
  // QCD
  //
  bool fakemu = false;
  bool fakee = false;
  for(unsigned i=0;i<Fakemuons.size();i++){
	  if(Fakemuons.at(i)==muidx){
		  fakemu=true;
		  fakeRateMu = FakeRateMuVec.at(i);
	  }
  }
  for(unsigned i=0;i<Fakeelectrons.size();i++){
	  if(Fakeelectrons.at(i)==eidx){
		  fakee=true;
		  fakeRateE = FakeRateEVec.at(i);
	  }
  }
  fakeRate = 1.;
  if(pass.at(charge)
		  && Ntp->isData()
		  ){
	  if(fakemu || fakee) fakeRate = 0.;
  }
  if(!pass.at(charge)
		  && Ntp->isData()
		  ){
	  if(fakemu && !fakee){
		  fakeRate = fakeRateMu;
		  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
		  pass.at(charge)=true;
	  }else if(fakee && !fakemu){
		  fakeRate = fakeRateE;
		  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
		  pass.at(charge)=true;
	  }else if(fakemu && fakee){
		  fakeRate = fakeRateMu*fakeRateE;
		  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
		  pass.at(charge) = true;
	  }
  }
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1),ww(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=34){
    w*=Ntp->EvtWeight3D();
    ww*=Ntp->EvtWeight3D();
    if(verbose)std::cout << "void  ZtoEMu_Skim::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1*fakeRate;wobs=1;ww=1*fakeRate;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  //NPV.at(t).Fill(Ntp->NVtx(),w);
  if(muidx!=999){
	  mupt.at(t).Fill(Ntp->Muons_p4(muidx).Pt(),w);
	  mueta.at(t).Fill(Ntp->Muons_p4(muidx).Eta(),w);
	  muptw.at(t).Fill(Ntp->Muons_p4(muidx).Pt(),ww);
	  muetaw.at(t).Fill(Ntp->Muons_p4(muidx).Eta(),ww);
  }
  if(eidx!=999){
	  ept.at(t).Fill(Ntp->Electron_p4(eidx).Pt(),w);
	  eeta.at(t).Fill(Ntp->Electron_supercluster_eta(eidx),w);
	  eptw.at(t).Fill(Ntp->Electron_p4(eidx).Pt(),ww);
	  eetaw.at(t).Fill(Ntp->Electron_supercluster_eta(eidx),ww);
  }
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(charge)
		  ){
	  NPV.at(t).Fill(Ntp->NVtx(),w);
	  mupt_n0.at(t).Fill(Ntp->Muons_p4(muidx).Pt(),w);
	  mueta_n0.at(t).Fill(Ntp->Muons_p4(muidx).Eta(),w);
	  ept_n0.at(t).Fill(Ntp->Electron_p4(eidx).Pt(),w);
	  eeta_n0.at(t).Fill(Ntp->Electron_supercluster_eta(eidx),w);
  }
  goodmuons.at(t).Fill(GoodMuons.size());
  fakemuons.at(t).Fill(Fakemuons.size());
  goodelectrons.at(t).Fill(GoodElectrons.size());
  fakeelectrons.at(t).Fill(Fakeelectrons.size());
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(!Ntp->Electron_HasMatchedConversions(i) &&
			  Ntp->Electron_numberOfMissedHits(i)<=1){
		  if(Ntp->Electron_p4(i).Pt()>10 && Ntp->Electron_p4(i).Pt()<20){
			  nontrigsm20.at(t).Fill(Ntp->Electron_MVA_NonTrig_discriminator(i));
			  trigsm20.at(t).Fill(Ntp->Electron_MVA_Trig_discriminator(i));
			  trignoipsm20.at(t).Fill(Ntp->Electron_MVA_TrigNoIP_discriminator(i));
		  }else if(Ntp->Electron_p4(i).Pt()>=20){
			  nontriggr20.at(t).Fill(Ntp->Electron_MVA_NonTrig_discriminator(i));
			  triggr20.at(t).Fill(Ntp->Electron_MVA_Trig_discriminator(i));
			  trignoipgr20.at(t).Fill(Ntp->Electron_MVA_TrigNoIP_discriminator(i));
		  }
	  }
  }
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_Skim::isGoodVtx(unsigned int i){
	if(fabs(Ntp->Vtx(i).z())>=24) return false;
	if(Ntp->Vtx(i).Perp()>=2) return false;
	if(Ntp->Vtx_ndof(i)<=4) return false;
	if(Ntp->Vtx_isFake(i)!=0) return false;
	return true;
}

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

double ZtoEMu_Skim::vertexSignificance(TVector3 vec, unsigned int vertex){
	if(vertex>=0 && vertex<Ntp->NVtx()){
		const float elm[3] = {(vec.X()-Ntp->Vtx(vertex).X()),(vec.Y()-Ntp->Vtx(vertex).Y()),(vec.Z()-Ntp->Vtx(vertex).Z())};
		TVectorF diff(3,elm);
		TMatrixF M(Ntp->Vtx_Cov(vertex));
		if(M.IsValid()){
			double mag = diff.Norm2Sqr();
			double sim = M.Similarity(diff);
			return mag/sqrt(sim);
		}
	}
	return 999;
}

bool ZtoEMu_Skim::matchTrigger(unsigned int i, double dr, std::string trigger, std::string object){
	unsigned int id = 0;
	TLorentzVector particle(0.,0.,0.,0.);
	TLorentzVector triggerObj(0.,0.,0.,0.);
	if(object=="electron"){
		id = 82;
		particle = Ntp->Electron_p4(i);
	}
	if(object=="muon"){
		id = 83;
		particle = Ntp->Muons_p4(i);
	}
	for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
		if(Ntp->HLTTrigger_objs_trigger(i).find(trigger) != string::npos){
			for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
				if(Ntp->HLTTrigger_objs_Id(i,j)==id){
					triggerObj.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
							Ntp->HLTTrigger_objs_Eta(i,j),
							Ntp->HLTTrigger_objs_Phi(i,j),
							Ntp->HLTTrigger_objs_E(i,j));
				}
				if(triggerObj.Pt()>0.
						&& particle.Pt()>0.
						&& particle.DeltaR(triggerObj)<dr) return true;
			}
		}
	}
	return false;
}

int ZtoEMu_Skim::matchTruth(TLorentzVector tvector){
	double testdr=1.;
	int pdgid = 0;
	for(unsigned i=0;i<Ntp->NMCParticles();i++){
		if(Ntp->MCParticle_p4(i).Pt()>0.){
			if(tvector.DeltaR(Ntp->MCParticle_p4(i))<testdr){
				testdr = tvector.DeltaR(Ntp->MCParticle_p4(i));
				pdgid = Ntp->MCParticle_pdgid(i);
			}
		}
	}
	return pdgid;
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_Skim::isTightMuon(unsigned int i){
	if(!Ntp->Muon_isGlobalMuon(i)) return false;
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(Ntp->Muon_normChi2(i)>=10.) return false;
	if(Ntp->Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
	if(Ntp->Muon_numberOfMatchedStations(i)<=1) return false;
	if(Ntp->Muon_numberofValidPixelHits(i)<=0) return false;
	if(Ntp->Muon_trackerLayersWithMeasurement(i)<=5) return false;
	return true;
}

bool ZtoEMu_Skim::isTightMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.5) return false;
	return true;
}

bool ZtoEMu_Skim::isLooseMuon(unsigned int i){
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(!(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))) return false;
	return true;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int i){
	if(!Ntp->Muon_isGlobalMuon(i)) return false;
	if(Ntp->Muons_p4(i).Pt()<=10) return false;
	if(fabs(Ntp->Muons_p4(i).Eta())>2.4) return false;
	if(Ntp->Muons_p4(i).Pt()<=20){
		if(Ntp->Muon_sumPt03(i)>=8.) return false;
		if(Ntp->Muon_emEt03(i)>=8.) return false;
		if(Ntp->Muon_hadEt03(i)>=8.) return false;
	}
	if(Ntp->Muons_p4(i).Pt()>20){
		if(Ntp->Muon_sumPt03(i)/Ntp->Muons_p4(i).Pt()>=0.4) return false;
		if(Ntp->Muon_emEt03(i)/Ntp->Muons_p4(i).Pt()>=0.4) return false;
		if(Ntp->Muon_hadEt03(i)/Ntp->Muons_p4(i).Pt()>=0.4) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeMuon(i)) return false;
	if(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
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

bool ZtoEMu_Skim::isMVATrigNoIPElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(mvapt<20){
		if(mvaeta<0.8){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=-0.5375) return false;
		}
		if(mvaeta>=0.8 && mvaeta<1.479){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=-0.375) return false;
		}
		if(mvaeta>=1.479 && mvaeta<2.5){
			if(Electron_RelIso(i)>=0.10) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=-0.025) return false;
		}
	}
	if(mvapt>=20){
		if(mvaeta<0.8){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=0.325) return false;
		}
		if(mvaeta>=0.8 && mvaeta<1.479){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
		}
		if(mvaeta>=1.479 && mvaeta<2.5){
			if(Electron_RelIso(i)>=0.10) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
		}
	}
	return true;
}

bool ZtoEMu_Skim::isMVANonTrigElectron(unsigned int i, unsigned int j){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)>1) return false;
	if(vertexSignificance(Ntp->Electron_Poca(i),j)>=4) return false;
	if(Electron_RelIso(i)>=0.4) return false;
	if(mvapt>7. && mvapt<10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.47) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.004) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.295) return false;
	}
	if(mvapt>=10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=-0.34) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=-0.65) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.6) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isMVATrigElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Electron_RelIso(i)>=0.15) return false;
	if(mvapt>10. && mvapt<20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.00) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.10) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.62) return false;
	}else if(mvapt>=20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.94) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.85) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.92) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isTightElectron(unsigned int i){
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_RelIso(i)>=0.1) return false;
	if(fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))>=0.05) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.004) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.03) return false;
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.01) return false;
		if(Ntp->Electron_hadronicOverEm(i)>=0.12) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(i))>1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.005) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.02) return false;
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.03) return false;
		if(Ntp->Electron_hadronicOverEm(i)>=0.10) return false;
		if(Ntp->Electron_p4(i).Pt()<20 && Electron_RelIso(i)>=0.07) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isTightElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightElectron(i)) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int i){
	if(Ntp->Electron_p4(i).Pt()<=10) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>1) return false;
	if(Ntp->Electron_tkSumPt03(i)/Ntp->Electron_p4(i).Pt()>=0.2) return false;
	if(std::max(Ntp->Electron_ecalRecHitSumEt03(i)-1.,0.)/Ntp->Electron_p4(i).Pt()>=0.2) return false;
	if((Ntp->Electron_hcalDepth1TowerSumEt03(i)+Ntp->Electron_hcalDepth2TowerSumEt03(i))/Ntp->Electron_p4(i).Pt()>=0.2) return false;
	//if(Electron_RelIso(i)>=0.2) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.01) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.15) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.007) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.03) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.10) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.009) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(i))>=2.5) return false;
	return true;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeElectron(i)) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	return true;
}

double ZtoEMu_Skim::Electron_RelIso(unsigned int i){
	return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff_R04(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();
}

double ZtoEMu_Skim::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.3) return 0.194;
	else if(eta>=2.4) return 0.261;
}

double ZtoEMu_Skim::Electron_Aeff_R03(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.13;
	else if(eta>=1. && eta<1.479) return 0.14;
	else if(eta>=1.479 && eta<2.) return 0.07;
	else if(eta>=2. && eta<2.2) return 0.09;
	else if(eta>=2.2 && eta<2.3) return 0.11;
	else if(eta>=2.3 && eta<2.3) return 0.11;
	else if(eta>=2.4) return 0.14;
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
	
	return fakerate/(1-fakerate);
}

//////////////////////////////
//
// Finish function
//

void ZtoEMu_Skim::Finish(){
	Selection::Finish();
}
