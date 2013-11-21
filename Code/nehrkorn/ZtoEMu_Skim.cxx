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
  ,mu_eta(2.4)
  ,e_eta(2.5)
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
  muptw=HConfig.GetTH1D(Name+"_muptw","muptw",40,0.,100.,"p_{T}^{#mu}");
  muetaw=HConfig.GetTH1D(Name+"_muetaw","muetaw",20,-2.5,2.5,"#eta_{#mu}");
  eptw=HConfig.GetTH1D(Name+"_eptw","eptw",40,0.,100.,"p_{T}^{e}");
  eetaw=HConfig.GetTH1D(Name+"_eetaw","eetaw",20,-2.5,2.5,"#eta_{e}");
  goodmuons=HConfig.GetTH1D(Name+"_goodmuons","goodmuons",20,0.,20.,"Number of tight muons");
  fakemuons=HConfig.GetTH1D(Name+"_fakemuons","fakemuons",20,0.,20.,"Number of fake muons");
  goodelectrons=HConfig.GetTH1D(Name+"_goodelectrons","goodelectrons",20,0.,20.,"Number of tight electrons");
  fakeelectrons=HConfig.GetTH1D(Name+"_fakeelectrons","fakeelectrons",20,0.,20.,"Number of fake electrons");
  discrgr20=HConfig.GetTH1D(Name+"_discrgr20","discrgr20",100,-1.,1.,"MVA discriminator (p_{t}#geq20 GeV)");
  discrsm20=HConfig.GetTH1D(Name+"_discrsm20","discrsm20",100,-1.,1.,"MVA discriminator (p_{t}<20 GeV)");
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
	Extradist1d.push_back(&goodmuons);
	Extradist1d.push_back(&fakemuons);
	Extradist1d.push_back(&goodelectrons);
	Extradist1d.push_back(&fakeelectrons);
	Extradist1d.push_back(&discrgr20);
	Extradist1d.push_back(&discrsm20);
}

void  ZtoEMu_Skim::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
  
	/*double embedd_weight = 1.;
	if(id==34){
		embedd_weight = Ntp->EmbeddedWeight();
		if(Ntp->EmbeddedWeight()!=Ntp->TauSpinnerWeight()*Ntp->SelEffWeight()*Ntp->RadiationCorrWeight()*Ntp->MinVisPtFilter()*Ntp->KinWeightPt()*Ntp->KinWeightEta()*Ntp->KinWeightMassPt()){
			std::cout << "Product = " << Ntp->TauSpinnerWeight()*Ntp->SelEffWeight()*Ntp->RadiationCorrWeight()*Ntp->MinVisPtFilter()*Ntp->KinWeightPt()*Ntp->KinWeightEta()*Ntp->KinWeightMassPt() << std::endl;
			std::cout << "Embedded weight = " << embedd_weight << std::endl;
		}
	}*/

  value.at(TriggerOk)=0;
  /*std::cout << "Number of HLT Triggers:" << Ntp->NHLTTriggers() << std::endl;
  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
	  std::cout << "######################################" << std::endl;
	  std::cout << "Trigger: " << Ntp->HTLTriggerName(i) << std::endl;
	  std::cout << "Acceptet (1) or not (0): " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << std::endl;
	  std::cout << "######################################" << std::endl;
  }*/
  //if(Ntp->TriggerAccept("HLT_IsoMu24_eta2p1"))value.at(TriggerOk)=1;
  if(Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloId")
		  && !Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloId")
		  ){
	  value.at(TriggerOk)=1;
	  mu_pt = 18;
	  e_pt = 9;
  }
  if(Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloId")
		  && !Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloId")
		  ){
	  value.at(TriggerOk)=1;
	  mu_pt = 9;
	  e_pt = 18;
  }
  if(Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloId")
		  && Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloId")
		  ){
	  value.at(TriggerOk)=1;
	  mu_pt = 18;
	  e_pt = 18;
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
  std::vector<int> Fakemuons;
  Fakemuons.clear();
  
  // muon ID cuts (including eta-dependent isolation)
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muons_p4(i).Pt()>=mu_pt &&
			  isTightMuon(i,vertex) &&
			  dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.02 &&
			  dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.5
			  ){
		  if(fabs(Ntp->Muons_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15){
			  GoodMuons.push_back(i);
		  }else if(fabs(Ntp->Muons_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.1){
			  GoodMuons.push_back(i);
		  }
	  }else if(isFakeMuon(i,vertex)
			  && Ntp->isData()
			  && MVA_ID
			  ){
		  Fakemuons.push_back(i);
		  GoodMuons.push_back(i);
		  fakeRateMu = Fakerate(Ntp->Muons_p4(i),MuonFakeRate,"muon");
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

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  std::vector<unsigned int> Fakeelectrons;
  Fakeelectrons.clear();
  
  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>=e_pt &&
			  muidx!=999 &&
			  vertex>=0 &&
			  Ntp->Electron_p4(i).DeltaR(Ntp->Muons_p4(muidx))>0.2 &&
			  !Ntp->Electron_HasMatchedConversions(i) &&
			  Ntp->Electron_numberOfMissedHits(i)==0 &&
			  dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.2 &&
			  dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.02
			  ){
		  if(MVA_ID){
			  if(isMVATrigNoIPElectron(i)){
				  GoodElectrons.push_back(i);
			  }else if(isFakeElectron(i,vertex) && Ntp->isData()){
				  Fakeelectrons.push_back(i);
				  GoodElectrons.push_back(i);
				  fakeRateE = Fakerate(Ntp->Electron_p4(i),ElectronFakeRate,"electron");
			  }
		  }else{
		  	if(isTightElectron(i,vertex))GoodElectrons.push_back(i);
		  }
	  }
  }
  
  // electron pt cut
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEPt).push_back(Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
    if(Ntp->Electron_p4(GoodElectrons.at(i)).Et()<e_pt){
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
  pass.at(NE)=(value.at(NE)>=cut.at(NE));
  
  unsigned int eidx(999);
  double hardeste(0);
  if(GoodElectrons.size()>1){
	  for(unsigned i=0;i<GoodElectrons.size();i++){
		  if(Ntp->Electron_p4(GoodElectrons.at(i)).Pt()>hardeste){
			  hardeste = Ntp->Electron_p4(GoodElectrons.at(i)).Pt();
			  eidx = GoodElectrons.at(i);
		  }
	  }
  }
  if(GoodElectrons.size()==1){eidx=GoodElectrons.at(0);}
  
  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=-5;
  if(eidx!=999 && muidx!=999){
	  value.at(charge)=Ntp->Electron_Charge(GoodElectrons.at(0))+Ntp->Muon_Charge(GoodMuons.at(0));
  }
  pass.at(charge)=(value.at(charge)==cut.at(charge));
  
  ////////////////////////////////////////////////
  //
  // QCD
  //
  bool fakemu = false;
  bool fakee = false;
  for(unsigned i=0;i<Fakemuons.size();i++){
	  if(Fakemuons.at(i)==muidx)fakemu=true;
  }
  for(unsigned i=0;i<Fakeelectrons.size();i++){
	  if(Fakeelectrons.at(i)==eidx)fakee=true;
  }
  if(MVA_ID){
	  fakeRate = 0.;
	  if(!pass.at(charge)
		  && Ntp->isData()
		  ){
		  if(fakemu && !fakee){
			  fakeRate = fakeRateMu;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
		  }else if(fakee && !fakemu){
			  fakeRate = fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
		  }else if(fakemu && fakee){
			  fakeRate = fakeRateMu*fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
		  }else{
			  fakeRate = 1.;
		  }
	  }else if(pass.at(charge)
		  && Ntp->isData()
		  && !fakemu
		  && !fakee
		  ){
		  fakeRate = 1.;
	  }
  }else{
	fakeRate = 1.;
  }
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1),ww(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=34){
    w*=Ntp->EvtWeight3D();
    ww*=Ntp->EvtWeight3D();
    if(pass.at(NE)){
		ww*=ElectronIDeff(eidx)*ElectronTriggerEff(eidx);
	}
	if(pass.at(NMu)){
		ww*=MuonIDeff(muidx)*MuonTriggerEff(muidx);
	}
    if(verbose)std::cout << "void  ZtoEMu_Skim::doEvent() k" << w << " " << wobs << std::endl;
  }/*else if(Ntp->GetMCID()==34){
		w*=embedd_weight;
		if(pass.at(NE))w*=ElectronEffRecHit(GoodElectrons.at(0));
		if(pass.at(NMu))w*=MuonDataSF(GoodMuons.at(0));
  }*/
  else{w=1*fakeRate;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  NPV.at(t).Fill(Ntp->NVtx(),w);
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
  goodmuons.at(t).Fill(GoodMuons.size());
  fakemuons.at(t).Fill(Fakemuons.size());
  goodelectrons.at(t).Fill(GoodElectrons.size());
  fakeelectrons.at(t).Fill(Fakeelectrons.size());
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(vertex>=0 &&
			  !Ntp->Electron_HasMatchedConversions(i) &&
			  Ntp->Electron_numberOfMissedHits(i)==0 &&
			  dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.2 &&
			  dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.02){
		  if(Ntp->Electron_p4(i).Pt()>10 && Ntp->Electron_p4(i).Pt()<20){
			  discrsm20.at(t).Fill(Ntp->Electron_MVA_discriminator(i));
		  }else if(Ntp->Electron_p4(i).Pt()>=20){
			  discrgr20.at(t).Fill(Ntp->Electron_MVA_discriminator(i));
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
	if(fabs(Ntp->Vtx(i).z())<24
			&& Ntp->Vtx(i).Perp()<2
			&& Ntp->Vtx_ndof(i)>4
			&& Ntp->Vtx_isFake(i)==0
			){
		return true;
	}
	return false;
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

double ZtoEMu_Fakerate::vertexSignificance(TVector3 vec, unsigned int vertex){
	TVectorF diff;
	if(vertex<Ntp->NVtx()){
		float elements[3] = {(vec.X()-Ntp->Vtx(vertex).X()),(vec.Y()-Ntp->Vtx(vertex).Y()),(vec.Z()-Ntp->Vtx(vertex).Z())};
		diff = TVectorF(3,elements);
		TVectorF v2 = diff;
		v2*=Ntp->Vtx_Cov(vertex);
		return diff.Norm2Sqr()/sqrt(diff*v2);
	}
	return 999;
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
	if(j<0 || j>=Ntp->NVtx()){
		return false;
	}
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
		if(Ntp->Muons_p4(i).Pt()>=20 &&
				Muon_RelIso(i)<0.4
				){
			return true;
		}else if(Ntp->Muons_p4(i).Pt()<20 &&
				Muon_AbsIso(i)<8.
				){
			return true;
		}
	}
	return false;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()){
		return false;
	}
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

bool ZtoEMu_Fakerate::isMVATrigNoIPElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(mvapt<20){
		if(mvaeta<=0.8 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_TrigNoIP_discriminator(i)>-0.5375) return true;
		else if(mvaeta>0.8 && mvaeta<=1.479 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_TrigNoIP_discriminator(i)>-0.375) return true;
		else if(mvaeta>1.479 && Electron_RelIso(i)<0.10 && Ntp->Electron_MVA_TrigNoIP_discriminator(i)>-0.025) return true;
	}else if(mvapt>=20){
		if(mvaeta<=0.8 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_TrigNoIP_discriminator(i)>0.325) return true;
		else if(mvaeta>0.8 && mvaeta<=1.479 && Electron_RelIso(i)<0.15 && Ntp->Electron_MVA_TrigNoIP_discriminator(i)>0.775) return true;
		else if(mvaeta>1.479 && Electron_RelIso(i)<0.10 && Ntp->Electron_MVA_TrigNoIP_discriminator(i)>0.775) return true;
	}
	return false;
}

bool ZtoEMu_Fakerate::isMVANonTrigElectron(unsigned int i, unsigned int j){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)<=1
			&& vertexSignificance(Ntp->Electron_Poca(i),j)<4
			&& Electron_RelIso(i)<0.4){
		if(mvapt>7. && mvapt<10.){
			if(mvaeta<=0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)>0.47) return true;
			else if(mvaeta>0.8 && mvaeta<=1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)>0.004) return true;
			else if(mvaeta>1.479 && mvaeta<=2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)>0.295) return true;
		}else if(mvapt>=10.){
			if(mvaeta<=0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)>-0.34) return true;
			else if(mvaeta>0.8 && mvaeta<=1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)>-0.65) return true;
			else if(mvaeta>1.479 && mvaeta<=2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)>0.6) return true;
		}
	}
	return false;
}

bool ZtoEMu_Fakerate::isMVATrigElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)==0
			&& !Ntp->Electron_HasMatchedConversions(i)
			&& Electron_RelIso(i)<0.15){
		if(mvapt>10. && mvapt<20.){
			if(mvaeta<=0.8 && Ntp->Electron_MVA_Trig_discriminator(i)>0.00) return true;
			else if(mvaeta>0.8 && mvaeta<=1.479 && Ntp->Electron_MVA_Trig_discriminator(i)>0.10) return true;
			else if(mvaeta>1.479 && mvaeta<=2.5 && Ntp->Electron_MVA_Trig_discriminator(i)>0.62) return true;
		}else if(mvapt>=20.){
			if(mvaeta<=0.8 && Ntp->Electron_MVA_Trig_discriminator(i)>0.94) return true;
			else if(mvaeta>0.8 && mvaeta<=1.479 && Ntp->Electron_MVA_Trig_discriminator(i)>0.85) return true;
			else if(mvaeta>1.479 && mvaeta<=2.5 && Ntp->Electron_MVA_Trig_discriminator(i)>0.92) return true;
		}
	}
	return false;
}

bool ZtoEMu_Fakerate::isTightElectron(unsigned int i){
	if(verbose)std::cout << "isTightElectron(unsigned int i)" << std::endl;
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){ //barrel
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.004 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.03 &&
				Ntp->Electron_sigmaIetaIeta(i)<0.01 &&
				Ntp->Electron_hadronicOverEm(i)<0.12 &&
				fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))<0.05 &&
				Electron_RelIso(i)<0.10 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_numberOfMissedHits(i)==0
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(i))>1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){ //endcaps
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

bool ZtoEMu_Fakerate::isTightElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()){
		return false;
	}
	if(isTightElectron(i)
			&& dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.02
			&& dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.1
			){
		return true;
	}
	return false;
}

bool ZtoEMu_Fakerate::isFakeElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){
		if(Ntp->Electron_p4(i).Pt()>10 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_sigmaIetaIeta(i)<0.01 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.15 &&
				Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.007 &&
				Electron_RelIso(i)<0.2
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(i))>1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_p4(i).Pt()>10 &&
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

bool ZtoEMu_Fakerate::isFakeElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()){
		return false;
	}
	if(isFakeElectron(i)
			&& dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.1
			&& dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.03
			){
		return true;
	}
	return false;
}

double ZtoEMu_Fakerate::Electron_RelIso(unsigned int i){
	return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff_R04(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();
}

double ZtoEMu_Fakerate::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.3) return 0.194;
	else if(eta>=2.4) return 0.261;
}

double ZtoEMu_Fakerate::Electron_Aeff_R03(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.130;
	else if(eta>=1. && eta<1.479) return 0.137;
	else if(eta>=1.479 && eta<2.) return 0.067;
	else if(eta>=2. && eta<2.2) return 0.089;
	else if(eta>=2.2 && eta<2.3) return 0.107;
	else if(eta>=2.3 && eta<2.3) return 0.110;
	else if(eta>=2.4) return 0.138;
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
// Trigger & ID efficiencies
//

double ZtoEMu_Skim::MuonIDeff(unsigned int i){
	double pt = Ntp->Muons_p4(i).Pt();
	double eta = fabs(Ntp->Muons_p4(i).Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9771;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9746;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9644;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9891;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9548;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9701;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9766;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9892;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9648;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9836;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9820;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9909;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9676;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9817;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9886;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9883;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9883;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9833;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9910;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9900;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9826;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9841;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9900;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9886;
		}
	}
}

double ZtoEMu_Skim::MuonIDerr(unsigned int i){
	double pt = Ntp->Muons_p4(i).Pt();
	double eta = fabs(Ntp->Muons_p4(i).Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.0107;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0091;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0078;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0080;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.0046;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0049;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0049;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0049;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.0023;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0030;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0029;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0030;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.0012;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0018;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0018;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0019;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.0008;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0012;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0013;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0016;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.0005;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0004;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0003;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0004;
		}
	}
}

double ZtoEMu_Skim::MuonTriggerEff(unsigned int i){
	double pt = Ntp->Muons_p4(i).Pt();
	double eta = fabs(Ntp->Muons_p4(i).Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9829;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9745;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9943;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9158;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9850;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9852;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9743;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9333;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9951;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9610;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9716;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9459;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9869;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9779;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9665;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9501;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9959;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9881;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9932;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9391;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9986;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9540;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9549;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9386;
		}
	}
}

double ZtoEMu_Skim::MuonTriggerErr(unsigned int i){
	double pt = Ntp->Muons_p4(i).Pt();
	double eta = fabs(Ntp->Muons_p4(i).Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.0058;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0124;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0164;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0176;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.0056;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0171;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0179;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0162;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.0060;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0116;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0141;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0159;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.0074;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0187;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0184;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0251;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.0085;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0227;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0271;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0307;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.0087;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0165;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0211;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0209;
		}
	}
}

double ZtoEMu_Skim::ElectronIDeff(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.7654;
		}else if(eta>=0.8 && eta<1.5){
			return 0.7693;
		}else if(eta>=1.5 && eta<2.3){
			return 0.5719;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.8394;
		}else if(eta>=0.8 && eta<1.5){
			return 0.8457;
		}else if(eta>=1.5 && eta<2.3){
			return 0.7024;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.8772;
		}else if(eta>=0.8 && eta<1.5){
			return 0.8530;
		}else if(eta>=1.5 && eta<2.3){
			return 0.7631;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9006;
		}else if(eta>=0.8 && eta<1.5){
			return 0.8874;
		}else if(eta>=1.5 && eta<2.3){
			return 0.8092;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9261;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9199;
		}else if(eta>=1.5 && eta<2.3){
			return 0.8469;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9514;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9445;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9078;
		}
	}
}

double ZtoEMu_Skim::ElectronIDerr(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.0149;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0164;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0131;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.0045;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0061;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0075;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.0023;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0039;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0061;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.0018;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0017;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0024;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.0007;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0010;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0027;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.0002;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0003;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0007;
		}
	}
}

double ZtoEMu_Skim::ElectronTriggerEff(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9548;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9015;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9017;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9830;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9672;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9463;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9707;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9731;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9691;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9768;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9870;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9727;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 1.0047;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9891;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9858;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 1.0063;
		}else if(eta>=0.8 && eta<1.5){
			return 1.0047;
		}else if(eta>=1.5 && eta<2.3){
			return 1.0015;
		}
	}
}

double ZtoEMu_Skim::ElectronTriggerErr(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.0197;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0205;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0470;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.0115;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0113;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0212;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.0087;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0083;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0149;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.0084;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0083;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0162;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.0100;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0111;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0112;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.0078;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0073;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0135;
		}
	}
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
