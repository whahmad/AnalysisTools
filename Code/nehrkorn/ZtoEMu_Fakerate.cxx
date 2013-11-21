#include "ZtoEMu_Fakerate.h"
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

ZtoEMu_Fakerate::ZtoEMu_Fakerate(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(10)
  ,e_pt(10)
  ,mu_eta(2.4)
  ,e_eta(2.5)
{
    //verbose=true;
	//outfile = new TFile("ownfakerates.root","RECREATE");
}

ZtoEMu_Fakerate::~ZtoEMu_Fakerate(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu_Fakerate::~ZtoEMu_Fakerate Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu_Fakerate::~ZtoEMu_Fakerate()" << std::endl;
}

void  ZtoEMu_Fakerate::Configure(){

  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
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
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  tightmu=HConfig.GetTH2D(Name+"_tightmu","tightmu",100,0.,100.,100,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  tightmu_rebin=HConfig.GetTH2D(Name+"_tightmu_rebin","tightmu_rebin",5,10.,35.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  tighte=HConfig.GetTH2D(Name+"_tighte","tighte",100,0.,100.,100,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  tighte_rebin=HConfig.GetTH2D(Name+"_tighte_rebin","tighte_rebin",5,10.,35.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  fakemu=HConfig.GetTH2D(Name+"_fakemu","fakemu",100,0.,100.,100,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  fakemu_rebin=HConfig.GetTH2D(Name+"_fakemu_rebin","fakemu_rebin",5,10.,35.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  fakee=HConfig.GetTH2D(Name+"_fakee","fakee",100,0.,100.,100,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  fakee_rebin=HConfig.GetTH2D(Name+"_fakee_rebin","fakee_rebin",5,10.,35.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  mueff=HConfig.GetTH2D(Name+"_mueff","mueff",5,10.,35.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eeff=HConfig.GetTH2D(Name+"_eeff","eeff",5,10.,35.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  mudr=HConfig.GetTH1D(Name+"_mudr","mudr",100,0.,0.5,"#DeltaR(#mu_{trig},#mu_{reco})");
  mupt=HConfig.GetTH1D(Name+"_mupt","mupt",100,0.,5.,"|p_{T}^{#mu,trig}-p_{T}^{#mu,reco}|");
  edr=HConfig.GetTH1D(Name+"_edr","edr",100,0.,0.5,"#DeltaR(e_{trig},e_{reco})");
  ept=HConfig.GetTH1D(Name+"_ept","ept",100,0.,5.,"|p_{T}^{e,trig}-p_{T}^{e,reco}|");

  muleg_numerator=HConfig.GetTH2D(Name+"_muleg_numerator","muleg_numerator",100,0.,100.,100,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  muleg_denominator=HConfig.GetTH2D(Name+"_muleg_denominator","muleg_denominator",100,0.,100.,100,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_numerator=HConfig.GetTH2D(Name+"_eleg_numerator","eleg_numerator",100,0.,100.,100,-2.4,2.4,"p_{T}^{e} / GeV","#eta_{e}");
  eleg_denominator=HConfig.GetTH2D(Name+"_eleg_denominator","eleg_denominator",100,0.,100.,100,-2.4,2.4,"p_{T}^{e} / GeV","#eta_{e}");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_Fakerate::Store_ExtraDist(){
	Extradist2d.push_back(&tightmu);
	Extradist2d.push_back(&tightmu_rebin);
	Extradist2d.push_back(&tighte);
	Extradist2d.push_back(&tighte_rebin);
	Extradist2d.push_back(&fakemu);
	Extradist2d.push_back(&fakemu_rebin);
	Extradist2d.push_back(&fakee);
	Extradist2d.push_back(&fakee_rebin);
	Extradist2d.push_back(&mueff);
	Extradist2d.push_back(&eeff);

	Extradist1d.push_back(&mudr);
	Extradist1d.push_back(&mupt);
	Extradist1d.push_back(&edr);
	Extradist1d.push_back(&ept);

	Extradist2d.push_back(&muleg_numerator);
	Extradist2d.push_back(&muleg_denominator);
	Extradist2d.push_back(&eleg_numerator);
	Extradist2d.push_back(&eleg_denominator);
}

void  ZtoEMu_Fakerate::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Fakerate::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
  
  value.at(TriggerOk)=0;

  if(Ntp->TriggerAccept("HLT_QuadJet80_v")
		  || Ntp->TriggerAccept("HLT_Mu17_v")
		  || Ntp->TriggerAccept("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")
		  || Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdL")
		  || Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL")
		  || Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL")
		  || Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdL")
		  || Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL")
		  || Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL")
		  ){
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
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1),ww(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=34){
    w*=Ntp->EvtWeight3D();
    ww*=Ntp->EvtWeight3D();
    if(verbose)std::cout << "void  ZtoEMu_Fakerate::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;

  ///////////////////////////////////////////////////////////
  //
  // Fakerate calculation
  //

  // muon fakerate
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muons_p4(i).Pt()>=mu_pt &&
			  fabs(Ntp->Muons_p4(i).Eta())<mu_eta &&
			  vertex>=0 &&
			  (Ntp->TriggerAccept("HLT_Mu8_v") || Ntp->TriggerAccept("HLT_QuadJet80_v")) &&
			  isFakeMuon(i,vertex) &&
			  isTightMuon(i,vertex) &&
			  dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.02 &&
			  dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.5 &&
			  Ntp->isData()
			  ){
		  if(fabs(Ntp->Muons_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15){
			  tightmu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
		  }else if(fabs(Ntp->Muons_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.1){
			  tightmu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
		  }
	  }else if(isFakeMuon(i,vertex)
			  && !isTightMuon(i,vertex)
			  && Ntp->isData()
			  ){
		  fakemu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
	  }
  }

  // electron fakerate
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>=e_pt &&
			  fabs(Ntp->Electron_supercluster_eta(i))<e_eta &&
			  vertex>=0 &&
			  (Ntp->TriggerAccept("HLT_QuadJet80_v") || Ntp->TriggerAccept("HLT_Ele8_CaloIdT_TrkIdVL_v") || Ntp->TriggerAccept("HLT_Ele8_CaloIdL_CaloIsoVL_v")) &&
			  !Ntp->Electron_HasMatchedConversions(i) &&
			  Ntp->Electron_numberOfMissedHits(i)==0 &&
			  dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.2 &&
			  dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.02 &&
			  Ntp->isData()
			  ){
		  if(isFakeElectron(i,vertex) && isMVANonTrigElectron(i,vertex)){
			  tighte.at(t).Fill(Ntp->Electron_p4(i).Pt(),Ntp->Electron_supercluster_eta(i));
		  }else if(isFakeElectron(i,vertex) && !isMVANonTrigElectron(i,vertex)){
			  fakee.at(t).Fill(Ntp->Electron_p4(i).Pt(),Ntp->Electron_supercluster_eta(i));
		  }
	  }
  }

  ///////////////////////////////////////////////////////////
  //
  // Trigger efficiency
  //

  // save all tight muons and loose MVA electrons in event
  if(verbose) std::cout << "saving objects passing nominal identification" << std::endl;
  std::vector<unsigned int> SingleMuons;
  std::vector<unsigned int> ProbeMuons;
  std::vector<unsigned int> SingleElectrons;
  std::vector<unsigned int> ProbeElectrons;
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muons_p4(i).Pt()>=mu_pt &&
			  fabs(Ntp->Muons_p4(i).Eta())<mu_eta &&
			  vertex>=0 &&
			  isFakeMuon(i,vertex) &&
			  isTightMuon(i,vertex) &&
			  dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.02 &&
			  dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.5
			  ){
		  if(fabs(Ntp->Muons_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15){
			  SingleMuons.push_back(i);
		  }else if(fabs(Ntp->Muons_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.1){
			  SingleMuons.push_back(i);
		  }
	  }
  }
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>=e_pt &&
			  fabs(Ntp->Electron_supercluster_eta(i))<e_eta &&
			  vertex>=0 &&
			  !Ntp->Electron_HasMatchedConversions(i) &&
			  Ntp->Electron_numberOfMissedHits(i)==0 &&
			  dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.2 &&
			  dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.02
			  ){
		  if(isMVANonTrigElectron(i,vertex)){
			  SingleElectrons.push_back(i);
		  }
	  }
  }

  // match muons to trigger muon for electron leg measurement
  if(verbose) std::cout << "examining electron leg" << std::endl;
  if(Ntp->TriggerAccept("HLT_Mu17_v")){
	  // find muon that triggered
	  unsigned int tagmu = 999;
	  float testdr = 0.1;
	  float testpt = 1.;
	  TLorentzVector testmu;
	  for(unsigned i=0;i<SingleMuons.size();i++){
		  for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
			  if(Ntp->HLTTrigger_objs_trigger(j).find("HLT_Mu17_v") != string::npos){
				  for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
					  if(Ntp->HLTTrigger_objs_Id(j,k)==83){
						  testmu.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
								  Ntp->HLTTrigger_objs_Eta(j,k),
								  Ntp->HLTTrigger_objs_Phi(j,k),
								  Ntp->HLTTrigger_objs_E(j,k));
						  mudr.at(t).Fill(testmu.DeltaR(Ntp->Muons_p4(SingleMuons.at(i))),w);
						  mupt.at(t).Fill(fabs(testmu.Pt()-Ntp->Muons_p4(SingleMuons.at(i)).Pt()),w);
						  if(testmu.DeltaR(Ntp->Muons_p4(SingleMuons.at(i)))<testdr
								  && fabs(testmu.Pt()-Ntp->Muons_p4(SingleMuons.at(i)).Pt())<testpt
								  ){
							  tagmu=SingleMuons.at(i);
						  }
					  }
				  }
				  break;
			  }
		  }
	  }
	  // create collection of probe electrons
	  if(tagmu!=999){
		  for(unsigned i=0;i<SingleElectrons.size();i++){
			  if(Ntp->Electron_Charge(SingleElectrons.at(i))*Ntp->Muon_Charge(tagmu)<0. &&
					  (Ntp->Electron_p4(SingleElectrons.at(i))+Ntp->Muons_p4(tagmu)).M()>30. &&
					  (Ntp->Electron_p4(SingleElectrons.at(i))+Ntp->Muons_p4(tagmu)).M()<90.
					  ){
				  ProbeElectrons.push_back(SingleElectrons.at(i));
			  }
		  }
	  }
	  std::cout << "number of probe electrons: " << ProbeElectrons.size() << std::endl;
	  // pick probe electron with highest pt
	  float probeept = 10.;
	  unsigned int probee = 999;
	  for(unsigned i=0;i<ProbeElectrons.size();i++){
		  if(Ntp->Electron_p4(ProbeElectrons.at(i)).Pt()>probeept){
			  probeept = Ntp->Electron_p4(ProbeElectrons.at(i)).Pt();
			  probee = ProbeElectrons.at(i);
		  }
	  }
	  // match probe electron to trigger electron of cross trigger
	  if(probee!=999){
		  TLorentzVector teste;
		  for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
			  if(Ntp->HLTTrigger_objs_trigger(i).find("HLT_Mu17_Ele8_") != string::npos){
				  for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
					  if(Ntp->HLTTrigger_objs_Id(i,j)==82){
						  teste.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
								  Ntp->HLTTrigger_objs_Eta(i,j),
								  Ntp->HLTTrigger_objs_Phi(i,j),
								  Ntp->HLTTrigger_objs_E(i,j));
					  }
				  }
				  break;
			  }
		  }
		  // plot numerator and denominator
		  if(teste.DeltaR(Ntp->Electron_p4(probee))<testdr &&
				  fabs(teste.Pt()-Ntp->Electron_p4(probee).Pt())<testpt
				  ){
			  eleg_numerator.at(t).Fill(Ntp->Electron_p4(probee).Pt(),Ntp->Electron_supercluster_eta(probee),w);
		  }else{
			  eleg_denominator.at(t).Fill(Ntp->Electron_p4(probee).Pt(),Ntp->Electron_supercluster_eta(probee),w);
		  }
	  }
  }


  // match electrons to trigger electron for muon leg measurement
  if(verbose) std::cout << "examining muon leg" << std::endl;
  if(Ntp->TriggerAccept("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")){
	  // find electron that triggered
	  unsigned int tage = 999;
	  float testdr = 0.1;
	  float testpt = 1.;
	  TLorentzVector teste;
	  for(unsigned i=0;i<SingleElectrons.size();i++){
		  for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
			  if(Ntp->HLTTrigger_objs_trigger(j).find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") != string::npos){
				  for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
					  if(Ntp->HLTTrigger_objs_Id(j,k)==82){
						  teste.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
								  Ntp->HLTTrigger_objs_Eta(j,k),
								  Ntp->HLTTrigger_objs_Phi(j,k),
								  Ntp->HLTTrigger_objs_E(j,k));
						  edr.at(t).Fill(teste.DeltaR(Ntp->Electron_p4(SingleElectrons.at(i))),w);
						  ept.at(t).Fill(fabs(teste.Pt()-Ntp->Electron_p4(SingleElectrons.at(i)).Pt()),w);
						  if(teste.DeltaR(Ntp->Electron_p4(SingleElectrons.at(i)))<testdr
								  && fabs(teste.Pt()-Ntp->Electron_p4(SingleElectrons.at(i)).Pt())<testpt
								  ){
							  tage=SingleElectrons.at(i);
						  }
					  }
				  }
				  break;
			  }
		  }
	  }
	  // create collection of probe muons
	  if(tage!=999){
		  for(unsigned i=0;i<SingleMuons.size();i++){
			  if(Ntp->Muon_Charge(SingleMuons.at(i))*Ntp->Electron_Charge(tage)<0. &&
					  (Ntp->Muons_p4(SingleMuons.at(i))+Ntp->Electron_p4(tage)).M()>30. &&
					  (Ntp->Muons_p4(SingleMuons.at(i))+Ntp->Electron_p4(tage)).M()<90.
					  ){
				  ProbeMuons.push_back(SingleMuons.at(i));
			  }
		  }
	  }
	  std::cout << "number of probe muons: " << ProbeMuons.size() << std::endl;
	  // pick probe muon with highest pt
	  float probemupt = 10.;
	  unsigned int probemu = 999;
	  for(unsigned i=0;i<ProbeMuons.size();i++){
		  if(Ntp->Muons_p4(ProbeMuons.at(i)).Pt()>probemupt){
			  probemupt = Ntp->Muons_p4(ProbeMuons.at(i)).Pt();
			  probemu = ProbeMuons.at(i);
		  }
	  }
	  // match probe muon to trigger muon of cross trigger
	  if(probemu!=999){
		  TLorentzVector testmu;
		  for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
			  if(Ntp->HLTTrigger_objs_trigger(i).find("HLT_Mu8_Ele17_") != string::npos){
				  for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
					  if(Ntp->HLTTrigger_objs_Id(i,j)==83){
						  testmu.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
								  Ntp->HLTTrigger_objs_Eta(i,j),
								  Ntp->HLTTrigger_objs_Phi(i,j),
								  Ntp->HLTTrigger_objs_E(i,j));
					  }
				  }
				  break;
			  }
		  }
		  // plot numerator and denominator
		  if(testmu.DeltaR(Ntp->Muons_p4(probemu))<testdr &&
				  fabs(testmu.Pt()-Ntp->Muons_p4(probemu).Pt())<testpt
				  ){
			  muleg_numerator.at(t).Fill(Ntp->Muons_p4(probemu).Pt(),Ntp->Muons_p4(probemu).Eta(),w);
		  }else{
			  muleg_denominator.at(t).Fill(Ntp->Muons_p4(probemu).Pt(),Ntp->Muons_p4(probemu).Eta(),w);
		  }
	  }
  }

  if(verbose)std::cout << "ZtoEMu_Fakerate::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_Fakerate::isGoodVtx(unsigned int i){
	if(fabs(Ntp->Vtx(i).z())<24
			&& Ntp->Vtx(i).Perp()<2
			&& Ntp->Vtx_ndof(i)>4
			&& Ntp->Vtx_isFake(i)==0
			){
		return true;
	}
	return false;
}

bool ZtoEMu_Fakerate::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu_Fakerate::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu_Fakerate::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu_Fakerate::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

double ZtoEMu_Fakerate::vertexSignificance(TVector3 vec, unsigned int vertex){
	if(vertex<Ntp->NVtx() && vertex>=0){
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

int ZtoEMu_Fakerate::getxbin(double pt){
	if(pt<15.) return 1;
	else if(pt>=15. && pt<20.) return 2;
	else if(pt>=20. && pt<25.) return 3;
	else if(pt>=25. && pt<30.) return 4;
	else if(pt>=30. && pt<35.) return 5;
	else return 0;
}

int ZtoEMu_Fakerate::getybin(double eta,std::string object){
	if(object=="muon"){
		if(eta<0){
			if(fabs(eta)>=2.1) return 1;
			else if(fabs(eta)<2.1 && fabs(eta)>=1.2) return 2;
			else if(fabs(eta)<1.2 && fabs(eta)>=0.8) return 3;
			else if(fabs(eta)<0.8 && fabs(eta)>=0.) return 4;
		}
		else if(eta>=0){
			if(eta>=0. && eta<0.8) return 5;
			else if(eta>=0.8 && eta<1.2) return 6;
			else if(eta>=1.2 && eta<2.1) return 7;
			else if(eta>=2.1) return 8;
		}
	}
	else if(object=="electron"){
		if(eta<0){
			if(fabs(eta)>=2.) return 1;
			else if(fabs(eta)<2. && fabs(eta)>=1.479) return 2;
			else if(fabs(eta)<1.479 && fabs(eta)>=0.8) return 3;
			else if(fabs(eta)<0.8 && fabs(eta)>=0.) return 4;
		}
		else if(eta>=0){
			if(eta>=0. && eta<0.8) return 5;
			else if(eta>=0.8 && eta<1.479) return 6;
			else if(eta>=1.479 && eta<2.) return 7;
			else if(eta>=2.) return 8;
		}
	}
	else{
		return 0;
	}
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_Fakerate::isTightMuon(unsigned int i){
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

bool ZtoEMu_Fakerate::isTightMuon(unsigned int i, unsigned int j){
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

bool ZtoEMu_Fakerate::isLooseMuon(unsigned int i){
	if(Ntp->Muon_isPFMuon(i) &&
			(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))
			){
		return true;
	}
	return false;
}

bool ZtoEMu_Fakerate::isFakeMuon(unsigned int i){
	if(Ntp->Muon_isGlobalMuon(i) &&
			Ntp->Muons_p4(i).Pt()>10 &&
			fabs(Ntp->Muons_p4(i).Eta())<2.5
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

bool ZtoEMu_Fakerate::isFakeMuon(unsigned int i, unsigned int j){
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

double ZtoEMu_Fakerate::Muon_RelIso(unsigned int i){
	return (Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i)))/Ntp->Muons_p4(i).Pt();
}

double ZtoEMu_Fakerate::Muon_AbsIso(unsigned int i){
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
// Finish function
//

void ZtoEMu_Fakerate::Finish(){
	for(unsigned int i=1;i<=fakemu.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<fakemu.at(0).GetNbinsY();j++){
			int xbin = getxbin(fakemu.at(0).GetXaxis()->GetBinLowEdge(i));
			int ybin = getybin(fakemu.at(0).GetYaxis()->GetBinLowEdge(j),"muon");
			fakemu_rebin.at(0).SetBinContent(xbin,ybin,fakemu.at(0).GetBinContent(i,j)+fakemu_rebin.at(0).GetBinContent(xbin,ybin));
			xbin = getxbin(tightmu.at(0).GetXaxis()->GetBinLowEdge(i));
			ybin = getybin(tightmu.at(0).GetYaxis()->GetBinLowEdge(j),"muon");
			tightmu_rebin.at(0).SetBinContent(xbin,ybin,tightmu.at(0).GetBinContent(i,j)+tightmu_rebin.at(0).GetBinContent(xbin,ybin));
			
			xbin = getxbin(fakee.at(0).GetXaxis()->GetBinLowEdge(i));
			ybin = getybin(fakee.at(0).GetYaxis()->GetBinLowEdge(j),"muon");
			fakee_rebin.at(0).SetBinContent(xbin,ybin,fakee.at(0).GetBinContent(i,j)+fakee_rebin.at(0).GetBinContent(xbin,ybin));
			xbin = getxbin(tighte.at(0).GetXaxis()->GetBinLowEdge(i));
			ybin = getybin(tighte.at(0).GetYaxis()->GetBinLowEdge(j),"muon");
			tighte_rebin.at(0).SetBinContent(xbin,ybin,tighte.at(0).GetBinContent(i,j)+tighte_rebin.at(0).GetBinContent(xbin,ybin));
		}
	}
	double muprob;
	double eprob;
	for(unsigned int i=1;i<=fakemu_rebin.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<fakemu_rebin.at(0).GetNbinsY();j++){
			if(fakemu_rebin.at(0).GetBinContent(i,j)!=0){
				muprob = tightmu_rebin.at(0).GetBinContent(i,j)/fakemu_rebin.at(0).GetBinContent(i,j);
				mueff.at(0).SetBinContent(i,j,muprob/(1.-muprob));
			}
		}
	}
	for(unsigned int i=1;i<=fakee_rebin.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<=fakee_rebin.at(0).GetNbinsY();j++){
			if(fakee_rebin.at(0).GetBinContent(i,j)!=0){
				eprob = tighte_rebin.at(0).GetBinContent(i,j)/fakee_rebin.at(0).GetBinContent(i,j);
				eeff.at(0).SetBinContent(i,j,eprob/(1.-eprob));
			}
		}
	}
	//mueff.at(0).Write();
	//eeff.at(0).Write();
	//outfile->Write();
	//outfile->Close();
	Selection::Finish();
}
