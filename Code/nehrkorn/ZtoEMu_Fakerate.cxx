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
  ,mu_pt(8)
  ,e_pt(8)
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
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  tightmu=HConfig.GetTH2D(Name+"_tightmu","tightmu",90,10.,100.,48,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  tightmu_rebin=HConfig.GetTH2D(Name+"_tightmu_rebin","tightmu_rebin",18,10.,100.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  tighte=HConfig.GetTH2D(Name+"_tighte","tighte",90,10.,100.,50,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  tighte_rebin=HConfig.GetTH2D(Name+"_tighte_rebin","tighte_rebin",18,10.,100.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  fakemu=HConfig.GetTH2D(Name+"_fakemu","fakemu",90,10.,100.,48,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  fakemu_rebin=HConfig.GetTH2D(Name+"_fakemu_rebin","fakemu_rebin",18,10.,100.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  fakee=HConfig.GetTH2D(Name+"_fakee","fakee",90,10.,100.,50,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  fakee_rebin=HConfig.GetTH2D(Name+"_fakee_rebin","fakee_rebin",18,10.,100.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  mueff=HConfig.GetTH2D(Name+"_mueff","mueff",18,10.,100.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eeff=HConfig.GetTH2D(Name+"_eeff","eeff",18,10.,100.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_numerator=HConfig.GetTH2D(Name+"_muleg_numerator","muleg_numerator",90,10.,100.,48,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  muleg_denominator=HConfig.GetTH2D(Name+"_muleg_denominator","muleg_denominator",90,10.,100.,48,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_numerator=HConfig.GetTH2D(Name+"_eleg_numerator","eleg_numerator",90,10.,100.,50,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  eleg_denominator=HConfig.GetTH2D(Name+"_eleg_denominator","eleg_denominator",90,10.,100.,50,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_numerator_rebin=HConfig.GetTH2D(Name+"_muleg_numerator_rebin","muleg_numerator_rebin",18,10.,100.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  muleg_denominator_rebin=HConfig.GetTH2D(Name+"_muleg_denominator_rebin","muleg_denominator_rebin",18,10.,100.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_numerator_rebin=HConfig.GetTH2D(Name+"_eleg_numerator_rebin","eleg_numerator_rebin",18,10.,100.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  eleg_denominator_rebin=HConfig.GetTH2D(Name+"_eleg_denominator_rebin","eleg_denominator_rebin",18,10.,100.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_eff=HConfig.GetTH2D(Name+"_muleg_eff","muleg_eff",18,10.,100.,8,-2.4,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_eff=HConfig.GetTH2D(Name+"_eleg_eff","eleg_eff",18,10.,100.,8,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  tagmupt=HConfig.GetTH1D(Name+"_tagmupt","tagmupt",40,0.,100.,"p_{T}^{#mu} / GeV");
  tagmueta=HConfig.GetTH1D(Name+"_tagmueta","tagmueta",25,-2.5,2.5,"#eta_{#mu}");
  tagept=HConfig.GetTH1D(Name+"_tagept","tagept",40,0.,100.,"p_{T}^{e} / GeV");
  tageeta=HConfig.GetTH1D(Name+"_tageeta","tageeta",25,-2.5,2.5,"#eta_{e}");
  probemupt=HConfig.GetTH1D(Name+"_probemupt","probemupt",40,0.,100.,"p_{T}^{#mu} / GeV");
  probemueta=HConfig.GetTH1D(Name+"_probemueta","probemueta",25,-2.5,2.5,"#eta_{#mu}");
  probeept=HConfig.GetTH1D(Name+"_probeept","probeept",40,0.,100.,"p_{T}^{e} / GeV");
  probeeeta=HConfig.GetTH1D(Name+"_probeeeta","probeeeta",25,-2.5,2.5,"#eta_{e}");
  drtagmuprobee=HConfig.GetTH1D(Name+"_drtagmuprobee","drtagmuprobee",20,0.,1.,"#DeltaR(e,#mu)");
  drtageprobemu=HConfig.GetTH1D(Name+"_drtageprobemu","drtageprobemu",20,0.,1.,"#DeltaR(e,#mu)");
  drmumu=HConfig.GetTH1D(Name+"_drmumu","drmumu",60,0.,3.,"#DeltaR(#mu,#mu)");
  dree=HConfig.GetTH1D(Name+"_dree","dree",60,0.,3.,"#DeltaR(e,e)");
  ptbaltagmuprobee=HConfig.GetTH1D(Name+"_ptbaltagmuprobee","ptbaltagmuprobee",40,0.,200.,"p_{T}^{e#mu} / GeV");
  ptbaltageprobemu=HConfig.GetTH1D(Name+"_ptbaltageprobemu","ptbaltageprobemu",40,0.,200.,"p_{T}^{e#mu} / GeV");
  ptbalmumu=HConfig.GetTH1D(Name+"_ptbalmumu","ptbalmumu",40,0.,200.,"p_{T}^{#mu#mu} / GeV");
  ptbalee=HConfig.GetTH1D(Name+"_ptbalee","ptbalee",40,0.,200.,"p_{T}^{ee} / GeV");
  mttagmu=HConfig.GetTH1D(Name+"_mttagmu","mttagmu",40,0.,160.,"m_{t}^{#mu} / GeV");
  mttage=HConfig.GetTH1D(Name+"_mttage","mttage",40,0.,160.,"m_{t}^{e} / GeV");
  mtprobemu=HConfig.GetTH1D(Name+"_mtprobemu","mtprobemu",40,0.,160.,"m_{t}^{#mu} / GeV");
  mtprobee=HConfig.GetTH1D(Name+"_mtprobee","mtprobee",40,0.,160.,"m_{t}^{e} / GeV");
  mtagmuprobee=HConfig.GetTH1D(Name+"_mtagmuprobee","mtagmuprobee",30,30.,90.,"m_{e,#mu} / GeV");
  mtageprobemu=HConfig.GetTH1D(Name+"_mtageprobemu","mtageprobemu",30,30.,90.,"m_{e,#mu} / GeV");
  mmumu=HConfig.GetTH1D(Name+"_mmumu","mmumu",30,60.,120.,"m_{#mu,#mu} / GeV");
  mee=HConfig.GetTH1D(Name+"_mee","mee",30,60.,120.,"m_{e,e} / GeV");

  sip=HConfig.GetTH1D(Name+"_sip","sip",10,0,10,"sip");
  nmu=HConfig.GetTH1D(Name+"_nmu","nmu",10,0,10,"nmu");
  ne=HConfig.GetTH1D(Name+"_ne","ne",10,0,10,"ne");
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

	Extradist2d.push_back(&muleg_numerator);
	Extradist2d.push_back(&muleg_denominator);
	Extradist2d.push_back(&eleg_numerator);
	Extradist2d.push_back(&eleg_denominator);

	Extradist2d.push_back(&muleg_numerator_rebin);
	Extradist2d.push_back(&muleg_denominator_rebin);
	Extradist2d.push_back(&eleg_numerator_rebin);
	Extradist2d.push_back(&eleg_denominator_rebin);

	Extradist2d.push_back(&muleg_eff);
	Extradist2d.push_back(&eleg_eff);

	Extradist1d.push_back(&tagmupt);
	Extradist1d.push_back(&tagmueta);
	Extradist1d.push_back(&tagept);
	Extradist1d.push_back(&tageeta);
	Extradist1d.push_back(&probemupt);
	Extradist1d.push_back(&probemueta);
	Extradist1d.push_back(&probeept);
	Extradist1d.push_back(&probeeeta);
	Extradist1d.push_back(&drtagmuprobee);
	Extradist1d.push_back(&drtageprobemu);
	Extradist1d.push_back(&drmumu);
	Extradist1d.push_back(&dree);
	Extradist1d.push_back(&ptbaltagmuprobee);
	Extradist1d.push_back(&ptbaltageprobemu);
	Extradist1d.push_back(&ptbalmumu);
	Extradist1d.push_back(&ptbalee);
	Extradist1d.push_back(&mttagmu);
	Extradist1d.push_back(&mttage);
	Extradist1d.push_back(&mtprobemu);
	Extradist1d.push_back(&mtprobee);
	Extradist1d.push_back(&mtagmuprobee);
	Extradist1d.push_back(&mtageprobemu);
	Extradist1d.push_back(&mmumu);
	Extradist1d.push_back(&mee);

	Extradist1d.push_back(&sip);
	Extradist1d.push_back(&nmu);
	Extradist1d.push_back(&ne);
}

void  ZtoEMu_Fakerate::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Fakerate::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
  
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
	  if(Ntp->Muons_p4(i).Pt()>mu_pt
			  && fabs(Ntp->Muons_p4(i).Eta())<mu_eta
			  && vertex>=0
			  && (Ntp->TriggerAccept("HLT_Mu8_v")
					  || Ntp->TriggerAccept("HLT_Mu12_v")
					  || Ntp->TriggerAccept("HLT_Mu17_v")
					  || Ntp->TriggerAccept("HLT_Mu40_v")
					  || Ntp->TriggerAccept("HLT_QuadJet80_v"))
					  && Ntp->isData()
			  ){
		  if(isFakeMuon(i,vertex)){
			  fakemu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
			  fakemu_rebin.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
			  if(isTightMuon(i,vertex)
					  && Muon_RelIso(i)<0.12){
				  /*if(fabs(Ntp->Muons_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15){
					  tightmu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
					  tightmu_rebin.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
				  }else if(fabs(Ntp->Muons_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.1){
					  tightmu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
					  tightmu_rebin.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
				  }*/
				  tightmu.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
				  tightmu_rebin.at(t).Fill(Ntp->Muons_p4(i).Pt(),Ntp->Muons_p4(i).Eta());
			  }
		  }
	  }
  }

  // electron fakerate
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>e_pt
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && vertex>=0
			  && (Ntp->TriggerAccept("HLT_Ele8_CaloIdT_TrkIdVL_v")
					  || Ntp->TriggerAccept("HLT_Ele8_CaloIdL_CaloIsoVL_v")
					  || Ntp->TriggerAccept("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")
					  || Ntp->TriggerAccept("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_v")
					  || Ntp->TriggerAccept("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_BTagIPIter_v")
					  || Ntp->TriggerAccept("HLT_Ele27_WP80_v")
					  || Ntp->TriggerAccept("HLT_QuadJet80_v"))
					  && Ntp->isData()
			  ){
		  if(isFakeElectron(i,vertex)){
			  fakee.at(t).Fill(Ntp->Electron_p4(i).Pt(),Ntp->Electron_supercluster_eta(i));
			  fakee_rebin.at(t).Fill(Ntp->Electron_p4(i).Pt(),Ntp->Electron_supercluster_eta(i));
			  if(isMVANonTrigElectron(i,vertex)// &&
					  //dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.1 &&
					  //dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.02
					  ){
				  tighte.at(t).Fill(Ntp->Electron_p4(i).Pt(),Ntp->Electron_supercluster_eta(i));
				  tighte_rebin.at(t).Fill(Ntp->Electron_p4(i).Pt(),Ntp->Electron_supercluster_eta(i));
			  }
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
	  if(Ntp->Muons_p4(i).Pt()>mu_pt
			  && fabs(Ntp->Muons_p4(i).Eta())<mu_eta
			  && vertex>=0
			  && isTightMuon(i,vertex)
			  && Muon_RelIso(i)<0.12
			  ){
		  /*if(fabs(Ntp->Muons_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15){
			  SingleMuons.push_back(i);
		  }else if(fabs(Ntp->Muons_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.1){
			  SingleMuons.push_back(i);
		  }*/
		  SingleMuons.push_back(i);
	  }
  }
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(vertex>=0)sip.at(t).Fill(vertexSignificance(Ntp->Electron_Poca(i),vertex),w);
	  if(Ntp->Electron_p4(i).Et()>e_pt
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && vertex>=0// &&
			  //dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.1 &&
			  //dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(vertex))<0.02
			  ){
		  if(isMVANonTrigElectron(i,vertex)){
			  SingleElectrons.push_back(i);
		  }
	  }
  }

  /*triggerMatch(SingleMuons,SingleElectrons,"HLT_Mu17_v","HLT_Mu17_Ele8_","muon",eleg_denominator.at(t),eleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_IsoMu24_eta2p1_v","HLT_Mu17_Ele8_","muon",eleg_denominator.at(t),eleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_IsoMu24_v","HLT_Mu17_Ele8_","muon",eleg_denominator.at(t),eleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
  		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_IsoMu30_v","HLT_Mu17_Ele8_","muon",eleg_denominator.at(t),eleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
  		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_Mu40_v","HLT_Mu17_Ele8_","muon",eleg_denominator.at(t),eleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
  		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_Mu17_v","HLT_Mu17_Ele8_","muon",eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_IsoMu24_eta2p1_v","HLT_Mu17_Ele8_","muon",eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_IsoMu24_v","HLT_Mu17_Ele8_","muon",eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_IsoMu30_v","HLT_Mu17_Ele8_","muon",eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),w);
  triggerMatch(SingleMuons,SingleElectrons,"HLT_Mu40_v","HLT_Mu17_Ele8_","muon",eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),w);

  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","HLT_Mu8_Ele17_","electron",muleg_denominator.at(t),muleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_v","HLT_Mu8_Ele17_","electron",muleg_denominator.at(t),muleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
  		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_BTagIPIter_v","HLT_Mu8_Ele17_","electron",muleg_denominator.at(t),muleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
  		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele27_WP80_v","HLT_Mu8_Ele17_","electron",muleg_denominator.at(t),muleg_numerator.at(t),tagmupt.at(t),tagmueta.at(t),tagept.at(t),tageeta.at(t),probemupt.at(t),probemueta.at(t),probeept.at(t),probeeeta.at(t),
    		  drtagmuprobee.at(t),drtageprobemu.at(t),ptbaltagmuprobee.at(t),ptbaltageprobemu.at(t),mttagmu.at(t),mttage.at(t),mtprobemu.at(t),mtprobee.at(t),mtagmuprobee.at(t),mtageprobemu.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","HLT_Mu8_Ele17_","electron",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_v","HLT_Mu8_Ele17_","electron",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_BTagIPIter_v","HLT_Mu8_Ele17_","electron",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),w);
  triggerMatch(SingleElectrons,SingleMuons,"HLT_Ele27_WP80_v","HLT_Mu8_Ele17_","electron",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),w);*/

  doubleMuE(SingleMuons,"muon","HLT_Mu17_v","HLT_Mu17_Mu8_v",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),tagmupt.at(t),tagmueta.at(t),probemupt.at(t),probemueta.at(t),
		  drmumu.at(t),ptbalmumu.at(t),mttagmu.at(t),mtprobemu.at(t),mmumu.at(t),w);
  doubleMuE(SingleMuons,"muon","HLT_Mu40_v","HLT_Mu17_Mu8_v",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),tagmupt.at(t),tagmueta.at(t),probemupt.at(t),probemueta.at(t),
  		  drmumu.at(t),ptbalmumu.at(t),mttagmu.at(t),mtprobemu.at(t),mmumu.at(t),w);
  doubleMuE(SingleMuons,"muon","HLT_IsoMu24_v","HLT_Mu17_Mu8_v",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),tagmupt.at(t),tagmueta.at(t),probemupt.at(t),probemueta.at(t),
  		  drmumu.at(t),ptbalmumu.at(t),mttagmu.at(t),mtprobemu.at(t),mmumu.at(t),w);
  doubleMuE(SingleMuons,"muon","HLT_IsoMu24_eta2p1_v","HLT_Mu17_Mu8_v",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),tagmupt.at(t),tagmueta.at(t),probemupt.at(t),probemueta.at(t),
    		  drmumu.at(t),ptbalmumu.at(t),mttagmu.at(t),mtprobemu.at(t),mmumu.at(t),w);
  doubleMuE(SingleMuons,"muon","HLT_IsoMu30_v","HLT_Mu17_Mu8_v",muleg_denominator_rebin.at(t),muleg_numerator_rebin.at(t),tagmupt.at(t),tagmueta.at(t),probemupt.at(t),probemueta.at(t),
    		  drmumu.at(t),ptbalmumu.at(t),mttagmu.at(t),mtprobemu.at(t),mmumu.at(t),w);

  doubleMuE(SingleElectrons,"electron","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
		  eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),tagept.at(t),tageeta.at(t),probeept.at(t),probeeeta.at(t),dree.at(t),ptbalee.at(t),mttage.at(t),mtprobee.at(t),mee.at(t),w);
  doubleMuE(SingleElectrons,"electron","HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_v","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
  		  eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),tagept.at(t),tageeta.at(t),probeept.at(t),probeeeta.at(t),dree.at(t),ptbalee.at(t),mttage.at(t),mtprobee.at(t),mee.at(t),w);
  doubleMuE(SingleElectrons,"electron","HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_BTagIPIter_v","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
  		  eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),tagept.at(t),tageeta.at(t),probeept.at(t),probeeeta.at(t),dree.at(t),ptbalee.at(t),mttage.at(t),mtprobee.at(t),mee.at(t),w);
  doubleMuE(SingleElectrons,"electron","HLT_Ele27_WP80_v","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",
  		  eleg_denominator_rebin.at(t),eleg_numerator_rebin.at(t),tagept.at(t),tageeta.at(t),probeept.at(t),probeeeta.at(t),dree.at(t),ptbalee.at(t),mttage.at(t),mtprobee.at(t),mee.at(t),w);

  if(verbose)std::cout << "ZtoEMu_Fakerate::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_Fakerate::isGoodVtx(unsigned int i){
	if(fabs(Ntp->Vtx(i).z())>=24) return false;
	if(Ntp->Vtx(i).Perp()>=2) return false;
	if(Ntp->Vtx_ndof(i)<=4) return false;
	if(Ntp->Vtx_isFake(i)!=0) return false;
	return true;
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
	return 0;
}

void ZtoEMu_Fakerate::doubleMuE(std::vector<unsigned int> objects, std::string particle,std::string tagtrigger, std::string probetrigger, TH2D &denominator, TH2D &numerator,
		  TH1D &tagpt, TH1D &tageta, TH1D &probept, TH1D &probeeta, TH1D &dr, TH1D &ptbal, TH1D &mttag, TH1D &mtprobe, TH1D &m, double w){
	unsigned int tag = 999;
	bool flag;
	float testdr;
	std::vector<unsigned int> probes;
	unsigned int probe = 999;
	// match tag to trigger for probe leg measurement
	if(Ntp->TriggerAccept(tagtrigger)){
		flag = false;
		// find particle that triggered (dR<0.5 and smallest dR of all trigger objects)
		TLorentzVector testtag(0.,0.,0.,0.);
		for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
			if(Ntp->HLTTrigger_objs_trigger(i).find(tagtrigger)!=string::npos){
				testdr = 0.5;
				for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
					if(particle=="muon"){
						if(Ntp->HLTTrigger_objs_Id(i,j)==83){
							testtag.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
									Ntp->HLTTrigger_objs_Eta(i,j),
									Ntp->HLTTrigger_objs_Phi(i,j),
									Ntp->HLTTrigger_objs_E(i,j));
							if(testtag.Pt()>0.){
								for(unsigned k=0;k<objects.size();k++){
									if(testtag.DeltaR(Ntp->Muons_p4(objects.at(k)))<testdr){
										tag = objects.at(k);
										testdr = testtag.DeltaR(Ntp->Muons_p4(objects.at(k)));
										flag = true;
									}
								}
							}
						}
					}
					if(particle=="electron"){
						if(Ntp->HLTTrigger_objs_Id(i,j)==82){
							testtag.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
									Ntp->HLTTrigger_objs_Eta(i,j),
									Ntp->HLTTrigger_objs_Phi(i,j),
									Ntp->HLTTrigger_objs_E(i,j));
							if(testtag.Pt()>0.){
								for(unsigned k=0;k<objects.size();k++){
									if(testtag.DeltaR(Ntp->Electron_p4(objects.at(k)))<testdr){
										tag = objects.at(k);
										testdr = testtag.DeltaR(Ntp->Electron_p4(objects.at(k)));
										flag = true;
									}
								}
							}
						}
					}
					if(flag) break;
				}
			}
			if(flag) break;
		}
		// create collection of probes
		if(tag!=999){
			if(particle=="muon"){
				tagpt.Fill(Ntp->Muons_p4(tag).Pt(),w);
				tageta.Fill(Ntp->Muons_p4(tag).Eta(),w);
				mttag.Fill(sqrt(2*Ntp->Muons_p4(tag).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muons_p4(tag).Px(),Ntp->Muons_p4(tag).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
				for(unsigned i=0;i<objects.size();i++){
					if(i!=tag &&
							Ntp->Muon_Charge(tag)+Ntp->Muon_Charge(i)==0
							&& (Ntp->Muons_p4(tag)+Ntp->Muons_p4(i)).M()>60.
							&& (Ntp->Muons_p4(tag)+Ntp->Muons_p4(i)).M()<120.){
						probes.push_back(objects.at(i));
					}
				}
			}
			if(particle=="electron"){
				tagpt.Fill(Ntp->Electron_p4(tag).Pt(),w);
				tageta.Fill(Ntp->Electron_supercluster_eta(tag),w);
				mttag.Fill(sqrt(2*Ntp->Electron_p4(tag).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(tag).Px(),Ntp->Electron_p4(tag).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
				for(unsigned i=0;i<objects.size();i++){
					if(i!=tag &&
							Ntp->Electron_Charge(tag)+Ntp->Electron_Charge(i)==0
							&& (Ntp->Electron_p4(tag)+Ntp->Electron_p4(i)).M()>60.
							&& (Ntp->Electron_p4(tag)+Ntp->Electron_p4(i)).M()<120.){
						probes.push_back(objects.at(i));
					}
				}
			}
		}
		// pick probe electron with highest pt
		float testpt = 8;
		unsigned int probe = 999;
		for(unsigned i=0;i<probes.size();i++){
			if(particle=="muon"){
				if(Ntp->Muons_p4(probes.at(i)).Pt()>testpt){
					testpt = Ntp->Muons_p4(probes.at(i)).Pt();
					probe = probes.at(i);
				}
			}
			if(particle=="electron"){
				if(Ntp->Electron_p4(probes.at(i)).Pt()>testpt){
					testpt = Ntp->Electron_p4(probes.at(i)).Pt();
					probe = probes.at(i);
				}
			}
		}
		// match probe(s) to trigger of probe trigger
		if(probe!=999){
			if(particle=="muon"){
				probept.Fill(Ntp->Muons_p4(probe).Pt(),w);
				probeeta.Fill(Ntp->Muons_p4(probe).Eta(),w);
				mtprobe.Fill(sqrt(2*Ntp->Muons_p4(probe).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muons_p4(probe).Px(),Ntp->Muons_p4(probe).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
				dr.Fill(Ntp->Muons_p4(tag).DeltaR(Ntp->Muons_p4(probe)),w);
				ptbal.Fill((Ntp->Muons_p4(tag)+Ntp->Muons_p4(probe)).Pt(),w);
				m.Fill((Ntp->Muons_p4(tag)+Ntp->Muons_p4(probe)).M(),w);
			}
			if(particle=="electron"){
				probept.Fill(Ntp->Electron_p4(probe).Pt(),w);
				probeeta.Fill(Ntp->Electron_supercluster_eta(probe),w);
				mtprobe.Fill(sqrt(2*Ntp->Electron_p4(probe).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(probe).Px(),Ntp->Electron_p4(probe).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
				dr.Fill(Ntp->Electron_p4(tag).DeltaR(Ntp->Electron_p4(probe)),w);
				ptbal.Fill((Ntp->Electron_p4(tag)+Ntp->Electron_p4(probe)).Pt(),w);
				m.Fill((Ntp->Electron_p4(tag)+Ntp->Electron_p4(probe)).M(),w);
			}
			TLorentzVector testprobe(0.,0.,0.,0.);
			flag = false;
			testdr = 0.5;
			for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
				if(Ntp->HLTTrigger_objs_trigger(i).find(probetrigger) != string::npos){
					for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
						if(particle=="muon"){
							if(Ntp->HLTTrigger_objs_Id(i,j)==83){
								testprobe.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
										Ntp->HLTTrigger_objs_Eta(i,j),
										Ntp->HLTTrigger_objs_Phi(i,j),
										Ntp->HLTTrigger_objs_E(i,j));
								//if(testprobe.DeltaR(testtag)>testdr)flag = true;
							}
						}
						if(particle=="electron"){
							if(Ntp->HLTTrigger_objs_Id(i,j)==82){
								testprobe.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
										Ntp->HLTTrigger_objs_Eta(i,j),
										Ntp->HLTTrigger_objs_Phi(i,j),
										Ntp->HLTTrigger_objs_E(i,j));
								//if(testprobe.DeltaR(testtag)>testdr)flag = true;
							}
						}
						if(flag) break;
					}
				}
				if(flag) break;
			}
			if(testprobe.Pt()>0.){
				if(particle=="muon"){
					denominator.Fill(Ntp->Muons_p4(probe).Pt(),Ntp->Muons_p4(probe).Eta(),w);
					if(testprobe.DeltaR(Ntp->Muons_p4(probe))<testdr){
						numerator.Fill(Ntp->Muons_p4(probe).Pt(),Ntp->Muons_p4(probe).Eta(),w);
					}
				}
				if(particle=="electron"){
					denominator.Fill(Ntp->Electron_p4(probe).Pt(),Ntp->Electron_supercluster_eta(probe),w);
					if(testprobe.DeltaR(Ntp->Electron_p4(probe))<testdr){
						numerator.Fill(Ntp->Electron_p4(probe).Pt(),Ntp->Electron_supercluster_eta(probe),w);
					}
				}
			}
		}
	}
}

void ZtoEMu_Fakerate::triggerMatch(std::vector<unsigned int> tags, std::vector<unsigned int> probes, std::string tagtrigger, std::string probetrigger, std::string tagname, TH2D &denominator, TH2D &numerator, double w){
	unsigned int finalprobe = 999;
	std::vector<unsigned int> Probe;
	bool flag;
	// match tag-particle to trigger for probe leg measurement
	if(Ntp->TriggerAccept(tagtrigger)){
	  flag = false;
	  // find muon that triggered
	  unsigned int tag = 999;
	  float testdr = 0.5;
	  TLorentzVector testtag(0.,0.,0.,0.);
	  for(unsigned i=0;i<tags.size();i++){
		  for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
			  if(Ntp->HLTTrigger_objs_trigger(j).find(tagtrigger) != string::npos){
				  for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
					  if(tagname=="muon"){
						  if(Ntp->HLTTrigger_objs_Id(j,k)==83){
							  testtag.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
								  Ntp->HLTTrigger_objs_Eta(j,k),
								  Ntp->HLTTrigger_objs_Phi(j,k),
								  Ntp->HLTTrigger_objs_E(j,k));
							  if(testtag.Pt()>0.
									  && testtag.DeltaR(Ntp->Muons_p4(tags.at(i)))<testdr
									  ){
								  tag = tags.at(i);
								  flag = true;
							  }
						  }
					  }
					  if(tagname=="electron"){
						  if(Ntp->HLTTrigger_objs_Id(j,k)==82){
							  testtag.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
								  Ntp->HLTTrigger_objs_Eta(j,k),
								  Ntp->HLTTrigger_objs_Phi(j,k),
								  Ntp->HLTTrigger_objs_E(j,k));
							  if(testtag.Pt()>0.
									  && testtag.DeltaR(Ntp->Electron_p4(tags.at(i)))<testdr
									  ){
								  tag = tags.at(i);
								  flag = true;
							  }
						  }
					  }
					  if(flag) break;
				  }
			  }
			  if(flag) break;
		  }
	  }
	  // create collection of probes
	  if(tag!=999){
		  if(tagname=="muon"){
			  for(unsigned i=0;i<probes.size();i++){
				  if(Ntp->Electron_Charge(probes.at(i))+Ntp->Muon_Charge(tag)==0
						  && (Ntp->Electron_p4(probes.at(i))+Ntp->Muons_p4(tag)).M()>30.
						  && (Ntp->Electron_p4(probes.at(i))+Ntp->Muons_p4(tag)).M()<90.
						  ){
							  Probe.push_back(probes.at(i));
						  }
			  }
		  }
		  if(tagname=="electron"){
			  for(unsigned i=0;i<probes.size();i++){
				  if(Ntp->Muon_Charge(probes.at(i))+Ntp->Electron_Charge(tag)==0
						  && (Ntp->Muons_p4(probes.at(i))+Ntp->Electron_p4(tag)).M()>30.
						  && (Ntp->Muons_p4(probes.at(i))+Ntp->Electron_p4(tag)).M()<90.
						  ){
					  Probe.push_back(probes.at(i));
				  }
			  }
		  }
	  }
	  // pick probe electron with highest pt
	  float probept = 8.;
	  unsigned int probe = 999;
	  for(unsigned i=0;i<Probe.size();i++){
		  if(tagname=="muon"){
			  if(Ntp->Electron_p4(Probe.at(i)).Pt()>probept){
				  probept = Ntp->Electron_p4(Probe.at(i)).Pt();
				  probe = Probe.at(i);
			  }
		  }
		  if(tagname=="electron"){
			  if(Ntp->Muons_p4(Probe.at(i)).Pt()>probept){
				  probept = Ntp->Muons_p4(Probe.at(i)).Pt();
				  probe = Probe.at(i);
			  }
		  }
	  }
	  // match probe particle to trigger of probe trigger
	  if(probe!=999){
		  TLorentzVector testprobe(0.,0.,0.,0.);
		  flag = false;
		  for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
			  if(Ntp->HLTTrigger_objs_trigger(i).find(probetrigger) != string::npos){
				  for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
					  if(tagname=="muon"){
						  if(Ntp->HLTTrigger_objs_Id(i,j)==82){
							  testprobe.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
									  Ntp->HLTTrigger_objs_Eta(i,j),
									  Ntp->HLTTrigger_objs_Phi(i,j),
									  Ntp->HLTTrigger_objs_E(i,j));
							  flag = true;
						  }
					  }
					  if(tagname=="electron"){
						  if(Ntp->HLTTrigger_objs_Id(i,j)==83){
							  testprobe.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
									  Ntp->HLTTrigger_objs_Eta(i,j),
									  Ntp->HLTTrigger_objs_Phi(i,j),
									  Ntp->HLTTrigger_objs_E(i,j));
							  flag = true;
						  }
					  }
					  if(flag) break;
				  }
			  }
			  if(flag) break;
		  }
		  if(tagname=="muon"){
			  denominator.Fill(Ntp->Electron_p4(probe).Pt(),Ntp->Electron_supercluster_eta(probe),w);
			  if(testprobe.Pt()>0 && testprobe.DeltaR(Ntp->Electron_p4(probe))<testdr){
				  numerator.Fill(Ntp->Electron_p4(probe).Pt(),Ntp->Electron_supercluster_eta(probe),w);
			  }
		  }
		  if(tagname=="electron"){
			  denominator.Fill(Ntp->Muons_p4(probe).Pt(),Ntp->Muons_p4(probe).Eta(),w);
			  if(testprobe.Pt()>0 && testprobe.DeltaR(Ntp->Muons_p4(probe))<testdr){
				  numerator.Fill(Ntp->Muons_p4(probe).Pt(),Ntp->Muons_p4(probe).Eta(),w);
			  }
		  }
	  }
	}
}

void ZtoEMu_Fakerate::triggerMatch(std::vector<unsigned int> tags, std::vector<unsigned int> probes, std::string tagtrigger, std::string probetrigger, std::string tagname, TH2D &denominator, TH2D &numerator,
		  TH1D &tagmupt, TH1D &tagmueta, TH1D &tagept, TH1D & tageeta, TH1D &probemupt, TH1D &probemueta, TH1D &probeept, TH1D &probeeeta,
		  TH1D &drtagmuprobee, TH1D &drtageprobemu, TH1D &ptbaltagmuprobee, TH1D &ptbaltageprobemu, TH1D &mttagmu, TH1D &mttage, TH1D &mtprobemu, TH1D &mtprobee, TH1D &mtagmuprobee, TH1D &mtageprobemu, double w){
	unsigned int finalprobe = 999;
	std::vector<unsigned int> Probe;
	bool flag;
	// match tag-particle to trigger for probe leg measurement
	if(Ntp->TriggerAccept(tagtrigger)){
	  flag = false;
	  // find muon that triggered
	  unsigned int tag = 999;
	  float testdr = 0.5;
	  TLorentzVector testtag(0.,0.,0.,0.);
	  for(unsigned i=0;i<tags.size();i++){
		  for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
			  if(Ntp->HLTTrigger_objs_trigger(j).find(tagtrigger) != string::npos){
				  for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
					  if(tagname=="muon"){
						  if(Ntp->HLTTrigger_objs_Id(j,k)==83){
							  testtag.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
								  Ntp->HLTTrigger_objs_Eta(j,k),
								  Ntp->HLTTrigger_objs_Phi(j,k),
								  Ntp->HLTTrigger_objs_E(j,k));
							  if(testtag.Pt()>0.
									  && testtag.DeltaR(Ntp->Muons_p4(tags.at(i)))<testdr
									  ){
								  tag = tags.at(i);
								  flag = true;
								  tagmupt.Fill(Ntp->Muons_p4(tag).Pt(),w);
								  tagmueta.Fill(Ntp->Muons_p4(tag).Eta(),w);
								  mttagmu.Fill(sqrt(2*Ntp->Muons_p4(tag).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muons_p4(tag).Px(),Ntp->Muons_p4(tag).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
							  }
						  }
					  }
					  if(tagname=="electron"){
						  if(Ntp->HLTTrigger_objs_Id(j,k)==82){
							  testtag.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
								  Ntp->HLTTrigger_objs_Eta(j,k),
								  Ntp->HLTTrigger_objs_Phi(j,k),
								  Ntp->HLTTrigger_objs_E(j,k));
							  if(testtag.Pt()>0.
									  && testtag.DeltaR(Ntp->Electron_p4(tags.at(i)))<testdr
									  ){
								  tag = tags.at(i);
								  flag = true;
								  tagept.Fill(Ntp->Electron_p4(tag).Pt(),w);
								  tageeta.Fill(Ntp->Electron_supercluster_eta(tag),w);
								  mttage.Fill(sqrt(2*Ntp->Electron_p4(tag).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(tag).Px(),Ntp->Electron_p4(tag).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
							  }
						  }
					  }
					  if(flag) break;
				  }
			  }
			  if(flag) break;
		  }
	  }
	  // create collection of probes
	  if(tag!=999){
		  if(tagname=="muon"){
			  for(unsigned i=0;i<probes.size();i++){
				  if(Ntp->Electron_Charge(probes.at(i))+Ntp->Muon_Charge(tag)==0
						  && (Ntp->Electron_p4(probes.at(i))+Ntp->Muons_p4(tag)).M()>30.
						  && (Ntp->Electron_p4(probes.at(i))+Ntp->Muons_p4(tag)).M()<90.
						  ){
							  Probe.push_back(probes.at(i));
						  }
			  }
		  }
		  if(tagname=="electron"){
			  for(unsigned i=0;i<probes.size();i++){
				  if(Ntp->Muon_Charge(probes.at(i))+Ntp->Electron_Charge(tag)==0
						  && (Ntp->Muons_p4(probes.at(i))+Ntp->Electron_p4(tag)).M()>30.
						  && (Ntp->Muons_p4(probes.at(i))+Ntp->Electron_p4(tag)).M()<90.
						  ){
					  Probe.push_back(probes.at(i));
				  }
			  }
		  }
	  }
	  // pick probe electron with highest pt
	  float probept = 8.;
	  unsigned int probe = 999;
	  for(unsigned i=0;i<Probe.size();i++){
		  if(tagname=="muon"){
			  if(Ntp->Electron_p4(Probe.at(i)).Pt()>probept){
				  probept = Ntp->Electron_p4(Probe.at(i)).Pt();
				  probe = Probe.at(i);
			  }
		  }
		  if(tagname=="electron"){
			  if(Ntp->Muons_p4(Probe.at(i)).Pt()>probept){
				  probept = Ntp->Muons_p4(Probe.at(i)).Pt();
				  probe = Probe.at(i);
			  }
		  }
	  }
	  // match probe particle to trigger of probe trigger
	  if(probe!=999){
		  TLorentzVector testprobe(0.,0.,0.,0.);
		  flag = false;
		  for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
			  if(Ntp->HLTTrigger_objs_trigger(i).find(probetrigger) != string::npos){
				  for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
					  if(tagname=="muon"){
						  if(Ntp->HLTTrigger_objs_Id(i,j)==82){
							  testprobe.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
									  Ntp->HLTTrigger_objs_Eta(i,j),
									  Ntp->HLTTrigger_objs_Phi(i,j),
									  Ntp->HLTTrigger_objs_E(i,j));
							  flag = true;
						  }
					  }
					  if(tagname=="electron"){
						  if(Ntp->HLTTrigger_objs_Id(i,j)==83){
							  testprobe.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
									  Ntp->HLTTrigger_objs_Eta(i,j),
									  Ntp->HLTTrigger_objs_Phi(i,j),
									  Ntp->HLTTrigger_objs_E(i,j));
							  flag = true;
						  }
					  }
					  if(flag) break;
				  }
			  }
			  if(flag) break;
		  }
		  if(tagname=="muon"){
			  denominator.Fill(Ntp->Electron_p4(probe).Pt(),Ntp->Electron_supercluster_eta(probe),w);
			  probeept.Fill(Ntp->Electron_p4(probe).Pt(),w);
			  probeeeta.Fill(Ntp->Electron_supercluster_eta(probe),w);
			  mtprobee.Fill(sqrt(2*Ntp->Electron_p4(probe).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(probe).Px(),Ntp->Electron_p4(probe).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
			  drtagmuprobee.Fill(Ntp->Muons_p4(tag).DeltaR(Ntp->Electron_p4(probe)),w);
			  ptbaltagmuprobee.Fill((Ntp->Muons_p4(tag)+Ntp->Electron_p4(probe)).Pt(),w);
			  mtagmuprobee.Fill((Ntp->Muons_p4(tag)+Ntp->Electron_p4(probe)).M(),w);
			  if(testprobe.Pt()>0. && testprobe.DeltaR(Ntp->Electron_p4(probe))<testdr){
				  numerator.Fill(Ntp->Electron_p4(probe).Pt(),Ntp->Electron_supercluster_eta(probe),w);
			  }
		  }
		  if(tagname=="electron"){
			  denominator.Fill(Ntp->Muons_p4(probe).Pt(),Ntp->Muons_p4(probe).Eta(),w);
			  probemupt.Fill(Ntp->Muons_p4(probe).Pt(),w);
			  probemueta.Fill(Ntp->Muons_p4(probe).Eta(),w);
			  mtprobemu.Fill(sqrt(2*Ntp->Muons_p4(probe).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muons_p4(probe).Px(),Ntp->Muons_p4(probe).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
			  drtageprobemu.Fill(Ntp->Electron_p4(tag).DeltaR(Ntp->Muons_p4(probe)),w);
			  ptbaltageprobemu.Fill((Ntp->Electron_p4(tag)+Ntp->Muons_p4(probe)).Pt(),w);
			  mtageprobemu.Fill((Ntp->Electron_p4(tag)+Ntp->Muons_p4(probe)).M(),w);
			  if(testprobe.Pt()>0. && testprobe.DeltaR(Ntp->Muons_p4(probe))<testdr){
				  numerator.Fill(Ntp->Muons_p4(probe).Pt(),Ntp->Muons_p4(probe).Eta(),w);
			  }
		  }
	  }
	}
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_Fakerate::isTightMuon(unsigned int i){
	if(!Ntp->Muon_isGlobalMuon(i)) return false;
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(Ntp->Muon_normChi2(i)>=10.) return false;
	if(Ntp->Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
	if(Ntp->Muon_numberOfMatchedStations(i)<=1) return false;
	if(Ntp->Muon_numberofValidPixelHits(i)<=0) return false;
	if(Ntp->Muon_trackerLayersWithMeasurement(i)<=5) return false;
	return true;
}

bool ZtoEMu_Fakerate::isTightMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.5) return false;
	return true;
}

bool ZtoEMu_Fakerate::isLooseMuon(unsigned int i){
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(!(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))) return false;
	return true;
}

bool ZtoEMu_Fakerate::isFakeMuon(unsigned int i){
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

bool ZtoEMu_Fakerate::isFakeMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeMuon(i)) return false;
	if(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
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

bool ZtoEMu_Fakerate::isMVANonTrigElectron(unsigned int i, unsigned int j){
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

bool ZtoEMu_Fakerate::isMVATrigElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isTightElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isTightElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightElectron(i)) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

bool ZtoEMu_Fakerate::isFakeElectron(unsigned int i){
	if(Ntp->Electron_p4(i).Pt()<=10) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
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

bool ZtoEMu_Fakerate::isFakeElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeElectron(i)) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	return true;
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
// Finish function
//

void ZtoEMu_Fakerate::Finish(){
	/*for(unsigned int i=1;i<=fakemu.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<=fakemu.at(0).GetNbinsY();j++){
			int xbin = getxbin(fakemu.at(0).GetXaxis()->GetBinLowEdge(i));
			int ybin = getybin(fakemu.at(0).GetYaxis()->GetBinLowEdge(j),"muon");
			fakemu_rebin.at(0).SetBinContent(xbin,ybin,fakemu.at(0).GetBinContent(i,j)+fakemu_rebin.at(0).GetBinContent(xbin,ybin));
			xbin = getxbin(tightmu.at(0).GetXaxis()->GetBinLowEdge(i));
			ybin = getybin(tightmu.at(0).GetYaxis()->GetBinLowEdge(j),"muon");
			tightmu_rebin.at(0).SetBinContent(xbin,ybin,tightmu.at(0).GetBinContent(i,j)+tightmu_rebin.at(0).GetBinContent(xbin,ybin));
			
			xbin = getxbin(fakee.at(0).GetXaxis()->GetBinLowEdge(i));
			ybin = getybin(fakee.at(0).GetYaxis()->GetBinLowEdge(j),"electron");
			fakee_rebin.at(0).SetBinContent(xbin,ybin,fakee.at(0).GetBinContent(i,j)+fakee_rebin.at(0).GetBinContent(xbin,ybin));
			xbin = getxbin(tighte.at(0).GetXaxis()->GetBinLowEdge(i));
			ybin = getybin(tighte.at(0).GetYaxis()->GetBinLowEdge(j),"electron");
			tighte_rebin.at(0).SetBinContent(xbin,ybin,tighte.at(0).GetBinContent(i,j)+tighte_rebin.at(0).GetBinContent(xbin,ybin));
		}
	}*/
	double muprob;
	double eprob;
	for(unsigned int i=1;i<=fakemu_rebin.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<=fakemu_rebin.at(0).GetNbinsY();j++){
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
	for(unsigned i=1;i<=muleg_numerator_rebin.at(0).GetNbinsX();i++){
		for(unsigned j=1;j<=muleg_numerator_rebin.at(0).GetNbinsY();j++){
			if(muleg_denominator_rebin.at(0).GetBinContent(i,j)!=0){
				muleg_eff.at(0).SetBinContent(i,j,muleg_numerator_rebin.at(0).GetBinContent(i,j)/muleg_denominator_rebin.at(0).GetBinContent(i,j));
			}
			if(muleg_denominator_rebin.at(11).GetBinContent(i,j)!=0){
				muleg_eff.at(11).SetBinContent(i,j,muleg_numerator_rebin.at(11).GetBinContent(i,j)/muleg_denominator_rebin.at(11).GetBinContent(i,j));
			}
			if(muleg_denominator_rebin.at(12).GetBinContent(i,j)!=0){
				muleg_eff.at(12).SetBinContent(i,j,muleg_numerator_rebin.at(12).GetBinContent(i,j)/muleg_denominator_rebin.at(12).GetBinContent(i,j));
			}
			if(muleg_denominator_rebin.at(13).GetBinContent(i,j)!=0){
				muleg_eff.at(13).SetBinContent(i,j,muleg_numerator_rebin.at(13).GetBinContent(i,j)/muleg_denominator_rebin.at(13).GetBinContent(i,j));
			}
		}
	}
	for(unsigned i=1;i<=eleg_numerator_rebin.at(0).GetNbinsX();i++){
		for(unsigned j=1;j<=eleg_numerator_rebin.at(0).GetNbinsY();j++){
			if(eleg_denominator_rebin.at(0).GetBinContent(i,j)!=0){
				eleg_eff.at(0).SetBinContent(i,j,eleg_numerator_rebin.at(0).GetBinContent(i,j)/eleg_denominator_rebin.at(0).GetBinContent(i,j));
			}
			if(eleg_denominator_rebin.at(11).GetBinContent(i,j)!=0){
				eleg_eff.at(11).SetBinContent(i,j,eleg_numerator_rebin.at(11).GetBinContent(i,j)/eleg_denominator_rebin.at(11).GetBinContent(i,j));
			}
			if(eleg_denominator_rebin.at(12).GetBinContent(i,j)!=0){
				eleg_eff.at(12).SetBinContent(i,j,eleg_numerator_rebin.at(12).GetBinContent(i,j)/eleg_denominator_rebin.at(12).GetBinContent(i,j));
			}
			if(eleg_denominator_rebin.at(13).GetBinContent(i,j)!=0){
				eleg_eff.at(13).SetBinContent(i,j,eleg_numerator_rebin.at(13).GetBinContent(i,j)/eleg_denominator_rebin.at(13).GetBinContent(i,j));
			}
		}
	}

	//mueff.at(0).Write();
	//eeff.at(0).Write();
	//outfile->Write();
	//outfile->Close();
	Selection::Finish();
}
