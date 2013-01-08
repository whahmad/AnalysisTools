#include "ZtoEMu.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

#include "Parameters.h"
#include "TMath.h"

ZtoEMu::ZtoEMu(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(24)
  ,e_pt(20)
  ,mu_eta(2.0)
  ,e_eta(2.0)
  ,jet_pt(30)
  ,jet_eta(2.4)
{
//    verbose=true;
}

ZtoEMu::~ZtoEMu(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu::~ZtoEMu Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu::~ZtoEMu()" << std::endl;
}

void  ZtoEMu::Configure(){

  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
//    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==MET)                cut.at(MET)=40;
    if(i==NMu)                cut.at(NMu)=1;
    if(i==NMuPt)              cut.at(NMuPt)=1;
    if(i==NMuEta)             cut.at(NMuEta)=1;
    if(i==MuIso)              cut.at(MuIso)=0.4;
    if(i==NE)                 cut.at(NE)=1;
    if(i==NEPt)               cut.at(NEPt)=1;
    if(i==NEEta)              cut.at(NEEta)=1;
    if(i==EIso)               cut.at(EIso)=0.1;
    if(i==deltaPhi)           cut.at(deltaPhi)=0.;
    if(i==cosdeltaPhi)        cut.at(cosdeltaPhi)=-0.2;
    if(i==charge)             cut.at(charge)=0;
    if(i==ptBalance)          cut.at(ptBalance)=15;
    if(i==ZMassmax)           cut.at(ZMassmax)=PDG_Var::Z_mass()+30;
    if(i==ZMassmin)           cut.at(ZMassmin)=80;//PDG_Var::Z_mass()-30;
	if(i==jetVeto)            cut.at(jetVeto)=0;
	if(i==Mt)                 cut.at(Mt)=500;
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
/*    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
    }
*/
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
      title.at(i)+=mu_eta;
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
      title.at(i)+=e_eta;
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
    else if(i==MET){
      title.at(i)="$E_{T}^{Miss} < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="E_{T}^{Miss} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==ptBalance){
      title.at(i)="$p_{t,e}-p_{t,\\mu} < $";
      title.at(i)+=cut.at(ptBalance);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="p_t balance (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ptBalance_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ptBalance_",htitle,40,0,200,hlabel,"Events"));
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
    else if(i==deltaPhi){
      title.at(i)="$\\phi_{e-\\mu} > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#phi_{e-#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,40,0.,2*TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,40,0.,2*TMath::Pi(),hlabel,"Events"));
    }
    else if(i==cosdeltaPhi){
      title.at(i)="$cos(\\phi_{e-\\mu}) < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(cosdeltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="cos(#phi_{e-#mu})";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_cosdeltaPhi_",htitle,40,-1,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_cosdeltaPhi_",htitle,40,-1,1,hlabel,"Events"));
    }
    else if(i==ZMassmax){
      title.at(i)="$M_{e,\\mu} < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmax));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{e,#mu} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,60,50,130,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,60,50,130,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_ZMassmax_","M_{e,#mu} (N-1 Distribution)",60,0,130,"M_{e,#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_ZMassmax_","M_{e,#mu} (Accumulative Distribution)",60,0,130,"M_{e,#mu} (GeV)","Events");
    }
    else if(i==ZMassmin){
      title.at(i)="$M_{e,\\mu} > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmin));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{e,#mu} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmin_",htitle,60,50,130,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmin_",htitle,60,50,130,hlabel,"Events"));
   } 
    else if(i==MuIso){
      title.at(i)="$dr(\\mu,Jet) > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(MuIso));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="dr(#mu,Jet) ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuIso_",htitle,40,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuIso_",htitle,40,0,20,hlabel,"Events"));
    }
    else if(i==EIso){
      title.at(i)="$dr(e,Jet) > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(EIso));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="dr(e,Jet) ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_EIso_",htitle,40,0,20,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_EIso_",htitle,40,0,20,hlabel,"Events"));
    }
	else if(i==jetVeto){
      title.at(i)="Number of Jets with high $P_{T} > $";
      title.at(i)+=jet_pt;
      title.at(i)+=" $=$";
      title.at(i)+=cut.at(jetVeto);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="high P_{T} jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_jetVeto_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_jetVeto_",htitle,40,0,200,hlabel,"Events"));
    }
	else if(i==Mt){
      title.at(i)="$m_{T}^{Miss} < $";
      title.at(i)+=cut.at(Mt);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="m_{T}^{Miss} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Mt_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Mt_",htitle,40,0,200,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  RelIsoE=HConfig.GetTH1D(Name+"_RelIsoE","RelIsoE",100,0.,10.,"Relative Isolation of Electron");
  RelIsoMu=HConfig.GetTH1D(Name+"_RelIsoMu","RelIsoMu",100,0.,10.,"Relative Isolation of Muon");
  EPt=HConfig.GetTH1D(Name+"_PtE","PtE",40,0.,200.,"Pt of electron");
  MuPt=HConfig.GetTH1D(Name+"_PtMu","PtMu",40,0.,200.,"Pt of muon");
  mtMu=HConfig.GetTH1D(Name+"_mtMu","mtMu",40,0.,200.,"Mt of muon and MET");
  phiMtMu=HConfig.GetTH1D(Name+"_phiMtMu","phiMtMu",40,0.,2*TMath::Pi(),"angle between muon and MET");
  mtE=HConfig.GetTH1D(Name+"_mtE","mtE",40,0.,200.,"Mt of electron and MET");
  phiMtE=HConfig.GetTH1D(Name+"_phiMtE","phiMtE",40,0.,2*TMath::Pi(),"angle between electron and MET");
  NJets=HConfig.GetTH1D(Name+"_NJets","NJets",20,0,20,"number of jets");
  pzeta=HConfig.GetTH1D(Name+"_pzeta","pzeta",40,-100.,100.,"pzeta");
  pzetaDQM=HConfig.GetTH1D(Name+"_pzetaDQM","pzetaDQM",40,-100.,100.,"pzetaDQM");
  //METvsMt=HConfig.GetTH2D(Name+"METvsMt","MET vs. m_t",40,0.,200.,40,0.,200.,"m_t","MET");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&RelIsoE);
 Extradist1d.push_back(&RelIsoMu);
 Extradist1d.push_back(&EPt);
 Extradist1d.push_back(&MuPt);
 Extradist1d.push_back(&mtMu);
 Extradist1d.push_back(&phiMtMu);
 Extradist1d.push_back(&mtE);
 Extradist1d.push_back(&phiMtE);
 Extradist1d.push_back(&NJets);
 Extradist1d.push_back(&pzeta);
 Extradist1d.push_back(&pzetaDQM);
 //Extradist2d.push_back(&METvsMt);

}

void  ZtoEMu::doEvent(){
  if(verbose)std::cout << "ZtoEMu::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
/*
  value.at(TriggerOk)=0;
  std::cout << "Number of HLT Triggers:" << Ntp->NHLTTriggers() << std::endl;
  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
	  std::cout << "######################################" << std::endl;
	  std::cout << "Trigger: " << Ntp->HTLTriggerName(i) << std::endl;
	  std::cout << "Acceptet (1) or not (0): " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << std::endl;
	  std::cout << "######################################" << std::endl;
  }
  if(Ntp->TriggerAccept("HLT_IsoMu24_eta2p1"))value.at(TriggerOk)=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);*/

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
  // Mu Cuts
  //
  if(verbose) std::cout << "Muon cuts" << std::endl;
  std::vector<unsigned int> GoodMuons;
  for(unsigned i=0;i<Ntp->NMuons();i++){
    if(Ntp->isGoodMuon(i))GoodMuons.push_back(i);
  }

  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuPt).push_back(Ntp->Muons_p4(GoodMuons.at(i)).Pt());
    if(Ntp->Muons_p4(GoodMuons.at(i)).Pt()<mu_pt){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  value.at(NMuPt)=GoodMuons.size();
  pass.at(NMuPt)=(value.at(NMuPt)>=cut.at(NMuPt));
  
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
  if(verbose)std::cout << "void  ZtoEMu::doEvent() E " << muidx1 <<" "<< muInfoidx1 << std::endl;

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  for(unsigned i=0;i<Ntp->NElectrons();i++){
    if(Ntp->isGoodElectron(i))GoodElectrons.push_back(i);
  }
  
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEPt).push_back(Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
    if(Ntp->Electron_p4(GoodElectrons.at(i)).Pt()<e_pt){
      GoodElectrons.erase(GoodElectrons.begin()+i);
      i--;
    }
  }
  value.at(NEPt)=GoodElectrons.size();
  pass.at(NEPt)=(value.at(NEPt)>=cut.at(NEPt));
  
  for(unsigned i=0;i<GoodElectrons.size();i++){
    dist.at(NEEta).push_back(Ntp->Electron_p4(GoodElectrons.at(i)).Eta());
    if(fabs(Ntp->Electron_p4(GoodElectrons.at(i)).Eta())>e_eta){
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
  if(verbose)std::cout << "void  ZtoEMu::doEvent() E " << eidx1 <<" "<< eInfoidx1 << std::endl;

  ///////////////////////////////////////////////
  //
  // Pt balance
  //
  if(verbose) std::cout << "pt balance" << std::endl;
  double balance(999);
  if(GoodMuons.size()==GoodElectrons.size()){
	  for(unsigned int i=0; i<GoodMuons.size(); i++){
		  balance = fabs(Ntp->Muons_p4(GoodMuons.at(i)).Pt()-Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
	  }
  }else if(GoodMuons.size()<GoodElectrons.size()){
	  for(unsigned int i=0; i<GoodMuons.size(); i++){
		  balance = fabs(Ntp->Muons_p4(GoodMuons.at(i)).Pt()-Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
	  }
  }else if(GoodMuons.size()>GoodElectrons.size()){
	  for(unsigned int i=0; i<GoodElectrons.size(); i++){
		  balance = fabs(Ntp->Muons_p4(GoodMuons.at(i)).Pt()-Ntp->Electron_p4(GoodElectrons.at(i)).Pt());
	  }
  }
  value.at(ptBalance)=balance;
  pass.at(ptBalance)=(value.at(ptBalance)<=cut.at(ptBalance));

  ///////////////////////////////////////////////
  //
  // Charge, deltaPhi, Isolation, mt
  //
  if(verbose) std::cout << "Charge, deltaPhi, Isolation, mt" << std::endl;
  value.at(charge)=-5;

  value.at(deltaPhi)=-1;
  value.at(cosdeltaPhi)=1;

  value.at(EIso)=5;
  value.at(MuIso)=5;
  value.at(ZMassmax)=8000;
  value.at(ZMassmin)=40;
  value.at(jetVeto)=0;
  int overlapE21(0), overlapE22(0);
  int overlapMu21(0), overlapMu22(0);
  int overlapE11(0), overlapMu11(0);
  double dphi(0);

  if(eInfoidx1!=999 && muInfoidx1!=999){
	value.at(charge)=Ntp->Electron_Charge(eInfoidx1)+Ntp->Muon_Charge(muInfoidx1);
  }
  if(eidx1!=999 && muidx1!=999){
	TLorentzVector E11=Ntp->Electron_p4(eidx1);
	TLorentzVector Mu11=Ntp->Muons_p4(muidx1);
	TLorentzVector ZEMu=E11+Mu11;
	dphi=Tools::DeltaPhi(E11.Phi(),Mu11.Phi());
	if(dphi<0)dphi+=2*TMath::Pi();
		value.at(deltaPhi)=dphi;
		value.at(cosdeltaPhi)=cos(Tools::DeltaPhi(E11.Phi(),Mu11.Phi()));
		value.at(ZMassmax)=ZEMu.M();
		value.at(ZMassmin)=ZEMu.M();
  }
  for(int i=0;i<Ntp->NPFJets();i++){
	double dREMu=5;
	if(eidx1!=999)dREMu=Tools::dr(Ntp->Electron_p4(eidx1),Ntp->PFJet_p4(i));
	if(dREMu<cut.at(EIso))overlapE11++;
	if(muidx1!=999)dREMu=Tools::dr(Ntp->Muons_p4(muidx1),Ntp->PFJet_p4(i));
	if(dREMu<cut.at(MuIso))overlapMu11++;
	if(sqrt(pow(Ntp->PFJet_Poca(i).X()-Ntp->Vtx(0).X(),2)+pow(Ntp->PFJet_Poca(i).Y()-Ntp->Vtx(0).Y(),2))<0.2 && fabs(Ntp->PFJet_Poca(i).Z()-Ntp->Vtx(0).Z())<0.5 && Ntp->PFJet_p4(i).Pt()>=jet_pt)value.at(jetVeto)=1;
  }
  
  value.at(Mt)=999;
  if(muidx1!=999){
    value.at(Mt)=sqrt(2*Ntp->Muons_p4(muidx1).Pt()*Ntp->MET_et()*(1-cos(Ntp->Muons_p4(muidx1).Phi()-Ntp->MET_phi())));
  }

  pass.at(jetVeto)=(value.at(jetVeto)<=cut.at(jetVeto));
  pass.at(charge)=(value.at(charge)==0);
  pass.at(deltaPhi)=(fabs(value.at(deltaPhi))>cut.at(deltaPhi));
  pass.at(cosdeltaPhi)=(value.at(cosdeltaPhi)<cut.at(cosdeltaPhi));
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));
  pass.at(ZMassmin)=(value.at(ZMassmin)>cut.at(ZMassmin));
  pass.at(EIso)=(1==overlapE11);
  pass.at(MuIso)=(1==overlapMu11);
  pass.at(Mt)=(value.at(Mt)<=cut.at(Mt));

  ///////////////////////////////////////////////
  //
  // Jet Cuts
  //
  if(verbose) std::cout << "Jet cuts" << std::endl;
  std::vector<unsigned int> GoodJets;
  for(int i=0;i<Ntp->NPFJets();i++){
    if(/*Ntp->isGoodJet(i) &&*/ Ntp->PFJet_p4(i).Pt()>jet_pt && fabs(Ntp->PFJet_p4(i).Eta())<jet_eta){
      bool overlap=false;
      for(unsigned j=0;j<GoodMuons.size();j++){
		if(Tools::dr(Ntp->Muons_p4(GoodMuons.at(j)),Ntp->PFJet_p4(i))<0.5) overlap=true;		
      }
	  for(unsigned j=0;j<GoodElectrons.size();j++){
		if(Tools::dr(Ntp->Electron_p4(GoodElectrons.at(j)),Ntp->PFJet_p4(i))<0.5) overlap=true;
	  }
      if(!overlap)GoodJets.push_back(i);
    }
  }
  //value.at(NJets)=GoodJets.size();
  //pass.at(NJets)=(value.at(NJets)>=cut.at(NJets));
  //pass.at(NJets)=true;
  ///////////////////////////////////////////////
  //
  // Event Shape and Energy Cuts
  //
  if(verbose) std::cout << "Event Shape and Energy Cuts" << std::endl;
  double MET_Ex(Ntp->MET_ex()),MET_Ey(Ntp->MET_ey());
  // correct for neutrino from mus
  value.at(MET)=sqrt(MET_Ex*MET_Ex+MET_Ey*MET_Ey);
  pass.at(MET)=(value.at(MET)<cut.at(MET));

  ///////////////////////////////////////////////////////////
  //Do QCD bkg
  /*if(!pass.at(charge)){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(charge)=true;
    }
    }*/
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData()){
    //w*=Ntp->EvtWeight3D();
    if(verbose)std::cout << "void  ZtoEMu::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);
  }
  if(verbose)std::cout << "filling NJets" << std::endl;
  //NJets.at(t).Fill(Ntp->NPFJets(),w);
  NJets.at(t).Fill(GoodJets.size(),w);
  if(verbose)std::cout << "filled NJtes" << std::endl;
  
  if(verbose)std::cout << "looping over good electrons" << std::endl;
  double phimte;
  for(unsigned int i=0;i<GoodElectrons.size();i++){
        EPt.at(t).Fill(Ntp->Electron_p4(GoodElectrons.at(i)).Pt(),w);
        mtE.at(t).Fill(sqrt(2*Ntp->Electron_p4(GoodElectrons.at(i)).Pt()*Ntp->MET_et()*(1-cos(Ntp->Electron_p4(GoodElectrons.at(i)).Phi()-Ntp->MET_phi()))),w);
        phimte=Ntp->Electron_p4(GoodElectrons.at(i)).Phi()-Ntp->MET_phi();
        if(phimte<0)phimte+=2*TMath::Pi();
        phiMtE.at(t).Fill(phimte,w);
        RelIsoE.at(t).Fill((Ntp->Electron_Gsf_dr03TkSumPt(GoodElectrons.at(i))+Ntp->Electron_Gsf_dr03HcalTowerSumEt(GoodElectrons.at(i))+Ntp->Electron_Gsf_dr03EcalRecHitSumE(GoodElectrons.at(i)))/Ntp->Electron_p4(GoodElectrons.at(i)).Pt(),w);
  }
  if(verbose)std::cout << "looping over good muons" << std::endl;
  double phimtmu;	
  for(unsigned int i=0;i<GoodMuons.size();i++){
      MuPt.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).Pt(),w);
      mtMu.at(t).Fill(sqrt(2*Ntp->Muons_p4(GoodMuons.at(i)).Pt()*Ntp->MET_et()*(1-cos(Ntp->Muons_p4(GoodMuons.at(i)).Phi()-Ntp->MET_phi()))),w);
      phimtmu=Ntp->Muons_p4(GoodMuons.at(i)).Phi()-Ntp->MET_phi();
      if(phimtmu<0)phimtmu+=2*TMath::Pi();
      phiMtMu.at(t).Fill(phimtmu,w);
      RelIsoMu.at(t).Fill((Ntp->Muon_sumPt03(GoodMuons.at(i))+Ntp->Muon_hadEt03(GoodMuons.at(i))+Ntp->Muon_emEt03(GoodMuons.at(i)))/Ntp->Muons_p4(GoodMuons.at(i)).Pt(),w);
  }
  
  if(verbose)std::cout << "Calculating pzeta" << std::endl;
  if(muidx1!=999 && eidx1!=999){
	  if(GoodMuons.size()==GoodElectrons.size()){
	  	  for(unsigned i=0;i<GoodMuons.size();i++){
	  		  pzeta.at(t).Fill(calculatePzeta(i,GoodElectrons,GoodMuons),w);
	  		  pzetaDQM.at(t).Fill(calculatePzetaDQM(i,GoodElectrons,GoodMuons),w);
	  	  }
	  }
	  else if(GoodMuons.size()<GoodElectrons.size()){
	  	  for(unsigned i=0;i<GoodMuons.size();i++){
	  		  pzeta.at(t).Fill(calculatePzeta(i,GoodElectrons,GoodMuons),w);
	  		  pzetaDQM.at(t).Fill(calculatePzetaDQM(i,GoodElectrons,GoodMuons),w);
	  	  }
	  }
	  else if(GoodMuons.size()>GoodElectrons.size()){
	  	  for(unsigned i=0;i<GoodElectrons.size();i++){
	  		  pzeta.at(t).Fill(calculatePzeta(i,GoodElectrons,GoodMuons),w);
	  		  pzetaDQM.at(t).Fill(calculatePzetaDQM(i,GoodElectrons,GoodMuons),w);
	  	  }
	  }
  }

  if(verbose)std::cout << "ZtoEMu::doEvent() doEvent END" << std::endl;
}

double ZtoEMu::calculatePzeta(int iterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2){
  pex=Ntp->Electron_p4(vec1.at(iterator)).Px();
  pey=Ntp->Electron_p4(vec1.at(iterator)).Py();
  pmux=Ntp->Muons_p4(vec2.at(iterator)).Px();
  pmuy=Ntp->Muons_p4(vec2.at(iterator)).Py();
  phie=Ntp->Electron_p4(vec1.at(iterator)).Phi();
  phimu=Ntp->Muons_p4(vec2.at(iterator)).Phi();
  combpt=TMath::Sqrt(pow(pex+pmux,2)+pow(pey+pmuy,2));
  aemu=TMath::ACos(pmux*pex+pmuy*pey/(Ntp->Muons_p4(vec2.at(iterator)).Pt()*Ntp->Electron_p4(vec1.at(iterator)).Pt()));
  if(phie<phimu && fabs(phie-phimu)<TMath::Pi())phi1=phie;
  else if(phimu<phie && fabs(phie-phimu)>TMath::Pi())phi1=phie;
  else if(phie<phimu && fabs(phie-phimu)>TMath::Pi())phi1=phimu;
  else if(phimu<phie && fabs(phie-phimu)<TMath::Pi())phi1=phimu;
  beta=TMath::ACos(((pex+pmux)*TMath::Cos(phi1+0.5*aemu)+(pey+pmuy)*TMath::Sin(phi1+0.5*aemu))/combpt);
  gamma=TMath::ACos((Ntp->MET_ex()*TMath::Cos(phi1+0.5*aemu)+Ntp->MET_ey()*TMath::Sin(phi1+0.5*aemu))/Ntp->MET_et());
  if(Ntp->MET_phi()>(phi1+0.5*aemu+0.5*TMath::Pi()) && Ntp->MET_phi()<(phi1+0.5*aemu+1.5*TMath::Pi()))gamma*=-1;
  pvis=TMath::Sin(beta)*combpt;
  pmiss=TMath::Sin(gamma)*Ntp->MET_et();
  return pmiss-pvis;
}

double ZtoEMu::calculatePzetaDQM(int iterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2){
	double cosPhi1 = TMath::Cos(Ntp->Electron_p4(vec1.at(iterator)).Phi());
	double sinPhi1 = TMath::Sin(Ntp->Electron_p4(vec1.at(iterator)).Phi());
	double cosPhi2 = TMath::Cos(Ntp->Muons_p4(vec2.at(iterator)).Phi());
	double sinPhi2 = TMath::Sin(Ntp->Muons_p4(vec2.at(iterator)).Phi());
	double zetaX = cosPhi1 + cosPhi2;
	double zetaY = sinPhi1 + sinPhi2;
	double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
	if(zetaR>0.){
		zetaX/=zetaR;
		zetaY/=zetaR;
	}
	double pxVis=Ntp->Electron_p4(vec1.at(iterator)).Px()+Ntp->Muons_p4(vec2.at(iterator)).Px();
	double pyVis=Ntp->Electron_p4(vec1.at(iterator)).Py()+Ntp->Muons_p4(vec2.at(iterator)).Py();
	double pZetaVis=pxVis*zetaX+pyVis*zetaY;
	double px=pxVis+Ntp->MET_ex();
	double py=pyVis+Ntp->MET_ey();
	double pZeta=px*zetaX+py*zetaY;
	return pZeta-2*pZetaVis;
}

