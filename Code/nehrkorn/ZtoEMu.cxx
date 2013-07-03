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

#include <TFile.h>

ZtoEMu::ZtoEMu(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(30)
  ,e_pt(30)
  ,mu_eta(2.1)
  ,e_eta(2.3)
  ,jet_pt(30)
  ,jet_eta(2.4)
  ,jet_sum(70)
  ,zmin(88)
  ,zmax(94)
  ,phimin(0.8)
  ,phimax(1.2)
  ,mtmu(40)
  ,ptbalance(20)
  ,dRmue(0.2)
{
    //verbose=true;
	MVA_ID = true;
	FRFile = new TFile("/net/scratch_cms/institut_3b/nehrkorn/FakeRates_2012_19ifb_rereco.root");
	ElectronFakeRate = (TH2D*)(FRFile->Get("ElectronFakeRateHist"));
	MuonFakeRate = (TH2D*)(FRFile->Get("MuonFakeRateHist"));
	twod = false;

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
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==NMu)                cut.at(NMu)=1;
    if(i==NMuPt)              cut.at(NMuPt)=1;
    if(i==NMuEta)             cut.at(NMuEta)=1;
    if(i==diMuonVeto)         cut.at(diMuonVeto)=0;
    if(i==triLeptonVeto)      cut.at(triLeptonVeto)=0;
    if(i==looseMuonVeto)      cut.at(looseMuonVeto)=0;
    if(i==NE)                 cut.at(NE)=1;
    if(i==NEPt)               cut.at(NEPt)=1;
    if(i==NEEta)              cut.at(NEEta)=1;
    if(i==SameVtx)            cut.at(SameVtx)=true;
    if(i==charge)             cut.at(charge)=0;
    if(i==ptBalance)          cut.at(ptBalance)=ptbalance;
    if(i==ZMassmax)           cut.at(ZMassmax)=zmax;
    if(i==ZMassmin)           cut.at(ZMassmin)=zmin;
	if(i==jetVeto)            cut.at(jetVeto)=jet_sum;
	if(i==MtMu)               cut.at(MtMu)=mtmu;
	if(i==drMuE)              cut.at(drMuE)=dRmue;
	if(i==qualitycuts)        cut.at(qualitycuts)=true;
	if(i==Phimin)             cut.at(Phimin)=phimin;
	if(i==Phimax)             cut.at(Phimax)=phimax;
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
    else if(i==ptBalance){
      title.at(i)="$p_{t,e+\\mu} < $";
      title.at(i)+=cut.at(ptBalance);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="p_t balance / GeV";
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
    else if(i==ZMassmax){
      title.at(i)="$M_{e,\\mu} < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmax));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{e,#mu} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,20,60,120,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,20,60,120,hlabel,"Events"));
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
      hlabel="M_{e,#mu} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmin_",htitle,20,60,120,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmin_",htitle,20,60,120,hlabel,"Events"));
    } 
	else if(i==looseMuonVeto){
		title.at(i)="Number of loose $\\mu = $";
		title.at(i)+=cut.at(looseMuonVeto);
		title.at(i)+="";
		htitle=title.at(i);
		htitle.ReplaceAll("$","");
		htitle.ReplaceAll("\\","#");
		hlabel="loose muon veto";
		Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_looseMuonVeto_",htitle,40,0,20,hlabel,"Events"));
		Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_looseMuonVeto_",htitle,40,0,20,hlabel,"Events"));
	}
	else if(i==jetVeto){
      title.at(i)="$P_{t}$ sum of 2 highest $P_{t}$-jets$ < $";
      title.at(i)+=cut.at(jetVeto);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="jet P_{T} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_jetVeto_",htitle,39,2,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_jetVeto_",htitle,39,2,200,hlabel,"Events"));
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
    else if(i==Phimin){
  	  title.at(i)="$\\varphi_{e,\\mu} / \\pi > $";
  	  char buffer[50];
  	  sprintf(buffer,"%5.2f",cut.at(Phimin));
  	  title.at(i)+=buffer;
  	  htitle=title.at(i);
  	  htitle.ReplaceAll("$","");
  	  htitle.ReplaceAll("\\","#");
  	  hlabel="#varphi_{e,#mu}";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Phimin_",htitle,40,0.,2.,hlabel,"Events"));
  	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Phimin_",htitle,40,0.,2.,hlabel,"Events"));
  	}
    else if(i==Phimax){
	  title.at(i)="$\\varphi_{e,\\mu} / \\pi < $";
	  char buffer[50];
	  sprintf(buffer,"%5.2f",cut.at(Phimax));
	  title.at(i)+=buffer;
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="#varphi_{e,#mu}";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Phimax_",htitle,40,0.,2.,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Phimax_",htitle,40,0.,2.,hlabel,"Events"));
	}
    else if(i==diMuonVeto){
	  title.at(i)="Number of $\\mu$ with $dR(\\mu,e)<=0.3 = $";
	  title.at(i)+=cut.at(diMuonVeto);
	  title.at(i)+="";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="dimuon veto";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_diMuonVeto_",htitle,40,0,20,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_diMuonVeto_",htitle,40,0,20,hlabel,"Events"));
	}
    else if(i==triLeptonVeto){
  	  title.at(i)="Number of additional $\\mu$ or e in event $ = $";
  	  title.at(i)+=cut.at(triLeptonVeto);
  	  title.at(i)+="";
  	  htitle=title.at(i);
  	  htitle.ReplaceAll("$","");
  	  htitle.ReplaceAll("\\","#");
  	  hlabel="trilepton veto";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_triLeptonVeto_",htitle,40,0,20,hlabel,"Events"));
  	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_triLeptonVeto_",htitle,40,0,20,hlabel,"Events"));
  	}
    
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  RelIsoE=HConfig.GetTH1D(Name+"_RelIsoE","RelIsoE",20,0.,1.,"Relative Isolation of Electron");
  RelIsoMu=HConfig.GetTH1D(Name+"_RelIsoMu","RelIsoMu",20,0.,1.,"Relative Isolation of Muon");
  EPt=HConfig.GetTH1D(Name+"_PtE","PtE",40,0.,200.,"P_{t}^{e} / GeV");
  MuPt=HConfig.GetTH1D(Name+"_PtMu","PtMu",40,0.,200.,"P_{t}^{#mu} / GeV");
  mtMu=HConfig.GetTH1D(Name+"_mtMu","mtMu",40,0.,200.,"m_{t}^{#mu} / GeV");
  mtE=HConfig.GetTH1D(Name+"_mtE","mtE",40,0.,200.,"m_{t}^{e} / GeV");
  NJets=HConfig.GetTH1D(Name+"_NJets","NJets",20,0,20,"number of jets");
  pzeta=HConfig.GetTH1D(Name+"_pzeta","pzeta",40,-100.,100.,"pzeta");
  pzetaDQM=HConfig.GetTH1D(Name+"_pzetaDQM","pzetaDQM",40,-100.,100.,"pzetaDQM");
  invmass=HConfig.GetTH1D(Name+"_invmass","invmass",40,0.,200.,"m_{e,#mu} / GeV");
  etaMu=HConfig.GetTH1D(Name+"_etaMu","etaMu",40,-3.5,3.5,"#eta_{#mu}");
  etaE=HConfig.GetTH1D(Name+"_etaE","etaE",40,-3.5,3.5,"#eta_{e}");
  jetsum=HConfig.GetTH1D(Name+"_jetsum","jetsum",80,0.,400,"P_{t}^{jets} / GeV");
  chargesum=HConfig.GetTH1D(Name+"_chargesum","chargesum",21,-5.5,5.5,"charge sum");
  chargesumID=HConfig.GetTH1D(Name+"_chargesumID","chargesumID",21,-5.5,5.5,"charge sum of tight e & tight #mu");
  drmue=HConfig.GetTH1D(Name+"_drmue","drmue",20,0.,1.,"dR(e,#mu)");
  drmueID=HConfig.GetTH1D(Name+"_drmueID","drmueID",20,0.,1.,"dR(e,#mu) (tight)");
  deltaphi=HConfig.GetTH1D(Name+"_deltaphi","deltaphi",40,0.,2,"#phi_{e,#mu}");
  deltaphiID=HConfig.GetTH1D(Name+"_deltaphiID","deltaphiID",40,0.,2,"#phi_{e,#mu} tight");
  ptbal=HConfig.GetTH1D(Name+"_ptbal","ptbal",40,0.,200.,"p_{t}^{e+#mu} / GeV");
  chargesumsigned=HConfig.GetTH1D(Name+"_chargesumsigned","chargesumsigned",21,-5.5,5.5,"charge sum");
  chargesumIDsigned=HConfig.GetTH1D(Name+"_chargesumIDsigned","chargesumIDsigned",21,-5.5,5.5,"charge sum of tight e & tight #mu");
  FirstJetPt=HConfig.GetTH1D(Name+"_FirstJetPt","FirstJetPt",40,0.,200.,"P_{t}^{1st jet} / GeV");
  SecondJetPt=HConfig.GetTH1D(Name+"_SecondJetPt","SecondJetPt",40,0.,200.,"P_{t}^{2nd jet} / GeV");
  ThirdJetPt=HConfig.GetTH1D(Name+"_ThirdJetPt","ThirdJetPt",40,0.,200.,"P_{t}^{3rd jet} / GeV");
  FourthJetPt=HConfig.GetTH1D(Name+"_FourthJetPt","FourthJetPt",40,0.,200.,"P_{t}^{4th jet} / GeV");
  
  invmass_zmass=HConfig.GetTH1D(Name+"_invmass_zmass","invmass_zmass",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_ptbalance=HConfig.GetTH1D(Name+"_invmass_ptbalance","invmass_ptbalance",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_mtmu=HConfig.GetTH1D(Name+"_invmass_mtmu","invmass_mtmu",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_jetveto=HConfig.GetTH1D(Name+"_invmass_jetveto","invmass_jetveto",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_charge=HConfig.GetTH1D(Name+"_invmass_charge","invmass_charge",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_loosemuonveto=HConfig.GetTH1D(Name+"_invmass_loosemuonveto","invmass_loosemuonveto",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_dremu=HConfig.GetTH1D(Name+"_invmass_dremu","invmass_dremu",20,60.,120.,"m_{e,#mu} / GeV");
  invmass_only_object_id=HConfig.GetTH1D(Name+"_invmass_only_object_id","invmass_only_object_id",20,60.,120.,"m_{e,#mu} / GeV");
  
  nm2_charge=HConfig.GetTH1D(Name+"_nm2_charge","nm2_charge",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_jetveto=HConfig.GetTH1D(Name+"_nm2_jetveto","nm2_jetveto",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_mtmu=HConfig.GetTH1D(Name+"_nm2_mtmu","nm2_mtmu",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_ptbalance=HConfig.GetTH1D(Name+"_nm2_ptbalance","nm2_ptbalance",20,60.,120.,"m_{e,#mu} / GeV");
  nm2_drmue=HConfig.GetTH1D(Name+"_nm2_drmue","nm2_drmue",20,60.,120.,"m_{e,#mu} / GeV");

  no_jetveto_mtmu=HConfig.GetTH1D(Name+"_no_jetveto_mtmu","no_jetveto_mtmu",20,60.,120.,"m_{e,#mu} / GeV");
  no_charge_ptbalance=HConfig.GetTH1D(Name+"_no_charge_ptbalance","no_charge_ptbalance",20,60.,120.,"m_{e,#mu} / GeV");
  no_charge=HConfig.GetTH1D(Name+"_no_charge","no_charge",20,60.,120.,"m_{e,#mu} / GeV");
  
  drmue_plus_charge=HConfig.GetTH1D(Name+"_drmue_plus_charge","drmue_plus_charge",20,60.,120.,"m_{e,#mu} / GeV");
  drmue_plus_jetveto=HConfig.GetTH1D(Name+"_drmue_plus_jetveto","drmue_plus_jetveto",20,60.,120.,"m_{e,#mu} / GeV");
  drmue_plus_mtmu=HConfig.GetTH1D(Name+"_drmue_plus_mtmu","drmue_plus_mtmu",20,60.,120.,"m_{e,#mu} / GeV");
  drmue_plus_ptbalance=HConfig.GetTH1D(Name+"_drmue_plus_ptbalance","drmue_plus_ptbalance",20,60.,120.,"m_{e,#mu} / GeV");
  
  mproj1=HConfig.GetTH1D(Name+"_mproj1","mproj1",7,60.,81.,"Upper left projected to mass");
  mproj2=HConfig.GetTH1D(Name+"_mproj2","mproj2",7,81.,102.,"Upper middle projected to mass");
  mproj3=HConfig.GetTH1D(Name+"_mproj3","mproj3",6,102.,120.,"Upper right projected to mass");
  mproj4=HConfig.GetTH1D(Name+"_mproj4","mproj4",7,60.,81.,"Middle left projected to mass");
  mproj5=HConfig.GetTH1D(Name+"_mproj5","mrpoj5",7,81.,102.,"Middle middle projected to mass");
  mproj6=HConfig.GetTH1D(Name+"_mproj6","mproj6",6,102.,120.,"Middle right projected to mass");
  mproj7=HConfig.GetTH1D(Name+"_mproj7","mproj7",7,60.,81.,"Lower left projected to mass");
  mproj8=HConfig.GetTH1D(Name+"_mproj8","mproj8",7,81.,102.,"Lower middle projected to mass");
  mproj9=HConfig.GetTH1D(Name+"_mproj9","mproj9",6,102.,120.,"Lower right projected to mass");
  phiproj1=HConfig.GetTH1D(Name+"_phiproj1","phiproj1",16,0.,2.5,"Upper left projected to phi");
  phiproj2=HConfig.GetTH1D(Name+"_phiproj2","phiproj2",16,0.,2.5,"Upper middle projected to phi");
  phiproj3=HConfig.GetTH1D(Name+"_phiproj3","phiproj3",16,0.,2.5,"Upper right projected to phi");
  phiproj4=HConfig.GetTH1D(Name+"_phiproj4","phiproj4",8,2.5,3.7,"Middle left projected to phi");
  phiproj5=HConfig.GetTH1D(Name+"_phiproj5","phiproj5",8,2.5,3.7,"Middle middle projected to phi");
  phiproj6=HConfig.GetTH1D(Name+"_phiproj6","phiproj6",8,2.5,3.7,"Middle right projected to phi");
  phiproj7=HConfig.GetTH1D(Name+"_phiproj7","phiproj7",17,3.7,2*TMath::Pi(),"Lower left projected to phi");
  phiproj8=HConfig.GetTH1D(Name+"_phiproj8","phiproj8",17,3.7,2*TMath::Pi(),"Lower middle projected to phi");
  phiproj9=HConfig.GetTH1D(Name+"_phiproj9","phiproj9",17,3.7,2*TMath::Pi(),"Lower right projected to phi");
  
  m1=HConfig.GetTH1D(Name+"_m1","m1",20,60.,120.,"top slice projected to mass");
  m2=HConfig.GetTH1D(Name+"_m2","m2",20,60.,120.,"middle slice projected to mass");
  m3=HConfig.GetTH1D(Name+"_m3","m3",20,60.,120.,"bottom slice projected to mass");
  phi1=HConfig.GetTH1D(Name+"_phi1","phi1",40,0.,2,"left slice projected to phi");
  phi2=HConfig.GetTH1D(Name+"_phi2","phi2",40,0.,2,"middle slice projected to phi");
  phi3=HConfig.GetTH1D(Name+"_phi3","phi3",40,0.,2,"right slice projected to phi");
  
  phi1_nopt=HConfig.GetTH1D(Name+"_phi1_nopt","phi1_nopt",40,0.,2,"left slice projected to phi");
  phi2_nopt=HConfig.GetTH1D(Name+"_phi2_nopt","phi2_nopt",40,0.,2,"middle slice projected to phi");
  phi3_nopt=HConfig.GetTH1D(Name+"_phi3_nopt","phi3_nopt",40,0.,2,"right slice projected to phi");
  
  phi1_ptdiff=HConfig.GetTH1D(Name+"_phi1_ptdiff","phi1_ptdiff",40,0.,2,"left slice projected to phi");
  phi2_ptdiff=HConfig.GetTH1D(Name+"_phi2_ptdiff","phi2_ptdiff",40,0.,2,"middle slice projected to phi");
  phi3_ptdiff=HConfig.GetTH1D(Name+"_phi3_ptdiff","phi3_ptdiff",40,0.,2,"right slice projected to phi");
  
  ptbal2=HConfig.GetTH1D(Name+"_ptbal2","ptbal2",40,0.,200.,"|p_{t}^{e} - p_{t}^{#mu}| / GeV");
  
  frMu=HConfig.GetTH1D(Name+"_frMu","frMu",100,0.,1.,"muon fakerate");
  frE=HConfig.GetTH1D(Name+"_frE","frE",100,0.,1.,"electron fakerate");
  
  pzetaCut=HConfig.GetTH1D(Name+"_pzetaCut","pzetaCut",40,-100.,100.,"pzetaCut");
  pzetaStatus=HConfig.GetTH1D(Name+"_pzetaStatus","pzetaStatus",40,-100.,100.,"pzetaStatus");
  metCut=HConfig.GetTH1D(Name+"_metCut","metCut",40,0.,200.,"MET / GeV");
  metStatus=HConfig.GetTH1D(Name+"_metStatus","metStatus",40,0.,200.,"MET / GeV");

  InvmassVsDeltaPhi=HConfig.GetTH2D(Name+"_InvmassVsDeltaPhi","m_{e,#mu} vs. #Delta#phi_{e,#mu}",20,60,120,40,0.,2.,"m_{e,#mu} / GeV","#Delta#phi_{e,#mu} / rad");
  PtDiffVsDeltaPhi=HConfig.GetTH2D(Name+"_PtDiffVsDeltaPhi","|p_{t}^{e} - p_{t}^{#mu}| vs. #Delta#phi_{e,#mu}",40,0.,200.,40,0.,2.,"|p_{t}^{e} - p_{t}^{#mu}| / GeV","#Delta#phi_{e,#mu} / rad");

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
 Extradist1d.push_back(&mtE);
 Extradist1d.push_back(&NJets);
 Extradist1d.push_back(&pzeta);
 Extradist1d.push_back(&pzetaDQM);
 Extradist1d.push_back(&invmass);
 Extradist1d.push_back(&jetsum);
 Extradist1d.push_back(&chargesum);
 Extradist1d.push_back(&chargesumID);
 Extradist1d.push_back(&drmue);
 Extradist1d.push_back(&drmueID);
 Extradist1d.push_back(&deltaphi);
 Extradist1d.push_back(&deltaphiID);
 Extradist1d.push_back(&ptbal);
 Extradist1d.push_back(&chargesumsigned);
 Extradist1d.push_back(&chargesumIDsigned);
 Extradist1d.push_back(&FirstJetPt);
 Extradist1d.push_back(&SecondJetPt);
 Extradist1d.push_back(&ThirdJetPt);
 Extradist1d.push_back(&FourthJetPt);
 
 Extradist1d.push_back(&invmass_zmass);
 Extradist1d.push_back(&invmass_ptbalance);
 Extradist1d.push_back(&invmass_mtmu);
 Extradist1d.push_back(&invmass_jetveto);
 Extradist1d.push_back(&invmass_charge);
 Extradist1d.push_back(&invmass_loosemuonveto);
 Extradist1d.push_back(&invmass_dremu);
 Extradist1d.push_back(&invmass_only_object_id);
 
 Extradist1d.push_back(&nm2_charge);
 Extradist1d.push_back(&nm2_jetveto);
 Extradist1d.push_back(&nm2_mtmu);
 Extradist1d.push_back(&nm2_ptbalance);
 Extradist1d.push_back(&nm2_drmue);

 Extradist1d.push_back(&no_jetveto_mtmu);
 Extradist1d.push_back(&no_charge_ptbalance);
 Extradist1d.push_back(&no_charge);
 
 Extradist1d.push_back(&drmue_plus_charge);
 Extradist1d.push_back(&drmue_plus_jetveto);
 Extradist1d.push_back(&drmue_plus_mtmu);
 Extradist1d.push_back(&drmue_plus_ptbalance);
 
 Extradist1d.push_back(&mproj1);
 Extradist1d.push_back(&mproj2);
 Extradist1d.push_back(&mproj3);
 Extradist1d.push_back(&mproj4);
 Extradist1d.push_back(&mproj5);
 Extradist1d.push_back(&mproj6);
 Extradist1d.push_back(&mproj7);
 Extradist1d.push_back(&mproj8);
 Extradist1d.push_back(&mproj9);
 Extradist1d.push_back(&phiproj1);
 Extradist1d.push_back(&phiproj2);
 Extradist1d.push_back(&phiproj3);
 Extradist1d.push_back(&phiproj4);
 Extradist1d.push_back(&phiproj5);
 Extradist1d.push_back(&phiproj6);
 Extradist1d.push_back(&phiproj7);
 Extradist1d.push_back(&phiproj8);
 Extradist1d.push_back(&phiproj9);
 
 Extradist1d.push_back(&m1);
 Extradist1d.push_back(&m2);
 Extradist1d.push_back(&m3);
 Extradist1d.push_back(&phi1);
 Extradist1d.push_back(&phi2);
 Extradist1d.push_back(&phi3);
 
 Extradist1d.push_back(&phi1_nopt);
 Extradist1d.push_back(&phi2_nopt);
 Extradist1d.push_back(&phi3_nopt);
 
 Extradist1d.push_back(&phi1_ptdiff);
 Extradist1d.push_back(&phi2_ptdiff);
 Extradist1d.push_back(&phi3_ptdiff);
 
 Extradist1d.push_back(&ptbal2);
 
 Extradist1d.push_back(&frMu);
 Extradist1d.push_back(&frE);
 
 Extradist1d.push_back(&pzetaCut);
 Extradist1d.push_back(&pzetaStatus);
 Extradist1d.push_back(&metCut);
 Extradist1d.push_back(&metStatus);
 
 Extradist2d.push_back(&InvmassVsDeltaPhi);
 Extradist2d.push_back(&PtDiffVsDeltaPhi);

}

void  ZtoEMu::doEvent(){
  if(verbose)std::cout << "ZtoEMu::doEvent() START" << std::endl;
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
  // Blinding
  //
  bool blinddataonly = true;
  bool blindingactive = false;
  int blind = 1;
  
  // if combination of muon & electron has mass within +-5 Z-widths, data (and or mc) will be weighted by 0
  if(pass.at(qualitycuts)){
	  double mass = (Ntp->Muons_p4(qualitymuons.at(0))+Ntp->Electron_p4(qualityelectrons.at(0))).M();
	  double delphi = Ntp->Muons_p4(qualitymuons.at(0)).DeltaPhi(Ntp->Electron_p4(qualityelectrons.at(0)))/TMath::Pi();
	  if(delphi<0)delphi+=2;
	  if(blinddataonly){
		  if(twod){
			  if(Ntp->isData() && mass>=81. && mass<=102. && delphi>=0.8 && delphi<=1.2){
				  blindingactive = true;
				  blind = 0;
			  }
		  }else{
			  if(Ntp->isData() && mass>=81. && mass<=102.){
				  blindingactive = true;
				  blind = 0;
			  }
		  }
	  }else{
		  if(twod){
			  if(mass>=81. && mass<=102. && delphi>=0.8 && delphi<=1.2){
				  blindingactive = true;
				  blind = 0;
			  }
		  }else{
			  if(mass>=81. && mass<=102.){
				  blindingactive = true;
				  blind = 0;
			  }
		  }
		}
  }
  
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
  if(verbose)std::cout << "void  ZtoEMu::doEvent() E " << muidx1 <<" "<< muInfoidx1 << std::endl;

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
  if(verbose)std::cout << "void  ZtoEMu::doEvent() E " << eidx1 <<" "<< eInfoidx1 << std::endl;
  
  ///////////////////////////////////////////////
  //
  // dR cleaning of e & mu
  //
  if(verbose)std::cout << "dR cleaning" << std::endl;
  
  value.at(drMuE)=10.;
  if(muidx1!=999 && eidx1!=999){
	  if(Ntp->Muon_Track_idx(GoodMuons.at(0))==Ntp->Electron_Track_idx(GoodElectrons.at(0))){
		  value.at(drMuE)=0.;
	  }else{
		  value.at(drMuE)=Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0)));
	  }
  }
  pass.at(drMuE)=(value.at(drMuE)>cut.at(drMuE));
  
  ///////////////////////////////////////////////
  //
  // Loose Muon Veto
  //
  if(verbose)std::cout << "Loose muon veto" << std::endl;
  int LooseMuons = 0;
  
  // loop over muons and check whether there are loose muons coming from selected vertex
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(isLooseMuon(i)
			  && !isFakeMuon(i,posmu)
			  && !isTightMuon(i,posmu)
			  && dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(posmu))<0.02
			  && dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(posmu))<0.2
			  ){
		  LooseMuons++;
	  }
  }
  value.at(looseMuonVeto)=LooseMuons;
  pass.at(looseMuonVeto)=(value.at(looseMuonVeto)==cut.at(looseMuonVeto));
  
  ///////////////////////////////////////////////
  //
  // Di Muon Veto
  //
  if(verbose)std::cout << "dimuon veto" << std::endl;
  int dimu(0);
  if(eidx1!=999 && muidx1!=999){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(i!=GoodMuons.at(0)
				  && Ntp->Muons_p4(i).Pt()>3
				  && fabs(Ntp->Muons_p4(i).Eta())<2.4
				  && Ntp->Electron_p4(GoodElectrons.at(0)).DeltaR(Ntp->Muons_p4(i))<=0.3
				  ){
			  dimu++;
		  }
	  }
  }
  value.at(diMuonVeto)=dimu;
  pass.at(diMuonVeto)=(value.at(diMuonVeto)==cut.at(diMuonVeto));
  
  ///////////////////////////////////////////////
  //
  // Tri Lepton Veto
  //
  if(verbose)std::cout << "trilepton veto" << std::endl;
  int trilep(0);
  if(eidx1!=999 && muidx1!=999){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(i!=GoodMuons.at(0)
				  && Ntp->Muons_p4(i).Pt()>10
				  && fabs(Ntp->Muons_p4(i).Eta())<2.4
				  && isTightMuon(i)
				  && dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(posmu))<0.045
				  && dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(posmu))<0.2
				  && Muon_RelIso(i)<0.3
				  ){
			  trilep++;
		  }
	  }
	  for(unsigned i=0;i<Ntp->NElectrons();i++){
		  if(i!=GoodElectrons.at(0)
				  && Ntp->Electron_p4(i).Pt()>10
				  && fabs(Ntp->Electron_supercluster_eta(i))<2.5
				  && isLooseElectron(i)
				  && dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(pose))<0.045
				  && dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(pose))<0.2
				  && Electron_RelIso(i)<0.3
				  ){
			  trilep++;
		  }
	  }
  }
  value.at(triLeptonVeto)=trilep;
  pass.at(triLeptonVeto)=(value.at(triLeptonVeto)==cut.at(triLeptonVeto));

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
  // jet veto
  //
  if(verbose)std::cout << "jet veto" << std::endl;
  if(verbose)std::cout << "Cleaning jets" << std::endl;
  std::vector<int> jetidx;
  bool etrackjet;
  bool mutrackjet;
  jetidx.clear();
  
  // loop over all jets & only save the ones that do not overlap with the selected muon/electron
  if(eidx1!=999 && muidx1!=999){
	  for(unsigned i=0;i<Ntp->NPFJets();i++){
		  etrackjet=false;
		  mutrackjet=false;
		  for(unsigned j=0;j<Ntp->PFJet_nTrk(i);j++){
			  if(Ntp->PFJet_TracksP4(i,j).DeltaR(Ntp->Muons_p4(GoodMuons.at(0)))<0.001){
				  mutrackjet=true;
			  }
			  if(Ntp->PFJet_TracksP4(i,j).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0)))<0.1){
				  etrackjet=true;
			  }
		  }
		  if(!mutrackjet && !etrackjet){
			  jetidx.push_back(i);
		  }
	  }
  }
  
  if(verbose)std::cout << "Looking for jets from Vtx" << std::endl;
  int leadingtrack;
  bool jetfromvtx;
  std::vector<int> jetsfromvtx;
  std::vector<int> leadingtracks;
  jetsfromvtx.clear();
  leadingtracks.clear();
  
  // loop over selected jets, find each leading track & check if it's assigned to selected vertex. save only jets from vertex
  for(unsigned i=0;i<jetidx.size();i++){
	  leadingtrack = 0;
	  jetfromvtx = false;
	  int counter = 0;
	  for(unsigned j=0;j<Ntp->PFJet_nTrk(jetidx.at(i));j++){
		  if(Ntp->PFJet_TracksP4(jetidx.at(i),j).Pt()>Ntp->PFJet_TracksP4(jetidx.at(i),leadingtrack).Pt()){
			  leadingtrack=j;
		  }
	  }
	  if(Ntp->PFJet_nTrk(jetidx.at(i))==0)leadingtrack = -1;
	  for(unsigned j=0;j<Ntp->Vtx_nTrk(pos);j++){
		  if(leadingtrack>=0 && Ntp->PFJet_TracksP4(jetidx.at(i),leadingtrack).DeltaR(Ntp->Vtx_TracksP4(pos,j))<0.00001){
			  jetfromvtx=true;
			  counter++;
		  }
	  }
	  if(counter>1)std::cout << "More than one vtx track associated to leading track from jet" << std::endl;
	  if(jetfromvtx){
		  jetsfromvtx.push_back(jetidx.at(i));
		  leadingtracks.push_back(leadingtrack);
	  }
  }
  
  if(verbose)std::cout<< "Find two highest pt jets" << std::endl;
  int firstjet_idx=-1;
  int secondjet_idx=-1;
  double initialpt=0.;
  
  // loop over jets from selected vertex & find the two jets with the highest pt
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()>initialpt){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt();
		  firstjet_idx=jetsfromvtx.at(i);
	  }
  }
  initialpt=0.;
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1 && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()>initialpt && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()<Ntp->PFJet_p4(firstjet_idx).Pt()){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt();
		  secondjet_idx=jetsfromvtx.at(i);
	  }
  }
  
  if(verbose)std::cout << "applying veto" << std::endl;
  value.at(jetVeto)=0;
  if(jetsfromvtx.size()>1){
	  value.at(jetVeto)=Ntp->PFJet_p4(firstjet_idx).Pt()+Ntp->PFJet_p4(secondjet_idx).Pt();
  }else if(jetsfromvtx.size()==1){
	  value.at(jetVeto)=Ntp->PFJet_p4(firstjet_idx).Pt();
	  cut.at(jetVeto)=40;
  }
  pass.at(jetVeto)=(value.at(jetVeto)<cut.at(jetVeto));
  
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
  
  ///////////////////////////////////////////////
  //
  // Pt balance cut
  //
  if(verbose) std::cout << "pt balance cut" << std::endl;
  value.at(ptBalance)=999.;
  for(unsigned i=0;i<GoodMuons.size();i++){
	  for(unsigned j=0;j<GoodElectrons.size();j++){
		  if(twod){
			  value.at(ptBalance) = fabs(Ntp->Muons_p4(GoodMuons.at(i)).Pt()-Ntp->Electron_p4(GoodElectrons.at(j)).Pt());
		  }else{
			  value.at(ptBalance) = (Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).Pt();
		  }
	  }
  }
  pass.at(ptBalance)=(value.at(ptBalance)<cut.at(ptBalance));
  
  ///////////////////////////////////////////////
  //
  // Invariant mass cut
  //
  if(verbose) std::cout << "Invariant mass cut" << std::endl;
  value.at(ZMassmax)=zmax+1.;
  value.at(ZMassmin)=zmin-1.;
  if(eidx1!=999 && muidx1!=999){
	  value.at(ZMassmax)=(Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M();
	  value.at(ZMassmin)=(Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M();
  }
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));
  pass.at(ZMassmin)=(value.at(ZMassmin)>cut.at(ZMassmin));
  
  ///////////////////////////////////////////////
  //
  // Delta phi cut
  //
  if(verbose) std::cout << "Delta phi cut" << std::endl;
  value.at(Phimin)=0.;
  value.at(Phimax)=2.;
  if(eidx1!=999 && muidx1!=999){
	  value.at(Phimin)=Ntp->Muons_p4(GoodMuons.at(0)).DeltaPhi(Ntp->Electron_p4(GoodElectrons.at(0)))/TMath::Pi();
	  value.at(Phimax)=Ntp->Muons_p4(GoodMuons.at(0)).DeltaPhi(Ntp->Electron_p4(GoodElectrons.at(0)))/TMath::Pi();
	  if(value.at(Phimin)<0)value.at(Phimin)+=2;
	  if(value.at(Phimax)<0)value.at(Phimax)+=2;
  }
  if(twod){
	  pass.at(Phimin)=(value.at(Phimin)>cut.at(Phimin));
	  pass.at(Phimax)=(value.at(Phimax)<cut.at(Phimax));
  }else{
	  pass.at(Phimin)=true;
	  pass.at(Phimax)=true;
  }
  
  ////////////////////////////////////////////////
  //
  // QCD
  //
  if(verbose) std::cout << "QCD" << std::endl;
  if(MVA_ID){
	  fakeRate = 0.;
	  if(!pass.at(charge)
		  && Ntp->isData()
		  ){
		  if(fakeMuon && !fakeElectron){
			  fakeRate = fakeRateMu;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
			  blind = 1;
		  }else if(fakeElectron && !fakeMuon){
			  fakeRate = fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge) = true;
			  blind = 1;
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
    w*=Ntp->EvtWeight3D()*blind;
    if(verbose)std::cout << "void  ZtoEMu::doEvent() k" << w << " " << wobs << std::endl;
  }else if(Ntp->GetMCID()==34){
		w*=blind*embedd_weight;
		if(pass.at(NE))w*=ElectronEffRecHit(GoodElectrons.at(0));
		if(pass.at(NMu))w*=MuonDataSF(GoodMuons.at(0));
	}
  else{w=1*fakeRate*blind;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  blindingactive = false;
  if(!blindingactive){
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
  
  NJets.at(t).Fill(Ntp->NPFJets(),w);
  
  if(jetsfromvtx.size()>1){
	  jetsum.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Pt()+Ntp->PFJet_p4(secondjet_idx).Pt(),w);
  }
  
  if(status && jetsfromvtx.size()>0)FirstJetPt.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Pt(),w);
  if(jetsfromvtx.size()>1)SecondJetPt.at(t).Fill(Ntp->PFJet_p4(secondjet_idx).Pt(),w);
  if(jetsfromvtx.size()>2)ThirdJetPt.at(t).Fill(Ntp->PFJet_p4(jetsfromvtx.at(2)).Pt(),w);
  if(jetsfromvtx.size()>3)FourthJetPt.at(t).Fill(Ntp->PFJet_p4(jetsfromvtx.at(3)).Pt(),w);
  
  
  if(verbose)std::cout << "looping over good electrons" << std::endl;
  for(unsigned int i=0;i<GoodElectrons.size();i++){
		EPt.at(t).Fill(Ntp->Electron_p4(GoodElectrons.at(i)).Pt(),w);
		etaE.at(t).Fill(Ntp->Electron_p4(GoodElectrons.at(i)).Eta(),w);
		mtE.at(t).Fill(sqrt(2*Ntp->Electron_p4(GoodElectrons.at(i)).Pt()*Ntp->MET_Corr_et()*(1-cosphi2d(Ntp->Electron_p4(GoodElectrons.at(i)).Px(),Ntp->Electron_p4(GoodElectrons.at(i)).Py(),Ntp->MET_Corr_ex(),Ntp->MET_Corr_ey()))),w);
		RelIsoE.at(t).Fill(Ntp->Electron_RelIso(GoodElectrons.at(i)),w);
  }
  if(verbose)std::cout << "looping over good muons" << std::endl;	
  for(unsigned int i=0;i<GoodMuons.size();i++){
	  MuPt.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).Pt(),w);
	  etaMu.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).Eta(),w);
	  mtMu.at(t).Fill(sqrt(2*Ntp->Muons_p4(GoodMuons.at(i)).Pt()*Ntp->MET_Corr_et()*(1-cosphi2d(Ntp->Muons_p4(GoodMuons.at(i)).Px(),Ntp->Muons_p4(GoodMuons.at(i)).Py(),Ntp->MET_Corr_ex(),Ntp->MET_Corr_ey()))),w);
	  RelIsoMu.at(t).Fill(Ntp->Muon_RelIso(GoodMuons.at(i)),w);
  }
  
  if(verbose)std::cout << "Calculating pzeta & poca difference" << std::endl;
  if(verbose)std::cout << "Calculating charges and dr(e,mu)" << std::endl;
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  for(unsigned j=0;j<Ntp->NElectrons();j++){
		  double dp = Ntp->Muons_p4(i).DeltaPhi(Ntp->Electron_p4(j))/TMath::Pi();
		  if(dp<0)dp+=2;
		  deltaphi.at(t).Fill(dp,w);
		  drmue.at(t).Fill(Ntp->Muons_p4(i).DeltaR(Ntp->Electron_p4(j)),w);
		  chargesum.at(t).Fill(fabs(Ntp->Muon_Charge(i)+Ntp->Electron_Charge(j)),w);
		  chargesumsigned.at(t).Fill(Ntp->Muon_Charge(i)+Ntp->Electron_Charge(j),w);
	  }
  }

  if(pass.at(qualitycuts)){
	  // projections from 2d plot invmass vs deltaphi
	  double mass = (Ntp->Muons_p4(qualitymuons.at(0))+Ntp->Electron_p4(qualityelectrons.at(0))).M();
	  double qualityphi = Ntp->Muons_p4(qualitymuons.at(0)).DeltaPhi(Ntp->Electron_p4(qualityelectrons.at(0)))/TMath::Pi();
	  if(qualityphi<0)qualityphi+=2;
	  if(pass.at(TriggerOk)
			  && pass.at(PrimeVtx)
			  && pass.at(qualitycuts)
			  && pass.at(SameVtx)
			  && pass.at(NMuPt)
			  && pass.at(NMuEta)
			  && pass.at(NEPt)
			  && pass.at(NEEta)
			  && pass.at(NMu)
			  && pass.at(NE)
			  && pass.at(drMuE)
			  && pass.at(charge)
			  && pass.at(jetVeto)
			  && pass.at(MtMu)
			  && pass.at(ptBalance)
			  ){
		  if(mass<81 && qualityphi<2.5){
			  mproj1.at(t).Fill(mass,w);
			  phiproj1.at(t).Fill(qualityphi,w);
		  }else if(mass>=81 && mass<=102 && qualityphi<2.5){
			  mproj2.at(t).Fill(mass,w);
			  phiproj2.at(t).Fill(qualityphi,w);
		  }else if(mass>102 && qualityphi<2.5){
			  mproj3.at(t).Fill(mass,w);
			  phiproj3.at(t).Fill(qualityphi,w);
		  }else if(mass<81 && qualityphi>=2.5 && qualityphi<=3.7){
			  mproj4.at(t).Fill(mass,w);
			  phiproj4.at(t).Fill(qualityphi,w);
		  }else if(mass>=81 && mass<=102 && qualityphi>=2.5 && qualityphi<=3.7){
			  mproj5.at(t).Fill(mass,w);
			  phiproj5.at(t).Fill(qualityphi,w);
		  }else if(mass>102 && qualityphi>=2.5 && qualityphi<=3.7){
			  mproj6.at(t).Fill(mass,w);
			  phiproj6.at(t).Fill(qualityphi,w);
		  }else if(mass<81 && qualityphi>3.7){
			  mproj7.at(t).Fill(mass,w);
			  phiproj7.at(t).Fill(qualityphi,w);
		  }else if(mass>=81 && mass<=102 && qualityphi>3.7){
			  mproj8.at(t).Fill(mass,w);
			  phiproj8.at(t).Fill(qualityphi,w);
		  }else if(mass>102 && qualityphi>3.7){
			  mproj9.at(t).Fill(mass,w);
			  phiproj9.at(t).Fill(qualityphi,w);
		  }
		  if(mass<81){
			  phi1.at(t).Fill(qualityphi,w);
		  }else if(mass>=81 && mass<=102){
			  phi2.at(t).Fill(qualityphi,w);
		  }else if(mass>102){
			  phi3.at(t).Fill(qualityphi,w);
		  }
		  if(qualityphi<(TMath::Pi()-0.6)){
			  m1.at(t).Fill(mass,w);
		  }else if(qualityphi>=(TMath::Pi()-0.6) && qualityphi<=(TMath::Pi()+0.6)){
			  m2.at(t).Fill(mass,w);
		  }else if(qualityphi>(TMath::Pi()+0.6)){
			  m3.at(t).Fill(mass,w);
		  }
	  }
	  if(pass.at(TriggerOk)
			  && pass.at(PrimeVtx)
			  && pass.at(qualitycuts)
			  && pass.at(SameVtx)
			  && pass.at(NMuPt)
			  && pass.at(NMuEta)
			  && pass.at(NEPt)
			  && pass.at(NEEta)
			  && pass.at(NMu)
			  && pass.at(NE)
			  && pass.at(drMuE)
			  && pass.at(charge)
			  && pass.at(jetVeto)
			  && pass.at(MtMu)
			  ){
		  if(mass<81){
			  if(fabs(Ntp->Muons_p4(qualitymuons.at(0)).Pt()-Ntp->Electron_p4(qualityelectrons.at(0)).Pt())<20.){
				  phi1_ptdiff.at(t).Fill(qualityphi,w);
			  }
			  phi1_nopt.at(t).Fill(qualityphi,w);
		  }else if(mass>=81 && mass<=102){
			  if(fabs(Ntp->Muons_p4(qualitymuons.at(0)).Pt()-Ntp->Electron_p4(qualityelectrons.at(0)).Pt())<20.){
				  phi2_ptdiff.at(t).Fill(qualityphi,w);
			  }
			  phi2_nopt.at(t).Fill(qualityphi,w);
		  }else if(mass>102){
			  if(fabs(Ntp->Muons_p4(qualitymuons.at(0)).Pt()-Ntp->Electron_p4(qualityelectrons.at(0)).Pt())<20.){
				  phi3_ptdiff.at(t).Fill(qualityphi,w);
			  }
			  phi3_nopt.at(t).Fill(qualityphi,w);
		  }
		  pzetaCut.at(t).Fill(calculatePzetaDQM(qualitymuons.at(0),qualityelectrons.at(0)),w);
		  metCut.at(t).Fill(Ntp->MET_Corr_et(),w);
		  if(pass.at(ptBalance)
				  && pass.at(ZMassmax)
				  && pass.at(ZMassmin)
				  ){
			  pzetaStatus.at(t).Fill(calculatePzetaDQM(qualitymuons.at(0),qualityelectrons.at(0)),w);
			  metStatus.at(t).Fill(Ntp->MET_Corr_et(),w);
		  }
	  }
  }

  for(unsigned i=0;i<GoodMuons.size();i++){
	  for(unsigned j=0;j<GoodElectrons.size();j++){
		  pzeta.at(t).Fill(calculatePzeta(i,j,GoodElectrons,GoodMuons),w);
		  pzetaDQM.at(t).Fill(calculatePzetaDQM(i,j),w);
		  ptbal.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).Pt(),w);
		  if(pass.at(qualitycuts))ptbal2.at(t).Fill(fabs(Ntp->Muons_p4(GoodMuons.at(i)).Pt()-Ntp->Electron_p4(GoodElectrons.at(j)).Pt()),w);
		  double dpID = Ntp->Muons_p4(GoodMuons.at(i)).DeltaPhi(Ntp->Electron_p4(GoodElectrons.at(j)))/TMath::Pi();
		  if(dpID<0)dpID+=2;
		  deltaphiID.at(t).Fill(dpID,w);
		  if(Ntp->Muons_p4(GoodMuons.at(i)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(j)))>0.2){
			if(pass.at(qualitycuts))InvmassVsDeltaPhi.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).M(),dpID,w);
			if(pass.at(qualitycuts))PtDiffVsDeltaPhi.at(t).Fill(fabs(Ntp->Muons_p4(GoodMuons.at(i)).Pt()-Ntp->Electron_p4(GoodElectrons.at(j)).Pt()),dpID,w);
			if(pass.at(qualitycuts))invmass.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(i))+Ntp->Electron_p4(GoodElectrons.at(j))).M(),w);
			chargesumID.at(t).Fill(fabs(Ntp->Muon_Charge(GoodMuons.at(i))+Ntp->Electron_Charge(GoodElectrons.at(j))),w);
			chargesumIDsigned.at(t).Fill(Ntp->Muon_Charge(GoodMuons.at(i))+Ntp->Electron_Charge(GoodElectrons.at(j)),w);
		  }
		  drmueID.at(t).Fill(Ntp->Muons_p4(GoodMuons.at(i)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(j))),w);
	  }
  }
  
  if(verbose)std::cout << "invariant mass with different cuts applied" << std::endl;
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(qualitycuts)
		  && pass.at(SameVtx)
		  && pass.at(NMuPt)
		  && pass.at(NMuEta)
		  && pass.at(NEPt)
		  && pass.at(NEEta)
		  && pass.at(NMu)
		  && pass.at(NE)
		  ){
	  invmass_only_object_id.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
	  if(pass.at(drMuE)){
		  invmass_dremu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  if(pass.at(looseMuonVeto)){
			  invmass_loosemuonveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
			  if(pass.at(charge)){
				  invmass_charge.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
				  if(pass.at(jetVeto)){
					  invmass_jetveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
					  if(pass.at(MtMu)){
						  invmass_mtmu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
						  if(pass.at(ptBalance)){
							  invmass_ptbalance.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
							  if(pass.at(ZMassmax)
									  && pass.at(ZMassmin)){
								  invmass_zmass.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
							  }
						  }
					  }
				  }
			  }
		  if(pass.at(charge)
				  && pass.at(jetVeto)
				  && pass.at(MtMu)
				  ){
			  nm2_ptbalance.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge)
				  && pass.at(jetVeto)
				  && pass.at(ptBalance)
				  ){
			  nm2_mtmu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge)
				  && pass.at(MtMu)
				  && pass.at(ptBalance)
				  ){
			  nm2_jetveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(jetVeto)
				  && pass.at(MtMu)
				  && pass.at(ptBalance)
				  ){
			  nm2_charge.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge)
				  && pass.at(jetVeto)
				  && pass.at(MtMu)
				  && pass.at(ptBalance)
				  ){
			  nm2_drmue.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge)
				  && pass.at(ptBalance)
				  ){
				  no_jetveto_mtmu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(looseMuonVeto)
				  && pass.at(jetVeto)
				  && pass.at(MtMu)
				  ){
				  no_charge_ptbalance.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(looseMuonVeto)
				  && pass.at(jetVeto)
				  && pass.at(MtMu)
				  && pass.at(ptBalance)
				  ){
			  no_charge.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  }
		  if(pass.at(charge))drmue_plus_charge.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  if(pass.at(jetVeto))drmue_plus_jetveto.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  if(pass.at(MtMu))drmue_plus_mtmu.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
		  if(pass.at(ptBalance))drmue_plus_ptbalance.at(t).Fill((Ntp->Muons_p4(GoodMuons.at(0))+Ntp->Electron_p4(GoodElectrons.at(0))).M(),w);
	  }
  }
  }
  
  if(fakeMuon)frMu.at(t).Fill(fakeRateMu);
  if(fakeElectron)frE.at(t).Fill(fakeRateE);

  }
  
  double mudr(999),edr(999);
  int muid(0),eid(0);
  
  if(status){
	  if(id==30){
		  std::cout << "##################" << std::endl;
		  for(unsigned i=0;i<Ntp->NMCParticles();i++){
			  if(fabs(Ntp->MC_pdgid(i))==11
					  && (Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->MC_p4(i))<=0.3 || Ntp->Electron_p4(GoodElectrons.at(0)).DeltaR(Ntp->MC_p4(i))<=0.3)
					  ){
				  std::cout << "------------------------------" << std::endl;
				  std::cout << "electron mother is " << Ntp->MC_pdgid(Ntp->MC_midx(i)) << std::endl;
				  std::cout << "electron has charge " << Ntp->MC_charge(i) << std::endl;
				  std::cout << "selected mu charge is " << Ntp->Muon_Charge(GoodMuons.at(0)) << std::endl;
				  std::cout << "selected e charge is " << Ntp->Electron_Charge(GoodElectrons.at(0)) << std::endl;
				  std::cout << "dR of this e & selected mu is " << Ntp->MC_p4(i).DeltaR(Ntp->Muons_p4(GoodMuons.at(0))) << std::endl;
				  std::cout << "dR of this e & selected e is " << Ntp->MC_p4(i).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0))) << std::endl;
				  std::cout << "dR of selected e & mu is " << Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0))) << std::endl;
			  }
			  if(fabs(Ntp->MC_pdgid(i))==13
					  && (Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->MC_p4(i))<=0.3 || Ntp->Electron_p4(GoodElectrons.at(0)).DeltaR(Ntp->MC_p4(i))<=0.3)
					  ){
				  std::cout << "------------------------------" << std::endl;
				  std::cout << "muon mother is " << Ntp->MC_pdgid(Ntp->MC_midx(i)) << std::endl;
				  std::cout << "muon has charge " << Ntp->MC_charge(i) << std::endl;
				  std::cout << "selected mu charge is " << Ntp->Muon_Charge(GoodMuons.at(0)) << std::endl;
				  std::cout << "selected e charge is " << Ntp->Electron_Charge(GoodElectrons.at(0)) << std::endl;
				  std::cout << "dR of this mu & selected mu is " << Ntp->MC_p4(i).DeltaR(Ntp->Muons_p4(GoodMuons.at(0))) << std::endl;
				  std::cout << "dR of this mu & selected e is " << Ntp->MC_p4(i).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0))) << std::endl;
				  std::cout << "dR of selected e & mu is " << Ntp->Muons_p4(GoodMuons.at(0)).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0))) << std::endl;
				  for(unsigned i=0;i<Ntp->NMuons();i++){
					  if(i!=GoodMuons.at(0) && Ntp->Muons_p4(i).DeltaR(Ntp->Electron_p4(GoodElectrons.at(0)))<=0.3){
						  std::cout << "!!!!!!!!!!!!!!!!!!" << std::endl;
						  std::cout << "muon in dR-cone <=0.3 of selected electron" << std::endl;
						  std::cout << "muon has pt of " << Ntp->Muons_p4(i).Pt() << std::endl;
						  
					  }
				  }
			  }
		  }
	  }
  }

  if(verbose)std::cout << "ZtoEMu::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

double ZtoEMu::calculatePzeta(int muiterator, int eiterator,std::vector<unsigned int> vec1,std::vector<unsigned int> vec2){
  pex=Ntp->Electron_p4(vec1.at(eiterator)).Px();
  pey=Ntp->Electron_p4(vec1.at(eiterator)).Py();
  pmux=Ntp->Muons_p4(vec2.at(muiterator)).Px();
  pmuy=Ntp->Muons_p4(vec2.at(muiterator)).Py();
  phie=Ntp->Electron_p4(vec1.at(eiterator)).Phi();
  phimu=Ntp->Muons_p4(vec2.at(muiterator)).Phi();
  combpt=TMath::Sqrt(pow(pex+pmux,2)+pow(pey+pmuy,2));
  aemu=TMath::ACos(pmux*pex+pmuy*pey/(Ntp->Muons_p4(vec2.at(muiterator)).Pt()*Ntp->Electron_p4(vec1.at(eiterator)).Pt()));
  if(phie<phimu && fabs(phie-phimu)<TMath::Pi())phismall=phie;
  else if(phimu<phie && fabs(phie-phimu)>TMath::Pi())phismall=phie;
  else if(phie<phimu && fabs(phie-phimu)>TMath::Pi())phismall=phimu;
  else if(phimu<phie && fabs(phie-phimu)<TMath::Pi())phismall=phimu;
  beta=TMath::ACos(((pex+pmux)*TMath::Cos(phismall+0.5*aemu)+(pey+pmuy)*TMath::Sin(phismall+0.5*aemu))/combpt);
  gamma=TMath::ACos((Ntp->MET_Corr_ex()*TMath::Cos(phismall+0.5*aemu)+Ntp->MET_Corr_ey()*TMath::Sin(phismall+0.5*aemu))/Ntp->MET_Corr_et());
  if(Ntp->MET_Corr_phi()>(phismall+0.5*aemu+0.5*TMath::Pi()) && Ntp->MET_Corr_phi()<(phismall+0.5*aemu+1.5*TMath::Pi()))gamma*=-1;
  pvis=TMath::Sin(beta)*combpt;
  pmiss=TMath::Sin(gamma)*Ntp->MET_Corr_et();
  return pmiss-pvis;
}

double ZtoEMu::calculatePzetaDQM(int muiterator, int eiterator){
	double cosPhi1 = TMath::Cos(Ntp->Electron_p4(eiterator).Phi());
	double sinPhi1 = TMath::Sin(Ntp->Electron_p4(eiterator).Phi());
	double cosPhi2 = TMath::Cos(Ntp->Muons_p4(muiterator).Phi());
	double sinPhi2 = TMath::Sin(Ntp->Muons_p4(muiterator).Phi());
	double zetaX = cosPhi1 + cosPhi2;
	double zetaY = sinPhi1 + sinPhi2;
	double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
	if(zetaR>0.){
		zetaX/=zetaR;
		zetaY/=zetaR;
	}
	double pxVis=Ntp->Electron_p4(eiterator).Px()+Ntp->Muons_p4(muiterator).Px();
	double pyVis=Ntp->Electron_p4(eiterator).Py()+Ntp->Muons_p4(muiterator).Py();
	double pZetaVis=pxVis*zetaX+pyVis*zetaY;
	double px=pxVis+Ntp->MET_Corr_ex();
	double py=pyVis+Ntp->MET_Corr_ey();
	double pZeta=px*zetaX+py*zetaY;
	return pZeta-1.5*pZetaVis;
}

double ZtoEMu::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu::cosphi3d(TVector3 vec1, TVector3 vec2){
	return (vec1.Dot(vec2))/vec1.Mag()/vec2.Mag();
}

bool ZtoEMu::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu::isTightMuon(unsigned int i){
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

bool ZtoEMu::isTightMuon(unsigned int i, unsigned int j){
	if(isTightMuon(i) &&
			dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))<0.2 &&
			dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))<0.5
			  ){
		return true;
	}
	return false;
}

bool ZtoEMu::isLooseMuon(unsigned int i){
	if(Ntp->Muon_isPFMuon(i) &&
			(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))
			){
		return true;
	}
	return false;
}

bool ZtoEMu::isFakeMuon(unsigned int i){
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

bool ZtoEMu::isFakeMuon(unsigned int i, unsigned int j){
	if(isFakeMuon(i)
			&& dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))<0.2
			){
		return true;
	}
	return false;
}

double ZtoEMu::Muon_RelIso(unsigned int i){
	return (Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i)))/Ntp->Muons_p4(i).Pt();
}

double ZtoEMu::Muon_AbsIso(unsigned int i){
	return Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i));
}

//////////////////////////////
//
// Electron related functions
//

bool ZtoEMu::isMVAElectron(unsigned int i){
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

bool ZtoEMu::isTightElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){ //barrel
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

bool ZtoEMu::isTightElectron(unsigned int i, unsigned int j){
	if(isTightElectron(i)
			&& dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.02
			&& dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.1
			){
		return true;
	}
	return false;
}

bool ZtoEMu::isLooseElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){ //barrel
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.007 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.15 &&
				Ntp->Electron_sigmaIetaIeta(i)<0.01 &&
				Ntp->Electron_hadronicOverEm(i)<0.12 &&
				fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))<0.05 &&
				Ntp->Electron_RelIso(i)<0.15 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_numberOfMissedHits(i)<=1
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(i))>1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){ //endcaps
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.009 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.10 &&
				Ntp->Electron_sigmaIetaIeta(i)<0.03 &&
				Ntp->Electron_hadronicOverEm(i)<0.10 &&
				fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))<0.05 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_numberOfMissedHits(i)<=1
				){
			if(Ntp->Electron_p4(i).Pt()>=20.0 && Electron_RelIso(i)<0.15){
				return true;
			}else if(Ntp->Electron_p4(i).Pt()<20.0 && Electron_RelIso(i)<0.10){
				return true;
			}
		}
	}
	return false;
}

bool ZtoEMu::isFakeElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){
		if(Ntp->Electron_p4(i).Pt()>20 &&
				!Ntp->Electron_HasMatchedConversions(i) &&
				Ntp->Electron_sigmaIetaIeta(i)<0.01 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)<0.15 &&
				Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)<0.007 &&
				Electron_RelIso(i)<0.2
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(i))>1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
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

bool ZtoEMu::isFakeElectron(unsigned int i, unsigned int j){
	if(isFakeElectron(i)
			&& dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.1
			&& dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))<0.03
			){
		return true;
	}
	return false;
}

double ZtoEMu::Electron_RelIso(unsigned int i){
	return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();
}

double ZtoEMu::Electron_Aeff(double Eta){
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

double ZtoEMu::MuonSF(unsigned int i){
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

double ZtoEMu::MuonDataSF(unsigned int i){
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

double ZtoEMu::ElectronSF(unsigned int i){
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

double ZtoEMu::ElectronDataSF(unsigned int i){
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

double ZtoEMu::ElectronEffRecHit(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = Ntp->Electron_supercluster_eta(i);
	
	Double_t xAxis1[10] = {10, 15, 20, 25, 30, 40, 55, 70, 100, 200}; 
	Double_t yAxis1[4] = {0, 0.8, 1.479, 2.5};
	   
	TH2D* hPtEtaSFL = new TH2D("hPtEtaSFL","",9,xAxis1,3,yAxis1);
	hPtEtaSFL->SetBinContent(12,0.81);
	hPtEtaSFL->SetBinContent(13,0.91);
	hPtEtaSFL->SetBinContent(14,0.95);
	hPtEtaSFL->SetBinContent(15,0.96);
	hPtEtaSFL->SetBinContent(16,0.97);
	hPtEtaSFL->SetBinContent(17,0.98);
	hPtEtaSFL->SetBinContent(18,0.99);
	hPtEtaSFL->SetBinContent(19,0.98);
	hPtEtaSFL->SetBinContent(20,0.99);
	hPtEtaSFL->SetBinContent(21,0.98);
	hPtEtaSFL->SetBinContent(23,0.78);
	hPtEtaSFL->SetBinContent(24,0.89);
	hPtEtaSFL->SetBinContent(25,0.92);
	hPtEtaSFL->SetBinContent(26,0.94);
	hPtEtaSFL->SetBinContent(27,0.94);
	hPtEtaSFL->SetBinContent(28,0.97);
	hPtEtaSFL->SetBinContent(29,0.97);
	hPtEtaSFL->SetBinContent(30,0.99);
	hPtEtaSFL->SetBinContent(31,1.00);
	hPtEtaSFL->SetBinContent(32,1.00);
	hPtEtaSFL->SetBinContent(34,0.46);
	hPtEtaSFL->SetBinContent(35,0.66);
	hPtEtaSFL->SetBinContent(36,0.73);
	hPtEtaSFL->SetBinContent(37,0.80);
	hPtEtaSFL->SetBinContent(38,0.83);
	hPtEtaSFL->SetBinContent(39,0.86);
	hPtEtaSFL->SetBinContent(40,0.88);
	hPtEtaSFL->SetBinContent(41,0.91);
	hPtEtaSFL->SetBinContent(42,0.93);
	hPtEtaSFL->SetBinContent(43,1.00);
	
	if(pt>199.99)pt=199.9;
	eta=fabs(eta);
	if(eta>2.49)eta=2.49;
	if(pt<10)return 0;
	
	Float_t eff=0;
	Int_t bin = hPtEtaSFL->FindFixBin(pt,eta);
	eff = hPtEtaSFL->GetBinContent(bin);
	//std::cout<<" pt="<<pt<<", eta="<<eta<<", eff="<<eff<<std::endl;
	
	return eff;
}

//////////////////////////////
//
// Calculate fakerate
//

double ZtoEMu::Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type){
	
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

void ZtoEMu::Finish(){
	unsigned int data=0;
	unsigned int qcd=1;
	unsigned int ww=2;
	unsigned int wz=3;
	unsigned int zz=4;
	unsigned int tt=5;
	unsigned int dytautau=6;
	unsigned int dyll=7;
	unsigned int wjets=8;
	unsigned int sig=7;//9;
	
	//if(Nminus0.at(0).at(data).Integral()!=0)if(HConfig.GetHisto(false,DataMCType::QCD,data))ScaleAllHistOfType(data,1./Nminus0.at(0).at(data).Integral());
	//if(Nminus0.at(0).at(dytautau).Integral()!=0)if(HConfig.GetHisto(false,DataMCType::DY_embedded,dytautau))ScaleAllHistOfType(dytautau,19789.302*1121.0/9770590.);
	//if(Nminus1.at(0).at(dytautau).Integral()!=0)if(HConfig.GetHisto(false,DataMCType::DY_tautau,dytautau))ScaleAllHistOfType(dytautau,19789.302*1168.0/869489./Nminus1.at(0).at(dytautau).Integral());

	std::cout << "#############################" << std::endl;
	std::cout << "Number of data events in sideband phi = " << phi1.at(data).Integral() << std::endl;
	std::cout << "Number of QCD events in sideband phi = " << phi1.at(qcd).Integral() << std::endl;
	std::cout << "Number of WW events in sideband phi = " << phi1.at(ww).Integral()*19789.302*5.824/1933235 << std::endl;
	std::cout << "Number of WZ events in sideband phi = " << phi1.at(wz).Integral()*19789.302*1.058/2979624 << std::endl;
	std::cout << "Number of ZZ events in sideband phi = " << phi1.at(zz).Integral()*19789.302*8.4/9799908 << std::endl;
	std::cout << "Number of ttbar events in sideband phi = " << phi1.at(tt).Integral()*19789.302*225.2/6923750 << std::endl;
	std::cout << "Number of dytautau events in sideband phi = " << phi1.at(dytautau).Integral()*19789.302*1168.0/10152445 << std::endl;
	std::cout << "Number of signal events in sideband phi = " << phi1.at(sig).Integral()*19789.302*0.057/100000 << std::endl;
	std::cout << "#############################" << std::endl;
	std::cout << "Number of data events in invmass = " << invmass_ptbalance.at(data).Integral() << std::endl;
	std::cout << "Number of QCD events in invmass = " << invmass_ptbalance.at(qcd).Integral() << std::endl;
	std::cout << "Number of WW events in invmass = " << invmass_ptbalance.at(ww).Integral()*19789.302*5.824/1933235 << std::endl;
	std::cout << "Number of WZ events in invmass = " << invmass_ptbalance.at(wz).Integral()*19789.302*1.058/2979624 << std::endl;
	std::cout << "Number of ZZ events in invmass = " << invmass_ptbalance.at(zz).Integral()*19789.302*8.4/9799908 << std::endl;
	std::cout << "Number of ttbar events in invmass = " << invmass_ptbalance.at(tt).Integral()*19789.302*225.2/6923750 << std::endl;
	std::cout << "Number of dytautau events in invmass = " << invmass_ptbalance.at(dytautau).Integral()*19789.302*1168.0/10152445 << std::endl;
	std::cout << "Number of signal events in invmass = " << invmass_ptbalance.at(sig).Integral()*19789.302*0.057/100000 << std::endl;
	std::cout << "#############################" << std::endl;
	std::cout << "Number of data events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(data).Integral() << std::endl;
	std::cout << "Number of QCD events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(qcd).Integral() << std::endl;
	std::cout << "Number of WW events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(ww).Integral()*19789.302*5.824/1933235 << std::endl;
	std::cout << "Number of WZ events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(wz).Integral()*19789.302*1.058/2979624 << std::endl;
	std::cout << "Number of ZZ events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(zz).Integral()*19789.302*8.4/9799908 << std::endl;
	std::cout << "Number of ttbar events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(tt).Integral()*19789.302*225.2/6923750 << std::endl;
	std::cout << "Number of dytautau events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(dytautau).Integral()*19789.302*1168.0/10152445 << std::endl;
	std::cout << "Number of signal events in Nminus1 = " << Nminus1.at(Nminus1.size()-1).at(sig).Integral()*19789.302*0.057/100000 << std::endl;
	std::cout << "#############################" << std::endl;
	
	TH1D* qcdhist = new TH1D("qcdhist","QCD",20,60.,120.);
	TH1D* wwhist = new TH1D("wwhist","WW",20,60.,120.);
	TH1D* wzhist = new TH1D("wzhist","WZ",20,60.,120.);
	TH1D* zzhist = new TH1D("zzhist","ZZ",20,60.,120.);
	TH1D* tthist = new TH1D("tthist","tt",20,60.,120.);
	TH1D* dyhist = new TH1D("dyhist","dy",20,60.,120.);
	
	qcdhist = (TH1D*)Nminus1.at(Nminus1.size()-1).at(1).Clone("qcdhist");
	wwhist = (TH1D*)Nminus1.at(Nminus1.size()-1).at(2).Clone("wwhist");
	wzhist = (TH1D*)Nminus1.at(Nminus1.size()-1).at(3).Clone("wzhist");
	zzhist = (TH1D*)Nminus1.at(Nminus1.size()-1).at(4).Clone("zzhist");
	tthist = (TH1D*)Nminus1.at(Nminus1.size()-1).at(5).Clone("tthist");
	dyhist = (TH1D*)Nminus1.at(Nminus1.size()-1).at(6).Clone("dyhist");
	
	if(wwhist->Integral()!=0)wwhist->Scale(93.7257/wwhist->Integral());
	if(wzhist->Integral()!=0)wzhist->Scale(2.44924/wzhist->Integral());
	if(zzhist->Integral()!=0)zzhist->Scale(0.383402/zzhist->Integral());
	if(tthist->Integral()!=0)tthist->Scale(35.1889/tthist->Integral());
	if(dyhist->Integral()!=0)dyhist->Scale(1119.08/dyhist->Integral());
	
	for(unsigned i=1;i<7;i++){
		std::cout << "integral of " << i << " = " << Nminus1.at(Nminus1.size()-1).at(i).Integral() << std::endl;
	}
	
	std::cout << "qcd integral = " << qcdhist->Integral() << std::endl;
	std::cout << "ww integral = " << wwhist->Integral() << std::endl;
	std::cout << "wz integral = " << wzhist->Integral() << std::endl;
	std::cout << "zz integral = " << zzhist->Integral() << std::endl;
	std::cout << "tt integral = " << tthist->Integral() << std::endl;
	std::cout << "dy integral = " << dyhist->Integral() << std::endl;
	
	TH1D* mchist = new TH1D("mchist","Data-MC comparison",7,60.,81.);
	TH1D* datahist = new TH1D("datahist","Data-MC comparison",7,60.,81.);
	TH1D* total = new TH1D("total","Data-MC comparison",7,60.,81.);
	//mchist->Sumw2();
	double mc_bincontent = 0.;
	double data_bincontent = 0.;
	for(unsigned i=1;i<8;i++){
		mc_bincontent = 0.;
		data_bincontent = Nminus1.at(Nminus1.size()-1).at(data).GetBinContent(i);
		mc_bincontent = qcdhist->GetBinContent(i) + wwhist->GetBinContent(i) + wzhist->GetBinContent(i) + zzhist->GetBinContent(i) + tthist->GetBinContent(i) + dyhist->GetBinContent(i);
		mchist->SetBinContent(i,mc_bincontent);
		datahist->SetBinContent(i,data_bincontent);
	}
	std::cout << "total mc = " << mchist->Integral() << std::endl;
	std::cout << "total data = " << datahist->Integral() << std::endl;
	std::cout << "#############################" << std::endl;
	std::cout << "ks result = " << datahist->KolmogorovTest(mchist,"N") << std::endl;
	std::cout << "chi2 test = " << datahist->Chi2Test(mchist,"UW") << std::endl;
	std::cout << "#############################" << std::endl;
	
	Selection::Finish();
}
