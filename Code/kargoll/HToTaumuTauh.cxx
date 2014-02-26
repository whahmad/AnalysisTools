#include "HToTaumuTauh.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

HToTaumuTauh::HToTaumuTauh(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_dxy(0.045),
  cMu_dz(0.2),
  cMu_relIso(0.1),
  cMu_pt(20.0),
  cMu_eta(2.1),
  cTau_pt(20.0),
  cTau_eta(2.3),
  cTau_rawIso(1.5),
  cMuTau_dR(0.3),
  cMuTriLep_pt(10.0),
  cMuTriLep_eta(2.4),
  cEleTriLep_pt(10.0),
  cEleTriLep_eta(2.5),
  cCat_jetPt(30.0),
  cCat_jetEta(4.7),
  cCat_bjetPt(20.0),
  cCat_bjetEta(2.4),
  cCat_btagDisc(0.679), // medium WP, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP#B_tagging_Operating_Points_for_5
  cCat_splitTauPt(40.0),
  cJetClean_dR(0.5)
{
	TString trigNames[] = {"HLT_IsoMu18_eta2p1_LooseIsoPFTau20","HLT_IsoMu17_eta2p1_LooseIsoPFTau20"};
	std::vector<TString> temp (trigNames, trigNames + sizeof(trigNames) / sizeof(TString) );
	cTriggerNames = temp;

	// implemented categories:
	// VBF, OneJetHigh, OneJetLow, ZeroJetHigh, ZeroJetLow, NoCategory //TODO: implement b-tagging categories
	// TODO: implement b-tag veto in all categories
	categoryFlag = "NoCategory";
}

HToTaumuTauh::~HToTaumuTauh(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "HToTaumuTauh::~HToTaumuTauh Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "HToTaumuTauh::~HToTaumuTauh()" << std::endl;
}

void  HToTaumuTauh::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    	cut.at(TriggerOk)=0;
    if(i==PrimeVtx)     	cut.at(PrimeVtx)=1;
    if(i==NMuId)			cut.at(NMuId)=1;
    if(i==NMuKin)			cut.at(NMuKin)=1;
    if(i==DiMuonVeto)		cut.at(DiMuonVeto)=false;
    if(i==NTauId)			cut.at(NTauId)=1;
    if(i==NTauIso)			cut.at(NTauIso)=1;
    if(i==NTauKin)			cut.at(NTauKin)=1;
    if(i==TriLeptonVeto)	cut.at(TriLeptonVeto)=0;
    if(i==OppCharge)		cut.at(OppCharge)=0;
    if(i==MT)				cut.at(MT)=20.0;
    //category-specific values are set in the corresponding configure function
    // set them to dummy value -10.0 here
    cut_VBF.push_back(-10.);
    cut_OneJet.push_back(-10.);
    cut_ZeroJet.push_back(-10.);
    if(i>=CatCut1){
    	cut.at(i)=-10.0;
    }
  }

  // Setup Category Cut Values
  for(unsigned i = CatCut1; i< NCuts; i++){
	  if(i==VbfNJet)		cut_VBF.at(VbfNJet)		 =2;
	  if(i==VbfDeltaEta)	cut_VBF.at(VbfDeltaEta)	 =3.5;
	  if(i==VbfNJetRapGap)	cut_VBF.at(VbfNJetRapGap)=0;
	  if(i==VbfJetInvM)		cut_VBF.at(VbfJetInvM)	 =500.0;

	  if(i==OneJetNJet) 	cut_OneJet.at(OneJetNJet) = 1;
	  if(i==OneJetNoVBF)	cut_OneJet.at(OneJetNoVBF)= true;
	  if(i==OneJetNBtagJets)cut_OneJet.at(OneJetNBtagJets) = 0;
	  if(i==OneJetTauPt)	cut_OneJet.at(OneJetTauPt) = cCat_splitTauPt;

	  if(i==ZeroJetNJet) 		cut_ZeroJet.at(ZeroJetNJet) = 5;// 0;
	  if(i==ZeroJetNBtagJets)	cut_ZeroJet.at(ZeroJetNBtagJets)= 0;
	  if(i==ZeroJetTauPt)		cut_ZeroJet.at(ZeroJetTauPt) = cCat_splitTauPt;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i_cut=0; i_cut<NCuts; i_cut++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i_cut;
  
    if(i_cut==PrimeVtx){
      title.at(i_cut)="Number of Prime Vertices $(N>$";
      title.at(i_cut)+=cut.at(PrimeVtx);
      title.at(i_cut)+=")";
      htitle=title.at(i_cut);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i_cut==TriggerOk){
      title.at(i_cut)="Trigger ";
      hlabel="Trigger ";

      std::vector<TH1D> Nm1Temp = HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,cTriggerNames.size()+2,-1.5,cTriggerNames.size()+0.5,hlabel,"Events");
      std::vector<TH1D> Nm0Temp = HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,cTriggerNames.size()+2,-1.5,cTriggerNames.size()+0.5,hlabel,"Events");
      for (unsigned i_hist = 0; i_hist < Nm1Temp.size(); i_hist++){
    	  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(1,"not fired");
    	  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(1,"not fired");
    	  for (unsigned i_bin = 2; i_bin < cTriggerNames.size()+2; i_bin++){
    		  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(i_bin,cTriggerNames.at(i_bin-2));
    		  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(i_bin,cTriggerNames.at(i_bin-2));
    	  }
    	  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(cTriggerNames.size()+2,"multiple fired");
    	  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(cTriggerNames.size()+2,"multiple fired");
      }
      Nminus1.push_back(Nm1Temp);
      Nminus0.push_back(Nm0Temp);
    }
    else if(i_cut==NMuId){
    	title.at(i_cut)="Number $\\mu_{ID} >=$";
    	title.at(i_cut)+=cut.at(NMuId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NMuKin){
    	title.at(i_cut)="Number $\\mu_{sel} >=$";
    	title.at(i_cut)+=cut.at(NMuKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{sel}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i_cut==DiMuonVeto){
    	title.at(i_cut)="Veto on $\\mu_{veto}$ pair";
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Veto on #mu_{veto} pair";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i_cut==NTauId){
    	title.at(i_cut)="Number $\\tau_{ID} >=$";
    	title.at(i_cut)+=cut.at(NTauId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauId_",htitle,26,-0.5,25.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauId_",htitle,26,-0.5,25.5,hlabel,"Events"));
    }
    else if(i_cut==NTauIso){
    	title.at(i_cut)="Number $\\tau_{Iso} >=$";
    	title.at(i_cut)+=cut.at(NTauIso);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauIso_",htitle,16,-0.5,15.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauIso_",htitle,16,-0.5,15.5,hlabel,"Events"));
    }
    else if(i_cut==NTauKin){
    	title.at(i_cut)="Number $\\tau_{sel} >=$";
    	title.at(i_cut)+=cut.at(NTauKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{sel}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==TriLeptonVeto){
    	title.at(i_cut)="3 lepton veto: $N(\\mu)+N(e) =$";
    	title.at(i_cut)+=cut.at(TriLeptonVeto);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of tri-lepton veto leptons";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriLeptonVeto_",htitle,5,-0.5,4.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriLeptonVeto_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i_cut==OppCharge){
    	title.at(i_cut)="$q(\\mu)+q(\\tau) =$";
    	title.at(i_cut)+=cut.at(OppCharge);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="q(#mu)+q(#tau)";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OppCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OppCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }
    else if(i_cut==MT){
    	title.at(i_cut)="$m_{T}(\\mu,E_{T}^{miss}) <$";
    	title.at(i_cut)+=cut.at(MT);
    	title.at(i_cut)+=" GeV";
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="m_{T}(#mu,E_{T}^{miss})/GeV";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_",htitle,50,0.,100.,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_",htitle,50,0.,100.,hlabel,"Events"));
    }
    else if(i_cut>=CatCut1){
    	// set histograms to dummy values
    	// will be overwritten in configure_{Category} method
    	title.at(i_cut) = "Category Dummy ";
    	title.at(i_cut)	+=(i_cut-CatCut1);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel = title.at(i_cut);
    	TString n = Name+c+"_Nminus1_CatDummy";
    	n += (i_cut-CatCut1);
    	n += "_";
    	Nminus1.push_back(HConfig.GetTH1D(n,htitle,50,-50.,50.,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(n,htitle,50,-50.,50.,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  // Setup Extra Histograms
  NCatFired=HConfig.GetTH1D(Name+"_NCatFired","NCatFired",6,-0.5,5.5,"Num. of passed categories");

  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"N(Vtx) before selection");
  NVtxFullSelection=HConfig.GetTH1D(Name+"_NVtxFullSelection","NVtxFullSelection",26,-0.5,25.5,"N(Vertex) after selection");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"N(good Vertex)");
  VtxZ=HConfig.GetTH1D(Name+"_VtxZ","VtxZ",50,-50.0,50.0,"z(Vtx)/cm");
  VtxRho=HConfig.GetTH1D(Name+"_VtxRho","VtxRho",50,0.0,2.0,"#rho(Vtx)/cm");
  VtxPhi=HConfig.GetTH1D(Name+"_VtxPhi","VtxPhi",50,0.0,3.2,"#phi(Vtx)");
  VtxNdof=HConfig.GetTH1D(Name+"_VtxNdof","VtxNdof",50,-0.5,49.5,"NDoF(Vtx Fit)");
  VtxIsfake=HConfig.GetTH1D(Name+"_VtxIsfake","VtxIsfake",2,-0.5,1.5,"IsFake(Vtx)");

  MuDxy=HConfig.GetTH1D(Name+"_MuDxy","MuDxy",60,-0.3,0.3,"d_{xy}(#mu,Vtx)/cm");
  MuDz=HConfig.GetTH1D(Name+"_MuDz","MuDz",60,-.6,.6,"d_{z}(#mu,Vtx)/cm");
  MuRelIso=HConfig.GetTH1D(Name+"_MuRelIso","MuRelIso",50,0.,1.,"relIso(#mu)");
  MuPt=HConfig.GetTH1D(Name+"_MuPt","MuPt",50,0.,200.,"p_{T}(#mu)/GeV");
  MuEta=HConfig.GetTH1D(Name+"_MuEta","MuEta",50,-2.5,2.5,"#eta(#mu)");
  MuPhi=HConfig.GetTH1D(Name+"_MuPhi","MuPhi",50,-3.14159,3.14159,"#phi(#mu)");

  MuSelPt=HConfig.GetTH1D(Name+"_MuSelPt","MuSelPt",50,0.,200.,"p_{T}(#mu_{sel})/GeV");
  MuSelEta=HConfig.GetTH1D(Name+"_MuSelEta","MuSelEta",50,-2.5,2.5,"#eta(#mu_{sel})");
  MuSelPhi=HConfig.GetTH1D(Name+"_MuSelPhi","MuSelPhi",50,-3.14159,3.14159,"#phi(#mu_{sel})");
  MuSelFakesTauID=HConfig.GetTH1D(Name+"_MuSelFakesTauID","MuSelFakesTauID",2,-0.5,1.5,"#mu_{sel} fakes #tau_{h}");

  TauPt=HConfig.GetTH1D(Name+"_TauPt","TauPt",50,0.,200.,"p_{T}(#tau)/GeV");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",50,-2.5,2.5,"#eta(#tau)");
  TauPhi=HConfig.GetTH1D(Name+"_TauPhi","TauPhi",50,-3.14159,3.14159,"#phi(#tau)");

  TauSelPt=HConfig.GetTH1D(Name+"_TauSelPt","TauSelPt",50,0.,200.,"p_{T}(#tau_{sel})/GeV");
  TauSelEta=HConfig.GetTH1D(Name+"_TauSelEta","TauSelEta",50,-2.5,2.5,"#eta(#tau_{sel})");
  TauSelPhi=HConfig.GetTH1D(Name+"_TauSelPhi","TauSelPhi",50,-3.14159,3.14159,"#phi(#tau_{sel})");

  MuVetoDPtSelMuon=HConfig.GetTH1D(Name+"_MuVetoDPtSelMuon","MuVetoDPtSelMuon",100,-100.,100.,"#Deltap_{T}(#mu_{veto},#mu)/GeV");
  MuVetoInvM=HConfig.GetTH1D(Name+"_MuVetoInvM","MuVetoInvM",100,0.,200,"m_{inv}(#mu_{veto}^{1},#mu_{veto}^{2})/GeV");
  MuVetoPtPositive=HConfig.GetTH1D(Name+"_MuVetoPtPositive","MuVetoPtPositive",50,0.,200.,"p_{T}(#mu_{veto}^{+})/GeV");
  MuVetoPtNegative=HConfig.GetTH1D(Name+"_MuVetoPtNegative","MuVetoPtNegative",50,0.,200.,"p_{T}(#mu_{veto}^{-})/GeV");
  MuVetoDRTau=HConfig.GetTH1D(Name+"_MuVetoDRTau","MuVetoDRTau",50,0.,5.,"#DeltaR(#mu_{veto},#tau_{h})");

  NMuonTriLepVeto=HConfig.GetTH1D(Name+"_NMuonTriLepVeto","NMuonTriLepVeto",5,-0.5,4.5,"N(#mu_{3l veto})");
  NElecTriLepVeto=HConfig.GetTH1D(Name+"_NElecTriLepVeto","NElecTriLepVeto",5,-0.5,4.5,"N(e_{3l veto})");

  MuCharge=HConfig.GetTH1D(Name+"_MuCharge","MuCharge",3,-1.5,1.5,"q(#mu)/e");
  TauCharge=HConfig.GetTH1D(Name+"_TauCharge","TauCharge",3,-1.5,1.5,"q(#tau)/e");

  MuTauDR=HConfig.GetTH1D(Name+"_MuTauDR","MuTauDR",50,0.,5.,"#DeltaR(#mu,#tau_{h})");
  MuTauDPhi=HConfig.GetTH1D(Name+"_MuTauDPhi","MuTauDPhi",50,0.,3.2,"#Delta#phi(#mu,#tau_{h})");
  MuTauDEta=HConfig.GetTH1D(Name+"_MuTauDEta","MuTauDEta",100,-6.,6.,"#Delta#eta(#mu,#tau_{h})");
  MuTauDPt=HConfig.GetTH1D(Name+"_MuTauDPt","MuTauDPt",100,-100.,100.,"#Deltap_{T}(#mu,#tau_{h})/GeV");
  MuTauRelDPt=HConfig.GetTH1D(Name+"_MuTauRelDPt","MuTauRelDPt",100,-2.,2.,"#Deltap_{T}(#mu,#tau_{h})/p_{T}(#mu)");

  MetPt  = HConfig.GetTH1D(Name+"_MetPt","MetPt",50,0.,200.,"E_{T}^{miss}/GeV");
  MetPhi = HConfig.GetTH1D(Name+"_MetPhi","MetPhi",50,-3.14159,3.14159,"#phi(E_{T}^{miss})");

  NJetsKin = HConfig.GetTH1D(Name+"_NJetsKin","NJetsKin",11,-0.5,10.5,"N(j_{kin})");
  JetKin1Pt = HConfig.GetTH1D(Name+"_JetKin1Pt","JetKin1Pt",50,0.,200.,"p_{T}(j_{kin}^{1})/GeV");
  JetKin1Eta = HConfig.GetTH1D(Name+"_JetKin1Eta","JetKin1Eta",100,-5.0,5.0,"#eta(j_{kin}^{1})");
  JetKin1Phi = HConfig.GetTH1D(Name+"_JetKin1Phi","JetKin1Phi",50,-3.14159,3.14159,"#phi(j_{kin}^{1})");
  JetKin1IsLooseId = HConfig.GetTH1D(Name+"_JetKin1IsLooseId","JetKin1IsLooseId",2,-0.5,1.5,"isLoosePUJetID(j_{kin}^{1}");
  JetKin2IsLooseId = HConfig.GetTH1D(Name+"_JetKin2IsLooseId","JetKin2IsLooseId",2,-0.5,1.5,"isLoosePUJetID(j_{kin}^{2}");
  JetKin2Pt = HConfig.GetTH1D(Name+"_JetKin2Pt","JetKin2Pt",50,0.,200.,"p_{T}(j_{kin}^{2})/GeV");
  JetKin2Eta = HConfig.GetTH1D(Name+"_JetKin2Eta","JetKin2Eta",100,-5.0,5.0,"#eta(j_{kin}^{2})");
  JetKin2Phi = HConfig.GetTH1D(Name+"_JetKin2Phi","JetKin2Phi",50,-3.14159,3.14159,"#phi(j_{kin}^{2})");
  NJetsId = HConfig.GetTH1D(Name+"_NJetsId","NJetsId",11,-0.5,10.5,"N(jets)");
  Jet1Pt = HConfig.GetTH1D(Name+"_Jet1Pt","Jet1Pt",50,0.,200.,"p_{T}(j^{1})/GeV");
  Jet1Eta = HConfig.GetTH1D(Name+"_Jet1Eta","Jet1Eta",100,-5.0,5.0,"#eta(j^{1})");
  Jet1Phi = HConfig.GetTH1D(Name+"_Jet1Phi","Jet1Phi",50,-3.14159,3.14159,"#phi(j^{1})");
  Jet1IsB = HConfig.GetTH1D(Name+"_Jet1IsB","Jet1IsB",2,-0.5,1.5,"isBJet(j^{1})");
  Jet2Pt = HConfig.GetTH1D(Name+"_Jet2Pt","Jet2Pt",50,0.,200.,"p_{T}(j^{2})/GeV");
  Jet2Eta = HConfig.GetTH1D(Name+"_Jet2Eta","Jet2Eta",100,-5.0,5.0,"#eta(j^{2})");
  Jet2Phi = HConfig.GetTH1D(Name+"_Jet2Phi","Jet2Phi",50,-3.14159,3.14159,"#phi(j^{2})");
  Jet2IsB = HConfig.GetTH1D(Name+"_Jet2IsB","Jet2IsB",2,-0.5,1.5,"isBJet(j^{2})");

  NBJets = HConfig.GetTH1D(Name+"_NBJets","NBJets",6,-0.5,5.5,"N(bjets)");
  BJet1Pt = HConfig.GetTH1D(Name+"_BJet1Pt","BJet1Pt",50,0.,200.,"p_{T}(b^{1})/GeV");
  BJet1Eta = HConfig.GetTH1D(Name+"_BJet1Eta","BJet1Eta",100,-5.0,5.0,"#eta(b^{1})");
  BJet1Phi = HConfig.GetTH1D(Name+"_BJet1Phi","BJet1Phi",50,-3.14159,3.14159,"#phi(b^{1})");

  // configure category
  if (categoryFlag == "VBF")		configure_VBF();
  if (categoryFlag == "OneJetHigh")	configure_OneJetHigh();
  if (categoryFlag == "OneJetLow")	configure_OneJetLow();
  if (categoryFlag == "ZeroJetHigh")configure_ZeroJetHigh();
  if (categoryFlag == "ZeroJetLow") configure_ZeroJetLow();
  if (categoryFlag == "NoCategory")	configure_NoCategory();

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}


void  HToTaumuTauh::Store_ExtraDist(){
 Extradist1d.push_back(&NCatFired);

 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NVtxFullSelection);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&VtxZ);
 Extradist1d.push_back(&VtxRho);
 Extradist1d.push_back(&VtxNdof);
 Extradist1d.push_back(&VtxIsfake);

 Extradist1d.push_back(&MuDxy);
 Extradist1d.push_back(&MuDz );
 Extradist1d.push_back(&MuRelIso);
 Extradist1d.push_back(&MuPt  );
 Extradist1d.push_back(&MuEta  );
 Extradist1d.push_back(&MuPhi  );

 Extradist1d.push_back(&MuSelPt  );
 Extradist1d.push_back(&MuSelEta  );
 Extradist1d.push_back(&MuSelPhi  );
 Extradist1d.push_back(&MuSelFakesTauID  );

 Extradist1d.push_back(&TauPt  );
 Extradist1d.push_back(&TauEta  );
 Extradist1d.push_back(&TauPhi  );

 Extradist1d.push_back(&TauSelPt  );
 Extradist1d.push_back(&TauSelEta  );
 Extradist1d.push_back(&TauSelPhi  );

 Extradist1d.push_back(&MuVetoDPtSelMuon);
 Extradist1d.push_back(&MuVetoInvM);
 Extradist1d.push_back(&MuVetoPtPositive);
 Extradist1d.push_back(&MuVetoPtNegative);
 Extradist1d.push_back(&MuVetoDRTau);

 Extradist1d.push_back(&NMuonTriLepVeto);
 Extradist1d.push_back(&NElecTriLepVeto);

 Extradist1d.push_back(&MuCharge  );
 Extradist1d.push_back(&TauCharge  );

 Extradist1d.push_back(&MuTauDR);
 Extradist1d.push_back(&MuTauDPhi);
 Extradist1d.push_back(&MuTauDEta);
 Extradist1d.push_back(&MuTauDPt);
 Extradist1d.push_back(&MuTauRelDPt);

 Extradist1d.push_back(&MetPt);
 Extradist1d.push_back(&MetPhi);

 Extradist1d.push_back(&NJetsKin);
 Extradist1d.push_back(&JetKin1Pt);
 Extradist1d.push_back(&JetKin1Eta);
 Extradist1d.push_back(&JetKin1Phi);
 Extradist1d.push_back(&JetKin1IsLooseId);
 Extradist1d.push_back(&JetKin2IsLooseId);
 Extradist1d.push_back(&JetKin2Pt);
 Extradist1d.push_back(&JetKin2Eta);
 Extradist1d.push_back(&JetKin2Phi);
 Extradist1d.push_back(&NJetsId);
 Extradist1d.push_back(&Jet1Pt);
 Extradist1d.push_back(&Jet1Eta);
 Extradist1d.push_back(&Jet1Phi);
 Extradist1d.push_back(&Jet1IsB);
 Extradist1d.push_back(&Jet2Pt);
 Extradist1d.push_back(&Jet2Eta);
 Extradist1d.push_back(&Jet2Phi);
 Extradist1d.push_back(&Jet2IsB);

 Extradist1d.push_back(&NBJets);
 Extradist1d.push_back(&BJet1Pt);
 Extradist1d.push_back(&BJet1Eta);
 Extradist1d.push_back(&BJet1Phi);
}

void  HToTaumuTauh::doEvent(){
  // set variables to hold selected objects to default values
  selVertex = -1;
  selMuon = -1;
  selTau = -1;
  selKinJets.clear();
  selBJets.clear();
  selMjj = -1;
  selJetdeta = -100;
  selNjetingap = -1;

  unsigned int t;
  int id(Ntp->GetMCID());
  //std::cout << "ID before = " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  double wobs=1;
  double w;
  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}

  // Apply Selection

  // Vertex
  unsigned int nGoodVtx=0;
  for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
    if(Ntp->isGoodVtx(i_vtx)){
    	if(selVertex == -1) selVertex = i_vtx; // selected vertex = first vertex (highest sum[pT^2]) to fulfill vertex requirements
    	nGoodVtx++;
    }
  }

  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  // Trigger
  value.at(TriggerOk) = -1;
  for (std::vector<TString>::iterator it_trig = cTriggerNames.begin(); it_trig != cTriggerNames.end(); ++it_trig){
	  if(Ntp->TriggerAccept(*it_trig)){
		  if ( value.at(TriggerOk) == -1 )
			  value.at(TriggerOk) = it_trig - cTriggerNames.begin();
		  else // more than 1 trigger fired, save this separately
			  value.at(TriggerOk) = cTriggerNames.size();
	  }
  }
  pass.at(TriggerOk) = (value.at(TriggerOk) >= cut.at(TriggerOk));
  
  // Muon cuts
  std::vector<int> selectedMuonsId;
  selectedMuonsId.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_Id(i_mu,selVertex) ) {
		  selectedMuonsId.push_back(i_mu);
	  }
  }
  value.at(NMuId)=selectedMuonsId.size();
  pass.at(NMuId)=(value.at(NMuId)>=cut.at(NMuId));

  std::vector<int> selectedMuons;	// full selection: ID and Kinematics
  selectedMuons.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsId.begin(); it_mu != selectedMuonsId.end(); ++it_mu){
	  if( selectMuon_Kinematics(*it_mu)) {
		  selectedMuons.push_back(*it_mu);
	  }
  }
  value.at(NMuKin)=selectedMuons.size();
  pass.at(NMuKin)=(value.at(NMuKin)>=cut.at(NMuKin));
  if (selectedMuons.size() > 0) selMuon = selectedMuons.at(0);

  std::vector<int> diMuonVetoMuonsPositive;	// muons selected for the dimuon veto
  diMuonVetoMuonsPositive.clear();
  std::vector<int> diMuonVetoMuonsNegative;	// muons selected for the dimuon veto
  diMuonVetoMuonsNegative.clear();
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if( selectMuon_diMuonVeto(i, selVertex) ) {
		  if (Ntp->Muon_Charge(i) == 1) {
			  diMuonVetoMuonsPositive.push_back(i);
		  }
		  else if (Ntp->Muon_Charge(i) == -1) {
			  diMuonVetoMuonsNegative.push_back(i);
		  }
	  }
  }
  value.at(DiMuonVeto) = (diMuonVetoMuonsPositive.size() >= 1 && diMuonVetoMuonsNegative.size() >= 1);
  pass.at(DiMuonVeto)=(value.at(DiMuonVeto) == cut.at(DiMuonVeto));

  // Tau cuts
  std::vector<int> selectedTausId;
  selectedTausId.clear();
  for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
	  if ( selectPFTau_Id(i_tau,selectedMuonsId) ){
		  selectedTausId.push_back(i_tau);
	  }
  }
  value.at(NTauId)=selectedTausId.size();
  pass.at(NTauId)=(value.at(NTauId)>=cut.at(NTauId));

  std::vector<int> selectedTausIso;
  selectedTausIso.clear();
  for(std::vector<int>::iterator it_tau = selectedTausId.begin(); it_tau != selectedTausId.end(); ++it_tau){
	  if ( selectPFTau_Iso(*it_tau) ){
		  selectedTausIso.push_back(*it_tau);
	  }
  }
  value.at(NTauIso)=selectedTausIso.size();
  pass.at(NTauIso)=(value.at(NTauIso)>=cut.at(NTauIso));

  std::vector<int> selectedTaus;
  selectedTaus.clear();
  for(std::vector<int>::iterator it_tau = selectedTausIso.begin(); it_tau != selectedTausIso.end(); ++it_tau){
	  if ( selectPFTau_Kinematics(*it_tau) ){
		  selectedTaus.push_back(*it_tau);
	  }
  }
  value.at(NTauKin)=selectedTaus.size();
  pass.at(NTauKin)=(value.at(NTauKin)>=cut.at(NTauKin));
  if(selectedTaus.size() > 0) selTau = selectedTaus.at(0);

  // Tri-lepton veto
  std::vector<int> triLepVetoMuons;
  triLepVetoMuons.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_triLeptonVeto(i_mu,selMuon,selVertex) ) {
		  triLepVetoMuons.push_back(i_mu);
	  }
  }
  std::vector<int> triLepVetoElecs;
  triLepVetoElecs.clear();
  for(unsigned i_el=0;i_el<Ntp->NElectrons();i_el++){
	  if( selectElectron_triLeptonVeto(i_el,selVertex,selectedMuonsId) ) {
		  triLepVetoElecs.push_back(i_el);
	  }
  }
  value.at(TriLeptonVeto) = triLepVetoMuons.size() + triLepVetoElecs.size();
  pass.at(TriLeptonVeto) = (value.at(TriLeptonVeto) <= cut.at(TriLeptonVeto));

  // Opposite charge
  if (selMuon != -1 && selTau != -1){
	  value.at(OppCharge) = Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau);
  }
  else value.at(OppCharge) = -9;
  if (cut.at(OppCharge) == 999) // set to 999 to disable oppcharge cut
	  pass.at(OppCharge) = true;
  else
	  pass.at(OppCharge) = (value.at(OppCharge) == cut.at(OppCharge));

  // Transverse mass
  if(selMuon == -1){ // no good muon in event: set MT to high value -> fail MT cut
	  value.at(MT) = 500.0;
  }
  else{
	  double pT 	= Ntp->Muon_p4(selMuon).Pt();
	  double phi	= Ntp->Muon_p4(selMuon).Phi();
	  double eTmiss = Ntp->MET_CorrMVA_et();
	  double eTmPhi = Ntp->MET_CorrMVA_phi();
	  value.at(MT)	= transverseMass(pT,phi,eTmiss,eTmPhi);
  }
  if (cut.at(MT) == 999) // set to 999 to disable mt cut
	  pass.at(MT) = true;
  else
	  pass.at(MT) = (value.at(MT) < cut.at(MT));

  // sort jets by corrected pt
  std::vector<int> sortedPFJets = sortPFjets();
//  bool print = false;
//  for (unsigned j = 0; j < Ntp->NPFJets(); j++){
//	  if (j != sortedPFJets.at(j)) print = true;
//  }
//  if (print){
//	  printf("There are %d unsorted and %d sorted jets.\n", Ntp->NPFJets(), sortedPFJets.size());
//	  for (unsigned j = 0; j < Ntp->NPFJets(); j++){
//		  printf("   jet %2d with pt= %4.2f       jet %2d with pt= %4.2f \n", j, Ntp->PFJet_p4(j).Pt(), sortedPFJets.at(j), Ntp->PFJet_p4(sortedPFJets.at(j)).Pt());
//	  }
//  }
  // select jets for categories
  // PFJet and bjet collections can have mutual elements!
  std::vector<int> selectedJetsKin;
  selectedJetsKin.clear();
  std::vector<int> selectedJets;
  selectedJets.clear();
  std::vector<int> selectedBJets;
  selectedBJets.clear();
  for (unsigned i_jet = 0; i_jet < Ntp->NPFJets(); i_jet++){
	  if ( selectPFJet_Kinematics(sortedPFJets.at(i_jet), selMuon, selTau) ) {
		  selectedJetsKin.push_back(sortedPFJets.at(i_jet));
		  if ( selectPFJet_Id(sortedPFJets.at(i_jet)) ){
			  selectedJets.push_back(sortedPFJets.at(i_jet));
		  }
	  }
	  if ( selectBJet(sortedPFJets.at(i_jet), selMuon, selTau) ) {
		  selectedBJets.push_back(sortedPFJets.at(i_jet));
	  }
  }
  selKinJets = selectedJetsKin;
  selBJets = selectedBJets;

  // run categories
  bool passed_VBF 			= category_VBF(selectedJetsKin, selectedBJets);
  bool passed_OneJetHigh	= category_OneJetHigh(selTau, selectedJetsKin, selectedBJets, passed_VBF);
  bool passed_OneJetLow		= category_OneJetLow(selTau, selectedJetsKin, selectedBJets, passed_VBF);
  bool passed_ZeroJetHigh	= category_ZeroJetHigh(selTau,selectedJetsKin, selectedBJets);
  bool passed_ZeroJetLow	= category_ZeroJetLow(selTau,selectedJetsKin, selectedBJets);
  bool passed_NoCategory	= category_NoCategory();

  // fill plot checking if multiple categories have passed, which should never happen
  unsigned nCat = 0;
  if (passed_VBF 		) nCat++;
  if (passed_OneJetHigh	) nCat++;
  if (passed_OneJetLow	) nCat++;
  if (passed_ZeroJetHigh) nCat++;
  if (passed_ZeroJetLow	) nCat++;

  NCatFired.at(t).Fill(nCat);



  bool status=AnalysisCuts(t,w,wobs); // true only if full selection passed

  ///////////////////////////////////////////////////////////
  // Add plots
  ///////////////////////////////////////////////////////////

  //////// plots filled before any cuts
  // Vertex plots
  NVtx.at(t).Fill(Ntp->NVtx(),w);
  for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
	VtxZ.at(t).Fill(Ntp->Vtx(i_vtx).z(),w);
	VtxRho.at(t).Fill(sqrt(Ntp->Vtx(i_vtx).x()*Ntp->Vtx(i_vtx).x() + Ntp->Vtx(i_vtx).y()*Ntp->Vtx(i_vtx).y()), w);
	VtxNdof.at(t).Fill(Ntp->Vtx_ndof(i_vtx), w);
	VtxIsfake.at(t).Fill(Ntp->Vtx_isFake(i_vtx), w);
  }
  NGoodVtx.at(t).Fill(nGoodVtx,w);

  //////// plots filled after Vertex selection: Object selection

  if(pass.at(TriggerOk) && pass.at(PrimeVtx)){
	  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
		  if(	Ntp->isTightMuon(i_mu,selVertex) ){
			  MuDxy.at(t).Fill(dxy(Ntp->Muon_p4(i_mu),Ntp->Muon_Poca(i_mu),Ntp->Vtx(selVertex)), w);
			  MuDz.at(t).Fill(dz(Ntp->Muon_p4(i_mu),Ntp->Muon_Poca(i_mu),Ntp->Vtx(selVertex)), w);
			  MuRelIso.at(t).Fill(Ntp->Muon_RelIso(i_mu), w);
		  }
	  }
  }

  //////// plots filled after muon ID selection: Muon Kinematics
  if(pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NMuId)){
	  for(std::vector<int>::iterator it_mu = selectedMuonsId.begin();it_mu != selectedMuonsId.end(); ++it_mu){
		  MuPt.at(t).Fill(Ntp->Muon_p4(*it_mu).Pt(), w);
		  MuEta.at(t).Fill(Ntp->Muon_p4(*it_mu).Eta(), w);
		  MuPhi.at(t).Fill(Ntp->Muon_p4(*it_mu).Phi(), w);
	  }

	  //////// plots filled only with selected muon
	  if(pass.at(NMuKin)){
		  MuSelPt.at(t).Fill(Ntp->Muon_p4(selMuon).Pt(), w);
		  MuSelEta.at(t).Fill(Ntp->Muon_p4(selMuon).Eta(), w);
		  MuSelPhi.at(t).Fill(Ntp->Muon_p4(selMuon).Phi(), w);

		  // Does the muon fake the tau_ID+Iso?
		  bool fakes = false;
		  for( unsigned  i_tau = 0; i_tau < Ntp->NPFTaus(); i_tau++){
			  if (	  selectPFTau_Id(i_tau) &&
					  selectPFTau_Iso(i_tau) &&
					  Ntp->Muon_p4(selMuon).DeltaR(Ntp->PFTau_p4(i_tau)) < cMuTau_dR){
				  fakes = true;
				  break;
			  }
		  }
		  MuSelFakesTauID.at(t).Fill(fakes, w);
	  }
  }

  //////// plots filled after tau ID + Iso selection: Tau Kinematics
  if(pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NTauId) && pass.at(NTauIso)){
	  for(std::vector<int>::iterator it_tau = selectedTausIso.begin(); it_tau != selectedTausIso.end(); ++it_tau){
		  TauPt.at(t).Fill(Ntp->PFTau_p4(*it_tau).Pt(), w);
		  TauEta.at(t).Fill(Ntp->PFTau_p4(*it_tau).Eta(), w);
		  TauPhi.at(t).Fill(Ntp->PFTau_p4(*it_tau).Phi(), w);
	  }

	  //////// plots filled only with selected tau
	  if(pass.at(NTauKin)){
		  TauSelPt.at(t).Fill(Ntp->PFTau_p4(selTau).Pt(), w);
		  TauSelEta.at(t).Fill(Ntp->PFTau_p4(selTau).Eta(), w);
		  TauSelPhi.at(t).Fill(Ntp->PFTau_p4(selTau).Phi(), w);
	  }
  }

  //////// plots filled after full muon and tau selection
  bool passedObjects = pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NMuId   ) && pass.at(NMuKin  ) && pass.at(NTauId  ) && pass.at(NTauIso ) && pass.at(NTauKin );
  if(passedObjects && !pass.at(DiMuonVeto)){
	  // Investigate events discarded by the DiMuon Veto
	  if (Ntp->Muon_Charge(selMuon) == 1){
		  MuVetoDPtSelMuon.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)).Pt() - Ntp->Muon_p4(selMuon).Pt(), w );
		  MuVetoDRTau.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)).DeltaR(Ntp->PFTau_p4(selTau)), w);
	  }
	  else if (Ntp->Muon_Charge(selMuon) == -1){
		  MuVetoDPtSelMuon.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).Pt() - Ntp->Muon_p4(selMuon).Pt(), w );
		  MuVetoDRTau.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).DeltaR(Ntp->PFTau_p4(selTau)), w);
	  }
	  MuVetoInvM.at(t).Fill( (Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)) + Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0))).M() , w);
	  MuVetoPtPositive.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).Pt(), w);
	  MuVetoPtNegative.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)).Pt(), w);
  }

  if(passedObjects && pass.at(DiMuonVeto)){
	  // Tri-lepton vetoes
	  NMuonTriLepVeto.at(t).Fill(triLepVetoMuons.size(), w);
	  NElecTriLepVeto.at(t).Fill(triLepVetoElecs.size(), w);
  }

  //////// plots filled after full selection (without categories)
  bool fullSel = passedObjects && pass.at(DiMuonVeto) && pass.at(TriLeptonVeto) && pass.at(OppCharge) && pass.at(MT);
  if (fullSel){
	  // Mu-Tau correlations
	  MuTauDR    .at(t).Fill( Ntp->Muon_p4(selMuon).DeltaR(Ntp->PFTau_p4(selTau)), w );
	  MuTauDPhi  .at(t).Fill( Ntp->Muon_p4(selMuon).DeltaPhi(Ntp->PFTau_p4(selTau)), w );
	  MuTauDEta  .at(t).Fill( Ntp->Muon_p4(selMuon).Eta() - Ntp->PFTau_p4(selTau).Eta(), w );
	  MuTauDPt   .at(t).Fill( Ntp->Muon_p4(selMuon).Pt() - Ntp->PFTau_p4(selTau).Pt(), w );
	  MuTauRelDPt.at(t).Fill( (Ntp->Muon_p4(selMuon).Pt() - Ntp->PFTau_p4(selTau).Pt()) / Ntp->Muon_p4(selMuon).Pt() , w);

	  // lepton charge
	  MuCharge.at(t).Fill( Ntp->Muon_Charge(selMuon), w);
	  TauCharge.at(t).Fill( Ntp->PFTau_Charge(selMuon), w);

	  // MET
	  MetPt.at(t).Fill( Ntp->MET_CorrMVA_et(), w);
	  MetPhi.at(t).Fill( Ntp->MET_CorrMVA_phi(), w);

	  // Jets
	  NJetsKin.at(t).Fill( selectedJetsKin.size(), w);
	  if (selectedJetsKin.size() > 0){
		  JetKin1Pt.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(0)).Pt(), w);
		  JetKin1Eta.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(0)).Eta(), w);
		  JetKin1Phi.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(0)).Phi(), w);
		  JetKin1IsLooseId.at(t).Fill( Ntp->PFJet_PUJetID_looseWP(selectedJetsKin.at(0)), w);
	  }
	  if (selectedJetsKin.size() > 1){
		  JetKin2IsLooseId.at(t).Fill( Ntp->PFJet_PUJetID_looseWP(selectedJetsKin.at(1)), w);
		  JetKin2Pt.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(1)).Pt(), w);
		  JetKin2Eta.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(1)).Eta(), w);
		  JetKin2Phi.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(1)).Phi(), w);
	  }
	  NJetsId.at(t).Fill( selectedJets.size(), w);
	  if (selectedJets.size() > 0){
		  Jet1Pt.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(0)).Pt(), w);
		  Jet1Eta.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(0)).Eta(), w);
		  Jet1Phi.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(0)).Phi(), w);
	  	  Jet1IsB.at(t).Fill( Ntp->PFJet_bDiscriminator(selectedJets.at(0)) > cCat_btagDisc , w);
	  }
	  if (selectedJets.size() > 1){
		  Jet2Pt.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(1)).Pt(), w);
		  Jet2Eta.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(1)).Eta(), w);
		  Jet2Phi.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(1)).Phi(), w);
		  Jet2IsB.at(t).Fill( Ntp->PFJet_bDiscriminator(selectedJets.at(1)) > cCat_btagDisc, w);
	  }

	  NBJets.at(t).Fill( selectedBJets.size(), w);
	  if (selectedBJets.size() > 0){
		  BJet1Pt.at(t).Fill( Ntp->PFJet_p4(selectedBJets.at(0)).Pt(), w);
		  BJet1Eta.at(t).Fill( Ntp->PFJet_p4(selectedBJets.at(0)).Eta(), w);
		  BJet1Phi.at(t).Fill( Ntp->PFJet_p4(selectedBJets.at(0)).Phi(), w);
	  }
  }

  //////// plots filled after full selection
  if(status){
    NVtxFullSelection.at(t).Fill(Ntp->NVtx(),w);
    //std::cout << "ID after = " << id << std::endl;
  }
}


void  HToTaumuTauh::Finish(){
  Selection::Finish();
}


/////////////////////////////////////////
// Definition of selection and helper functions
/////////////////////////////////////////

///////// Helper functions
double HToTaumuTauh::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return (-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt();
}

double HToTaumuTauh::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2);
}

///////// Muons

bool HToTaumuTauh::selectMuon_Id(unsigned i, unsigned vertex){
	if(	Ntp->isSelectedMuon(i,vertex,cMu_dxy,cMu_dz) &&
		Ntp->Muon_RelIso(i) < cMu_relIso
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_Kinematics(unsigned i){
	if(	Ntp->Muon_p4(i).Pt() >= cMu_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) <= cMu_eta
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_diMuonVeto(unsigned i, unsigned i_vtx){
	if(	Ntp->Muon_p4(i).Pt() > 15.0 &&
		fabs(Ntp->Muon_p4(i).Eta()) < 2.4 &&
		Ntp->Muon_isPFMuon(i) &&
		Ntp->Muon_isGlobalMuon(i) &&
		Ntp->Muon_isTrackerMuon(i) &&
		Ntp->Muon_RelIso(i) < 0.3 &&
		fabs(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < 0.2
		) {
	  return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_triLeptonVeto(unsigned i, int selectedMuon, unsigned i_vtx){
	if(	i != selectedMuon &&
		Ntp->isTightMuon(i,i_vtx) &&
		fabs(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < cMu_dxy &&
		fabs(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < cMu_dz &&
		Ntp->Muon_RelIso(i) < 0.3 &&
		Ntp->Muon_p4(i).Pt() > cMuTriLep_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) < cMuTriLep_eta
			){
			return true;
	}
	return false;
}


///////// Electrons
bool HToTaumuTauh::selectElectron_triLeptonVeto(unsigned i, unsigned i_vtx, std::vector<int> muonCollection){
	// check if elec is matched to a muon, if so this is not a good elec (should become obsolete when using top projections)
	for(std::vector<int>::iterator it_mu = muonCollection.begin(); it_mu != muonCollection.end(); ++it_mu){
	  if( Ntp->Electron_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR ) {
		  return false;
	  }
	}

	if ( 	Ntp->isSelectedElectron(i,i_vtx,0.045,0.2) &&
			//TODO: This electron isolation is using rho corrections, but should use deltaBeta corrections
			//documentation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
			Ntp->Electron_RelIso04(i) < 0.3 &&
			Ntp->Electron_p4(i).Pt() > 10.0 &&
			fabs(Ntp->Electron_p4(i).Eta()) < 2.5
			){
		return true;
	}
	return false;
}

///////// Taus

bool HToTaumuTauh::selectPFTau_Id(unsigned i){
	if ( 	Ntp->PFTau_isHPSByDecayModeFinding(i) &&
			Ntp->PFTau_isHPSAgainstElectronsLoose(i) &&
			Ntp->PFTau_isHPSAgainstMuonTight(i)
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Id(unsigned i, std::vector<int> muonCollection){
	// check if tau is matched to a muon, if so this is not a good tau
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Sync_Issues
	for(std::vector<int>::iterator it_mu = muonCollection.begin(); it_mu != muonCollection.end(); ++it_mu){
	  if( Ntp->PFTau_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR ) {
		  return false;
	  }
	}
	if ( 	selectPFTau_Id(i) ){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Iso(unsigned i){
	if ( 	Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(i) < cTau_rawIso
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Kinematics(unsigned i){
	if ( 	Ntp->PFTau_p4(i).Pt() >= cTau_pt &&
			fabs(Ntp->PFTau_p4(i).Eta()) <= cTau_eta
			){
		return true;
	}
	return false;
}

std::vector<int> HToTaumuTauh::sortPFjets(){
	// create vector of pairs to allow for sorting by jet pt
	std::vector< std::pair<int, double> > jetIdxPt;
	for (unsigned i = 0; i<Ntp->NPFJets(); i++ ){
		jetIdxPt.push_back( std::make_pair(i,Ntp->PFJet_p4(i).Pt()) );
	}
	// sort vector of pairs
	std::sort(jetIdxPt.begin(), jetIdxPt.end(), sortIdxByValue());
	// create vector of indices in correct order
	std::vector<int> sortedJets;
	for (unsigned i = 0; i<jetIdxPt.size(); i++){
		sortedJets.push_back(jetIdxPt.at(i).first);
	}
	return sortedJets;
}

bool HToTaumuTauh::selectPFJet_Kinematics(unsigned i, int selectedMuon, int selectedTau){
	// clean against selected muon and tau
	if (selectedMuon >= 0) {
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(selectedMuon)) < cJetClean_dR) return false;
	}
	if (selectedTau >= 0){
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->PFTau_p4(selectedTau)) < cJetClean_dR) return false;
	}

	if ( 	fabs(Ntp->PFJet_p4(i).Eta()) < cCat_jetEta &&
			Ntp->PFJet_p4(i).Pt() > cCat_jetPt){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFJet_Id(unsigned i){
	if (	Ntp->PFJet_PUJetID_looseWP(i) ){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectBJet(unsigned i, int selectedMuon, int selectedTau){
	// clean against selected muon and tau
	if (selectedMuon >= 0) {
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(selectedMuon)) < cJetClean_dR) return false;
	}
	if (selectedTau >= 0){
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->PFTau_p4(selectedTau)) < cJetClean_dR) return false;
	}

	if (	fabs(Ntp->PFJet_p4(i).Eta()) < cCat_bjetEta &&
			Ntp->PFJet_p4(i).Pt() > cCat_bjetPt &&
			Ntp->PFJet_bDiscriminator(i) > cCat_btagDisc){
		return true;
	}
	return false;
}


//// *****functions defining the categories*****

void HToTaumuTauh::configure_VBF(){
	// to be called only if VBF is chosen category

	// set cut values to be the cut values of this category
	cut.at(VbfNJet) = cut_VBF.at(VbfNJet);
	cut.at(VbfDeltaEta) = cut_VBF.at(VbfDeltaEta);
	cut.at(VbfNJetRapGap) = cut_VBF.at(VbfNJetRapGap);
	cut.at(VbfJetInvM) = cut_VBF.at(VbfJetInvM);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(VbfNJet)="Number VBF Jets $>=$";
	title.at(VbfNJet)+=cut.at(VbfNJet);
	htitle=title.at(VbfNJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF}";
	c="_Cut_";c+=VbfNJet;
	Nminus1.at(VbfNJet) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfNJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(VbfNJet) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfNJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(VbfDeltaEta)="$\\Delta\\eta(jj) >$";
	title.at(VbfDeltaEta)+=cut.at(VbfDeltaEta);
	htitle=title.at(VbfDeltaEta);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="#Delta#eta(Jet_{VBF}^{1},Jet_{VBF}^{2})";
	c="_Cut_";c+=VbfDeltaEta;
	Nminus1.at(VbfDeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfDeltaEta_",htitle,50,-10.,10.,hlabel,"Events");
	Nminus0.at(VbfDeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus0__VbfDeltaEta",htitle,50,-10.,10.,hlabel,"Events");

	title.at(VbfNJetRapGap)="Number Jets in $\\eta$ gap $<=$";
	title.at(VbfNJetRapGap)+=cut.at(VbfNJetRapGap);
	htitle=title.at(VbfNJetRapGap);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF} in rapidity gap";
	c="_Cut_";c+=VbfNJetRapGap;
	Nminus1.at(VbfNJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfNJetRapGap_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(VbfNJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfNJetRapGap_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(VbfJetInvM)="$m_{jj}(VBF) >$";
	title.at(VbfJetInvM)+=cut.at(VbfJetInvM);
	title.at(VbfJetInvM)+=" GeV";
	htitle=title.at(VbfJetInvM);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="m_{inv}(jj) of VBF-jets";
	c="_Cut_";c+=VbfJetInvM;
	Nminus1.at(VbfJetInvM) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfJetInvM_",htitle,50,0.,2000.,hlabel,"Events");
	Nminus0.at(VbfJetInvM) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfJetInvM_",htitle,50,0.,2000.,hlabel,"Events");
}
bool HToTaumuTauh::category_VBF(std::vector<int> jetCollection, std::vector<int> bJetCollection){
	std::vector<float> value_VBF;
	std::vector<float> pass_VBF;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_VBF.push_back(-10.);
	pass_VBF.push_back(false);
	}

	value_VBF.at(VbfNJet) = jetCollection.size();
	pass_VBF.at(VbfNJet) = (value_VBF.at(VbfNJet) >= cut_VBF.at(VbfNJet));

	if(pass_VBF.at(VbfNJet)){
		double vbfJetEta1 = Ntp->PFJet_p4(jetCollection.at(0)).Eta();
		double vbfJetEta2 = Ntp->PFJet_p4(jetCollection.at(1)).Eta();

		value_VBF.at(VbfDeltaEta) = vbfJetEta1 - vbfJetEta2;
		selJetdeta = value_VBF.at(VbfDeltaEta);
		pass_VBF.at(VbfDeltaEta) = (fabs(value_VBF.at(VbfDeltaEta)) > cut_VBF.at(VbfDeltaEta));

		int jetsInRapidityGap = 0;
		for(std::vector<int>::iterator it_jet = jetCollection.begin()+2; it_jet != jetCollection.end(); ++it_jet){
			double etaPos = ( value_VBF.at(VbfDeltaEta) >= 0) ? vbfJetEta1 : vbfJetEta2;
			double etaNeg = ( value_VBF.at(VbfDeltaEta) >= 0) ? vbfJetEta2 : vbfJetEta1;
			  if (	Ntp->PFJet_p4(*it_jet).Eta() > etaNeg &&
					Ntp->PFJet_p4(*it_jet).Eta() < etaPos){
				  jetsInRapidityGap++;
			  }
		}
		value_VBF.at(VbfNJetRapGap) = jetsInRapidityGap;
		selNjetingap = int(value_VBF.at(VbfNJetRapGap));
		pass_VBF.at(VbfNJetRapGap) = (value_VBF.at(VbfNJetRapGap) <= cut_VBF.at(VbfNJetRapGap));

		double invM = (Ntp->PFJet_p4(jetCollection.at(0)) + Ntp->PFJet_p4(jetCollection.at(1))).M();
		value_VBF.at(VbfJetInvM) = invM;
		selMjj = invM;
		pass_VBF.at(VbfJetInvM) = (value_VBF.at(VbfJetInvM) > cut_VBF.at(VbfJetInvM));
	}

	// migrate into main analysis if this is chosen category
	bool catPassed = true;
	for (unsigned i_cut = VbfNJet; i_cut < NCuts; i_cut++){
		if (categoryFlag == "VBF"){
			if (i_cut < VbfNCuts){
				value.at(i_cut) = value_VBF.at(i_cut);
				pass.at(i_cut) = pass_VBF.at(i_cut);
			}
			else{
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		catPassed = catPassed && pass_VBF.at(i_cut);
	}

	return catPassed;
}

void HToTaumuTauh::configure_OneJetHigh(){
	// to be called only if OneJetHigh is chosen category

	// set cut values to be the cut values of this category
	cut.at(OneJetNJet) = cut_OneJet.at(OneJetNJet);
	cut.at(OneJetNoVBF) = cut_OneJet.at(OneJetNoVBF);
	cut.at(OneJetNBtagJets) = cut_OneJet.at(OneJetNBtagJets);
	cut.at(OneJetTauPt) = cut_OneJet.at(OneJetTauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(OneJetNJet)="Number Jets $>=$";
	title.at(OneJetNJet)+=cut.at(OneJetNJet);
	htitle=title.at(OneJetNJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=OneJetNJet;
	Nminus1.at(OneJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetNoVBF)="No VBF $==$";
	title.at(OneJetNoVBF)+=cut.at(OneJetNoVBF);
	htitle=title.at(OneJetNoVBF);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=OneJetNoVBF;
	Nminus1.at(OneJetNoVBF) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetNoVBF_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(OneJetNoVBF) = HConfig.GetTH1D(Name+c+"_Nminus0__OneJetNoVBF",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(OneJetNBtagJets)="Number btag jets $<=$";
	title.at(OneJetNBtagJets)+=cut.at(OneJetNBtagJets);
	htitle=title.at(OneJetNBtagJets);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of btag jets";
	c="_Cut_";c+=OneJetNBtagJets;
	Nminus1.at(OneJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetTauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(OneJetTauPt)+=cut.at(OneJetTauPt);
	title.at(OneJetTauPt)+=" GeV";
	htitle=title.at(OneJetTauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=OneJetTauPt;
	Nminus1.at(OneJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_OneJetHigh(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection, bool passedVBF){
	std::vector<float> value_OneJetHigh;
	std::vector<float> pass_OneJetHigh;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_OneJetHigh.push_back(-10.);
	pass_OneJetHigh.push_back(false);
	}

	value_OneJetHigh.at(OneJetNJet) = jetCollection.size();
	pass_OneJetHigh.at(OneJetNJet) = ( value_OneJetHigh.at(OneJetNJet) >= cut_OneJet.at(OneJetNJet) );

	value_OneJetHigh.at(OneJetNoVBF) = !passedVBF;
	pass_OneJetHigh.at(OneJetNoVBF) = ( value_OneJetHigh.at(OneJetNoVBF) == cut_OneJet.at(OneJetNoVBF) );

	value_OneJetHigh.at(OneJetNBtagJets) = bJetCollection.size();
	pass_OneJetHigh.at(OneJetNBtagJets) = ( value_OneJetHigh.at(OneJetNBtagJets) <= cut_OneJet.at(OneJetNBtagJets) );

	if (selTau == -1){
		value_OneJetHigh.at(OneJetTauPt) = -10.;
		pass_OneJetHigh.at(OneJetTauPt) = false;
	}
	else{
		value_OneJetHigh.at(OneJetTauPt) = Ntp->PFTau_p4(selTau).Pt();
		pass_OneJetHigh.at(OneJetTauPt) = ( value_OneJetHigh.at(OneJetTauPt) >= cut_OneJet.at(OneJetTauPt) );
	}


	// migrate into main analysis if this is chosen category
	bool catPassed = true;
	for (unsigned i_cut = OneJetNJet; i_cut < NCuts; i_cut++){
		if (categoryFlag == "OneJetHigh"){
			if (i_cut < OneJetNCuts){
				value.at(i_cut) = value_OneJetHigh.at(i_cut);
				pass.at(i_cut) = pass_OneJetHigh.at(i_cut);
			}
			else{
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		catPassed = catPassed && pass_OneJetHigh.at(i_cut);
	}

	return catPassed;
}

void HToTaumuTauh::configure_OneJetLow(){
	// to be called only if OneJetLow is chosen category

	// set cut values to be the cut values of this category
	cut.at(OneJetNJet) = cut_OneJet.at(OneJetNJet);
	cut.at(OneJetNoVBF) = cut_OneJet.at(OneJetNoVBF);
	cut.at(OneJetNBtagJets) = cut_OneJet.at(OneJetNBtagJets);
	cut.at(OneJetTauPt) = cut_OneJet.at(OneJetTauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(OneJetNJet)="Number Jets $>=$";
	title.at(OneJetNJet)+=cut.at(OneJetNJet);
	htitle=title.at(OneJetNJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=OneJetNJet;
	Nminus1.at(OneJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetNoVBF)="No VBF $==$";
	title.at(OneJetNoVBF)+=cut.at(OneJetNoVBF);
	htitle=title.at(OneJetNoVBF);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=OneJetNoVBF;
	Nminus1.at(OneJetNoVBF) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetNoVBF_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(OneJetNoVBF) = HConfig.GetTH1D(Name+c+"_Nminus0__OneJetNoVBF",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(OneJetNBtagJets)="Number btag jets $<=$";
	title.at(OneJetNBtagJets)+=cut.at(OneJetNBtagJets);
	htitle=title.at(OneJetNBtagJets);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of btag jets";
	c="_Cut_";c+=OneJetNBtagJets;
	Nminus1.at(OneJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetTauPt)="$p_{T}(\\tau_{h}) <$";
	title.at(OneJetTauPt)+=cut.at(OneJetTauPt);
	title.at(OneJetTauPt)+=" GeV";
	htitle=title.at(OneJetTauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=OneJetTauPt;
	Nminus1.at(OneJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_OneJetLow(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection, bool passedVBF){
	std::vector<float> value_OneJetLow;
	std::vector<float> pass_OneJetLow;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_OneJetLow.push_back(-10.);
	pass_OneJetLow.push_back(false);
	}

	value_OneJetLow.at(OneJetNJet) = jetCollection.size();
	pass_OneJetLow.at(OneJetNJet) = ( value_OneJetLow.at(OneJetNJet) >= cut_OneJet.at(OneJetNJet) );

	value_OneJetLow.at(OneJetNoVBF) = !passedVBF;
	pass_OneJetLow.at(OneJetNoVBF) = ( value_OneJetLow.at(OneJetNoVBF) == cut_OneJet.at(OneJetNoVBF) );

	value_OneJetLow.at(OneJetNBtagJets) = bJetCollection.size();
	pass_OneJetLow.at(OneJetNBtagJets) = ( value_OneJetLow.at(OneJetNBtagJets) <= cut_OneJet.at(OneJetNBtagJets) );

	if (selTau == -1){
		value_OneJetLow.at(OneJetTauPt) = -10.;
		pass_OneJetLow.at(OneJetTauPt) = false;
	}
	else{
		value_OneJetLow.at(OneJetTauPt) = Ntp->PFTau_p4(selTau).Pt();
		pass_OneJetLow.at(OneJetTauPt) = ( value_OneJetLow.at(OneJetTauPt) >= cut_OneJet.at(OneJetTauPt) );
	}


	// migrate into main analysis if this is chosen category
	bool catPassed = true;
	for (unsigned i_cut = OneJetNJet; i_cut < NCuts; i_cut++){
		if (categoryFlag == "OneJetLow"){
			if (i_cut < OneJetNCuts){
				value.at(i_cut) = value_OneJetLow.at(i_cut);
				pass.at(i_cut) = pass_OneJetLow.at(i_cut);
			}
			else{
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		catPassed = catPassed && pass_OneJetLow.at(i_cut);
	}

	return catPassed;
}

void HToTaumuTauh::configure_ZeroJetHigh(){
	// to be called only if ZeroJetHigh is chosen category

	// set cut values to be the cut values of this category
	cut.at(ZeroJetNJet) = cut_ZeroJet.at(ZeroJetNJet);
	cut.at(ZeroJetNBtagJets) = cut_ZeroJet.at(ZeroJetNBtagJets);
	cut.at(ZeroJetTauPt) = cut_ZeroJet.at(ZeroJetTauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(ZeroJetNJet)="Number Jets $<=$";
	title.at(ZeroJetNJet)+=cut.at(ZeroJetNJet);
	htitle=title.at(ZeroJetNJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=ZeroJetNJet;
	Nminus1.at(ZeroJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(ZeroJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(ZeroJetNBtagJets)="Number btag jets $<=$";
	title.at(ZeroJetNBtagJets)+=cut.at(ZeroJetNBtagJets);
	htitle=title.at(ZeroJetNBtagJets);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of btag jets";
	c="_Cut_";c+=ZeroJetNBtagJets;
	Nminus1.at(ZeroJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(ZeroJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(ZeroJetTauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(ZeroJetTauPt)+=cut.at(ZeroJetTauPt);
	title.at(ZeroJetTauPt)+=" GeV";
	htitle=title.at(ZeroJetTauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=ZeroJetTauPt;
	Nminus1.at(ZeroJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(ZeroJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_ZeroJetHigh(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection){
	std::vector<float> value_ZeroJetHigh;
	std::vector<float> pass_ZeroJetHigh;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_ZeroJetHigh.push_back(-10.);
	pass_ZeroJetHigh.push_back(false);
	}

	value_ZeroJetHigh.at(ZeroJetNJet) = jetCollection.size();
	pass_ZeroJetHigh.at(ZeroJetNJet) = ( value_ZeroJetHigh.at(ZeroJetNJet) <= cut_ZeroJet.at(ZeroJetNJet) );

	value_ZeroJetHigh.at(ZeroJetNBtagJets) = bJetCollection.size();
	pass_ZeroJetHigh.at(ZeroJetNBtagJets) = ( value_ZeroJetHigh.at(ZeroJetNBtagJets) <= cut_ZeroJet.at(ZeroJetNBtagJets) );

	if (selTau == -1){
		value_ZeroJetHigh.at(ZeroJetTauPt) = -10.;
		pass_ZeroJetHigh.at(ZeroJetTauPt) = false;
	}
	else{
		value_ZeroJetHigh.at(ZeroJetTauPt) = Ntp->PFTau_p4(selTau).Pt();
		pass_ZeroJetHigh.at(ZeroJetTauPt) = ( value_ZeroJetHigh.at(ZeroJetTauPt) >= cut_ZeroJet.at(ZeroJetTauPt) );
	}


	// migrate into main analysis if this is chosen category
	bool catPassed = true;
	for (unsigned i_cut = ZeroJetNJet; i_cut < NCuts; i_cut++){
		if (categoryFlag == "ZeroJetHigh"){
			if (i_cut < ZeroJetNCuts){
				value.at(i_cut) = value_ZeroJetHigh.at(i_cut);
				pass.at(i_cut) = pass_ZeroJetHigh.at(i_cut);
			}
			else{
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		catPassed = catPassed && pass_ZeroJetHigh.at(i_cut);
	}

	return catPassed;
}

void HToTaumuTauh::configure_ZeroJetLow(){
	// to be called only if ZeroJetLow is chosen category

	// set cut values to be the cut values of this category
	cut.at(ZeroJetNJet) = cut_ZeroJet.at(ZeroJetNJet);
	cut.at(ZeroJetNBtagJets) = cut_ZeroJet.at(ZeroJetNBtagJets);
	cut.at(ZeroJetTauPt) = cut_ZeroJet.at(ZeroJetTauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(ZeroJetNJet)="Number Jets $<=$";
	title.at(ZeroJetNJet)+=cut.at(ZeroJetNJet);
	htitle=title.at(ZeroJetNJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=ZeroJetNJet;
	Nminus1.at(ZeroJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(ZeroJetNJet) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetNJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(ZeroJetNBtagJets)="Number btag jets $<=$";
	title.at(ZeroJetNBtagJets)+=cut.at(ZeroJetNBtagJets);
	htitle=title.at(ZeroJetNBtagJets);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of btag jets";
	c="_Cut_";c+=ZeroJetNBtagJets;
	Nminus1.at(ZeroJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(ZeroJetNBtagJets) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetNBtagJets_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(ZeroJetTauPt)="$p_{T}(\\tau_{h}) <$";
	title.at(ZeroJetTauPt)+=cut.at(ZeroJetTauPt);
	title.at(ZeroJetTauPt)+=" GeV";
	htitle=title.at(ZeroJetTauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=ZeroJetTauPt;
	Nminus1.at(ZeroJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(ZeroJetTauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetTauPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_ZeroJetLow(int selTau, std::vector<int> jetCollection, std::vector<int> bJetCollection){
	std::vector<float> value_ZeroJetLow;
	std::vector<float> pass_ZeroJetLow;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_ZeroJetLow.push_back(-10.);
	pass_ZeroJetLow.push_back(false);
	}

	value_ZeroJetLow.at(ZeroJetNJet) = jetCollection.size();
	pass_ZeroJetLow.at(ZeroJetNJet) = ( value_ZeroJetLow.at(ZeroJetNJet) <= cut_ZeroJet.at(ZeroJetNJet) );

	value_ZeroJetLow.at(ZeroJetNBtagJets) = bJetCollection.size();
	pass_ZeroJetLow.at(ZeroJetNBtagJets) = ( value_ZeroJetLow.at(ZeroJetNBtagJets) <= cut_ZeroJet.at(ZeroJetNBtagJets) );

	if (selTau == -1){
		value_ZeroJetLow.at(ZeroJetTauPt) = -10.;
		pass_ZeroJetLow.at(ZeroJetTauPt) = false;
	}
	else{
		value_ZeroJetLow.at(ZeroJetTauPt) = Ntp->PFTau_p4(selTau).Pt();
		pass_ZeroJetLow.at(ZeroJetTauPt) = ( value_ZeroJetLow.at(ZeroJetTauPt) < cut_ZeroJet.at(ZeroJetTauPt) );
	}


	// migrate into main analysis if this is chosen category
	bool catPassed = true;
	for (unsigned i_cut = ZeroJetNJet; i_cut < NCuts; i_cut++){
		if (categoryFlag == "ZeroJetLow"){
			if (i_cut < ZeroJetNCuts){
				value.at(i_cut) = value_ZeroJetLow.at(i_cut);
				pass.at(i_cut) = pass_ZeroJetLow.at(i_cut);
			}
			else{
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		catPassed = catPassed && pass_ZeroJetLow.at(i_cut);
	}
	return catPassed;
}

void HToTaumuTauh::configure_NoCategory(){
	// to be called only if No Category shall be run

	// set cut values to be the cut values of this category
	// nothing to do

	// set histograms of category cuts
	// no histograms to be set
}
bool HToTaumuTauh::category_NoCategory(){
	std::vector<float> value_NoCategory;
	std::vector<float> pass_NoCategory;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_NoCategory.push_back(-10.);
	pass_NoCategory.push_back(false);
	}

	// not cuts to compute

	// migrate into main analysis if this is chosen category
	bool catPassed = true;
	for (unsigned i_cut = NoCategoryNCuts; i_cut < NCuts; i_cut++){
		if (categoryFlag == "NoCategory"){
			if(i_cut < NoCategoryNCuts){
				value.at(i_cut) = value_NoCategory.at(i_cut);
				pass.at(i_cut) = pass_NoCategory.at(i_cut);
			}
			else{
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		catPassed = catPassed && pass_NoCategory.at(i_cut);
	}

	return catPassed;
}

