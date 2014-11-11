#include "ZToTaumuTauh.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SkimConfig.h"

ZToTaumuTauh::ZToTaumuTauh(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_dxy(0.045),// 100 is dummy value
  cMu_dz(0.2),
  cMu_relIso(0.1),
  cMu_pt(20.),
  cMu_eta(2.1),
  cMu_dRHltMatch(0.5),
  cTau_pt(20.),//Recommended by Tau POG
  cTau_eta(2.3),
  cMuTau_dR(0.3),
  cTau_IsoRaw(1.5)
{
	TString trigNames[] = {"HLT_IsoMu18_eta2p1_LooseIsoPFTau20","HLT_IsoMu17_eta2p1_LooseIsoPFTau20"};
	std::vector<TString> temp (trigNames, trigNames + sizeof(trigNames) / sizeof(TString) );
	cTriggerNames = temp;

	//WJets BG Method
	SB_lowerLimit = 70; //in GeV, Region > SB_lowerLimit dominated by W+Jets
	SB_upperLimit = 140; //in GeV, Region < SB_upperLimit
	Scaleby_Counting = true; // == false --> Scale by Integral

	//Set verbose boolean
	verbose = false;
}

ZToTaumuTauh::~ZToTaumuTauh(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZToTaumuTauh::~ZToTaumuTauh Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZToTaumuTauh::~ZToTaumuTauh()" << std::endl;
}

void  ZToTaumuTauh::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    cut.at(TriggerOk)=0;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=1;
    if(i==NMuId)		cut.at(NMuId)=1;
    if(i==NMuKin)		cut.at(NMuKin)=1;
    if(i==NMuIso)		cut.at(NMuIso)=1;
    if(i==NTauId)		cut.at(NTauId)=1;
    if(i==NTauKin)		cut.at(NTauKin)=1;
    if(i==NTauIso)		cut.at(NTauIso)=1;
    if(i==ChargeSum)	cut.at(ChargeSum)=0;
    if(i==MT_MuMET)		cut.at(MT_MuMET)=44;
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,41,-0.5,40.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,41,-0.5,40.5,hlabel,"Events"));
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
    	title.at(i_cut)="Number $\\mu_{kin} >=$";
    	title.at(i_cut)+=cut.at(NMuKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{Kin}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i_cut==NMuIso){
    	title.at(i_cut)="Number $\\mu_{Iso} >=$";
    	title.at(i_cut)+=cut.at(NMuIso);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i_cut==NTauId){
    	title.at(i_cut)="Number $\\tau_{ID} >=$";
    	title.at(i_cut)+=cut.at(NTauId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NTauKin){
    	title.at(i_cut)="Number $\\tau_{Kin} >=$";
    	title.at(i_cut)+=cut.at(NTauKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Kin}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NTauIso){
    	title.at(i_cut)="Number $\\tau_{Iso} >=$";
    	title.at(i_cut)+=cut.at(NTauIso);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauIso_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauIso_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==ChargeSum){
    	title.at(i_cut)="Sum of Charges = ";
    	title.at(i_cut)+=cut.at(ChargeSum);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Sum of Charges";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ChargeSum_",htitle,5,-2.5,2.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ChargeSum_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }
    else if(i_cut==MT_MuMET){
    	title.at(i_cut)="$m_{T}(\\mu,MET) <$";
    	title.at(i_cut)+=cut.at(MT_MuMET);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="m_{T}(#mu,MET)";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_MuMET_",htitle,75,0,150.,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_MuMET_",htitle,75,0,150.,hlabel,"Events"));
    }
  }
  // Setup NPassed Histograms
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",41,-0.5,40.5,"Number of Vertices","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of good Vertices","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  NMtTauMET=HConfig.GetTH1D(Name+"_NMtTauMET","NMtTauMET",75,-0.5,149.5,"m_{T}(#tau,MET)","Events");
  NMvis=HConfig.GetTH1D(Name+"_NMvis","NMvis",75,-0.5,149.5,"m_{vis}(#mu,#tau)","Events");
  NSB=HConfig.GetTH1D(Name+"_NSB","NSB",2,-0.5,1.5,"SB Region (all Samples)","Events");

  Mu_pt=HConfig.GetTH1D(Name+"_Mu_pt","Mu_pt",100,0,200,"Muon p_{t}","Events");
  Mu_phi=HConfig.GetTH1D(Name+"_Mu_phi","Mu_phi",10,0,3.14159265359,"Muon #phi","Events");
  Mu_eta=HConfig.GetTH1D(Name+"_Mu_phi","Mu_eta",50,-2.5,2.5,"Muon #eta","Events");

  Tau_pt=HConfig.GetTH1D(Name+"_Tau_pt","Tau_pt",100,0,200,"Tau p_{t}","Events");
  Tau_phi=HConfig.GetTH1D(Name+"_Tau_phi","Tau_phi",10,0,3.14159265359,"Tau #phi","Events");
  Tau_eta=HConfig.GetTH1D(Name+"_Tau_phi","Tau_eta",50,-2.5,2.5,"Tau #eta","Events");

  NQCD_MT_MuMET_A=HConfig.GetTH1D(Name+"_NQCD_MT_MuMET_A","NQCD_MT_MuMET_A",75,0,150.,"m_{T}(#mu,MET) in A","Events");
  NQCD_MT_MuMET_B=HConfig.GetTH1D(Name+"_NQCD_MT_MuMET_B","NQCD_MT_MuMET_B",75,0,150.,"m_{T}(#mu,MET) in B","Events");
  NQCD_MT_MuMET_C=HConfig.GetTH1D(Name+"_NQCD_MT_MuMET_C","NQCD_MT_MuMET_C",75,0,150.,"m_{T}(#mu,MET) in C","Events");
  NQCD_MT_MuMET_D=HConfig.GetTH1D(Name+"_NQCD_MT_MuMET_D","NQCD_MT_MuMET_D",75,0,150.,"m_{T}(#mu,MET) in D","Events");

  NQCD_MET_A=HConfig.GetTH1D(Name+"_NQCD_MET_A","NQCD_MET_A",75,0,150.,"MVA MET in A","Events");
  NQCD_MET_B=HConfig.GetTH1D(Name+"_NQCD_MET_B","NQCD_MET_B",75,0,150.,"MVA MET in B","Events");
  NQCD_MET_C=HConfig.GetTH1D(Name+"_NQCD_MET_C","NQCD_MET_C",75,0,150.,"MVA MET in C","Events");
  NQCD_MET_D=HConfig.GetTH1D(Name+"_NQCD_MET_D","NQCD_MET_D",75,0,150.,"MVA MET in D","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  ZToTaumuTauh::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&NMtTauMET);
 Extradist1d.push_back(&NMvis);
 Extradist1d.push_back(&NSB);
 Extradist1d.push_back(&Mu_pt);
 Extradist1d.push_back(&Mu_phi);
 Extradist1d.push_back(&Mu_eta);
 Extradist1d.push_back(&Tau_pt);
 Extradist1d.push_back(&Tau_phi);
 Extradist1d.push_back(&Tau_eta);
 Extradist1d.push_back(&NQCD_MT_MuMET_A);
 Extradist1d.push_back(&NQCD_MT_MuMET_B);
 Extradist1d.push_back(&NQCD_MT_MuMET_C);
 Extradist1d.push_back(&NQCD_MT_MuMET_D);
 Extradist1d.push_back(&NQCD_MET_A);
 Extradist1d.push_back(&NQCD_MET_B);
 Extradist1d.push_back(&NQCD_MET_C);
 Extradist1d.push_back(&NQCD_MET_D);
}

void  ZToTaumuTauh::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  int selVertex = -1;
  int selMuon_Iso = -1;
  int selMuon_AntiIso = -1;
  int selTau = -1;

  TString tau_corr = "";
  if(Ntp->GetMCID() == DataMCType::DY_tautau || (Ntp->GetMCID()>=10 && Ntp->GetMCID()<= 13)) tau_corr = "scalecorr";

  // Apply Selection
  if(verbose) std::cout << "Cut on good vertex" << std::endl;
  unsigned int nGoodVtx=0;
  for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
    if(Ntp->isVtxGood(i_vtx)){
    	if(selVertex == -1) selVertex = i_vtx; // selected vertex = first vertex (highest sum[pT^2]) to fulfill vertex requirements
    	nGoodVtx++;
    }
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  // Trigger
  if(verbose) std::cout << "Cut on Trigger" << std::endl;
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
  if(verbose) std::cout << "Cut on MuonID" << std::endl;
  std::vector<int> selectedMuonsId;
  selectedMuonsId.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_Id(i_mu,selVertex) ) {
		  selectedMuonsId.push_back(i_mu);
	  }
  }
  value.at(NMuId)=selectedMuonsId.size();
  pass.at(NMuId)=(value.at(NMuId)>=cut.at(NMuId));

  if(verbose) std::cout << "Cut on Muon Kinematics" << std::endl;
  std::vector<int> selectedMuonsKin;
  selectedMuonsKin.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsId.begin(); it_mu != selectedMuonsId.end(); ++it_mu){
	  if( selectMuon_Kinematics(*it_mu) ) {
		  selectedMuonsKin.push_back(*it_mu);
	  }
  }
  value.at(NMuKin)=selectedMuonsKin.size();
  pass.at(NMuKin)=(value.at(NMuKin)>=cut.at(NMuKin));

  if(verbose) std::cout << "Cut on Muon Isolation (Iso)" << std::endl;
  std::vector<int> selectedMuonsIso;
  selectedMuonsIso.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsKin.begin(); it_mu != selectedMuonsKin.end(); ++it_mu){
	  if( selectMuon_Isolation(*it_mu) ) {
		  if(selMuon_Iso == -1) selMuon_Iso = *it_mu;
		  selectedMuonsIso.push_back(*it_mu);
	  }
  }
  value.at(NMuIso)=selectedMuonsIso.size();
  pass.at(NMuIso)=(value.at(NMuIso)>=cut.at(NMuIso));

  if(verbose) std::cout << "Cut on Muon Isolation (Anti Iso)" << std::endl;
  std::vector<int> selectedMuonsAntiIso;
  selectedMuonsAntiIso.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsKin.begin(); it_mu != selectedMuonsKin.end(); ++it_mu){
	  if( selectMuon_AntiIsolation(*it_mu) ) {
		  if(selMuon_AntiIso == -1 && selMuon_Iso == -1) selMuon_AntiIso = *it_mu;
		  selectedMuonsAntiIso.push_back(*it_mu);
	  }
  }

  // Tau cuts
  if(verbose) std::cout << "Cut on TauID" << std::endl;
  std::vector<int> selectedTausId;
  selectedTausId.clear();
  for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
	  if ( selectPFTau_Id(i_tau,selectedMuonsId) ){
		  selectedTausId.push_back(i_tau);
	  }
  }
  value.at(NTauId)=selectedTausId.size();
  pass.at(NTauId)=(value.at(NTauId)>=cut.at(NTauId));

  if(verbose) std::cout << "Cut on Tau Kinematics" << std::endl;
  std::vector<int> selectedTausKin;
  selectedTausKin.clear();
  for(std::vector<int>::iterator it_tau = selectedTausId.begin(); it_tau != selectedTausId.end(); ++it_tau){
	  if ( selectPFTau_Kinematics(*it_tau) ){
		  selectedTausKin.push_back(*it_tau);
	  }
  }
  value.at(NTauKin)=selectedTausKin.size();
  pass.at(NTauKin)=(value.at(NTauKin)>=cut.at(NTauKin));

  if(verbose) std::cout << "Cut on Tau Isolation" << std::endl;
  std::vector<int> selectedTausIso;
  selectedTausIso.clear();
  for(std::vector<int>::iterator it_tau = selectedTausKin.begin(); it_tau != selectedTausKin.end(); ++it_tau){
	  if ( selectPFTau_Isolation(*it_tau) ){
		  if(selTau == -1) selTau = *it_tau;
		  selectedTausIso.push_back(*it_tau);
	  }
  }
  value.at(NTauIso)=selectedTausIso.size();
  pass.at(NTauIso)=(value.at(NTauIso)>=cut.at(NTauIso));

  // Charge of MuTau
  if(verbose) std::cout << "Cut on Charge of MuTau System" << std::endl;
  value.at(ChargeSum) = 666;
  if(selTau != -1){
	  if(selMuon_Iso != -1 && selMuon_AntiIso == -1){
		  Charge = Ntp->Muon_Charge(selMuon_Iso) + Ntp->PFTau_Charge(selTau);
		  value.at(ChargeSum) = Charge;
  	  }
  	  else if(selMuon_Iso == -1 && selMuon_AntiIso != -1){
  		  Charge = Ntp->Muon_Charge(selMuon_AntiIso) + Ntp->PFTau_Charge(selTau);
  		  value.at(ChargeSum) = Charge;
  	  }
  	  else{
  		  Charge = 666;
  	  }
  }
  else{
	  Charge = 666;
  }
  pass.at(ChargeSum)=(value.at(ChargeSum)==cut.at(ChargeSum));


  // MT calculation
  if(verbose) std::cout << "Calculation and Cut on MT distribution" << std::endl;
  double pT,phi,eTmiss,eTmPhi;
  double MT_TauMET, MT_MuMET_AntiIso;

  if(selMuon_Iso == -1 && selMuon_AntiIso == -1){
	  value.at(MT_MuMET) = -10;
	  if(verbose) std::cout << "No Muon selected: neither isolated or anti isolated" << std::endl;
  }
  else if(selMuon_Iso != -1 && selMuon_AntiIso != -1){
	  value.at(MT_MuMET) = -10;
	  std::cout << "CRITICAL: SELECTED MUON PASSED ISOLATION AND ANTI-ISOLATION CUT --> FIX" << std::endl;
  }
  else if(selMuon_Iso != -1 && selMuon_AntiIso == -1){
	  //std::cout << "selMuon_Iso is " << selMuon_Iso << std::endl;
	  eTmiss				 	= Ntp->MET_CorrMVAMuTau_et();
	  eTmPhi 					= Ntp->MET_CorrMVAMuTau_phi();
	  pT 						= Ntp->Muon_p4(selMuon_Iso).Pt();
	  phi						= Ntp->Muon_p4(selMuon_Iso).Phi();
	  value.at(MT_MuMET)		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	  //std::cout << "MT_MuMET is " << MT_MuMET << std::endl;
  }
  else if(selMuon_AntiIso != -1 && selMuon_Iso == -1){
	  eTmiss				 	= Ntp->MET_CorrMVAMuTau_et();
	  eTmPhi 					= Ntp->MET_CorrMVAMuTau_phi();
	  pT 						= Ntp->Muon_p4(selMuon_AntiIso).Pt();
	  phi						= Ntp->Muon_p4(selMuon_AntiIso).Phi();
	  MT_MuMET_AntiIso			= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	  value.at(MT_MuMET) 		= MT_MuMET_AntiIso;
  }
  if(value.at(MT_MuMET) != -10) pass.at(MT_MuMET)=(value.at(MT_MuMET)<cut.at(MT_MuMET));

  if(selTau == -1){
	  //std::cout << "selTau is " << selTau << std::endl;
	  MT_TauMET = -10;
  }
  else{
	  //std::cout << "selTau is " << selTau << std::endl;
	  eTmiss 		= Ntp->MET_CorrMVAMuTau_et();
	  eTmPhi 		= Ntp->MET_CorrMVAMuTau_phi();
	  pT 			= Ntp->PFTau_p4(selTau, tau_corr).Pt();
	  phi			= Ntp->PFTau_p4(selTau, tau_corr).Phi();
	  MT_TauMET		= Ntp->transverseMass(pT,phi,eTmiss,eTmPhi);
	  //std::cout << "MT_TauMET is " << MT_TauMET << std::endl;
  }

  // Mvis
  if(verbose) std::cout << "Calculation of Mvis" << std::endl;
  double Mvis;

  if(selTau != -1){
	  if(selMuon_Iso != -1 && selMuon_AntiIso == -1){
		  Mvis = (Ntp->PFTau_p4(selTau, tau_corr) + Ntp->Muon_p4(selMuon_Iso)).M();
	  }
	  else if(selMuon_Iso == -1 && selMuon_AntiIso != -1){
		  Mvis = (Ntp->PFTau_p4(selTau, tau_corr) + Ntp->Muon_p4(selMuon_AntiIso)).M();
	  }
	  else{
		  Mvis = -10;
	  }
  }
  else{
	  Mvis = -10;
  }

  // Weights

  double wobs=1;
  double w;
  if(!Ntp->isData()){
	  w = Ntp->PUWeight();
	  if(selMuon_Iso != -1){
		  w *= RSF->HiggsTauTau_MuTau_Id_Mu(Ntp->Muon_p4(selMuon_Iso));
		  w *= RSF->HiggsTauTau_MuTau_Iso_Mu(Ntp->Muon_p4(selMuon_Iso));
		  w *= RSF->HiggsTauTau_MuTau_Trigger_Mu(Ntp->Muon_p4(selMuon_Iso));
	  }
	  if(selTau != -1){
		  w *= RSF->HiggsTauTau_MuTau_Trigger_Tau(Ntp->PFTau_p4(selTau, tau_corr));
	  }
  }
  else{w=1;}

  // W+Jets BG Method
  if(verbose) std::cout << "W+Jets BG Method" << std::endl;
  std::vector<int> ChargeSumMT_MuMET;
  ChargeSumMT_MuMET.clear();
  ChargeSumMT_MuMET.push_back(ChargeSum);
  ChargeSumMT_MuMET.push_back(MT_MuMET);
  if(passAllBut(ChargeSumMT_MuMET)){
  	  if(pass.at(ChargeSum)){ //Opposite Sign WJets yield (bin #1 with value 0)
  		  if(value.at(MT_MuMET) > SB_lowerLimit && value.at(MT_MuMET) < SB_upperLimit){
  			  NSB.at(t).Fill(0.,w);
  		  }
  	  }
  	  if(!pass.at(ChargeSum) && Charge != 666){ //Same Sign WJets yield (bin #2 with value 1)
  		  if(value.at(MT_MuMET) > SB_lowerLimit && value.at(MT_MuMET) < SB_upperLimit){
  			  NSB.at(t).Fill(1.,w);
  		  }
  	  }
  }

  // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  SS
   *        |        *       S
   * ----------------*------ i
   *        |        *       g
   *    A   |    B   *  OS   n
   *        |        *
   *******************
   *  Iso   | AntiIso
   *
   *     relIso(mu)
   */
  if(verbose) std::cout << "QCD ABCD BG Method" << std::endl;
  std::vector<int> MuIsoChargeSumMT_MuMET;
  MuIsoChargeSumMT_MuMET.clear();
  MuIsoChargeSumMT_MuMET.push_back(NMuIso);
  MuIsoChargeSumMT_MuMET.push_back(ChargeSum);
  MuIsoChargeSumMT_MuMET.push_back(MT_MuMET);
  if(passAllBut(MuIsoChargeSumMT_MuMET)){
	  if(pass.at(NMuIso) && selMuon_Iso != -1){
		  if(pass.at(ChargeSum)){ //A --> Signal-Selection (pass all cuts)
			  NQCD_MT_MuMET_A.at(t).Fill(value.at(MT_MuMET),w);
	 	 	  NQCD_MET_A.at(t).Fill(Ntp->MET_CorrMVAMuTau_et());
	 	  }
	 	  if(!pass.at(ChargeSum) && value.at(ChargeSum) != 666){ //C
	 	 	  NQCD_MT_MuMET_C.at(t).Fill(value.at(MT_MuMET),w);
	 	      NQCD_MET_C.at(t).Fill(Ntp->MET_CorrMVAMuTau_et());
	 	  }
	  }
	  if(!pass.at(NMuIso) && selMuon_AntiIso != -1){
		  if(pass.at(ChargeSum)){ //B
	 		  NQCD_MT_MuMET_B.at(t).Fill(MT_MuMET_AntiIso,w);
	 		  NQCD_MET_B.at(t).Fill(Ntp->MET_CorrMVAMuTau_et());
	 	  }
	 	  if(!pass.at(ChargeSum) && Charge != 666){ //D
	 	 	  NQCD_MT_MuMET_D.at(t).Fill(MT_MuMET_AntiIso,w);
	 	 	  NQCD_MET_D.at(t).Fill(Ntp->MET_CorrMVAMuTau_et());
	 	 	  if(id == DataMCType::Data){
		 	 	NQCD_MT_MuMET_D.at(HConfig.GetType(DataMCType::QCD)).Fill(MT_MuMET_AntiIso,w);
		 	 	NQCD_MET_D.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->MET_CorrMVAMuTau_et());
		 	 	NQCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::QCD)).Fill(MT_MuMET_AntiIso,w);
		 	 	NQCD_MET_A.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->MET_CorrMVAMuTau_et());
	 	 		pass.at(ChargeSum) = true;
	 	 		pass.at(NMuIso) = true;
	 	 		bool QCD_status = AnalysisCuts(HConfig.GetType(DataMCType::QCD),w,wobs);
	 	 		if(QCD_status){
	 	 			NVtx.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->NVtx(),w);
	 	 			unsigned int nGoodVtx=0;
	 	 			for(unsigned int i=0;i<Ntp->NVtx();i++){
	 	 				NTrackperVtx.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Vtx_Track_idx(i).size(),w);
	 	 				if(Ntp->isVtxGood(i))nGoodVtx++;
	 	 			}
	 	 			NGoodVtx.at(HConfig.GetType(DataMCType::QCD)).Fill(nGoodVtx,w);
	 	 			NMtTauMET.at(HConfig.GetType(DataMCType::QCD)).Fill(MT_TauMET,w);
	 	 			NMvis.at(HConfig.GetType(DataMCType::QCD)).Fill(Mvis,w);
	 	 		    Mu_pt.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Muon_p4(selMuon_Iso).Pt(),w);
	 	 		    Mu_phi.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Muon_p4(selMuon_Iso).Phi(),w);
	 	 		    Mu_eta.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->Muon_p4(selMuon_Iso).Eta(),w);
	 	 		    Tau_pt.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).Pt(),w);
	 	 		    Tau_phi.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).Phi(),w);
	 	 		    Tau_eta.at(HConfig.GetType(DataMCType::QCD)).Fill(Ntp->PFTau_p4(selTau, tau_corr).Eta(),w);
	 	 		}
	 	 		pass.at(ChargeSum) = false;
	 	 		pass.at(NMuIso) = false;
	 	 	  }
	 	  }
	  }
  }

  bool status=AnalysisCuts(t,w,wobs);
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);
    NMtTauMET.at(t).Fill(MT_TauMET,w);
    NMvis.at(t).Fill(Mvis,w);
    Mu_pt.at(t).Fill(Ntp->Muon_p4(selMuon_Iso).Pt(),w);
    Mu_phi.at(t).Fill(Ntp->Muon_p4(selMuon_Iso).Phi(),w);
    Mu_eta.at(t).Fill(Ntp->Muon_p4(selMuon_Iso).Eta(),w);
    Tau_pt.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).Pt(),w);
    Tau_phi.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).Phi(),w);
    Tau_eta.at(t).Fill(Ntp->PFTau_p4(selTau, tau_corr).Eta(),w);
  }

}//final bracket of DoEvent


void  ZToTaumuTauh::Finish(){
  std::cout << "Starting Finish" << std::endl;
  if(mode==RECONSTRUCT){
	  std::cout << "Enter mode==RECONSTRUCT" << std::endl;
	  SkimConfig SC;
	  SC.ApplySkimEfficiency(types,Npassed,Npassed_noweight);
		for(int i=0; i<Npassed.size();i++){
			nevents_noweight_default.push_back(Npassed_noweight.at(i).GetBinContent(1));
		}
  	  std::vector<double> SB_Integral;
  	  double SB_Integral_Data_minus_MC = 0;
  	  double SB_Integral_WJets = 0;

  	  std::vector<double> SB_Counting_OS;
  	  double SB_Counting_Data_minus_MC_OS = 0;
  	  double SB_Counting_WJets_OS = 0;

  	  std::vector<double> SB_Counting_SS;
  	  double SB_Counting_Data_minus_MC_SS = 0;
  	  double SB_Counting_WJets_SS = 0;

  	  std::vector<double> QCD_Integral_B;
  	  double QCD_Integral_B_Data_minus_MC = 0;
  	  std::vector<double> QCD_Integral_C;
  	  double QCD_Integral_C_Data_minus_MC = 0;
  	  std::vector<double> QCD_Integral_D;
  	  double QCD_Integral_D_Data_minus_MC = 0;

  	  std::cout << "W+Jets Background Method " << std::endl;
  	  for(int i=0;i<CrossSectionandAcceptance.size();i++){
  	  	  int Bin_SB_lowerLimit = Nminus1.at(MT_MuMET).at(i).FindFixBin(SB_lowerLimit);
  	  	  int Bin_SB_upperLimit = Nminus1.at(MT_MuMET).at(i).FindFixBin(SB_upperLimit) -1;
  		  SB_Integral.push_back(Nminus1.at(MT_MuMET).at(i).Integral(Bin_SB_lowerLimit, Bin_SB_upperLimit));
  		  SB_Counting_OS.push_back(NSB.at(i).GetBinContent(1));
  		  SB_Counting_SS.push_back(NSB.at(i).GetBinContent(2));
  		  if(CrossSectionandAcceptance.at(i)>0){
			  SB_Integral.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
			  SB_Counting_OS.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
			  SB_Counting_SS.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
		  }
		  std::cout << "SB integral 70-140GeV from Nminus1 plot for Sample " << i << " : " << SB_Integral.at(i) << std::endl;
		  std::cout << "SB integral 70-140GeV from Counting for OS for Sample " << i << " : " << SB_Counting_OS.at(i) << std::endl;
		  std::cout << "SB integral 70-140GeV from Counting for SS for Sample " << i << " : " << SB_Counting_SS.at(i) << std::endl;
  	  }
  	  for(int i=0;i<CrossSectionandAcceptance.size();i++){
  		  if(HConfig.GetID(i) == DataMCType::W_lnu || HConfig.GetID(i) == DataMCType::W_taunu){
  			  SB_Integral_WJets += SB_Integral.at(i);
  			  SB_Counting_WJets_OS += SB_Counting_OS.at(i);
  			  SB_Counting_WJets_SS += SB_Counting_SS.at(i);
  		  }
  		  else if(HConfig.GetID(i) == DataMCType::Data){
  			  SB_Integral_Data_minus_MC += SB_Integral.at(i);
  			  SB_Counting_Data_minus_MC_OS += SB_Counting_OS.at(i);
  			  SB_Counting_Data_minus_MC_SS += SB_Counting_SS.at(i);
  		  }
  		  else if(HConfig.GetID(i) != DataMCType::Data){
  			  SB_Integral_Data_minus_MC -= SB_Integral.at(i);
  			  SB_Counting_Data_minus_MC_OS -= SB_Counting_OS.at(i);
  			  SB_Counting_Data_minus_MC_SS -= SB_Counting_SS.at(i);
  		  }
  	  }
  	  if(SB_Integral_Data_minus_MC > 0 && SB_Integral_WJets > 0){
  		  std::cout << "Scale Factor for W+Jets Sample with Histo.Integral Method: " << SB_Integral_Data_minus_MC/SB_Integral_WJets << std::endl;
  		  if(!Scaleby_Counting){
  	  		  std::cout << "Scaling by Histo.Integral Method" << std::endl;
  			  ScaleAllHistOfType(HConfig.GetType(DataMCType::W_lnu), SB_Integral_Data_minus_MC/SB_Integral_WJets);
  			  ScaleAllHistOfType(HConfig.GetType(DataMCType::W_taunu), SB_Integral_Data_minus_MC/SB_Integral_WJets);
  		  }
  	  }
  	  else{
		  std::cout << "SB_Integral_WJets is: " << SB_Integral_WJets << std::endl;
		  std::cout << "SB_Integral_Data_minus_MC is: " << SB_Integral_Data_minus_MC << std::endl;
  	  }
  	  if(SB_Counting_Data_minus_MC_OS > 0 && SB_Counting_WJets_OS > 0){
		  std::cout << "Scale Factor for W+Jets Sample with Counting Method for OS: " << SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS << std::endl;
		  std::cout << "Scaleby_Counting boolean is set to " << Scaleby_Counting << std::endl;
  		  if(Scaleby_Counting){
  	  		  std::cout << "Scaling by Counting Method" << std::endl;
  			  //ScaleAllHistOfType(HConfig.GetType(DataMCType::W_lnu), SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  			  //ScaleAllHistOfType(HConfig.GetType(DataMCType::W_taunu), SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MT_MuMET_B.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MT_MuMET_B.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MET_A.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MET_A.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MET_B.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  	  		  NQCD_MET_B.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  	  for(int i=0; i<Nminus1.size(); i++){
  		  		  Nminus0.at(i).at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Nminus0.at(i).at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Nminus1.at(i).at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  		  Nminus1.at(i).at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  	  }
  		  	  NVtx.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  	  NVtx.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  NMtTauMET.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  NMtTauMET.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  NMvis.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  NMvis.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	    	  Mu_pt.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	    	  Mu_pt.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Mu_phi.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Mu_phi.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Mu_eta.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Mu_eta.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Tau_pt.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Tau_pt.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Tau_phi.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 	  	  Tau_phi.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Tau_eta.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
	 		  Tau_eta.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_OS/SB_Counting_WJets_OS);
  		  }
  	  }
  	  else{
		  std::cout << "SB_Counting_WJets_OS is: " << SB_Counting_WJets_OS << std::endl;
		  std::cout << "SB_Counting_Data_minus_MC_OS is: " << SB_Counting_Data_minus_MC_OS << std::endl;
  	  }
  	  if(SB_Counting_Data_minus_MC_SS > 0 && SB_Counting_WJets_SS > 0){
		  std::cout << "Scale Factor for W+Jets Sample with Counting Method for SS: " << SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS << std::endl;
		  std::cout << "Scaleby_Counting boolean is set to " << Scaleby_Counting << std::endl;
		  NQCD_MT_MuMET_C.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
		  NQCD_MT_MuMET_C.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  NQCD_MT_MuMET_D.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  NQCD_MT_MuMET_D.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  NQCD_MET_C.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  NQCD_MET_C.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  NQCD_MET_D.at(HConfig.GetType(DataMCType::W_lnu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
	  	  NQCD_MET_D.at(HConfig.GetType(DataMCType::W_taunu)).Scale(SB_Counting_Data_minus_MC_SS/SB_Counting_WJets_SS);
  	  }
  	  else{
		  std::cout << "SB_Counting_WJets_SS is: " << SB_Counting_WJets_SS << std::endl;
		  std::cout << "SB_Counting_Data_minus_MC_SS is: " << SB_Counting_Data_minus_MC_SS << std::endl;
  	  }
  	  for(int i=0;i<CrossSectionandAcceptance.size();i++){
  		  QCD_Integral_B.push_back(NQCD_MT_MuMET_B.at(i).Integral());
  		  QCD_Integral_C.push_back(NQCD_MT_MuMET_C.at(i).Integral());
  		  QCD_Integral_D.push_back(NQCD_MT_MuMET_D.at(i).Integral());
  		  if(CrossSectionandAcceptance.at(i)>0){
	  		  QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	  		  QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	  		  QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
		  }
  	  }
  	  for(int i=0;i<CrossSectionandAcceptance.size();i++){
  		  if(HConfig.GetID(i) == DataMCType::Data){
  			  QCD_Integral_B_Data_minus_MC += QCD_Integral_B.at(i);
  			  QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
  			  QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
  		  }
  		  else if(HConfig.GetID(i) != DataMCType::QCD){
  			  QCD_Integral_B_Data_minus_MC -= QCD_Integral_B.at(i);
  			  QCD_Integral_C_Data_minus_MC -= QCD_Integral_C.at(i);
  			  QCD_Integral_D_Data_minus_MC -= QCD_Integral_D.at(i);
  		  }
  	  }
  	  if(QCD_Integral_B_Data_minus_MC > 0 && QCD_Integral_C_Data_minus_MC > 0 && QCD_Integral_D_Data_minus_MC > 0){
		  std::cout << "Factor AntiIso OS/SS QCD Sample: " << QCD_Integral_B_Data_minus_MC/QCD_Integral_D_Data_minus_MC << std::endl;
		  std::cout << "Scale Factor for QCD Sample: " << QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC << std::endl;
  		  //ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD), QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
		  NQCD_MT_MuMET_A.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
		  NQCD_MET_A.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
		  for(int i=0; i<Nminus1.size(); i++){
			  Nminus0.at(i).at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
		  	  Nminus1.at(i).at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
		  }
		  NVtx.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  NMtTauMET.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  NMvis.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
    	  Mu_pt.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  Mu_phi.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  Mu_eta.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  Tau_pt.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  Tau_phi.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
 		  Tau_eta.at(HConfig.GetType(DataMCType::QCD)).Scale(QCD_Integral_B_Data_minus_MC*QCD_Integral_C_Data_minus_MC/QCD_Integral_D_Data_minus_MC/QCD_Integral_D_Data_minus_MC);
  	  }
  	  else{
		  std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
		  std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
		  std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
  	  }
  }
  // weight all Histograms
  std::cout << "Starting Selection::Finish" << std::endl;
  Selection::Finish();
}


/////////////////////////////////////////
// Definition of selection and helper functions
/////////////////////////////////////////

///////// Muons

bool ZToTaumuTauh::selectMuon_Id(unsigned i, unsigned vertex){
	if(	Ntp->isSelectedMuon(i,vertex,cMu_dxy,cMu_dz) &&
		Ntp->matchTrigger(Ntp->Muon_p4(i),cTriggerNames,"muon") < cMu_dRHltMatch
			){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectMuon_Kinematics(unsigned i){
	if(	Ntp->Muon_p4(i).Pt() >= cMu_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) <= cMu_eta
			){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectMuon_Isolation(unsigned i){
	if(Ntp->Muon_RelIso(i) < cMu_relIso){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectMuon_AntiIsolation(unsigned i){
	if(		Ntp->Muon_RelIso(i) < 0.5 &&
			Ntp->Muon_RelIso(i) > 0.2){
		return true;
	}
	return false;
}

///////// Taus
bool ZToTaumuTauh::selectPFTau_Id(unsigned i){
	if ( 	Ntp->PFTau_isHPSByDecayModeFinding(i) &&
			Ntp->PFTau_isHPSAgainstElectronsLoose(i) &&
			Ntp->PFTau_isHPSAgainstMuonTight(i)
			){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectPFTau_Id(unsigned i, std::vector<int> muonCollection){
	// check if tau is matched to a muon, if so this is not a good tau
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Sync_Issues
	for(std::vector<int>::iterator it_mu = muonCollection.begin(); it_mu != muonCollection.end(); ++it_mu){
	  if( Ntp->PFTau_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR ) {
		  return false;
	  }
	}
	// trigger matching
	if (Ntp->matchTrigger(Ntp->PFTau_p4(i),cTriggerNames,"tau") > cMu_dRHltMatch) {
		return false;
	}

	if ( 	selectPFTau_Id(i) ){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectPFTau_Isolation(unsigned i){
	if(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(i) < cTau_IsoRaw){
		return true;
	}
	return false;
}
bool ZToTaumuTauh::selectPFTau_Kinematics(unsigned i){
	if(	Ntp->PFTau_p4(i).Pt() >= cTau_pt &&
		fabs(Ntp->PFTau_p4(i).Eta()) <= cTau_eta
			){
		return true;
	}
	return false;
}
