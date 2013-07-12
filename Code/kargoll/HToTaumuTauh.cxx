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
  cTau_eta(2.3)
{
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
    if(i==TriggerOk)    	cut.at(TriggerOk)=1;
    if(i==PrimeVtx)     	cut.at(PrimeVtx)=1;
    if(i==NMuId)			cut.at(NMuId)=1;
    if(i==NMuKin)			cut.at(NMuKin)=1;
    if(i==DiMuonVeto)		cut.at(DiMuonVeto)=false;
    if(i==NTauId)			cut.at(NTauId)=1;
    if(i==NTauIso)			cut.at(NTauIso)=1;
    if(i==NTauKin)			cut.at(NTauKin)=1;
    if(i==OppCharge)		cut.at(OppCharge)=0;
    if(i==TriLeptonVeto)	cut.at(TriLeptonVeto)=0;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NMuId){
    	title.at(i)="Number $\\mu_{ID} =$";
    	title.at(i)+=cut.at(NMuId);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==NMuKin){
    	title.at(i)="Number $\\mu_{sel} =$";
    	title.at(i)+=cut.at(NMuKin);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{sel}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==DiMuonVeto){
    	title.at(i)="Veto on $\\mu_{veto}$ pair";
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Veto on #mu_{veto} pair";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NTauId){
    	title.at(i)="Number $\\tau_{ID} =$";
    	title.at(i)+=cut.at(NTauId);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauId_",htitle,26,-0.5,25.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauId_",htitle,26,-0.5,25.5,hlabel,"Events"));
    }
    else if(i==NTauIso){
    	title.at(i)="Number $\\tau_{Iso} =$";
    	title.at(i)+=cut.at(NTauIso);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauIso_",htitle,16,-0.5,15.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauIso_",htitle,16,-0.5,15.5,hlabel,"Events"));
    }
    else if(i==NTauKin){
    	title.at(i)="Number $\\tau_{sel} =$";
    	title.at(i)+=cut.at(NTauKin);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{sel}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==OppCharge){
    	title.at(i)="$q(\\mu)+q(\\tau) =$";
    	title.at(i)+=cut.at(OppCharge);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="q(#mu)+q(#tau)";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OppCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OppCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }
    else if(i==TriLeptonVeto){
    	title.at(i)="3 lepton veto: $N(\\mu)+N(e) =$";
    	title.at(i)+=cut.at(TriLeptonVeto);
    	htitle=title.at(i);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of tri-lepton veto leptons";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriLeptonVeto_",htitle,5,-0.5,4.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriLeptonVeto_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  // Setup Extra Histograms
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

  TauPt=HConfig.GetTH1D(Name+"_TauPt","TauPt",50,0.,200.,"p_{T}(#tau)/GeV");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",50,-2.5,2.5,"#eta(#tau)");

  MuVetoDPtSelMuon=HConfig.GetTH1D(Name+"_MuVetoDPtSelMuon","MuVetoDPtSelMuon",100,-100.,100.,"#Deltap_{T}(#mu_{veto},#mu)/GeV");
  MuVetoInvM=HConfig.GetTH1D(Name+"_MuVetoInvM","MuVetoInvM",100,0.,200,"m_{inv}(#mu_{veto}^1,#mu_{veto}^2)/GeV");
  MuVetoPtPositive=HConfig.GetTH1D(Name+"_MuVetoPtPositive","MuVetoPtPositive",50,0.,200.,"p_{T}(#mu_{veto}^{+})/GeV");
  MuVetoPtNegative=HConfig.GetTH1D(Name+"_MuVetoPtNegative","MuVetoPtNegative",50,0.,200.,"p_{T}(#mu_{veto}^{-})/GeV");
  MuVetoDRTau=HConfig.GetTH1D(Name+"_MuVetoDRTau","MuVetoDRTau",50,0.,5.,"#DeltaR(#mu_{veto},#tau_{h})");

  NMuonTriLepVeto=HConfig.GetTH1D(Name+"_NMuonTriLepVeto","NMuonTriLepVeto",5,-0.5,4.5,"N(#mu_{3l veto})");
  NElecTriLepVeto=HConfig.GetTH1D(Name+"_NElecTriLepVeto","NElecTriLepVeto",5,-0.5,4.5,"N(e_{3l veto})");

  MuCharge=HConfig.GetTH1D(Name+"_MuCharge","MuCharge",3,-1.5,1.5,"q(#mu)/e");
  TauCharge=HConfig.GetTH1D(Name+"_TauCharge","TauCharge",3,-1.5,1.5,"q(#tau)/e");

  MuTauDR=HConfig.GetTH1D(Name+"_MuTauDR","MuTauDR",50,0.,5.,"#DeltaR(#mu,#tau_{h})");
  MuTauDPhi=HConfig.GetTH1D(Name+"_MuTauDPhi","MuTauDPhi",50,0.,3.2,"#Delta#phi(#mu,#tau_{h})");
  MuTauDEta=HConfig.GetTH1D(Name+"_MuTauDEta","MuTauDEta",50,-5.,5.,"#Delta#eta(#mu,#tau_{h})");
  MuTauDPt=HConfig.GetTH1D(Name+"_MuTauDPt","MuTauDPt",100,-100.,100.,"#Deltap_{T}(#mu,#tau_{h})/GeV");
  MuTauRelDPt=HConfig.GetTH1D(Name+"_MuTauRelDPt","MuTauRelDPt",100,-2.,2.,"#Deltap_{T}(#mu,#tau_{h})/p_{T}(#mu)");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}


void  HToTaumuTauh::Store_ExtraDist(){
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

 Extradist1d.push_back(&TauPt  );
 Extradist1d.push_back(&TauEta  );

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
}

void  HToTaumuTauh::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  double wobs=1;
  double w;
  if(!Ntp->isData()){w = Ntp->EvtWeight3D();}
  else{w=1;}

  // Apply Selection

  // Vertex
  unsigned int nGoodVtx=0;
  int selVertex = -1;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(selectVertex(i)){
    	if(selVertex == -1) selVertex = i; // selected vertex = first vertex (highest sum[pT^2]) to fulfill vertex requirements
    	nGoodVtx++;
    }
  }

  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  // Trigger
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;
  
  // Muon cuts
  std::vector<int> selectedMuonsId;
  selectedMuonsId.clear();
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if( selectMuon_Id(i,selVertex) ) {
		  selectedMuonsId.push_back(i);
	  }
  }
  value.at(NMuId)=selectedMuonsId.size();
  pass.at(NMuId)=(value.at(NMuId)>=cut.at(NMuId));

  std::vector<int> selectedMuons;	// full selection: ID and Kinematics
  selectedMuons.clear();
  for(unsigned i=0;i<selectedMuonsId.size();i++){
	  if( selectMuon_Kinematics(i)) {
		  selectedMuons.push_back(i);
	  }
  }
  value.at(NMuKin)=selectedMuons.size();
  pass.at(NMuKin)=(value.at(NMuKin)>=cut.at(NMuKin));
  int selMuon = (selectedMuons.size() > 0) ? selectedMuons.at(0) : -1;

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
  for(unsigned i=0; i < Ntp->NPFTaus(); i++){
	  if ( selectPFTau_Id(i) ){
		  selectedTausId.push_back(i);
	  }
  }
  value.at(NTauId)=selectedTausId.size();
  pass.at(NTauId)=(value.at(NTauId)>=cut.at(NTauId));

  std::vector<int> selectedTausIso;
  selectedTausIso.clear();
  for(unsigned i=0; i < selectedTausId.size(); i++){
	  if ( selectPFTau_Iso(i) ){
		  selectedTausIso.push_back(i);
	  }
  }
  value.at(NTauIso)=selectedTausIso.size();
  pass.at(NTauIso)=(value.at(NTauIso)>=cut.at(NTauIso));

  std::vector<int> selectedTaus;
  selectedTaus.clear();
  for(unsigned i=0; i < selectedTausIso.size(); i++){
	  if ( selectPFTau_Kinematics(i) ){
		  selectedTaus.push_back(i);
	  }
  }
  value.at(NTauKin)=selectedTaus.size();
  pass.at(NTauKin)=(value.at(NTauKin)>=cut.at(NTauKin));
  int selTau = (selectedTaus.size() > 0) ? selectedTaus.at(0) : -1;

  // Opposite charge
  if (selMuon != -1 && selTau != -1){
	  value.at(OppCharge) = Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau);
  }
  else value.at(OppCharge) = -9;
  pass.at(OppCharge) = (value.at(OppCharge) == cut.at(OppCharge));

  // Tri-lepton veto
  std::vector<int> triLepVetoMuons;
  triLepVetoMuons.clear();
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if( selectMuon_triLeptonVeto(i,selMuon,selVertex) ) {
		  triLepVetoMuons.push_back(i);
	  }
  }
  std::vector<int> triLepVetoElecs;
  triLepVetoElecs.clear();
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if( selectElectron_triLeptonVeto(i,selVertex) ) {
		  triLepVetoElecs.push_back(i);
	  }
  }
  value.at(TriLeptonVeto) = triLepVetoMuons.size() + triLepVetoElecs.size();
  pass.at(TriLeptonVeto) = (value.at(TriLeptonVeto) <= cut.at(TriLeptonVeto));



  bool status=AnalysisCuts(t,w,wobs); // true only if full selection passed

  ///////////////////////////////////////////////////////////
  // Add plots
  ///////////////////////////////////////////////////////////

  //////// plots filled before any cuts
  // Vertex plots
  NVtx.at(t).Fill(Ntp->NVtx(),w);
  for(unsigned int i=0;i<Ntp->NVtx();i++){
	VtxZ.at(t).Fill(Ntp->Vtx(i).z(),w);
	VtxRho.at(t).Fill(sqrt(Ntp->Vtx(i).x()*Ntp->Vtx(i).x() + Ntp->Vtx(i).y()*Ntp->Vtx(i).y()), w);
	VtxNdof.at(t).Fill(Ntp->Vtx_ndof(i), w);
	VtxIsfake.at(t).Fill(Ntp->Vtx_isFake(i), w);
  }
  NGoodVtx.at(t).Fill(nGoodVtx,w);

  //////// plots filled after Vertex selection: Object selection

  if(pass.at(TriggerOk) && pass.at(PrimeVtx)){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(	isTightMuon(i,selVertex) ){
			  MuDxy.at(t).Fill(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(selVertex)), w);
			  MuDz.at(t).Fill(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(selVertex)), w);
			  MuRelIso.at(t).Fill(Muon_RelIso(i), w);
		  }
	  }
  }

  //////// plots filled after muon ID selection: Muon Kinematics
  if(pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NMuId)){
	  for(unsigned i=0;i<selectedMuonsId.size();i++){
		  MuPt.at(t).Fill(Ntp->Muons_p4(i).Pt(), w);
		  MuEta.at(t).Fill(Ntp->Muons_p4(i).Eta(), w);
	  }
  }

  //////// plots filled after tau ID + Iso selection: Tau Kinematics
  if(pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NTauId) && pass.at(NTauIso)){
	  for(unsigned i=0;i<selectedTausIso.size();i++){
		  TauPt.at(t).Fill(Ntp->PFTau_p4(i).Pt(), w);
		  TauEta.at(t).Fill(Ntp->PFTau_p4(i).Eta(), w);
	  }
  }

  //////// plots filled after full muon and tau selection
  bool passedObjects = pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NMuId   ) && pass.at(NMuKin  ) && pass.at(NTauId  ) && pass.at(NTauIso ) && pass.at(NTauKin );
  if(passedObjects){
	  if(!pass.at(DiMuonVeto)){
		  // Investigate events discarded by the DiMuon Veto
		  if (Ntp->Muon_Charge(selMuon) == 1){
			  MuVetoDPtSelMuon.at(t).Fill( Ntp->Muons_p4(diMuonVetoMuonsNegative.at(0)).Pt() - Ntp->Muons_p4(selMuon).Pt(), w );
			  MuVetoDRTau.at(t).Fill( Ntp->Muons_p4(diMuonVetoMuonsNegative.at(0)).DeltaR(Ntp->PFTau_p4(selTau)), w);
		  }
		  else if (Ntp->Muon_Charge(selMuon) == -1){
			  MuVetoDPtSelMuon.at(t).Fill( Ntp->Muons_p4(diMuonVetoMuonsPositive.at(0)).Pt() - Ntp->Muons_p4(selMuon).Pt(), w );
			  MuVetoDRTau.at(t).Fill( Ntp->Muons_p4(diMuonVetoMuonsPositive.at(0)).DeltaR(Ntp->PFTau_p4(selTau)), w);
		  }
		  MuVetoInvM.at(t).Fill( (Ntp->Muons_p4(diMuonVetoMuonsPositive.at(0)) + Ntp->Muons_p4(diMuonVetoMuonsNegative.at(0))).M() , w);
		  MuVetoPtPositive.at(t).Fill( Ntp->Muons_p4(diMuonVetoMuonsPositive.at(0)).Pt(), w);
		  MuVetoPtNegative.at(t).Fill( Ntp->Muons_p4(diMuonVetoMuonsNegative.at(0)).Pt(), w);
	  }

	  if(pass.at(DiMuonVeto)){
		  // Mu-Tau correlations
		  MuTauDR    .at(t).Fill( Ntp->Muons_p4(selMuon).DeltaR(Ntp->PFTau_p4(selTau)), w );
		  MuTauDPhi  .at(t).Fill( Ntp->Muons_p4(selMuon).DeltaPhi(Ntp->PFTau_p4(selTau)), w );
		  MuTauDEta  .at(t).Fill( Ntp->Muons_p4(selMuon).Eta() - Ntp->PFTau_p4(selTau).Eta(), w );
		  MuTauDPt   .at(t).Fill( Ntp->Muons_p4(selMuon).Pt() - Ntp->PFTau_p4(selTau).Pt(), w );
		  MuTauRelDPt.at(t).Fill( (Ntp->Muons_p4(selMuon).Pt() - Ntp->PFTau_p4(selTau).Pt()) / Ntp->Muons_p4(selMuon).Pt() , w);
		  // Tri-lepton vetoes
		  NMuonTriLepVeto.at(t).Fill(triLepVetoMuons.size(), w);
		  NElecTriLepVeto.at(t).Fill(triLepVetoElecs.size(), w);
		  // lepton charge
		  MuCharge.at(t).Fill( Ntp->Muon_Charge(selMuon), w);
		  TauCharge.at(t).Fill( Ntp->PFTau_Charge(selMuon), w);
	  }
  }

  //////// plots filled after full selection
  if(status){
    NVtxFullSelection.at(t).Fill(Ntp->NVtx(),w);
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

///////// Vertices
bool HToTaumuTauh::selectVertex(unsigned i){
	if(verbose)std::cout << "selectVertex(unsigned i)" << std::endl;
	if( 	(fabs(Ntp->Vtx(i).z()) < 24)
			&& (Ntp->Vtx(i).Perp() < 2)
			&& (Ntp->Vtx_ndof(i) > 4)
			&& (Ntp->Vtx_isFake(i) == 0)
			) return true;
	return false;
}

///////// Muons

// isolation
double HToTaumuTauh::Muon_AbsIso(unsigned int i){
	return Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i));
}

double HToTaumuTauh::Muon_RelIso(unsigned int i){
	return (Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i)))/Ntp->Muons_p4(i).Pt();
}


// tightMuon without vertex constrains
bool HToTaumuTauh::isTightMuon(unsigned i){
	if(verbose)std::cout << "isTightMuon(unsigned i, unsigned i_vtx)" << std::endl;
	if(		Ntp->Muon_isGlobalMuon(i) &&
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

// tightMuon with vertex constrains
bool HToTaumuTauh::isTightMuon(unsigned i, unsigned i_vtx){
	if(			isTightMuon(i) &&
				fabs(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx)))<0.2 &&
				fabs(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx)))<0.5
				  ){
			return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_Id(unsigned i, unsigned vertex){
	if(	isTightMuon(i,vertex) &&
		fabs(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))) < cMu_dxy &&
		fabs(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))) < cMu_dz &&
		Muon_RelIso(i) < cMu_relIso
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_Kinematics(unsigned i){
	if(	Ntp->Muons_p4(i).Pt() >= cMu_pt &&
		fabs(Ntp->Muons_p4(i).Eta()) <= cMu_eta
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_diMuonVeto(unsigned i, unsigned i_vtx){
	if(	Ntp->Muons_p4(i).Pt() > 15.0 &&
		fabs(Ntp->Muons_p4(i).Eta()) < 2.4 &&
		Ntp->Muon_isPFMuon(i) &&
		Ntp->Muon_isGlobalMuon(i) &&
		Ntp->Muon_isTrackerMuon(i) &&
		Muon_RelIso(i) < 0.3 &&
		fabs(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < 0.2
		) {
	  return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_triLeptonVeto(unsigned i, int selectedMuon, unsigned i_vtx){
	if(	i != selectedMuon &&
		isTightMuon(i,i_vtx) &&
		fabs(dxy(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < cMu_dxy &&
		fabs(dz(Ntp->Muons_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < cMu_dz &&
		Muon_RelIso(i) < 0.3 &&
		Ntp->Muons_p4(i).Pt() > 10.0 &&
		fabs(Ntp->Muons_p4(i).Eta()) < 2.4
			){
			return true;
	}
	return false;
}


///////// Electrons

bool HToTaumuTauh::isLooseMVAElectron(unsigned i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));

//	TODO: uncomment as soon as Electron_MVA_discriminator is in Ntuple and Ntuple_Controller (was implemented already by Alex)
//	if(mvapt<20){
//		if(mvaeta<0.8 && Ntp->Electron_MVA_discriminator(i)>0.925){ // Cat. 1
//			return true;
//		}else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_discriminator(i)>0.915){ // Cat. 2
//			return true;
//		}else if(mvaeta>=1.479 && Ntp->Electron_MVA_discriminator(i)>0.965){ // Cat. 3
//			return true;
//		}
//	}else if(mvapt>=20){
//		if(mvaeta<0.8 && Ntp->Electron_MVA_discriminator(i)>0.905){ // Cat. 4
//			return true;
//		}else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_discriminator(i)>0.955){ // Cat. 5
//			return true;
//		}else if(mvaeta>=1.479 && Ntp->Electron_MVA_discriminator(i)>0.975){ // Cat. 6
//			return true;
//		}
//	}
//	return false;
	return false;
}

double HToTaumuTauh::Electron_RelIso(unsigned i){
	//	TODO: uncomment as soon as Electron_MVA_discriminator is in Ntuple and Ntuple_Controller (was implemented already by Alex)
	//return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();

	return 10.0;
}

bool HToTaumuTauh::selectElectron_triLeptonVeto(unsigned i, unsigned i_vtx){
	if ( 	Ntp->Electron_numberOfMissedHits(i) < 0.1 && // no missing hits
			!Ntp->Electron_HasMatchedConversions(i) &&
			fabs( dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(i_vtx)) ) < 0.2 &&
			fabs( dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(i_vtx)) ) < 0.045 &&
			isLooseMVAElectron(i) &&
			Electron_RelIso(i) < 0.3 &&
			Ntp->Electron_p4(i).Pt() > 10.0 &&
			fabs(Ntp->Electron_p4(i).Eta()) < 2.5
			){
		return true;
	}
	return false;
}

///////// Taus

bool HToTaumuTauh::selectPFTau_Id(unsigned i, unsigned i_muon){
	if ( 	Ntp->PFTau_p4(i).DeltaR(Ntp->Muons_p4(i_muon))  &&
			Ntp->PFTau_isHPSByDecayModeFinding(i) &&
			Ntp->PFTau_isHPSAgainstElectronsLoose(i) &&
			Ntp->PFTau_isHPSAgainstMuonTight(i)
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Iso(unsigned i){
	if ( 	Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(i)
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
