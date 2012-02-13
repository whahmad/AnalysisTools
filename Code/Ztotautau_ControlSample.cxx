#include "Ztotautau_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"

Ztotautau_ControlSample::Ztotautau_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
  ,channel(muontag)
{
  if(Get_Name().Contains("muontag")) channel=muontag;
  if(Get_Name().Contains("electrontag")) channel=muontag;  // not implemented yet
  if(Get_Name().Contains("rhotag")) channel=muontag;       // not implemented yet
  if(Get_Name().Contains("threepiontag")) channel=muontag; // not implemented yet

}

Ztotautau_ControlSample::~Ztotautau_ControlSample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Ztotautau_ControlSample::~Ztotautau_ControlSample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Ztotautau_ControlSample::~Ztotautau_ControlSample()" << std::endl;
}

void  Ztotautau_ControlSample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==hasTag)             cut.at(hasTag)=1;
    if(i==TagPt)              cut.at(TagPt)=18;
    if(i==TagIso)             cut.at(TagIso)=0.2;
    if(i==NJets)              cut.at(NJets)=1;
    if(i==JetPt)              cut.at(JetPt)=10;
    if(i==deltaPhi)           cut.at(deltaPhi)==TMath::Pi()*3.0/4.0;
    if(i==MET)                cut.at(MET)=40;
    if(i==MT)                 cut.at(MT)=30;
    if(i==PInBalance)         cut.at(PInBalance)=0.1;
    if(i==TauAvgMETPhi)       cut.at(TauAvgMETPhi)=0.75;
    if(i==ZMassV)             cut.at(ZMassV)=120;
    if(i==charge)             cut.at(charge)=0;
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //-----------
   else if(i==hasTag){
      title.at(i)="has Tag";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="has Tag (bool)";
      if(channel==muontag){
	title.at(i)="Number of Muons $=$";
	title.at(i)+=cut.at(hasTag);
	htitle=title.at(i);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="N_{Muons}";
      }
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
   else if(i==TagPt){
      title.at(i)="$P_{T}^{Tag}>$";
      title.at(i)+=cut.at(TagPt);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="P_{T}^{Tag} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPt_",htitle,20,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPt_",htitle,20,0,100,hlabel,"Events"));
    }
   else if(i==TagIso){
      title.at(i)="Relative Isolation (Tag) $<$";
      title.at(i)+=cut.at(TagIso);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Relative Isolation (Tag)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagIso_",htitle,40,0,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagIso_",htitle,40,0,1,hlabel,"Events"));
    }
   else if(i==NJets){
     title.at(i)="$Number of Jet>$";
     title.at(i)+=cut.at(NJets);
     title.at(i)+="";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="N_{Jets}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NJets_",htitle,11,-0.0,10.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NJets_",htitle,11,-0.5,10.5,hlabel,"Events"));
   }
   else if(i==JetPt){
     title.at(i)="$P_{T}^{Jet}>$";
     title.at(i)+=cut.at(JetPt);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="P_{T}^{Jet} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_JetPt_",htitle,20,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_JetPt_",htitle,20,0,200,hlabel,"Events"));
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
   else if(i==MT){
     title.at(i)="$M_{T} < $";
     title.at(i)+=cut.at(MT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_",htitle,40,0,200,hlabel,"Events"));
   }
   else if(i==PInBalance){
     title.at(i)="$|P_{T}^{\\mu-Jet}-P_{T}^{Jets}|/|P_{T}^{\\mu-Jet}+P_{T}^{Jets}|$";
     title.at(i)+=cut.at(PInBalance);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="|P_{T}^{#mu-Jet}-P_{T}^{Jets}|/|P_{T}^{#mu-Jet}+P_{T}^{Jets}|";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PInBalance_",htitle,40,0,2.0,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PInBalance_",htitle,40,0,2.0,hlabel,"Events"));
   }
   else if(i==TauAvgMETPhi){
     title.at(i)="$sin(\\phi_{\\mu}-\\phi_{E_{T}^{Miss}})$";
     title.at(i)+=cut.at(TauAvgMETPhi);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="sin(#phi_{#mu}-#phi_{E^{T}^{miss}})";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauAvgMETPhi_",htitle,20,-1,1,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauAvgMETPhi_",htitle,20,-1,1,hlabel,"Events"));
   }
   else if(i==deltaPhi){
      title.at(i)="$\\Delta\\phi(Tag,Jet) < $";
      title.at(i)+=cut.at(deltaPhi);
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta#phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,32,-TMath::Pi(),TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,32,-TMath::Pi(),TMath::Pi(),hlabel,"Events"));
    } 
   else if(i==ZMassV){
      title.at(i)="$M_{Z} $";
      title.at(i)+=cut.at(ZMassV);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{Z} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassV_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassV_",htitle,40,0,200,hlabel,"Events"));
    } 
   else if(i==charge){
      title.at(i)="$C_{\\mu}\\times 4\\sum_{i=0}^{ntracks}P_{i}C_{i}/\\sum_{i=0}^{ntracks}P_{i}$";
      title.at(i)+=cut.at(charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="C_{#mu}#times 4#sum_{i=0}^{ntracks}P_{i}C_{i}/#sum_{i=0}^{ntracks}P_{i}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,110,-10.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,110,-10.5,10.5,hlabel,"Events"));
    } 

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  Ztotautau_ControlSample::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
}

void  Ztotautau_ControlSample::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() A" << std::endl;
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;
  
  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() B" << std::endl;

  unsigned mu_idx(999),nmus(0);
  double mu_pt(0);
  if(channel==muontag){
    for(unsigned int i=0;i<Ntp->NMuons();i++){
      if(Ntp->isGoodMuon(i)){
	if(mu_pt<Ntp->Muons_p4(i).Pt()){mu_idx=i;mu_pt=Ntp->Muons_p4(i).Pt();}
	nmus++;
      }
    }
    if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() C" << std::endl;
    value.at(hasTag)=nmus;
    pass.at(hasTag)=(value.at(hasTag)==cut.at(hasTag));
    std::cout << nmus << std::endl;    

    value.at(TagPt)=mu_pt;
    pass.at(TagPt)=(value.at(TagPt)>=cut.at(TagPt));
    if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() C2 - " << mu_idx << " " << Ntp->NMuons() << std::endl;
    if(mu_idx!=999){value.at(TagIso) = (Ntp->Muon_emEt05(mu_idx) + Ntp->Muon_hadEt05(mu_idx) + Ntp->Muon_sumPt05(mu_idx))/Ntp->Muons_p4(mu_idx).Pt();}
    else{value.at(TagIso)=999;}
    pass.at(TagIso)=(value.at(TagIso)<=cut.at(TagIso));
    
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() d" << std::endl;
  unsigned int jet_idx(999),njets(0);
  double jet_pt(0);
  for(int i=0;i<Ntp->NPFJets();i++){
    if(verbose)std::cout << "jet loop " << i << " " << Ntp->NPFJets() << std::endl;
    if(Ntp->isGoodJet(i)){
      if(verbose)std::cout << "jet loop is good" << std::endl;
      if(jet_pt<Ntp->PFJet_p4(i).Pt()){jet_idx=i;jet_pt=Ntp->PFJet_p4(i).Pt();}
      njets++;
      if(verbose)std::cout << "jet loop is good end" << std::endl;
    }
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() e" << std::endl;
  value.at(NJets)=njets;
  pass.at(NJets)=(value.at(NJets)>=cut.at(NJets));

  value.at(JetPt)=jet_pt;
  pass.at(JetPt)=(value.at(JetPt)>=cut.at(JetPt));
    

  if(mu_idx!=999 && jet_idx!=999){
    value.at(deltaPhi)=Tools::DeltaPhi(Ntp->Muons_p4(mu_idx),Ntp->PFJet_p4(jet_idx));
  }
  else { 
    value.at(deltaPhi)=0;
  }
  pass.at(deltaPhi)=(fabs(value.at(deltaPhi))>=cut.at(deltaPhi));
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() f" << std::endl;

  value.at(MET)=Ntp->MET_et();
  pass.at(MET)=(value.at(MET)<cut.at(MET));

  if(mu_idx!=999){
    value.at(MT)=sqrt(2*(Ntp->MET_et())*Ntp->Muons_p4(mu_idx).Pt()*fabs(1-cos(Ntp->Muons_p4(mu_idx).Phi()-Ntp->MET_phi())));
  }
  else{
    value.at(MT)=999;
  }
  pass.at(MT)=(value.at(MT)<=cut.at(MT));

  TLorentzVector Z_lv(0,0,0,0);
  if(mu_idx!=999)Z_lv+=Ntp->Muons_p4(mu_idx);
  if(jet_idx!=999)Z_lv+=Ntp->PFJet_p4(jet_idx);
  value.at(ZMassV)=Z_lv.M();
  pass.at(ZMassV)=(value.at(ZMassV)<cut.at(ZMassV));
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() g" << std::endl;

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() h" << std::endl;
  if(jet_idx!=999 && mu_idx!=999){
    Ntp->Muons_p4(mu_idx);
    Ntp->PFJet_p4(jet_idx);
    double N=(((Ntp->PFJet_p4(jet_idx).Px()-Ntp->Muons_p4(mu_idx).Px())*Ntp->MET_et()*sin(Ntp->MET_phi()))-
	      ((Ntp->PFJet_p4(jet_idx).Py()-Ntp->Muons_p4(mu_idx).Py())*Ntp->MET_et()*cos(Ntp->MET_phi())));
    double D=sqrt(pow(Ntp->Muons_p4(mu_idx).Px()-Ntp->PFJet_p4(jet_idx).Px(),2.0)+pow(Ntp->Muons_p4(mu_idx).Py()-Ntp->PFJet_p4(jet_idx).Py(),2.0))*fabs(Ntp->MET_et());
    value.at(TauAvgMETPhi)=N/D;
    pass.at(TauAvgMETPhi)=(fabs(value.at(TauAvgMETPhi))>cut.at(TauAvgMETPhi));
    std::cout << "TauAvgMETPhi" << value.at(TauAvgMETPhi) << std::endl;
  }
  else{
    value.at(TauAvgMETPhi)=0;
    pass.at(TauAvgMETPhi)=false;
  }
  pass.at(TauAvgMETPhi)=true;

  if(jet_idx!=999 && mu_idx!=999){
    unsigned int jetmatch_idx=0;
    bool hasmatch=Ntp->muonhasJetMatch(mu_idx,jetmatch_idx);
    TLorentzVector LV_diff=Ntp->PFJet_p4(jet_idx);
    TLorentzVector LV_sum=Ntp->PFJet_p4(jet_idx);
    double pt_sum(Ntp->PFJet_p4(jet_idx).Pt()), pt_diff(Ntp->PFJet_p4(jet_idx).Pt());
    if(hasmatch){
      LV_diff+=Ntp->PFJet_p4(jetmatch_idx);
      LV_sum-=Ntp->PFJet_p4(jetmatch_idx);
      pt_diff-=Ntp->PFJet_p4(jetmatch_idx).Pt();
      pt_sum+=Ntp->PFJet_p4(jetmatch_idx).Pt();
    }
    else{
      LV_diff+=Ntp->Muons_p4(mu_idx);
      LV_sum-=Ntp->Muons_p4(mu_idx);
      pt_diff-=Ntp->Muons_p4(mu_idx).Pt();
      pt_sum+=Ntp->Muons_p4(mu_idx).Pt();
    }

    value.at(PInBalance)=fabs(LV_diff.Pt())/LV_sum.Pt();
    pass.at(PInBalance)=value.at(PInBalance)>cut.at(PInBalance);
  }
  else{
    value.at(PInBalance)=0;
    pass.at(PInBalance)=false;
  }
  pass.at(PInBalance)=true;
  // Momentum weighted charge is slow only do if nminus1
  // must be last cut
  bool charge_nminus1(true);
  for(unsigned int i=0;i<NCuts;i++){
    if(!pass.at(i) && i!=charge) charge_nminus1=false;
  }
  double PWeightedCharge(0),PTracks(0);
  if(charge_nminus1 && jet_idx!=999 && mu_idx!=999){
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    for(int i=0; i<PFJet_Track_idx.size();i++){
      unsigned int idx=PFJet_Track_idx.at(i);
      if(idx<Ntp->NTracks()){
	double pt=Ntp->Track_p4(idx).Pt();
	PWeightedCharge+=pt*pt*Ntp->Track_charge(idx);
	PTracks+=pt*pt;
      }
    }
    value.at(charge)=Ntp->Muon_Charge(mu_idx)*4.0*PWeightedCharge/PTracks;
    pass.at(charge)=(value.at(charge)>cut.at(charge));
  }
  else{
    value.at(charge)=0;
    pass.at(charge)=false;
  }
  pass.at(charge)=true;
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() i" << std::endl;
  /*
  if( Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex)!=-1 && Ntp->NTracks() >= Ntp->Muon_Track_idx(HighestPtMuonIndex)){
    value.at(charge) =Ntp->KFTau_Fit_charge(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex))*Ntp->Track_charge(Ntp->Muon_Track_idx(HighestPtMuonIndex));
    
    //  pass.at(TauPt)=(value.at(TauPt)>=cut.at(TauPt));
    pass.at(charge)=(Ntp->KFTau_Fit_charge(Ntp->KFTau_indexOfFitInfo(HighestPtTauIndex))*Ntp->Track_charge(Ntp->Muon_Track_idx(HighestPtMuonIndex)) == -1);
  }
  */

  double wobs(1),w(1);
  /*
  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
  }
  else{w=1;}
  std::cout << "w=" << w << std::endl;
  */
  bool status=AnalysisCuts(t,w,wobs); 
 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
  }
}




