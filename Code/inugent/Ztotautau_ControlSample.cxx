#include "Ztotautau_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"

Ztotautau_ControlSample::Ztotautau_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
  ,channel(muontag)
  ,jeteta(2.4)
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
    if(i==hasTag)             cut.at(hasTag)=1;
    if(i==TagPtmin)           cut.at(TagPtmin)=25;
    if(i==TagPtmax)           cut.at(TagPtmax)=38;
    if(i==TagIso)             cut.at(TagIso)=0.2;
    if(i==NJets)              cut.at(NJets)=1;
    if(i==JetPt)              cut.at(JetPt)=10;
    if(i==deltaPhi)           cut.at(deltaPhi)=TMath::Pi()*7.0/8.0;
    if(i==MET)                cut.at(MET)=40;
    if(i==MT)                 cut.at(MT)=30;
    if(i==etaq)               cut.at(etaq)=0.1;
    if(i==sumcosdeltaphi)     cut.at(sumcosdeltaphi)=0.75;
    if(i==ZMassmin)           cut.at(ZMassmin)=70;
    if(i==ZMassmax)           cut.at(ZMassmax)=100;
    if(i==HT)                 cut.at(HT)=200;
    if(i==charge)             cut.at(charge)=-1.0;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";
    if(i<10)c+="0";
    c+=i;
  
    if(i==TriggerOk){
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
   else if(i==TagPtmin){
      title.at(i)="$P_{T}^{Tag}>$";
      title.at(i)+=cut.at(TagPtmin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="P_{T}^{Tag} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtmin_",htitle,50,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtmin_",htitle,50,0,100,hlabel,"Events"));
    }
   else if(i==TagPtmax){
     title.at(i)="$P_{T}^{Tag}<$";
     title.at(i)+=cut.at(TagPtmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="P_{T}^{Tag} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtmax_",htitle,50,0,100,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtmax_",htitle,50,0,100,hlabel,"Events"));
   }
   else if(i==TagIso){
      title.at(i)="Relative Isolation (Tag) $<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(TagIso));
      title.at(i)+=buffer;
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
   else if(i==HT){
     title.at(i)="$H_{T} < $";
     title.at(i)+=cut.at(HT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="H_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HT_",htitle,50,0,500,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HT_",htitle,50,0,500,hlabel,"Events"));
   }
   else if(i==etaq){
     title.at(i)="$q_{\\mu}|\\eta|$";
     char buffer[50];
     sprintf(buffer,"%5.2f",cut.at(etaq));
     title.at(i)+=buffer;
     title.at(i)+="";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="q_{#mu}|#eta|";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_etaq_",htitle,40,-2.0,2.0,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_etaq_",htitle,40,-2.0,2.0,hlabel,"Events"));
   }
   else if(i==sumcosdeltaphi){
     title.at(i)="$\\sum cos(\\delta\\phi)$";
     char buffer[50];
     sprintf(buffer,"%5.2f",cut.at(sumcosdeltaphi));
     title.at(i)+=buffer;
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="#sum cos(#delta#phi)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_sumcosdeltaphi_",htitle,20,-2,2,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_sumcosdeltaphi_",htitle,20,-2,2,hlabel,"Events"));
   }
   else if(i==deltaPhi){
      title.at(i)="$\\Delta\\phi(Tag,Jet) < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta#phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
    } 
   else if(i==ZMassmin){
      title.at(i)="$M_{Z}< $";
      title.at(i)+=cut.at(ZMassmin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{Z} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmin_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmin_",htitle,40,0,200,hlabel,"Events"));
    } 
   else if(i==ZMassmax){
     title.at(i)="$M_{Z}< $";
     title.at(i)+=cut.at(ZMassmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{Z} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
   }
   else if(i==charge){
      title.at(i)="$C_{\\mu}\\times\\sum_{i=0}^{3}C_{i}|_{P_{T}^{max}}$";
      title.at(i)+=cut.at(charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="C_{#mu}#times #sum_{i=0}^{3}C_{i}|_{P_{T}^{max}}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,11,-10.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,11,-10.5,10.5,hlabel,"Events"));
    } 

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  TagEtaPT=HConfig.GetTH2D(Name+"_TagEtaPT","TagEtaPT",25,0,2.5,50,0,50,"#eta","P_{T}^{Tag}");

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
 Extradist2d.push_back(&TagEtaPT);
}

void  Ztotautau_ControlSample::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() A" << std::endl;

  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("IsoMu24"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);
  
  // Apply Selection
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

    value.at(TagPtmin)=mu_pt;
    pass.at(TagPtmin)=(value.at(TagPtmin)>cut.at(TagPtmin));
    value.at(TagPtmax)=mu_pt;
    pass.at(TagPtmax)=(value.at(TagPtmax)<cut.at(TagPtmax));

    if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() C2 - " << mu_idx << " " << Ntp->NMuons() << std::endl;
    if(mu_idx!=999){value.at(TagIso) = (Ntp->Muon_emEt05(mu_idx) + Ntp->Muon_hadEt05(mu_idx) + Ntp->Muon_sumPt05(mu_idx))/Ntp->Muons_p4(mu_idx).Pt();}
    else{value.at(TagIso)=999;}
    pass.at(TagIso)=(value.at(TagIso)<=cut.at(TagIso));
    
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() d" << std::endl;
  unsigned int jet_idx(999),njets(0);
  double jet_pt(0);
  if(mu_idx!=999){
    for(int i=0;i<Ntp->NPFJets();i++){
      if(Tools::dr(Ntp->Muons_p4(mu_idx),Ntp->PFJet_p4(i))>0.4){
	if(jeteta>fabs(Ntp->PFJet_p4(i).Eta())){
	  if(jet_pt<Ntp->PFJet_p4(i).Pt()){jet_idx=i;jet_pt=Ntp->PFJet_p4(i).Pt();}
	  njets++;
	}
      }
    }
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() e" << std::endl;
  value.at(NJets)=njets;
  pass.at(NJets)=(value.at(NJets)>=cut.at(NJets));

  value.at(JetPt)=jet_pt;
  pass.at(JetPt)=(value.at(JetPt)>=cut.at(JetPt));
    

  if(mu_idx!=999 && jet_idx!=999){
    value.at(deltaPhi)=fabs(Tools::DeltaPhi(Ntp->Muons_p4(mu_idx),Ntp->PFJet_p4(jet_idx)));
  }
  else { 
    value.at(deltaPhi)=0;
  }
  pass.at(deltaPhi)=(fabs(value.at(deltaPhi))>=cut.at(deltaPhi));


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
  value.at(ZMassmin)=Z_lv.M();
  pass.at(ZMassmin)=(value.at(ZMassmin)>cut.at(ZMassmin));
  value.at(ZMassmax)=Z_lv.M();
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() g" << std::endl;

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() h" << std::endl;
  if(jet_idx!=999 && mu_idx!=999){
    value.at(sumcosdeltaphi)=cos(Tools::DeltaPhi(Ntp->Muons_p4(mu_idx),Ntp->PFJet_p4(jet_idx)))+cos(Tools::DeltaPhi(Ntp->Muons_p4(mu_idx),Ntp->MET_phi()));
    pass.at(sumcosdeltaphi)=(fabs(value.at(sumcosdeltaphi))>cut.at(sumcosdeltaphi));
  }
  else{
    value.at(sumcosdeltaphi)=0;
    pass.at(sumcosdeltaphi)=false;
  }
  pass.at(sumcosdeltaphi)=true;

  if(mu_idx!=999){
    value.at(etaq)=Ntp->Muon_Charge(mu_idx)*fabs(Ntp->Muons_p4(mu_idx).Eta());
    pass.at(etaq)=value.at(etaq)>cut.at(etaq);
  }
  else{
    value.at(etaq)=0;
    pass.at(etaq)=false;
  }
  pass.at(etaq)=true;

  if(jet_idx!=999 && mu_idx!=999){
    value.at(HT)=0;
    value.at(HT)+=Ntp->PFJet_p4(jet_idx).Pt();
    value.at(HT)+=Ntp->Muons_p4(mu_idx).Pt();
    value.at(HT)+=Ntp->MET_et();
    for(int i=0;i<Ntp->NPFJets();i++){
      if(i!=jet_idx){
	if(Ntp->isGoodJet(i))value.at(HT)+=Ntp->PFJet_p4(jet_idx).Pt();
      }
    }
  }
  pass.at(HT)=value.at(HT)<cut.at(HT);
  // Momentum weighted charge is slow only do if nminus1
  // must be last cut
  int charge_nminus1(0);
  for(unsigned int i=0;i<NCuts;i++){
    if(!pass.at(i) && i!=charge) charge_nminus1++;
  }
  value.at(charge)=0;
  if(charge_nminus1<=1 && jet_idx!=999 && mu_idx!=999){
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    std::vector<TLorentzVector> track_p4;
    std::vector<float> track_charge;
    for(int i=0; i<PFJet_Track_idx.size();i++){ 
      if(PFJet_Track_idx.at(i)<Ntp->NTracks() && PFJet_Track_idx.at(i)>=0){
	track_p4.push_back(Ntp->Track_p4(PFJet_Track_idx.at(i)));
	track_charge.push_back(Ntp->Track_charge(PFJet_Track_idx.at(i)));
      }
    }
    for(int i=0; i<track_p4.size();i++){
      double JetCharge(0);
      int ntracks(0);
      for(int j=0; j<track_p4.size();j++){
	if(track_p4.at(i).Pt()<=track_p4.at(j).Pt()){
	  ntracks++;
	  JetCharge+=track_charge.at(j);
	}
      }
      if(ntracks==3){
	value.at(charge)=JetCharge*Ntp->Muon_Charge(mu_idx);
	break;
      }
    }
    pass.at(charge)=(fabs(value.at(charge)-cut.at(charge))<0.5);
  }  
  else{
    value.at(charge)=0;
    pass.at(charge)=false;
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() i" << std::endl;


  double wobs(1),w(1);
  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
  }
  else{w=1;}
  std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
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
    if(mu_idx!=999) TagEtaPT.at(t).Fill(fabs(Ntp->Muons_p4(mu_idx).Eta()),Ntp->Muons_p4(mu_idx).Pt(),w);

  }
}




