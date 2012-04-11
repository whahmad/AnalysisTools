#include "Ztotautau_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

Ztotautau_ControlSample::Ztotautau_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
  ,channel(muontag)
  ,jeteta(2.4)
{
  if(Get_Name().Contains("muontag")) channel=muontag;
  if(Get_Name().Contains("electrontag")) channel=muontag;  // not implemented yet
  if(Get_Name().Contains("rhotag")) channel=muontag;       // not implemented yet
  if(Get_Name().Contains("threepiontag")) channel=muontag; // not implemented yet
  //  verbose=true;
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
    if(i==ZMassmin)           cut.at(ZMassmin)=60;
    if(i==ZMassmax)           cut.at(ZMassmax)=90;
    if(i==HT)                 cut.at(HT)=200;
    if(i==charge)             cut.at(charge)=0.0;
    if(i==MaxTracksinJet)     cut.at(MaxTracksinJet)=5.0;
    if(i==MinTracksinJet)     cut.at(MinTracksinJet)=1.0;
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
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NJets_",htitle,21,-0.0,20.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NJets_",htitle,21,-0.5,20.5,hlabel,"Events"));
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,21,-10.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,21,-10.5,10.5,hlabel,"Events"));
    } 
   else if(i==MaxTracksinJet){
     title.at(i)="$N_{Tracks}^{Jet}<$";
     title.at(i)+=cut.at(MaxTracksinJet);
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="N_{Tracks}^{Jet}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MaxTracksinJet_",htitle,11,-0.5,10.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MaxTracksinJet_",htitle,11,-0.5,10.5,hlabel,"Events"));
   }
   else if(i==MinTracksinJet){
     title.at(i)="$N_{Tracks}^{Jet}>$";
     title.at(i)+=cut.at(MinTracksinJet);
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="N_{Tracks}^{Jet}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MinTracksinJet_",htitle,11,-0.5,10.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MinTracksinJet_",htitle,11,-0.5,10.5,hlabel,"Events"));
     distindx.at(i)=true;
     Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_MinTracksinJet_","P_{T,Track} (N-1 Distribution)",100,0,200,"P_{T,Track} (GeV)","Events");
     Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_MinTracksinJet_","P_{T,Track} (Accumulative Distribution)",100,0,200,"P_{T,Track} (GeV)","Events");
   }

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  TagEtaPT=HConfig.GetTH2D(Name+"_TagEtaPT","TagEtaPT",25,0,2.5,50,0,50,"|#eta|","P_{T}^{Tag}");

  TauCandFound =HConfig.GetTH1D(Name+"_TauCandFound","TauCandFound",2,-0.5,1.5,"Found #tau Candidate (bool)","Events");
  TauCandEtaPhi =HConfig.GetTH2D(Name+"_TauCandPtEta","TauCandPtEta",25,0,2.5,40,0,200,"|#eta|","#phi (rad)");

  TauCandPhi =HConfig.GetTH1D(Name+"_TauCandPhi","TauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  TauCandPhiRes =HConfig.GetTH1D(Name+"_TauCandPhiRes","TauCandPhiRes",100,-0.2,0.2,"#sigma(#phi_{KF-#tau}) (rad)","Events");
  TauCandEta =HConfig.GetTH1D(Name+"_TauCandEta","TauCandEta",25,-2.5,2.5,"#eta_{KF-#tau}","Events");
  TauCandEtaRes =HConfig.GetTH1D(Name+"_TauCandEtaRes","TauCandEtaRes",100,-0.2,0.2,"#sigma(#eta_{KF-#tau}) ","Events");
  TauCandE =HConfig.GetTH1D(Name+"_TauCandE","TauCandE",40,0,200,"E_{KF-#tau} (GeV)","Events");
  TauCandERes =HConfig.GetTH1D(Name+"_TauCandERes","TauCandERes",50,-50,50,"#sigma(E_{KF-#tau}) (GeV)","Events");

  TauCandMass =HConfig.GetTH1D(Name+"_TauCandMass","TauCandMass",100,0,2,"M(KF-#tau) (GeV)","Events");
  TauCandPiMass =HConfig.GetTH1D(Name+"_TauCandPiMass","TauCandPiMass",100,0,1,"M(#pi) (GeV)","Events");
  TauCandNuMass =HConfig.GetTH1D(Name+"_TauCandNuMass","TauCandNuMass",110,-0.1,1,"M(#nu) (GeV)","Events");

  TauSolutionResult = HConfig.GetTH1D(Name+"_TauSolutionResult","TauSolutionResult",3,-1.5,1.5,"Solution: -:R<0:+","Events");
  EstimatedTauE =HConfig.GetTH1D(Name+"_TauCandE","TauCandE",40,0,200,"E_{Est-#tau} (GeV)","Events");
  EstimatedTauPhi =HConfig.GetTH1D(Name+"_TauCandPhi","TauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{Est-#tau} (rad)","Events");
  EstimatedTauEta  =HConfig.GetTH1D(Name+"_TauCandEta","TauCandEta",25,-2.5,2.5,"#eta_{Est-#tau}","Events");

  EstimatedTauERes =HConfig.GetTH1D(Name+"_TauEstimatedERes","TauEstimatedERes",40,-50,50,"#sigma(E_{Est-#tau}) (GeV)","Events");
  EstimatedTauPhiRes =HConfig.GetTH1D(Name+"_TauEstimatedPhiRes","TauEstimatedPhiRes",100,-0.2,0.2,"#sigma(#phi_{Est-#tau}) (rad)","Events");
  EstimatedTauEtaRes =HConfig.GetTH1D(Name+"_TauEstimatedEtaRes","TauEstimatedEtaRes",100,-0.2,0.2,"#sigma(#eta_{Est-#tau}) ","Events");

  EstimatedTauDirPhiRes=HConfig.GetTH1D(Name+"_TauCandDirPhiRes","TauCandDirPhiRes",100,-0.2,0.2,"#sigma(#phi_{Direction-#tau}) (rad)","Events");
  EstimatedTauDirEtaRes=HConfig.GetTH1D(Name+"_TauCandDirEtaRes","TauCandDirEtaRes",100,-0.2,0.2,"#sigma(#eta_{Direction-#tau}) (rad)","Events");

  KFTau_Fit_chiprob=HConfig.GetTH1D(Name+"_KFTau_Fit_prob","KFTau_Fit_prob",25,0,1,"Kinematic Fit Probability","Events");
  KFTau_Fit_a1mass=HConfig.GetTH1D(Name+"_KFTau_Fit_a1mass","KFTau_Fit_a1mass",25,0,2.5,"Kinematic Fit a_{1} Mass","Events");
  KFTau_Fit_chi2=HConfig.GetTH1D(Name+"_KFTau_Fit_chi2","KFTau_Fit_chi",25,0,25,"Kinematic Fit #chi^{2}","Events");
  KFTau_Fit_ndf=HConfig.GetTH1D(Name+"_KFTau_Fit_ndf","KFTau_Fit_ndf",25,0,25,"Kinematic Fit ndf","Events");
  KFTau_Fit_ambiguity=HConfig.GetTH1D(Name+"_KFTau_Fit_ambiguity","KFTau_Fit_ambiguity",5,-2.5,2.5,"Kinematic Fit Ambiguity","Events");
  KFTau_Fit_csum=HConfig.GetTH1D(Name+"_KFTau_Fit_csum","KFTau_Fit_csum",100,0,1.0,"Kinematic Fit csum","Events");
  KFTau_Fit_iterations=HConfig.GetTH1D(Name+"_KFTau_Fit_iterations","KFTau_Fit_iterations",25,0,25,"Kinematic Fit Iterations","Events");
  KFTau_Fit_TauEnergyFraction=HConfig.GetTH1D(Name+"_KFTau_Fit_TauEnergyFraction","KFTau_Fit_TauEnergyFraction",50,0,5.0,"Transverse Energy Fration,","Events");
  KFTau_Fit_PV_PV_significance=HConfig.GetTH1D(Name+"_ KFTau_Fit_PV_PV_significance"," KFTau_Fit_PV_PV_significance",50,0,10,"Kinematic Fit PV-PV_{Rotated} Significance","Events");
  KFTau_Fit_SV_PV_significance=HConfig.GetTH1D(Name+"_ KFTau_Fit_SV_PV_significance","KFTau_Fit_SV_PV_significance",50,0,10,"Kinematic Fit SV-PV_{Rotated} Significance","Events");


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

 Extradist1d.push_back(&TauCandFound);
 Extradist2d.push_back(&TauCandEtaPhi); 

 Extradist1d.push_back(&TauCandPhi);
 Extradist1d.push_back(&TauCandPhiRes);
 Extradist1d.push_back(&TauCandEta);
 Extradist1d.push_back(&TauCandEtaRes);
 Extradist1d.push_back(&TauCandE);
 Extradist1d.push_back(&TauCandERes);

 Extradist1d.push_back(&TauCandMass);
 Extradist1d.push_back(&TauCandPiMass);
 Extradist1d.push_back(&TauCandNuMass);

 Extradist1d.push_back(&TauSolutionResult);
 Extradist1d.push_back(&EstimatedTauE);
 Extradist1d.push_back(&EstimatedTauPhi);
 Extradist1d.push_back(&EstimatedTauEta);

 Extradist1d.push_back(&EstimatedTauERes);
 Extradist1d.push_back(&EstimatedTauPhiRes);
 Extradist1d.push_back(&EstimatedTauEtaRes);

 Extradist1d.push_back(&EstimatedTauDirEtaRes);
 Extradist1d.push_back(&EstimatedTauDirPhiRes);

 Extradist1d.push_back(&KFTau_Fit_chiprob); 
 Extradist1d.push_back(&KFTau_Fit_a1mass);
 Extradist1d.push_back(&KFTau_Fit_chi2);
 Extradist1d.push_back(&KFTau_Fit_ndf);
 Extradist1d.push_back(&KFTau_Fit_ambiguity);
 Extradist1d.push_back(&KFTau_Fit_csum);
 Extradist1d.push_back(&KFTau_Fit_iterations);
 Extradist1d.push_back(&KFTau_Fit_TauEnergyFraction);
 Extradist1d.push_back(&KFTau_Fit_PV_PV_significance);
 Extradist1d.push_back(&KFTau_Fit_SV_PV_significance);

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
  ///////////////////////////////////////
  if(jet_idx!=999){
    value.at(MaxTracksinJet)=Ntp->PFJet_Track_idx(jet_idx).size();
    pass.at(MaxTracksinJet)=true;//cut.at(MaxTracksinJet)>value.at(MaxTracksinJet); 
    value.at(MinTracksinJet)=Ntp->PFJet_Track_idx(jet_idx).size();
    pass.at(MinTracksinJet)=cut.at(MinTracksinJet)<value.at(MinTracksinJet);
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    for(int i=0; i<PFJet_Track_idx.size();i++){
      dist.at(MinTracksinJet).push_back(Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt());
    }
  }

  // Momentum weighted charge is slow only do if nminus1
  // must be last cut
  value.at(charge)=0;
  if(jet_idx!=999 && mu_idx!=999){
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
      if((ntracks==3 && track_p4.size()>=3) || (ntracks==track_p4.size() && track_p4.size()<3)){
	value.at(charge)=JetCharge+Ntp->Muon_Charge(mu_idx);
	break;
      }
    }
    pass.at(charge)=(fabs(value.at(charge)-cut.at(charge))<0.5);
  }  
  else{
    value.at(charge)=-10.0;
    pass.at(charge)=false;
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() i" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////
  // QCD Control sample
  if(!pass.at(charge) && fabs(value.at(charge))==1){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(charge)=true;
    }
  }

  double wobs(1),w(1);
  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
  }
  else{w=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
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

 
    /////////////////////////////////////////
    //
    // Analyze Tau
    //
    unsigned int tau_idx(999);
    if(jet_idx!=999){
      double mydr=0.5;
      for(unsigned i=0;i<Ntp->NKFTau();i++){
        if(Ntp->isGoodKFTau(i)){
          if(Tools::dr(Ntp->KFTau_TauVis_p4(i),Ntp->PFJet_p4(jet_idx))<mydr){
            tau_idx=i;
	    mydr=Tools::dr(Ntp->KFTau_TauVis_p4(i),Ntp->PFJet_p4(jet_idx));
          }
        }
      }
      
      if(tau_idx!=999 ){
        TauCandFound.at(t).Fill(1,w);
        TauCandEtaPhi.at(t).Fill(fabs(Ntp->KFTau_TauFit_p4(tau_idx).Eta()),Ntp->KFTau_TauFit_p4(tau_idx).Pt(),w);

        TauCandPhi.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).Phi(),w);
        TauCandEta.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).Eta(),w);
        TauCandE.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).E(),w);
	
	KFTau_Fit_chiprob.at(t).Fill(Ntp->KFTau_Fit_Chi2Prob(tau_idx),w);
	KFTau_Fit_a1mass.at(t).Fill(Ntp->KFTau_Fit_RefitVisibleMass(tau_idx),w);
	KFTau_Fit_chi2.at(t).Fill(Ntp->KFTau_Fit_chi2(tau_idx),w);
	KFTau_Fit_ndf.at(t).Fill(Ntp->KFTau_Fit_ndf(tau_idx),w);
	KFTau_Fit_ambiguity.at(t).Fill(Ntp->KFTau_Fit_ambiguity(tau_idx),w);
	KFTau_Fit_csum.at(t).Fill(Ntp->KFTau_Fit_csum(tau_idx),w);
	KFTau_Fit_iterations.at(t).Fill(Ntp->KFTau_Fit_iterations(tau_idx),w);
	KFTau_Fit_TauEnergyFraction.at(t).Fill(Ntp->KFTau_Fit_TauEnergyFraction(tau_idx),w);
	KFTau_Fit_PV_PV_significance.at(t).Fill(Ntp->KFTau_Fit_PV_PV_significance(tau_idx),w);
	KFTau_Fit_SV_PV_significance.at(t).Fill(Ntp-> KFTau_Fit_SV_PV_significance(tau_idx),w);
      }
      else{
	TauCandFound.at(t).Fill(0.0,w);
      }
    }
  }
  if(id==DataMCType::Signal){
    unsigned int mcBoson_idx,mctau_idx;
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,mcBoson_idx,TauDecay::JAK_A1_3PI,mctau_idx)){
      TLorentzVector MCTau_LV(0,0,0,0),X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(mctau_idx);i++){
	if(abs(Ntp->MCTauandProd_pdgid(mctau_idx,i))==abs(PdtPdgMini::tau_minus) && Ntp->MCTau_JAK(mctau_idx)==TauDecay::JAK_A1_3PI){
	  MCTau_LV=Ntp->MCTauandProd_p4(mctau_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(mctau_idx,i))==abs(PdtPdgMini::pi_plus) ||
		abs(Ntp->MCTauandProd_pdgid(mctau_idx,i))==abs(PdtPdgMini::pi0)
		){
	  X_LV+=Ntp->MCTauandProd_p4(mctau_idx,i);
	  if(verbose)std::cout << "MC pi: " << Ntp->MCTauandProd_pdgid(mctau_idx,i) << " 4-vec: " << Ntp->MCTauandProd_p4(mctau_idx,i).E()
			       << " " << Ntp->MCTauandProd_p4(mctau_idx,i).Phi() 
			       << " " << Ntp->MCTauandProd_p4(mctau_idx,i).Theta() 
			       << " " << Ntp->MCTauandProd_p4(mctau_idx,i).M() << std::endl;
	}
      
      }
      unsigned int tau_idx(999);
      double mydr=4.0;
      for(unsigned i=0;i<Ntp->NKFTau();i++){
	//if(Ntp->isGoodKFTau(i)){
	  if(Tools::dr(Ntp->KFTau_TauVis_p4(i),MCTau_LV)<mydr){
	    tau_idx=i;
	    mydr=Tools::dr(Ntp->KFTau_TauVis_p4(i),MCTau_LV);
	  }
	  //}
      }
      if(tau_idx!=999 ){
        TLorentzVector TauVisInput(0,0,0,0);
        TLorentzVector TauSolution(0,0,0,0);
	TLorentzVector TauSolution1;
	TLorentzVector TauSolution2;

	double Enu1(-999),Enu2(-999),R(0);
        TVector3 TauDirection=Ntp->KFTau_InitialSecondaryVtx(tau_idx)-Ntp->KFTau_ReducedVtx(); 
	unsigned int npi=0;
	
        for(unsigned int i=0;i<Ntp->KFTau_NDaughter(tau_idx);i++){
          if(abs(Ntp->KFTau_Daughter_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_plus)){
            TauCandMass.at(t).Fill(Ntp->KFTau_Daughter_par(tau_idx,i,Ntuple_Controller::KFTau_m),w);
	  }
	  else if(abs(Ntp->KFTau_Daughter_pdgid(tau_idx,i))==abs(PdtPdgMini::nu_tau)){
            TauCandNuMass.at(t).Fill(Ntp->KFTau_Daughter_par(tau_idx,i,Ntuple_Controller::KFTau_m),w);
          }
          else if(abs(Ntp->KFTau_Daughter_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus)){
            TauCandPiMass.at(t).Fill(Ntp->KFTau_Daughter_par(tau_idx,i,Ntuple_Controller::KFTau_m),w);
            TLorentzVector pi(Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_px),
			      Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_py),
			      Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_pz),
			      sqrt(pow(Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_m),2.0)+
				   pow(Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_px),2.0)+
				   pow(Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_py),2.0)+
				   pow(Ntp->KFTau_Daughter_inputpar(tau_idx,i,Ntuple_Controller::KFTau_pz),2.0)));
	    if(verbose)std::cout << "pi 4-vec E " << pi.E() << " Phi " <<pi.Phi() << " Theta " << pi.Theta() << " M " << pi.M() << std::endl;
            TauVisInput+=pi;
	    npi++;
          }
        }
	if(verbose)std::cout << "npi: " << npi << std::endl;
        if(npi==3){
	  TauVisInput=X_LV;
	  TauDirection=MCTau_LV.Vect();
	  double phi(TauDirection.Phi()),theta(TauDirection.Theta());
	  TLorentzVector TauVis1=TauVisInput;
	  TLorentzVector TauVis2=TauVisInput;
          TauVisInput.RotateZ(-phi);
          TauVisInput.RotateY(-theta);
	  double Enu(0), Ea1(TauVisInput.E()),Pz(TauVisInput.Pz()),mtau(PDG_Var::Tau_mass()),pt(TauVisInput.Pt()),ma1(TauVisInput.M());
          unsigned loop(0);
          bool Enuok=false;
          for(unsigned int loop=0;loop<10;loop++){
            if(loop>0){
              double factor=0.999;
              Ea1=Ea1*Ea1-pt*pt+factor*factor*pt*pt;
              pt*=0.999;
            }
	    
	  	    
	    std::cout << "Tau E " <<  MCTau_LV.E() << "Tau phi " << MCTau_LV.Phi() << " Tau Theta " <<MCTau_LV.Theta() << std::endl; 
	    std::cout << "Tau Dir  phi " << phi << " theta " << theta << std::endl;
	    std::cout << "TauVis E " <<  TauVisInput.E() << "Tau phi " << TauVis1.Phi() << " Tau Theta " <<TauVis1.Theta() << std::endl;
	    std::cout << "True Vis " << X_LV.E() << " phi " << X_LV.Phi() << " Theta " << X_LV.Theta() << std::endl;  

	    unsigned int method=2;
	    // Method 1: Solve for neutrino Pz by rotation
	    if(method==1){
	      double a=1-Pz*Pz/(Ea1*Ea1);
	      double K=(mtau*mtau-ma1*ma1-2.0*pt*pt);
	      double b=-K*Pz/(Ea1*Ea1);
	      double c=pt*pt-K*K/(4.0*Ea1*Ea1);
	      R=b*b-4.0*a*c;
	      std::cout << "R: " << R << std::endl;
	      if(R<0){if(verbose)std::cout << "R is negative " << R << "setting R to 0." << std::endl;R=0;}
	      double Pnuz1=(-b-sqrt(R))/(2.0*a);
	      double Pnuz2=(-b+sqrt(R))/(2.0*a);
	      Enu1=sqrt(Pnuz1*Pnuz1+pt*pt);
	      Enu2=sqrt(Pnuz2*Pnuz2+pt*pt);
	      if(Pnuz1<0)Enu1=pt;
	      if(Pnuz2<0)Enu2=pt;
	      std::cout << " Ea1 " << Ea1 << " Pz " << Pz << " ma1 " << ma1 << " pt " << pt << " K " << K << " mtau " << mtau << std::endl;
	      std::cout << "a " << a << " b " << b << " c " << c << std::endl;
	      std::cout << "R: " << R << " -b/2a " << (-b)/(2*a) << std::endl;
	      TLorentzVector Neutrino1(-TauVisInput.Px(),-TauVisInput.Py(),sqrt(Enu1*Enu1-pt*pt),Enu1);
	      Neutrino1.RotateY(theta);
	      Neutrino1.RotateZ(phi);
	      TauSolution1=TauVis1+Neutrino1;
              TLorentzVector Neutrino2(-TauVisInput.Px(),-TauVisInput.Py(),sqrt(Enu2*Enu2-pt*pt),Enu2);
              Neutrino2.RotateY(theta);
              Neutrino2.RotateZ(phi);
              TauSolution2=TauVis2+Neutrino2;
	    }

	    // Method 2: Solve for neutrino E by rotation
	    if(method==2){
	      double a=1-Ea1*Ea1/(Pz*Pz);
	      double K=(mtau*mtau-ma1*ma1-2*pt*pt);
	      double b=K*Ea1/(Pz*Pz);
	      double c=-(pt*pt+K*K/(4*Pz*Pz));
	      R=b*b-4*a*c;
	      std::cout << "R: " << R << std::endl;
	      if(R<0){if(verbose)std::cout << "R is negative " << R << "setting R to 0." << std::endl;R=0;}
	      Enu1=(-b-sqrt(R))/(2.0*a);
	      Enu2=(-b+sqrt(R))/(2.0*a);
	      std::cout << " Ea1 " << Ea1 << " Pz " << Pz << " ma1 " << ma1 << " pt " << pt << " K " << K << " mtau " << mtau << std::endl;
	      std::cout << "a " << a << " b " << b << " c " << c << std::endl;
	      std::cout << "R: " << R << " -b/2a " << (-b)/(2*a) << std::endl;
              TLorentzVector Neutrino1(-TauVisInput.Px(),-TauVisInput.Py(),sqrt(Enu1*Enu1-pt*pt),Enu1);
              Neutrino1.RotateY(theta);
              Neutrino1.RotateZ(phi);
	      TLorentzVector Neutrino2(-TauVisInput.Px(),-TauVisInput.Py(),sqrt(Enu2*Enu2-pt*pt),Enu2);
	      Neutrino2.RotateY(theta);
	      Neutrino2.RotateZ(phi);
	      TauSolution2=TauVis2+Neutrino2;
	    }
	    // Method 3: Solve for neutrino Pz without rotating frames
            if(method==3){
	      TVector3 ntau=TauDirection;
	      if(ntau.Mag()!=0)ntau*=1/ntau.Mag();
	      Pz=ntau.Dot(TauVis1.Vect());
	      TVector3 Pz_vec=ntau;
	      Pz_vec*=Pz;
	      TVector3 Pt_vec=TauVis1.Vect();
	      Pt_vec-=Pz_vec;
	      pt=Pt_vec.Mag();
	      ///////////////////////////////////////////////////////////////////////////////////////////////////////
              double a=1-Pz*Pz/(Ea1*Ea1);
              double K=(mtau*mtau-ma1*ma1-2.0*pt*pt);
              double b=-K*Pz/(Ea1*Ea1);
              double c=pt*pt-K*K/(4.0*Ea1*Ea1);
              R=b*b-4.0*a*c;
	      std::cout << "R: " << R << std::endl;
              if(R<0){if(verbose)std::cout << "R is negative " << R << "setting R to 0." << std::endl;R=0;}
              double Pnuz1=(-b-sqrt(R))/(2.0*a);
              double Pnuz2=(-b+sqrt(R))/(2.0*a);
              Enu1=sqrt(Pnuz1*Pnuz1+pt*pt);
              Enu2=sqrt(Pnuz2*Pnuz2+pt*pt);
	      std::cout << " Ea1 " << Ea1 << " Pz " << Pz << " ma1 " << ma1 << " pt " << pt << " K " << K << " mtau " << mtau << std::endl;
	      std::cout << "a " << a << " b " << b << " c " << c << std::endl;
	      std::cout << "R: " << R << " -b/2a " << (-b)/(2*a) << std::endl;
	      ///////////////////////////////////////////////////////////////////////////////////////////////////////
              TLorentzVector Neutrino1(Pnuz1*ntau.Px()-Pt_vec.Px(),Pnuz1*ntau.Py()-Pt_vec.Py(),Pnuz1*ntau.Pz()-Pt_vec.Pz(),Enu1);
              TauSolution1=TauVis1+Neutrino1;
	      TLorentzVector Neutrino2(Pnuz2*ntau.Px()-Pt_vec.Px(),Pnuz2*ntau.Py()-Pt_vec.Py(),Pnuz2*ntau.Pz()-Pt_vec.Pz(),Enu2);
              TauSolution1=TauVis2+Neutrino2;
            }
            // Method 4: Solve for neutrino E without rotating frames
            if(method==2){
              TVector3 ntau=TauDirection;
              if(ntau.Mag()!=0)ntau*=1/ntau.Mag();
              Pz=ntau.Dot(TauVis1.Vect());
              TVector3 Pz_vec=ntau;
              Pz_vec*=Pz;
              TVector3 Pt_vec=TauVis1.Vect();
              Pt_vec-=Pz_vec;
              pt=Pt_vec.Mag();
              ///////////////////////////////////////////////////////////////////////////////////////////////////////
              double a=1-Ea1*Ea1/(Pz*Pz);
              double K=(mtau*mtau-ma1*ma1-2*pt*pt);
              double b=K*Ea1/(Pz*Pz);
              double c=-(pt*pt+K*K/(4*Pz*Pz));
              R=b*b-4*a*c;
	      std::cout << "R: " << R << std::endl;
              if(R<0){if(verbose)std::cout << "R is negative " << R << "setting R to 0." << std::endl;R=0;}
              Enu1=(-b-sqrt(R))/(2.0*a);
              Enu2=(-b+sqrt(R))/(2.0*a);
	      double Pnuz1=0;if(Enu1>pt)Pnuz1=sqrt(Enu1*Enu1-pt*pt);
	      double Pnuz2=0;if(Enu2>pt)Pnuz2=sqrt(Enu2*Enu2-pt*pt);
	      std::cout << " Ea1 " << Ea1 << " Pz " << Pz << " ma1 " << ma1 << " pt " << pt << " K " << K << " mtau " << mtau << std::endl;
	      std::cout << "a " << a << " b " << b << " c " << c << std::endl;
	      std::cout << "R: " << R << " -b/2a " << (-b)/(2*a) << std::endl;
              ///////////////////////////////////////////////////////////////////////////////////////////////////////
              TLorentzVector Neutrino1(Pnuz1*ntau.Px()-Pt_vec.Px(),Pnuz1*ntau.Py()-Pt_vec.Py(),Pnuz1*ntau.Pz()-Pt_vec.Pz(),Enu1);
              TauSolution1=TauVis1+Neutrino1;
              TLorentzVector Neutrino2(Pnuz2*ntau.Px()-Pt_vec.Px(),Pnuz2*ntau.Py()-Pt_vec.Py(),Pnuz2*ntau.Pz()-Pt_vec.Pz(),Enu2);
              TauSolution2=TauVis2+Neutrino2;
            }



	    // Summary output
	    /*	    std::cout << "Sol1 4-vector: " 
		      << -TauVisInput.Px() << " " << -TauVisInput.Py() << " " << sqrt(Enu1*Enu1-pt*pt) << " " << Enu1 <<std::endl;	    
	    std::cout << "Sol1 lab 4-vector: "
                      << Neutrino1.Px() << " " << Neutrino1.Py() << " " << Neutrino1.Pz() << " " << Neutrino1.E() <<std::endl;
	    std::cout << "Sol2 4-vector: "
                      << -TauVisInput.Px() << " " << -TauVisInput.Py() << " " << sqrt(Enu2*Enu2-pt*pt) << " " << Enu2 <<std::endl;
	    std::cout << "Sol2 lab 4-vector: "
                      << Neutrino2.Px() << " " << Neutrino2.Py() << " " << Neutrino2.Pz() << " " << Neutrino2.E() <<std::endl;
	    */
	    if(fabs(MCTau_LV.E()-TauSolution1.E())<fabs(MCTau_LV.E()-TauSolution2.E()) && TauSolution1.E()>0){
	      Enu=Enu1;
	      TauSolution=TauSolution1;
	      TauSolutionResult.at(t).Fill(-1.0,w);
	      break;
	    }
	    else{
	      Enu=Enu2;
	      TauSolution=TauSolution2;
	      TauSolutionResult.at(t).Fill(1.0,w);
	      break;
	    }

	}
	  if(verbose)std::cout << "Tau 4-vector: Px=" << Ntp->KFTau_TauFit_p4(tau_idx).Px()
			       << " Py=" << Ntp->KFTau_TauFit_p4(tau_idx).Py()
			       << " Pz=" << Ntp->KFTau_TauFit_p4(tau_idx).Pz()
			       << " M=" << Ntp->KFTau_TauFit_p4(tau_idx).M() <<std::endl;         
       
          EstimatedTauE.at(t).Fill(TauSolution.E(),w);
          EstimatedTauPhi.at(t).Fill(TauSolution.Phi(),w);
          EstimatedTauEta.at(t).Fill(TauSolution.Eta(),w);
	  
	  if(verbose)std::cout << "Tau direction: phi= " << TauDirection.Phi() << " theta=" << TauDirection.Theta() << std::endl;
	  if(verbose)std::cout << "MC Tau direction: phi= " << MCTau_LV.Phi() << " theta=" << MCTau_LV.Theta() << std::endl; 
	  if(verbose)std::cout << "MC Tau 4-vector: Px=" << MCTau_LV.Px() << " Py=" << MCTau_LV.Py() <<" Pz=" << MCTau_LV.Pz() 
			       <<" Pm=" << MCTau_LV.M() << std::endl; 
	  if(verbose)std::cout << "Tau KF Resolution dphi=" 
			       << Tools::DeltaPhi(Ntp->KFTau_TauFit_p4(tau_idx).Phi(),MCTau_LV.Phi()) 
			       << " deta" << Ntp->KFTau_TauFit_p4(tau_idx).Eta()-MCTau_LV.Eta() 
			       << " dE=" << Ntp->KFTau_TauFit_p4(tau_idx).E()-MCTau_LV.E() << std::endl;
	  //	  if(verbose)
          std::cout << "Tau Solution Energy  Resolution Solution dE="
		    << TauSolution.E()-MCTau_LV.E()
		    << " Enu1 dE=" << TauSolution1.E()-MCTau_LV.E() 
		    << " Enu2 dE=" << TauSolution2.E()-MCTau_LV.E() << std::endl;
	  //      if(verbose)
	  std::cout << "Tau Estimate vs direction check dphi="
				<< Tools::DeltaPhi(TauSolution.Phi(),TauDirection.Phi())
				<< " deta" << TauSolution.Eta()- TauDirection.Eta() << std::endl;

	  TauCandERes.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).E()-MCTau_LV.E(),w);
          TauCandPhiRes.at(t).Fill(Tools::DeltaPhi(Ntp->KFTau_TauFit_p4(tau_idx).Phi(),MCTau_LV.Phi()),w);
          TauCandEtaRes.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).Eta()-MCTau_LV.Eta(),w);
	  
          EstimatedTauERes.at(t).Fill(TauSolution.E()-MCTau_LV.E(),w);
          EstimatedTauPhiRes.at(t).Fill(Tools::DeltaPhi(TauSolution.Phi(),MCTau_LV.Phi()),w);
          EstimatedTauEtaRes.at(t).Fill(TauSolution.Eta()-MCTau_LV.Eta(),w);

	  EstimatedTauDirPhiRes.at(t).Fill(Tools::DeltaPhi(TauDirection.Phi(),MCTau_LV.Phi()),w);
	  EstimatedTauDirEtaRes.at(t).Fill(TauDirection.Eta()-MCTau_LV.Eta(),w);
        }
      }
    }
  }
}
      
