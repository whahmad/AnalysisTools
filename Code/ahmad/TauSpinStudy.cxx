#include "TauSpinStudy.h"  
#include "TLorentzVector.h"   
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TF1.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSpinerInterface.h"
//#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

                 
TauSpinStudy::TauSpinStudy(TString Name_, TString id_):
  Selection(Name_,id_)
{   
  //verbose=true;
}

TauSpinStudy::~TauSpinStudy() {
  for(unsigned int j=0; j<Npassed.size(); j++) {
    std::cout << "TauSpinStudy::~TauSpinStudy Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TauSpinStudy::~TauSpinStudy()" << std::endl;
}
      
void  TauSpinStudy::Configure() {
  // Setup Cut Values
  for(int i=0; i<NCuts;i++) {
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==isZtautautopimu)          cut.at(isZtautautopimu)=1;
  }

  TString hlabel;
  TString htitle; for(int i=0; i<NCuts; i++) {
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==isZtautautopimu) {
      title.at(i)="Is $Z\\rightarrow\\tau\\tau\\rightarrow\\mu\\pi$ MC (bool)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Is Z#rightarrow#tau#tau#rightarrow#mu#pi MC (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isZtautautopimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isZtautautopimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  //mu
  mu_ExoverEtau=HConfig.GetTH1D(Name+"_mu_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{#mu}/E_{#tau}","Events");
  mu_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{+}}","Events");
  mu_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{-}}","Events");
  mu_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_mu_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{Spin}","Events");
  mu_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_mu_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{UnSpin}","Events");
  mu_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_mu_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{FlipSpin}","Events");
  mu_WT_Spin=HConfig.GetTH1D(Name+"_mu_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#mu}","Events");
  mu_WT_UnSpin=HConfig.GetTH1D(Name+"_mu_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#mu}","Events");
  mu_WT_FlipSpin=HConfig.GetTH1D(Name+"_mu_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#mu}","Events");
  mu_PtRatio_hplus=HConfig.GetTH1D(Name+"_mu_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#mu}/P_{t,#tau}|_{h^{+}}","Events");
  mu_PtRatio_hminus=HConfig.GetTH1D(Name+"_mu_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#mu}/P_{t,#tau}|_{h^{-}}","Events");

  //murec
  murec_ExoverEtaurec=HConfig.GetTH1D(Name+"_murec_ExoverEtaurec","ExoverEtaurec",20,0.0,1.0,"E_{#mu_{rec}}/E_{#tau_{rec}}","Events");
  murec_ExoverEtaurec_hplus=HConfig.GetTH1D(Name+"_murec_ExoverEtaurec_hplus","ExoverEtaurec_hplus",20,0.0,1.0,"E_{#mu_{rec}}/E_{#tau_{rec}}|_{h^{+}}","Events");
  murec_ExoverEtaurec_hminus=HConfig.GetTH1D(Name+"_murec_ExoverEtaurec_hminus","ExoverEtaurec_hminus",20,0.0,1.0,"E_{#mu_{rec}}/E_{#tau_{rec}}|_{h^{-}}","Events");
  murec_ExoverEtaurec_Spin=HConfig.GetTH1D(Name+"_murec_ExoverEtaurec_Spin","ExoverEtaurec_Spin",20,0.0,1.0,"E_{#mu_{rec}}/E_{#tau_{rec}}|_{Spin}","Events");
  murec_ExoverEtaurec_UnSpin=HConfig.GetTH1D(Name+"_murec_ExoverEtaurec_UnSpin","ExoverEtaurec_UnSpin",20,0.0,1.0,"E_{#mu_{rec}}/E_{#tau_{rec}}|_{UnSpin}","Events");
  murec_ExoverEtaurec_FlipSpin=HConfig.GetTH1D(Name+"_murec_ExoverEtaurec_FlipSpin","ExoverEtaurec_FlipSpin",20,0.0,1.0,"E_{#mu_{rec}}/E_{#tau_{rec}}|_{FlipSpin}","Events");
  murec_WT_Spin=HConfig.GetTH1D(Name+"_murec_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#mu_{rec}}","Events");
  murec_WT_UnSpin=HConfig.GetTH1D(Name+"_murec_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#mu_{rec}}","Events");
  murec_WT_FlipSpin=HConfig.GetTH1D(Name+"_murec_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#mu_{rec}}","Events");
  murec_PtRatio_hplus=HConfig.GetTH1D(Name+"_murec_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#mu_{rec}}/P_{t,#tau_{rec}}|_{h^{+}}","Events");
  murec_PtRatio_hminus=HConfig.GetTH1D(Name+"_murec_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#mu_{rec}}/P_{t,#tau_{rec}}|_{h^{-}}","Events");
  //pinu
  pi_ExoverEtau=HConfig.GetTH1D(Name+"_pi_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{#pi}/E_{#tau}","Events");
  pi_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}","Events");
  pi_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}","Events");
  pi_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_pi_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{Spin}","Events");
  pi_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_pi_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{UnSpin}","Events");
  pi_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_pi_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{FlipSpin}","Events");
  pi_WT_Spin=HConfig.GetTH1D(Name+"_pi_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#pi}","Events");
  pi_WT_UnSpin=HConfig.GetTH1D(Name+"_pi_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#pi}","Events");
  pi_WT_FlipSpin=HConfig.GetTH1D(Name+"_pi_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#pi}","Events");
  pi_PtRatio_hplus=HConfig.GetTH1D(Name+"_pi_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{+}}","Events");
  pi_PtRatio_hminus=HConfig.GetTH1D(Name+"_pi_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#pi}/P_t,{#tau}|_{h^{-}}","Events");

  pi_Mvis=HConfig.GetTH1D(Name+"_pi_Mvis","Mvis",20,0.0,1.0,"M_{Vis[#pi#pi]}}","Events");
  pi_Mvis_hplus=HConfig.GetTH1D(Name+"_pi_Mvis_hplus","Mvis_hplus",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{h^{+}}","Events");
  pi_Mvis_hminus=HConfig.GetTH1D(Name+"_pi_Mvis_hminus","Mvis_hminus",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{h^{-}}","Events");
  pi_Mvis_Spin=HConfig.GetTH1D(Name+"_pi_Mvis_Spin","Mvis_Spin",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{Spin}","Events");
  pi_Mvis_UnSpin=HConfig.GetTH1D(Name+"_pi_Mvis_UnSpin","Mvis_UnSpin",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{UnSpin}","Events");
  pi_Mvis_FlipSpin=HConfig.GetTH1D(Name+"_pi_Mvis_FlipSpin","Mvis_FlipSpin",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{FlipSpin}","Events");

 //pinu rec
  pirec_ExoverEtaurec=HConfig.GetTH1D(Name+"_pirec_ExoverEtaurec","ExoverEtaurec",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau}","Events");
  pirec_ExoverEtaurec_hplus=HConfig.GetTH1D(Name+"_pirec_ExoverEtaurec_hplus","ExoverEtaurec_hplus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau}|_{h^{+}}","Events");
  pirec_ExoverEtaurec_hminus=HConfig.GetTH1D(Name+"_pirec_ExoverEtaurec_hminus","ExoverEtaurec_hminus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau}|_{h^{-}}","Events");
  pirec_ExoverEtaurec_Spin=HConfig.GetTH1D(Name+"_pirec_ExoverEtaurec_Spin","ExoverEtaurec_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau}|_{Spin}","Events");
  pirec_ExoverEtaurec_UnSpin=HConfig.GetTH1D(Name+"_pirec_ExoverEtaurec_UnSpin","ExoverEtaurec_UnSpin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau}|_{UnSpin}","Events");
  pirec_ExoverEtaurec_FlipSpin=HConfig.GetTH1D(Name+"_pirec_ExoverEtaurec_FlipSpin","ExoverEtaurec_FlipSpin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau}|_{FlipSpin}","Events");
  pirec_WT_Spin=HConfig.GetTH1D(Name+"_pirec_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#pi_{rec}}","Events");
  pirec_WT_UnSpin=HConfig.GetTH1D(Name+"_pirec_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#pi_{rec}}","Events");
  pirec_WT_FlipSpin=HConfig.GetTH1D(Name+"_pirec_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#pi_{rec}}","Events");
  pirec_PtRatio_hplus=HConfig.GetTH1D(Name+"_pirec_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{+}}","Events");
  pirec_PtRatio_hminus=HConfig.GetTH1D(Name+"_pirec_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#pi_{rec}}/P_t,{#tau_{rec}}|_{h^{-}}","Events");
  

  // Setup Extra Histograms
  // Reco Objects
  Pion_pt=HConfig.GetTH1D(Name+"_Pion_pt","Pion_pt",50,0,100,"Pion p_{t}","Events");
  Pion_phi=HConfig.GetTH1D(Name+"_Pion_phi","Pion_phi",30,-TMath::Pi(),TMath::Pi(),"Pion #phi","Events");
  Pion_eta=HConfig.GetTH1D(Name+"_Pion_eta","Pion_eta",50,-2.5,2.5,"Pion #eta","Events");
  Pion_theta=HConfig.GetTH1D(Name+"_Pion_theta","Pion_theta",50,-2.5,2.5,"Pion #theta","Events");
  Mu_pt=HConfig.GetTH1D(Name+"_Mu_pt","Mu_pt",50,0,100,"Muon p_{t}","Events");
  Mu_phi=HConfig.GetTH1D(Name+"_Mu_phi","Mu_phi",30,-TMath::Pi(),TMath::Pi(),"Muon #phi","Events");
  Mu_eta=HConfig.GetTH1D(Name+"_Mu_eta","Mu_eta",50,-2.5,2.5,"Muon #eta","Events");
  Mu_theta=HConfig.GetTH1D(Name+"_Mu_theta","Mu_theta",50,-2.5,2.5,"Muon #theta","Events");
  MET_phi=HConfig.GetTH1D(Name+"_MET_phi","MET_phi",30,-TMath::Pi(),TMath::Pi(),"MET #phi","Events");
  MET_eta=HConfig.GetTH1D(Name+"_MET_eta","MET_eta",50,-2.5,2.5,"MET #eta","Events");
  // Generator Objects
  //Z particle
  M_Z_gen=HConfig.GetTH1D(Name+"_M_Z_gen","M_Z_gen",180,50.,160.,"M_Z_gen","Events");
  Pt_Z_gen=HConfig.GetTH1D(Name+"_Pt_Z_gen","Pt_Z_gen",30,0.,80.,"Pt_Z_gen","Events");
  Eta_Z_gen=HConfig.GetTH1D(Name+"_Eta_Z_gen","Eta_Z_gen",50,-2.5,2.5,"Eta_Z_gen","Events");
  Phi_Z_gen=HConfig.GetTH1D(Name+"_Phi_Z_gen","Phi_Z_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Z_gen","Events");
  // Z1 particle
  M_Z1_gen=HConfig.GetTH1D(Name+"_M_Z1_gen","M_Z1_gen",180,50.,160.,"M_Z1_gen","Events");
  Pt_Z1_gen=HConfig.GetTH1D(Name+"_Pt_Z1_gen","Pt_Z1_gen",30,0.,80.,"Pt_Z1_gen","Events");
  Eta_Z1_gen=HConfig.GetTH1D(Name+"_Eta_Z1_gen","Eta_Z1_gen",50,-2.5,2.5,"Eta_Z1_gen","Events");
  Phi_Z1_gen=HConfig.GetTH1D(Name+"_Phi_Z1_gen","Phi_Z1_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Z1_gen","Events");
  //Z---> TauPi
  M_TauPi_gen=HConfig.GetTH1D(Name+"_M_TauPi_gen","M_TauPi_gen",180,0.,4.,"M_TauPi_gen","Events");
  Pt_TauPi_gen=HConfig.GetTH1D(Name+"_Pt_TauPi_gen","Pt_TauPi_gen",30,0.,70.,"Pt_TauPi_gen","Events");
  Eta_TauPi_gen=HConfig.GetTH1D(Name+"_Eta_TauPi_gen","Eta_TauPi_gen",50,-2.5,2.5,"Eta_TauPi_gen","Events");
  Phi_TauPi_gen=HConfig.GetTH1D(Name+"_Phi_TauPi_gen","Phi_TauPi_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_gen","Events");
  //Z---> TauMu
  M_TauMu_gen=HConfig.GetTH1D(Name+"_M_TauMu_gen","M_TauMu_gen",180,0.,4.,"M_TauMu_gen","Events");
  Pt_TauMu_gen=HConfig.GetTH1D(Name+"_Pt_TauMu_gen","Pt_TauMu_gen",30,0.,70.,"Pt_TauMu_gen","Events");
  Eta_TauMu_gen=HConfig.GetTH1D(Name+"_Eta_TauMu_gen","Eta_TauMu_gen",50,-2.5,2.5,"Eta_TauMu_gen","Events");
  Phi_TauMu_gen=HConfig.GetTH1D(Name+"_Phi_TauMu_gen","Phi_TauMu_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_gen","Events");
  //TauPi--->Pi
  M_Pi_gen=HConfig.GetTH1D(Name+"_M_Pi_gen","M_Pi_gen",180,0.,0.3,"M_Pi_gen","Events");
  Pt_Pi_gen=HConfig.GetTH1D(Name+"_Pt_Pi_gen","Pt_Pi_gen",30,0.,70.,"Pt_Pi_gen","Events");
  Eta_Pi_gen=HConfig.GetTH1D(Name+"_Eta_Pi_gen","Eta_Pi_gen",50,-2.5,2.5,"Eta_Pi_gen","Events");
  Phi_Pi_gen=HConfig.GetTH1D(Name+"_Phi_Pi_gen","Phi_Pi_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Pi_gen","Events");
  //TauMu--->Mu
  M_Mu_gen=HConfig.GetTH1D(Name+"_M_Mu_gen","M_Mu_gen",180,0.,0.21,"M_Mu_gen","Events");
  Pt_Mu_gen=HConfig.GetTH1D(Name+"_Pt_Mu_gen","Pt_Mu_gen",30,0.,70.,"Pt_Mu_gen","Events");
  Eta_Mu_gen=HConfig.GetTH1D(Name+"_Eta_Mu_gen","Eta_Mu_gen",50,-2.5,2.5,"Eta_Mu_gen","Events");
  Phi_Mu_gen=HConfig.GetTH1D(Name+"_Phi_Mu_gen","Phi_Mu_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Mu_gen","Events");
  //TauPi---> Nutp
  M_Nutp_gen=HConfig.GetTH1D(Name+"_M_Nutp_gen","M_Nutp_gen",180,0.,0.000000005,"M_Nutp_gen","Events");
  Pt_Nutp_gen=HConfig.GetTH1D(Name+"_Pt_Nutp_gen","Pt_Nutp_gen",30,0.,70.,"Pt_Nutp_gen","Events");
  Eta_Nutp_gen=HConfig.GetTH1D(Name+"_Eta_Nutp_gen","Eta_Nutp_gen",50,-2.5,2.5,"Eta_Nutp_gen","Events");
  Phi_Nutp_gen=HConfig.GetTH1D(Name+"_Phi_Nutp_gen","Phi_Nutp_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Nutp_gen","Events");
  //TauMu-> Nutm
  M_Nutm_gen=HConfig.GetTH1D(Name+"_M_Nutm_gen","M_Nutm_gen",180,0.,0.000000005,"M_Nutm_gen","Events");
  Pt_Nutm_gen=HConfig.GetTH1D(Name+"_Pt_Nutm_gen","Pt_Nutm_gen",30,0.,70.,"Pt_Nutm_gen","Events");
  Eta_Nutm_gen=HConfig.GetTH1D(Name+"_Eta_Nutm_gen","Eta_Nutm_gen",50,-2.5,2.5,"Eta_Nutm_gen","Events");
  Phi_Nutm_gen=HConfig.GetTH1D(Name+"_Phi_Nutm_gen","Phi_Nutm_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Nutm_gen","Events");
  //TauMu-> Num
  M_Num_gen=HConfig.GetTH1D(Name+"_M_Num_gen","M_Num_gen",180,0.,0.000000005,"M_Num_gen","Events");
  Pt_Num_gen=HConfig.GetTH1D(Name+"_Pt_Num_gen","Pt_Num_gen",30,0.,70.,"Pt_Num_gen","Events");
  Eta_Num_gen=HConfig.GetTH1D(Name+"_Eta_Num_gen","Eta_Num_gen",50,-2.5,2.5,"Eta_Num_gen","Events");
  Phi_Num_gen=HConfig.GetTH1D(Name+"_Phi_Num_gen","Phi_Num_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_Num_gen","Events");
  //Delta_R(Mu_rec,Mu_gen), Delta_R(Pi_rec,Pi_gen)
  dR_recMu_genMu=HConfig.GetTH1D(Name+"_dR_recMu_genMu","dR_recMu_genMu",100,0.,8.,"dR_recMu_genMu","Events");
  mindR_recMu_genMu=HConfig.GetTH1D(Name+"_mindR_recMu_genMu","mindR_recMu_genMu",100,0.,1.,"mindR_recMu_genMu","Events");
  dR_recPi_genPi=HConfig.GetTH1D(Name+"_dR_recPi_genPi","dR_recPi_genPi",100,0.,8.,"dR_recPi_genPi","Events");
  mindR_recPi_genPi=HConfig.GetTH1D(Name+"_mindR_recPi_genPi","mindR_recPi_genPi",100,0.,1.,"mindR_recPi_genPi","Events");

  // reco midR Mu Pi
  PionRec_pt=HConfig.GetTH1D(Name+"_PionRec_pt","PionRec_pt",50,0,100,"PionRec p_{t}","Events");
  PionRec_phi=HConfig.GetTH1D(Name+"_PionRec_phi","PionRec_phi",30,-TMath::Pi(),TMath::Pi(),"PionRec #phi","Events");
  PionRec_eta=HConfig.GetTH1D(Name+"_PionRec_eta","PionRec_eta",50,-2.5,2.5,"PionRec #eta","Events");
  Pion_mdR_pt=HConfig.GetTH1D(Name+"_Pion_mdR_pt","Pion_mdR_pt",50,0,100,"Pion_mdR p_{t}","Events");
  Pion_mdR_phi=HConfig.GetTH1D(Name+"_Pion_mdR_phi","Pion_mdR_phi",30,-TMath::Pi(),TMath::Pi(),"Pion_mdR #phi","Events");
  Pion_mdR_eta=HConfig.GetTH1D(Name+"_Pion_mdR_eta","Pion_mdR_eta",50,-2.5,2.5,"Pion_mdR #eta","Events");
  
  MuonRec_pt=HConfig.GetTH1D(Name+"_MuonRec_pt","MuonRec_pt",50,0,100,"MuonRec p_{t}","Events");
  MuonRec_phi=HConfig.GetTH1D(Name+"_MuonRec_phi","MuonRec_phi",30,-TMath::Pi(),TMath::Pi(),"MuonRec #phi","Events");
  MuonRec_eta=HConfig.GetTH1D(Name+"_MuonRec_eta","MuonRec_eta",50,-2.5,2.5,"MuonRec #eta","Events");
  Muon_mdR_pt=HConfig.GetTH1D(Name+"_Muon_mdR_pt","Muon_mdR_pt",50,0,100,"Muon_mdR p_{t}","Events");
  Muon_mdR_phi=HConfig.GetTH1D(Name+"_Muon_mdR_phi","Muon_mdR_phi",30,-TMath::Pi(),TMath::Pi(),"Muon_mdR #phi","Events");
  Muon_mdR_eta=HConfig.GetTH1D(Name+"_Muon_mdR_eta","Muon_mdR_eta",50,-2.5,2.5,"Muon_mdR #eta","Events");
  //Z---> TauPirec
  M_TauPi_rec=HConfig.GetTH1D(Name+"_M_TauPi_rec","M_TauPi_rec",180,0.,4.,"M_TauPi_rec","Events");
  Pt_TauPi_rec=HConfig.GetTH1D(Name+"_Pt_TauPi_rec","Pt_TauPi_rec",30,0.,70.,"Pt_TauPi_rec","Events");
  Eta_TauPi_rec=HConfig.GetTH1D(Name+"_Eta_TauPi_rec","Eta_TauPi_rec",50,-2.5,2.5,"Eta_TauPi_rec","Events");
  Phi_TauPi_rec=HConfig.GetTH1D(Name+"_Phi_TauPi_rec","Phi_TauPi_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_rec","Events");
  //Z---> TauMurec
  M_TauMu_rec=HConfig.GetTH1D(Name+"_M_TauMu_rec","M_TauMu_rec",180,0.,4.,"M_TauMu_rec","Events");
  Pt_TauMu_rec=HConfig.GetTH1D(Name+"_Pt_TauMu_rec","Pt_TauMu_rec",30,0.,70.,"Pt_TauMu_rec","Events");
  Eta_TauMu_rec=HConfig.GetTH1D(Name+"_Eta_TauMu_rec","Eta_TauMu_rec",50,-2.5,2.5,"Eta_TauMu_rec","Events");
  Phi_TauMu_rec=HConfig.GetTH1D(Name+"_Phi_TauMu_rec","Phi_TauMu_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_rec","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  TauSpinStudy::Store_ExtraDist() {
  
  Extradist1d.push_back(&mu_ExoverEtau);
  Extradist1d.push_back(&mu_ExoverEtau_hplus);
  Extradist1d.push_back(&mu_ExoverEtau_hminus);
  Extradist1d.push_back(&mu_ExoverEtau_Spin);
  Extradist1d.push_back(&mu_ExoverEtau_UnSpin);
  Extradist1d.push_back(&mu_ExoverEtau_FlipSpin);
  Extradist1d.push_back(&mu_WT_Spin);
  Extradist1d.push_back(&mu_WT_UnSpin);
  Extradist1d.push_back(&mu_WT_FlipSpin);
  Extradist1d.push_back(&mu_PtRatio_hplus);
  Extradist1d.push_back(&mu_PtRatio_hminus);

  Extradist1d.push_back(&murec_ExoverEtaurec);
  Extradist1d.push_back(&murec_ExoverEtaurec_hplus);
  Extradist1d.push_back(&murec_ExoverEtaurec_hminus);
  Extradist1d.push_back(&murec_ExoverEtaurec_Spin);
  Extradist1d.push_back(&murec_ExoverEtaurec_UnSpin);
  Extradist1d.push_back(&murec_ExoverEtaurec_FlipSpin);
  Extradist1d.push_back(&murec_WT_Spin);
  Extradist1d.push_back(&murec_WT_UnSpin);
  Extradist1d.push_back(&murec_WT_FlipSpin);
  Extradist1d.push_back(&murec_PtRatio_hplus);
  Extradist1d.push_back(&murec_PtRatio_hminus);

  Extradist1d.push_back(&pi_ExoverEtau);
  Extradist1d.push_back(&pi_ExoverEtau_hplus);
  Extradist1d.push_back(&pi_ExoverEtau_hminus);
  Extradist1d.push_back(&pi_ExoverEtau_Spin);
  Extradist1d.push_back(&pi_ExoverEtau_UnSpin);
  Extradist1d.push_back(&pi_ExoverEtau_FlipSpin);
  Extradist1d.push_back(&pi_WT_Spin);
  Extradist1d.push_back(&pi_WT_UnSpin);
  Extradist1d.push_back(&pi_WT_FlipSpin);
  Extradist1d.push_back(&pi_PtRatio_hplus);
  Extradist1d.push_back(&pi_PtRatio_hminus);
  
  Extradist1d.push_back(&pirec_ExoverEtaurec);
  Extradist1d.push_back(&pirec_ExoverEtaurec_hplus);
  Extradist1d.push_back(&pirec_ExoverEtaurec_hminus);
  Extradist1d.push_back(&pirec_ExoverEtaurec_Spin);
  Extradist1d.push_back(&pirec_ExoverEtaurec_UnSpin);
  Extradist1d.push_back(&pirec_ExoverEtaurec_FlipSpin);
  Extradist1d.push_back(&pirec_WT_Spin);
  Extradist1d.push_back(&pirec_WT_UnSpin);
  Extradist1d.push_back(&pirec_WT_FlipSpin);
  Extradist1d.push_back(&pirec_PtRatio_hplus);
  Extradist1d.push_back(&pirec_PtRatio_hminus);

  Extradist1d.push_back(&pi_Mvis);
  Extradist1d.push_back(&pi_Mvis_hplus);
  Extradist1d.push_back(&pi_Mvis_hminus);
  Extradist1d.push_back(&pi_Mvis_Spin);
  Extradist1d.push_back(&pi_Mvis_UnSpin);
  Extradist1d.push_back(&pi_Mvis_FlipSpin);

  //Reco objects
  Extradist1d.push_back(&Pion_pt);            Extradist1d.push_back(&Pion_eta);
  Extradist1d.push_back(&Pion_phi);           Extradist1d.push_back(&Mu_pt);
  Extradist1d.push_back(&Mu_eta);             Extradist1d.push_back(&Mu_phi);           
  Extradist1d.push_back(&MET_phi);            Extradist1d.push_back(&MET_eta);  

  Extradist1d.push_back(&PionRec_pt);         Extradist1d.push_back(&PionRec_eta);
  Extradist1d.push_back(&PionRec_phi);        Extradist1d.push_back(&Pion_mdR_pt);
  Extradist1d.push_back(&Pion_mdR_eta);       Extradist1d.push_back(&Pion_mdR_phi); 
  
  Extradist1d.push_back(&MuonRec_pt);         Extradist1d.push_back(&MuonRec_eta);
  Extradist1d.push_back(&MuonRec_phi);        Extradist1d.push_back(&Muon_mdR_pt);
  Extradist1d.push_back(&Muon_mdR_eta);       Extradist1d.push_back(&Muon_mdR_phi);

  Extradist1d.push_back(&M_TauPi_rec);        Extradist1d.push_back(&Pt_TauPi_rec);
  Extradist1d.push_back(&Eta_TauPi_rec);      Extradist1d.push_back(&Phi_TauPi_rec);
  Extradist1d.push_back(&M_TauMu_rec);        Extradist1d.push_back(&Pt_TauMu_rec);
  Extradist1d.push_back(&Eta_TauMu_rec);      Extradist1d.push_back(&Phi_TauMu_rec);
        
  //Generator objets   
  Extradist1d.push_back(&M_Z_gen);            Extradist1d.push_back(&Pt_Z_gen);
  Extradist1d.push_back(&Eta_Z_gen);          Extradist1d.push_back(&Phi_Z_gen);
  Extradist1d.push_back(&M_Z1_gen);           Extradist1d.push_back(&Pt_Z1_gen);
  Extradist1d.push_back(&Eta_Z1_gen);         Extradist1d.push_back(&Phi_Z1_gen);
  Extradist1d.push_back(&M_TauPi_gen);        Extradist1d.push_back(&Pt_TauPi_gen);
  Extradist1d.push_back(&Eta_TauPi_gen);      Extradist1d.push_back(&Phi_TauPi_gen);
  Extradist1d.push_back(&M_TauMu_gen);        Extradist1d.push_back(&Pt_TauMu_gen);
  Extradist1d.push_back(&Eta_TauMu_gen);      Extradist1d.push_back(&Phi_TauMu_gen);
  Extradist1d.push_back(&M_Pi_gen);           Extradist1d.push_back(&Pt_Pi_gen);
  Extradist1d.push_back(&Eta_Pi_gen);         Extradist1d.push_back(&Phi_Pi_gen);   
  Extradist1d.push_back(&M_Mu_gen);           Extradist1d.push_back(&Pt_Mu_gen);
  Extradist1d.push_back(&Eta_Mu_gen);         Extradist1d.push_back(&Phi_Mu_gen);   
  Extradist1d.push_back(&M_Nutp_gen);         Extradist1d.push_back(&Pt_Nutp_gen);       
  Extradist1d.push_back(&Eta_Nutp_gen);       Extradist1d.push_back(&Phi_Nutp_gen); 
  Extradist1d.push_back(&M_Nutm_gen);         Extradist1d.push_back(&Pt_Nutm_gen);     
  Extradist1d.push_back(&Eta_Nutm_gen);       Extradist1d.push_back(&Phi_Nutm_gen); 
  Extradist1d.push_back(&M_Num_gen);          Extradist1d.push_back(&Pt_Num_gen);
  Extradist1d.push_back(&Eta_Num_gen);        Extradist1d.push_back(&Phi_Num_gen);    
  Extradist1d.push_back(&dR_recMu_genMu);     Extradist1d.push_back(&dR_recPi_genPi);
  Extradist1d.push_back(&mindR_recMu_genMu);  Extradist1d.push_back(&mindR_recPi_genPi);
  
}

void  TauSpinStudy::doEvent() {
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)) {t=2;}// std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,Boson3_idx,Boson4_idx,tau1_idx(0),tau2_idx(0),tau3_idx(0),tau4_idx(0);
  value.at(isZtautautopimu)=0;
  pass.at(isZtautautopimu) = true;
  if(pass.at(isZtautautopimu))value.at(isZtautautopimu)=1;
  
  double wobs=1;
  double w=1;

  //  if(verbose)  std::cout << "MCID=="<< Ntp->GetMCID() << "||Size==" << Npassed.size() << "||t== " << t << "||BosonIdx== " << Boson_idx << "||Tau1== " << tau1_idx << "||Tau2== " << tau1_idx << "||Ntaus== " << Ntp->NMCTaus() << std::endl;
  
  bool status=AnalysisCuts(t,w,wobs); 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status) {
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    // reject unsupported modes
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,tau1_idx,tau2_idx)) {
      bool ImpTau=true;
      double jakid1=Ntp->MCTau_JAK(tau1_idx);
      double jakid2=Ntp->MCTau_JAK(tau2_idx);
      if(!(jakid1==TauDecay::JAK_PION || jakid1==TauDecay::JAK_MUON)) ImpTau=false;
      if(!(jakid2==TauDecay::JAK_PION || jakid2==TauDecay::JAK_MUON)) ImpTau=false;
      //if(jakid1!=jakid2) return;
      if(!ImpTau) { std::cout << "Decay modes not implemented in TauSpinner " << std::endl; return;}
    }
    else {
      std::cout << "Not a valid Z0->tau+tau- decay" <<std::endl;
      return;
    }
    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    std::cout << "Spin_WT ==" << Spin_WT << std::endl;

    double UnSpin_WT=1;//Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
    double FlipSpin_WT=1;//Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
    double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    double hminus=1-hplus;//Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
    std::cout << "hplus " << hplus << " hminus " << hminus << std::endl;
    ////if(Spin_WT<=0 || 2<=Spin_WT) {Spin_WT=0.000001;UnSpin_WT=0.000001;FlipSpin_WT=0.000001;}
    if(verbose)std::cout<< "A " <<std::endl;
    ////////////////////////////////////////////////
    //
    // Spin Validation
    //  
    ////////////////////////////////////////////////
    
    ///////////////////////////////////////////    
    ////////////   Muon && Pion    ////////////
    ///////////////////////////////////////////
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson3_idx,TauDecay::JAK_MUON,tau3_idx) && Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson4_idx,TauDecay::JAK_PION,tau4_idx)) { //JAK_MUON && JAK_PION
      if(verbose)std::cout<< "pi-mu " <<std::endl;
      
      TLorentzVector GenZ=Ntp->MCSignalParticle_p4(Boson3_idx);
      TLorentzVector GenTauMu(0,0,0,0);
      TLorentzVector GenMu(0,0,0,0);  
      TLorentzVector GenNum(0,0,0,0);     
      TLorentzVector GenNutm(0,0,0,0);
      TLorentzVector Muon(0,0,0,0);
      TLorentzVector MuonRec(0,0,0,0);
      TLorentzVector Muon_mdR(0,0,0,0);
      
      TLorentzVector GenZ1=Ntp->MCSignalParticle_p4(Boson4_idx);
      TLorentzVector GenTauPi(0,0,0,0);
      TLorentzVector GenPi(0,0,0,0);
      TLorentzVector GenNutp(0,0,0,0);
      TLorentzVector Pion(0,0,0,0);
      TLorentzVector PionRec(0,0,0,0);
      TLorentzVector Pion_mdR(0,0,0,0);
      
      bool hasPirec=false;  bool hasPigen=false; 
      bool hasMurec=false;  bool hasMugen=false; 
      double dR_Muons=-1;   double dR_Pions=-1;
      double mindRMu=999;   double mindRPi=999; //dR_Muon < 0.1 && dR_Pion<0.4
      int mindRMu_Index=-1; int mindRPi_Index=-1;
      //unsigned int Nmuon=-1; unsigned int Npion=-1;
      
      //Gen Infos      
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau3_idx);i++) {
        if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::tau_minus)) {
          GenTauMu=Ntp->MCTauandProd_p4(tau3_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::mu_plus)) {
          GenMu+=Ntp->MCTauandProd_p4(tau3_idx,i);
	  hasMugen=true;
	}
	else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::nu_mu)) {
          GenNum+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::nu_tau)) {
          GenNutm+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
      } //end for loop 
       if(hasMugen) {
   	M_Z_gen.at(t).Fill(GenZ.M(),w);
	Pt_Z_gen.at(t).Fill(GenZ.Pt(),w);            
	Eta_Z_gen.at(t).Fill(GenZ.Eta(),w);
     	Phi_Z_gen.at(t).Fill(GenZ.Phi(),w);
	M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);
	Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);   
	Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);
     	Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w);
	M_Mu_gen.at(t).Fill(GenMu.M(),w); 
	Pt_Mu_gen.at(t).Fill(GenMu.Pt(),w);
	Eta_Mu_gen.at(t).Fill(GenMu.Eta(),w); 
     	Phi_Mu_gen.at(t).Fill(GenMu.Phi(),w); 
	M_Nutm_gen.at(t).Fill(GenNutm.M(),w); 
	Pt_Nutm_gen.at(t).Fill(GenNutm.Pt(),w);  
	Eta_Nutm_gen.at(t).Fill(GenNutm.Eta(),w); 
     	Phi_Nutm_gen.at(t).Fill(GenNutm.Phi(),w);
	M_Num_gen.at(t).Fill(GenNum.M(),w);
	Pt_Num_gen.at(t).Fill(GenNum.Pt(),w);    
	Eta_Num_gen.at(t).Fill(GenNum.Eta(),w); 
     	Phi_Num_gen.at(t).Fill(GenNum.Phi(),w);  
      }
      //Reco Muon
      for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++) {
	//std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
	Muon = Ntp->Muon_p4(iMuon);
	hasMurec = true;
      }
      if(Muon.Pt()>0) {
	Mu_pt.at(t).Fill(Muon.Pt(),w);
	Mu_phi.at(t).Fill(Muon.Phi(),w);
	Mu_eta.at(t).Fill(Muon.Eta(),w);
      }	
      for(unsigned int i=0; i < Ntp->NMuons(); i++) {
	//std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
	dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu) ;
	//std::cout<< "DR_Muon==" << dR_Muons <<std::endl;
	dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
	
	if(dR_Muons < mindRMu) { mindRMu = dR_Muons;  mindRMu_Index = i; }
	//std::cout<< "Moun_minDRindex==" << mindRMu_Index <<std::endl;
	MuonRec = Ntp->Muon_p4(mindRMu_Index);
      }
      if(MuonRec.Pt()>0 && mindRMu_Index>-1) {
	MuonRec_pt.at(t).Fill(MuonRec.Pt(),w);
	MuonRec_phi.at(t).Fill(MuonRec.Phi(),w);
	MuonRec_eta.at(t).Fill(MuonRec.Eta(),w);
      }
      if (mindRMu<0.1 && mindRMu_Index>-1) { 
	mindR_recMu_genMu.at(t).Fill(mindRMu,w);   
	std::cout<< "Moun_minDR==" << mindRMu <<std::endl;
	Muon_mdR=MuonRec;
	if(Muon_mdR.Pt()>0) {
	  Muon_mdR_pt.at(t).Fill(Muon_mdR.Pt(),w);
	  Muon_mdR_phi.at(t).Fill(Muon_mdR.Phi(),w);
	  Muon_mdR_eta.at(t).Fill(Muon_mdR.Eta(),w);
	}
      }
      if(GenTauMu.E()>0) {
        mu_PtRatio_hplus.at(t).Fill(GenMu.Pt()/GenTauMu.Pt(),w*hplus*UnSpin_WT);
        mu_PtRatio_hminus.at(t).Fill(GenMu.Pt()/GenTauMu.Pt(),w*hminus*UnSpin_WT);
      	GenTauMu.Boost(-GenZ.BoostVector());
        GenMu.Boost(-GenZ.BoostVector()); 
        GenZ.Boost(-GenZ.BoostVector());
        //Now fill results                                                                                        
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau3_idx));
      	mu_ExoverEtau.at(t).Fill(GenMu.E()/GenTauMu.E(),w);
      	mu_ExoverEtau_Spin.at(t).Fill(GenMu.E()/GenTauMu.E(),w*Spin_WT);
      	mu_ExoverEtau_UnSpin.at(t).Fill(GenMu.E()/GenTauMu.E(),w*UnSpin_WT);
      	mu_ExoverEtau_FlipSpin.at(t).Fill(GenMu.E()/GenTauMu.E(),w*FlipSpin_WT);
      	mu_ExoverEtau_hplus.at(t).Fill(GenMu.E()/GenTauMu.E(),w*hplus*UnSpin_WT);
      	mu_ExoverEtau_hminus.at(t).Fill(GenMu.E()/GenTauMu.E(),w*hminus*UnSpin_WT);
      	mu_WT_Spin.at(t).Fill(Spin_WT,w);
      	mu_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
      	mu_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
	
      	if(verbose) {
      	  for(int i=0;i<Ntp->NMCTauDecayProducts(tau3_idx);i++) {
      	    TLorentzVector LV=Ntp->MCTauandProd_p4(tau3_idx,i);
      	    //std::cout <<"MC_ID="<<Ntp->MCTauandProd_pdgid(tau3_idx,i)<< "||Px="<<LV.Px()<<"||Py="<<LV.Py()<<"||Pz="<<LV.Pz()<<"||E="<<LV.E()<< std::endl;
      	  }
      	}
      	//if(verbose)std::cout <<"mu:x="<<GenMu.E()/GenTauMu.E()<<"||w="<<w<<"||hplus="<<hplus<<"||hminus="<<hminus <<std::endl;
      } 
      
      //if(verbose)std::cout<< "muon done " <<std::endl;  
      //Gen Pion
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau4_idx);i++) { 
        if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::tau_minus)) {
          GenTauPi=Ntp->MCTauandProd_p4(tau4_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::pi_plus) ) {
          GenPi+=Ntp->MCTauandProd_p4(tau4_idx,i);
	  hasPigen = true;
	}
        else if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::nu_tau) ) {
	  GenNutp+=Ntp->MCTauandProd_p4(tau4_idx,i);
	}
      }      
      if(hasPigen) {
	M_Z1_gen.at(t).Fill(GenZ1.M(),w);
	Pt_Z1_gen.at(t).Fill(GenZ1.Pt(),w);            
	Eta_Z1_gen.at(t).Fill(GenZ1.Eta(),w);
     	Phi_Z1_gen.at(t).Fill(GenZ1.Phi(),w);

	M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);
	Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);
	Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w); 
     	Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w);
	
     	M_Pi_gen.at(t).Fill(GenPi.M(),w); 
	Pt_Pi_gen.at(t).Fill(GenPi.Pt(),w);
	Eta_Pi_gen.at(t).Fill(GenPi.Eta(),w); 
     	Phi_Pi_gen.at(t).Fill(GenPi.Phi(),w); 
	
     	M_Nutp_gen.at(t).Fill(GenNutp.M(),w);
        Pt_Nutp_gen.at(t).Fill(GenNutp.Pt(),w); 
	Eta_Nutp_gen.at(t).Fill(GenNutp.Eta(),w); 
     	Phi_Nutp_gen.at(t).Fill(GenNutp.Phi(),w); 
      }
      
      //Reco Pion
      for(unsigned int iPion=0; iPion < Ntp->NPFTaus(); iPion++) {
	//std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
	Pion = Ntp->PFTau_p4(iPion);
	hasPirec  = true;
      }
      if(Pion.Pt()>0) {
	Pion_pt.at(t).Fill(Pion.Pt(),w);  
	Pion_phi.at(t).Fill(Pion.Phi(),w);
	Pion_eta.at(t).Fill(Pion.Eta(),w);
      }
      for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
	//std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
	dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi) ;
	//std::cout<< "DR_Pion==" << dR_Pions <<std::endl;
	dR_recPi_genPi.at(t).Fill(dR_Pions,w);
	
	if(dR_Pions < mindRPi) { mindRPi = dR_Pions;  mindRPi_Index = i; }
	//std::cout<< "Pion_minDRindex==" << mindRPi_Index <<std::endl;
	PionRec = Ntp->PFTau_p4(mindRPi_Index);
      }
      if(PionRec.Pt()>0 && mindRPi_Index>-1) {
	PionRec_pt.at(t).Fill(PionRec.Pt(),w);
	PionRec_phi.at(t).Fill(PionRec.Phi(),w);
	PionRec_eta.at(t).Fill(PionRec.Eta(),w);
      }      
      if (mindRPi<0.4 && mindRPi_Index>-1) { 
	mindR_recPi_genPi.at(t).Fill(mindRPi,w);
	//std::cout<< "Pion_minDR==" << mindRPi <<std::endl;
	Pion_mdR=PionRec;
	if(Pion_mdR.Pt()>0) {
	  Pion_mdR_pt.at(t).Fill(Pion_mdR.Pt(),w);
	  Pion_mdR_phi.at(t).Fill(Pion_mdR.Phi(),w);
	  Pion_mdR_eta.at(t).Fill(Pion_mdR.Eta(),w);
	}      
      }
      if(GenTauPi.E()>0) {
      	pi_PtRatio_hplus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hplus*UnSpin_WT);
        pi_PtRatio_hminus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hminus*UnSpin_WT);
      	GenTauPi.Boost(-GenZ1.BoostVector());
        GenPi.Boost(-GenZ1.BoostVector());
        GenZ1.Boost(-GenZ1.BoostVector());
        // Now fill results                                   
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau4_idx));
      	pi_ExoverEtau.at(t).Fill(GenPi.E()/GenTauPi.E(),w);
      	pi_ExoverEtau_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*Spin_WT);
      	pi_ExoverEtau_UnSpin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*UnSpin_WT);
      	pi_ExoverEtau_FlipSpin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*FlipSpin_WT);
      	pi_ExoverEtau_hplus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hplus*UnSpin_WT);
      	pi_ExoverEtau_hminus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hminus*UnSpin_WT);
      	pi_WT_Spin.at(t).Fill(Spin_WT,w);
      	pi_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
      	pi_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
	
      	//if(verbose)std::cout << "pi : x " << GenPi.E()/GenTauPi.E() << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << std::endl;
      }
      //if(verbose)std::cout<< "pion done " <<std::endl; //if(verbose)std::cout<< "///////////////////// " <<std::endl;
      //if(verbose)std::cout<< "muon && pion done " <<std::endl;
      
      /////////////////////////////////
      //// Tau_rec  reconstruction ////
      /////////////////////////////////
      TLorentzVector TauMuRec(0,0,0,0);
      TLorentzVector TauPiRec(0,0,0,0);
      if(mindRMu<0.1 && mindRMu_Index>-1 && mindRPi<0.4 && mindRPi_Index>-1) {
	TauMuRec = Muon_mdR + GenNum + GenNutm;
	TauPiRec = Pion_mdR + GenNutp;
	if(TauMuRec.E()>0) {
	  M_TauMu_rec.at(t).Fill(TauMuRec.M(),w);
	  Pt_TauMu_rec.at(t).Fill(TauMuRec.Pt(),w);   
	  Eta_TauMu_rec.at(t).Fill(TauMuRec.Eta(),w);
	  Phi_TauMu_rec.at(t).Fill(TauMuRec.Phi(),w);
	}
	if(TauPiRec.E()>0) {	
	  M_TauPi_rec.at(t).Fill(TauPiRec.M(),w);
	  Pt_TauPi_rec.at(t).Fill(TauPiRec.Pt(),w);
	  Eta_TauPi_rec.at(t).Fill(TauPiRec.Eta(),w); 
	  Phi_TauPi_rec.at(t).Fill(TauPiRec.Phi(),w);
	}
	//}
      if(TauMuRec.E()>0) {
	murec_PtRatio_hplus.at(t).Fill(Muon_mdR.Pt()/TauMuRec.Pt(),w*hplus*UnSpin_WT);
	murec_PtRatio_hminus.at(t).Fill(Muon_mdR.Pt()/TauMuRec.Pt(),w*hminus*UnSpin_WT);
	
	TauMuRec.Boost(-GenZ.BoostVector());
	Muon_mdR.Boost(-GenZ.BoostVector());
	GenZ.Boost(-GenZ.BoostVector());
	// Now fill results                               
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau3_idx));
	murec_ExoverEtaurec.at(t).Fill(Muon_mdR.E()/TauMuRec.E(),w);
	murec_ExoverEtaurec_Spin.at(t).Fill(Muon_mdR.E()/TauMuRec.E(),w*Spin_WT);
	murec_ExoverEtaurec_UnSpin.at(t).Fill(Muon_mdR.E()/TauMuRec.E(),w*UnSpin_WT);
	murec_ExoverEtaurec_FlipSpin.at(t).Fill(Muon_mdR.E()/TauMuRec.E(),w*FlipSpin_WT);
	murec_ExoverEtaurec_hplus.at(t).Fill(Muon_mdR.E()/TauMuRec.E(),w*hplus*UnSpin_WT);
	murec_ExoverEtaurec_hminus.at(t).Fill(Muon_mdR.E()/TauMuRec.E(),w*hminus*UnSpin_WT);
	murec_WT_Spin.at(t).Fill(Spin_WT,w);
	murec_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	murec_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
      if(TauPiRec.E()>0) {
	pirec_PtRatio_hplus.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hplus*UnSpin_WT);
	pirec_PtRatio_hminus.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hminus*UnSpin_WT);
	
	TauPiRec.Boost(-GenZ1.BoostVector());
	Pion_mdR.Boost(-GenZ1.BoostVector());
	GenZ1.Boost(-GenZ1.BoostVector());
	// Now fill results 
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau4_idx));
	pirec_ExoverEtaurec.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w);
	pirec_ExoverEtaurec_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*Spin_WT);
	pirec_ExoverEtaurec_UnSpin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*UnSpin_WT);
	pirec_ExoverEtaurec_FlipSpin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*FlipSpin_WT);
	pirec_ExoverEtaurec_hplus.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hplus*UnSpin_WT);
	pirec_ExoverEtaurec_hminus.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hminus*UnSpin_WT);
	pirec_WT_Spin.at(t).Fill(Spin_WT,w);
	pirec_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	pirec_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
     }   
      
    } //end JAK_MUON && JAK_PION 
      ////////////   end Muon && Pion    ////////////    
    
    
  } //end if(Status)
  
  if(verbose)std::cout << "done" << std::endl;
} //end doEvent()




void  TauSpinStudy::Finish() {
  unsigned int tdata=0;
  TF1 fmu_ExoverEtau_hplus( "fmu_ExoverEtau_hplus", "5.0/3.0-3.0*(x^2.0)+(4.0/3.0)*x^3.0-1*(-1/3+3*x^2-8/3*x^3)",0,1);
  TF1 fmu_ExoverEtau_hminus("fmu_ExoverEtau_hminus","5.0/3.0-3.0*(x^2.0)+(4.0/3.0)*x^3.0+1*(-1/3+3*x^2-8/3*x^3)",0,1);
  TF1 fpi_ExoverEtau_hplus( "fpi_ExoverEtau_hplus", "x",0,1);
  TF1 fpi_ExoverEtau_hminus("fpi_ExoverEtau_hminus","1-x",0,1);
  for(int i=1;i<mu_ExoverEtau_hplus.at(tdata).GetNbinsX();i++) {
    double min=mu_ExoverEtau_hplus.at(tdata).GetBinLowEdge(i);
    double max=min+mu_ExoverEtau_hplus.at(tdata).GetBinWidth(i);
    mu_ExoverEtau_hplus.at(tdata).SetBinContent(i,fmu_ExoverEtau_hplus.Integral(min,max));
    mu_ExoverEtau_hplus.at(tdata).SetBinError(i,0.00001);
    mu_ExoverEtau_hminus.at(tdata).SetBinContent(i,fmu_ExoverEtau_hminus.Integral(min,max));
    mu_ExoverEtau_hminus.at(tdata).SetBinError(i,0.00001);
    pi_ExoverEtau_hplus.at(tdata).SetBinContent(i,fpi_ExoverEtau_hplus.Integral(min,max));
    pi_ExoverEtau_hplus.at(tdata).SetBinError(i,0.00001);
    pi_ExoverEtau_hminus.at(tdata).SetBinContent(i,fpi_ExoverEtau_hminus.Integral(min,max));
    pi_ExoverEtau_hminus.at(tdata).SetBinError(i,0.00001);
  }











  Selection::Finish();
}



     



//////++++++++++++++++++++++++++++//////
////// Z ---> Tau Tau ----> Pi Mu ////// 
//////++++++++++++++++++++++++++++//////

// if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,tau3_idx,tau4_idx) ) {
//   if((Ntp->MCTau_JAK(tau3_idx)==TauDecay::JAK_PION && Ntp->MCTau_JAK(tau4_idx)==TauDecay::JAK_MUON) || (Ntp->MCTau_JAK(tau3_idx)==TauDecay::JAK_MUON && Ntp->MCTau_JAK(tau4_idx)==TauDecay::JAK_PION)) {
// 	if(verbose)std::cout<< "pi-mu " <<std::endl;
// 	TLorentzVector Zboson=Ntp->MCSignalParticle_p4(Boson_idx);
// 	TLorentzVector LVTau=Ntp->MCTau_p4(tau4_idx);
// 	LVTau.Boost(-1*Zboson.BoostVector());

/////////////////// end of Z ---> Tau Tau ----> Pi Mu /////////////////////////////
///////////////////++++++++++++++++++++++++++++++++++ /////////////////////////////


    // ////////////   Muon    ////////////
    // if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,TauDecay::JAK_MUON,tau_idx)) {
    //   if(verbose)std::cout << "muon" << std::endl;
    //   TLorentzVector GenZ=Ntp->MCSignalParticle_p4(Boson_idx);
    //   TLorentzVector GenTauMu(0,0,0,0);
    //   TLorentzVector GenMu(0,0,0,0);  
    //   TLorentzVector GenNum(0,0,0,0);     
    //   TLorentzVector GenNutm(0,0,0,0);
    //   TLorentzVector Muon(0,0,0,0); 
    //   bool hasMurec=false; bool hasMugen=false; 
    //   double dR_Muons=-1;   
            
    //   //Gen Infos      
    //   for(int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++) {
    //     if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::tau_minus)) {
    //       GenTauMu=Ntp->MCTauandProd_p4(tau_idx,i);
    //     }
    //     else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::mu_plus)) {
    //       GenMu+=Ntp->MCTauandProd_p4(tau_idx,i);
    // 	  hasMugen=true;
    // 	}
    // 	else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::nu_mu)) {
    //       GenNum+=Ntp->MCTauandProd_p4(tau_idx,i);
    // 	}
    // 	else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::nu_tau)) {
    //       GenNutm+=Ntp->MCTauandProd_p4(tau_idx,i);
    // 	}
    //   } //end for loop 
    //   if(GenTauMu.E()>0) {
    //     mu_PtRatio_hplus.at(t).Fill(GenMu.Pt()/GenTauMu.Pt(),w*hplus*UnSpin_WT);
    //     mu_PtRatio_hminus.at(t).Fill(GenMu.Pt()/GenTauMu.Pt(),w*hminus*UnSpin_WT);
    // 	//GenTauMu.Boost(-GenZ.BoostVector());
    //     //GenMu.Boost(-GenZ.BoostVector());
    //     //GenZ.Boost(-GenZ.BoostVector());
    //     // Now fill results                                                                                        
    //     Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
    // 	mu_ExoverEtau.at(t).Fill(GenMu.E()/GenTauMu.E(),w);
    // 	mu_ExoverEtau_Spin.at(t).Fill(GenMu.E()/GenTauMu.E(),w*Spin_WT);
    // 	mu_ExoverEtau_UnSpin.at(t).Fill(GenMu.E()/GenTauMu.E(),w*UnSpin_WT);
    // 	mu_ExoverEtau_FlipSpin.at(t).Fill(GenMu.E()/GenTauMu.E(),w*FlipSpin_WT);
    // 	mu_ExoverEtau_hplus.at(t).Fill(GenMu.E()/GenTauMu.E(),w*hplus*UnSpin_WT);
    // 	mu_ExoverEtau_hminus.at(t).Fill(GenMu.E()/GenTauMu.E(),w*hminus*UnSpin_WT);
    // 	mu_WT_Spin.at(t).Fill(Spin_WT,w);
    // 	mu_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
    // 	mu_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
	
    // 	if(verbose) {
    // 	  for(int i=0;i<Ntp->NMCTauDecayProducts(tau_idx);i++) {
    // 	    TLorentzVector LV=Ntp->MCTauandProd_p4(tau_idx,i);
    // 	    std::cout <<"MC_ID="<<Ntp->MCTauandProd_pdgid(tau_idx,i)<< "||Px="<<LV.Px()<<"||Py="<<LV.Py()<<"||Pz="<<LV.Pz()<<"||E="<<LV.E()<< std::endl;
    // 	  }
    // 	}
    // 	if(verbose)std::cout <<"mu:x="<<GenMu.E()/GenTauMu.E()<<"||w="<<w<<"||hplus="<<hplus<<"||hminus="<<hminus <<std::endl;	
    //   }
    //   if(hasMugen) {
    // 	M_Z_gen.at(t).Fill(GenZ.M(),w);
    // 	Pt_Z_gen.at(t).Fill(GenZ.Pt(),w);            
    // 	Eta_Z_gen.at(t).Fill(GenZ.Eta(),w);
    //  	Phi_Z_gen.at(t).Fill(GenZ.Phi(),w);
    // 	M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);
    // 	Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);   
    // 	Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);
    //  	Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w);
    // 	M_Mu_gen.at(t).Fill(GenMu.M(),w); 
    // 	Pt_Mu_gen.at(t).Fill(GenMu.Pt(),w);
    // 	Eta_Mu_gen.at(t).Fill(GenMu.Eta(),w); 
    //  	Phi_Mu_gen.at(t).Fill(GenMu.Phi(),w); 
    // 	M_Nutm_gen.at(t).Fill(GenNutm.M(),w); 
    // 	Pt_Nutm_gen.at(t).Fill(GenNutm.Pt(),w);  
    // 	Eta_Nutm_gen.at(t).Fill(GenNutm.Eta(),w); 
    //  	Phi_Nutm_gen.at(t).Fill(GenNutm.Phi(),w);
    // 	M_Num_gen.at(t).Fill(GenNum.M(),w);
    // 	Pt_Num_gen.at(t).Fill(GenNum.Pt(),w);    
    // 	Eta_Num_gen.at(t).Fill(GenNum.Eta(),w); 
    //  	Phi_Num_gen.at(t).Fill(GenNum.Phi(),w);  
    //   }
    //   //Reco Muon
    //   for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++) {
    // 	//std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
    // 	Muon = Ntp->Muon_p4(iMuon);
    // 	hasMurec = true;
    //   }
    //   if(hasMurec) {
    // 	Mu_pt.at(t).Fill(Muon.Pt(),w);
    // 	Mu_phi.at(t).Fill(Muon.Phi(),w);
    // 	Mu_eta.at(t).Fill(Muon.Eta(),w);
    //   }	
    //   if(hasMugen) {
    // 	for(unsigned int i=0; i < Ntp->NMuons(); i++) {
    // 	  std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
    // 	  dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu) ;
    // 	  std::cout<< "DR_Muon==" << dR_Muons <<std::endl;
    // 	  dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
    // 	}  
    //   }  
      
    //   if(verbose)std::cout<< "muon done " <<std::endl;
    // } //end JAK_MUON
    // ////////////  end Muon  ////////////
    
    // ////////////   Pion    ////////////
    // if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,TauDecay::JAK_PION,tau_idx)) {
    //   if(verbose)std::cout << "pion" << std::endl;
    //   TLorentzVector GenZ1=Ntp->MCSignalParticle_p4(Boson_idx);
    //   TLorentzVector GenTauPi(0,0,0,0);
    //   TLorentzVector GenPi(0,0,0,0);
    //   TLorentzVector GenNutp(0,0,0,0);
    //   TLorentzVector Pion(0,0,0,0);
    //   bool hasPirec=false; bool hasPigen=false; 
    //   double dR_Pions=-1;       
      
    //   //Gen Pion
    //   for(int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++) {
    //     if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::tau_minus)) {
    //       GenTauPi=Ntp->MCTauandProd_p4(tau_idx,i);
    //     }
    //     else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) ) {
    //       GenPi+=Ntp->MCTauandProd_p4(tau_idx,i);
    // 	  hasPigen = true;
    // 	}
    //     else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::nu_tau) ) {
    // 	  GenNutp+=Ntp->MCTauandProd_p4(tau_idx,i);
    // 	}
    //   }      
    //   if(GenTauPi.E()>0) {
    // 	pi_PtRatio_hplus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hplus*UnSpin_WT);
    //     pi_PtRatio_hminus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hminus*UnSpin_WT);
    // 	//GenTauPi.Boost(-GenZ1.BoostVector());
    //     //GenPi.Boost(-GenZ1.BoostVector());
    //     //Boson_LV.Boost(-GenZ1.BoostVector());
    //     // Now fill results                                   
    //     Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
    // 	pi_ExoverEtau.at(t).Fill(GenPi.E()/GenTauPi.E(),w);
    // 	pi_ExoverEtau_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*Spin_WT);
    // 	pi_ExoverEtau_UnSpin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*UnSpin_WT);
    // 	pi_ExoverEtau_FlipSpin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*FlipSpin_WT);
    // 	pi_ExoverEtau_hplus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hplus*UnSpin_WT);
    // 	pi_ExoverEtau_hminus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hminus*UnSpin_WT);
    // 	pi_WT_Spin.at(t).Fill(Spin_WT,w);
    // 	pi_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
    // 	pi_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
	
    // 	if(verbose)std::cout << "pi : x " << GenPi.E()/GenTauPi.E() << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << std::endl;
    //   }
    //   if(hasPigen) {
    // 	M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);
    // 	Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);
    // 	Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w); 
    //  	Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w);
	
    //  	M_Pi_gen.at(t).Fill(GenPi.M(),w); 
    // 	Pt_Pi_gen.at(t).Fill(GenPi.Pt(),w);
    // 	Eta_Pi_gen.at(t).Fill(GenPi.Eta(),w); 
    //  	Phi_Pi_gen.at(t).Fill(GenPi.Phi(),w); 
	
    //  	M_Nutp_gen.at(t).Fill(GenNutp.M(),w);
    //     Pt_Nutp_gen.at(t).Fill(GenNutp.Pt(),w); 
    // 	Eta_Nutp_gen.at(t).Fill(GenNutp.Eta(),w); 
    //  	Phi_Nutp_gen.at(t).Fill(GenNutp.Phi(),w); 
    //   }
      
    //   //Reco Pion
    //   for(unsigned int iPion=0; iPion < Ntp->NPFTaus(); iPion++) {
    // 	//std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
    // 	Pion = Ntp->PFTau_p4(iPion);
    // 	hasPirec  = true;
    //   }
    //   if (hasPirec) {
    // 	Pion_pt.at(t).Fill(Pion.Pt(),w);  
    // 	Pion_phi.at(t).Fill(Pion.Phi(),w);
    // 	Pion_eta.at(t).Fill(Pion.Eta(),w);
    //   }
    //   if(hasPigen) {
    // 	for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
    //  	  std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
    //  	  dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi) ;
    //  	  std::cout<< "DR_Pion==" << dR_Pions <<std::endl;
    //  	  dR_recPi_genPi.at(t).Fill(dR_Pions,w); 
    //  	}
    //   }
      
    //   if(verbose)std::cout<< "pion done " <<std::endl;
    // } // end JAK_PION
    // ////////////////// end Poin ////////////////////////////



    // ////////////////// Muon && Poin ////////////////////////////    
    // if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,tau3_idx,tau4_idx) ) {
    // 	 if((Ntp->MCTau_JAK(tau3_idx)==TauDecay::JAK_PION && Ntp->MCTau_JAK(tau4_idx)==TauDecay::JAK_MUON) || (Ntp->MCTau_JAK(tau3_idx)==TauDecay::JAK_MUON && Ntp->MCTau_JAK(tau4_idx)==TauDecay::JAK_PION)) {
    // 	   if(verbose)std::cout<< "pi-mu " <<std::endl;
    // 	   TLorentzVector Zboson=Ntp->MCSignalParticle_p4(Boson_idx);
    // 	   TLorentzVector LVTau=Ntp->MCTau_p4(tau4_idx);
    // 	   LVTau.Boost(-1*Zboson.BoostVector());
    
    // 	 }
    // }    
 
    // ////////////////// end Muon && Poin ////////////////////////////    


















///////////////////////////////////////////////////////////////////////
// if(status) { 
// TLorentzVector Pion;     TLorentzVector Muon; 
// bool hasPionrec=false;   bool hasMurec=false; 
// double dR_Pions=-1;      double dR_Muons=-1;     
// TLorentzVector GenTauMu, GenTauPi, GenZ, GenPi, GenMu, GenNutp, GenNutm, GenNum;
// TLorentzVector DiTau, MuPi;
// double Delta_M_DiTau, dEta_genTaus, dPhi_genTaus, dTheta_genTaus, dR_genTaus;
// double dEta_genMuPi, dPhi_genMuPi, dTheta_genMuPi, dR_genMu_genPi;
// bool hasZ    = false;  bool hasTaumu = false;  bool hasTaupi = false; 
// bool hasMu   = false;  bool hasPi    = false; 
    
// if(id != DataMCType::Data && (id == DataMCType::DY_tautau )) { //DY_tautau
//   if(Ntp->NMCTaus() == 2) {  //  NMCTaus() == 2       
	
// 	////////////////////////////
// 	////     Reco Study     ////    
// 	/// Pion and Muon not classified 
// 	for(unsigned int iPion=0; iPion < Ntp->NPFTaus(); iPion++) {
// 	  //std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
// 	  Pion = Ntp->PFTau_p4(iPion);
// 	  hasPionrec  = true;
// 	}
// 	for(unsigned int iMuon=0; iMuon < Ntp->NMuons(); iMuon++) {
// 	  //std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
// 	  Muon = Ntp->Muon_p4(iMuon);
// 	  hasMurec = true;
// 	}
// 	if (hasPionrec && hasMurec) {
// 	  Pion_pt.at(t).Fill(Pion.Pt(),w);  
// 	  Pion_phi.at(t).Fill(Pion.Phi(),w);
// 	  Pion_eta.at(t).Fill(Pion.Eta(),w);
// 	  Pion_theta.at(t).Fill(Pion.Theta(),w);
	  
// 	  Mu_pt.at(t).Fill(Muon.Pt(),w);
// 	  Mu_phi.at(t).Fill(Muon.Phi(),w);
// 	  Mu_eta.at(t).Fill(Muon.Eta(),w);
// 	  Mu_theta.at(t).Fill(Muon.Theta(),w);
// 	}
// 	////////////////////////////
// 	//// Generator Study my ////
// 	// if(id != DataMCType::Data && (id == DataMCType::DY_tautau )) { //DY_tautau
// 	//if(Ntp->NMCTaus() == 2) {  //  NMCTaus() == 2  
	
// 	for(unsigned int i=0; i<Ntp->NMCParticles(); i++) { // loop over MC Particles 
// 	  if(Ntp->MCParticle_pdgid(i) == PDGInfo::Z0) {     // Z bososn 
// 	    GenZ = Ntp->MCParticle_p4(i);
// 	    hasZ = true;
// 	  }	
// 	}   
// 	for(int i=0; i<Ntp->NMCTaus(); i++) {  //loop over NMCTaus
// 	  if(Ntp->MCTau_JAK(i) == TauDecay::JAK_MUON) {         //Tau->Muon   //JAK_MUON=2
// 	    GenTauMu = Ntp->MCTau_p4(i);
// 	    hasTaumu=true;
	    
// 	    for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++) {
// 	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::mu_minus) {
// 		GenMu = Ntp->MCTauandProd_p4(i,j); 
// 		hasMu = true; 
// 	      } 
// 	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_mu) {
// 		GenNum = Ntp->MCTauandProd_p4(i,j);
// 	      } 
// 	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_tau) {
// 		GenNutm = Ntp->MCTauandProd_p4(i,j);
// 	      }
// 	    }
// 	  } //end Tau->Muon   //JAK_MUON=2
	  
// 	  if(Ntp->MCTau_JAK(i) == TauDecay::JAK_PION) {//Tau->Pion  //JAK_PION=3
// 	    GenTauPi = Ntp->MCTau_p4(i);
// 	    hasTaupi = true;
	    
// 	    for(int j=0; j<Ntp->NMCTauDecayProducts(i); j++) {
// 	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::pi_plus) {
// 		GenPi = Ntp->MCTauandProd_p4(i,j); 
// 		hasPi = true;
// 	      }
// 	      if(fabs(Ntp->MCTauandProd_pdgid(i,j)) == PDGInfo::nu_tau) {
// 		GenNutp = Ntp->MCTauandProd_p4(i,j);
// 	      }
// 	    }
// 	  }  //end Tau->Pion  //JAK_PION=3
// 	} //end loop over NMCTaus
//   }//end Taus && NMCTaus() == 2
      
//   if(hasZ && hasMu && hasPi ) { // fill variables
// 	M_Z_gen.at(t).Fill(GenZ.M(),w);              Pt_Z_gen.at(t).Fill(GenZ.Pt(),w);                  Eta_Z_gen.at(t).Fill(GenZ.Eta(),w);
// 	Phi_Z_gen.at(t).Fill(GenZ.Phi(),w);          Theta_Z_gen.at(t).Fill(GenZ.Theta(),w);
	
// 	M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);      Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);          Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);
// 	Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w);  Theta_TauMu_gen.at(t).Fill(GenTauMu.Theta(),w);
// 	M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);      Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);          Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w); 
// 	Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w);  Theta_TauPi_gen.at(t).Fill(GenTauPi.Theta(),w); 
	
// 	M_Mu_gen.at(t).Fill(GenMu.M(),w);            Pt_Mu_gen.at(t).Fill(GenMu.Pt(),w);                Eta_Mu_gen.at(t).Fill(GenMu.Eta(),w); 
// 	Phi_Mu_gen.at(t).Fill(GenMu.Phi(),w);        Theta_Mu_gen.at(t).Fill(GenMu.Theta(),w); 
// 	M_Pi_gen.at(t).Fill(GenPi.M(),w);   	     Pt_Pi_gen.at(t).Fill(GenPi.Pt(),w);                Eta_Pi_gen.at(t).Fill(GenPi.Eta(),w); 
// 	Phi_Pi_gen.at(t).Fill(GenPi.Phi(),w);        Theta_Pi_gen.at(t).Fill(GenPi.Theta(),w);    	
	
// 	M_Nutp_gen.at(t).Fill(GenNutp.M(),w);        Pt_Nutp_gen.at(t).Fill(GenNutp.Pt(),w);            Eta_Nutp_gen.at(t).Fill(GenNutp.Eta(),w); 
// 	Phi_Nutp_gen.at(t).Fill(GenNutp.Phi(),w);    Theta_Nutp_gen.at(t).Fill(GenNutp.Theta(),w); 
// 	M_Nutm_gen.at(t).Fill(GenNutm.M(),w);        Pt_Nutm_gen.at(t).Fill(GenNutm.Pt(),w);            Eta_Nutm_gen.at(t).Fill(GenNutm.Eta(),w); 
// 	Phi_Nutm_gen.at(t).Fill(GenNutm.Phi(),w);    Theta_Nutm_gen.at(t).Fill(GenNutm.Theta(),w); 
// 	M_Num_gen.at(t).Fill(GenNum.M(),w);          Pt_Num_gen.at(t).Fill(GenNum.Pt(),w);             	Eta_Num_gen.at(t).Fill(GenNum.Eta(),w); 
// 	Phi_Num_gen.at(t).Fill(GenNum.Phi(),w);      Theta_Num_gen.at(t).Fill(GenNum.Theta(),w); 
// 	//DiTau_gen
// 	DiTau = GenTauPi + GenTauMu;
// 	M_DiTau_gen.at(t).Fill(DiTau.M(),w);                   Pt_DiTau_gen.at(t).Fill(DiTau.Pt(),w);
// 	Delta_M_DiTau = DiTau.M() - GenZ.M();                  dM_DiTau_gen.at(t).Fill(Delta_M_DiTau,w);	
// 	dEta_genTaus = GenTauPi.Eta() - GenTauMu.Eta();        dEta_DiTau_gen.at(t).Fill(dEta_genTaus,w);
// 	dPhi_genTaus = GenTauPi.Phi() - GenTauMu.Phi();        dPhi_DiTau_gen.at(t).Fill(dPhi_genTaus,w);
// 	dTheta_genTaus = GenTauPi.Theta() - GenTauMu.Theta();  dTheta_DiTau_gen.at(t).Fill(dTheta_genTaus,w);
// 	dR_genTaus = sqrt(pow(dEta_genTaus,2) + pow(dPhi_genTaus,2));      
// 	dR_DiTau_gen.at(t).Fill(dR_genTaus,w); -
// 	//MuPi_gen
// 	MuPi = GenMu + GenPi;
// 	M_MuPi_gen.at(t).Fill(GenMu.M()+GenPi.M(),w);          Pt_MuPi_gen.at(t).Fill(GenMu.Pt()+GenPi.Pt(),w);
// 	dEta_genMuPi = GenMu.Eta() - GenPi.Eta();              dEta_MuPi_gen.at(t).Fill(dEta_genMuPi,w);
// 	dPhi_genMuPi = GenMu.Phi() - GenPi.Phi();              dPhi_MuPi_gen.at(t).Fill(dPhi_genMuPi,w); 
// 	dTheta_genMuPi = GenMu.Theta() - GenPi.Theta();        dTheta_MuPi_gen.at(t).Fill(dTheta_genMuPi,w);
// 	dR_genMu_genPi = sqrt(pow(GenMu.Eta() - GenPi.Eta(),2) + pow( GenMu.Phi() - GenPi.Phi(),2));      
// 	//dR_genMu_genPi = sqrt(pow(dEta_genMuPi,2) + pow(dPhi_genMuPi,2));      
// 	dR_MuPi_gen.at(t).Fill(dR_genMu_genPi,w);
      
//   } // end fill variables
      
//   /////////////////////////////////////////
//   // DR(Pirec , Pigen)   , DR(Murec ,Mugen)
//   ////////////////////////////////////////
//   if (hasMu && hasPi) {
// 	for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
// 	  std::cout<< "NPion==" << Ntp->NPFTaus() <<std::endl;
// 	  dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi) ;
// 	  std::cout<< "DR_Pion==" << dR_Pions <<std::endl;
// 	  dR_recPi_genPi.at(t).Fill(dR_Pions,w); 
// 	}
// 	for(unsigned int i=0; i < Ntp->NMuons(); i++) {
// 	  std::cout<< "NMuons==" << Ntp->NMuons() <<std::endl;
// 	  dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu) ;
// 	  std::cout<< "DR_Muon==" << dR_Muons <<std::endl;
// 	  dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
// 	}
//   }
// } //end DY_tautau
// } // End if(status)


































