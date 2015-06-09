#include "PolarizedTau.h"  
#include "TLorentzVector.h"   
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TF1.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSpinerInterface.h"
                 
PolarizedTau::PolarizedTau(TString Name_, TString id_):
  Selection(Name_,id_)
{   
  //verbose=true;
}

PolarizedTau::~PolarizedTau() {
  for(unsigned int j=0; j<Npassed.size(); j++) {
    std::cout << "PolarizedTau::~PolarizedTau Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "PolarizedTau::~PolarizedTau()" << std::endl;
}
      
void  PolarizedTau::Configure() {
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

  /////// Tau_pi --> Pi + Nu //////  
  spin_WT=HConfig.GetTH1D(Name+"_spin_WT","spin_WT",20,0.0,2.1,"spin_WT","Events"); 

  //////////////////////////////////////
  // Gen_Infos /////// Lab_Frame ///////    
  // X1 = Pi_Pt / Tau_Pt
  pi_PtRatio=HConfig.GetTH1D(Name+"_pi_PtRatio","pi_PtRatio",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}","Events");
  pi_PtRatio_hplus=HConfig.GetTH1D(Name+"_pi_PtRatio_hplus","pi_PtRatio_hplus",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{+}}","Events");
  pi_PtRatio_hminus=HConfig.GetTH1D(Name+"_pi_PtRatio_hminus","pi_PtRatio_hminus",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{-}}","Events");

  pi_PtRatio_Spin=HConfig.GetTH1D(Name+"_pi_PtRatio_Spin","pi_PtRatio_Spin",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_spin","Events");
  pi_PtRatio_hplus_Spin=HConfig.GetTH1D(Name+"_pi_PtRatio_hplus_Spin","pi_PtRatio_hplus_Spin",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{+}}|_spin","Events");
  pi_PtRatio_hminus_Spin=HConfig.GetTH1D(Name+"_pi_PtRatio_hminus_Spin","pi_PtRatio_hminus_Spin",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{-}}|_spin","Events");
  // X2 = Pi_Pt / Z_Mt
  pi_PtoverZmt=HConfig.GetTH1D(Name+"_pi_PtoverZmt","pi_PtoverZmt",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}","Events");
  pi_PtoverZmt_hplus=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hplus","pi_PtoverZmt_hplus",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{+}}","Events");
  pi_PtoverZmt_hminus=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hminus","pi_PtoverZmt_hminus",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{-}}","Events");
  pi_PtoverZmt_Spin=HConfig.GetTH1D(Name+"_pi_PtoverZmt_Spin","pi_PtoverZmt_Spin",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_spin","Events");
  pi_PtoverZmt_hplus_Spin=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hplus_Spin","pi_PtoverZmt_hplus_Spin",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{+}}|_spin","Events");
  pi_PtoverZmt_hminus_Spin=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hminus_Spin","pi_PtoverZmt_hminus_Spin",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{-}}|_spin","Events");
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  pi_EoverEtau=HConfig.GetTH1D(Name+"_pi_EoverEtau","pi_EoverEtau",20,0.0,1.0,"E_{#pi}/E_{#tau}","Events");
  pi_EoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_EoverEtau_hplus","pi_EoverEtau_hplus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}","Events");
  pi_EoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_EoverEtau_hminus","pi_EoverEtau_hminus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}","Events");
  pi_EoverEtau_Spin=HConfig.GetTH1D(Name+"_pi_EoverEtau_Spin","pi_EoverEtau_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_spin","Events");
  pi_EoverEtau_hplus_Spin=HConfig.GetTH1D(Name+"_pi_EoverEtau_hplus_Spin","pi_EoverEtau_hplus_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}|_spin","Events");
  pi_EoverEtau_hminus_Spin=HConfig.GetTH1D(Name+"_pi_EoverEtau_hminus_Spin","pi_EoverEtau_hminus_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}|_spin","Events");
  // X4 = Pi_E / Z_M
  pi_EoverZm=HConfig.GetTH1D(Name+"_pi_EoverZm","pi_EoverZm",20,0.0,1.0,"E_{#pi}/M_{Z}","Events");
  pi_EoverZm_hplus=HConfig.GetTH1D(Name+"_pi_EoverZm_hplus","pi_EoverZm_hplus",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{+}}","Events");
  pi_EoverZm_hminus=HConfig.GetTH1D(Name+"_pi_EoverZm_hminus","pi_EoverZm_hminus",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{-}}","Events");
  pi_EoverZm_Spin=HConfig.GetTH1D(Name+"_pi_EoverZm_Spin","pi_EoverZm_Spin",20,0.0,1.0,"E_{#pi}/M_{Z}|_spin","Events");
  pi_EoverZm_hplus_Spin=HConfig.GetTH1D(Name+"_pi_EoverZm_hplus_Spin","pi_EoverZm_hplus_Spin",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{+}}|_spin","Events");
  pi_EoverZm_hminus_Spin=HConfig.GetTH1D(Name+"_pi_EoverZm_hminus_Spin","pi_EoverZm_hminus_Spin",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{-}}|_spin","Events");
  
  ////////////////////////////////////// 
  // Rec_Infos /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  pirec_PtRatio=HConfig.GetTH1D(Name+"_pirec_PtRatio","pirec_PtRatio",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}","Events");
  pirec_PtRatio_hplus=HConfig.GetTH1D(Name+"_pirec_PtRatio_hplus","pirec_PtRatio_hplus",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{+}}","Events");
  pirec_PtRatio_hminus=HConfig.GetTH1D(Name+"_pirec_PtRatio_hminus","pirec_PtRatio_hminus",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{-}}","Events");
  pirec_PtRatio_Spin=HConfig.GetTH1D(Name+"_pirec_PtRatio_Spin","pirec_PtRatio_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_spin","Events");
  pirec_PtRatio_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_PtRatio_hplus_Spin","pirec_PtRatio_hplus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{+}}|_spin","Events");
  pirec_PtRatio_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_PtRatio_hminus_Spin","pirec_PtRatio_hminus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{-}}|_spin","Events");
  
  // X2 = Pirec_Pt / Zrec_Mt
  pirec_PtoverZmt=HConfig.GetTH1D(Name+"_pirec_PtoverZmt","pirec_PtoverZmt",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}","Events");
  pirec_PtoverZmt_hplus=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hplus","pirec_PtoverZmt_hplus",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{+}}","Events");
  pirec_PtoverZmt_hminus=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hminus","pirec_PtoverZmt_hminus",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{-}}","Events");
  pirec_PtoverZmt_Spin=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_Spin","pirec_PtoverZmt_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_spin","Events");
  pirec_PtoverZmt_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hplus_Spin","pirec_PtoverZmt_hplus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{+}}|_spin","Events");
  pirec_PtoverZmt_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hminus_Spin","pirec_PtoverZmt_hminus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{-}}|_spin","Events");
  
  ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  pirec_EoverEtau=HConfig.GetTH1D(Name+"_pirec_EoverEtau","pirec_EoverEtau",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}","Events");
  pirec_EoverEtau_hplus=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hplus","pirec_EoverEtau_hplus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}","Events");
  pirec_EoverEtau_hminus=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hminus","pirec_EoverEtau_hminus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}","Events");
  pirec_EoverEtau_Spin=HConfig.GetTH1D(Name+"_pirec_EoverEtau_Spin","pirec_EoverEtau_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_spin","Events");
  pirec_EoverEtau_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hplus_Spin","pirec_EoverEtau_hplus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}|_spin","Events");
  pirec_EoverEtau_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hminus_Spin","pirec_EoverEtau_hminus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}|_spin","Events");
  
  // X4 = Pirec_E / Zrec_M
  pirec_EoverZm=HConfig.GetTH1D(Name+"_pirec_EoverZm","pirec_EoverZm",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}","Events");
  pirec_EoverZm_hplus=HConfig.GetTH1D(Name+"_pirec_EoverZm_hplus","pirec_EoverZm_hplus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}","Events");
  pirec_EoverZm_hminus=HConfig.GetTH1D(Name+"_pirec_EoverZm_hminus","pirec_EoverZm_hminus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}","Events");
  pirec_EoverZm_Spin=HConfig.GetTH1D(Name+"_pirec_EoverZm_Spin","pirec_EoverZm_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_spin","Events");
  pirec_EoverZm_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverZm_hplus_Spin","pirec_EoverZm_hplus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}|_spin","Events");
  pirec_EoverZm_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverZm_hminus_Spin","pirec_EoverZm_hminus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}|_spin","Events");
  
    
  // Setup Extra Histograms
  // Reco Objects
  //Delta_R(Mu_rec,Mu_gen), Delta_R(Pi_rec,Pi_gen)
  dR_recMu_genMu=HConfig.GetTH1D(Name+"_dR_recMu_genMu","dR_recMu_genMu",100,0.,8.,"dR_recMu_genMu","Events");
  mindR_recMu_genMu=HConfig.GetTH1D(Name+"_mindR_recMu_genMu","mindR_recMu_genMu",100,0.,1.,"mindR_recMu_genMu","Events");
  dR_recPi_genPi=HConfig.GetTH1D(Name+"_dR_recPi_genPi","dR_recPi_genPi",100,0.,8.,"dR_recPi_genPi","Events");
  mindR_recPi_genPi=HConfig.GetTH1D(Name+"_mindR_recPi_genPi","mindR_recPi_genPi",100,0.,1.,"mindR_recPi_genPi","Events");

  // DR(vis,mis)
  dR_Mu_MisNum_10_15=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_10_15","dR_Mu_MisNum_10_15",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_15_20=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_15_20","dR_Mu_MisNum_15_20",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_20_25=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_20_25","dR_Mu_MisNum_20_25",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_25_30=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_25_30","dR_Mu_MisNum_25_30",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_30_35=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_30_35","dR_Mu_MisNum_30_35",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_35_40=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_35_40","dR_Mu_MisNum_35_40",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_40_45=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_40_45","dR_Mu_MisNum_40_45",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_45_50=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_45_50","dR_Mu_MisNum_45_50",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_50_55=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_50_55","dR_Mu_MisNum_50_55",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_55_60=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_55_60","dR_Mu_MisNum_55_60",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_60_65=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_60_65","dR_Mu_MisNum_60_65",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_65_70=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_65_70","dR_Mu_MisNum_65_70",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_70_75=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_70_75","dR_Mu_MisNum_70_75",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_75_80=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_75_80","dR_Mu_MisNum_75_80",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_80_85=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_80_85","dR_Mu_MisNum_80_85",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_85_90=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_85_90","dR_Mu_MisNum_85_90",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_90_95=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_90_95","dR_Mu_MisNum_90_95",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Mu_MisNum_95_100=HConfig.GetTH1D(Name+"_dR_Mu_MisNum_95_100","dR_Mu_MisNum_95_100",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");

  dR_Pi_MisNup_10_15=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_10_15","dR_Pi_MisNup_10_15",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_15_20=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_15_20","dR_Pi_MisNup_15_20",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_20_25=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_20_25","dR_Pi_MisNup_20_25",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_25_30=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_25_30","dR_Pi_MisNup_25_30",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_30_35=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_30_35","dR_Pi_MisNup_30_35",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_35_40=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_35_40","dR_Pi_MisNup_35_40",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_40_45=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_40_45","dR_Pi_MisNup_40_45",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_45_50=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_45_50","dR_Pi_MisNup_45_50",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_50_55=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_50_55","dR_Pi_MisNup_50_55",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_55_60=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_55_60","dR_Pi_MisNup_55_60",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_60_65=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_60_65","dR_Pi_MisNup_60_65",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_65_70=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_65_70","dR_Pi_MisNup_65_70",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_70_75=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_70_75","dR_Pi_MisNup_70_75",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_75_80=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_75_80","dR_Pi_MisNup_75_80",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_80_85=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_80_85","dR_Pi_MisNup_80_85",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_85_90=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_85_90","dR_Pi_MisNup_85_90",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_90_95=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_90_95","dR_Pi_MisNup_90_95",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");
  dR_Pi_MisNup_95_100=HConfig.GetTH1D(Name+"_dR_Pi_MisNup_95_100","dR_Pi_MisNup_95_100",60,0.,0.8,"#DeltaR(#vec{P}_{vis},#vec{P}_{mis})","Events");

  //Z---> TauPigen                                                                                                                                           
  M_TauPi_gen=HConfig.GetTH1D(Name+"_M_TauPi_gen","M_TauPi_gen",180,0.,4.,"M_TauPi_gen","Events");
  Pt_TauPi_gen=HConfig.GetTH1D(Name+"_Pt_TauPi_gen","Pt_TauPi_gen",30,0.,120.,"Pt_TauPi_gen","Events");
  Eta_TauPi_gen=HConfig.GetTH1D(Name+"_Eta_TauPi_gen","Eta_TauPi_gen",50,-2.5,2.5,"Eta_TauPi_gen","Events");
  Phi_TauPi_gen=HConfig.GetTH1D(Name+"_Phi_TauPi_gen","Phi_TauPi_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_gen","Events");
  //Z---> TauMugen                                                                                                                                           
  M_TauMu_gen=HConfig.GetTH1D(Name+"_M_TauMu_gen","M_TauMu_gen",180,0.,4.,"M_TauMu_gen","Events");
  Pt_TauMu_gen=HConfig.GetTH1D(Name+"_Pt_TauMu_gen","Pt_TauMu_gen",30,0.,120.,"Pt_TauMu_gen","Events");
  Eta_TauMu_gen=HConfig.GetTH1D(Name+"_Eta_TauMu_gen","Eta_TauMu_gen",50,-2.5,2.5,"Eta_TauMu_gen","Events");
  Phi_TauMu_gen=HConfig.GetTH1D(Name+"_Phi_TauMu_gen","Phi_TauMu_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_gen","Events");

  //Z---> TauPirec
  M_TauPi_rec=HConfig.GetTH1D(Name+"_M_TauPi_rec","M_TauPi_rec",180,0.,4.,"M_TauPi_rec","Events");
  Pt_TauPi_rec=HConfig.GetTH1D(Name+"_Pt_TauPi_rec","Pt_TauPi_rec",30,0.,120.,"Pt_TauPi_rec","Events");
  Eta_TauPi_rec=HConfig.GetTH1D(Name+"_Eta_TauPi_rec","Eta_TauPi_rec",50,-2.5,2.5,"Eta_TauPi_rec","Events");
  Phi_TauPi_rec=HConfig.GetTH1D(Name+"_Phi_TauPi_rec","Phi_TauPi_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_rec","Events");
  //Z---> TauMurec
  M_TauMu_rec=HConfig.GetTH1D(Name+"_M_TauMu_rec","M_TauMu_rec",180,0.,4.,"M_TauMu_rec","Events");
  Pt_TauMu_rec=HConfig.GetTH1D(Name+"_Pt_TauMu_rec","Pt_TauMu_rec",30,0.,120.,"Pt_TauMu_rec","Events");
  Eta_TauMu_rec=HConfig.GetTH1D(Name+"_Eta_TauMu_rec","Eta_TauMu_rec",50,-2.5,2.5,"Eta_TauMu_rec","Events");
  Phi_TauMu_rec=HConfig.GetTH1D(Name+"_Phi_TauMu_rec","Phi_TauMu_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_rec","Events");

  //////...Delta Phi and Theta...///////
  dPhi_genTauMu_genTauPi=HConfig.GetTH1D(Name+"_dPhi_genTauMu_genTauPi","dPhi_genTauMu_genTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_gen,#tau_{#pi}_gen)","Events");
  dPhi_genTauMu_genMu=HConfig.GetTH1D(Name+"_dPhi_genTauMu_genMu","dPhi_genTauMu_genMu",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_gen,#mu_gen)","Events");
  dPhi_genTauPi_genPi=HConfig.GetTH1D(Name+"_dPhi_genTauPi_genPi","dPhi_genTauPi_genPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#pi}_gen,#pi_gen)","Events");
  dPhi_genTauMu_genNum=HConfig.GetTH1D(Name+"_dPhi_genTauMu_genNum","dPhi_genTauMu_genNum",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_gen,#nu_{#mu}_gen)","Events");
  dPhi_genTauMu_genNutm=HConfig.GetTH1D(Name+"_dPhi_genTauMu_genTauPi","dPhi_genTauMu_genTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_gen,#nu_{#tau}_gen)","Events");
  dPhi_genTauPi_genNutp=HConfig.GetTH1D(Name+"_dPhi_genTauMu_genTauPi","dPhi_genTauMu_genTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#pi}_gen,#nu_{#tau}_gen)","Events");

  dPhi_recTauMu_recTauPi=HConfig.GetTH1D(Name+"_dPhi_recTauMu_recTauPi","dPhi_recTauMu_recTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_rec,#tau_{#pi}_rec)","Events");
  dPhi_recTauMu_recMu=HConfig.GetTH1D(Name+"_dPhi_recTauMu_recMu","dPhi_recTauMu_recMu",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_rec,#mu_rec)","Events");
  dPhi_recTauPi_recPi=HConfig.GetTH1D(Name+"_dPhi_recTauPi_recPi","dPhi_recTauPi_recPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#pi}_rec,#pi_rec)","Events");
  dPhi_recTauMu_recNum=HConfig.GetTH1D(Name+"_dPhi_recTauMu_recNum","dPhi_recTauMu_recNum",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_rec,#nu_{#mu}_rec)","Events");
  dPhi_recTauMu_recNutm=HConfig.GetTH1D(Name+"_dPhi_recTauMu_recTauPi","dPhi_recTauMu_recTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#mu}_rec,#nu_{#tau}_rec)","Events");
  dPhi_recTauPi_recNutp=HConfig.GetTH1D(Name+"_dPhi_recnTauMu_recTauPi","dPhi_recTauMu_recTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#phi(#tau_{#pi}_rec,#nu_{#tau}_rec)","Events");


  dTheta_genTauMu_genTauPi=HConfig.GetTH1D(Name+"_dTheta_genTauMu_genTauPi","dTheta_genTauMu_genTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_gen,#tau_{#pi}_gen)","Events");
  dTheta_genTauMu_genMu=HConfig.GetTH1D(Name+"_dTheta_genTauMu_genMu","dTheta_genTauMu_genMu",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_gen,#mu_gen)","Events");
  dTheta_genTauPi_genPi=HConfig.GetTH1D(Name+"_dTheta_genTauPi_genPi","dTheta_genTauPi_genPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#pi}_gen,#pi_gen)","Events");
  dTheta_genTauMu_genNum=HConfig.GetTH1D(Name+"_dTheta_genTauMu_genNum","dTheta_genTauMu_genNum",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_gen,#nu_{#mu}_gen)","Events");
  dTheta_genTauMu_genNutm=HConfig.GetTH1D(Name+"_dTheta_genTauMu_genTauPi","dTheta_genTauMu_genTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_gen,#nu_{#tau}_gen)","Events");
  dTheta_genTauPi_genNutp=HConfig.GetTH1D(Name+"_dTheta_genTauMu_genTauPi","dTheta_genTauMu_genTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#pi}_gen,#nu_{#tau}_gen)","Events");

  dTheta_recTauMu_recTauPi=HConfig.GetTH1D(Name+"_dTheta_recTauMu_recTauPi","dTheta_recTauMu_recTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_rec,#tau_{#pi}_rec)","Events");
  dTheta_recTauMu_recMu=HConfig.GetTH1D(Name+"_dTheta_recTauMu_recMu","dTheta_recTauMu_recMu",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_rec,#mu_rec)","Events");
  dTheta_recTauPi_recPi=HConfig.GetTH1D(Name+"_dTheta_recTauPi_recPi","dTheta_recTauPi_recPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#pi}_rec,#pi_rec)","Events");
  dTheta_recTauMu_recNum=HConfig.GetTH1D(Name+"_dTheta_recTauMu_recNum","dTheta_recTauMu_recNum",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_rec,#nu_{#mu}_rec)","Events");
  dTheta_recTauMu_recNutm=HConfig.GetTH1D(Name+"_dTheta_recTauMu_recTauPi","dTheta_recTauMu_recTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#theta(#tau_{#mu}_rec,#nu_{#tau}_rec)","Events");
  dTheta_recTauPi_recNutp=HConfig.GetTH1D(Name+"_dTheta_recnTauMu_recTauPi","dTheta_recTauMu_recTauPi",32,-TMath::Pi(),TMath::Pi(),"#Delta#thetaTheta(#tau_{#pi}_rec,#nu_{#tau}_rec)","Events");

  // extra Histos
  mu_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{+}}","Events");
  mu_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{-}}","Events");
  pi_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}","Events");
  pi_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}","Events");
  


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  PolarizedTau::Store_ExtraDist() {
  ////////////////////////////////////// 

  Extradist1d.push_back(&spin_WT);
  // Gen_Infos /////// Lab_Frame ///////
  // X1 = Pi_Pt / Tau_Pt     
  Extradist1d.push_back(&pi_PtRatio);
  Extradist1d.push_back(&pi_PtRatio_hplus);
  Extradist1d.push_back(&pi_PtRatio_hminus);
  Extradist1d.push_back(&pi_PtRatio_Spin);
  Extradist1d.push_back(&pi_PtRatio_hplus_Spin);
  Extradist1d.push_back(&pi_PtRatio_hminus_Spin);
  // X2 = Pi_Pt / Z_Mt
  Extradist1d.push_back(&pi_PtoverZmt);
  Extradist1d.push_back(&pi_PtoverZmt_hplus);
  Extradist1d.push_back(&pi_PtoverZmt_hminus);
  Extradist1d.push_back(&pi_PtoverZmt_Spin);
  Extradist1d.push_back(&pi_PtoverZmt_hplus_Spin);
  Extradist1d.push_back(&pi_PtoverZmt_hminus_Spin);
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  Extradist1d.push_back(&pi_EoverEtau);
  Extradist1d.push_back(&pi_EoverEtau_hplus);
  Extradist1d.push_back(&pi_EoverEtau_hminus);
  Extradist1d.push_back(&pi_EoverEtau_Spin);
  Extradist1d.push_back(&pi_EoverEtau_hplus_Spin);
  Extradist1d.push_back(&pi_EoverEtau_hminus_Spin);
  // X4 = Pi_E / Z_M
  Extradist1d.push_back(&pi_EoverZm);
  Extradist1d.push_back(&pi_EoverZm_hplus);
  Extradist1d.push_back(&pi_EoverZm_hminus);
  Extradist1d.push_back(&pi_EoverZm_Spin);
  Extradist1d.push_back(&pi_EoverZm_hplus_Spin);
  Extradist1d.push_back(&pi_EoverZm_hminus_Spin);
  ////////////////////////////////////// 
  // Rec_Infos /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  Extradist1d.push_back(&pirec_PtRatio);
  Extradist1d.push_back(&pirec_PtRatio_hplus);
  Extradist1d.push_back(&pirec_PtRatio_hminus);
  Extradist1d.push_back(&pirec_PtRatio_Spin);
  Extradist1d.push_back(&pirec_PtRatio_hplus_Spin);
  Extradist1d.push_back(&pirec_PtRatio_hminus_Spin);
  // X2 = Extradist1d.Push_Back(&Pirec_Pt / Zrec_Mt
  Extradist1d.push_back(&pirec_PtoverZmt);
  Extradist1d.push_back(&pirec_PtoverZmt_hplus);
  Extradist1d.push_back(&pirec_PtoverZmt_hminus);
  Extradist1d.push_back(&pirec_PtoverZmt_Spin);
  Extradist1d.push_back(&pirec_PtoverZmt_hplus_Spin);
  Extradist1d.push_back(&pirec_PtoverZmt_hminus_Spin);
    ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  Extradist1d.push_back(&pirec_EoverEtau);
  Extradist1d.push_back(&pirec_EoverEtau_hplus);
  Extradist1d.push_back(&pirec_EoverEtau_hminus);
  Extradist1d.push_back(&pirec_EoverEtau_Spin);
  Extradist1d.push_back(&pirec_EoverEtau_hplus_Spin);
  Extradist1d.push_back(&pirec_EoverEtau_hminus_Spin);
  // X4 = Extradist1d.Push_Back(&Pirec_E / Zrec_M
  Extradist1d.push_back(&pirec_EoverZm);
  Extradist1d.push_back(&pirec_EoverZm_hplus);
  Extradist1d.push_back(&pirec_EoverZm_hminus);
  Extradist1d.push_back(&pirec_EoverZm_Spin);
  Extradist1d.push_back(&pirec_EoverZm_hplus_Spin);
  Extradist1d.push_back(&pirec_EoverZm_hminus_Spin);
  //Gen objects
  Extradist1d.push_back(&M_TauPi_gen);        Extradist1d.push_back(&Pt_TauPi_gen);
  Extradist1d.push_back(&Eta_TauPi_gen);      Extradist1d.push_back(&Phi_TauPi_gen);
  Extradist1d.push_back(&M_TauMu_gen);        Extradist1d.push_back(&Pt_TauMu_gen);
  Extradist1d.push_back(&Eta_TauMu_gen);      Extradist1d.push_back(&Phi_TauMu_gen);        
  //Reco objects
  Extradist1d.push_back(&M_TauPi_rec);        Extradist1d.push_back(&Pt_TauPi_rec);
  Extradist1d.push_back(&Eta_TauPi_rec);      Extradist1d.push_back(&Phi_TauPi_rec);
  Extradist1d.push_back(&M_TauMu_rec);        Extradist1d.push_back(&Pt_TauMu_rec);
  Extradist1d.push_back(&Eta_TauMu_rec);      Extradist1d.push_back(&Phi_TauMu_rec);        
  // Delta_R (,)
  Extradist1d.push_back(&dR_recMu_genMu);     Extradist1d.push_back(&dR_recPi_genPi);
  Extradist1d.push_back(&mindR_recMu_genMu);  Extradist1d.push_back(&mindR_recPi_genPi);
  // Delta_R(vis,mis)
  Extradist1d.push_back(&dR_Mu_MisNum_10_15); Extradist1d.push_back(&dR_Mu_MisNum_15_20);
  Extradist1d.push_back(&dR_Mu_MisNum_20_25); Extradist1d.push_back(&dR_Mu_MisNum_25_30);
  Extradist1d.push_back(&dR_Mu_MisNum_30_35); Extradist1d.push_back(&dR_Mu_MisNum_35_40);
  Extradist1d.push_back(&dR_Mu_MisNum_40_45); Extradist1d.push_back(&dR_Mu_MisNum_45_50);
  Extradist1d.push_back(&dR_Mu_MisNum_50_55); Extradist1d.push_back(&dR_Mu_MisNum_55_60);
  Extradist1d.push_back(&dR_Mu_MisNum_60_65); Extradist1d.push_back(&dR_Mu_MisNum_65_70);
  Extradist1d.push_back(&dR_Mu_MisNum_70_75); Extradist1d.push_back(&dR_Mu_MisNum_75_80);
  Extradist1d.push_back(&dR_Mu_MisNum_80_85); Extradist1d.push_back(&dR_Mu_MisNum_85_90);
  Extradist1d.push_back(&dR_Mu_MisNum_90_95); Extradist1d.push_back(&dR_Mu_MisNum_95_100);
  Extradist1d.push_back(&dR_Pi_MisNup_10_15); Extradist1d.push_back(&dR_Pi_MisNup_15_20);
  Extradist1d.push_back(&dR_Pi_MisNup_20_25); Extradist1d.push_back(&dR_Pi_MisNup_25_30);
  Extradist1d.push_back(&dR_Pi_MisNup_30_35); Extradist1d.push_back(&dR_Pi_MisNup_35_40);
  Extradist1d.push_back(&dR_Pi_MisNup_40_45); Extradist1d.push_back(&dR_Pi_MisNup_45_50);
  Extradist1d.push_back(&dR_Pi_MisNup_50_55); Extradist1d.push_back(&dR_Pi_MisNup_55_60);
  Extradist1d.push_back(&dR_Pi_MisNup_60_65); Extradist1d.push_back(&dR_Pi_MisNup_65_70);
  Extradist1d.push_back(&dR_Pi_MisNup_70_75); Extradist1d.push_back(&dR_Pi_MisNup_75_80);
  Extradist1d.push_back(&dR_Pi_MisNup_80_85); Extradist1d.push_back(&dR_Pi_MisNup_85_90);
  Extradist1d.push_back(&dR_Pi_MisNup_90_95); Extradist1d.push_back(&dR_Pi_MisNup_95_100);
  //////...Delta Phi and Theta...///////   
  Extradist1d.push_back(&dPhi_genTauMu_genTauPi);   Extradist1d.push_back(&dPhi_genTauMu_genMu);
  Extradist1d.push_back(&dPhi_genTauPi_genPi);      Extradist1d.push_back(&dPhi_genTauMu_genNum);
  Extradist1d.push_back(&dPhi_genTauMu_genNutm);    Extradist1d.push_back(&dPhi_genTauPi_genNutp);
  Extradist1d.push_back(&dPhi_recTauMu_recTauPi);   Extradist1d.push_back(&dPhi_recTauMu_recMu);
  Extradist1d.push_back(&dPhi_recTauPi_recPi);      Extradist1d.push_back(&dPhi_recTauMu_recNum);
  Extradist1d.push_back(&dPhi_recTauMu_recNutm);    Extradist1d.push_back(&dPhi_recTauPi_recNutp);
  Extradist1d.push_back(&dTheta_genTauMu_genTauPi); Extradist1d.push_back(&dTheta_genTauMu_genMu);
  Extradist1d.push_back(&dTheta_genTauPi_genPi);    Extradist1d.push_back(&dTheta_genTauMu_genNum);
  Extradist1d.push_back(&dTheta_genTauMu_genNutm);  Extradist1d.push_back(&dTheta_genTauPi_genNutp);
  Extradist1d.push_back(&dTheta_recTauMu_recTauPi); Extradist1d.push_back(&dTheta_recTauMu_recMu);
  Extradist1d.push_back(&dTheta_recTauPi_recPi);    Extradist1d.push_back(&dTheta_recTauMu_recNum);
  Extradist1d.push_back(&dTheta_recTauMu_recNutm);  Extradist1d.push_back(&dTheta_recTauPi_recNutp);
  // Extra Histos
  Extradist1d.push_back(&mu_ExoverEtau_hplus);Extradist1d.push_back(&mu_ExoverEtau_hminus);
  Extradist1d.push_back(&pi_ExoverEtau_hplus);Extradist1d.push_back(&pi_ExoverEtau_hminus);
  
}

void  PolarizedTau::doEvent() {
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
    double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    double hminus=1-hplus;//Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
    std::cout << "hplus " << hplus << " hminus " << hminus << std::endl;

    //spin_WT.at(t).Fill(Spin_WT,w);


    
    if(verbose)std::cout<< "A " <<std::endl;
    ////////////////////////////////////////////////
    //
    // Spin Validation
    //  
    ////////////////////////////////////////////////
    
    /////////////////////////////////    
    //
    //  Muon && Pion  Channel
    //
    /////////////////////////////////
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson3_idx,TauDecay::JAK_MUON,tau3_idx) && Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson4_idx,TauDecay::JAK_PION,tau4_idx)) { //JAK_MUON && JAK_PION
      if(verbose)std::cout<< "pi-mu " <<std::endl;
      
      TLorentzVector GenZ=Ntp->MCSignalParticle_p4(Boson3_idx);
      TLorentzVector GenTauMu(0,0,0,0);       TLorentzVector GenMu(0,0,0,0);  
      TLorentzVector GenNum(0,0,0,0);         TLorentzVector GenNutm(0,0,0,0);
      TLorentzVector Muon(0,0,0,0);           TLorentzVector TauMuRec(0,0,0,0);
      TLorentzVector MuonRec(0,0,0,0);        TLorentzVector Muon_mdR(0,0,0,0);
      TLorentzVector GenNumtm(0,0,0,0);

      TLorentzVector GenZ1=Ntp->MCSignalParticle_p4(Boson4_idx);
      TLorentzVector GenTauPi(0,0,0,0);      TLorentzVector GenPi(0,0,0,0);
      TLorentzVector GenNutp(0,0,0,0);       TLorentzVector Pion(0,0,0,0);
      TLorentzVector Zrec(0,0,0,0);          TLorentzVector TauPiRec(0,0,0,0);
      TLorentzVector PionRec(0,0,0,0);       TLorentzVector Pion_mdR(0,0,0,0);  
      
      double dR_Muons = -1;   double dR_Pions = -1;
      double mindRMu = 999;   double mindRPi=999; //dR_Muon < 0.1 && dR_Pion<0.4
      int mindRMu_Index = -1; int mindRPi_Index = -1;
      double Zrec_Mt = -1;    double Zgen_Mt = -1;

      double dphi_genTaus    =-1000; double dtheta_genTaus    =-1000;
      double dphi_genTaumu   =-1000; double dtheta_genTaumu   =-1000;
      double dphi_genTaupi   =-1000; double dtheta_genTaupi   =-1000;
      double dphi_genTauNum  =-1000; double dtheta_genTauNum  =-1000;
      double dphi_genTauNutm =-1000; double dtheta_genTauNutm =-1000;
      double dphi_genTauNutp =-1000; double dtheta_genTauNutp =-1000;
      
      double dphi_recTaus    =-1000; double dtheta_recTaus    =-1000;
      double dphi_recTaumu   =-1000; double dtheta_recTaumu   =-1000;
      double dphi_recTaupi   =-1000; double dtheta_recTaupi   =-1000;
      double dphi_recTauNum  =-1000; double dtheta_recTauNum  =-1000;
      double dphi_recTauNutm =-1000; double dtheta_recTauNutm =-1000;
      double dphi_recTauNutp =-1000; double dtheta_recTauNutp =-1000;

      double dR_m1 = -1; double dR_m2 = -1; double dR_m3 = -1;
      double dR_m4 = -1; double dR_m5 = -1; double dR_m6 = -1; 
      double dR_m7 = -1; double dR_m8 = -1; double dR_m9 = -1;
      double dR_m10 = -1; double dR_m11 = -1; double dR_m12 = -1;
      double dR_m13 = -1; double dR_m14 = -1; double dR_m15 = -1;
      double dR_m16 = -1; double dR_m17 = -1; double dR_m18 = -1;

      double dR_p1 = -1; double dR_p2 = -1; double dR_p3 = -1;
      double dR_p4 = -1; double dR_p5 = -1; double dR_p6 = -1;
      double dR_p7 = -1; double dR_p8 = -1; double dR_p9 = -1;
      double dR_p10 = -1; double dR_p11 = -1; double dR_p12 = -1;
      double dR_p13 = -1; double dR_p14 = -1; double dR_p15 = -1;
      double dR_p16 = -1; double dR_p17 = -1; double dR_p18 = -1;
      
      
      //Gen Infos      
      //Gen Muon//
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau3_idx);i++) {
        if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::tau_minus)) {
          GenTauMu=Ntp->MCTauandProd_p4(tau3_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::mu_plus)) {
          GenMu+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::nu_mu)) {
          GenNum+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::nu_tau)) {
          GenNutm+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
      } //end for loop 
      //Reco Muon
      for(unsigned int i=0; i < Ntp->NMuons(); i++) {
	dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu);
	dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
	if(dR_Muons < mindRMu) { mindRMu = dR_Muons;  mindRMu_Index = i; }
	MuonRec = Ntp->Muon_p4(mindRMu_Index);
      }
      if(mindRMu<0.1 && mindRMu_Index>-1) { 
	mindR_recMu_genMu.at(t).Fill(mindRMu,w);   
	std::cout<< "Moun_minDR==" << mindRMu <<std::endl;
	Muon_mdR=MuonRec;
      }
      //Gen Pion//
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau4_idx);i++) { 
        if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::tau_minus)) {
          GenTauPi=Ntp->MCTauandProd_p4(tau4_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::pi_plus) ) {
          GenPi+=Ntp->MCTauandProd_p4(tau4_idx,i);
	}
        else if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::nu_tau) ) {
	  GenNutp+=Ntp->MCTauandProd_p4(tau4_idx,i);
	}
      }      
      //Reco Pion
      for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
	dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi);
	dR_recPi_genPi.at(t).Fill(dR_Pions,w);
	if(dR_Pions < mindRPi) { mindRPi = dR_Pions;  mindRPi_Index = i; }
	PionRec = Ntp->PFTau_p4(mindRPi_Index);
      }
      if(mindRPi<0.4 && mindRPi_Index>-1) { 
	mindR_recPi_genPi.at(t).Fill(mindRPi,w);
	Pion_mdR=PionRec;
      }     
      
      // Transvers mass of Z boson
      // make sure that Spin_WT !=NaN
      if(Spin_WT > 0. &&  Spin_WT < 2.) {
	spin_WT.at(t).Fill(Spin_WT,w);
	
	/////////////////////////////////
	//// Tau_gen-XPlots ////
	/////////////////////////////////
	/////////////////////////////////
	//// Tau_gen  infos ////
	/////////////////////////////////                                                                                                                                                   
	if(GenTauMu.E()>0) {
	  M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);	  Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);
	  Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);	  Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w); }
	if(GenTauPi.E()>0) {
	  M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);	  Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);
	  Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w);	  Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w); }
	
	if(GenTauMu.E()>0 && GenTauPi.E()>0) {
	  dphi_genTaus      = GenTauPi.Phi()   - GenTauMu.Phi();        dPhi_genTauMu_genTauPi.at(t).Fill(dphi_genTaus,w);
	  dtheta_genTaus    = GenTauPi.Theta() - GenTauMu.Theta();      dTheta_genTauMu_genTauPi.at(t).Fill(dtheta_genTaus,w);
	  dphi_genTaumu     = GenTauMu.Phi()   - GenMu.Phi();           dPhi_genTauMu_genMu.at(t).Fill(dphi_genTaumu,w); 
	  dtheta_genTaumu   = GenTauMu.Theta() - GenMu.Theta();         dTheta_genTauMu_genMu.at(t).Fill(dtheta_genTaumu,w); 
	  dphi_genTaupi     = GenTauPi.Phi()   - GenPi.Phi();           dPhi_genTauPi_genPi.at(t).Fill(dphi_genTaupi,w);
	  dtheta_genTaupi   = GenTauPi.Theta() - GenPi.Theta();         dTheta_genTauPi_genPi.at(t).Fill(dtheta_genTaupi,w);
	  dphi_genTauNum    = GenTauMu.Phi()   - GenNum.Phi();          dPhi_genTauMu_genNum.at(t).Fill(dphi_genTauNum,w);
	  dtheta_genTauNum  = GenTauMu.Theta() - GenNum.Theta();        dTheta_genTauMu_genNum.at(t).Fill(dtheta_genTauNum,w);
	  dphi_genTauNutm   = GenTauMu.Phi()   - GenNutm.Phi();         dPhi_genTauMu_genNutm.at(t).Fill(dphi_genTauNutm,w);
	  dtheta_genTauNutm = GenTauMu.Theta() - GenNutm.Theta();       dTheta_genTauMu_genNutm.at(t).Fill(dtheta_genTauNutm,w);
	  dphi_genTauNutp   = GenTauPi.Phi()   - GenNutp.Phi();         dPhi_genTauPi_genNutp.at(t).Fill(dphi_genTauNutp,w);
	  dtheta_genTauNutp = GenTauPi.Theta() - GenNutp.Theta();       dTheta_genTauPi_genNutp.at(t).Fill(dtheta_genTauNutp,w);
	}
	  
	//////// Delta_R(vis,mis) ////////
	GenNumtm = GenNum + GenNutm; // sum of neutrinos.
	// dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu) ;
	// dR_recMu_genMu.at(t).Fill(dR_Muons,w);   

	if(10 <=GenTauMu.Pt() && GenTauMu.Pt()<= 15) {dR_m1=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_10_15.at(t).Fill(dR_m1,w);}
	if(15 <=GenTauMu.Pt() && GenTauMu.Pt()<= 20) {dR_m2=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_15_20.at(t).Fill(dR_m2,w);}
	if(20 <=GenTauMu.Pt() && GenTauMu.Pt()<= 25) {dR_m3=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_20_25.at(t).Fill(dR_m3,w);}
	if(25 <=GenTauMu.Pt() && GenTauMu.Pt()<= 30) {dR_m4=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_25_30.at(t).Fill(dR_m4,w);}
	if(30 <=GenTauMu.Pt() && GenTauMu.Pt()<= 35) {dR_m5=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_30_35.at(t).Fill(dR_m5,w);}
	if(35 <=GenTauMu.Pt() && GenTauMu.Pt()<= 40) {dR_m6=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_35_40.at(t).Fill(dR_m6,w);}
	if(40 <=GenTauMu.Pt() && GenTauMu.Pt()<= 45) {dR_m7=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_40_45.at(t).Fill(dR_m7,w);}
	if(45 <=GenTauMu.Pt() && GenTauMu.Pt()<= 50) {dR_m8=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_45_50.at(t).Fill(dR_m8,w);}
	if(50 <=GenTauMu.Pt() && GenTauMu.Pt()<= 55) {dR_m9=GenNumtm.DeltaR(GenMu);   dR_Mu_MisNum_50_55.at(t).Fill(dR_m9,w);}
	if(55 <=GenTauMu.Pt() && GenTauMu.Pt()<= 60) {dR_m10=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_55_60.at(t).Fill(dR_m10,w);}
	if(60 <=GenTauMu.Pt() && GenTauMu.Pt()<= 65) {dR_m11=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_60_65.at(t).Fill(dR_m11,w);}
	if(65 <=GenTauMu.Pt() && GenTauMu.Pt()<= 70) {dR_m12=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_65_70.at(t).Fill(dR_m12,w);}
	if(70 <=GenTauMu.Pt() && GenTauMu.Pt()<= 75) {dR_m13=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_70_75.at(t).Fill(dR_m13,w);}
	if(75 <=GenTauMu.Pt() && GenTauMu.Pt()<= 80) {dR_m14=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_75_80.at(t).Fill(dR_m14,w);}	  
	if(80 <=GenTauMu.Pt() && GenTauMu.Pt()<= 85) {dR_m15=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_80_85.at(t).Fill(dR_m15,w);}
	if(85 <=GenTauMu.Pt() && GenTauMu.Pt()<= 90) {dR_m16=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_85_90.at(t).Fill(dR_m16,w);}
	if(90 <=GenTauMu.Pt() && GenTauMu.Pt()<= 95) {dR_m17=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_90_95.at(t).Fill(dR_m17,w);}
	if(95 <=GenTauMu.Pt() && GenTauMu.Pt()<= 100){dR_m18=GenNumtm.DeltaR(GenMu);  dR_Mu_MisNum_95_100.at(t).Fill(dR_m18,w);}

	if(10 <=GenTauPi.Pt() && GenTauPi.Pt()<= 15) {dR_p1=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_10_15.at(t).Fill(dR_p1,w);}
	if(15 <=GenTauPi.Pt() && GenTauPi.Pt()<= 20) {dR_p2=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_15_20.at(t).Fill(dR_p2,w);}
	if(20 <=GenTauPi.Pt() && GenTauPi.Pt()<= 25) {dR_p3=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_20_25.at(t).Fill(dR_p3,w);}
	if(25 <=GenTauPi.Pt() && GenTauPi.Pt()<= 30) {dR_p4=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_25_30.at(t).Fill(dR_p4,w);}
	if(30 <=GenTauPi.Pt() && GenTauPi.Pt()<= 35) {dR_p5=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_30_35.at(t).Fill(dR_p5,w);}
	if(35 <=GenTauPi.Pt() && GenTauPi.Pt()<= 40) {dR_p6=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_35_40.at(t).Fill(dR_p6,w);}
	if(40 <=GenTauPi.Pt() && GenTauPi.Pt()<= 45) {dR_p7=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_40_45.at(t).Fill(dR_p7,w);}
	if(45 <=GenTauPi.Pt() && GenTauPi.Pt()<= 50) {dR_p8=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_45_50.at(t).Fill(dR_p8,w);}
	if(50 <=GenTauPi.Pt() && GenTauPi.Pt()<= 55) {dR_p9=GenNutp.DeltaR(GenPi);    dR_Pi_MisNup_50_55.at(t).Fill(dR_p9,w);}
	if(55 <=GenTauPi.Pt() && GenTauPi.Pt()<= 60) {dR_p10=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_55_60.at(t).Fill(dR_p10,w);}
	if(60 <=GenTauPi.Pt() && GenTauPi.Pt()<= 65) {dR_p11=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_60_65.at(t).Fill(dR_p11,w);}
	if(65 <=GenTauPi.Pt() && GenTauPi.Pt()<= 70) {dR_p12=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_65_70.at(t).Fill(dR_p12,w);}
	if(70 <=GenTauPi.Pt() && GenTauPi.Pt()<= 75) {dR_p13=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_70_75.at(t).Fill(dR_p13,w);}
	if(75 <=GenTauPi.Pt() && GenTauPi.Pt()<= 80) {dR_p14=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_75_80.at(t).Fill(dR_p14,w);}	  
	if(80 <=GenTauPi.Pt() && GenTauPi.Pt()<= 85) {dR_p15=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_80_85.at(t).Fill(dR_p15,w);}
	if(85 <=GenTauPi.Pt() && GenTauPi.Pt()<= 90) {dR_p16=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_85_90.at(t).Fill(dR_p16,w);}
	if(90 <=GenTauPi.Pt() && GenTauPi.Pt()<= 95) {dR_p17=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_90_95.at(t).Fill(dR_p17,w);}
	if(95 <=GenTauPi.Pt() && GenTauPi.Pt()<= 100){dR_p18=GenNutp.DeltaR(GenPi);   dR_Pi_MisNup_95_100.at(t).Fill(dR_p18,w);}

	  
	//Lab-Frame
	// X1 = Pi_Pt / Tau_Pt  
	std::cout<<" SpinWt=="<< Spin_WT<< std::endl;
	if(GenTauPi.Pt()>0 ) {
	  pi_PtRatio.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w);
	  pi_PtRatio_hplus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hplus);
	  pi_PtRatio_hminus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hminus);
	
	  pi_PtRatio_Spin.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*Spin_WT);
	  pi_PtRatio_hplus_Spin.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hplus*Spin_WT);
	  pi_PtRatio_hminus_Spin.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hminus*Spin_WT);
	}
	// X2 = Pi_Pt / Z_Mt 
	Zgen_Mt = sqrt(2*GenTauMu.Pt()*GenTauPi.Pt()*(1-cos(GenTauMu.Phi() - GenTauPi.Phi())));	
	std::cout<<"Z_Mt=="<< Zgen_Mt<<std::endl;
	//std::cout<<"Z_Mt1=="<< GenZ1.Mt()<<std::endl;
	if(Zgen_Mt>0) {
	  pi_PtoverZmt.at(t).Fill(GenPi.Pt()/Zgen_Mt,w);
	  pi_PtoverZmt_hplus.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hplus);
	  pi_PtoverZmt_hminus.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hminus); 
	
	  pi_PtoverZmt_Spin.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*Spin_WT);
	  pi_PtoverZmt_hplus_Spin.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hplus*Spin_WT);
	  pi_PtoverZmt_hminus_Spin.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hminus*Spin_WT);
	}
	// Z-Rest-Frame
	GenTauPi.Boost(-GenZ1.BoostVector());
	GenPi.Boost(-GenZ1.BoostVector());
	GenZ1.Boost(-GenZ1.BoostVector());
	// Now fill results                                   
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau4_idx));
	// X3 = Pi_E / Tau_E 
	if(GenTauPi.E()>0) {
	  pi_EoverEtau.at(t).Fill(GenPi.E()/GenTauPi.E(),w);
	  pi_EoverEtau_hplus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hplus);
	  pi_EoverEtau_hminus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hminus);
	  pi_EoverEtau_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*Spin_WT);
	  pi_EoverEtau_hplus_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hplus*Spin_WT);
	  pi_EoverEtau_hminus_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hminus*Spin_WT);
	}
	// X4 = Pi_E / Z_M
	std::cout<<"Z_M =="<< GenZ1.M()<<std::endl;
	if(GenZ1.M()>0) {
	  pi_EoverZm.at(t).Fill(2.*GenPi.E()/GenZ1.M(),w);
	  pi_EoverZm_hplus.at(t).Fill(2.*GenPi.E()/GenZ1.M(),w*hplus);
	  pi_EoverZm_hminus.at(t).Fill(2.*GenPi.E()/GenZ1.M(),w*hminus); 
	  pi_EoverZm_Spin.at(t).Fill(2.*GenPi.E()/GenZ1.M(),w*Spin_WT);
	  pi_EoverZm_hplus_Spin.at(t).Fill(2.*GenPi.E()/GenZ1.M(),w*hplus*Spin_WT);
	  pi_EoverZm_hminus_Spin.at(t).Fill(2.*GenPi.E()/GenZ1.M(),w*hminus*Spin_WT);    
	}
      
      
	/////////////////////////////////
	//// Tau_rec  reconstruction ////
	/////////////////////////////////
	if(mindRMu < 0.1 && mindRPi < 0.4 && mindRMu_Index >-1 && mindRPi_Index >-1) {
	  TauMuRec = Muon_mdR + GenNum + GenNutm;
	  TauPiRec = Pion_mdR + GenNutp;
	  Zrec = TauMuRec + TauPiRec;
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
	}
      
	if(TauMuRec.E()>0 && TauPiRec.E()>0) {
	
	  dphi_recTaus      = TauPiRec.Phi()   - TauMuRec.Phi();            dPhi_recTauMu_recTauPi.at(t).Fill(dphi_recTaus,w);
	  dtheta_recTaus    = TauPiRec.Theta() - TauMuRec.Theta();          dTheta_recTauMu_recTauPi.at(t).Fill(dtheta_recTaus,w);
	  dphi_recTaumu     = TauMuRec.Phi()   - Muon_mdR.Phi();            dPhi_recTauMu_recMu.at(t).Fill(dphi_recTaumu,w); 
	  dtheta_recTaumu   = TauMuRec.Theta() - Muon_mdR.Theta();          dTheta_recTauMu_recMu.at(t).Fill(dtheta_recTaumu,w); 
	  dphi_recTaupi     = TauPiRec.Phi()   - Pion_mdR.Phi();            dPhi_recTauPi_recPi.at(t).Fill(dphi_recTaupi,w);
	  dtheta_recTaupi   = TauPiRec.Theta() - Pion_mdR.Theta();          dTheta_recTauPi_recPi.at(t).Fill(dtheta_recTaupi,w);
	  dphi_recTauNum    = TauMuRec.Phi()   - GenNum.Phi();              dPhi_recTauMu_recNum.at(t).Fill(dphi_recTauNum,w);
	  dtheta_recTauNum  = TauMuRec.Theta() - GenNum.Theta();            dTheta_recTauMu_recNum.at(t).Fill(dtheta_recTauNum,w);
	  dphi_recTauNutm   = TauMuRec.Phi()   - GenNutm.Phi();             dPhi_recTauMu_recNutm.at(t).Fill(dphi_recTauNutm,w);
	  dtheta_recTauNutm = TauMuRec.Theta() - GenNutm.Theta();           dTheta_recTauMu_recNutm.at(t).Fill(dtheta_recTauNutm,w);
	  dphi_recTauNutp   = TauPiRec.Phi()   - GenNutp.Phi();             dPhi_recTauPi_recNutp.at(t).Fill(dphi_recTauNutp,w);
	  dtheta_recTauNutp = TauPiRec.Theta() - GenNutp.Theta();           dTheta_recTauPi_recNutp.at(t).Fill(dtheta_recTauNutp,w);
	}
      
        //Lab-Frame
	// X1 = Pirec_Pt / Tau_Pt  
	std::cout<<" SpinWt(rec)=="<< Spin_WT<< std::endl;
	if(TauPiRec.Pt()>0) {  
	  pirec_PtRatio.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w);
	  pirec_PtRatio_hplus.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hplus);
	  pirec_PtRatio_hminus.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hminus);
	  pirec_PtRatio_Spin.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*Spin_WT);
	  pirec_PtRatio_hplus_Spin.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hplus*Spin_WT);
	  pirec_PtRatio_hminus_Spin.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hminus*Spin_WT);
	}
	// X2 = Pirec_Pt / Z_Mt 
	Zrec_Mt = sqrt(2.*TauMuRec.Pt()*TauPiRec.Pt()*(1-cos(TauMuRec.Phi() - TauPiRec.Phi())));
	std::cout<<"Zrec_Mt=="<< Zrec_Mt<< std::endl;
	//std::cout<<"Zrec_Mt=="<< Zrec.Mt()<<std::endl;
	if(Zrec_Mt>0) {  
	  pirec_PtoverZmt.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w);
	  pirec_PtoverZmt_hplus.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hplus);
	  pirec_PtoverZmt_hminus.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hminus); 
	  pirec_PtoverZmt_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*Spin_WT);
	  pirec_PtoverZmt_hplus_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hplus*Spin_WT);
	  pirec_PtoverZmt_hminus_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hminus*Spin_WT);
	}
	// Z-Rest-Frame
	TauPiRec.Boost(-Zrec.BoostVector());
	Pion_mdR.Boost(-Zrec.BoostVector());
	Zrec.Boost(-Zrec.BoostVector());
	// Now fill results                                   
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau4_idx));
	// X3 = Pirec_E / Tau_E 
	if(TauPiRec.E()>0) {  
	  pirec_EoverEtau.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w);
	  pirec_EoverEtau_hplus.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hplus);
	  pirec_EoverEtau_hminus.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hminus);
	  pirec_EoverEtau_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*Spin_WT);
	  pirec_EoverEtau_hplus_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hplus*Spin_WT);
	  pirec_EoverEtau_hminus_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hminus*Spin_WT);
	}
	// X4 = Pirec_E / Z_M
	//double Z1rec_M = TauMuRec.M() + TauPiRec.M();
	//std::cout<<"Z1rec_M =="<< Z1rec_M<<std::endl;
	std::cout<<"Zrec_M =="<< Zrec.M()<< std::endl;
	if(Zrec.M()>0) {   
	  pirec_EoverZm.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w);
	  pirec_EoverZm_hplus.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hplus);
	  pirec_EoverZm_hminus.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hminus); 
	  pirec_EoverZm_Spin.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*Spin_WT);
	  pirec_EoverZm_hplus_Spin.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hplus*Spin_WT);
	  pirec_EoverZm_hminus_Spin.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hminus*Spin_WT);    
	}
      } //end if( 0.< Spin_WT <2.)  
      
    } //end JAK_MUON && JAK_PION   ////////////   end Muon && Pion    ////////////    
  } //end if(Status)   
} //end doEvent()


void  PolarizedTau::Finish() {
  
  Selection::Finish();
}







     






















