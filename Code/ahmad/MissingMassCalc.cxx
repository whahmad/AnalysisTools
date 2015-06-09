#include "MissingMassCalc.h"  
#include "TLorentzVector.h"   
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>  
#include "TF1.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSpinerInterface.h"
 
#include <iostream>
#include "TError.h"
#include "TH1D.h"
#include <fstream>
#include <string>  
#include <sstream>  
#include <iomanip>

using namespace TMath;
using namespace std;

bool custom_isnan(double var)
{
  volatile double d = var;
  return d != d;
}  

                 
MissingMassCalc::MissingMassCalc(TString Name_, TString id_):
  Selection(Name_,id_)
{   
  //verbose=true;
}     

MissingMassCalc::~MissingMassCalc() {
  for(unsigned int j=0; j<Npassed.size(); j++) {
    std::cout << "MissingMassCalc::~MissingMassCalc Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "MissingMassCalc::~MissingMassCalc()" << std::endl;
}
      
void  MissingMassCalc::Configure() {
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
 
  // MET infos
  MetX=HConfig.GetTH1D(Name+"_MetX","MetX",30,-100.,120.,"E_{X}^{miss}","Events");
  MetY=HConfig.GetTH1D(Name+"_MetY","MetY",30,-100.,120.,"E_{Y}^{miss}","Events");
  MetX_gen=HConfig.GetTH1D(Name+"_MetX_gen","MetX_gen",30,-100.,120.,"E_{X_{gen}}^{miss}","Events");
  MetY_gen=HConfig.GetTH1D(Name+"_MetY_gen","MetY_gen",30,-100.,120.,"E_{Y_{gen}}^{miss}","Events");   

  // Di-Tau mass
  mtautau11=HConfig.GetTH1D(Name+"_ditau11", "ditau11", 100.,0.,150.,"M_tautau_11","Events");
  mtautau12=HConfig.GetTH1D(Name+"_ditau12", "ditau12", 100.,0.,150.,"M_tautau_12","Events");
  mtautau21=HConfig.GetTH1D(Name+"_ditau21", "ditau21", 100.,0.,150.,"M_tautau_21","Events");
  mtautau22=HConfig.GetTH1D(Name+"_ditau22", "ditau22", 100.,0.,150.,"M_tautau_22","Events");
  mtautau3=HConfig.GetTH1D(Name+"_ditau3", "ditau3", 100.,0.,150.,"M_tautau_3","Events");
  mtautau4=HConfig.GetTH1D(Name+"_ditau4", "ditau4", 100.,0.,150.,"M_tautau_4","Events");
  mtautau5=HConfig.GetTH1D(Name+"-ditau5", "ditau5", 100.,0.,150.,"M_tautau_5","Events");

  Lw=HConfig.GetTH1D(Name+"-Lw", "Lw", 100.,-5.,2.,"Likelihood_weight","Events");
  mtautau=HConfig.GetTH1D(Name+"-ditau", "ditau", 100.,0.,150.,"M_tautau_Final","Events");
  

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  MissingMassCalc::Store_ExtraDist() {
  ////////////////////////////////////// 

  Extradist1d.push_back(&spin_WT);
  // Gen_Infos /////// Lab_Frame ///////
  // X1 = Pi_Pt / Tau_Pt     
  Extradist1d.push_back(&pi_PtRatio);                Extradist1d.push_back(&pi_PtRatio_hplus);
  Extradist1d.push_back(&pi_PtRatio_hminus);         Extradist1d.push_back(&pi_PtRatio_Spin);
  Extradist1d.push_back(&pi_PtRatio_hplus_Spin);     Extradist1d.push_back(&pi_PtRatio_hminus_Spin);
  // X2 = Pi_Pt / Z_Mt
  Extradist1d.push_back(&pi_PtoverZmt);              Extradist1d.push_back(&pi_PtoverZmt_hplus);
  Extradist1d.push_back(&pi_PtoverZmt_hminus);       Extradist1d.push_back(&pi_PtoverZmt_Spin);
  Extradist1d.push_back(&pi_PtoverZmt_hplus_Spin);   Extradist1d.push_back(&pi_PtoverZmt_hminus_Spin);
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  Extradist1d.push_back(&pi_EoverEtau);              Extradist1d.push_back(&pi_EoverEtau_hplus);
  Extradist1d.push_back(&pi_EoverEtau_hminus);       Extradist1d.push_back(&pi_EoverEtau_Spin);
  Extradist1d.push_back(&pi_EoverEtau_hplus_Spin);   Extradist1d.push_back(&pi_EoverEtau_hminus_Spin);
  // X4 = Pi_E / Z_M
  Extradist1d.push_back(&pi_EoverZm);                Extradist1d.push_back(&pi_EoverZm_hplus);
  Extradist1d.push_back(&pi_EoverZm_hminus);         Extradist1d.push_back(&pi_EoverZm_Spin);
  Extradist1d.push_back(&pi_EoverZm_hplus_Spin);     Extradist1d.push_back(&pi_EoverZm_hminus_Spin);
  ////////////////////////////////////// 
  // Rec_Infos /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  Extradist1d.push_back(&pirec_PtRatio);             Extradist1d.push_back(&pirec_PtRatio_hplus);
  Extradist1d.push_back(&pirec_PtRatio_hminus);      Extradist1d.push_back(&pirec_PtRatio_Spin);
  Extradist1d.push_back(&pirec_PtRatio_hplus_Spin);  Extradist1d.push_back(&pirec_PtRatio_hminus_Spin);
  // X2 = Extradist1d.Push_Back(&Pirec_Pt / Zrec_Mt
  Extradist1d.push_back(&pirec_PtoverZmt);           Extradist1d.push_back(&pirec_PtoverZmt_hplus);
  Extradist1d.push_back(&pirec_PtoverZmt_hminus);    Extradist1d.push_back(&pirec_PtoverZmt_Spin);
  Extradist1d.push_back(&pirec_PtoverZmt_hplus_Spin);Extradist1d.push_back(&pirec_PtoverZmt_hminus_Spin);
  ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  Extradist1d.push_back(&pirec_EoverEtau);           Extradist1d.push_back(&pirec_EoverEtau_hplus);
  Extradist1d.push_back(&pirec_EoverEtau_hminus);    Extradist1d.push_back(&pirec_EoverEtau_Spin);
  Extradist1d.push_back(&pirec_EoverEtau_hplus_Spin);Extradist1d.push_back(&pirec_EoverEtau_hminus_Spin);
  // X4 = Extradist1d.Push_Back(&Pirec_E / Zrec_M
  Extradist1d.push_back(&pirec_EoverZm);             Extradist1d.push_back(&pirec_EoverZm_hplus);
  Extradist1d.push_back(&pirec_EoverZm_hminus);      Extradist1d.push_back(&pirec_EoverZm_Spin);
  Extradist1d.push_back(&pirec_EoverZm_hplus_Spin);  Extradist1d.push_back(&pirec_EoverZm_hminus_Spin);
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
  // MET infos
  Extradist1d.push_back(&MetX);               Extradist1d.push_back(&MetY);
  Extradist1d.push_back(&MetX_gen);           Extradist1d.push_back(&MetY_gen);
  // Di-Tau mass
  Extradist1d.push_back(&mtautau11);          Extradist1d.push_back(&mtautau12);
  Extradist1d.push_back(&mtautau21);          Extradist1d.push_back(&mtautau22);
  Extradist1d.push_back(&mtautau3);           Extradist1d.push_back(&mtautau4); 
  Extradist1d.push_back(&mtautau5);           Extradist1d.push_back(&mtautau); 
  Extradist1d.push_back(&Lw); 
  
}

TF1  MissingMassCalc::DeltaR_had(double pt) {
  
  TF1 cons2("cons2","pol1",1,66.88);                 // landau 1st parameter
  cons2.SetParameter(0,-0.04697871);    
  cons2.SetParameter(1,0.03170594);
  TF1 MPV("MPV","expo+pol1(2)",1,92.8);              // landau 1st parameter
  MPV.SetParameter(0,-0.5269114);
  MPV.SetParameter(1,-0.09367247);
  MPV.SetParameter(2,0.112788); 
  MPV.SetParameter(3,-0.0008607203);
  TF1 sigma2("sigma2","expo+pol1(2)",10.72,77.68);   // landau 1st parameter
  sigma2.SetParameter(0,-2.376518);
  sigma2.SetParameter(1,-0.1253568);
  sigma2.SetParameter(2,0.00586322); 
  sigma2.SetParameter(3,-3.839789e-005);
  
  TF1 total("DeltaR","gaus(0)+landau(3)",0,0.4);
  double par[6];          //double *err;
  par[0]=0;     	  par[1]=0;
  par[2]=0;     	  par[3]=cons2.Eval(pt);
  par[4]=MPV.Eval(pt);	  par[5]=sigma2.Eval(pt);
  
  total.SetParameters(par);
  //total.Draw();
  return total; 
}

TF1  MissingMassCalc::DeltaR_lep(double pt) {
  ////**** Gaus Fit ****////
  TF1 cons1 ("cons1","pol1",1,109); 
  cons1.SetParameter(0,-0.0004049933);	 cons1.SetParameter(1,0.001609134);// gaus 1st parameter
  TF1 mean("mean","expo(0)+pol1(2)",10.72,85.24);// gaus 2nd parameter
  mean.SetParameter(0,-1.319604);        mean.SetParameter(1,-0.0698018);
  mean.SetParameter(2,0.05926357);       mean.SetParameter(3,-0.0004089469);
  TF1 sigma1("sigma1","expo+pol1(2)",10.72,109);// gaus 3rd parameter   
  sigma1.SetParameter(0,-2.227225);	 sigma1.SetParameter(1,-0.04167413);
  sigma1.SetParameter(2,6.679525e-005);  sigma1.SetParameter(3,0.0001051946);
  ////**** End Gaus Fit ****////
  ////**** Landau Fit ****////
  TF1 cons2 ("cons2","pol1",1,109);// landau 1st parameter
  cons2.SetParameter(0,-0.03423635);	 cons2.SetParameter(1,0.008789224);
  TF1 MPV("MPV","expo+pol1(2)",10.72,96.04); // landau 2nd parameter
  MPV.SetParameter(0,-0.8407024);        MPV.SetParameter(1,-0.06564579);
  MPV.SetParameter(2,0.07128014);	 MPV.SetParameter(3,-0.0004138105);
  TF1 sigma2("sigma2","expo+pol1(2)",11.8,92.8);// landau 3rd parameter
  sigma2.SetParameter(0,-2.364371);      sigma2.SetParameter(1,-0.09803685);
  sigma2.SetParameter(2,0.01046975);     sigma2.SetParameter(3,-8.072633e-005);
  ////**** End Landau Fit ****////
  
  TF1 total("DeltaR","gaus(0)+landau(3)",0,1);
  Double_t par[6];             // Double_t *err;
  par[0]=cons1.Eval(pt);       // par[0] = 1st gaus parameter (norm factor)
  par[1]=mean.Eval(pt);        // par[1] = 2nd gaus parameter (mean)
  par[2]=sigma1.Eval(pt);      // par[2] = 3rd gaus parameter (sigma
  par[3]=cons2.Eval(pt);       // par[3] = 1st landau parameter
  par[4]=MPV.Eval(pt);         // par[4] = 2nd landau parameter
  par[5]=sigma2.Eval(pt);      // par[5] = 3rd landau parameter
  
  total.SetParameters(par);
  //total.Draw();
  return total;  
}


void  MissingMassCalc::doEvent() {
 
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
    //double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double Spin_WT=1;
    std::cout << "Spin_WT ==" << Spin_WT << std::endl;
    //double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    //double hminus=1-hplus;//Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
    double hplus=1;
    double hminus=1;
    std::cout << "hplus " << hplus << " hminus " << hminus << std::endl;
    
    spin_WT.at(t).Fill(Spin_WT,w);
  
    if(verbose)std::cout<< "A " <<std::endl;
    ////////////////////////////////////////////////
    // Spin Validation
    ////////////////////////////////////////////////
    
    /////////////////////////////////    
    //  Muon && Pion  Channel
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
      
      double METX= -1000;     double METY= -1000;
      double dR_Muons = -1;   double dR_Pions = -1;
      double mindRMu = 999;   double mindRPi=999; //dR_Muon < 0.1 && dR_Pion<0.4
      int mindRMu_Index = -1; int mindRPi_Index = -1;
      double Zrec_Mt = -1;    double Zgen_Mt = -1;

      double dR_m1 = -1;  double dR_m2 = -1;  double dR_m3 = -1;
      double dR_m4 = -1;  double dR_m5 = -1;  double dR_m6 = -1; 
      double dR_m7 = -1;  double dR_m8 = -1;  double dR_m9 = -1;
      double dR_m10 = -1; double dR_m11 = -1; double dR_m12 = -1;
      double dR_m13 = -1; double dR_m14 = -1; double dR_m15 = -1;
      double dR_m16 = -1; double dR_m17 = -1; double dR_m18 = -1;

      double dR_p1 = -1;  double dR_p2 = -1;  double dR_p3 = -1;
      double dR_p4 = -1;  double dR_p5 = -1;  double dR_p6 = -1;
      double dR_p7 = -1;  double dR_p8 = -1;  double dR_p9 = -1;
      double dR_p10 = -1; double dR_p11 = -1; double dR_p12 = -1;
      double dR_p13 = -1; double dR_p14 = -1; double dR_p15 = -1;
      double dR_p16 = -1; double dR_p17 = -1; double dR_p18 = -1;
      
      TLorentzVector p4miss11(0,0,0,0),p4miss12(0,0,0,0),p4miss21(0,0,0,0),p4miss22(0,0,0,0),p4vis1(0,0,0,0),p4vis2(0,0,0,0);
      TLorentzVector p41(0,0,0,0),p42(0,0,0,0), p43(0,0,0,0),miss(0,0,0,0);
      
      double m_tau= 1.777;     
      double p_vis1= -1000;    double theta_vis1= -1000;  double phi_vis1= -1000;  double m_vis1= -1000; 
      double p_vis2= -1000;    double theta_vis2= -1000;  double phi_vis2= -1000;  double m_vis2= -1000; 
      double p_miss1= -1000;   double theta_miss1= -1000; double phi_miss1= -1000; double m_miss1= 0;
      double p_miss2= -1000;   double theta_miss2= -1000; double phi_miss2= -1000; double m_miss2= -1000; 
      double pz11= -1000;      double pz12= -1000;        double pz21= -1000;      double pz22= -1000;
      double pt_miss1= -1000;  double pt_miss2= -1000;    double METX_gen= -1000;  double METY_gen= -1000; 
      double W3 = 1.;  double W31 = 1.;  double W32 = 1.; //weight
      double W4 = 1.;  double W41 = 1.;  double W42 = 1.; //weight
      double W5 = 1.;
      ////**** Gen- Reco Infos ****////        
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
      } 
      //Reco Muon
      for(unsigned int i=0; i < Ntp->NMuons(); i++) {
	dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu);
	dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
	if(dR_Muons < mindRMu) { mindRMu = dR_Muons;  mindRMu_Index = i; }
	MuonRec = Ntp->Muon_p4(mindRMu_Index);
      }
      if(mindRMu<0.1 && mindRMu_Index>-1) { 
	mindR_recMu_genMu.at(t).Fill(mindRMu,w);   
	//	std::cout<< "Moun_minDR==" << mindRMu <<std::endl;
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

      //Gen & Reco MET            
      GenNumtm = GenNum + GenNutm;       // sum of neutrinos from leptonic side.
      miss = GenNutp + GenNum + GenNutm; // sum of neutrinos from both sides lep & had.
      
      METX_gen= miss.Px(); 
      METY_gen= miss.Py(); 
      // std::cout<<"METX_gen="<<METX_gen<<std::endl;      std::cout<<"METY_gen="<<METY_gen<<std::endl;
      //METX = Ntp->MET_Uncorr_ex(); //MET_CorrMVA_ex();      //METY = Ntp->MET_Uncorr_ey(); //MET_CorrMVA_ey();
      METX = Ntp->MET_CorrMVAMuTau_ex();
      METY = Ntp->MET_CorrMVAMuTau_ey();
      // std::cout<<"metx="<<METX<<std::endl;      std::cout<<"mety="<<METY<<std::endl;

      ////////////////////////////
      ////**** MMC infos **** ////
      ////////////////////////////

      // // Reco
      // p_vis1= Pion_mdR.P(); theta_vis1= Pion_mdR.Theta(); phi_vis1= Pion_mdR.Phi(); m_vis1= Pion_mdR.M(); 
      // p_vis2= Muon_mdR.P(); theta_vis2= Muon_mdR.Theta(); phi_vis2= Muon_mdR.Phi(); m_vis2= Muon_mdR.M();     
       
      // //Gen
        p_vis1= GenPi.P(); theta_vis1= GenPi.Theta(); phi_vis1= GenPi.Phi(); m_vis1= GenPi.M(); 
        p_vis2= GenMu.P(); theta_vis2= GenMu.Theta(); phi_vis2= GenMu.Phi(); m_vis2= GenMu.M();  
      
      // p_miss1= GenNutp.P(); theta_miss1= GenNutp.Theta(); phi_miss1= GenNutp.Phi(); m_miss1= GenNutp.M(); 
      // p_miss2= GenNumtm.P(); theta_miss2= GenNumtm.Theta(); phi_miss2= GenNumtm.Phi(); m_miss2= GenNumtm.M(); 
      
      ////**** end MMC infos **** ////
      

      // make sure that Spin_WT !=NaN
	if(Spin_WT > 0. &&  Spin_WT < 2.) {
	spin_WT.at(t).Fill(Spin_WT,w);
	
	/////////////////////////////////
	//// MET-Plots ////
	/////////////////////////////////
	MetX.at(t).Fill(METX,w);
	MetY.at(t).Fill(METY,w);
	MetX_gen.at(t).Fill(METX_gen,w);
	MetY_gen.at(t).Fill(METY_gen,w);

	/////////////////////////////////
	//// Tau_gen  infos ////
	/////////////////////////////////                                                                                                                                                   
	if(GenTauMu.E()>0) {
	  M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);	  Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);
	  Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);	  Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w); }
	if(GenTauPi.E()>0) {
	  M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);	  Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);
	  Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w);	  Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w); }
	

	////********************************////
	////**** paramete space wheight ****////
	////********************************////
	//////// Delta_R(vis,mis) ////////
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
	// for(int i=10; i<105.; i+=5.) {}
	
 
	////************************************////	  
	////**** Spin-Observables Generated ****////
	////************************************////
	//Lab-Frame
	// X1 = Pi_Pt / Tau_Pt  
	//	std::cout<<" SpinWt=="<< Spin_WT<< std::endl;
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
	//	std::cout<<"Z_Mt=="<< Zgen_Mt<<std::endl;
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
	//	std::cout<<"Z_M =="<< GenZ1.M()<<std::endl;
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

     	////****************************************////
	////**** Spin-Observables Reconstructed ****////     
	////****************************************////
        //Lab-Frame
	// X1 = Pirec_Pt / Tau_Pt  
	//	std::cout<<" SpinWt(rec)=="<< Spin_WT<< std::endl;
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
	//	std::cout<<"Zrec_Mt=="<< Zrec_Mt<< std::endl;
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
	//	std::cout<<"Zrec_M =="<< Zrec.M()<< std::endl;
	if(Zrec.M()>0) {   
	  pirec_EoverZm.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w);
	  pirec_EoverZm_hplus.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hplus);
	  pirec_EoverZm_hminus.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hminus); 
	  pirec_EoverZm_Spin.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*Spin_WT);
	  pirec_EoverZm_hplus_Spin.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hplus*Spin_WT);
	  pirec_EoverZm_hminus_Spin.at(t).Fill(2.*Pion_mdR.E()/Zrec.M(),w*hminus*Spin_WT);    
	}


      	////***********************////
      	////      MMC-Method       ////
      	////***********************////
    	// double pt1 = (GenTauPi+GenNutp).Pt(); double pt2 = (GenTauMu+GenNumtm).Pt();
	// std::cout<<"Pt_had=="<<pt1<<std::endl; std::cout<<"Pt_lep=="<<pt2<<std::endl;
		
	////***********************////
      	//// Solving of 4 equtions ////
      	////***********************////
	TH1D mtt_est("mtt_est", "mtt_est2", 100.,0.,150.);

	//if(mindRMu < 0.1 && mindRPi < 0.4 && mindRMu_Index >-1 && mindRPi_Index >-1) {
	  
	  //phace space (phi_miss1, phi_miss2, m_miss2)
	  m_miss1=0;//hadronic- sdie
	  for (double ii=-TMath::Pi();ii<TMath::Pi();ii+=0.1 ) {
	    phi_miss1=ii;
	   
	    for (double jj=TMath::Pi();jj>-TMath::Pi();jj-=0.1 ) {
	      phi_miss2=jj;
	     
	      for (double kk=0;kk<2.;kk+=0.1 ) {
		m_miss2=kk;//leptonic-side
	       
	       
	      
		// // calculate Pt form the first two equations.
		pt_miss1=(-sin(phi_miss2)*METX_gen+cos(phi_miss2)*METY_gen)/sin(phi_miss1-phi_miss2); //gen_MET
		pt_miss2=(sin(phi_miss1)*METX_gen-cos(phi_miss1)*METY_gen)/sin(phi_miss1-phi_miss2);  //gen_MEt
		
		//	pt_miss1=(-sin(phi_miss2)*METX+cos(phi_miss2)*METY)/sin(phi_miss1-phi_miss2); //rec_MET
		//	pt_miss2=(sin(phi_miss1)*METX-cos(phi_miss1)*METY)/sin(phi_miss1-phi_miss2); //rec_MET
		
    
		pz11=(1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))))*(-cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1+cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1+sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*	cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1)))));
    
		pz12=(-1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))))*(cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1-cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1+sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*	cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1)))));
    
		pz21=(1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))))*(-cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2+cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2+sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*	cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2)))));
    
		pz22=(-1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))))*(cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2-cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2+sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*	cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2)))));
    
		////Be sure the solutions are not a nan. 
		//if ((custom_isnan(pz11) && custom_isnan(pz12)) || (custom_isnan(pz21) && custom_isnan(pz22)  )) continue;
		// if ((custom_isnan(pz11) && custom_isnan(pz12)) || (custom_isnan(pz21) && custom_isnan(pz22)  )) {
		if (!((custom_isnan(pz11) && custom_isnan(pz12)) || (custom_isnan(pz21) && custom_isnan(pz22) ))) {
	      
		  p4miss11.SetXYZT(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz11,sqrt(pow(pz11,2)+pow(pt_miss1,2)+pow(m_miss1,2)));
		  p4miss12.SetXYZT(pt_miss1*cos(phi_miss1),pt_miss1*sin(phi_miss1),pz12,sqrt(pow(pz12,2)+pow(pt_miss1,2)+pow(m_miss1,2)));
		  p4miss21.SetXYZT(pt_miss2*cos(phi_miss2),pt_miss2*sin(phi_miss2),pz21,sqrt(pow(pz21,2)+pow(pt_miss2,2)+pow(m_miss2,2)));
		  p4miss22.SetXYZT(pt_miss2*cos(phi_miss2),pt_miss2*sin(phi_miss2),pz22,sqrt(pow(pz22,2)+pow(pt_miss2,2)+pow(m_miss2,2)));
    
		  p4vis1.SetXYZT(p_vis1*sin(theta_vis1)*cos(phi_vis1),p_vis1*sin(theta_vis1)*sin(phi_vis1), p_vis1*cos(theta_vis1), sqrt(pow(p_vis1,2)+pow(m_vis1,2)));
		  p4vis2.SetXYZT(p_vis2*sin(theta_vis2)*cos(phi_vis2),p_vis2*sin(theta_vis2)*sin(phi_vis2), p_vis2*cos(theta_vis2), sqrt(pow(p_vis2,2)+pow(m_vis2,2)));
	      
		
		  // likelihood..... 
		  W31= DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11)) * DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)); //probability ...as in paper
		  W32= DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11)) * DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22)); //probability ...as in paper
		  W41= DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12)) * DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)); //probability ...as in paper
		  W42= DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12)) * DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22)); //probability ...as in paper

		  // W31= -log(DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11)) * DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21)));
		  // W32= -log(DeltaR_had((p4vis1+p4miss11).Pt()).Eval(p4vis1.DeltaR(p4miss11)) * DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22))); 
		  // W41= -log(DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12)) * DeltaR_lep((p4vis2+p4miss21).Pt()).Eval(p4vis2.DeltaR(p4miss21))); 
		  // W42= -log(DeltaR_had((p4vis1+p4miss12).Pt()).Eval(p4vis1.DeltaR(p4miss12)) * DeltaR_lep((p4vis2+p4miss22).Pt()).Eval(p4vis2.DeltaR(p4miss22))); 
		  
		  // std::cout<<"W31=="<<W31<<std::endl;
		  // std::cout<<"W32=="<<W32<<std::endl;
		  // std::cout<<"W41=="<<W41<<std::endl;
		  // std::cout<<"W42=="<<W42<<std::endl;  
		
		  std::cout<<"MMC_Check="<<std::endl;
		  //cout<<"//////////////////// Begin New Event/////////////////////"<<endl;
		  		  std::cout<<"Root values     = x[0] "<<p4miss11.P()<<" "<<p4miss12.P()<<"   x[1] ="<<setw(13)<<p4miss11.Theta()<<" "<<p4miss12.Theta()<<"   x[2] ="
		  	   <<setw(13)<<p4miss21.P()<<" "<<p4miss22.P()<<"   x[3] ="<<setw(13)<<p4miss21.Theta()<<" "<<p4miss22.Theta()<<std::endl;
		  //	std::cout<<"Actual values   = x[0] ="<<setw(13)<<p_miss1<<"   x[1] ="<< setw(13)<<theta_miss1<<"   x[2] ="<<setw(13)<<p_miss2<<"   x[3] ="
		  //	 <<setw(13)<<theta_miss2<<std::endl;
		
		  //// *************************** ////
		  //// select the correct solution //// 
		  //// *************************** ////
		
		  if(fabs((p4vis1+p4miss11 + p4vis2+p4miss21).M() - 91.18) < fabs((p4vis1+p4miss11 + p4vis2+p4miss22).M() - 91.18))
		    { p41 = p4vis1+p4miss11 + p4vis2+p4miss21;
		      W3 = W31; 
		      mtautau3.at(t).Fill(p41.M(),W3);
		      //mtautau3.at(t).Scale(1/mtautau3.at(t).Integral());
		    }
		  else {
		    p41 = p4vis1+p4miss11 + p4vis2+p4miss22; 
		    W3 = W32;
		    mtautau3.at(t).Fill(p41.M(),W3);
		    //mtautau3.at(t).Scale(1/mtautau3.at(t).Integral());
		  }
		
		  //std::cout<<"diff(11-21)=="<<fabs((p4vis1+p4miss11 + p4vis2+p4miss21).M() - 91.18)<<std::endl;
		  //std::cout<<"diff(11-22)=="<<fabs((p4vis1+p4miss11 + p4vis2+p4miss22).M() - 91.18)<<std::endl;     
		
		  if(fabs((p4vis1+p4miss12 + p4vis2+p4miss21).M() - 91.18) < fabs((p4vis1+p4miss12 + p4vis2+p4miss22).M() - 91.18))
		    { p42=p4vis1+p4miss12 + p4vis2+p4miss21;
		      W4 = W41;
		      mtautau4.at(t).Fill(p42.M(),W4);
		      //mtautau4.at(t).Scale(1/mtautau4.at(t).Integral());
		    }
		  else { p42=p4vis1+p4miss12 + p4vis2+p4miss22;
		    W4 = W42;
		    mtautau4.at(t).Fill(p42.M(),W4); 
		    //mtautau4.at(t).Scale(1/mtautau4.at(t).Integral());
		  }
		
		  //std::cout<<"diff(12-21)=="<<fabs((p4vis1+p4miss12 + p4vis2+p4miss21).M() - 91.18)<<std::endl;
		  //std::cout<<"diff(12-22)=="<<fabs((p4vis1+p4miss12 + p4vis2+p4miss22).M() - 91.18)<<std::endl;
		
		  if(fabs(p41.M() - 91.18) < fabs(p42.M() - 91.18))
		    { p43=p41;
		      W5 = W3;
		      mtautau5.at(t).Fill(p43.M(),W5);
		      //mtautau5.at(t).Scale(1/mtautau5.at(t).Integral());
		      mtt_est.Fill(p43.M(),W5);	    
		    }
		  else { p43=p42;
		    W5 =W4;
		    mtautau5.at(t).Fill(p43.M(),W5);
		    //mtautau5.at(t).Scale(1/mtautau5.at(t).Integral());
		    mtt_est.Fill(p43.M(),W5);
		  }
		
		  Lw.at(t).Fill(W5);

		  // //std::cout<<"diff(1)=="<<fabs(p41.M() - 91.18)<<std::endl;
		  // //std::cout<<"diff(2)=="<<fabs(p42.M() - 91.18)<<std::endl;
		  // std::cout<<"W3=="<<W3<<std::endl;
		  // std::cout<<"W4=="<<W4<<std::endl;
		  // std::cout<<"W5=="<<W5<<std::endl;
		
		  // finish......
		
		  mtautau11.at(t).Fill((p4vis1+p4miss11 +p4vis2+ p4miss21).M(),W31);
		  //mtautau11.at(t).Scale(1/mtautau11.at(t).Integral());
		  mtautau12.at(t).Fill((p4vis1+p4miss11 +p4vis2+ p4miss22).M(),W32);
		  //mtautau12.at(t).Scale(1/mtautau12.at(t).Integral());
		  mtautau21.at(t).Fill((p4vis1+p4miss12 +p4vis2+ p4miss21).M(),W41);
		  //mtautau21.at(t).Scale(1/mtautau21.at(t).Integral());
		  mtautau22.at(t).Fill((p4vis1+p4miss12 +p4vis2+ p4miss22).M(),W42);
		  //mtautau22.at(t).Scale(1/mtautau22.at(t).Integral());
		
		} //end if custom_isnan    
	      } // end for  kk
	    } // end for  jj
	  } // end for  ii
 
	   if(mtt_est.Integral()>0) {	
	     int binmax = mtt_est.GetMaximumBin();
	     double x = mtt_est.GetBinCenter(binmax);
	     mtautau.at(t).Fill(x);
	     //mtautau.at(t).Scale(1/mtautau.at(t).Integral());
	  

	     //t t =0;
	  //std::cout<<"t=="<<t<<std::endl;
	     
	  // for(int i= mtt_est.at(t).FindBin(-5. + mtt_est.at(t).GetBinCenter(mtt_est.at(t).GetMaximumBin())); i< mtt_est.at(t).FindBin(6. + mtt_est.at(t).GetBinCenter(mtt_est.at(t).GetMaximumBin())); i++) {
	  //   std::cout<<"i=="<<i<<std::endl;
	  //   t += mtt_est.at(t).GetBinCenter(i); 
	  //   std::cout<<"t=="<<t<<std::endl;
	  // }
	  
	  // if(mtt_est.at(t).Integral()>0) {	
	  //   int binmax = t/11.;//mtt_est.at(t).GetMaximumBin();
	  //   std::cout<<"Average=="<<binmax<<std::endl;
	    
	  //   double x = mtt_est.at(t).GetBinCenter(binmax);
	  //   mtautau.at(t).Fill(x);
	    
	    

	  // for(int i=-5. + mtt_est.at(t).GetBinCenter(mtt_est.at(t).GetMaximumBin()); i< 6. + mtt_est.at(t).GetBinCenter(mtt_est.at(t).GetMaximumBin()); i++) {
	  //   std::cout<<"i=="<<i<<std::endl;
	  //   t += i;  
	  //   std::cout<<"t=="<<t<<std::endl;
	  // }
	  
	  // if(mtt_est.at(t).Integral()>0) {	
	  //   int binmax = t/11.;//mtt_est.at(t).GetMaximumBin();
	  //   std::cout<<"Average=="<<binmax<<std::endl; 
	  //   double x = mtt_est.at(t).GetBinCenter(binmax);
	  //   mtautau.at(t).Fill(x);   
	    //mtt_est.at(t).Reset();
	   }	
	   //	}
      } //end if( 0.< Spin_WT <2.)  
    } //end JAK_MUON && JAK_PION   ////////////   end Muon && Pion    ////////////    
  } //end if(Status)
} //end doEvent()


void  MissingMassCalc::Finish() {
  
  Selection::Finish();
}







     






















