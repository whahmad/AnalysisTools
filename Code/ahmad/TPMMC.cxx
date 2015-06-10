#include "TPMMC.h"  
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
#include <TH2D.h>
#include <fstream>
#include <string>  
#include <sstream>  
#include <iomanip>

using namespace TMath;
using namespace std;


bool Isnan(double var)
{
  volatile double d = var;
  return d != d;
}  

bool NaN(double l)
{
  volatile double b = l;
  return b != b;
}  
                 
TPMMC::TPMMC(TString Name_, TString id_):
  Selection(Name_,id_)
{   
  //verbose=true;
}     

TPMMC::~TPMMC() {
  for(unsigned int j=0; j<Npassed.size(); j++) {
    std::cout << "TPMMC::~TPMMC Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TPMMC::~TPMMC()" << std::endl;
}
      
void  TPMMC::Configure() {
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
    
  // Setup Extra Histograms
  // Reco Objects
  //Delta_R(Mu_rec,Mu_gen), Delta_R(Pi_rec,Pi_gen)
  dR_recMu_genMu=HConfig.GetTH1D(Name+"_dR_recMu_genMu","dR_recMu_genMu",100,0.,8.,"dR_recMu_genMu","Events");
  mindR_recMu_genMu=HConfig.GetTH1D(Name+"_mindR_recMu_genMu","mindR_recMu_genMu",100,0.,1.,"mindR_recMu_genMu","Events");
  dR_recPi_genPi=HConfig.GetTH1D(Name+"_dR_recPi_genPi","dR_recPi_genPi",100,0.,8.,"dR_recPi_genPi","Events");
  mindR_recPi_genPi=HConfig.GetTH1D(Name+"_mindR_recPi_genPi","mindR_recPi_genPi",100,0.,1.,"mindR_recPi_genPi","Events");

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

  Mz_gen=HConfig.GetTH1D(Name+"_Mz_gen", "Mz_gen", 150.,0.,150.,"Mz_{gen}","Events");
  Mz_gen1=HConfig.GetTH1D(Name+"_Mz_gen1", "Mz_gen1", 150.,0.,150.,"Mz_{gen1}","Events");
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
  mtautau11_gen=HConfig.GetTH1D(Name+"_ditau11_gen", "ditau11_gen", 150.,0.,150.,"M_tautau_11_gen","Events");
  mtautau12_gen=HConfig.GetTH1D(Name+"_ditau12_gen", "ditau12_gen", 150.,0.,150.,"M_tautau_12_gen","Events");
  mtautau21_gen=HConfig.GetTH1D(Name+"_ditau21_gen", "ditau21_gen", 150.,0.,150.,"M_tautau_21_gen","Events");
  mtautau22_gen=HConfig.GetTH1D(Name+"_ditau22_gen", "ditau22_gen", 150.,0.,150.,"M_tautau_22_gen","Events");
  mtautau3_gen=HConfig.GetTH1D(Name+"_ditau3_gen", "ditau3_gen", 150.,0.,150.,"M_tautau_3_gen","Events");
  mtautau4_gen=HConfig.GetTH1D(Name+"_ditau4_gen", "ditau4_gen", 150.,0.,150.,"M_tautau_4_gen","Events");
  mtautau5_gen=HConfig.GetTH1D(Name+"-ditau5_gen", "ditau5_gen", 150.,0.,150.,"M_tautau_5_gen","Events");

  Lw_gen=HConfig.GetTH1D(Name+"-Lw_gen", "Lw_gen", 150.,-5.,2.,"Likelihood_weight_gen","Events");
  mtautau_gen=HConfig.GetTH1D(Name+"-ditau_gen", "ditau_gen", 150.,0.,150.,"M_tautau_Final_gen","Events");

  // Di-Tau mass
  mtautau11=HConfig.GetTH1D(Name+"_ditau11", "ditau11", 150.,0.,150.,"M_tautau_11","Events");
  mtautau12=HConfig.GetTH1D(Name+"_ditau12", "ditau12", 150.,0.,150.,"M_tautau_12","Events");
  mtautau21=HConfig.GetTH1D(Name+"_ditau21", "ditau21", 150.,0.,150.,"M_tautau_21","Events");
  mtautau22=HConfig.GetTH1D(Name+"_ditau22", "ditau22", 150.,0.,150.,"M_tautau_22","Events");
  mtautau3=HConfig.GetTH1D(Name+"_ditau3", "ditau3", 150.,0.,150.,"M_tautau_3","Events");
  mtautau4=HConfig.GetTH1D(Name+"_ditau4", "ditau4", 150.,0.,150.,"M_tautau_4","Events");
  mtautau5=HConfig.GetTH1D(Name+"-ditau5", "ditau5", 150.,0.,150.,"M_tautau_5","Events");

  Lw=HConfig.GetTH1D(Name+"-Lw", "Lw", 150.,-5.,2.,"Likelihood_weight","Events");
  mtautau=HConfig.GetTH1D(Name+"-ditau", "ditau", 150.,0.,150.,"M_tautau_Final","Events");

  Mzgenmmc_gen=HConfig.GetTH2D(Name+"-Mzgenmmc_gen", "Mzgenmmc_gen", 150.,0.,150.,150.,0.,150.);
  Mzrecmmc_gen=HConfig.GetTH2D(Name+"-Mzrecmmc_gen", "Mzrecmmc_gen", 150.,0.,150.,150.,0.,150.);
  Mzrecmmc_genmmc=HConfig.GetTH2D(Name+"-Mzrecmmc_genmmc", "Mzrecmmc_genmmc", 150.,0.,150.,150.,0.,150.);

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  TPMMC::Store_ExtraDist() {
  ////////////////////////////////////// 

  Extradist1d.push_back(&spin_WT);
  //Gen objects
  Extradist1d.push_back(&M_TauPi_gen);        Extradist1d.push_back(&Pt_TauPi_gen);
  Extradist1d.push_back(&Eta_TauPi_gen);      Extradist1d.push_back(&Phi_TauPi_gen);
  Extradist1d.push_back(&M_TauMu_gen);        Extradist1d.push_back(&Pt_TauMu_gen);
  Extradist1d.push_back(&Eta_TauMu_gen);      Extradist1d.push_back(&Phi_TauMu_gen);        
  Extradist1d.push_back(&Mz_gen);             Extradist1d.push_back(&Mz_gen1);
  //Reco objects
  Extradist1d.push_back(&M_TauPi_rec);        Extradist1d.push_back(&Pt_TauPi_rec);
  Extradist1d.push_back(&Eta_TauPi_rec);      Extradist1d.push_back(&Phi_TauPi_rec);
  Extradist1d.push_back(&M_TauMu_rec);        Extradist1d.push_back(&Pt_TauMu_rec);
  Extradist1d.push_back(&Eta_TauMu_rec);      Extradist1d.push_back(&Phi_TauMu_rec);        
  // Delta_R (,)
  Extradist1d.push_back(&dR_recMu_genMu);     Extradist1d.push_back(&dR_recPi_genPi);
  Extradist1d.push_back(&mindR_recMu_genMu);  Extradist1d.push_back(&mindR_recPi_genPi);
  // MET infos
  Extradist1d.push_back(&MetX);               Extradist1d.push_back(&MetY);
  Extradist1d.push_back(&MetX_gen);           Extradist1d.push_back(&MetY_gen);

  // Di-Tau mass
  Extradist1d.push_back(&mtautau11_gen);      Extradist1d.push_back(&mtautau12_gen);
  Extradist1d.push_back(&mtautau21_gen);      Extradist1d.push_back(&mtautau22_gen);
  Extradist1d.push_back(&mtautau3_gen);       Extradist1d.push_back(&mtautau4_gen); 
  Extradist1d.push_back(&mtautau5_gen);       Extradist1d.push_back(&mtautau_gen); 
  Extradist1d.push_back(&Lw_gen); 

  // Di-Tau mass
  Extradist1d.push_back(&mtautau11);          Extradist1d.push_back(&mtautau12);
  Extradist1d.push_back(&mtautau21);          Extradist1d.push_back(&mtautau22);
  Extradist1d.push_back(&mtautau3);           Extradist1d.push_back(&mtautau4); 
  Extradist1d.push_back(&mtautau5);           Extradist1d.push_back(&mtautau); 
  Extradist1d.push_back(&Lw); 

  Extradist2d.push_back(&Mzgenmmc_gen);
  Extradist2d.push_back(&Mzrecmmc_gen);
  Extradist2d.push_back(&Mzrecmmc_genmmc);
  
  
}

TF1  TPMMC::DeltaR_had(double pt) {
  TF1 cons2("cons2","pol1",1,66.88);      // landau 1st parameter
  cons2.SetParameter(0,-0.04697871);      cons2.SetParameter(1,0.03170594);
  TF1 MPV("MPV","expo+pol1(2)",1,92.8);   // landau 1st parameter
  MPV.SetParameter(0,-0.5269114);         MPV.SetParameter(1,-0.09367247);
  MPV.SetParameter(2,0.112788);           MPV.SetParameter(3,-0.0008607203);
  TF1 sigma2("sigma2","expo+pol1(2)",10.72,77.68);   // landau 1st parameter
  sigma2.SetParameter(0,-2.376518);       sigma2.SetParameter(1,-0.1253568);
  sigma2.SetParameter(2,0.00586322);      sigma2.SetParameter(3,-3.839789e-005);
  TF1 total("DeltaR","gaus(0)+landau(3)",0,0.4);
  double par[6];          //double *err;
  par[0]=0; par[1]=0; par[2]=0; par[3]=cons2.Eval(pt); par[4]=MPV.Eval(pt); par[5]=sigma2.Eval(pt);
  total.SetParameters(par);
  return total; }
TF1  TPMMC::DeltaR_lep(double pt) {
  TF1 cons1 ("cons1","pol1",1,109); 
  cons1.SetParameter(0,-0.0004049933);	 cons1.SetParameter(1,0.001609134);// gaus 1st parameter
  TF1 mean("mean","expo(0)+pol1(2)",10.72,85.24);// gaus 2nd parameter
  mean.SetParameter(0,-1.319604);        mean.SetParameter(1,-0.0698018);
  mean.SetParameter(2,0.05926357);       mean.SetParameter(3,-0.0004089469);
  TF1 sigma1("sigma1","expo+pol1(2)",10.72,109);// gaus 3rd parameter   
  sigma1.SetParameter(0,-2.227225);	 sigma1.SetParameter(1,-0.04167413);
  sigma1.SetParameter(2,6.679525e-005);  sigma1.SetParameter(3,0.0001051946);
  TF1 cons2 ("cons2","pol1",1,109);// landau 1st parameter
  cons2.SetParameter(0,-0.03423635);	 cons2.SetParameter(1,0.008789224);
  TF1 MPV("MPV","expo+pol1(2)",10.72,96.04); // landau 2nd parameter
  MPV.SetParameter(0,-0.8407024);        MPV.SetParameter(1,-0.06564579);
  MPV.SetParameter(2,0.07128014);	 MPV.SetParameter(3,-0.0004138105);
  TF1 sigma2("sigma2","expo+pol1(2)",11.8,92.8);// landau 3rd parameter
  sigma2.SetParameter(0,-2.364371);      sigma2.SetParameter(1,-0.09803685);
  sigma2.SetParameter(2,0.01046975);     sigma2.SetParameter(3,-8.072633e-005);
  TF1 total("DeltaR","gaus(0)+landau(3)",0,1);
  Double_t par[6];             
  par[0]=cons1.Eval(pt); par[1]=mean.Eval(pt); par[2]=sigma1.Eval(pt);
  par[3]=cons2.Eval(pt); par[4]=MPV.Eval(pt);  par[5]=sigma2.Eval(pt);      
  total.SetParameters(par);
  return total;  }

void  TPMMC::doEvent() {
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
      TLorentzVector Z_gen(0,0,0,0);         TLorentzVector Z_gen1(0,0,0,0);  
      double mzgen= -100;                    double mzgen1= -100;

      double dR_Muons = -1;   double dR_Pions = -1;   double mindRMu = 999; double mindRPi=999; //dR_Muon < 0.1 && dR_Pion<0.4
      int mindRMu_Index = -1; int mindRPi_Index = -1; double Zrec_Mt = -1;  double Zgen_Mt = -1;
  
      double m_tau= 1.777;     
      TLorentzVector p4miss11_gen(0,0,0,0),p4miss12_gen(0,0,0,0),p4miss21_gen(0,0,0,0),p4miss22_gen(0,0,0,0),p4vis1_gen(0,0,0,0),p4vis2_gen(0,0,0,0);
      TLorentzVector p41_gen(0,0,0,0),p42_gen(0,0,0,0), p43_gen(0,0,0,0),miss(0,0,0,0);
      double METX_gen= -1000;     double METY_gen= -1000;         double pt_miss1_gen= -1000;   double pt_miss2_gen= -1000; 
      double p_vis1_gen= -1000;   double theta_vis1_gen= -1000;   double phi_vis1_gen= -1000;   double m_vis1_gen= -1000; 
      double p_vis2_gen= -1000;   double theta_vis2_gen= -1000;   double phi_vis2_gen= -1000;   double m_vis2_gen= -1000; 
      // double p_miss1_gen= -1000;  double theta_miss1_gen= -1000;  double phi_miss1_gen= -1000;  double m_miss1_gen= 0;
      // double p_miss2_gen= -1000;  double theta_miss2_gen= -1000;  double phi_miss2_gen= -1000;  double m_miss2_gen= -1000; 
      double pz11_gen= -1000;     double pz12_gen= -1000;         double pz21_gen= -1000;       double pz22_gen= -1000;
      double W3_gen = 1.; double W31_gen = 1.; double W32_gen = 1.; double W4_gen = 1.; double W41_gen = 1.; double W42_gen = 1.; double W5_gen = 1.; //weight

      TLorentzVector p4miss11(0,0,0,0),p4miss12(0,0,0,0),p4miss21(0,0,0,0),p4miss22(0,0,0,0),p4vis1(0,0,0,0),p4vis2(0,0,0,0);
      TLorentzVector p41(0,0,0,0),p42(0,0,0,0), p43(0,0,0,0);
      double METX= -1000;         double METY= -1000;             double pt_miss1= -1000;       double pt_miss2= -1000;  
      double p_vis1= -1000;       double theta_vis1= -1000;       double phi_vis1= -1000;       double m_vis1= -1000; 
      double p_vis2= -1000;       double theta_vis2= -1000;       double phi_vis2= -1000;       double m_vis2= -1000; 
      double p_miss1= -1000;      double theta_miss1= -1000;      double phi_miss1= -1000;      double m_miss1= 0;
      double p_miss2= -1000;      double theta_miss2= -1000;      double phi_miss2= -1000;      double m_miss2= -1000; 
      double pz11= -1000;         double pz12= -1000;             double pz21= -1000;           double pz22= -1000;
      double W3 = 1.; double W31 = 1.; double W32 = 1.; double W4 = 1.; double W41 = 1.; double W42 = 1.; double W5 = 1.; //weight

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
      miss     = GenNutp + GenNum + GenNutm; // sum of neutrinos from both sides lep & had.
      Z_gen    = GenTauMu + GenTauPi;
      Z_gen1   = GenMu + GenPi + miss;


      /////////////////////////
      ////////   MET   ////////
      /////////////////////////

      // METX_gen= miss.Px(); 
      // METY_gen= miss.Py();
      METX_gen = miss.E()*sin(miss.Theta())*cos(miss.Phi());
      METY_gen = miss.E()*sin(miss.Theta())*sin(miss.Phi());


// enrec=sqrt(pow(pxnrec,2.)+pow(pynrec,2.)+pow(pznrec,2.));

      // METX_gen= miss.Px();
      // std::cout<<"METX_gen="<<METX_gen<<std::endl;      std::cout<<"METY_gen="<<METY_gen<<std::endl;

      //METX = Ntp->MET_Uncorr_ex(); //MET_CorrMVA_ex();      //METY = Ntp->MET_Uncorr_ey(); //MET_CorrMVA_ey();
      METX = Ntp->MET_CorrMVAMuTau_ex();
      METY = Ntp->MET_CorrMVAMuTau_ey();
      // std::cout<<"metx="<<METX<<std::endl;      std::cout<<"mety="<<METY<<std::endl;

      
      ////////////////////////////
      ////**** MMC infos **** ////
      ////////////////////////////
      //Gen
      p_vis1_gen= GenPi.P(); theta_vis1_gen= GenPi.Theta(); phi_vis1_gen= GenPi.Phi(); m_vis1_gen= GenPi.M(); 
      p_vis2_gen= GenMu.P(); theta_vis2_gen= GenMu.Theta(); phi_vis2_gen= GenMu.Phi(); m_vis2_gen= GenMu.M();  
      // Reco
      p_vis1= Pion_mdR.P(); theta_vis1= Pion_mdR.Theta(); phi_vis1= Pion_mdR.Phi(); m_vis1= Pion_mdR.M(); 
      p_vis2= Muon_mdR.P(); theta_vis2= Muon_mdR.Theta(); phi_vis2= Muon_mdR.Phi(); m_vis2= Muon_mdR.M();     
      
      // p_miss1= GenNutp.P(); theta_miss1= GenNutp.Theta(); phi_miss1= GenNutp.Phi(); m_miss1= GenNutp.M(); 
      // p_miss2= GenNumtm.P(); theta_miss2= GenNumtm.Theta(); phi_miss2= GenNumtm.Phi(); m_miss2= GenNumtm.M(); 
      
      ////**** end MMC infos **** ////
      
      // make sure that Spin_WT !=NaN
      if(Spin_WT > 0. &&  Spin_WT < 2.) {
	spin_WT.at(t).Fill(Spin_WT,w);

	//// MET-Plots ////
	/////////////////////////////////
	MetX.at(t).Fill(METX,w);	 MetY.at(t).Fill(METY,w);
	MetX_gen.at(t).Fill(METX_gen,w); MetY_gen.at(t).Fill(METY_gen,w);

	/////////////////////////////////
	////     Z_gen  infos      ////
	/////////////////////////////////                                                                                                                                                   
	mzgen= Z_gen.M();
	mzgen1= Z_gen1.M();
	if(Z_gen.E()>0) { Mz_gen.at(t).Fill(mzgen,w); }
	if(Z_gen1.E()>0) { Mz_gen1.at(t).Fill(mzgen1,w); }

	/////////////////////////////////
	////     Tau_gen  infos      ////
	/////////////////////////////////                                                                                                                                                   
	if(GenTauMu.E()>0) {
	  M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);	  Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);
	  Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);	  Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w); }
	if(GenTauPi.E()>0) {
	  M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);	  Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);
	  Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w);	  Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w); }

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

      	////***********************////
      	////      MMC-Method       ////
      	////***********************////
    	// double pt1 = (GenTauPi+GenNutp).Pt(); double pt2 = (GenTauMu+GenNumtm).Pt();
	// std::cout<<"Pt_had=="<<pt1<<std::endl; std::cout<<"Pt_lep=="<<pt2<<std::endl;
		
	////***********************////
      	//// Solving of 4 equtions ////
      	////***********************////
	TH1D mtt_est("mtt_est", "mtt_est2", 100.,0.,150.);
	TH1D mtt_est_gen("mtt_est_gen", "mtt_est2_gen", 100.,0.,150.);

	if(mindRMu < 0.1 && mindRPi < 0.4 && mindRMu_Index >-1 && mindRPi_Index >-1) {
	  
	  ////****  phace space (phi_miss1, phi_miss2, m_miss2)  ****////
	  m_miss1=0;//hadronic- sdie
	  for (double ii=-TMath::Pi();ii<TMath::Pi();ii+=0.1 ) {   
	    phi_miss1=ii;
	    for (double jj=TMath::Pi();jj>-TMath::Pi();jj-=0.1 ) {   
	      phi_miss2=jj;
	      for (double kk=0;kk<2.;kk+=0.1 ) {
		m_miss2=kk;//leptonic-side
		
		// Gen //
		// // calculate Pt form the first two equations.
		pt_miss1_gen=(-sin(phi_miss2)*METX_gen+cos(phi_miss2)*METY_gen)/sin(phi_miss1-phi_miss2); //gen_MET
		pt_miss2_gen=(sin(phi_miss1)*METX_gen-cos(phi_miss1)*METY_gen)/sin(phi_miss1-phi_miss2);  //gen_MEt
		
		pz11_gen=(1./(2*(pow(m_vis1_gen,2)+pow(sin(theta_vis1_gen),2)*pow(p_vis1_gen,2))))*(-cos(theta_vis1_gen)*(pow(m_miss1,2)+pow(m_vis1_gen,2)-pow(m_tau,2))*p_vis1_gen+cos(phi_miss1-phi_vis1_gen)*sin(2*theta_vis1_gen)*pow(p_vis1_gen,2)*pt_miss1_gen+sqrt((pow(m_vis1_gen,2)+pow(p_vis1_gen,2))*(pow(m_miss1,4)+pow(pow(m_vis1_gen,2)-pow(m_tau,2),2)+4*	cos(phi_miss1-phi_vis1_gen)*sin(theta_vis1_gen)*(-pow(m_vis1_gen,2)+pow(m_tau,2))*p_vis1_gen*pt_miss1_gen-4*(pow(m_vis1_gen,2)+pow(sin(theta_vis1_gen),2)*pow(sin(phi_miss1-phi_vis1_gen),2)*pow(p_vis1_gen,2))*pow(pt_miss1_gen,2)-2*pow(m_miss1,2)*(pow(m_vis1_gen,2)+pow(m_tau,2)+2*sin(theta_vis1_gen)*p_vis1_gen*(sin(theta_vis1_gen)*p_vis1_gen+cos(phi_miss1-phi_vis1_gen)*pt_miss1_gen)))));
		
		pz12_gen=(-1./(2*(pow(m_vis1_gen,2)+pow(sin(theta_vis1_gen),2)*pow(p_vis1_gen,2))))*(cos(theta_vis1_gen)*(pow(m_miss1,2)+pow(m_vis1_gen,2)-pow(m_tau,2))*p_vis1_gen-cos(phi_miss1-phi_vis1_gen)*sin(2*theta_vis1_gen)*pow(p_vis1_gen,2)*pt_miss1_gen+sqrt((pow(m_vis1_gen,2)+pow(p_vis1_gen,2))*(pow(m_miss1,4)+pow(pow(m_vis1_gen,2)-pow(m_tau,2),2)+4*	cos(phi_miss1-phi_vis1_gen)*sin(theta_vis1_gen)*(-pow(m_vis1_gen,2)+pow(m_tau,2))*p_vis1_gen*pt_miss1_gen-4*(pow(m_vis1_gen,2)+pow(sin(theta_vis1_gen),2)*pow(sin(phi_miss1-phi_vis1_gen),2)*pow(p_vis1_gen,2))*pow(pt_miss1_gen,2)-2*pow(m_miss1,2)*(pow(m_vis1_gen,2)+pow(m_tau,2)+2*sin(theta_vis1_gen)*p_vis1_gen*(sin(theta_vis1_gen)*p_vis1_gen+cos(phi_miss1-phi_vis1_gen)*pt_miss1_gen)))));
    
		pz21_gen=(1./(2*(pow(m_vis2_gen,2)+pow(sin(theta_vis2_gen),2)*pow(p_vis2_gen,2))))*(-cos(theta_vis2_gen)*(pow(m_miss2,2)+pow(m_vis2_gen,2)-pow(m_tau,2))*p_vis2_gen+cos(phi_miss2-phi_vis2_gen)*sin(2*theta_vis2_gen)*pow(p_vis2_gen,2)*pt_miss2_gen+sqrt((pow(m_vis2_gen,2)+pow(p_vis2_gen,2))*(pow(m_miss2,4)+pow(pow(m_vis2_gen,2)-pow(m_tau,2),2)+4*	cos(phi_miss2-phi_vis2_gen)*sin(theta_vis2_gen)*(-pow(m_vis2_gen,2)+pow(m_tau,2))*p_vis2_gen*pt_miss2_gen-4*(pow(m_vis2_gen,2)+pow(sin(theta_vis2_gen),2)*pow(sin(phi_miss2-phi_vis2_gen),2)*pow(p_vis2_gen,2))*pow(pt_miss2_gen,2)-2*pow(m_miss2,2)*(pow(m_vis2_gen,2)+pow(m_tau,2)+2*sin(theta_vis2_gen)*p_vis2_gen*(sin(theta_vis2_gen)*p_vis2_gen+cos(phi_miss2-phi_vis2_gen)*pt_miss2_gen)))));
    
		pz22_gen=(-1./(2*(pow(m_vis2_gen,2)+pow(sin(theta_vis2_gen),2)*pow(p_vis2_gen,2))))*(cos(theta_vis2_gen)*(pow(m_miss2,2)+pow(m_vis2_gen,2)-pow(m_tau,2))*p_vis2_gen-cos(phi_miss2-phi_vis2_gen)*sin(2*theta_vis2_gen)*pow(p_vis2_gen,2)*pt_miss2_gen+sqrt((pow(m_vis2_gen,2)+pow(p_vis2_gen,2))*(pow(m_miss2,4)+pow(pow(m_vis2_gen,2)-pow(m_tau,2),2)+4*	cos(phi_miss2-phi_vis2_gen)*sin(theta_vis2_gen)*(-pow(m_vis2_gen,2)+pow(m_tau,2))*p_vis2_gen*pt_miss2_gen-4*(pow(m_vis2_gen,2)+pow(sin(theta_vis2_gen),2)*pow(sin(phi_miss2-phi_vis2_gen),2)*pow(p_vis2_gen,2))*pow(pt_miss2_gen,2)-2*pow(m_miss2,2)*(pow(m_vis2_gen,2)+pow(m_tau,2)+2*sin(theta_vis2_gen)*p_vis2_gen*(sin(theta_vis2_gen)*p_vis2_gen+cos(phi_miss2-phi_vis2_gen)*pt_miss2_gen)))));
    
		////Be sure the solutions are not a nan. 
		if (!((Isnan(pz11_gen) && Isnan(pz12_gen)) || (Isnan(pz21_gen) && Isnan(pz22_gen) ))) {
		  
		  p4miss11_gen.SetXYZT(pt_miss1_gen*cos(phi_miss1),pt_miss1_gen*sin(phi_miss1),pz11,sqrt(pow(pz11,2)+pow(pt_miss1_gen,2)+pow(m_miss1,2)));
		  p4miss12_gen.SetXYZT(pt_miss1_gen*cos(phi_miss1),pt_miss1_gen*sin(phi_miss1),pz12,sqrt(pow(pz12,2)+pow(pt_miss1_gen,2)+pow(m_miss1,2)));
		  p4miss21_gen.SetXYZT(pt_miss2_gen*cos(phi_miss2),pt_miss2_gen*sin(phi_miss2),pz21,sqrt(pow(pz21,2)+pow(pt_miss2_gen,2)+pow(m_miss2,2)));
		  p4miss22_gen.SetXYZT(pt_miss2_gen*cos(phi_miss2),pt_miss2_gen*sin(phi_miss2),pz22,sqrt(pow(pz22,2)+pow(pt_miss2_gen,2)+pow(m_miss2,2)));
    
		  p4vis1_gen.SetXYZT(p_vis1_gen*sin(theta_vis1_gen)*cos(phi_vis1_gen),p_vis1_gen*sin(theta_vis1_gen)*sin(phi_vis1_gen), p_vis1_gen*cos(theta_vis1_gen), sqrt(pow(p_vis1_gen,2)+pow(m_vis1_gen,2)));
		  p4vis2_gen.SetXYZT(p_vis2_gen*sin(theta_vis2_gen)*cos(phi_vis2_gen),p_vis2_gen*sin(theta_vis2_gen)*sin(phi_vis2_gen), p_vis2_gen*cos(theta_vis2_gen), sqrt(pow(p_vis2_gen,2)+pow(m_vis2_gen,2)));
	      
		
		  // likelihood..... 
		  W31_gen= DeltaR_had((p4vis1_gen+p4miss11_gen).Pt()).Eval(p4vis1_gen.DeltaR(p4miss11_gen)) * DeltaR_lep((p4vis2_gen+p4miss21_gen).Pt()).Eval(p4vis2_gen.DeltaR(p4miss21_gen)); //probability ...as in paper
		  W32_gen= DeltaR_had((p4vis1_gen+p4miss11_gen).Pt()).Eval(p4vis1_gen.DeltaR(p4miss11_gen)) * DeltaR_lep((p4vis2_gen+p4miss22_gen).Pt()).Eval(p4vis2_gen.DeltaR(p4miss22_gen)); //probability ...as in paper
		  W41_gen= DeltaR_had((p4vis1_gen+p4miss12_gen).Pt()).Eval(p4vis1_gen.DeltaR(p4miss12_gen)) * DeltaR_lep((p4vis2_gen+p4miss21_gen).Pt()).Eval(p4vis2_gen.DeltaR(p4miss21_gen)); //probability ...as in paper
		  W42_gen= DeltaR_had((p4vis1_gen+p4miss12_gen).Pt()).Eval(p4vis1_gen.DeltaR(p4miss12_gen)) * DeltaR_lep((p4vis2_gen+p4miss22_gen).Pt()).Eval(p4vis2_gen.DeltaR(p4miss22_gen)); //probability ...as in paper

		
		  // //// *************************** ////
		  // //// select the correct solution //// 
		  // //// *************************** ////
		
		  // if(fabs((p4vis1_gen+p4miss11_gen + p4vis2_gen+p4miss21_gen).M() - 91.18) < fabs((p4vis1_gen+p4miss11_gen + p4vis2_gen+p4miss22_gen).M() - 91.18))
		  //   { p41_gen = p4vis1_gen+p4miss11_gen + p4vis2_gen+p4miss21_gen;
		  //     W3_gen = W31_gen; 
		  //     mtautau3_gen.at(t).Fill(p41_gen.M(),W3_gen);
		  //     //mtautau3_gen.at(t).Scale(1/mtautau3_gen.at(t).Integral()_gen);
		  //   }
		  // else {
		  //   p41_gen = p4vis1_gen+p4miss11_gen + p4vis2_gen+p4miss22_gen; 
		  //   W3_gen = W32_gen;
		  //   mtautau3_gen.at(t).Fill(p41_gen.M(),W3_gen);
		  //   //mtautau3_gen.at(t).Scale(1/mtautau3_gen.at(t).Integral()_gen);
		  // }
		
		  // //std::cout<<"diff(11-21)=="<<fabs((p4vis1_gen+p4miss11_gen + p4vis2_gen+p4miss21_gen).M() - 91.18)<<std::endl;
		  // //std::cout<<"diff(11-22)=="<<fabs((p4vis1_gen+p4miss11_gen + p4vis2_gen+p4miss22_gen).M() - 91.18)<<std::endl;     
		  
		  // if(fabs((p4vis1_gen+p4miss12_gen + p4vis2_gen+p4miss21_gen).M() - 91.18) < fabs((p4vis1_gen+p4miss12_gen + p4vis2_gen+p4miss22_gen).M() - 91.18))
		  //   { p42=p4vis1_gen+p4miss12_gen + p4vis2_gen+p4miss21_gen;
		  //     W4_gen = W41_gen;
		  //     mtautau4_gen.at(t).Fill(p42.M(),W4_gen);
		  //     //mtautau4_gen.at(t).Scale(1/mtautau4_gen.at(t).Integral()_gen);
		  //   }
		  // else { p42=p4vis1_gen+p4miss12_gen + p4vis2_gen+p4miss22_gen;
		  //   W4_gen = W42_gen;
		  //   mtautau4_gen.at(t).Fill(p42.M(),W4_gen); 
		  //   //mtautau4_gen.at(t).Scale(1/mtautau4_gen.at(t).Integral()_gen);
		  // }
		  
		  // //std::cout<<"diff(12-21)=="<<fabs((p4vis1_gen+p4miss12_gen + p4vis2_gen+p4miss21_gen).M() - 91.18)<<std::endl;
		  // //std::cout<<"diff(12-22)=="<<fabs((p4vis1_gen+p4miss12_gen + p4vis2_gen+p4miss22_gen).M() - 91.18)<<std::endl;
		  
		  // if(fabs(p41_gen.M() - 91.18) < fabs(p42.M() - 91.18))
		  //   { p43=p41_gen;
		  //     W5_gen = W3_gen;
		  //     mtautau5_gen.at(t).Fill(p43.M(),W5_gen);
		  //     mtt_est_gen.Fill(p43.M(),W5_gen);	    
		  //   }
		  // else { p43=p42;
		  //   W5_gen =W4_gen;
		  //   mtautau5_gen.at(t).Fill(p43.M(),W5_gen);
		  //   mtt_est_gen.Fill(p43.M(),W5_gen);
		  // }
		
		  Lw_gen.at(t).Fill(W5_gen);

		
		  mtautau11_gen.at(t).Fill((p4vis1_gen+p4miss11_gen +p4vis2_gen+ p4miss21_gen).M(),W31_gen);
		  mtt_est_gen.Fill((p4vis1_gen+p4miss11_gen +p4vis2_gen+ p4miss21_gen).M(),W31_gen);
		  mtautau12_gen.at(t).Fill((p4vis1_gen+p4miss11_gen +p4vis2_gen+ p4miss22_gen).M(),W32_gen);
		  mtt_est_gen.Fill((p4vis1_gen+p4miss11_gen +p4vis2_gen+ p4miss22_gen).M(),W32_gen);
		  mtautau21_gen.at(t).Fill((p4vis1_gen+p4miss12_gen +p4vis2_gen+ p4miss21_gen).M(),W41_gen);
		  mtt_est_gen.Fill((p4vis1_gen+p4miss12_gen +p4vis2_gen+ p4miss21_gen).M(),W41_gen);
		  mtautau22_gen.at(t).Fill((p4vis1_gen+p4miss12_gen +p4vis2_gen+ p4miss22_gen).M(),W42_gen);
		  mtt_est_gen.Fill((p4vis1_gen+p4miss12_gen +p4vis2_gen+ p4miss22_gen).M(),W42_gen);

		}
		  
		  // Rec //
		  pt_miss1=(-sin(phi_miss2)*METX+cos(phi_miss2)*METY)/sin(phi_miss1-phi_miss2); //rec_MET
		  pt_miss2=(sin(phi_miss1)*METX-cos(phi_miss1)*METY)/sin(phi_miss1-phi_miss2); //rec_MET
		
    
		  pz11=(1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))))*(-cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1+cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1+sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*	cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1)))));
    
		  pz12=(-1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))))*(cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1-cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1+sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*	cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1)))));
    
		  pz21=(1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))))*(-cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2+cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2+sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*	cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2)))));
    
		  pz22=(-1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))))*(cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2-cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2+sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*	cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2)))));
    
		  ////Be sure the solutions are not a nan. 
		  //if ((custom_isnan(pz11) && custom_isnan(pz12)) || (custom_isnan(pz21) && custom_isnan(pz22)  )) continue;
		  if (!((NaN(pz11) && NaN(pz12)) || (NaN(pz21) && NaN(pz22) ))) {
	      
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

		    // std::cout<<"W31=="<<W31<<std::endl;
		    // std::cout<<"W32=="<<W32<<std::endl;
		    // std::cout<<"W41=="<<W41<<std::endl;
		    // std::cout<<"W42=="<<W42<<std::endl;  
		
		    //		    std::cout<<"MMC_Check="<<std::endl;
		    //cout<<"//////////////////// Begin New Event/////////////////////"<<endl;
		    // std::cout<<"Root values     = x[0] "<<p4miss11.P()<<" "<<p4miss12.P()<<"   x[1] ="<<setw(13)<<p4miss11.Theta()<<" "<<p4miss12.Theta()<<"   x[2] ="
		    // 	     <<setw(13)<<p4miss21.P()<<" "<<p4miss22.P()<<"   x[3] ="<<setw(13)<<p4miss21.Theta()<<" "<<p4miss22.Theta()<<std::endl;
		    //	std::cout<<"Actual values   = x[0] ="<<setw(13)<<p_miss1<<"   x[1] ="<< setw(13)<<theta_miss1<<"   x[2] ="<<setw(13)<<p_miss2<<"   x[3] ="
		    //	 <<setw(13)<<theta_miss2<<std::endl;
		
		    // //// *************************** ////
		    // //// select the correct solution //// 
		    // //// *************************** ////
		
		    // if(fabs((p4vis1+p4miss11 + p4vis2+p4miss21).M() - 91.18) < fabs((p4vis1+p4miss11 + p4vis2+p4miss22).M() - 91.18))
		    //   { p41 = p4vis1+p4miss11 + p4vis2+p4miss21;
		    // 	W3 = W31; 
		    // 	mtautau3.at(t).Fill(p41.M(),W3);
		    // 	//mtautau3.at(t).Scale(1/mtautau3.at(t).Integral());
		    //   }
		    // else {
		    //   p41 = p4vis1+p4miss11 + p4vis2+p4miss22; 
		    //   W3 = W32;
		    //   mtautau3.at(t).Fill(p41.M(),W3);
		    //   //mtautau3.at(t).Scale(1/mtautau3.at(t).Integral());
		    // }
		
		    // //std::cout<<"diff(11-21)=="<<fabs((p4vis1+p4miss11 + p4vis2+p4miss21).M() - 91.18)<<std::endl;
		    // //std::cout<<"diff(11-22)=="<<fabs((p4vis1+p4miss11 + p4vis2+p4miss22).M() - 91.18)<<std::endl;     
		
		    // if(fabs((p4vis1+p4miss12 + p4vis2+p4miss21).M() - 91.18) < fabs((p4vis1+p4miss12 + p4vis2+p4miss22).M() - 91.18))
		    //   { p42=p4vis1+p4miss12 + p4vis2+p4miss21;
		    // 	W4 = W41;
		    // 	mtautau4.at(t).Fill(p42.M(),W4);
		    // 	//mtautau4.at(t).Scale(1/mtautau4.at(t).Integral());
		    //   }
		    // else { p42=p4vis1+p4miss12 + p4vis2+p4miss22;
		    //   W4 = W42;
		    //   mtautau4.at(t).Fill(p42.M(),W4); 
		    //   //mtautau4.at(t).Scale(1/mtautau4.at(t).Integral());
		    // }
		
		    // //std::cout<<"diff(12-21)=="<<fabs((p4vis1+p4miss12 + p4vis2+p4miss21).M() - 91.18)<<std::endl;
		    // //std::cout<<"diff(12-22)=="<<fabs((p4vis1+p4miss12 + p4vis2+p4miss22).M() - 91.18)<<std::endl;
		
		    // if(fabs(p41.M() - 91.18) < fabs(p42.M() - 91.18))
		    //   { p43=p41;
		    // 	W5 = W3;
		    // 	mtautau5.at(t).Fill(p43.M(),W5);
		    // 	//mtautau5.at(t).Scale(1/mtautau5.at(t).Integral());
		    // 	mtt_est.Fill(p43.M(),W5);	    
		    //   }
		    // else { p43=p42;
		    //   W5 =W4;
		    //   mtautau5.at(t).Fill(p43.M(),W5);
		    //   //mtautau5.at(t).Scale(1/mtautau5.at(t).Integral());
		    //   mtt_est.Fill(p43.M(),W5);
		    // }
		
		    Lw.at(t).Fill(W5);
		    // //std::cout<<"diff(1)=="<<fabs(p41.M() - 91.18)<<std::endl;
		    // //std::cout<<"diff(2)=="<<fabs(p42.M() - 91.18)<<std::endl;
		    // std::cout<<"W3=="<<W3<<std::endl;    std::cout<<"W4=="<<W4<<std::endl;    std::cout<<"W5=="<<W5<<std::endl;
		    // finish......
		
		    mtautau11.at(t).Fill((p4vis1+p4miss11 +p4vis2+ p4miss21).M(),W31);
		    //mtautau11.at(t).Scale(1/mtautau11.at(t).Integral());
		    mtt_est.Fill((p4vis1+p4miss11 +p4vis2+ p4miss21).M(),W31);
		    mtautau12.at(t).Fill((p4vis1+p4miss11 +p4vis2+ p4miss22).M(),W32);
		    //mtautau12.at(t).Scale(1/mtautau12.at(t).Integral());
		    mtt_est.Fill((p4vis1+p4miss11 +p4vis2+ p4miss22).M(),W32);
		    mtautau21.at(t).Fill((p4vis1+p4miss12 +p4vis2+ p4miss21).M(),W41);
		    //mtautau21.at(t).Scale(1/mtautau21.at(t).Integral());
		    mtt_est.Fill((p4vis1+p4miss12 +p4vis2+ p4miss21).M(),W41);
		    mtautau22.at(t).Fill((p4vis1+p4miss12 +p4vis2+ p4miss22).M(),W42);
		    //mtautau22.at(t).Scale(1/mtautau22.at(t).Integral());
		    mtt_est.Fill((p4vis1+p4miss12 +p4vis2+ p4miss22).M(),W42);
		
		  } //end if custom_isnan    

	      } // end for  kk
	    } // end for  jj
	  } // end for  ii
 
	    if(mtt_est_gen.Integral()>0) {	
	      int binmax_gen = mtt_est_gen.GetMaximumBin();
	      double x_gen = mtt_est_gen.GetBinCenter(binmax_gen);
	      mtautau_gen.at(t).Fill(x_gen);
	      //Mzgenmmc_gen.at(t).Fill(x_gen,mzgen);
	    }

	    if(mtt_est.Integral()>0) {	
	      int binmax = mtt_est.GetMaximumBin();
	      double x = mtt_est.GetBinCenter(binmax);
	      mtautau.at(t).Fill(x);
	      //Mzrecmmc_gen.at(t).Fill(x,mzgen);
	    }
      
	    if(mtt_est_gen.Integral()>0 && mtt_est.Integral()>0) {
	      int binmax_est_gen = mtt_est_gen.GetMaximumBin();
              double M_gen = mtt_est_gen.GetBinCenter(binmax_est_gen);
	      int binmax_rec = mtt_est.GetMaximumBin();
              double M_rec = mtt_est.GetBinCenter(binmax_rec);

	      Mzgenmmc_gen.at(t).Fill(M_gen,mzgen);
	      Mzrecmmc_gen.at(t).Fill(M_rec,mzgen); 
	      Mzrecmmc_genmmc.at(t).Fill(M_rec,M_gen);
	    }




	} // if(midRmu........	
	
      } //end if( 0.< Spin_WT <2.)  
    } //end JAK_MUON && JAK_PION   ////////////   end Muon && Pion    ////////////    
  } //end if(Status)
} //end doEvent()


void  TPMMC::Finish() {
  
  Selection::Finish();
}







     






















