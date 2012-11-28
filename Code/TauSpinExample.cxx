#include "TauSpinExample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TF1.h"
#include "TauSolver.h"

TauSpinExample::TauSpinExample(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

TauSpinExample::~TauSpinExample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TauSpinExample::~TauSpinExample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TauSpinExample::~TauSpinExample()" << std::endl;
}

void  TauSpinExample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==isZtautauto3pimu)          cut.at(isZtautauto3pimu)=1;
  }

  TString hlabel;
  TString htitle; for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==isZtautauto3pimu){
      title.at(i)="Is $Z\\rightarrow\\tau\\tau\\rightarrow\\mu\\pi\\pi\\pi$ MC (bool)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Is Z#rightarrow#tau#tau#rightarrow#mu#pi#pi#pi MC (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isZtautauto3pimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isZtautauto3pimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
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

  //a1nu
  a1_ExoverEtau=HConfig.GetTH1D(Name+"_a1_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}","Events");
  a1_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_a1_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{h^{+}}","Events");
  a1_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_a1_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{h^{-}}","Events");
  a1_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_a1_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{Spin}","Events");
  a1_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_a1_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{UnSpin}","Events");
  a1_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_a1_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{FlipSpin}","Events");
  a1_WT_Spin=HConfig.GetTH1D(Name+"_a1_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{a_{1}(1260)}","Events");
  a1_WT_UnSpin=HConfig.GetTH1D(Name+"_a1_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{a_{1}(1260)}","Events");
  a1_WT_FlipSpin=HConfig.GetTH1D(Name+"_a1_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{a_{1}(1260)}","Events");
  a1_PtRatio_hplus=HConfig.GetTH1D(Name+"_a1_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,a_{1}(1260)}/P_{t,#tau}|_{h^{+}}","Events");
  a1_PtRatio_hminus=HConfig.GetTH1D(Name+"_a1_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,a_{1}(1260)}/P_{t,#tau}|_{h^{-}}","Events");

  a1_cosbeta_hminus=HConfig.GetTH1D(Name+"_a1_cosbeta_hminus","cosbeta_hminus",20,0.0,1.0,"cos#beta|_{h^{-}}","Events");
  a1_cosbeta_hplus=HConfig.GetTH1D(Name+"_a1_cosbeta_hplus","cosbeta_hplus",20,0.0,1.0,"cos#beta|_{h^{+}}","Events");

  //rhonu
  rho_ExoverEtau=HConfig.GetTH1D(Name+"_rho_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{#rho}/E_{#tau}","Events");
  rho_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_rho_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{h^{+}}","Events");
  rho_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_rho_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{h^{-}}","Events");
  rho_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_rho_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{Spin}","Events");
  rho_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_rho_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{UnSpin}","Events");
  rho_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_rho_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{FlipSpin}","Events");
  rho_WT_Spin=HConfig.GetTH1D(Name+"_rho_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#rho}","Events");
  rho_WT_UnSpin=HConfig.GetTH1D(Name+"_rho_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#rho}","Events");
  rho_WT_FlipSpin=HConfig.GetTH1D(Name+"_rho_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#rho}","Events");
  rho_PtRatio_hplus=HConfig.GetTH1D(Name+"_rho_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#rho}/P_{t,#tau}|_{h^{+}}","Events");
  rho_PtRatio_hminus=HConfig.GetTH1D(Name+"_rho_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#rho}/P_{t,#tau}|_{h^{-}}","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  TauSpinExample::Store_ExtraDist(){

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
 
 Extradist1d.push_back(&a1_ExoverEtau);
 Extradist1d.push_back(&a1_ExoverEtau_hplus);
 Extradist1d.push_back(&a1_ExoverEtau_hminus);
 Extradist1d.push_back(&a1_ExoverEtau_Spin);
 Extradist1d.push_back(&a1_ExoverEtau_UnSpin);
 Extradist1d.push_back(&a1_ExoverEtau_FlipSpin);
 Extradist1d.push_back(&a1_WT_Spin);
 Extradist1d.push_back(&a1_WT_UnSpin);
 Extradist1d.push_back(&a1_WT_FlipSpin);
 Extradist1d.push_back(&a1_PtRatio_hplus);
 Extradist1d.push_back(&a1_PtRatio_hminus);

 Extradist1d.push_back(&a1_cosbeta_hplus);
 Extradist1d.push_back(&a1_cosbeta_hminus);
 
 Extradist1d.push_back(&rho_ExoverEtau);
 Extradist1d.push_back(&rho_ExoverEtau_hplus);
 Extradist1d.push_back(&rho_ExoverEtau_hminus);
 Extradist1d.push_back(&rho_ExoverEtau_Spin);
 Extradist1d.push_back(&rho_ExoverEtau_UnSpin);
 Extradist1d.push_back(&rho_ExoverEtau_FlipSpin);
 Extradist1d.push_back(&rho_WT_Spin);
 Extradist1d.push_back(&rho_WT_UnSpin);
 Extradist1d.push_back(&rho_WT_FlipSpin);
 Extradist1d.push_back(&rho_PtRatio_hplus);
 Extradist1d.push_back(&rho_PtRatio_hminus);

}

void  TauSpinExample::doEvent(){
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){t=2;}// std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,tau_idx;
  value.at(isZtautauto3pimu)=0;
  pass.at(isZtautauto3pimu) = true;
  if(pass.at(isZtautauto3pimu))value.at(isZtautauto3pimu)=1;

  double wobs=1;
  double w=1;
  if(verbose)  std::cout << Ntp->GetMCID() << " " << Npassed.size() << " " << t << " " << Boson_idx << " " << " " << tau_idx 
	    << " " << Ntp->NMCTaus() << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;

    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
    double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
  
      ////////////////////////////////////////////////
      //
      // Spin Validation
      //
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_MUON,tau_idx)){
      if(verbose)std::cout << "muon" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::mu_plus)){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
        mu_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
        mu_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);

        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	mu_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	mu_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	mu_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	mu_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	mu_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
	mu_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);
	mu_WT_Spin.at(t).Fill(Spin_WT,w);
	mu_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	mu_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_PION,tau_idx)){
      if(verbose)std::cout << "pion" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
	pi_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
        pi_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);

        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	pi_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	pi_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	pi_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	pi_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	pi_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
	pi_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);
	pi_WT_Spin.at(t).Fill(Spin_WT,w);
	pi_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	pi_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_A1_3PI,tau_idx)){
      if(verbose)std::cout << "a1to3pi" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      std::vector<TLorentzVector> pions;
      std::vector<float> pions_charge;
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi0)
                ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
	  pions.push_back(Ntp->MCTauandProd_p4(tau_idx,i));
	  pions_charge.push_back(Ntp->MCTauandProd_charge(tau_idx,i));
        }
      }
      if(Tau_LV.E()>0){
        a1_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
        a1_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);

        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	a1_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	a1_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	a1_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	a1_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	a1_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
	a1_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);
	a1_WT_Spin.at(t).Fill(Spin_WT,w);
	a1_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	a1_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
	TauSolver TS(Tau_LV.Vect(),X_LV);
	float cosbeta, gamma;
	TS.EulerAnglesfor3prong(pions,pions_charge,cosbeta,gamma,true,true);
	a1_cosbeta_hplus.at(t).Fill(cosbeta,w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
	a1_cosbeta_hminus.at(t).Fill(cosbeta,w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);
      }
    }
    if(Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_RHO_PIPI0,tau_idx)){
      if(verbose)std::cout << "rho" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
	  if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)){
	    Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
	  }
	  else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi_plus) ||
		  abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::pi0)
		  ){
	    X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
	  }
      }
      if(Tau_LV.E()>0){
        rho_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
        rho_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);

	Tau_LV.Boost(-Boson_LV.BoostVector());
	X_LV.Boost(-Boson_LV.BoostVector());
	Boson_LV.Boost(-Boson_LV.BoostVector());
	// Now fill results                                                                                                                                                                                                                  
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	rho_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	rho_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	rho_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	rho_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	rho_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus)*UnSpin_WT);
	rho_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus)*UnSpin_WT);
	rho_WT_Spin.at(t).Fill(Spin_WT,w);
	rho_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	rho_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
  }
  if(verbose)std::cout << "done" << std::endl;
}




void  TauSpinExample::Finish(){
  unsigned int tdata=0;
  TF1 fmu_ExoverEtau_hplus( "fmu_ExoverEtau_hplus", "5.0/3.0-3.0*(x^2.0)+(4.0/3.0)*x^3.0-1*(-1/3+3*x^2-8/3*x^3)",0,1);
  TF1 fmu_ExoverEtau_hminus("fmu_ExoverEtau_hminus","5.0/3.0-3.0*(x^2.0)+(4.0/3.0)*x^3.0+1*(-1/3+3*x^2-8/3*x^3)",0,1);
  TF1 fpi_ExoverEtau_hplus( "fpi_ExoverEtau_hplus", "x",0,1);
  TF1 fpi_ExoverEtau_hminus("fpi_ExoverEtau_hminus","1-x",0,1);
  for(int i=1;i<mu_ExoverEtau_hplus.at(tdata).GetNbinsX();i++){
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

