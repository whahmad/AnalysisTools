#include "TauSpinExample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

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
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
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
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Verrtex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  PmuoverEtau=HConfig.GetTH1D(Name+"_PmuoverEtau","PmuoverEtau",20,0.0,1.0,"P_{#mu}/P_{#tau}","Events");
  PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_PmuoverEtau_hplus","PmuoverEtau_hplus",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{h^{+}}","Events");
  PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_PmuoverEtau_hminus","PmuoverEtau_hminus",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{h^{-}}","Events");
  PmuoverEtau_Spin=HConfig.GetTH1D(Name+"_PmuoverEtau_Spin","PmuoverEtau_Spin",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{Spin}","Events");
  PmuoverEtau_UnSpin=HConfig.GetTH1D(Name+"_PmuoverEtau_UnSpin","PmuoverEtau_UnSpin",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{UnSpin}","Events");
  PmuoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_PmuoverEtau_FlipSpin","PmuoverEtau_FlipSpin",20,0.0,1.0,"P_{#mu}/P_{#tau}|_{FlipSpin}","Events");
  WT_Spin=HConfig.GetTH1D(Name+"_WT_Spin","WT_Spin",40,0.0,4.0,"WT","Events");
  WT_UnSpin=HConfig.GetTH1D(Name+"_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT","Events");
  WT_FlipSpin=HConfig.GetTH1D(Name+"_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)","Events");
  LongitudinalPolarization=HConfig.GetTH1D(Name+"_LongitudinalPolarization","LongitudinalPolarization",20,-2.0,2.0,"#rho_{l}","Events");
  LongitudinalPolarization_Spin=HConfig.GetTH1D(Name+"_LongitudinalPolarization_Spin","LongitudinalPolarization_Spin",20,-2.0,2.0,"#rho_{l}|_{Spin}","Events");
  LongitudinalPolarization_UnSpin=HConfig.GetTH1D(Name+"_LongitudinalPolarization_UnSpin","LongitudinalPolarization_UnSpin",20,-2.0,2.0,"#rho_{l}|_{UnSpin}","Events");
  LongitudinalPolarization_FlipSpin=HConfig.GetTH1D(Name+"_LongitudinalPolarization_FlipSpin","LongitudinalPolarization_FlipSpin",20,-2.0,2.0,"#rho_{l}|_{FlipSpin}","Events");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  TauSpinExample::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist1d.push_back(&PmuoverEtau);
 Extradist1d.push_back(&PmuoverEtau_hplus);
 Extradist1d.push_back(&PmuoverEtau_hminus);
 Extradist1d.push_back(&PmuoverEtau_Spin);
 Extradist1d.push_back(&PmuoverEtau_UnSpin);
 Extradist1d.push_back(&PmuoverEtau_FlipSpin);
 Extradist1d.push_back(&WT_Spin);
 Extradist1d.push_back(&WT_UnSpin);
 Extradist1d.push_back(&WT_FlipSpin);
 Extradist1d.push_back(&LongitudinalPolarization);
 Extradist1d.push_back(&LongitudinalPolarization_Spin);
 Extradist1d.push_back(&LongitudinalPolarization_UnSpin);
 Extradist1d.push_back(&LongitudinalPolarization_FlipSpin);
}

void  TauSpinExample::doEvent(){
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,tau_idx;
  value.at(isZtautauto3pimu)=0;
  pass.at(isZtautauto3pimu) =  Ntp->hasSignalTauDecay(PdtPdgMini::Z0,Boson_idx,TauDecay::JAK_MUON,tau_idx);
  if(pass.at(isZtautauto3pimu))value.at(isZtautauto3pimu)=1;
  double wobs=1;
  double w=1;
  /*  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
    }*/
  if(verbose)  std::cout << Ntp->GetMCID() << " " << Npassed.size() << " " << t << " " << Boson_idx << " " << " " << tau_idx 
	    << " " << Ntp->NMCTaus() << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  std::cout << "Status" << std::endl; 
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
    NGoodVtx.at(t).Fill(nGoodVtx,w);
    ////////////////////////////////////////////////
    //
    // Spin Validation
    //
    double WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
    if(verbose)std::cout <<  "TauSpiner WT: " << WT << std::endl;

    WT_Spin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::Spin),w);
    WT_UnSpin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::UnSpin),w);
    WT_FlipSpin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin),w);
    LongitudinalPolarization.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w);
    LongitudinalPolarization_Spin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w*Ntp->TauSpinerGet(TauSpinerInterface::Spin));
    LongitudinalPolarization_UnSpin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w*Ntp->TauSpinerGet(TauSpinerInterface::UnSpin));
    LongitudinalPolarization_FlipSpin.at(t).Fill(Ntp->TauSpinerGet(TauSpinerInterface::LPolarization),w*Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin));
    TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
    TLorentzVector Tau_LV(0,0,0,0);
    TLorentzVector Mu_LV(0,0,0,0);
    for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
      if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::tau_minus)) Tau_LV = Ntp->MCTauandProd_p4(tau_idx,i);
      if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PdtPdgMini::mu_minus)) Mu_LV  = Ntp->MCTauandProd_p4(tau_idx,i);
    }
    if(Tau_LV.E()>0){
      Tau_LV.Boost(-Boson_LV.BoostVector());
      Mu_LV.Boost(-Boson_LV.BoostVector());
      Boson_LV.Boost(-Boson_LV.BoostVector());
      PmuoverEtau.at(t).Fill(Mu_LV.E()/Tau_LV.E(),w);
      PmuoverEtau_Spin.at(t).Fill(Mu_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::Spin));
      PmuoverEtau_UnSpin.at(t).Fill(Mu_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::UnSpin));
      PmuoverEtau_FlipSpin.at(t).Fill(Mu_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin));
      PmuoverEtau_hplus.at(t).Fill(Mu_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hplus));
      PmuoverEtau_hminus.at(t).Fill(Mu_LV.E()/Tau_LV.E(),w*Ntp->TauSpinerGet(TauSpinerInterface::hminus));
    }
  }
  if(verbose)std::cout << "done" << std::endl;
}




