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
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Accumulative Cuts Passed","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  PmuoverEtau=HConfig.GetTH1D(Name+"_PmuoverEtau","PmuoverEtau",100,0.0,4.0,"P_{#mu}/P_{#tau}","Events");
  PmuoverEtau_hplus=HConfig.GetTH1D(Name+"_PmuoverEtau_hplus","PmuoverEtau_hplus",100,0.0,100.0,"P_{#mu}/P_{#tau}","Events");
  PmuoverEtau_hminus=HConfig.GetTH1D(Name+"_PmuoverEtau_hminus","PmuoverEtau_hminus",100,0.0,100,"P_{#mu}/P_{#tau}","Events");

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
}

void  TauSpinExample::doEvent(){
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  unsigned int idx;
  value.at(isZtautauto3pimu)=0;
  pass.at(isZtautauto3pimu) =  Ntp->hasSignalTauDecay(PdtPdgMini::Z0,TauDecay::JAK_MUON,idx);
  if(pass.at(isZtautauto3pimu))value.at(isZtautauto3pimu)=1;
  double wobs=1;
  double w=1;
  /*  if(!Ntp->isData()){
    w*=Ntp->EvtWeight3D();
    }*/

  std::cout << Ntp->GetMCID() << " " << Npassed.size() << " " << t << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  std::cout << "Status" << std::endl; 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
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
    double WT=Ntp->TauSpinerWeight(TauSpinerInterface::FlipSpin);
    std::cout <<  "TauSpiner WT: " << WT << std::endl;

    PmuoverEtau.at(t).Fill(Ntp->TauSpinerWeight(TauSpinerInterface::Spin));
    PmuoverEtau_hplus.at(t).Fill(Ntp->TauSpinerWeight(TauSpinerInterface::UnSpin));
    PmuoverEtau_hminus.at(t).Fill(Ntp->TauSpinerWeight(TauSpinerInterface::FlipSpin));


  }
  std::cout << "done" << std::endl;
}




