#ifndef TriggerStudyMC_h
#define TriggerStudyMC_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TriggerStudyMC : public Selection {

 public:
  TriggerStudyMC(TString Name_, TString id_);
  virtual ~TriggerStudyMC();

  virtual void  Configure();

  enum cuts {isSignal=0,NCuts};

  virtual void Finish();

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<std::vector<TH1D> > N_pt;
  std::vector<std::vector<TH1D> > n_pt;
  std::vector<std::vector<TH1D> > Eff_pt;
  std::vector<std::vector<TH1D> > N_eta;
  std::vector<std::vector<TH1D> > n_eta;
  std::vector<std::vector<TH1D> > Eff_eta;
  std::vector<std::vector<TH1D> > N_abs;
  std::vector<std::vector<TH1D> > n_abs;
  std::vector<std::vector<TH1D> > Eff_abs;

  std::vector<std::vector<TH1D> > N_pt_vis;
  std::vector<std::vector<TH1D> > n_pt_vis;
  std::vector<std::vector<TH1D> > Eff_pt_vis;
  std::vector<std::vector<TH1D> > N_eta_vis;
  std::vector<std::vector<TH1D> > n_eta_vis;
  std::vector<std::vector<TH1D> > Eff_eta_vis;

  std::vector<std::vector<TH1D> > N_pt_visreco;
  std::vector<std::vector<TH1D> > n_pt_visreco;
  std::vector<std::vector<TH1D> > Eff_pt_visreco;
  std::vector<std::vector<TH1D> > N_eta_visreco;
  std::vector<std::vector<TH1D> > n_eta_visreco;
  std::vector<std::vector<TH1D> > Eff_eta_visreco;


  std::vector<std::vector<TH2D> > N_2d_vis;
  std::vector<std::vector<TH2D> > n_2d_vis;
  std::vector<std::vector<TH2D> > Eff_2d_vis;

  std::vector<std::vector<TH2D> > n_2d_vis_Triggershift5GeV;
  std::vector<std::vector<TH2D> > n_2d_vis_Triggershift10GeV;

  bool doKFTau, doHPSLooseTau, doHPSMediumTau, doHPSTightTau;

};
#endif
