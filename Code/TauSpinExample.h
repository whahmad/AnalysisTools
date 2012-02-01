#ifndef TauSpinExample_h
#define TauSpinExample_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TauSpinExample : public Selection {

 public:
  TauSpinExample(TString Name_, TString id_);
  virtual ~TauSpinExample();

  virtual void  Configure();

  enum cuts {isZtautauto3pimu=0,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> LongitudinalPolarization;
  std::vector<TH1D> LongitudinalPolarization_Spin;
  std::vector<TH1D> LongitudinalPolarization_UnSpin;
  std::vector<TH1D> LongitudinalPolarization_FlipSpin;

  std::vector<TH1D> mu_PmuoverEtau;
  std::vector<TH1D> mu_PmuoverEtau_hplus;
  std::vector<TH1D> mu_PmuoverEtau_hminus;
  std::vector<TH1D> mu_PmuoverEtau_Spin;
  std::vector<TH1D> mu_PmuoverEtau_UnSpin;
  std::vector<TH1D> mu_PmuoverEtau_FlipSpin;
  std::vector<TH1D> mu_WT_Spin;
  std::vector<TH1D> mu_WT_UnSpin;
  std::vector<TH1D> mu_WT_FlipSpin;

  std::vector<TH1D> pi_PmuoverEtau;
  std::vector<TH1D> pi_PmuoverEtau_hplus;
  std::vector<TH1D> pi_PmuoverEtau_hminus;
  std::vector<TH1D> pi_PmuoverEtau_Spin;
  std::vector<TH1D> pi_PmuoverEtau_UnSpin;
  std::vector<TH1D> pi_PmuoverEtau_FlipSpin;
  std::vector<TH1D> pi_WT_Spin;
  std::vector<TH1D> pi_WT_UnSpin;
  std::vector<TH1D> pi_WT_FlipSpin;

  std::vector<TH1D> a1_PmuoverEtau;
  std::vector<TH1D> a1_PmuoverEtau_hplus;
  std::vector<TH1D> a1_PmuoverEtau_hminus;
  std::vector<TH1D> a1_PmuoverEtau_Spin;
  std::vector<TH1D> a1_PmuoverEtau_UnSpin;
  std::vector<TH1D> a1_PmuoverEtau_FlipSpin;
  std::vector<TH1D> a1_WT_Spin;
  std::vector<TH1D> a1_WT_UnSpin;
  std::vector<TH1D> a1_WT_FlipSpin;

  std::vector<TH1D> rho_PmuoverEtau;
  std::vector<TH1D> rho_PmuoverEtau_hplus;
  std::vector<TH1D> rho_PmuoverEtau_hminus;
  std::vector<TH1D> rho_PmuoverEtau_Spin;
  std::vector<TH1D> rho_PmuoverEtau_UnSpin;
  std::vector<TH1D> rho_PmuoverEtau_FlipSpin;
  std::vector<TH1D> rho_WT_Spin;
  std::vector<TH1D> rho_WT_UnSpin;
  std::vector<TH1D> rho_WT_FlipSpin;

};
#endif
