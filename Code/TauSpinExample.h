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
  std::vector<TH1D> PmuoverEtau;
  std::vector<TH1D> PmuoverEtau_hplus;
  std::vector<TH1D> PmuoverEtau_hminus;
  std::vector<TH1D> WT_Spin;
  std::vector<TH1D> WT_UnSpin;
  std::vector<TH1D> WT_FlipSpin;
  std::vector<TH1D> LongitudinalPolarization;
};
#endif
