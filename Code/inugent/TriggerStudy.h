#ifndef TriggerStudy_h
#define TriggerStudy_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class TriggerStudy : public Selection {

 public:
  TriggerStudy(TString Name_, TString id_);
  virtual ~TriggerStudy();

  virtual void  Configure();

  enum cuts {TriggerOk=0,NTau,TriggerMatch,NCuts};

  virtual void Finish();

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<std::vector<TH1D> > N_pt;
  std::vector<std::vector<TH1D> > N_eta;
  std::vector<std::vector<TH1D> > HpsMode;\
  std::vector<std::vector<TH2D> > N_2d;
  std::vector<TH1D> Mtautau;

};
#endif
