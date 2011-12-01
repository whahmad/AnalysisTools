#ifndef SkimConfig_h
#define SkimConfig_h

#include <vector>
#include "TString.h"
#include "TH1D.h"

class SkimConfig {

 public:
  SkimConfig();
  virtual ~SkimConfig();

  bool Load();
  bool Load(TString Name_);
  void LoadSkimEff(TString Name_);
  void CorrectNEvents(std::vector<TString> ids, std::vector<float> nevts); 
  void SaveEfficiency(TString Name,std::vector<TString> ids,std::vector<TH1D> NPassed, std::vector<TH1D> NPassed_noweight);
  void ApplySkimEfficiency(std::vector<TString> ids,std::vector<TH1D> &NPassed, std::vector<TH1D> &NPassed_noweight);


 private:
  static  std::vector<float>   nevents;
  static  std::vector<TString> dsname;
  static  std::vector<TString> SkimIDs;
  static  std::vector<float>   EventsPassed;
  static  std::vector<float>   EventsPassedErr;
  static  std::vector<float>   EventEff_noweight;

  bool loaded;
};
#endif
