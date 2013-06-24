#include "Selection_Factory.h"

#include "Example.h"
#include "cherepanov/Validation.h"
#include "cherepanov/Ztotautau_hadmu_ControlSample.h"
#include "cherepanov/Tau_momentum_calculation.h"
#include "inugent/Ztotautau_ControlSample.h"
#include "TauSpinExample.h"
#include "inugent/ChargedHiggs_dilepontic.h"
#include "inugent/ChargedHiggs_tauplusjet.h"
#include "inugent/ZDouble3prong.h"
#include "inugent/Ztomumu_ControlSample.h"
#include "inugent/TriggerStudy.h"
#include "inugent/TriggerStudyMC.h"
#include "nehrkorn/ZtoEMu.h"
#include "nehrkorn/ZtoEMu_Skim.h"

Selection_Factory::Selection_Factory(){
}

Selection_Factory::~Selection_Factory(){
}

Selection_Base* Selection_Factory::Factory(TString Analysis, TString UncertType,int mode, int runtype, double lumi){
  Selection_Base* s;
  Analysis.ToLower();

  /*if(Analysis.Contains("validation"))s=new Validation(Analysis,UncertType);
  else if(Analysis.Contains("tau_momentum_calculation"))s=new Tau_momentum_calculation(Analysis,UncertType);
  else if(Analysis.Contains("tauspinexample"))s=new TauSpinExample(Analysis,UncertType);
  else if(Analysis.Contains("triggerstudymc"))s=new TriggerStudyMC(Analysis,UncertType);
  else if(Analysis.Contains("triggerstudy"))s=new TriggerStudy(Analysis,UncertType);
  else */if(Analysis.Contains("ztoemu_mcsample"))s=new ZtoEMu(Analysis,UncertType);
  else if(Analysis.Contains("ztoemu_skim"))s=new ZtoEMu_Skim(Analysis,UncertType);
  else{
    std::cout << "WARNING: Selection_Factory::Factory INVALID ANALYSIS TYPE.... USING DEFAULT <Example.h> " << std::endl;
    s=new Example(Analysis,UncertType);
  }
  s->Configure();
  s->SetMode(mode);
  s->SetRunType(runtype);
  s->SetLumi(lumi);
  return s;
}
