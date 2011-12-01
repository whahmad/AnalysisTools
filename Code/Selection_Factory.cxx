#include "Selection_Factory.h"

#include "Example.h"

Selection_Factory::Selection_Factory(){
}

Selection_Factory::~Selection_Factory(){
}

Selection_Base* Selection_Factory::Factory(TString Analysis, TString UncertType,int mode, int runtype){
  Selection_Base* s;
  Analysis.ToLower();
  if(Analysis.Contains("example"))s=new Example(Analysis,UncertType);
  else{
    std::cout << "WARNING: Selection_Factory::Factory INVALID ANALYSIS TYPE.... USING DEFAULT <Example.h> " << std::endl;
    s=new Example(Analysis,UncertType);
  }
  s->Configure();
  s->SetMode(mode);
  s->SetRunType(runtype);
  return s;
}
