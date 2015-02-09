#ifndef SVFitStorage_h
#define SVFitStorage_h

#include <vector>
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTreeIndex.h"
#include "TFile.h"
#include "SVFitObject.h"
#include "DataStorage.h"
#include "TBranch.h"

class SVFitStorage : public DataStorage {
 public:
  SVFitStorage();
  ~SVFitStorage();

  void LoadTree();
  void SaveTree();
  void SaveEvent(Int_t RunNumber,Int_t EventNumber,SVFitObject svfit);
  const SVFitObject& GetEvent(Int_t RunNumber,Int_t EventNumber);
  
 private:
  TFile *infile_, *outfile_;
  TTree *intree_,*outtree_;
  TTreeIndex *index_;
  
  UInt_t RunNumber_;
  UInt_t EventNumber_;
  SVFitObject svfit_; 

  TBranch *b_RunNumber_;
  TBranch *b_EventNumber_;
  TBranch *b_svfit_;

  bool intree_loaded;
};
#endif
