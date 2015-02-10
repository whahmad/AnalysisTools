#ifndef DataStorage_h
#define DataStorage_h

#include <vector>
#include "TString.h"

class DataStorage {
 public:
  DataStorage();
  virtual ~DataStorage();

  int  GetFile(TString InFile, TString key);
  void StoreFile(TString File, TString savedFile);

 private:
  TString mydir;
};
#endif
