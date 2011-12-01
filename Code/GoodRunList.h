#ifndef GoodRunList_h
#define GoodRunList_h

#include <vector>
#include <map>
#include "TString.h"

class GoodRunList{

 public:
  GoodRunList();
  GoodRunList(TString file);
  bool isGoodLumiBlock(int Run, int LumiBlock);
  void ConfigureGRL(TString file);
  inline TString GetGRLFileName(){return GRLfilename;}

  ~GoodRunList();

 private:
  bool Get_Line(ifstream &input_file, TString &Item);
  bool GetLumiRange(TString &l,int &start, int &end);

  bool state,Configured;
  int run;
  int lumiblock;

  std::vector<std::vector<std::pair<int,int> > > LumiBlockRange; // [Runidx][LumiBlocks] 
  std::map<int,int> RunMap; //First=Run Second=Runidx
  TString GRLfilename;
};
#endif
