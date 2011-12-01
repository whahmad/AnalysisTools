#define GoodRunList_cxx

#include "GoodRunList.h"

#include <cstdlib>
#include <algorithm>                                                                                                                                
#include <iostream>
#include <fstream>
#include <cstdlib>


GoodRunList::GoodRunList():
  state(false),
  run(-1),
  lumiblock(-1),
  Configured(false)
{
}

GoodRunList::GoodRunList(TString file):
  state(false),
  run(-1),
  lumiblock(-1),
  Configured(false)
{
  ConfigureGRL(file);
}

GoodRunList::~GoodRunList(){

}


void GoodRunList::ConfigureGRL(TString file){
  LumiBlockRange.clear();
  RunMap.clear();
  GRLfilename=file;
  ifstream input_file;
  char *file_=(char*)file.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening xml file "<< file <<" for GoodRunList has failed.\n" << std::endl;
    exit(1);
  }
  ////////////////////////
  // Read in the xml file
  std::cout << "\nOpened GoodRunList xml file: "<< file <<".\n" << std::endl;
  bool flag(true);
  while(flag){
    TString line;
    flag=Get_Line(input_file,line);
    if(flag && line!=""){
      if(line.Contains("<LumiBlockCollection>")){
	flag=Get_Line(input_file,line);
	if(line.Contains("<Run>")){
	  line.ReplaceAll(" ","");
	  line.ReplaceAll("<Run>","");
	  line.ReplaceAll("<Run>","</Run>");
	  int therun=line.Atoi();
	  RunMap.insert(std::pair<int,int>(therun,LumiBlockRange.size()));
	  LumiBlockRange.push_back(std::vector<std::pair<int,int> >());
	  bool lumi=true;
	  //std::cout << "Found Run: " << therun << std::endl;
	  while(lumi && flag){
	    flag=Get_Line(input_file,line);
	    int start,end;
	    if(GetLumiRange(line,start,end)){
	      LumiBlockRange[LumiBlockRange.size()-1].push_back(std::pair<int,int>(start,end));
	      //std::cout << "Found LumiRange: " << start << " to " << end << std::endl;
	    }
	    else{
	      lumi=false;
	      break;
	    }
	  }
	}
      }
    }
  }
  input_file.close();
  Configured=true;
}

bool GoodRunList::isGoodLumiBlock(int Run, int LumiBlock ) {
  if(run==Run && lumiblock==LumiBlock && Configured){
    return state;
  }
  state=false;
  lumiblock=LumiBlock;
  run=Run;
  if(RunMap.count(run)){
    int runidx=RunMap[run];
    if(runidx<LumiBlockRange.size()){
      for(int i=0; i<LumiBlockRange[runidx].size();i++){
	if(LumiBlockRange[runidx][i].first<=lumiblock && lumiblock<=LumiBlockRange[runidx][i].second){
	  state=true;
	  return state;
	}
      }
    }
  }
  return state;
}




bool GoodRunList::Get_Line(ifstream &input_file, TString &Item){
  Item="";
  char buffer[5000];
  bool flag=false;
  if(!input_file.eof()){
    input_file.getline(buffer,5000);
    Item=buffer;
    flag=true;
    //std::cout << Item << std::endl;
  }
  if(Item==""){
    flag=false;
  }
  return flag;
}

bool GoodRunList::GetLumiRange(TString &l,int &start, int &end){
  if(!l.Contains("LBRange")) return false;
  int size=l.Length();
  std::vector<TString> v;
  TString s="";
  start=-1;
  end=-1;
  for(int i=0; i<size;i++){
    TString c=l.Data()[i];
    if(s.Length()>0 && c==" "){
      v.push_back(s);
      s="";
    }
    else if(c!=" "){
      s+=c;
    }
  } 
  v.push_back(s);
  for(int i=0; i<v.size();i++){
    if(v[i].Contains("Start")){
      v[i].ReplaceAll("Start=","");
      v[i].ReplaceAll("\"","");
      start=v[i].Atoi();
    }
    if(v[i].Contains("End")){
      v[i].ReplaceAll("End=","");
      v[i].ReplaceAll("\"","");
      end=v[i].Atoi();
    }
  }
  return true;
}
