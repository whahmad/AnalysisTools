#include "DataStorage.h"
#include "Parameters.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <fstream>
#include "TString.h"
#include <iostream>

DataStorage::DataStorage(){
  mydir= getenv ("PWD");
}

DataStorage::~DataStorage(){

}

int DataStorage::GetFile(TString InFile, TString key){
  Parameters Par; // assumes configured in Analysis.cxx
  TString gridsite="none";
  std::vector<TString> Files;
  Par.GetVectorString(key,Files);
  Par.GetString("GRIDSite:",gridsite);

  if(gridsite=="none" || Files.size()==0) return 0;
  
  for(unsigned int i=0;i<Files.size();i++){
    TString inFile=InFile; inFile+=i; inFile+=".root";
    TString cmd1= "srmcp srm://" + gridsite + ":8443/" + Files.at(i) + " file:////" + mydir + "/" + inFile;
    std::cout << "cmd:" << cmd1 << std::endl;
    system(cmd1.Data());
    ifstream f(inFile);
    if(!f) return 0;
  }
  return Files.size();
}

void DataStorage::StoreFile(TString File, TString savedFile){
  TString gridsite;
  Parameters Par; // assumes configured in Analysis.cxx
  Par.GetString("GRIDSite:",gridsite);
  TString cmd1 = "srmcp file:////" + mydir + "/" + File + " srm://" + gridsite + ":8443/" + savedFile;
  std::cout << "cmd:" << cmd1 << std::endl;
  system(cmd1.Data());
}
