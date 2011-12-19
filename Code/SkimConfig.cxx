#include "SkimConfig.h"

#include <cstdlib>
#include <algorithm>                                     
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "Riostream.h"
#include <fstream>
#include <istream>
#include <strstream>
#include <cstdlib>
#include "TH1D.h"
#include <math.h>

std::vector<int> SkimConfig::SkimIDs;
std::vector<float>   SkimConfig::NEvents;
std::vector<float>   SkimConfig::NEventsErr;
std::vector<float>   SkimConfig::NEvents_sel;
std::vector<float>   SkimConfig::NEventsErr_sel;
std::vector<float>   SkimConfig::NEvents_noweight;
std::vector<float>   SkimConfig::NEvents_noweight_sel;
bool SkimConfig::loaded=false;

SkimConfig::SkimConfig()
{
}

bool SkimConfig::Load()
{
  return SkimConfig::Load("Tools/DataSets.dat");
}

bool SkimConfig::Load(TString Name_)
{
  std::cout << "SkimConfig::LoadSkimEff("<< Name_ <<")" << std::endl;
  if(loaded) return false;
  SkimIDs.clear();
  NEvents.clear();
  NEventsErr.clear();
  NEvents_sel.clear();
  NEventsErr_sel.clear();
  NEvents_noweight.clear();
  NEvents_noweight_sel.clear();
  ifstream input_file;
  char *file_=(char*)Name_.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening SkimEff file "<< Name_ <<" for SkimConfig has failed.\n" << std::endl;
    return false;
  }
  std::cout << "\nOpened SkimConfig SkimEff file: "<< Name_ <<".\n" << std::endl;
  loaded =true;
  std::string s;
  int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>1000) break;
    std::stringstream line(s);
    int id;
    float nevents;
    float neventserr;
    float nevents_sel;
    float neventserr_sel;
    float noweight;
    float noweight_sel;
    line >> id >> nevents >> neventserr >> nevents_sel >> neventserr_sel >> noweight >> noweight_sel;
    SkimIDs.push_back(id);
    NEvents.push_back(nevents);
    NEventsErr.push_back(neventserr);
    NEvents_sel.push_back(nevents_sel);
    NEventsErr_sel.push_back(neventserr_sel);
    NEvents_noweight.push_back(noweight);
    NEvents_noweight_sel.push_back(noweight_sel);
  }
  input_file.close();
  return true;
}

SkimConfig::~SkimConfig(){

}


double SkimConfig::GetNEvents(int id){
  for(int i=0; i<SkimIDs.size();i++){
    if(id==SkimIDs.at(i))return NEvents.at(i);
  }
  return 0;
}

void SkimConfig::SaveEfficiency(TString Name,std::vector<int> ids,std::vector<TH1D> NPassed, std::vector<TH1D> NPassed_noweight){
  ofstream output;
  output.open(Name+"SkimEff.dat", std::ios::out);
  int nbins=NPassed.at(0).GetNbinsX();
  for(int i=0; i<ids.size();i++){
    (output) << ids.at(i)  << setprecision (15)
	     << " " << NPassed.at(i).GetBinContent(1)           << " " << NPassed.at(i).GetBinError(1)
	     << " " << NPassed.at(i).GetBinContent(nbins)       << " " << NPassed.at(i).GetBinError(nbins) 
	     << " " << NPassed_noweight.at(i).GetBinContent(1)  << " " << NPassed_noweight.at(i).GetBinContent(1) 
	     << std::endl;
  }
}

void SkimConfig::ApplySkimEfficiency(std::vector<int> ids,std::vector<TH1D> &NPassed, std::vector<TH1D> &NPassed_noweight){
  std::vector<bool> hasid(SkimIDs.size(),false);
  for(unsigned int i=0; i<ids.size();i++){
    for(unsigned int j=0; j<SkimIDs.size();j++){
      if(SkimIDs.at(i)==ids.at(i))hasid.at(j)=true;
    }
  }
  for(unsigned int i=0; i<ids.size();i++){
    NPassed.at(i).SetBinContent(1,0);
    NPassed.at(i).SetBinError(1,0);
    NPassed.at(i).SetBinContent(0,0);
    NPassed.at(i).SetBinError(0,0);
    NPassed_noweight.at(i).SetBinContent(1,0);
    NPassed_noweight.at(i).SetBinError(1,0);
  }
  for(unsigned int j=0; j<SkimIDs.size();j++){
    for(unsigned int i=0; i<ids.size();i++){
      if((hasid.at(j) && ids.at(i)==SkimIDs.at(i)) || (!hasid.at(j) && ids.at(i)==SkimIDs.at(i)%100)){
	std::cout << "ApplySkimEfficiency" << ids.at(i) << " " << SkimIDs.size() << std::endl;
	float tmp=NPassed.at(i).GetBinContent(1);
	NPassed.at(i).SetBinContent(1,tmp+NEvents.at(j));
	tmp=NPassed.at(i).GetBinError(1);
	NPassed.at(i).SetBinError(1,sqrt(tmp*tmp+NEventsErr.at(j)*NEventsErr.at(j)));
	NPassed.at(i).SetBinContent(0,NPassed.at(i).GetBinContent(1));
	NPassed.at(i).SetBinError(0,NPassed.at(i).GetBinError(1));
	tmp=NPassed_noweight.at(i).GetBinContent(1);
	NPassed_noweight.at(i).SetBinContent(1,tmp+NEvents_noweight.at(j));
	tmp=NPassed_noweight.at(i).GetBinError(1);
	NPassed_noweight.at(i).SetBinError(1,sqrt(tmp*tmp+NEvents_noweight.at(j)*NEvents_noweight.at(j)));
      }
    }
  }
}




void SkimConfig::CheckNEvents(std::vector<int> ids, std::vector<float> nevts){
  if(!loaded){
    std::cout << "Error SkimConfig::CorrectEvents input not loaded -> no skim summary " << std::endl;  
    return;
  }
  std::vector<bool> hasid(SkimIDs.size(),false); 
  std::vector<float> skimnevts(nevts.size(),0);
  for(unsigned int i=0; i<ids.size();i++){
    for(unsigned int j=0; j<SkimIDs.size();j++){
      if(SkimIDs.at(i)==ids.at(i))hasid.at(j)=true;
    }
  }
  for(unsigned int j=0; j<SkimIDs.size();j++){
    for(unsigned int i=0; i<ids.size();i++){
      if((hasid.at(j) && ids.at(i)==SkimIDs.at(j)) || (!hasid.at(j) && ids.at(i)==SkimIDs.at(j)%100)) skimnevts.at(i)+=nevts.at(i);
    }
  }
  for(unsigned int i=0; i<ids.size();i++){
    if(fabs(nevts.at(i)-skimnevts.at(i))>0.01){
      std::cout << "Failed " << ids.at(i) << " incorrect number of events "
		<< "Found: " << nevts.at(i) << " Expected: " << skimnevts.at(i) 
		<< " Ratio " <<   nevts.at(i)/skimnevts.at(i)
		<< std::endl;
    }
    else if(fabs(nevts.at(i)-skimnevts.at(i))<0.01){
      std::cout << "Passed " <<  ids.at(i) << " all events analysed "
		<< "Found: " << nevts.at(i) << " Expected: " << skimnevts.at(i)
		<< std::endl;
    }
  }
  return;
}




