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
std::vector<float>   SkimConfig::nevents;
std::vector<TString> SkimConfig::dsname;

std::vector<TString> SkimConfig::SkimIDs;
std::vector<float>   SkimConfig::EventsPassed;
std::vector<float>   SkimConfig::EventsPassedErr;
std::vector<float>   SkimConfig::EventEff_noweight;


SkimConfig::SkimConfig():
  loaded(false)
{
}

bool SkimConfig::Load()
{
  return SkimConfig::Load("Tools/DataSets.dat");
}

bool SkimConfig::Load(TString Name_)  
{
  if(loaded) return true;
  std::cout << "SkimConfig::Load("<< Name_ <<")" << std::endl;
  nevents.clear();
  dsname.clear();
  if(Name_==""){
    loaded=true;
    return true;
  }
  ifstream input_file;
  char *file_=(char*)Name_.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening xml file "<< Name_ <<" for SkimConfig has failed.\n" << std::endl;
    return false;
  }
  loaded=true;
  std::cout << "\nOpened SkimConfig xml file: "<< Name_ <<".\n" << std::endl;

  std::string s;
  unsigned int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>1000) break;
    std::stringstream line(s); 
    TString ds;
    int nfiles;
    int nevt;
    line >> ds >> nfiles >> nevt;
    dsname.push_back(ds);
    nevents.push_back(nevt);
  }
  input_file.close();
  for(int i=0; i<dsname.size();i++){
    std::cout << "SkimConfig::Load " << i << " "  << dsname.at(i) << " NEvents=" << nevents.at(i) << std::endl;  
  }
  std::cout << "SkimConfig::Load("<< Name_ <<") complete" << std::endl;
}

void SkimConfig::LoadSkimEff(TString Name_)
{
  std::cout << "SkimConfig::LoadSkimEff("<< Name_ <<")" << std::endl;
  bool first=true;
  //if(SkimIDs.size()>0)first=false;

  SkimIDs.clear();
  EventsPassed.clear();
  EventsPassedErr.clear();
  EventEff_noweight.clear();
  ifstream input_file;
  char *file_=(char*)Name_.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening SkimEff file "<< Name_ <<" for SkimConfig has failed.\n" << std::endl;
    return;
  }
  std::cout << "\nOpened SkimConfig SkimEff file: "<< Name_ <<".\n" << std::endl;

  std::string s;
  unsigned int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>1000) break;
    std::stringstream line(s);
    TString id;
    float nevents;
    float neventserr;
    float c;
    float noweighteff;
    line >> id >> nevents >> neventserr >> c >> noweighteff;
    if(first){
      std::cout << "First" << std::endl;
      SkimIDs.push_back(id);
      EventsPassed.push_back(nevents);
      EventsPassedErr.push_back(neventserr);
      EventEff_noweight.push_back(noweighteff);
    }
    else{
      for(int j=0; j<SkimIDs.size();j++){
	std::cout << j << std::endl;
	if(SkimIDs.at(j)==id){
	  EventEff_noweight.at(j)+=noweighteff;
	}
	std::cout << j << " " << "done" << std::endl;
      }
    }

  }
  input_file.close();
  std::cout << "SkimConfig::LoadSkimEff " << SkimIDs.size() << " " << EventsPassed.size() << " " << EventsPassedErr.size() 
	    << " " << EventEff_noweight.size() << std::endl;
  for(int i=0; i<SkimIDs.size();i++){
    std::cout << "SkimConfig::LoadSkimEff " << i << " "  << SkimIDs.at(i) << " " << EventsPassed.at(i) << " " 
	      << EventsPassedErr.at(i) << " " << EventEff_noweight.at(i) << std::endl;
  }
  std::cout << "SkimConfig::LoadSkimEff("<< Name_ <<") complete" << std::endl;
}

SkimConfig::~SkimConfig(){

}


void SkimConfig::SaveEfficiency(TString Name,std::vector<TString> ids,std::vector<TH1D> NPassed, std::vector<TH1D> NPassed_noweight){
  ofstream output;
  output.open(Name+"SkimEff.dat", std::ios::out);
  int nbins=NPassed.at(0).GetNbinsX();
  for(int i=0; i<ids.size();i++){
    ids.at(i).ReplaceAll(" ","");
    //if(NPassed_noweight.at(i).GetBinContent(1)!=0 && NPassed.at(i).GetBinContent(1)!=0){
    (output) << ids.at(i)  << setprecision (15)
	     << " " << NPassed.at(i).GetBinContent(1) << " " << NPassed.at(i).GetBinError(1) 
	     << " " << NPassed_noweight.at(i).GetBinContent(1) 
	     << " " << NPassed_noweight.at(i).GetBinContent(nbins)/NPassed_noweight.at(i).GetBinContent(1)
	     << " " << NPassed.at(i).GetBinContent(nbins) 
	     << " " << NPassed_noweight.at(i).GetBinContent(nbins)
	     << std::endl;
      cout << ids.at(i)
	   << " " << NPassed.at(i).GetBinContent(1) << " " << NPassed.at(i).GetBinError(1)
	   << " " << NPassed_noweight.at(i).GetBinContent(1)
	   << " " << NPassed_noweight.at(i).GetBinContent(nbins)/NPassed_noweight.at(i).GetBinContent(1)
	   << std::endl;

  }
}

void SkimConfig::ApplySkimEfficiency(std::vector<TString> ids,std::vector<TH1D> &NPassed, std::vector<TH1D> &NPassed_noweight){
  for(int i=0; i<ids.size();i++){
    std::cout << "ApplySkimEfficiency" << ids.at(i) << " " << SkimIDs.size() << std::endl;
    for(int j=0; j<SkimIDs.size();j++){
      ids.at(i).ReplaceAll(" ","");
      if(SkimIDs.at(j).Contains(ids.at(i))){
	NPassed.at(i).SetBinContent(1,EventsPassed.at(j));
	NPassed.at(i).SetBinError(1,EventsPassedErr.at(j));
        NPassed.at(i).SetBinContent(0,EventsPassed.at(j));
        NPassed.at(i).SetBinError(0,EventsPassedErr.at(j));
	float value=NPassed_noweight.at(i).GetBinContent(1);
	if(EventEff_noweight.at(j)!=0){
	  NPassed_noweight.at(i).SetBinContent(1,value/EventEff_noweight.at(j));
	  NPassed_noweight.at(i).SetBinError(1,sqrt(value/EventEff_noweight.at(j)));
	}
      }
    }
  }
}



void SkimConfig::CorrectNEvents(std::vector<TString> ids, std::vector<float> nevts){
  if(!loaded){
    std::cout << "Error SkimConfig::CorrectEvents input not loaded -> no skim summary " << std::endl;  
  }

  if(ids.size()!=nevts.size()){
    std::cout << "Error SkimConfig::CorrectEvents invalid sizes "  << ids.size() << " " << nevts.size() << std::endl;
    return;
  }
  for(int i=0; i<ids.size();i++){
    float sumevts=0; 
    for(int j=0; j<dsname.size();j++){
      if(dsname[j].Contains(ids[i])){
	sumevts+=nevents[j];
      }
    }	
    if(nevts[i]==0){
      std::cout << "Warning " << ids[i] << " not used " << std::endl;
    }
    else if(fabs(nevts[i]-sumevts)>0.01){
      std::cout << "Failed " << ids[i] << " incorrect number of events "
		<< "Found: " << nevts[i] << " Expected: " << sumevts 
		<< " Ratio " <<   nevts[i]/sumevts
		<< std::endl;
    }
    else if(fabs(nevts[i]-sumevts)<0.01){
      std::cout << "Passed " << ids[i] << " all events analysed "
		<< "Found: " << nevts[i] << " Expected: " << sumevts
		<< std::endl;
    }
  }
  return;
}




