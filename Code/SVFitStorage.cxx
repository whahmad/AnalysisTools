#include "SVFitStorage.h"
#include "Parameters.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTreeIndex.h"
#include "TFile.h"
#include "SVFitObject.h"
#include <iostream>

SVFitStorage::SVFitStorage():intree_loaded(false){
  gSystem->Load("libMy.so");

  outtree_= new TTree("SVFitTree", "SVFitTree");
  outtree_->Branch("RunNumber", &RunNumber_);
  outtree_->Branch("EventNumber", &EventNumber_);
  outtree_->Branch("svfit", &svfit_);

  LoadTree();
}

SVFitStorage::~SVFitStorage(){
  SaveTree();
  if(intree_loaded) infile_->Close();
}

void SVFitStorage::LoadTree(){
  TString key="InputAuxiliaryFile:";
  TString InFile="SVFIT";
  int nfiles=GetFile(InFile,key);
  if(nfiles>0){
    std::cout << "Error: unble to load: " << key << std::endl;
  }
  else{
    TDirectory *gdirectory_save= gDirectory;
    TChain *chain = new TChain("SVFitTree");
    for(int i=0;i<nfiles;i++){
      TString name=InFile;name+=i;name+="0.root";
      chain->Add(name);
    }
    intree_ = (TTree*)chain;
    
    //Set branches
    intree_->SetBranchAddress("RunNumber",&RunNumber_, &b_RunNumber_);
    intree_->SetBranchAddress("EventNumber",&EventNumber_,&b_EventNumber_);
    intree_->SetBranchAddress("svfit",&svfit_,&b_svfit_);
    
    //Build index tree
    intree_->BuildIndex("RunNumber","EventNumber");
    index_ = (TTreeIndex*)intree_->GetTreeIndex();
    gDirectory=gdirectory_save;
    gDirectory->cd();
    intree_loaded=true;
  }
}


void SVFitStorage::SaveTree(){
  //Load File name
  Parameters Par; // assumes configured in Analysis.cxx
  TString File;
  Par.GetString("OutputAuxiliaryFile:",File);
  
  //Save output
  TDirectory *gdirectory_save= gDirectory;
  TFile *outfile_=TFile::Open(File,"RECREATE");
  if (!outfile_) {
    std::cout << "ERROR:  " << File << " not saved" << std::endl;
    return;
  }
  outfile_->cd();
  outtree_->Write("SVFit_Tree");
  outfile_->Close();
  gDirectory = gdirectory_save;
  gDirectory->cd();
  
  //Store file on the grid
  StoreFile(File);
}

void SVFitStorage::SaveEvent(Int_t RunNumber,Int_t EventNumber,SVFitObject svfit){
  //Fill event
  RunNumber_=RunNumber;
  EventNumber_=EventNumber;
  svfit_=svfit;
  outtree_->Fill();
}

const SVFitObject& SVFitStorage::GetEvent(Int_t RunNumber,Int_t EventNumber){
  //Get tree entry using index and then get svfit
  Long64_t local = index_->GetEntryNumberWithIndex(RunNumber,EventNumber);
  b_svfit_->GetEntry(local);
  return svfit_;
}


