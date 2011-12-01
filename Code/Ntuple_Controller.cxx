//Ntuple_Controller.cxx IMPLEMENTATION FILE

#ifndef Ntuple_Controller_cxx
#define Ntuple_Controller_cxx

#include "Ntuple_Controller.h"


///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,ObjEvent(-1)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("t");
  std::cout << "Ntuple_Controller" << RootFiles.size() << std::endl;
  for(int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;  
  if(chain==0){
    std::cout << "failed" << std::endl;
  }
  std::cout << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  std::cout << "Ntuple Configured" << std::endl;
}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  std::cout << Ntp->fChain->GetEntries() << std::endl;
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Int_t ientry = Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  std::cout << "Ntuple_Controller::~Ntuple_Controller()" << std::endl;
  delete Ntp;
  std::cout << "Ntuple_Controller::~Ntuple_Controller() complete" << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
    std::cout << "Starting D3PD cloning " << std::endl;  
    newfile = new TFile(n+".root","recreate");
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  std::cout << "Ntuple_Controller::SaveCloneTree()"<< std::endl;
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  std::cout << "Ntuple_Controller::SaveCloneTree() Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  std::cout << "Ntuple_Controller::ThinTree" << std::endl;

  std::cout << "Ntuple_Controller::ThinTree complete" << std::endl;  
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}


//Physics get Functions
int Ntuple_Controller::GetMCID(){
  if(isData())return Data;
  return MC;
}

TMatrixF     Ntuple_Controller::Vtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      if(j*dim+k<Ntp->Vtx_Cov->at(i).size()){
	M[j][k]=Ntp->Vtx_Cov->at(i).at(j*dim+k);
	M[k][j]=Ntp->Vtx_Cov->at(i).at(j*dim+k);
      }
    }
  }
}

bool Ntuple_Controller::isVtxGood(unsigned int i){
  if(0<=i && i<NVtx()){
    if(Vtx_Track_idx(i).size()>4)return true;
  }
  return false;
}


#endif



