#include "MuTauSync.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

MuTauSync::MuTauSync(TString Name_, TString id_):
  HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class MuTauSync" << std::endl;
	// always run without category for sync exercise
	categoryFlag = "NoCategory";
}

MuTauSync::~MuTauSync(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "MuTauSync::~MuTauSync Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "MuTauSync::~MuTauSync()" << std::endl;
}

void  MuTauSync::Configure(){
	HToTaumuTauh::Configure();

	// implement changes between main analysis and sync exercise
	cut.at(OppCharge) = 999; // set to 999 to disable opp. charge cut
	title.at(OppCharge) += " DISABLED";
	cut.at(MT) = 999; // set to 999 to disable mt cut
	title.at(MT) += " DISABLED";
}
