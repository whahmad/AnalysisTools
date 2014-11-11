#include "ZtoEE.h"

ZtoEE::ZtoEE(TString Name_, TString id_):
	ZtoEMu(Name_,id_)
{
	std::cout << "Setting up the class ZtoEE" << std::endl;
	// run ZtoEE category
	categoryFlag = "ZtoEE";
}

ZtoEE::~ZtoEE()
{
	for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "ZtoEE::~ZtoEE Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	}
	std::cout << "ZtoEE::~ZtoEE()" << std::endl;
}
