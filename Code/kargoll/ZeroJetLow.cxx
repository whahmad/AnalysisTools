/*
 * ZeroJetLow.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetLow.h"

ZeroJetLow::ZeroJetLow(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class ZeroJetLow" << std::endl;
	// run ZeroJetLow category
	categoryFlag = "ZeroJetLow";

	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

ZeroJetLow::~ZeroJetLow() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "ZeroJetLow::~ZeroJetLow Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "ZeroJetLow::~ZeroJetLow()" << std::endl;
}

