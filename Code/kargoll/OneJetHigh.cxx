/*
 * OneJetHigh.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "OneJetHigh.h"

OneJetHigh::OneJetHigh(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class OneJetHigh" << std::endl;
	// run OneJetHigh category
	categoryFlag = "OneJetHigh";

	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

OneJetHigh::~OneJetHigh() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "OneJetHigh::~OneJetHigh Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "OneJetHigh::~OneJetHigh()" << std::endl;
}

