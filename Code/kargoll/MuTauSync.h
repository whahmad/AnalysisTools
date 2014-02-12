#ifndef MuTauSync_h
#define MuTauSync_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "HToTaumuTauh.h"

class MuTauSync: public HToTaumuTauh {

public:
	MuTauSync(TString Name_, TString id_);
	virtual ~MuTauSync();

	virtual void Configure();

};
#endif
