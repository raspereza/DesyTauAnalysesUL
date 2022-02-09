// Abstract class for Systematics evaluation
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef Systematics_h
#define Systematics_h

#include <iostream>

#include "TTree.h"
#include "TString.h"

#include "DesyTauAnalyses/Common/interface/functions.h"
#include "DesyTauAnalyses/BBHTT/interface/SynchTree.h"

class Systematics {
public:
    
  Systematics(){};
  virtual ~Systematics(){};
  
  virtual void Eval(utils::channel ch = utils::UNKNOWN) = 0;
  virtual void Write(const char *name="", Int_t option=0) = 0;
protected:
  virtual void Init(SynchTree*) = 0;
  
  TString label;
  SynchTree* cenTree;
  
  void Add(void* pvar, const char* name, const char* shift, const char* type){
    cenTree->fChain->Branch(TString(name)+"_"+label+shift, pvar, TString(name)+"_"+label+shift+"/"+type);
  };
};
  
#endif //!endif Systematics_h
