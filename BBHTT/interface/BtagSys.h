// Lapton Scale Systematics evaluator
// Author: Francesco Costanza (francesco.costanza@cern.ch)

#ifndef BtagSys_h
#define BtagSys_h

#define addvar(name, value, key, type) name[key] = value; this->Add( &this->name[key], #name, key, #type)
#include "DesyTauAnalyses/Common/interface/Config.h"
#include "DesyTauAnalyses/Common/interface/AC1B.h"
#include <DesyTauAnalyses/Common/interface/functions.h>
#include <DesyTauAnalyses/BBHTT/interface/jets_bbh.h>
#include <DesyTauAnalyses/BBHTT/interface/Systematics.h>

using namespace utils;

class BtagSys : public Systematics {
public:
  
  BtagSys(){};
  
  BtagSys(SynchTree* c, TString name){
    cenTree = c;
    label = "CMS_eff_b_13TeV";    
    if (name=="Mistag")
      label = "CMS_mistag_b_13TeV";
    this->SetUncertaintyName(name);
    this->Init(cenTree);
  };
  
  virtual ~BtagSys(){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      delete it->second;
  };
  
  virtual void Eval(utils::channel ch = utils::UNKNOWN){
    this->Central();
    this->ScaleUp();
    this->ScaleDown();
  };

  virtual void Write(const char *name="", Int_t option=0){
    for (std::map<std::string,TTree*>::iterator it=outTree.begin(); it!=outTree.end(); ++it)
      it->second->Write(name,option);
  };

void SetAC1B(const AC1B * tree){ 
  analysisTree = tree;
}; 

void SetConfig(Config * config){
  cfg = config;
};

void SetBtagScaling(const btag_scaling_inputs * _InputsBtagScaling){
	inputs_btag_scaling = _InputsBtagScaling;
};

void SetUncertaintyName(TString name){
	uncertainty_name = name;
};

protected:

  virtual void Init(SynchTree* c){
    cenTree = c;

    this->InitTree("Up");
    this->InitTree("Down");
  };

  virtual void InitTree(const char* shift){
    std::cout<<label+shift<<std::endl;
    outTree[shift] = cenTree->fChain->CloneTree(0);
    outTree[shift]->SetName(cenTree->fChain->GetName()+TString("_")+label+shift);
    outTree[shift]->SetTitle(cenTree->fChain->GetName()+TString("_")+label+shift);
    outTree[shift]->SetDirectory(cenTree->fChain->GetDirectory());
  };
  
  virtual void Central(){
  };

  virtual void ScaleUp(){
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Up");
    this->Fill("Up");
  };
  
  virtual void ScaleDown(){    
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, uncertainty_name, "Down");
    this->Fill("Down");
  };
  
  virtual void Fill(const char* shift){
    outTree[shift]->Fill();
    jets::counting_jets(analysisTree, cenTree, cfg, inputs_btag_scaling, "central");
  }

  const AC1B * analysisTree;
  Config * cfg;
  const btag_scaling_inputs * inputs_btag_scaling;
  std::map< std::string, TTree* >  outTree;
  TString uncertainty_name;
  
};

#undef addvar

#endif //!endif TopPtWeightSys_h
