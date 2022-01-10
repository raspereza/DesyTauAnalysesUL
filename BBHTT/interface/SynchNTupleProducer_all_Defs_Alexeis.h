#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TSystem.h"

#include "TVector3.h"
#include "TMatrix.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/Common/interface/Config.h"
#include "DesyTauAnalyses/Common/interface/AC1B.h"

#include "DesyTauAnalyses/BBHTT/interface/SynchTree.h"
#include "DesyTauAnalyses/BBHTT/interface/SynchGenTree.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "DesyTauAnalyses/Common/interface/functions.h"
#include "DesyTauAnalyses/BBHTT/interface/functionsSynch2017.h"
#include "DesyTauAnalyses/BBHTT/interface/functionsCP.h"

#include "DesyTauAnalyses/BBHTT/interface/leptau_jets.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/Common/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
#include "HiggsCPinTauDecays/IpCorrection/interface/IpCorrection.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/BBH/interface/Systematics.h"
#include "DesyTauAnalyses/BBH/interface/LeptonScaleSys.h"
#include "DesyTauAnalyses/BBH/interface/ZPtWeightSys.h"
#include "DesyTauAnalyses/BBH/interface/TopPtWeightSys.h"
#include "DesyTauAnalyses/BBH/interface/JetEnergyScaleSys.h"
#include "DesyTauAnalyses/BBH/interface/PuppiMETSys.h"
#include "DesyTauAnalyses/BBH/interface/PFMETSys.h"
#include "DesyTauAnalyses/BBH/interface/BtagSys.h"


#include "HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "RooFunctor.h"
