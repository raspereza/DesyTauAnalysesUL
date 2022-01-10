export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_26
cd CMSSW_10_6_26/src
cmsenv

git cms-init

cd ${CMSSW_BASE}/src
git clone https://github.com/DesyTau/DesyTauAnalysesUL DesyTauAnalyses

cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/HiggsCPinTauDecays.git

scram b -j 8

cd ${CMSSW_BASE}/src
git cms-addpkg RecoEgamma/EgammaTools
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools

scram b -j 8

cd ${CMSSW_BASE}/src

git clone https://github.com/svfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_19_02_2019
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF

cd ${CMSSW_BASE}/src

scram b -j 16
scram b -j 16

