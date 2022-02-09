export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_26
cd CMSSW_10_6_26/src
cmsenv

git cms-init

cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/HiggsCPinTauDecays.git

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
git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections

cd ${CMSSW_BASE}/src
git clone https://github.com/marmeyer/RecoilCorrections.git HTT-utilities/RecoilCorrections_KIT

cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/CorrectionsWorkspace HTT-utilities/CorrectionsWorkspace
cd ${CMSSW_BASE}/src/HTT-utilities/CorrectionsWorkspace
root -l -q CrystalBallEfficiency.cxx++
cd ${CMSSW_BASE}/src

git clone https://github.com/veelken/SVfit_standalone.git ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
cd ${CMSSW_BASE}/src/TauAnalysis/SVfitStandalone
git checkout HIG-16-006
cd ${CMSSW_BASE}/src/


cd ${CMSSW_BASE}/src
git clone https://github.com/cardinia/DesyTauAnalysesUL DesyTauAnalyses -b dev


cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.h TauAnalysis/SVfitStandalone/interface/
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/SVfitStandaloneAlgorithm.cc TauAnalysis/SVfitStandalone/src
cp ${CMSSW_BASE}/src/DesyTauAnalyses/patch/SVFit/testSVfitStandalone.cc TauAnalysis/SVfitStandalone/bin
rm TauAnalysis/SVfitStandalone/interface/SVfitStandaloneQuantities.h
rm TauAnalysis/SVfitStandalone/src/SVfitStandaloneQuantities.cc

cd ${CMSSW_BASE}/src

cp -r /nfs/dust/cms/user/cardinia/public/DesyTau_data/data ${CMSSW_BASE}/src/DesyTauAnalyses/Common/.


scram b -j 16
scram b -j 16

