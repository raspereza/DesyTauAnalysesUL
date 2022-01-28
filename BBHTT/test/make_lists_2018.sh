#!/bin/sh
dirData=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/data
dirEmbedded=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/emb
dirMC=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/mc
dirMC2=/pnfs/desy.de/cms/tier2/store/user/rasp/ntuples_Dec2020/2018/mc_2
dirMC_UL=/pnfs/desy.de/cms/tier2/store/user/acardini/ntuples/Oktoberfest21/2018/mc

OUTDIR=./2018
if [ ! -d "$OUTDIR" ]; then
  echo "Path does not exist: ${OUTDIR}"
  echo "Please create it"
  exit
fi

ls $dirMC_UL/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_ext/*root > $OUTDIR/DYJetsToLL_M-50
ls $dirMC_UL/DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY1JetsToLL_M-50
ls $dirMC_UL/DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY2JetsToLL_M-50
ls $dirMC_UL/DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY3JetsToLL_M-50
ls $dirMC_UL/DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/DY4JetsToLL_M-50 

ls $dirMC_UL/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/WJetsToLNu
ls $dirMC_UL/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W1JetsToLNu
ls $dirMC_UL/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W2JetsToLNu
ls $dirMC_UL/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W3JetsToLNu
ls $dirMC_UL/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/*root > $OUTDIR/W4JetsToLNu

ls $dirMC_UL/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTTo2L2Nu

# list is too long (need splitting) -> 
ls $dirMC_UL/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*0.root > $OUTDIR/TTToSemiLeptonic
for index in {1..9}
do
    ls $dirMC_UL/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/*${index}.root >> $OUTDIR/TTToSemiLeptonic
done


ls $dirMC_UL/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/TTToHadronic

ls $dirMC_UL/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > $OUTDIR/ST_t-channel_antitop_4f
ls $dirMC_UL/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/*root > $OUTDIR/ST_t-channel_top_4f
ls $dirMC2/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_antitop_5f
ls $dirMC2/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/*root > $OUTDIR/ST_tW_top_5f

ls $dirMC2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/VVTo2L2Nu
ls $dirMC/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/WZTo2L2Q
ls $dirMC_UL/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/*root > $OUTDIR/WZTo3LNu
ls $dirMC2/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/*root > $OUTDIR/ZZTo2L2Q
ls $dirMC_UL/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/*root > $OUTDIR/ZZTo4L

ls $dirMC_UL/GluGluHToTauTau/*root > $OUTDIR/GluGluHToTauTau_M125
ls $dirMC_UL/VBFHToTauTau/*.root > $OUTDIR/VBFHToTauTau_M125
ls $dirMC_UL/WplusHToTauTau/*.root > $OUTDIR/WplusHToTauTau_M125
ls $dirMC_UL/WminusHToTauTau/*.root > $OUTDIR/WminusHToTauTau_M125
ls $dirMC/ZHToTauTau_M125_13TeV_powheg_pythia8/*root > $OUTDIR/ZHToTauTau_M125_13TeV

ls $dirMC/GluGluHToWWTo2L2Nu/*root > $OUTDIR/GluGluHToWWTo2L2Nu_M125
ls $dirMC/VBFHToWWTo2L2Nu/*root > $OUTDIR/VBFHToWWTo2L2Nu_M125
ls $dirMC/HWminusJ_HToWW/*root > $OUTDIR/HWminusJ_HToWW_M125
ls $dirMC/HWplusJ_HToWW/*root > $OUTDIR/HWplusJ_HToWW_M125
ls $dirMC/HZJ_HToWW/*root > $OUTDIR/ZHJ_HToWW_M125

ls $dirData/MuonEG-Run2018A-UL2018/*.root > $OUTDIR/MuonEG_Run2018A
ls $dirData/MuonEG-Run2018B-UL2018/*.root > $OUTDIR/MuonEG_Run2018B
ls $dirData/MuonEG-Run2018C-UL2018/*.root > $OUTDIR/MuonEG_Run2018C
ls $dirData/MuonEG-Run2018D-UL2018/*.root > $OUTDIR/MuonEG_Run2018D

ls $dirEmbedded/EmbeddingRun2018A_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018A
ls $dirEmbedded/EmbeddingRun2018B_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018B
ls $dirEmbedded/EmbeddingRun2018C_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018C
ls $dirEmbedded/EmbeddingRun2018D_ElMu/*root > $OUTDIR/EmbeddedElMu_Run2018D

ls $dirMC_UL/BBHToTauTauYbYt/*.root > $OUTDIR/BBHToTauTauYbYt_M125_13TeV

ls $dirMC_UL/GluGluToBBHToTauTau/*.root > $OUTDIR/GluGluToBBHToTauTau_M125_13TeV
ls $dirMC_UL/GluGluToBBHToTauTau-ext1/*.root >> $OUTDIR/GluGluToBBHToTauTau_M125_13TeV
ls $dirMC_UL/GluGluToBBHToTauTau-ext2/*.root >> $OUTDIR/GluGluToBBHToTauTau_M125_13TeV

ls $dirMC_UL/BBHToTauTau/*.root > $OUTDIR/BBHToTauTau_M125_13TeV
ls $dirMC_UL/BBHToTauTau-ext1/*.root >> $OUTDIR/BBHToTauTau_M125_13TeV
ls $dirMC_UL/BBHToTauTau-ext2/*.root >> $OUTDIR/BBHToTauTau_M125_13TeV
