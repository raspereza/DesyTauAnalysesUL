#!/bin/bash

CHANNEL=$1
echo "CONFIGFILE,FILELIST" > parameters.txt

if [[ $CHANNEL == "em" ]]; then
    # data
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018A 20
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018B 20
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018C 20
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf MuonEG_Run2018D 20

    # Embedded
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018A 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018B 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018C 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedElMu_Run2018D 4
elif [[ $CHANNEL == "tt" ]]; then
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf Tau_Run2018A 20
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf Tau_Run2018B 20
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf Tau_Run2018C 20
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_data.conf Tau_Run2018D 20

    # Embedded
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedTauTau_Run2018A 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedTauTau_Run2018B 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedTauTau_Run2018C 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_embedded.conf EmbeddedTauTau_Run2018D 4
fi
    
# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DYJetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY1JetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY2JetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY3JetsToLL_M-50 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf DY4JetsToLL_M-50 20

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf W4JetsToLNu 20

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZZTo4L 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf TTToSemiLeptonic 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_t-channel_antitop_4f 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_t-channel_top_4f 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_tW_antitop_5f 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ST_tW_top_5f 10

# H->WW
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf GluGluHToWWTo2L2Nu_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VBFHToWWTo2L2Nu_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf HWminusJ_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf HWplusJ_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZHJ_HToWW_M125 10

# H->tautau
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf GluGluHToTauTau_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf VBFHToTauTau_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WplusHToTauTau_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf WminusHToTauTau_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf ZHToTauTau_M125_13TeV 10


# BBH
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf BBHToTauTau_M125_13TeV 5


# GGH+bb
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf GluGluToBBHToTauTau_M125_13TeV 5

# Interference
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf BBHToTauTauYbYt_M125_13TeV 5


# BBHToWW

./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf BBHToWW_M125_13TeV 5
./split_filelist.sh analysisMacroSynch_${CHANNEL}_18_MC.conf BBHToWWYbYt_M125_13TeV 5
