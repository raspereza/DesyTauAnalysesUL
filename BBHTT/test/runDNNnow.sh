#!/bin/bash

YEAR=$1
CHANNEL=$2

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters: Run with ./runDNNnow.sh YEAR={2016,2017,2018} CHANNEL={mt,et,emu}"
else
    if [ "$CHANNEL" == "mt" ]; then
	echo "mutau channel"
	./NohupParser.sh SingleMuon $YEAR $CHANNEL 
	./NohupParser.sh EmbeddedMuTau $YEAR $CHANNEL 
    elif [ "$CHANNEL" == "et" ]; then
	echo "etau channel"
	./NohupParser.sh SingleElectron $YEAR $CHANNEL 
	./NohupParser.sh EmbeddedElTau $YEAR $CHANNEL 
    else
	echo "emu channel"
	./NohupParser.sh MuonEG $YEAR $CHANNEL 
	./NohupParser.sh Embedded $YEAR $CHANNEL 
    fi
    ./NohupParser.sh DYJets $YEAR $CHANNEL 
    ./NohupParser.sh WJets $YEAR $CHANNEL 
    ./NohupParser.sh TTbar $YEAR $CHANNEL 
    ./NohupParser.sh Diboson $YEAR $CHANNEL 
    ./NohupParser.sh GluGluHToTauTau $YEAR $CHANNEL 
    ./NohupParser.sh SingleTop $YEAR $CHANNEL
    ./NohupParser.sh VBFHToTauTau $YEAR $CHANNEL 
    ./NohupParser.sh BBH $YEAR $CHANNEL
    ./NohupParser.sh ZHToTauTau $YEAR $CHANNEL
    ./NohupParser.sh WHToTauTau $YEAR $CHANNEL
    ./NohupParser.sh ggHToWW $YEAR $CHANNEL
    ./NohupParser.sh VBFHToWW $YEAR $CHANNEL
    ./NohupParser.sh ZHToWW $YEAR $CHANNEL
    ./NohupParser.sh WHToWW $YEAR $CHANNEL
fi
