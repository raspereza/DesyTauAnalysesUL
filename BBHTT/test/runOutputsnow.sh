#!/bin/bash

OUTDIR=$1

SAMPLES="ggHww VBFHww TTbar BBH Diboson DYJets ZHTT VBFHTT ggHTT ZHww QCD ST WHww WJets WHTT"


if [ "$#" -eq 1 ]; then
    for SAMPLE in $SAMPLES; do
	echo $SAMPLE
	./NohupParserOutput.sh $SAMPLE $OUTDIR
    done
elif [ "$#" -eq 0 ]; then
    for SAMPLE in $SAMPLES; do
	echo $SAMPLE
	./NohupParserOutput.sh $SAMPLE
    done
else
    echo "Illegal number of arguments: run either with default output location with no arguments, or with the output directory as ./runOutputsnow.sh /nfs/dust/cms/user/cardinia/test_output"
fi
    
