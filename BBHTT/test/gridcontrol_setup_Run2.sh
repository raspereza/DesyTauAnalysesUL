#!/bin/bash

YEAR=20${1}
CHANNEL=$2
if [[ $CHANNEL == "em" ]]; then
    OUTDIR=./emu/$YEAR
else  
    if [[ $CHANNEL == "tt" ]]; then
	OUTDIR=./tautau/$YEAR
    else
	echo "ERROR: please run the script with ./gridcontrol_setup_mt_Run2.sh <year={16,17,18}> <channel={tt,em}>"
	exit
    fi
fi

echo $OUTDIR

if [ ! -d "$OUTDIR" ]; then
  mkdir ${OUTDIR}
  cp ./run_${CHANNEL}_synchntuples.sh $OUTDIR/run_synchntuples.sh
  cp ./split_filelist.sh $OUTDIR/.
  cp ./gc_synch.conf $OUTDIR/. 
  cp ./make_parameter_file_$YEAR.sh $OUTDIR/.
  cp ./add_samples.sh $OUTDIR/.
fi

./make_config_Run2.sh $1 $CHANNEL MC 
./make_config_Run2.sh $1 $CHANNEL data 
./make_config_Run2.sh $1 $CHANNEL embedded 
./make_lists_${YEAR}.sh $CHANNEL

cd ./$OUTDIR
rm parameters.txt
./make_parameter_file_${YEAR}.sh $CHANNEL
cd -

echo "-----------------------------------------"
echo "DONT FORGET TO EDIT gc_synch.conf file!!!"
echo "-----------------------------------------"
