#!/bin/sh 
# $1 - sample


cat > $1_OutputRun.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
cd /afs/desy.de/user/m/makou/NN_BBH_Analysis
source venv/bin/activate
cd -
if [ "$#" -ne 2]; then
    python XGB_predict.py -s $1
else
    python XGB_predict.py -s $1 -o $2
fi
EOF1

chmod u+x $1_OutputRun.sh
nohup ./$1_OutputRun.sh > nohup_$1_OutputRun.out &
