#!/bin/bash
#Script to run RECO and DQM sequences on existing files using cmsDriver.py
#More background information: 
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCmsDriver
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookDataFormats

hostname
env
voms-proxy-info

#abort on error
set -e
set -x

#number of events to process per job
#passed through condor, but here is a default value
if [ -z "$PERJOB" ]; then
    PERJOB=200
fi

#
#set default conditions - run3 2021
#CONDITIONS=auto:phase1_2022_realistic ERA=Run3 GEOM=DB.Extended CUSTOM=
#
#conditions - 2018
#CONDITIONS=auto:phase1_2018_realistic ERA=Run2_2018 GEOM=DB.Extended CUSTOM=
#
#conditions - phase2
CONDITIONS=auto:phase2_realistic_T21 ERA=Phase2C17I13M9 GEOM=Extended2026D95

#Running with 2 threads allows to use more memory on grid
NTHREADS=20

#Argument parsing
if [ "$#" -ne 3 ]; then
    echo "Must pass exactly 3 arguments: run_relval.sh [Zee|TTbar|conf] [njob] [cluster]"
    exit 0
fi

#index of the job is used to keep track of which events / files to process in the reco step
NJOB=$(($3 + 1))

#set CMSSW environment and go to condor work dir
LAUNCHDIR=`pwd`
source /cvmfs/cms.cern.ch/cmsset_default.sh

#this environment variable comes from the condor submit script
cd $CMSSW_BASE
eval `scram runtime -sh`

#define HOME if not defined.
if [ -z "$HOME" ]; then
    export HOME=/tmp
fi

#if the _CONDOR_SCRATCH_DIR is not defined, we are not inside a condor batch job
if [ -z "$_CONDOR_SCRATCH_DIR" ]; then
    cd $LAUNCHDIR
else
    cd $_CONDOR_SCRATCH_DIR
fi

##RelVal samples
if [ "$1" == "Zee" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/MtdValidation/test/tmp/das_cache/Zee.txt
    OUTPUT_DIR=root://eosuser.cern.ch//eos/user/a/aperego/Timing/root_files/Zee
    NAME=Zee
elif [ "$1" == "TTbar" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/MtdValidation/test/tmp/das_cache/TTbar.txt
    OUTPUT_DIR=root://eosuser.cern.ch//eos/user/a/aperego/Timing/root_files/TTbar
    NAME=TTbar	
elif [ "$1" == "conf" ]; then
    INPUT_FILELIST=${CMSSW_BASE}/src/Validation/MtdValidation/test/tmp/das_cache/TTbar.txt
    NAME=conf
else
    echo "Argument 1 must be [Zee|TTbar|conf] but was $1"
    exit 1
fi

#skip njob*perjob events
SKIPEVENTS=$(($NJOB * $PERJOB))

#Just print out environment last time for debugging
echo $INPUT_FILELIST $NAME $STEP $SKIPEVENTS
#env

if [ $NAME == "conf" ]; then
    mkdir -p $NAME
    cd $NAME

    FILENAME=`sed -n "${NJOB}p" $INPUT_FILELIST`
    echo "FILENAME="$FILENAME

    cmsDriver.py step3 --conditions $CONDITIONS -s RAW2DIGI,L1Reco,RECOSIM,PAT --datatier MINIAODSIM --nThreads $NTHREADS -n -1 --era $ERA --eventcontent MINIAODSIM --geometry=$GEOM --filein step2.root --fileout file:step3_inMINIAODSIM.root --no_exec --python_filename=step3.py $CUSTOM
    
else
    
    #Start of workflow
    echo "Making subdirectory $NAME"

    if [ -e $NAME ]; then
        echo "directory $NAME exists, aborting"
        exit 1
    fi

    mkdir $NAME
    cd $NAME

    FILENAME=`sed -n "${NJOB}p" $INPUT_FILELIST`
    echo "FILENAME="$FILENAME
    #Run the actual CMS reco 
    echo "Running step RECO" 
    cmsDriver.py step3 -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions $CONDITIONS --datatier MINIAODSIM,DQMIO  --nThreads $NTHREADS -n -1 --era $ERA --eventcontent MINIAODSIM,DQM --no_exec --geometry=$GEOM --filein $FILENAME --fileout file:step3_inMINIAODSIM.root | tee step3.log  2>&1
    sed -i '$a process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root"))' step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py
    cmsRun -n $NTHREADS step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py
    xrdcp step3_inMINIAODSIM.root $OUTPUT_DIR/step3_inMINIAODSIM_$4_$3.root
    xrdcp step3_inMINIAODSIM_inDQM.root $OUTPUT_DIR/step3_inMINIAODSIM_inDQM_$4_$3.root
    xrdcp step3.log $OUTPUT_DIR/step3_$4_$3.log
    xrdcp histo.root $OUTPUT_DIR/histo_$4_$3.root
    # hadd ttree_Zee.root histo_$4_*.root
 
    cd ..
fi
