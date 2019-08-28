#!/bin/bash
source logincms_cvmfs_slc6.sh

export PATH=path:$PATH
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
export PATH=${ROOTSYS}/bin:${PATH}
export SCRAM_ARCH=slc7_amd64_gcc820

#echo $PWD

cmsswlibdir=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw-patch/CMSSW_9_4_13_patch4/lib/slc7_amd64_gcc820
cmsswincdir=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw-patch/CMSSW_9_4_13_patch4/src

#echo $PWD

macrosdir=/lustre/home/taliercio/SL7/CMSSW_9_4_13_patch4/src/condor/new_skimming
#tarname=CMSSW_9_4_13_patch4

cd ${macrosdir}

eval `scramv1 runtime -sh`

mkdir -p dir_skimmed
mkdir -p dir_jobs
mkdir -p dir_list

#savedir=$macrosdir/dir_skimmed
#savedir_jobs=$macrosdir/dir_jobs
#savedir_list=$macrosdir/dir_list

g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${cmsswincdir} CloneTree_new.C `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib -L ${cmsswlibdir} -o CloneTree_new.exe

./CloneTree_new.exe input_root

mv skimmed_sample*.root dir_skimmed #${savedir} #dir_skimmed
mv skimmer_sample*.sh dir_jobs #${savedir_jobs} #dir_jobs
mv condor_ex_new_sample*.cfg  dir_jobs #${savedir_jobs} #dir_jobs
mv list_sample*.txt dir_list #dir_list #${savedir_list} #dir_list
