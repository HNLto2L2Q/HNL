#!/bin/bash

savedir=${PWD}

echo $savedir

directory=/lustre/cms/store/user/atalierc/HNL_2017_mc_good/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_MLML_ext/190308_131226/0000/
output=DYJetsToLL_M-50

if [ -f $output.txt ]; then
rm $output.txt
fi

cd $directory

if [ -f $output.txt ]; then
rm $output.txt
fi

ls -1 *.root > fake.txt 

while read line 
do
echo $directory$line >> $output.txt
done < fake.txt

rm fake.txt

mv $output.txt $savedir

