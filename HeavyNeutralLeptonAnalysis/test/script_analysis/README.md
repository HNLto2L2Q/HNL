## Submit skim

Compile `CloneTree.C`

```bash
g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
```

### Interactive submit

`CloneTree.C` gets external arguments:

```bash
./CloneTree.exe [path/to/ntuples] [input_filename]  [isMC]
```

Where `[path/to/ntuples]` is the full path to the ntuples location, 
`[input_filename]` is the ROOT file name used as input and 
`[isMC]` is a boolen for switch on/off some part of the code


### __qsub__ submit
Instruction on how to submit HNL2L2Q skim on __VUB T2__.

1. Compile `CloneTree.C`

```bash
g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
```

2. Modify `sample.yml` where you specify the path to the different samples and the file name of each samples

```python
path_sig: "/user/moanwar/heavyneutrino/CMSSW_9_4_13/src/HNL/HeavyNeutralLeptonAnalysis/test/signal_samples_mu"
path_bkg: "/pnfs/iihe/cms/store/user/moanwar/SamplesToSkimm_run5_muon"
path_data: "/pnfs/iihe/cms/store/user/moanwar/SamplesToSkimm_run5_muon"
```

3. Run `python submit_skim.py`:
  - it creates `FARM/inputs`, `FARM/logs` and `FARM/outputs` that will contain respectively scripts to be submitted, logs and the ROOT output files.
  - it loops on the __samples__ list in `samples.yml` and create and store, in `FARM/inputs` the _.sh_ for each samples. The template bash script is:

```bash
bash_script_tmp = """#!/bin/bash
set -e

FOLDER=\"{folder}\"

source $VO_CMS_SW_DIR/cmsset_default.sh
pushd $FOLDER
eval `scramv1 runtime -sh`
popd

# g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit $FOLDER/CloneTree.C -o CloneTree.exe
# cp $FOLDER/CloneTree.exe .

{command}

mv *.root $FOLDER/{output}
"""
```

  - it creates `skimming_submit.sh` in `FARM/inputs`, that is a list of `qsub` commands
4. Run `bash FARM/inputs/skimming_submit.sh`
