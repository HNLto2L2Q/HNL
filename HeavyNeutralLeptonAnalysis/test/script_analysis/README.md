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


### Submit on batch system (multiple input files)
Instruction on how to submit HNL2L2Q skim on a batch system using as an input a list of file obtained from crab.
The script is `submit_skimming.py`

1. Compile `CloneTree.C`

```bash
g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
```

2. Modify `sample_skimming_2016.yml` where you specify the path to the different samples. It has three different categories: `data`, `backgrounds`, `signals`.
You can add the sample in the right category and specify the CRAB output path in the field `path`

```yml
data:
  Run2016E_SingleMuon:
    path: '/pnfs/iihe/cms/store/user/moanwar/SingleMuon/crab_Run15_Run2016E-17Jul2018_SingleMuon_trigAfilter_corr/200418_150505/0000/'
    group: data
backgrounds:
  ttbar:
    path: 'PUT_THE_RIGHT_PATH_HERE'
    weight: 0.000005408
    color: 'b'
signals:
  M_mu_M1_74.22:
    path: 'PUT_THE_RIGHT_PATH_HERE'
    ctau: 74.22
    mass: 1
    weight: 0.000333641
```

3. Run `python submit_skim.py` plus options:

Example

```bash
python submit_skimming.py --ymls samples_skimming_2016.yml --groups data --outputPath skimming_2016_RUNE --sourcefile ./CloneTree.exe
```

  - it creates `FARM/inputs`, `FARM/logs` and `skimming_2016_RUNE/` that will contain respectively scripts to be submitted, logs and the ROOT output files.
  - it loops on the __samples__ list in `yml` file that you parse and create and store, in `FARM/inputs` the _.sh_ for each file of the samples.

Options:

```bash
usage: submit_skimming.py [-h] --ymls [YMLS [YMLS ...]] --sourcefile
                          SOURCEFILE --outputPath OUTPUTPATH
                          [--onlyBackground] [--skipSignals]
                          [--groups [{dyjets,singletop,ttV,VVV,VV,data} [{dyjets,singletop,ttV,VVV,VV,data} ...]]]

optional arguments:
-h, --help            show this help message and exit
--ymls [YMLS [YMLS ...]]
yml files
--sourcefile SOURCEFILE
path to source file to be used
--outputPath OUTPUTPATH
output path
--onlyBackground      Submit skimming only for background samples
--skipSignals         Submit skimming skipping signals samples
--groups [{dyjets,singletop,ttV,VVV,VV,data} [{dyjets,singletop,ttV,VVV,VV,data} ...]]
Group of sample to process [dyjets, singletop, ttV,
VVV, VV, data]

```



### Submit on batch system (hadd input files)
Instruction on how to submit HNL2L2Q skim on a batch system.
The script `submit_skim.py` will find if is running on __lxplus__ or __VUB-T2__ and it will create scripts for HTCondor or PBS.

1. Compile `CloneTree.C`

```bash
g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
```

2. Modify `sample.yml` or `sample_lxplus.yml` where you specify the path to the different samples and the file name of each samples

```python
path_sig: "/user/moanwar/heavyneutrino/CMSSW_9_4_13/src/HNL/HeavyNeutralLeptonAnalysis/test/signal_samples_mu"
path_bkg: "/pnfs/iihe/cms/store/user/moanwar/SamplesToSkimm_run5_muon"
path_data: "/pnfs/iihe/cms/store/user/moanwar/SamplesToSkimm_run5_muon"
```

3. Run `python submit_skim.py` plus options:

  - it creates `FARM/inputs`, `FARM/logs` and `FARM/outputs` that will contain respectively scripts to be submitted, logs and the ROOT output files.
  - it loops on the __samples__ list in `yml` file that you parse and create and store, in `FARM/inputs` the _.sh_ for each samples.

- `submit_skim.py` options:

```bash
usage: submit_skim.py [-h] --ymls [YMLS [YMLS ...]] [--onlyBackground]
                      [--skipSignals]
                      [--groups [{dyjets,singletop,ttV,VVV,VV} [{dyjets,singletop,ttV,VVV,VV} ...]]]

optional arguments:
  -h, --help            show this help message and exit
  --ymls [YMLS [YMLS ...]]
                        yml files
  --onlyBackground      Submit skimming only for background samples
  --skipSignals         Submit skimming skipping signals samples
  --groups [{dyjets,singletop,ttV,VVV,VV} [{dyjets,singletop,ttV,VVV,VV} ...]]
                        Group of sample to process [dyjets, singletop, ttV,
                        VVV, VV]
```

- The template bash script is:


```bash
bash_script_tmp = """#!/bin/bash
set -e

FOLDER=\"{folder}\"

{source_command}
pushd {cmssw_src}
eval `scramv1 runtime -sh`
popd

{command}

mv *.root $FOLDER/{output}
"""
```

  - it creates `skimming_submit.sh` (qsub) or `skimming_submit.cmd` (HTCondor) in `FARM/inputs`
  - condor template:

```python
  condor_cfg = """Universe                = vanilla
  Environment             = CONDORJOBID=$(Process)
  requirements            = (OpSysAndVer =?= "SLCern6")
  notification            = Error
  when_to_transfer_output = ON_EXIT
  transfer_output_files   = ""
  should_transfer_files   = YES
  executable = $(filename)
  output = {logsfolder}/$Fn(filename).$(ClusterId).$(ProcId).out
  error  = {logsfolder}/$Fn(filename).$(ClusterId).$(ProcId).err
  log    = {logsfolder}/$Fn(filename).$(ClusterId).$(ProcId).log
  +JobFlavour = "testmatch"
  queue filename matching files {inputsfolder}/{script_wildcard}
  """
```
4. Run
  - __qsub__: `bash FARM/inputs/skimming_submit.sh`
  - __HTCondor__: `condor_submit FARM/inputs/skimming_submit.cmd`

#### Example

```bash
python submit_skim.py --ymls samples_lxplus.yml --onlyBackground --groups singletop
```

where _groups_ are defined for each sample in the yml file and should also be used to aggregate samples under the same macro process. If `--groups` is skipped, the script will loop on all the samples.

```yml
ST_s-channel_4f_leptonDecays:
  filename: 'ST_s-channel_4f_leptonDecays.root'
  weight: 0.000003681
  group: singletop
ST_tW_antitop_5f_inclusiveDecays:
  filename: 'ST_tW_antitop_5f_inclusiveDecays.root'
  weight: 0.000005262
  group: singletop
ST_tW_top_5f_inclusiveDecays:
  filename: 'ST_tW_top_5f_inclusiveDecays.root'
  weight: 0.000005156
  group: singletop
ST_t-channel_antitop_4f_inclusiveDecays:
  filename: 'ST_t-channel_antitop_4f_inclusiveDecays.root'
  weight: 0.000002101
  group: singletop
ST_t-channel_top_4f_inclusiveDecays:
  filename: 'ST_t-channel_top_4f_inclusiveDecays.root'
  weight: 0.000002027
  group: singletop
```
