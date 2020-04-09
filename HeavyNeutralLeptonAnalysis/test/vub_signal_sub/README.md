## How to submit on VUB-T2 batch system

You need a yml file containing infos about HNL signal samples.

If you don't have one, you could use `createYMLfile.sh` to generate it:

```bash
#! /bin/bash
# The script generate a .yml file from a list a HNL displaced samples
# Example:
# bash createYMLfile.sh  /pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/ \
#                        HeavyNeutrino_lljj*_e_Dirac_massiveAndCKM_LO \
#                        Fall17_lljj_e

if [ $# -lt 3 ]
then
  echo "Missing arguments!!"
  echo "Usage: $0 [path] [regex] [output]"
  echo "[path]   : path containing the different signal folders"
  echo "[regex]  : regular expression used for matching the signal wanted"
  echo "[output] : output file name (without .yml)"
  exit
fi
```

The script will generate yml.


Then run `python submit_signal_vub.py` that will generate all the bash scripts for the batch system

```bash
usage: submit_signal_vub.py [-h] --ymls [YMLS [YMLS ...]] --outputPath
                            OUTPUTPATH --cfgfile CFGFILE
                            [--masses [{1,2,3,4,5,6,8,10,15,20} [{1,2,3,4,5,6,8,10,15,20} ...]]]
                            [--newIVF]
```
