# How to skim ntuples

## SLURM 
You should have an account on _ingrid_. 

Load `SLURMUTILS`

```bash
module load slurm/slurm_utils
```

1. Generate a text file containing the absoulute path to your datasets
2. Compile `CloneTree_new.C`:
```cpp
g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree_new.C -o CloneTree_new.exe
```
3. Modify `slurm_SKIMMER_config.py` specifing 
    * the path to text files `inputFiles = glob.glob` (wildcards friendly)
    * the number of files to be processed per job `interval = 30`
    * the `config.stageoutDir` 
4. Submit with `slurm_submit -s slurm_SKIMMER_config.py`

